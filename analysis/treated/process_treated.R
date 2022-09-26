######################################

# This script:
# imports data extracted by the cohort extractor (or dummy data)
# fills in unknown ethnicity from GP records with ethnicity from SUS (secondary care)
# tidies missing values
# standardises some variables (eg convert to factor) and derives some new ones
# organises vaccination date data to "vax X type", "vax X date" (rather than "pfizer X date", "az X date", ...)
######################################

# Preliminaries ----

## Import libraries ----
library('tidyverse')
library('lubridate')
library('arrow')
library('here')
library('glue')

## import local functions and parameters ---

source(here("analysis", "design.R"))

source(here("lib", "functions", "utility.R"))


## get cohort-specific parameters study dates and parameters ---- 
lens <- sapply(study_dates, length)
dates_general <- map(study_dates[lens==1], as.Date)
dates_cohort <- map(study_dates[lens==3], ~map(.x, as.Date))
study_dates <- splice(dates_general, dates_cohort)[names(study_dates)]

## create output directory ----
fs::dir_create(here("output", "treated", "eligible"))
fs::dir_create(here("output", "pfizer", "treated"))
fs::dir_create(here("output", "moderna", "treated"))


# import data ----

# use externally created dummy data if not running in the server
# check variables are as they should be
if(Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")){

  # ideally in future this will check column existence and types from metadata,
  # rather than from a cohort-extractor-generated dummy data

  data_studydef_dummy <- read_feather(here("output", "treated", "extract", "input_treated.feather")) %>%
    # because date types are not returned consistently by cohort extractor
    mutate(across(ends_with("_date"), ~ as.Date(.))) %>%
    # because of a bug in cohort extractor -- remove once pulled new version
    mutate(patient_id = as.integer(patient_id))

  data_custom_dummy <- read_feather(ghere("lib", "dummydata", "dummy_treated.feather")) %>%
    mutate(
      msoa = sample(factor(c("1", "2")), size=n(), replace=TRUE) # override msoa so matching success more likely
    )


  not_in_studydef <- names(data_custom_dummy)[!( names(data_custom_dummy) %in% names(data_studydef_dummy) )]
  not_in_custom  <- names(data_studydef_dummy)[!( names(data_studydef_dummy) %in% names(data_custom_dummy) )]


  if(length(not_in_custom)!=0) stop(
    paste(
      "These variables are in studydef but not in custom: ",
      paste(not_in_custom, collapse=", ")
    )
  )

  if(length(not_in_studydef)!=0) stop(
    paste(
      "These variables are in custom but not in studydef: ",
      paste(not_in_studydef, collapse=", ")
    )
  )

  # reorder columns
  data_studydef_dummy <- data_studydef_dummy[,names(data_custom_dummy)]

  unmatched_types <- cbind(
    map_chr(data_studydef_dummy, ~paste(class(.), collapse=", ")),
    map_chr(data_custom_dummy, ~paste(class(.), collapse=", "))
  )[ (map_chr(data_studydef_dummy, ~paste(class(.), collapse=", ")) != map_chr(data_custom_dummy, ~paste(class(.), collapse=", ")) ), ] %>%
    as.data.frame() %>% rownames_to_column()


  if(nrow(unmatched_types)>0) stop(
    #unmatched_types
    "inconsistent typing in studydef : dummy dataset\n",
    apply(unmatched_types, 1, function(row) paste(paste(row, collapse=" : "), "\n"))
  )

  data_extract <- data_custom_dummy 
} else {
  data_extract <- read_feather(ghere("output", "treated", "extract", "input_treated.feather")) %>%
    #because date types are not returned consistently by cohort extractor
    mutate(across(ends_with("_date"),  as.Date))
}


# process data -----

## patient-level info ----

data_processed <- data_extract %>%
  mutate(
    
    ageband = cut(
      age,
      breaks=c(-Inf, 18, 40, 50, 60, 70, 80, 90, Inf),
      labels=c("under 18", "18-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90+"),
      right=FALSE
    ),
    
    sex = fct_case_when(
      sex == "F" ~ "Female",
      sex == "M" ~ "Male",
      #sex == "I" ~ "Inter-sex",
      #sex == "U" ~ "Unknown",
      TRUE ~ NA_character_
    ),
    
    ethnicity_combined = if_else(is.na(ethnicity), ethnicity_6_sus, ethnicity),
    
    ethnicity_combined = fct_case_when(
      ethnicity_combined == "1" ~ "White",
      ethnicity_combined == "4" ~ "Black",
      ethnicity_combined == "3" ~ "South Asian",
      ethnicity_combined == "2" ~ "Mixed",
      ethnicity_combined == "5" ~ "Other",
      TRUE ~ NA_character_
      
    ),
    
    region = fct_collapse(
      region,
      `East of England` = "East",
      `London` = "London",
      `Midlands` = c("West Midlands", "East Midlands"),
      `North East and Yorkshire` = c("Yorkshire and The Humber", "North East"),
      `North West` = "North West",
      `South East` = "South East",
      `South West` = "South West"
    ),
    
    imd_Q5 = factor(imd_Q5, levels = c("1 (most deprived)", "2", "3", "4", "5 (least deprived)", "Unknown")),
    
    rural_urban_group = fct_case_when(
      rural_urban %in% c(1,2) ~ "Urban conurbation",
      rural_urban %in% c(3,4) ~ "Urban city or town",
      rural_urban %in% c(5,6,7,8) ~ "Rural town or village",
      TRUE ~ NA_character_
    ),
    
    care_home_combined = care_home_tpp | care_home_code, # any carehome flag
    
    # clinically at-risk group
    cv = immunosuppressed | chronic_kidney_disease | chronic_resp_disease | diabetes | chronic_liver_disease |
      chronic_neuro_disease | chronic_heart_disease | asplenia | learndis | sev_mental,
    
    cev_cv = fct_case_when(
      cev ~ "Clinically extremely vulnerable",
      cv ~ "Clinically at-risk",
      TRUE ~ "Not clinically at-risk"
    ) %>% fct_rev(),
    
    multimorb =
      (sev_obesity) +
      (chronic_heart_disease) +
      (chronic_kidney_disease)+
      (diabetes) +
      (chronic_liver_disease)+
      (chronic_resp_disease | asthma)+
      (chronic_neuro_disease)#+
    #(learndis)+
    #(sev_mental),
    ,
    multimorb = cut(multimorb, breaks = c(0, 1, 2, Inf), labels=c("0", "1", "2+"), right=FALSE),
    immuno = immunosuppressed | asplenia,
    
    
    # original priority groups https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/1007737/Greenbook_chapter_14a_30July2021.pdf#page=15
    # new priority groups https://www.england.nhs.uk/coronavirus/wp-content/uploads/sites/52/2021/07/C1327-covid-19-vaccination-autumn-winter-phase-3-planning.pdf
    # group 10 split into 16-39 and 40-49 because of earlier roll-out in 40+ from 15 Nov https://www.gov.uk/government/news/jcvi-issues-advice-on-covid-19-booster-vaccines-for-those-aged-40-to-49-and-second-doses-for-16-to-17-year-olds
    
    jcvi_ageband = cut(
      age_aug2021,
      breaks=c(-Inf, 18, 40, 50, 55, 60, 65, 70, 75, 80, Inf),
      labels=c("under 18", "18-39", "40-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80+"),
      right=FALSE
    ),
    
    
    jcvi_group = fct_case_when(
      care_home_combined | hscworker  ~ "1",
      age_aug2021>=80 ~ "2",
      age_aug2021>=75 ~ "3",
      age_aug2021>=70 | (cev & (age_aug2021>=16)) ~ "4",
      age_aug2021>=65 ~ "5",
      between(age_aug2021, 16, 64.999) & cv ~ "6",
      age_aug2021>=60 ~ "7",
      age_aug2021>=55 ~ "8",
      age_aug2021>=50 ~ "9",
      age_aug2021>=40 ~ "10a",
      TRUE ~ "10b"
    ),
    
    jcvi_group_descr = fct_recode(
      jcvi_group,
      "Care home residents and health and social care workers"="1",
      "80+ years"="2",
      "75-79 years"="3",
      "70-74 years or clinically extremely vulnerable"="4",
      "65-69 years"="5",
      "16-64 years or clinically at-risk"="6",
      "60-64 years"="7",
      "55-59 years"="8",
      "50-54 years"="9",
      "40-49 years"="10a",
      "16-39 years"="10b"
    ),
    
    
    prior_tests_cat = cut(prior_covid_test_frequency, breaks=c(0, 1, 2, 3, Inf), labels=c("0", "1", "2", "3+"), right=FALSE),
    
    prior_covid_infection0 = (!is.na(positive_test_0_date)) | (!is.na(admitted_covid_0_date)) | (!is.na(primary_care_covid_case_0_date)),
    
    
    #covidemergency_1_date = pmin(covidemergency_1_date, covidadmitted_1_date, na.rm=TRUE),
    
    # because this value is returned as a factor by the study definition
    admitted_covid_ccdays_1 = as.numeric(as.character(admitted_covid_ccdays_1)),
    admitted_covid_ccdays_2 = as.numeric(as.character(admitted_covid_ccdays_2)),
    admitted_covid_ccdays_3 = as.numeric(as.character(admitted_covid_ccdays_3)),
    admitted_covid_ccdays_4 = as.numeric(as.character(admitted_covid_ccdays_4)),
    
    
    covidcc_1_date = case_when(
      admitted_covid_ccdays_1 > 0 ~ admitted_covid_1_date,
      admitted_covid_ccdays_2 > 0 ~ admitted_covid_2_date,
      admitted_covid_ccdays_3 > 0 ~ admitted_covid_3_date,
      admitted_covid_ccdays_4 > 0 ~ admitted_covid_4_date,
      TRUE ~ as.Date(NA_character_)
    ),
    
    covidcc_2_date = case_when(
      (admitted_covid_ccdays_2 > 0) & (admitted_covid_2_date > covidcc_1_date) ~ admitted_covid_2_date,
      (admitted_covid_ccdays_3 > 0) & (admitted_covid_3_date > covidcc_1_date) ~ admitted_covid_3_date,
      (admitted_covid_ccdays_4 > 0) & (admitted_covid_4_date > covidcc_1_date) ~ admitted_covid_4_date,
      TRUE ~ as.Date(NA_character_)
    ),
    
    covidcc_3_date = case_when(
      (admitted_covid_ccdays_3 > 0) & (admitted_covid_3_date > covidcc_2_date) ~ admitted_covid_3_date,
      (admitted_covid_ccdays_4 > 0) & (admitted_covid_4_date > covidcc_2_date) ~ admitted_covid_4_date,
      TRUE ~ as.Date(NA_character_)
    ),
    
    covidcc_4_date = case_when(
      (admitted_covid_ccdays_4 > 0) & (admitted_covid_4_date > covidcc_3_date) ~ admitted_covid_4_date,
      TRUE ~ as.Date(NA_character_)
    ),
    
    # latest covid event before study start
    anycovid_0_date = pmax(positive_test_0_date, covidemergency_0_date, admitted_covid_0_date, na.rm=TRUE),
    
    # earliest covid event after study start
    anycovid_1_date = pmin(positive_test_1_date, covidemergency_1_date, admitted_covid_1_date, covidcc_1_date, coviddeath_date, na.rm=TRUE),
    
    noncoviddeath_date = if_else(!is.na(death_date) & is.na(coviddeath_date), death_date, as.Date(NA_character_)),
    
    
    cause_of_death = fct_case_when(
      !is.na(coviddeath_date) ~ "covid-related",
      !is.na(death_date) ~ "not covid-related",
      TRUE ~ NA_character_
    ),
    
    
  )

## reshape vaccination data ----

data_vax <- local({
  
  data_vax_pfizer <- data_processed %>%
    select(patient_id, matches("covid\\_vax\\_pfizer\\_\\d+\\_date")) %>%
    pivot_longer(
      cols = -patient_id,
      names_to = c(NA, "vax_pfizer_index"),
      names_pattern = "^(.*)_(\\d+)_date",
      values_to = "date",
      values_drop_na = TRUE
    ) %>%
    arrange(patient_id, date)
  
  data_vax_az <- data_processed %>%
    select(patient_id, matches("covid\\_vax\\_az\\_\\d+\\_date")) %>%
    pivot_longer(
      cols = -patient_id,
      names_to = c(NA, "vax_az_index"),
      names_pattern = "^(.*)_(\\d+)_date",
      values_to = "date",
      values_drop_na = TRUE
    ) %>%
    arrange(patient_id, date)
  
  data_vax_moderna <- data_processed %>%
    select(patient_id, matches("covid\\_vax\\_moderna\\_\\d+\\_date")) %>%
    pivot_longer(
      cols = -patient_id,
      names_to = c(NA, "vax_moderna_index"),
      names_pattern = "^(.*)_(\\d+)_date",
      values_to = "date",
      values_drop_na = TRUE
    ) %>%
    arrange(patient_id, date)
  
  
  data_vax <-
    data_vax_pfizer %>%
    full_join(data_vax_az, by=c("patient_id", "date")) %>%
    full_join(data_vax_moderna, by=c("patient_id", "date")) %>%
    mutate(
      type = fct_case_when(
        (!is.na(vax_az_index)) & is.na(vax_pfizer_index) & is.na(vax_moderna_index) ~ "az",
        is.na(vax_az_index) & (!is.na(vax_pfizer_index)) & is.na(vax_moderna_index) ~ "pfizer",
        is.na(vax_az_index) & is.na(vax_pfizer_index) & (!is.na(vax_moderna_index)) ~ "moderna",
        (!is.na(vax_az_index)) + (!is.na(vax_pfizer_index)) + (!is.na(vax_moderna_index)) > 1 ~ "duplicate",
        TRUE ~ NA_character_
      )
    ) %>%
    arrange(patient_id, date) %>%
    group_by(patient_id) %>%
    mutate(
      vax_index=row_number()
    ) %>%
    ungroup()
  
  data_vax
  
})

# write_rds(data_vax, here("output", "vax", "data_vaxlong.rds"), compress="gz")

data_vax_wide = data_vax %>%
  pivot_wider(
    id_cols= patient_id,
    names_from = c("vax_index"),
    values_from = c("date", "type"),
    names_glue = "covid_vax_{vax_index}_{.value}"
  )

data_processed <- data_processed %>%
  left_join(data_vax_wide, by ="patient_id") %>%
  mutate(
    vax1_type = covid_vax_1_type,
    vax2_type = covid_vax_2_type,
    vax3_type = covid_vax_3_type,
    vax4_type = covid_vax_4_type,
    
    vax12_type = paste0(vax1_type, "-", vax2_type),
    
    
    
    vax1_type_descr = fct_case_when(
      vax1_type == "pfizer" ~ "BNT162b2",
      vax1_type == "az" ~ "ChAdOx1",
      vax1_type == "moderna" ~ "Moderna",
      TRUE ~ NA_character_
    ),
    vax2_type_descr = fct_case_when(
      vax2_type == "pfizer" ~ "BNT162b2",
      vax2_type == "az" ~ "ChAdOx1",
      vax2_type == "moderna" ~ "Moderna",
      TRUE ~ NA_character_
    ),
    vax3_type_descr = fct_case_when(
      vax3_type == "pfizer" ~ "BNT162b2",
      vax3_type == "az" ~ "ChAdOx1",
      vax3_type == "moderna" ~ "Moderna",
      TRUE ~ NA_character_
    ),
    vax4_type_descr = fct_case_when(
      vax4_type == "pfizer" ~ "BNT162b2",
      vax4_type == "az" ~ "ChAdOx1",
      vax4_type == "moderna" ~ "Moderna",
      TRUE ~ NA_character_
    ),
    
    vax12_type_descr = paste0(vax1_type_descr, "-", vax2_type_descr),
    
    vax1_date = covid_vax_1_date,
    vax2_date = covid_vax_2_date,
    vax3_date = covid_vax_3_date,
    vax4_date = covid_vax_4_date,
    vax1_day = as.integer(floor((vax1_date - study_dates$index_date))+1), # day 0 is the day before "start_date"
    vax2_day = as.integer(floor((vax2_date - study_dates$index_date))+1), # day 0 is the day before "start_date"
    vax3_day = as.integer(floor((vax3_date - study_dates$index_date))+1), # day 0 is the day before "start_date"
    vax4_day = as.integer(floor((vax4_date - study_dates$index_date))+1), # day 0 is the day before "start_date"
    vax1_week = as.integer(floor((vax1_date - study_dates$index_date)/7)+1), # week 1 is days 1-7.
    vax2_week = as.integer(floor((vax2_date - study_dates$index_date)/7)+1), # week 1 is days 1-7.
    vax3_week = as.integer(floor((vax3_date - study_dates$index_date)/7)+1), # week 1 is days 1-7.
    vax4_week = as.integer(floor((vax4_date - study_dates$index_date)/7)+1), # week 1 is days 1-7.
  ) %>%
  select(
    -starts_with("covid_vax_"),
  )

####################################################################################

# Define selection criteria ----
data_criteria <- data_processed %>%
  transmute(
    patient_id,
    has_age = !is.na(age),
    has_sex = !is.na(sex),
    has_imd = imd_Q5 != "Unknown",
    has_ethnicity = !is.na(ethnicity_combined),
    has_region = !is.na(region),
    #has_msoa = !is.na(msoa),
    isnot_hscworker = !hscworker,
    isnot_carehomeresident = !care_home_combined,
    isnot_endoflife = !endoflife,
    isnot_housebound = !housebound,
    vax1_afterfirstvaxdate = case_when(
      (vax1_type=="pfizer") & (vax1_date >= study_dates$firstpfizer_date) ~ TRUE,
      (vax1_type=="az") & (vax1_date >= study_dates$firstaz_date) ~ TRUE,
      (vax1_type=="moderna") & (vax1_date >= study_dates$firstmoderna_date) ~ TRUE,
      TRUE ~ FALSE
    ),
    vax2_beforelastvaxdate = !is.na(vax2_date) & (vax2_date <= study_dates$lastvax2_date),
    vax3_notbeforestartdate = case_when(
      (vax3_type=="pfizer") & (vax3_date < study_dates$pfizer$start_date) ~ FALSE,
      #(vax3_type=="az") & (vax1_date >= study_dates$azstart_date) ~ TRUE,
      (vax3_type=="moderna") & (vax3_date < study_dates$moderna$start_date) ~ FALSE,
      TRUE ~ TRUE
    ),
    vax3_beforeenddate = case_when(
      (vax3_type=="pfizer") & (vax3_date <= study_dates$pfizer$end_date) & !is.na(vax3_date) ~ TRUE,
      (vax3_type=="moderna") & (vax3_date <= study_dates$moderna$end_date) & !is.na(vax3_date) ~ TRUE,
      TRUE ~ FALSE
    ),
    vax12_homologous = vax1_type==vax2_type,
    has_vaxgap12 = vax2_date >= (vax1_date+17), # at least 17 days between first two vaccinations
    has_vaxgap23 = vax3_date >= (vax2_date+17) | is.na(vax3_date), # at least 17 days between second and third vaccinations
    has_knownvax1 = vax1_type %in% c("pfizer", "az"),
    has_knownvax2 = vax2_type %in% c("pfizer", "az"),
    has_expectedvax3type = vax3_type %in% c("pfizer", "moderna"),
    
    jcvi_group_6orhigher = jcvi_group %in% as.character(1:6),
    
    include = (
      #jcvi_group_6orhigher & # temporary until more data available
      vax1_afterfirstvaxdate &
        vax2_beforelastvaxdate &
        vax3_notbeforestartdate &
        has_age & has_sex & has_imd & has_ethnicity & has_region &
        has_vaxgap12 & has_vaxgap23 & has_knownvax1 & has_knownvax2 & vax12_homologous &
        isnot_hscworker &
        isnot_carehomeresident & isnot_endoflife &
        isnot_housebound
    ),
  )

data_treated_eligible <- data_criteria %>%
  filter(include) %>%
  select(patient_id) %>%
  left_join(data_processed, by="patient_id") %>%
  droplevels()

# save cohort-specific datasets ----
write_rds(data_treated_eligible %>% filter(vax3_type == "pfizer"), 
          here("output", "pfizer", "treated", "data_treatedeligible.rds"), compress="gz")

write_rds(data_treated_eligible %>% filter(vax3_type == "moderna"), 
          here("output", "moderna", "treated", "data_treatedeligible.rds"), compress="gz")

# create flowchart ----

data_flowchart <- data_criteria %>%
  transmute(
    c0 = vax1_afterfirstvaxdate & vax2_beforelastvaxdate & vax3_notbeforestartdate,
    #c1_1yearfup = c0_all & (has_follow_up_previous_year),
    c1 = c0 & (has_age & has_sex & has_imd & has_ethnicity & has_region),
    c2 = c1 & (has_vaxgap12 & has_vaxgap23 & has_knownvax1 & has_knownvax2 & vax12_homologous),
    c3 = c2 & (isnot_hscworker ),
    c4 = c3 & (isnot_carehomeresident & isnot_endoflife & isnot_housebound),
    c5 = c4 & vax3_beforeenddate & has_expectedvax3type
  ) %>%
  summarise(
    across(.fns=sum)
  ) %>%
  pivot_longer(
    cols=everything(),
    names_to="criteria",
    values_to="n"
  ) %>%
  mutate(
    n_exclude = lag(n) - n,
    pct_exclude = n_exclude/lag(n),
    pct_all = n / first(n),
    pct_step = n / lag(n),
    crit = str_extract(criteria, "^c\\d+"),
    criteria = fct_case_when(
      crit == "c0" ~ "Aged 18+ with 2nd dose on or before 31 Aug 2021", # paste0("Aged 18+\n with 2 doses on or before ", format(study_dates$lastvax2_date, "%d %b %Y")),
      crit == "c1" ~ "  with no missing demographic information",
      crit == "c2" ~ "  with homologous primary vaccination course of pfizer or AZ",
      crit == "c3" ~ "  and not a HSC worker",
      crit == "c4" ~ "  and not a care/nursing home resident, end-of-life or housebound",
      crit == "c5" ~ "  and third dose of pfizer or moderna within the study period",
      TRUE ~ NA_character_
    )
  )

write_rds(data_flowchart, here("output", "treated", "eligible", "flowchart_treatedeligible.rds"))
