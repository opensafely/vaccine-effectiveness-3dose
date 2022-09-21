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


## import command-line arguments ----

args <- commandArgs(trailingOnly=TRUE)


if(length(args)==0){
  # use for interactive testing
  removeobjects <- FALSE
  cohort <- "over12"
} else {
  #FIXME replace with actual eventual action variables
  removeobjects <- TRUE
  cohort <- args[[1]]
}

## get cohort-specific parameters study dates and parameters ----

dates <- map(study_dates[[cohort]], as.Date)
params <- study_params[[cohort]]

## create output directory ----
fs::dir_create(here("output", cohort, "treated"))


# import data ----

# use externally created dummy data if not running in the server
# check variables are as they should be
if(Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")){

  # ideally in future this will check column existence and types from metadata,
  # rather than from a cohort-extractor-generated dummy data

  data_studydef_dummy <- read_feather(ghere("output", cohort, "extract", "input_treated.feather")) %>%
    # because date types are not returned consistently by cohort extractor
    mutate(across(ends_with("_date"), ~ as.Date(.))) %>%
    # because of a bug in cohort extractor -- remove once pulled new version
    mutate(patient_id = as.integer(patient_id))

  data_custom_dummy <- read_feather(ghere("lib", "dummydata", "dummy_treated_{cohort}.feather")) %>%
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
  data_extract <- read_feather(ghere("output", cohort, "extract", "input_treated.feather")) %>%
    #because date types are not returned consistently by cohort extractor
    mutate(across(ends_with("_date"),  as.Date))
}


# process data -----

## patient-level info ----

data_processed <- data_extract %>%
  mutate(

    sex = fct_case_when(
      sex == "F" ~ "Female",
      sex == "M" ~ "Male",
      #sex == "I" ~ "Inter-sex",
      #sex == "U" ~ "Unknown",
      TRUE ~ NA_character_
    ),

    # ethnicity_combined = if_else(is.na(ethnicity), ethnicity_6_sus, ethnicity),
    # 
    # ethnicity_combined = fct_case_when(
    #   ethnicity_combined == "1" ~ "White",
    #   ethnicity_combined == "4" ~ "Black",
    #   ethnicity_combined == "3" ~ "South Asian",
    #   ethnicity_combined == "2" ~ "Mixed",
    #   ethnicity_combined == "5" ~ "Other",
    #   #TRUE ~ "Unknown",
    #   TRUE ~ NA_character_
    # 
    # ),

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

    # prior_tests_cat = cut(prior_covid_test_frequency, breaks=c(0, 1, 2, 3, Inf), labels=c("0", "1", "2", "3+"), right=FALSE),

    prior_covid_infection = (!is.na(postest_0_date)) | (!is.na(covidadmitted_0_date)) | (!is.na(primary_care_covid_case_0_date)),

    # latest covid event before study start
    anycovid_0_date = pmax(postest_0_date, covidemergency_0_date, covidadmitted_0_date, na.rm=TRUE),
    
    # # earliest covid event after study start
    # anycovid_1_date = pmin(postest_1_date, covidemergency_1_date, covidadmitted_1_date, covidcc_1_date, coviddeath_date, na.rm=TRUE),
    # 
    # noncoviddeath_date = if_else(!is.na(death_date) & is.na(coviddeath_date), death_date, as.Date(NA_character_)),
    # 
    # cause_of_death = fct_case_when(
    #   !is.na(coviddeath_date) ~ "covid-related",
    #   !is.na(death_date) ~ "not covid-related",
    #   TRUE ~ NA_character_
    # ),

  )

## reshape vaccination data ----

data_vax <- local({

  data_vax_any <- data_processed %>%
    select(patient_id, matches("covid\\_vax\\_any\\_\\d+\\_date")) %>%
    pivot_longer(
      cols = -patient_id,
      names_to = c(NA, "vax_any_index"),
      names_pattern = "^(.*)_(\\d+)_date",
      values_to = "date",
      values_drop_na = TRUE
    ) %>%
    arrange(patient_id, date)
  
  data_vax_pfizerA <- data_processed %>%
    select(patient_id, matches("covid\\_vax\\_pfizerA\\_\\d+\\_date")) %>%
    pivot_longer(
      cols = -patient_id,
      names_to = c(NA, "vax_pfizerA_index"),
      names_pattern = "^(.*)_(\\d+)_date",
      values_to = "date",
      values_drop_na = TRUE
    ) %>%
    arrange(patient_id, date)

  data_vax_pfizerC <- data_processed %>%
    select(patient_id, matches("covid\\_vax\\_pfizerC\\_\\d+\\_date")) %>%
    pivot_longer(
      cols = -patient_id,
      names_to = c(NA, "vax_pfizerC_index"),
      names_pattern = "^(.*)_(\\d+)_date",
      values_to = "date",
      values_drop_na = TRUE
    ) %>%
    arrange(patient_id, date)
  
 
  data_vax <-
    data_vax_any %>%
    full_join(data_vax_pfizerA, by=c("patient_id", "date")) %>%
    full_join(data_vax_pfizerC, by=c("patient_id", "date")) %>%
    mutate(
      type = fct_case_when(
        is.na(vax_pfizerC_index) & (!is.na(vax_pfizerA_index)) ~ "pfizerA",
        (!is.na(vax_pfizerC_index)) & is.na(vax_pfizerA_index)  ~ "pfizerC",
        !is.na(vax_any_index) ~ "other",
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
    
    vax1_type_descr = fct_case_when(
      vax1_type == "pfizerA" ~ "BNT162b2 30micrograms/0.3ml",
      vax1_type == "pfizerC" ~ "BNT162b2 10mcg/0.2ml",
      vax1_type == "any" ~ "Other",
      TRUE ~ NA_character_
    ),
    vax2_type_descr = fct_case_when(
      vax2_type == "pfizerA" ~ "BNT162b2 30micrograms/0.3ml",
      vax2_type == "pfizerC" ~ "BNT162b2 10mcg/0.2ml",
      vax2_type == "any" ~ "Other",
      TRUE ~ NA_character_
    ),
    
    vax1_date = covid_vax_1_date,
    vax2_date = covid_vax_2_date,
   
) %>%
select(
  -starts_with("covid_vax_"),
)



# apply eligibility criteria ----

## define criteria ----

# Define selection criteria ----
data_criteria <- data_processed %>%
  transmute(
    patient_id,
    has_age = !is.na(age),
    has_sex = !is.na(sex) & !(sex %in% c("I", "U")),
    has_imd = imd_Q5 != "Unknown",
    #has_ethnicity = !is.na(ethnicity_combined),
    has_region = !is.na(region),
    vax1_betweenentrydates = case_when(
      (vax1_type==params$treatment) & 
        (vax1_date >= dates$start_date) & 
        (vax1_date <= dates$end_date) ~ TRUE,
      TRUE ~ FALSE
    ),
    has_vaxgap12 = vax2_date >= (vax1_date+17) | is.na(vax2_date), # at least 17 days between first two vaccinations
    no_recentcovid90 = is.na(anycovid_0_date) |  ((vax1_date - anycovid_0_date)>90),

    include = (
      vax1_betweenentrydates & has_vaxgap12  &
        has_age & has_sex & has_imd & # has_ethnicity &
        has_region &
        no_recentcovid90
    ),
  )

## filter and export ----

data_treated_eligible <- 
  data_criteria %>%
  filter(include) %>%
  select(patient_id) %>%
  left_join(data_processed, by="patient_id") %>%
  droplevels()

write_rds(data_treated_eligible, ghere("output", cohort, "treated", "data_treatedeligible.rds"), compress="gz")


# create flowchart ----


data_flowchart <- data_criteria %>%
  transmute(
    c0 = vax1_betweenentrydates & has_vaxgap12,
    c1 = c0 & (has_age & has_sex & has_imd & has_region),
    c2 = c1 + no_recentcovid90
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
      crit == "c0" ~ "Received age-correct vaccine within study entry dates", 
      crit == "c1" ~ "  with no missing demographic information",
      crit == "c2" ~ "  with no COVID-19 90 days prior",
      TRUE ~ NA_character_
    )
  )
write_rds(data_flowchart, ghere("output", cohort, "treated", "flowchart_treatedeligible.rds"))


