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

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  # use for interactive testing
  group <- "treated"
  # group <- "control"
  # cohort <- "pfizer"
  # matching_round <- as.integer("1")
} else {
  group <- args[[1]]
  
  if (group == "treated") {
    if (length(args) > 1) 
      stop("No additional args to be specified when `group=\"treated\"")
  } else if (group == "control") {
    if (length(args) == 1) {
      stop("`cohort` and `matching_round` must be specified when `group=\"control\"`")
    }
    
    cohort <- args[[2]] # NULL if treated
    matching_round <- as.integer(args[[3]]) # NULL if treated    
    
  }
} 

## get cohort-specific parameters study dates and parameters ---- 
lens <- sapply(study_dates, length)
dates_general <- map(study_dates[lens==1], as.Date)
dates_cohort <- map(study_dates[lens==3], ~map(.x, as.Date))
study_dates <- splice(dates_general, dates_cohort)[names(study_dates)]

if (group == "control") {
  matching_round_date <- study_dates[[cohort]]$control_extract_dates[matching_round]
}

## create output directory ----
if (group == "treated") {
  fs::dir_create(here("output", "pfizer", "treated"))
  fs::dir_create(here("output", "moderna", "treated"))
  fs::dir_create(here("output", "treated", "eligible"))
} else if (group == "control") {
  fs::dir_create(ghere("output", cohort, "matchround{matching_round}", "process"))
} 


# import data ----

# use externally created dummy data if not running in the server
# check variables are as they should be
if(Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")){
  
  # ideally in future this will check column existence and types from metadata,
  # rather than from a cohort-extractor-generated dummy data
  
  if (group == "treated") {
    studydef_path <- here("output", "treated", "extract", "input_treated.feather")
    custom_path <- here("lib", "dummydata", "dummy_treated.feather")
  } else if (group == "control") {
    studydef_path <- ghere("output", cohort, "matchround{matching_round}", "extract", "input_controlpotential.feather")
    custom_path <- here("lib", "dummydata", "dummy_control_potential1.feather")
  }
  
  data_studydef_dummy <- read_feather(studydef_path) %>%
    # because date types are not returned consistently by cohort extractor
    mutate(across(ends_with("_date"), ~ as.Date(.))) %>%
    # because of a bug in cohort extractor -- remove once pulled new version
    mutate(patient_id = as.integer(patient_id))
  
  data_custom_dummy <- read_feather(custom_path) %>%
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

source(here("analysis", "process_functions.R"))

data_processed <- data_extract %>%
  process_jcvi() %>%
  process_demographic() %>%
  process_pre() 

if (group == "treated") {
  data_processed <- data_processed %>%
    process_post()
}

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

data_vax_wide = data_vax %>%
  pivot_wider(
    id_cols= patient_id,
    names_from = c("vax_index"),
    values_from = c("date", "type"),
    names_glue = "covid_vax_{vax_index}_{.value}"
  )

# only add variables corresponding to 4th dose if group = treated
if (group == "treated") {
  vax4_vars <- rlang::quos(
    
    vax4_type = covid_vax_4_type,
    
    vax4_type_descr = fct_case_when(
      vax4_type == "pfizer" ~ "BNT162b2",
      vax4_type == "az" ~ "ChAdOx1",
      vax4_type == "moderna" ~ "Moderna",
      TRUE ~ NA_character_
    ),
    
    vax4_date = covid_vax_4_date,
    
    vax4_day = as.integer(floor((vax4_date - study_dates$index_date))+1), # day 0 is the day before "start_date"
    
    vax4_week = as.integer(floor((vax4_date - study_dates$index_date)/7)+1), # week 1 is days 1-7.
  
  )
} else if (group == "control") {
  vax4_vars <- rlang::quos(
    NULL
  )
}

data_processed <- data_processed %>%
  left_join(data_vax_wide, by ="patient_id") %>%
  mutate(
    vax1_type = covid_vax_1_type,
    vax2_type = covid_vax_2_type,
    vax3_type = covid_vax_3_type,
    
    
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
    
    
    vax12_type_descr = paste0(vax1_type_descr, "-", vax2_type_descr),
    
    vax1_date = covid_vax_1_date,
    vax2_date = covid_vax_2_date,
    vax3_date = covid_vax_3_date,
    
    vax1_day = as.integer(floor((vax1_date - study_dates$index_date))+1), # day 0 is the day before "start_date"
    vax2_day = as.integer(floor((vax2_date - study_dates$index_date))+1), # day 0 is the day before "start_date"
    vax3_day = as.integer(floor((vax3_date - study_dates$index_date))+1), # day 0 is the day before "start_date"
    
    vax1_week = as.integer(floor((vax1_date - study_dates$index_date)/7)+1), # week 1 is days 1-7.
    vax2_week = as.integer(floor((vax2_date - study_dates$index_date)/7)+1), # week 1 is days 1-7.
    vax3_week = as.integer(floor((vax3_date - study_dates$index_date)/7)+1), # week 1 is days 1-7.
    
    !!! vax4_vars
    
  ) %>%
  select(
    -starts_with("covid_vax_"),
  ) 


####################################################################################

if (group == "treated") {
  selection_group <- rlang::quos(
    
    has_expectedvax3type = vax3_type %in% c("pfizer", "moderna"),
    
    has_vaxgap23 = vax3_date >= (vax2_date+17) | is.na(vax3_date), # at least 17 days between second and third vaccinations
    
    vax3_notbeforestartdate = case_when(
      (vax3_type=="pfizer") & (vax3_date < study_dates$pfizer$start_date) ~ FALSE,
      (vax3_type=="moderna") & (vax3_date < study_dates$moderna$start_date) ~ FALSE,
      TRUE ~ TRUE
    ),
    vax3_beforeenddate = case_when(
      (vax3_type=="pfizer") & (vax3_date <= study_dates$pfizer$end_date) & !is.na(vax3_date) ~ TRUE,
      (vax3_type=="moderna") & (vax3_date <= study_dates$moderna$end_date) & !is.na(vax3_date) ~ TRUE,
      TRUE ~ FALSE
    ),
    
    no_recentcovid30 = is.na(anycovid_0_date) | ((vax3_date - anycovid_0_date) > 30),
    
    c5 = c4 & vax3_notbeforestartdate & vax3_beforeenddate & has_expectedvax3type & has_vaxgap23,
    
  )
  
} else if (group == "control") {
  
  selection_group <- rlang::quos(
    
    vax3_notbeforematchingrounddate = case_when(
      is.na(vax3_date) | (vax3_date > matching_round_date) ~ TRUE,
      TRUE ~ FALSE
    ),
    
    no_recentcovid30 = is.na(anycovid_0_date) | ((matching_round_date - anycovid_0_date) > 30),
    
    c5 = c4 & vax3_notbeforematchingrounddate,
    
  )
  
}

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
    
    has_knownvax1 = vax1_type %in% c("pfizer", "az"),
    has_knownvax2 = vax2_type %in% c("pfizer", "az"),
    
    vax12_homologous = vax1_type==vax2_type,
    has_vaxgap12 = vax2_date >= (vax1_date+17), # at least 17 days between first two vaccinations
    
    c0 = vax1_afterfirstvaxdate & vax2_beforelastvaxdate,
    c1 = c0 & has_age & has_sex & has_imd & has_ethnicity & has_region,
    c2 = c1 & has_vaxgap12 & has_knownvax1 & has_knownvax2 & vax12_homologous,
    c3 = c2 & isnot_hscworker,
    c4 = c3 & isnot_carehomeresident & isnot_endoflife & isnot_housebound,
    
    !!! selection_group,
    
    c6 = c5 & no_recentcovid30,
    
    include = (c0 & c1 & c2 & c3 & c4 & c5 & c6),
    
  )

data_eligible <- data_criteria %>%
  filter(include) %>%
  select(patient_id) %>%
  left_join(data_processed, by="patient_id") %>%
  droplevels()

# save cohort-specific datasets ----
if (group == "treated") {
  
  write_rds(data_eligible %>% filter(vax3_type == "pfizer"), 
            here("output", "pfizer", "treated", "data_treatedeligible.rds"),
            compress="gz")
  
  write_rds(data_eligible %>% filter(vax3_type == "moderna"), 
            here("output", "moderna", "treated", "data_treatedeligible.rds"), 
            compress="gz")
  
}

if (group == "control") {
  
  write_rds(data_eligible, 
            ghere("output", cohort, "matchround{matching_round}", "process", "data_controlpotential.rds"),
            compress = "gz")
  
}


# create flowchart (only for treated) ----

if (group == "treated") {
  
  data_flowchart <- data_criteria %>%
    summarise(
      across(matches("^c\\d"), .fns=sum)
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
        crit == "c0" ~ "Aged 18+ with 2nd dose on or before 31 Aug 2021", 
        crit == "c1" ~ "  with no missing demographic information",
        crit == "c2" ~ "  with homologous primary vaccination course of pfizer or AZ and at least 17 days between doses",
        crit == "c3" ~ "  and not a HSC worker",
        crit == "c4" ~ "  and not a care/nursing home resident, end-of-life or housebound",
        crit == "c5" ~ "  and third dose of pfizer or moderna within the study period and at least 17 days after second dose",
        crit == "c6" ~ "  and no evidence of covid in 30 days before third dose",
        TRUE ~ NA_character_
      )
    )
  
  write_rds(data_flowchart, here("output", "treated", "eligible", "flowchart_treatedeligible.rds"))
  
}