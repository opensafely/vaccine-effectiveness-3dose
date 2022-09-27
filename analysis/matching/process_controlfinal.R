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
library('here')
library('glue')
library('arrow')

## import local functions and parameters ---

source(here("analysis", "design.R"))

source(here("lib", "functions", "utility.R"))



# import command-line arguments ----

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
fs::dir_create(ghere("output", cohort, "match"))


# import ----

# use externally created dummy data if not running in the server
# check variables are as they should be
if(Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")){
  
  source(here("analysis", "dummy", "dummydata_controlfinal.R"))
  
  data_studydef_dummy <- read_feather(ghere("output", cohort, "extract", "input_controlfinal.feather")) %>%
    #because date types are not returned consistently by cohort extractor
    mutate(across(ends_with("_date"),  as.Date))
  
  data_custom_dummy <- read_feather(fs::path("lib", "dummydata", glue("dummy_controlfinal_{cohort}.feather")))
  
  
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
  
  data_outcomes <- data_custom_dummy 
  

} else {
  data_outcomes <- read_feather(ghere("output", cohort, "extract", "input_controlfinal.feather")) %>%
    #because date types are not returned consistently by cohort extractor
    mutate(across(ends_with("_date"),  as.Date))
}


data_matchstatus <- read_rds(ghere("output", cohort, "matchround{n_matching_rounds}", "actual", "data_matchstatus_allrounds.rds"))

# import data for treated group and select those who were successfully matched

data_treatedeligible <- read_rds(ghere("output", cohort, "treated", "data_treatedeligible.rds"))

data_treated <- 
  left_join(
    data_matchstatus %>% filter(treated==1L),
    data_treatedeligible,
    by="patient_id"
  ) 


# import final dataset of matched controls, including matching variables
# alternative to this is re-extracting everything in the study definition
data_control <- 
  data_matchstatus %>% filter(treated==0L) %>%
  left_join(
    map_dfr(
      seq_len(n_matching_rounds), 
      ~{read_rds(ghere("output", cohort, glue("matchround", .x), "actual", "data_successful_matchedcontrols.rds"))}
    ) %>% select(-match_id, -trial_date, -treated, -controlistreated_date), # remove to avoid clash with already-stored variables
    by=c("patient_id")
  ) %>%
  # merge with outcomes data
  left_join(
    data_outcomes,
    by=c("patient_id", "match_id", "trial_date")
  ) %>%
  mutate(
    treated=0L
  )

# check final data agrees with matching status

all(data_control$patient_id %in% (data_matchstatus %>% filter(treated==0L) %>% pull(patient_id)))
all((data_matchstatus %>% filter(treated==0L) %>% pull(patient_id)) %in% data_control$patient_id)


# merge treated and control groups
data_matched <-
  bind_rows(
    data_treated,
    data_control
  ) %>%
  # derive some variables
  mutate(

    # earliest covid event after study start
    anycovid_date = pmin(postest_date, covidemergency_date, covidadmitted_date, covidcritcare_date, coviddeath_date, na.rm=TRUE),
    
    noncoviddeath_date = if_else(!is.na(death_date) & is.na(coviddeath_date), death_date, as.Date(NA_character_)),
    
    cause_of_death = fct_case_when(
      !is.na(coviddeath_date) ~ "covid-related",
      is.na(death_date) ~ "not covid-related",
      TRUE ~ NA_character_
    ),
    
  ) 


write_rds(data_matched, here("output", cohort, "match", glue("data_matched.rds")), compress="gz")

# matching status of all treated, eligible people ----

data_treatedeligible_matchstatus <- 
  left_join(
    data_treatedeligible %>% select(patient_id, vax1_date),
    data_matchstatus %>% filter(treated==1L),
    by="patient_id"
  ) %>%
  mutate(
    matched = if_else(is.na(match_id), 0L, 1L),
    treated = if_else(is.na(match_id), 1L, treated),
  )

print(
  glue(
    "all trial dates match vaccination dates for matched, treated people: ",
    data_treatedeligible_matchstatus %>% 
    filter(matched==1L) %>%
    mutate(
      agree = trial_date==vax1_date
    ) %>% pull(agree) %>% all()
  )
)

write_rds(data_treatedeligible_matchstatus, here("output", cohort, "match", "data_treatedeligible_matchstatus.rds"), compress="gz")

    

# Flowchart ----

## FIXME -- to add flowchart entry for all treated people who ended up with a matched control, and all treated people who were also used as a control in an earlier trial


