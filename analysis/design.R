# # # # # # # # # # # # # # # # # # # # #
# This script:
# creates metadata for aspects of the study design
# # # # # # # # # # # # # # # # # # # # #

# Preliminaries ----

## Import libraries ----
library('tidyverse')
library('here')

## create output directories ----
fs::dir_create(here("lib", "design"))



# define key dates ----

study_dates <- lst(
  
  pfizer = lst( # pfizer dose 3
   start_date = "2021-09-16", #start of recruitment thursday 16 september first pfizer booster jabs administered in england
  ),
  
  moderna = lst( # moderna dose 3
    start_date = "2021-10-29", #start of recruitment friday 29 october first moderna booster jabs administered in england
  ),
  
  # see page 21 of https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/1072064/Vaccine-surveillance-report-week-17.pdf
  # vaccine uptake leveled off in all groups by end of Feb (must earlier in older groups)
  recruitmentend_date = "2022-02-28",
  
  testend_date = "2022-03-31", # last day of public testing in England
  studyend_date = "2022-06-30", # end of available hospitalization data
  
  lastvax2_date = "2021-12-01", # don't recruit anyone with second vaccination after this date
  
  # dose 1 dates
  firstpfizer_date = "2020-12-08", # first pfizer vaccination in national roll-out
  firstaz_date = "2021-01-04", # first az vaccination in national roll-out
  firstmoderna_date = "2021-04-13", # first moderna vaccination in national roll-out
  firstpossiblevax_date = "2020-06-01", # used to catch "real" vaccination dates (eg not 1900-01-01)
)

study_dates <- rapply(study_dates, as.Date, how = "list")

study_dates$mrna$start_date <- min(study_dates$pfizer$start_date, study_dates$moderna$start_date)

cohorts <- c("pfizer", "moderna", "mrna")

extract_increment <- 14

for (c in cohorts) {
  study_dates[[c]]$control_extract_dates = seq(study_dates[[c]]$start_date, study_dates$recruitmentend_date, extract_increment)
}

jsonlite::write_json(study_dates, path = here("lib", "design", "study-dates.json"), auto_unbox=TRUE, pretty =TRUE)

# number of matching rounds to perform for each cohort

n_matching_rounds_list <- sapply(
  cohorts,
  function(x) length(study_dates[[x]][["control_extract_dates"]]),
  USE.NAMES = TRUE
  ) 

# define outcomes ----

events_lookup <- tribble(
  ~event, ~event_var, ~event_descr,
  
  # other
  "anytest", "anytest_date", "SARS-CoV-2 test",
  "symptest", "symptest_date", "SARS-CoV-2 test (symptomatic)",
  "pcrtest", "pcrtest_date", "SARS-CoV-2 test (PCR)",
  "lftest", "lftest_date", "SARS-CoV-2 test (lateral flow)",
  
  # effectiveness
  "postest", "positive_test_date", "Positive SARS-CoV-2 test",
  "covidadmitted", "covidadmitted_date", "COVID-19 hospitalisation",
  "covidcritcare", "covidcritcare_date", "COVID-19 critical care",
  "coviddeath", "coviddeath_date", "COVID-19 death",
  "covidcritcareordeath", "covidcritcareordeath_date", "COVID-19 critical care or death",
  
  # other
  "emergency", "emergency_date", "A&E attendance",
  "covidemergency", "covidemergency_date", "COVID-19 A&E attendance",
  "noncoviddeath", "noncoviddeath_date", "Non-COVID-19 death",
  "death", "death_date", "Any death",
  
)

# define treatments ----

treatement_lookup <-
  tribble(
    ~dose, ~treatment, ~treatment_descr,
    "3","pfizer", "BNT162b2",
    "3", "az", "ChAdOx1-S",
    "3", "moderna", "mRNA-1273",
    "primary", "pfizer-pfizer", "BNT162b2",
    "primary", "az-az", "ChAdOx1-S",
    "primary", "moderna-moderna", "mRNA-1273"
  )

## lookups to convert coded variables to full, descriptive variables ----

recoder <-
  lst(
    subgroups = c(
      `Main` = "all",
      `Third dose brand` = "vax3_type",
      `Prior SARS-CoV-2 infection` = "prior_covid_infection",
      `Primary course vaccine brand` = "vax12_type",
      `Age` = "ageband"
    ),
    status = c(
      `Unmatched`= "unmatched",
      `Matched` = "matched"
    ),
    treated = c(
      `Two doses` = "0",
      `Three doses` = "1"
    ),
    outcome = set_names(events_lookup$event, events_lookup$event_descr),
    all = c(` ` = "all"),
    prior_covid_infection = c(
      `No prior SARS-CoV-2 infection` = "FALSE",
      `Prior SARS-CoV-2 infection` = "TRUE"
    ),
    vax12_type = c(
      `BNT162b2` = "pfizer-pfizer",
      `ChAdOx1-S` = "az-az"
    ),
    ageband = c(
      `18-39 years` = "18-39",
      `40-64 years` = "40-64",
      `65-79 years` = "65-79",
      `80+ years` = "85+"
    )
  )


## follow-up time ----

# period width
postbaselinedays <- 28

# where to split follow-up time after recruitment
postbaselinecuts <- c(0, 14, 14 + (1:6)*postbaselinedays)

# maximum follow-up
maxfup <- max(postbaselinecuts)

# matching variables ----

# exact variables
exact_variables <- c(
  
  "jcvi_ageband",
  "cev_cv",
  "vax12_type",
  #"vax2_week",
  "region",
  #"sex",
  #"cev_cv",
  
  #"multimorb",
  "prior_covid_infection",
  #"immunosuppressed",
  #"status_hospplanned"
  NULL
)

# caliper variables
caliper_variables <- c(
  age = 3,
  vax2_day = 7,
  NULL
)
matching_variables <- c(exact_variables, names(caliper_variables))

# covariates ----

covariates <- c(
  NULL
)
