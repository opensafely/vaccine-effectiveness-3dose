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



# number of matching rounds to perform

n_matching_rounds <- 2


# define key dates ----

study_dates <- lst(
  over12 = lst(
   start_date = "2021-09-20", #start of recruitment monday 20 september pfizer licensed for for over 12yo in england
   end_date = "2021-12-19", # end of recruitment (13 weeks later)
   followupend_date = "2022-01-02", # end of follow-up
  ),
  under12 = lst(
    start_date = "2022-04-04", #start of recruitment monday 4 april moderna licensed for for under 12yo in england
    end_date = "2022-07-03", # end of recruitment (13 weeks later)
    followupend_date = "2022-07-10" # end of follow-up
  )
)

extract_increment <- 14

study_dates$over12$control_extract_dates = as.Date(study_dates$over12$start_date) + (0:26)*extract_increment
study_dates$under12$control_extract_dates = as.Date(study_dates$under12$start_date) + (0:26)*extract_increment

jsonlite::write_json(study_dates, path = here("lib", "design", "study-dates.json"), auto_unbox=TRUE, pretty =TRUE)

study_params <- lst(
  
  # over 12 params
  over12 = lst(
   minage = 12,
   maxage= 15,
  
   treatment = "pfizerA",
  ),
  
  # under 12 params
  under12 = lst(
    minage = 5,
    maxage= 11,
    treatment = "pfizerC",
  )
  
)

jsonlite::write_json(study_params, path = here("lib", "design", "study-params.json"), auto_unbox=TRUE, pretty =TRUE)

# define outcomes ----

events_lookup <- tribble(
  ~event, ~event_var, ~event_descr,

  # other
  "test", "covid_test_date", "SARS-CoV-2 test",
  "dereg", "dereg_date", "Deregistration",
  "primary_care_covid_case", "primary_care_covid_case_date", "Primary care COVID-19",
  
  # effectiveness
  "postest", "postest_date", "Positive SARS-CoV-2 test",
  "covidemergency", "covidemergency_date", "COVID-19 A&E attendance",
  "covidemergencyhosp", "covidemergencyhosp_date", "Admission from COVID-19 A&E",
  "covidadmitted", "covidadmitted_date", "COVID-19 hospitalisation",
  "noncovidadmitted", "noncovidadmitted_date", "Non-COVID-19 hospitalisation",
  
  "covidadmittedproxy1", "covidadmittedproxy1_date", "COVID-19 hospitalisation (A&E proxy)",
  "covidadmittedproxy2", "covidadmittedproxy2_date", "COVID-19 hospitalisation (A&E proxy v2)",
  "covidcritcare", "covidcritcare_date", "COVID-19 critical care",
  "coviddeath", "coviddeath_date", "COVID-19 death",
  "noncoviddeath", "noncoviddeath_date", "Non-COVID-19 death",
  "death", "death_date", "Any death",

  # safety
  "emergencyhosp", "emergencyhosp_date", "Admission from A&E",
  "admitted", "admitted_unplanned_date", "Unplanned hospitalisation",
  "emergency", "emergency_date", "A&E attendance",
)


## lookups to convert coded variables to full, descriptive variables ----

recoder <-
  lst(
    subgroups = c(
      `Main` = "all",
      `Prior SARS-CoV-2 infection` = "prior_covid_infection"
    ),
    status = c(
      `Unmatched`= "unmatched",
      `Matched` = "matched"
    ),
    treated = c(
      `Unvaccinated` = "0",
      `Vaccinated` = "1"
    ),
    outcome = set_names(events_lookup$event, events_lookup$event_descr),
    all = c(` ` = "all"),
    prior_covid_infection = c(
      `No prior SARS-CoV-2 infection` = "FALSE",
      `Prior SARS-CoV-2 infection` = "TRUE"
    ),
  )


## follow-up time ----

# period width
postbaselinedays <- 14

# where to split follow-up time after recruitment
postbaselinecuts <- (0:10)*postbaselinedays

# maximum follow-up
maxfup <- max(postbaselinecuts)

# matching variables ----

# exact variables
exact_variables <- c(
  "age_aug21",
  "region",
  "sex",
  "prior_covid_infection",
  NULL
)

# caliper variables
caliper_variables <- c(
  #age = 1,
  NULL
)
matching_variables <- c(exact_variables, names(caliper_variables))


