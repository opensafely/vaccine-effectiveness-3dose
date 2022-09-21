
library('tidyr')#, lib.loc = "C:/Program Files/R/R-4.1.3/library")
library('tidyverse')#, lib.loc = "C:/Program Files/R/R-4.1.3/library")
#library('tidyverse')
library('arrow')
library('here')
library('glue')

#source(here("analysis", "lib", "utility_functions.R"))

# remotes::install_github("https://github.com/wjchulme/dd4d")
library('dd4d')

source(here("lib", "functions", "utility.R"))

population_size <- 20000

# get nth largest value from list
nthmax <- function(x, n=1){
  dplyr::nth(sort(x, decreasing=TRUE), n)
}

source(here("analysis", "design.R"))

cohort <- "over12"
n_matching_rounds <- 2

dates <- map(study_dates[[cohort]], as.Date)
params <- study_params[[cohort]]


start_date <- as.Date(dates$start_date)
end_date <- as.Date(dates$end_date)
followupend_date <- as.Date(dates$followupend_date)
index_date <- as.Date(dates$start_date)

first_pfizerA_date <- as.Date(dates$start_date)
first_pfizerC_date <- as.Date(dates$start_date)

index_day <- 0L
start_day <- as.integer(start_date - index_date)
end_day <- as.integer(end_date - index_date)
first_pfizerA_day <- as.integer(first_pfizerA_date - index_date)
first_pfizerC_day <- as.integer(first_pfizerC_date - index_date)

known_variables <- c(
  "index_date", "start_date", "end_date", "first_pfizerA_date", "first_pfizerC_date",
  "index_day",  "start_day", "end_day", "first_pfizerA_day", "first_pfizerC_day"
)


data_matchstatus <- read_rds(ghere("output", cohort, "matchround{n_matching_rounds}", "actual", "data_matchstatus_allrounds.rds")) %>% filter(treated==0L)


# import all datasets of matched controls, including matching variables
data_matchedcontrols <- 
  map_dfr(
    seq_len(n_matching_rounds), 
    ~{read_rds(ghere("output", cohort, glue("matchround", .x), "actual", glue("data_successful_matchedcontrols.rds")))},
    .id="matching_round"
  ) %>%
  mutate(
    trial_day = as.integer(trial_date - start_date)
  ) %>%
  select(
    # see study_definition_finalmatched.py for variables to include
    
    # select variables with_value_from_file
    patient_id, trial_day, match_id,
    
    ## select variables in `variables_matching.py`
    ## or not, if they are saved in the "data_Successful_match" output in `matching_filter.R`
    # sex,
    # ethnicity,
    # practice_id,
    # msoa,
    # stp,
    # region,
    # imd_Q5,
    # primary_care_covid_case_0_date,
    # postest_0_date,
    # covidadmitted_0_date,
    # covidemergency_0_date,
    
    # variables in `variables_outcomes.py` are simulated below
    
  )
  

sim_list = lst(
  
  dereg_day = bn_node(
    ~as.integer(runif(n=..n, trial_day, trial_day+120)),
    missing_rate = ~0.99
  ),

  primary_care_covid_case_day = bn_node(
    ~as.integer(runif(n=..n, trial_day, trial_day+100)),
    missing_rate = ~0.7
  ),

  covid_test_day = bn_node(
    ~as.integer(runif(n=..n, trial_day, trial_day+100)),
    missing_rate = ~0.7
  ),

  postest_day = bn_node(
    ~as.integer(runif(n=..n, trial_day, trial_day+100)),
    missing_rate = ~0.7
  ),

  emergency_day = bn_node(
    ~as.integer(runif(n=..n, trial_day, trial_day+200)),
    missing_rate = ~0.8
  ),
  emergencyhosp_day = bn_node(
    ~as.integer(runif(n=..n, trial_day, trial_day+200)),
    missing_rate = ~0.85
  ),


  covidemergency_day = bn_node(
    ~as.integer(runif(n=..n, trial_day, trial_day+200)),
    missing_rate = ~0.8
  ),

  covidemergencyhosp_day = bn_node(
    ~as.integer(runif(n=..n, trial_day, trial_day+200)),
    missing_rate = ~0.85
  ),

  # respemergency_day = bn_node(
  #   ~as.integer(runif(n=..n, trial_day, trial_day+100)),
  #   missing_rate = ~0.95
  # ),
  # 
  # respemergencyhosp_day = bn_node(
  #   ~as.integer(runif(n=..n, trial_day, trial_day+100)),
  #   missing_rate = ~0.95
  # ),

  covidadmitted_day = bn_node(
    ~as.integer(runif(n=..n, trial_day, trial_day+100)),
    missing_rate = ~0.7
  ),

  # placeholder for single criticalcare variable ---
  covidcritcare_day = bn_node(
    ~as.integer(runif(n=..n, trial_day, trial_day+100)),
    missing_rate = ~0.8
  ),
  
  admitted_unplanned_day = bn_node(
    ~as.integer(runif(n=..n, trial_day, trial_day+100)),
    missing_rate = ~0.7
  ),
  # admitted_planned_day = bn_node(
  #   ~as.integer(runif(n=..n, trial_day, trial_day+100)),
  #   missing_rate = ~0.7
  # ),

  coviddeath_day = bn_node(
    ~death_day,
    missing_rate = ~0.7,
    needs = "death_day"
  ),

  death_day = bn_node(
    ~as.integer(runif(n=..n, trial_day, trial_day+100)),
    missing_rate = ~0.90
  ),

)
bn <- bn_create(sim_list, known_variables = c(known_variables, names(data_matchedcontrols)))

bn_plot(bn)
bn_plot(bn, connected_only=TRUE)

set.seed(10)

dummydata <-bn_simulate(bn, known_df = data_matchedcontrols, keep_all = FALSE, .id="patient_id")


dummydata_processed <- dummydata %>%
  #convert logical to integer as study defs output 0/1 not TRUE/FALSE
  # mutate(across(where(is.logical), ~ as.integer(.))) %>%
  #convert integer days to dates since index date and rename vars
  mutate(across(ends_with("_day"), ~ as.Date(as.character(index_date + .)))) %>%
  rename_with(~str_replace(., "_day", "_date"), ends_with("_day"))


fs::dir_create(here("lib", "dummydata"))
write_feather(dummydata_processed, sink = here("lib", "dummydata", "dummy_controlfinal.feather"))
