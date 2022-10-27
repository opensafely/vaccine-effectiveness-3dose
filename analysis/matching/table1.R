# # # # # # # # # # # # # # # # # # # # #
# Purpose: describe matching results
# creates "table 1"
# # # # # # # # # # # # # # # # # # # # #


# Preliminaries ----


## Import libraries ----
library('tidyverse')
library('lubridate')
library('here')
library('glue')
library('arrow')
library('gt')
library('gtsummary')

## import local functions and parameters ---

source(here("analysis", "design.R"))
source(here("lib", "functions", "utility.R"))
source(here("lib", "functions", "redaction.R"))

# import command-line arguments ----

args <- commandArgs(trailingOnly=TRUE)


if(length(args)==0){
  # use for interactive testing
  removeobjects <- FALSE
  cohort <- "mrna"
} else {
  #FIXME replace with actual eventual action variables
  removeobjects <- TRUE
  cohort <- args[[1]]
}


## get cohort-specific parameters study dates and parameters ----
dates <- study_dates[[cohort]]


## create output directories ----

output_dir <- here("output", cohort, "table1")
fs::dir_create(output_dir)

## Import data and derive some variables ----

data_matched <- read_rds(ghere("output", cohort, "match", "data_matched.rds")) %>%
  mutate(dayssincevax2 = as.integer(trial_date - vax2_date))

data_treatedeligible_matchstatus <- read_rds(here("output", cohort, "match", "data_treatedeligible_matchstatus.rds"))

# table 1 style baseline characteristics ----

var_labels <- list(
  N  ~ "Total N",
  
  treated ~ "Status",
  age ~ "Age",
  jcvi_ageband ~ "JCVI ageband",
  sex ~ "Sex",
  ethnicity ~ "Ethnicity",
  imd_Q5 ~ "Deprivation",
  region ~ "Region",
  rural_urban_group ~ "Rurality",
  
  bmi ~ "Body mass index",
  
  learndis ~ "Learning disability",
  sev_mental ~ "Severe mental illness",
  immunosuppressed ~ "Immunouppressed",
  
  multimorb ~ "Multimorbidity score",
  
  pregnancy ~ "Pregnancy",
  
  cv ~ "Clinically vulnerable",
  cev ~ "Clinically extremely vulnerable",
  
  flu_vaccine ~ "Flu vaccine",
  
  vax12_type ~ "Primary course vaccine type",
  dayssincevax2 ~ "Days since second dose",
  
  prior_test_cat ~ "Number of SARS-CoV-2 tests during unvaccinated period",
  prior_covid_infection ~ "Prior documented SARS-CoV-2 infection",
  time_since_infection ~ "Time since last evidence of SARS-CoV-2 infection"
  
  
  
) %>%
set_names(., map_chr(., all.vars))

map_chr(var_labels[-c(1,2)], ~last(as.character(.)))


# use gtsummary to obtain stnadardised table 1 data
tab_summary_baseline <-
  data_matched %>%
  mutate(
    N = 1L,
    #treated_descr = fct_recoderelevel(as.character(treated), recoder$treated),
    age = factor(age, levels=sort(unique(age)))
  ) %>%
  select(
    treated,
    all_of(names(var_labels)),
  ) %>%
  tbl_summary(
    by = treated,
    label = unname(var_labels[names(.)]),
    statistic = list(N = "{N}")
  ) 

raw_stats <- tab_summary_baseline$meta_data %>%
  select(var_label, df_stats) %>%
  unnest(df_stats)


raw_stats_redacted <- raw_stats %>%
  mutate(
    n=roundmid_any(n, 6),
    N=roundmid_any(N, 6),
    p=n/N,
    N_miss = roundmid_any(N_miss, 6),
    N_obs = roundmid_any(N_obs, 6),
    p_miss = N_miss/N_obs,
    N_nonmiss = roundmid_any(N_nonmiss, 6),
    p_nonmiss = N_nonmiss/N_obs,
    var_label = factor(var_label, levels=map_chr(var_labels[-c(1,2)], ~last(as.character(.)))),
    variable_levels = replace_na(as.character(variable_levels), "")
  ) 

write_csv(raw_stats_redacted, fs::path(output_dir, "table1.csv"))

# cross tab cv and cev ----
table(data_matched$cev_cv, data_matched$cev)
