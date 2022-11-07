
library(tidyverse)
library(here)
library(glue)

## import local functions and parameters ---

source(here("analysis", "design.R"))

source(here("lib", "functions", "utility.R"))


## import command-line arguments ----

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  cohort <- "mrna"
} else {
  cohort <- args[[1]]
} 


## import data ---
data_covidtests <- read_rds(ghere("output", cohort, "covidtests", "extract", "data_extract.rds"))

data_matched <- read_rds(here("output", cohort, "match", "data_matched.rds")) %>%
  select(patient_id, trial_date, death_date, dereg_date, controlistreated_date)

## import baseline data, restrict to matched individuals and derive time-to-event variables
data_covidtests <- data_matched %>%
  left_join(data_covidtests, by = c("patient_id", "trial_date")) %>%
  mutate(
    censor_date = pmin(
      dereg_date,
      death_date,
      study_dates$testend_date,
      trial_date + maxfup,
      controlistreated_date,
      na.rm = TRUE
    ),
    censor_date = trial_date + maxfup # use this to overwrite above definition until issue with `patients.minimum_of()` and date arithmetic is fixed
  )

# report number of tests ----


data_counts <- data_matched %>%
  group_by(treated, !!subgroup_sym) %>%
  summarise(
    n = roundmid_any(n(), threshold),
    persontime = sum(as.numeric(censor_date - (trial_date - 1))),
    test_rate = sum(test_count) / persontime,
    postest_rate = sum(postest_count) / persontime,
  )

write_rds(data_counts, fs::path(output_dir, "testcounts.rds"))
