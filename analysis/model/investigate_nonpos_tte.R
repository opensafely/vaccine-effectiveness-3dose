library(tidyverse)
library(glue)
library(here)

## import local functions and parameters ---

source(here("analysis", "design.R"))
source(here("lib", "functions", "utility.R"))
source(here("lib", "functions", "survival.R"))

## import command-line arguments ----

args <- commandArgs(trailingOnly=TRUE)


if(length(args)==0){
  # use for interactive testing
  cohort <- "mrna"
  
} else {
  cohort <- args[[1]]
}

output_dir <- ghere("output", cohort, "models")
fs::dir_create(output_dir)

# read final matched data
data_matched <- read_rds(ghere("output", cohort, "match", "data_matched.rds"))

# identify the patients who died or deregistered before their trial date
data_diedordereg <- data_matched %>%
  filter(death_date < trial_date | dereg_date < trial_date) %>%
  select(patient_id, treated, death_date, dereg_date)

cat("Number of deaths and deregistrations before trial_date in:\n")
for (m in 1:n_matching_rounds_list[[cohort]]) {
  
  cat(glue("Matching round {m}"), "\n")
  # read in the patients who were successfully matched in this round
  data_matchround <- read_rds(
    ghere("output", cohort, "matchround{m}", "actual", "data_successful_matchedcontrols.rds")
    ) %>%
    select(patient_id) %>%
    # keep only those who had died or deregistered before their trial date in the final matched dataset
    inner_join(data_diedordereg, by = "patient_id") %>%
    group_by(treated) %>%
    # summarise number of deaths and deregistrations by treatment group
    summarise(
      death = sum(death_date < trial_date), 
      dereg = sum(dereg_date < trial_date),
      .groups="keep"
    ) %>% 
    print()
  
}

# empty output so that action succeeds
tibble() %>% write_csv(file.path(output_dir, "tmp.csv"))
