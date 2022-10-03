# create final dummy data for control population ----

library('tidyverse')
library('arrow')
library('here')
library('glue')

# import command-line arguments ----

args <- commandArgs(trailingOnly=TRUE)


if(length(args)==0){
  cohort <- "pfizer"
} else {
  cohort <- args[[1]]
}

fs::dir_create(here("output", cohort,  "dummydata"))


if(Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")){
  
  source(here("lib", "functions", "utility.R"))
  
  source(here("analysis", "design.R"))
  
  dummy_post <- read_feather(here("lib", "dummydata", "dummy_post.feather")) 
  
  vars_post <- names(dummy_post)
  vars_post <- vars_post[str_detect(vars_post, "_day$")]
  
  # dummy_control_final
  
  # import all datasets of matched controls, including matching variables
  data_matchedcontrols <- 
    map_dfr(
      seq_len(n_matching_rounds), 
      ~{read_rds(ghere("output", cohort, glue("matchround", .x), "actual", glue("data_successful_matchedcontrols.rds")))},
      .id="matching_round"
    ) %>%
    select(
      # select variables with_value_from_file
      patient_id, trial_date, match_id,
    )
  
  data_matchedcontrols %>%
    left_join(dummy_post , by = "patient_id") %>%
    mutate(across(all_of(vars_post), ~ trial_date + lubridate::days(.))) %>% # index on trial_date
    rename_with(~str_replace(., "_day", "_date"), all_of(vars_post)) %>%
    write_feather(sink = here("output", cohort, "dummydata", "dummy_control_final.feather"))
  
  
} else {
  
  # save empty output to save space if running on real data
  tibble() %>%
    write_feather(sink = here("output", cohort, "dummydata", "dummy_control_final.feather"))
  
}
