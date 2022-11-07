# create dummy data for the tests data ----

library('tidyverse')
library('arrow')
library('here')
library('glue')

# import command-line arguments ----

args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  cohort <- "mrna"
} else {
  cohort <- args[[1]]
}

fs::dir_create(here("output", cohort,  "dummydata"))

if(Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")){ 
  
  source(here("lib", "functions", "utility.R"))
  
  source(here("analysis", "design.R"))
  
  index_date <- as.Date("2020-01-01") # doesn't matter what this is, just need some constant date
  
  set.seed(10)
  
  data_tests <- read_csv(here("output", cohort, "match", "data_matched.csv.gz")) 
  
  input_tests <- read_feather(here("output", cohort, "tests", "extract", "input_tests.feather"))
  
  case_category <- c("LFT_Only"=0.4, "PCR_Only"=0.5, "LFT_WithPCR"=0.1)
  symptomatic <- c("N"=0.4, "Y"=0.6)
  
  # add number of tests variables
  data_tests <- data_tests %>%
    bind_cols(
      map_dfc(
        .x = 0:6,
        .f = ~tibble(!! sym(str_c("anytest_", .x, "_n")) := rpois(n=nrow(data_tests), lambda = 2))
      )
    )
  
  
  data_tests <- data_tests %>%
    mutate(
      anytest_missing = rbernoulli(n = nrow(data_tests), p=0.3),
      anytest_1_day = as.integer(runif(n=nrow(data_tests), 0, 100)),
      anytest_1_symptomatic = as.integer(rbernoulli(n=nrow(data_tests), p=symptomatic[["Y"]]))
      ) %>%
    mutate(across(c(anytest_1_day, anytest_1_symptomatic), ~if_else(anytest_missing, NA_integer_, .x))) %>%
    mutate(
      firstpostestmissing = anytest_missing | rbernoulli(n = nrow(data_tests), p=0.3),
      firstpostest_day = if_else(firstpostestmissing, NA_integer_, anytest_1_day),
      firstpostest_category = factor(
        if_else(
        firstpostestmissing,
        NA_character_, 
        sample(x=names(case_category), size = nrow(data_tests), prob = unname(case_category), replace=TRUE)
        ),
        levels = names(case_category)
    )
    ) %>%
    select(-anytest_missing, -firstpostestmissing)
  
for (i in 2:10) {
    
    data_tests <- data_tests %>%
      mutate(
        anytest_missing = rbernoulli(n = nrow(data_tests), p=0.3),
        !! sym(glue("anytest_{i}_day")) := !! sym(glue("anytest_{i-1}_day")) + as.integer(runif(n=nrow(data_tests), 0, 50)),
        !! sym(glue("anytest_{i}_symptomatic")) := as.integer(rbernoulli(n=nrow(data_tests), p=symptomatic[["Y"]]))
      ) %>%
      mutate(across(c(glue("anytest_{i}_day"), glue("anytest_{i}_symptomatic")), ~if_else(anytest_missing, NA_integer_, .x))) %>%
      mutate(across(glue("anytest_{i}_symptomatic"), factor, levels = c(0,1), labels = names(symptomatic))) %>%
      select(-anytest_missing)
    
  }
    
  
  data_tests <- data_tests %>%
    mutate(across(ends_with("_day$"), ~trial_date+.x)) %>%
    rename_with(.f = ~str_replace(.x, "_day$", "_date$")) 
    
  
} else {
  
  # save empty output to save space if running on real data
  tibble() %>%
    write_feather(sink = here("output", cohort, "dummydata", "dummy_tests.feather"))
  
}
