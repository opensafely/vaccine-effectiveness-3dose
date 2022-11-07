######################################

# This script:

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
  cohort <- "mrna"
} else {
  cohort <- args[[1]]
} 

## create output directory ----
outdir <- ghere("output", cohort, "covidtests", "extract")
fs::dir_create(outdir)

# import data ----

studydef_path <- file.path(outdir, "input_covidtests.feather")
data_studydef_dummy <- read_feather(studydef_path) %>%
  # because date types are not returned consistently by cohort extractor
  mutate(across(ends_with("_date"), ~ as.Date(.))) %>%
  # because of a bug in cohort extractor -- remove once fixed
  mutate(patient_id = as.integer(patient_id))

# use externally created dummy data if not running in the server
# check variables are as they should be
if(Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")){
  
  # generate custom dummy data
  data_custom_dummy <- local({
    
    # set seed so that dummy data results are reproducible
    set.seed(10)
    
    case_category <- c("LFT_Only"=0.4, "PCR_Only"=0.5, "LFT_WithPCR"=0.1)
    symptomatic <- c("N"=0.4, "Y"=0.6)
    
    data_tests <- read_csv(here("output", cohort, "match", "data_matched.csv.gz")) 
    
    # add number of tests variables
    data_tests <- data_tests %>%
      bind_cols(
        map_dfc(
          .x = 1:3,
          .f = ~tibble(!! sym(str_c("anytestpre_", .x, "_n")) := rpois(n=nrow(data_tests), lambda = 2))
        ),
        map_dfc(
          .x = 0:6,
          .f = ~tibble(!! sym(str_c("anytestpost_", .x, "_n")) := rpois(n=nrow(data_tests), lambda = 2))
        )
      )
    
    data_tests <- data_tests %>%
      mutate(
        anytest_missing = rbernoulli(n = nrow(data_tests), p=0.3),
        anytest_1_day = as.integer(runif(n=nrow(data_tests), -3*28, 14+6*28)),
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
        select(-anytest_missing)
      
    }
    
    data_tests <- data_tests %>%
      mutate(across(matches("anytest_\\d+_symptomatic"), factor, levels = c(0,1), labels = names(symptomatic))) %>%
      mutate(across(ends_with("_day"), ~ as.Date(as.character(trial_date + .)))) %>%
      rename_with(~str_replace(., "_day", "_date"), ends_with("_day")) 
    
    return(data_tests)
    
    
  })
  
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
  
  
  # if(nrow(unmatched_types)>0) stop(
  #   #unmatched_types
  #   "inconsistent typing in studydef : dummy dataset\n",
  #   apply(unmatched_types, 1, function(row) paste(paste(row, collapse=" : "), "\n"))
  # )
  
  data_extract <- data_custom_dummy 
  
} else {
  
  data_extract <- data_studydef_dummy
  rm(data_studydef_dummy)
  
}

# summarise and save data_extract
my_skim(data_extract, path = file.path(outdir, "input_treated_skim.txt"))
write_rds(data_extract, file.path(outdir, "data_extract.rds"), compress = "gz")

# check that the sums of the anytest_*_date variables match the anytest*_*_n variables
# rerun study definition with larger n if not
check_n <- data_extract %>%
  select(patient_id, trial_date, matches("anytest_\\d+_date")) %>%
  # mutate(across(matches("anytest_\\d+_date"))) %>%
  pivot_longer(
    cols = matches("anytest_\\d+_date"),
    values_drop_na = TRUE
    ) %>%
  select(-name) %>%
  arrange(patient_id, trial_date, value) %>%
  mutate(anytest_cut = cut(
    as.integer(value-trial_date),
    breaks = c(seq(-3*28, 0, 28), seq(14, 14+6*28, 28)),
    right=FALSE
    )
  ) %>%
  mutate(
    name = factor(
      as.integer(anytest_cut), 
      labels = c(str_c("pre_", 1:3), str_c("post_", 0:6)))
    ) %>%
  filter(!is.na(anytest_cut)) %>%
  group_by(patient_id, anytest_cut, name) %>%
  count() %>%
  ungroup() %>%
  left_join(
    data_extract %>%
      select(patient_id, matches("anytest\\w+_n")) %>%
      pivot_longer(
        cols = -patient_id
      ) %>%
      mutate(across(name, ~str_extract(.x, "p\\w+_\\d+"))),
    by = c("patient_id", "name")
  ) %>%
  mutate(percent = 100*n/value) 

# print:
cat ("summarise number of tests missing per period:\n")
check_n %>%
  mutate(n_missing = value - n) %>%
  group_by(anytest_cut, name) %>%
  summarise(across(n_missing, list(min=min, max=max, mean=mean, median=median))) %>%
  ungroup() %>%
  print(n=Inf)

# plot distribution of sum anytest*_date as a percent of anytest_*_n per period
p <- check_n %>%
  ggplot(aes(x=percent)) +
  geom_freqpoly(binwidth=1) +
  facet_wrap(~anytest_cut, scales = "free_y", nrow=2) +
  theme_bw()
ggsave(filename = file.path(outdir, "check_anytest.png"),
       plot = p, width = 20, height = 15, units = "cm")
  
