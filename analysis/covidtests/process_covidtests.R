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
source(here("lib", "functions", "survival.R"))

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
    
    data_tests <- read_rds(here("output", cohort, "match", "data_matched.rds")) %>%
      select(patient_id, trial_date)
    
    # add number of tests variables
    data_tests <- data_tests %>%
      bind_cols(
        map_dfc(
          .x = 1:prebaselineperiods,
          .f = ~tibble(
            !! sym(str_c("anytestpre_", .x, "_n")) := rpois(n=nrow(data_tests), lambda = 2),
            !! sym(str_c("postestpre_", .x, "_n")) := as.integer(max(
              !! sym(str_c("anytestpre_", .x, "_n")) - rpois(n=nrow(data_tests), lambda = 1),
              0
              ))
            )
        ),
        map_dfc(
          .x = 0:postbaselineperiods,
          .f = ~tibble(
            !! sym(str_c("anytestpost_", .x, "_n")) := rpois(n=nrow(data_tests), lambda = 2),
            !! sym(str_c("postestpost_", .x, "_n")) := as.integer(max(
              !! sym(str_c("anytestpost_", .x, "_n")) - rpois(n=nrow(data_tests), lambda = 1),
              0
              ))
            )
        )
      )
    
    data_tests <- data_tests %>%
      mutate(
        # first anytest during the period
        anytest_1_day = if_else(
          rbernoulli(n = nrow(data_tests), p=0.3),
          NA_integer_,
          as.integer(runif(
            n=nrow(data_tests),
            -prebaselineperiods*postbaselinedays, 
            baselinedays+postbaselineperiods*postbaselinedays
          ))
        ),
        # symptom category of anytest_1_day
        anytest_1_symptomatic = if_else(
          is.na(anytest_1_day),
          NA_integer_,
          as.integer(rbernoulli(n=nrow(data_tests), p=symptomatic[["Y"]]))
          ),
        # first positive test in the period (anytest_1_day with added missingness for negative tests)
        postest_1_day = if_else(
          rbernoulli(n = nrow(data_tests), p=0.5),
          NA_integer_,
          anytest_1_day
        ),
        # first positive test in the period and ever (postest_1_day with added missingness for those that have had a positive test before the period)
        firstpostest_day = if_else(
          is.na(postest_1_day) | rbernoulli(n = nrow(data_tests), p=0.3),
          NA_integer_, 
          postest_1_day
        ),
        # type of test on firstpostest_day
        firstpostest_category = factor(
          if_else(
            is.na(firstpostest_day),
            NA_character_, 
            sample(x=names(case_category), size = nrow(data_tests), prob = unname(case_category), replace=TRUE)
          ),
          levels = names(case_category)
        )
      ) 
    
    for (i in 2:n_any) {
      
      # derive subsequent anytest_*_day and anytest_*_symptomatic
      data_tests <- data_tests %>%
        mutate(
          !! sym(glue("anytest_{i}_day")) := if_else(
            rbernoulli(n = nrow(data_tests), p=0.3),
            NA_integer_,
            !! sym(glue("anytest_{i-1}_day")) + as.integer(runif(n=nrow(data_tests), 0, 50))
          ),
          !! sym(glue("anytest_{i}_symptomatic")) := if_else(
            is.na(!! sym(glue("anytest_{i}_day"))),
            NA_integer_,
            as.integer(rbernoulli(n=nrow(data_tests), p=symptomatic[["Y"]]))
          )
        )
      
      # derive subsequent postest_*_day
      if (i <= n_pos) {
        data_tests <- data_tests %>%
          mutate(
            !! sym(glue("postest_{i}_day")) := if_else(
              rbernoulli(n = nrow(data_tests), p=0.5),
              NA_integer_,
              !! sym(glue("anytest_{i}_day"))
            )
          ) 
      }
      
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
  
  
  if(nrow(unmatched_types)>0) stop(
    #unmatched_types
    "inconsistent typing in studydef : dummy dataset\n",
    apply(unmatched_types, 1, function(row) paste(paste(row, collapse=" : "), "\n"))
  )
  
  data_extract <- data_custom_dummy 
  
} else {
  
  data_extract <- data_studydef_dummy
  rm(data_studydef_dummy)
  
}

# summarise and save data_extract
my_skim(data_extract, path = file.path(outdir, "input_treated_skim.txt"))

data_split <- local({
  
  # derive censor date
  data_matched <- read_rds(here("output", cohort, "match", "data_matched.rds")) %>%
    select(patient_id, trial_date, death_date, dereg_date, controlistreated_date) %>%
    mutate(
      censor_date = pmin(
        dereg_date,
        death_date,
        study_dates$testend_date,
        trial_date - 1 + maxfup,
        controlistreated_date - 1,
        na.rm = TRUE
      ),
      tte_censor = as.integer(censor_date-(trial_date-1)),
      ind_outcome = 0
      # censor_date = trial_date + maxfup # use this to overwrite above definition until issue with `patients.minimum_of()` and date arithmetic is fixed
    ) %>%
    select(patient_id, trial_date, censor_date, tte_censor) %>%
    group_by(patient_id, trial_date) %>%
    mutate(new_id = cur_group_id()) %>% 
    ungroup()
  
  # derive fup_split (extra processing required when variant_option %in% c("split", "restrict"))
  fup_split <-
    data_matched %>%
    select(new_id) %>%
    uncount(weights = length(postbaselinecuts)-1, .id="period_id") %>%
    mutate(
      fupstart_time = postbaselinecuts[period_id]#,
      # fupend_time = postbaselinecuts[period_id+1]-1,
    ) %>%
    droplevels() %>%
    select(
      new_id, period_id, fupstart_time#, fupend_time
    ) 
  
  data_split <-
    tmerge(
      data1 = data_matched,
      data2 = data_matched,
      id = new_id,
      tstart = 0,
      tstop = tte_censor
    ) %>%
    # add post-treatment periods
    tmerge(
      data1 = .,
      data2 = fup_split,
      id = new_id,
      period_id = tdc(fupstart_time, period_id)
    ) %>%
    mutate(
      fup_cut = cut(
        tstart,
        breaks = c(seq(-3*28, 0, 28), seq(14, 14+6*28, 28)),
        right=FALSE
      )
    ) %>%
    transmute(
      patient_id, trial_date, censor_date, fup_cut,
      persondays = as.integer(tstop-tstart)
    )
  
  return(as_tibble(data_split))
  
})


# reshape data_extract
data_anytest_long <- data_extract %>%
  select(patient_id, trial_date, matches("\\w+test_\\d+_\\w+")) %>%
  # rename to make it easier to reshape
  rename_with(
    .fn = ~str_c(str_extract(.x, "\\d+_"), str_remove(.x, "\\d+_")), 
    .cols = matches("\\w+test_\\d+_\\w+")
    ) %>%
  pivot_longer(
    cols = matches("\\d+_\\w+test_\\w+"),
    names_to = c("index", ".value"),
    names_pattern = "(.*)_(.*_.*)",
    values_drop_na = TRUE
  ) %>%
  select(-index) %>% 
  left_join(
    data_split %>%
      distinct(patient_id, trial_date, censor_date), 
    by = c("patient_id", "trial_date")
    ) %>%
  # duplicate variables for transforming
  mutate(
    anytestcensor_date=anytest_date,
    postestcensor_date=postest_date,
    anytest_cut=anytest_date,
    postest_cut=postest_date,
  ) %>%
  # create censored versions of anytest_date and postest_date
  mutate(across(
    c(anytestcensor_date, postestcensor_date),
    ~if_else(
      .x <= censor_date,
      .x,
      as.Date(NA_character_)
    )
    )) %>%
  # create binned versions of anytest_date and postest_date
  mutate(across(
    c(anytest_cut, postest_cut),
    ~ cut(
      as.integer(.x-trial_date),
      breaks = c(seq(-3*28, 0, 28), seq(14, 14+6*28, 28)),
      right=FALSE
    )
  )) %>%
  arrange(patient_id, trial_date, anytest_date) %>%
  mutate(
    # create a variable with same labelling as the \\w+test\\w+_n variables for joining
    name = factor(
      as.integer(anytest_cut), 
      labels = c(str_c("pre_", 1:3), str_c("post_", 0:6)))
  ) %>%
  # remove any that are outside the time periods of interest
  filter(!is.na(anytest_cut)) %>%
  # sum the number of tests per period
  group_by(patient_id, trial_date, anytest_cut, name) %>%
  summarise(
    # sum all dates (this is just used to check value of n in study definition if correct)
    sum_anytest=n(), 
    sum_postest=n(), 
    # sum censored dates (this will be used to calculate testing rates)
    sum_anytestcensor=sum(!is.na(anytestcensor_date)),  
    sum_postestcensor=sum(!is.na(anytestcensor_date)),  
    .groups="keep"
    ) %>%
  ungroup() %>%
  # join the total number of tests per period with returning=xxx in study definition
  left_join(
    data_extract %>%
      select(patient_id, trial_date, matches("\\w+test\\w+_n")) %>%
      pivot_longer(
        cols = matches("\\w+test\\w+_n"),
        names_pattern = "(.*test)(.*)_n",
        names_to = c(".value", "name")
      ),
    by = c("patient_id", "trial_date", "name")
  ) %>%
  # join data_split for persondays of follow-up
  left_join(
    data_split %>% select(-censor_date),
    by = c("patient_id", "trial_date", "anytest_cut" = "fup_cut")
  ) %>%
  # fill in persondays for prebaseline periods 
  # (must be postbaselinedays, otherwise they would have been censored before trialdate)
  mutate(across(
    persondays,
    ~if_else(
      str_detect(name, "^pre"),
      as.integer(postbaselinedays),
      persondays
    )))

# save firstpostest variables 
data_extract %>%
  select(patient_id, trial_date, starts_with("firstpostest")) %>%
  filter(!is.na(firstpostest_date)) %>%
  write_rds(file.path(outdir, "data_firstpostest.rds"), compress = "gz")
# save long variables
data_anytest_long %>%
  write_rds(file.path(outdir, "data_anytest_long.rds"), compress = "gz")

# checks ----

# sense check
cat("-----------------")
cat("Sense checks ----\n")
cat("When persondays=NA, check sum_*testcensor=0:\n")
data_anytest_long %>%
  filter(is.na(persondays)) %>%
  group_by(anytest_cut) %>%
  summarise(across(
    ends_with("censor"),
    list(min=min, max=max)
  ))

# check that the sums of the anytest_*_date variables match the anytest*_*_n variables
# if not, it's a flag that we need to increase n in the study definition
cat("------------------------------------------")
cat("Check `n_any` and `n_pos` appropriate ----\n")

cat ("Summarise number of tests missing per person per period when summing dates:\n")
data_anytest_long %>%
  mutate(
    n_missing_anytest = anytest - sum_anytest,
    n_missing_postest = postest - sum_postest
    ) %>%
  group_by(anytest_cut, name) %>%
  summarise(across(
    starts_with("n_missing"), 
    list(min=min, max=max, mean=mean, median=median)
    ), .groups = "keep") %>%
  ungroup() %>%
  pivot_longer(
    cols = starts_with("n_missing"),
    names_pattern = "n_missing_(.*)_(.*)",
    names_to = c("result", ".value")
  ) %>%
  arrange(result, anytest_cut) %>%
  group_split(result) %>% as.list() 

cat(glue("see {file.path(outdir, \"check_*.png\")} for distribution of sum *test*_date as a percent of *test_*_n per period"), "\n")

plot_function <- function(result) {
  p <- data_anytest_long %>%
    mutate(percent = 100*!!sym(glue("sum_{result}"))/!!sym(result)) %>%
    ggplot(aes(x=percent)) +
    geom_freqpoly(binwidth=1) +
    facet_wrap(~anytest_cut, scales = "free_y", nrow=2) +
    theme_bw()
  ggsave(filename = file.path(outdir, glue("check_{result}.png")),
         plot = p, width = 20, height = 15, units = "cm")
  return(p)
}

plot_function("anytest")
plot_function("postest")


