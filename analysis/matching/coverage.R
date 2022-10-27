# # # # # # # # # # # # # # # # # # # # #
# Purpose: describe matching results
# reports on matching coverage
# # # # # # # # # # # # # # # # # # # # #


# Preliminaries ----


## Import libraries ----
library('tidyverse')
library('lubridate')
library('here')
library('glue')
library('arrow')

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

data_matched <- read_rds(ghere("output", cohort, "match", "data_matched.rds")) 

data_treatedeligible_matchstatus <- read_rds(here("output", cohort, "match", "data_treatedeligible_matchstatus.rds"))


# matching coverage on each day of recruitment period ----


# matching coverage for boosted people
data_coverage <-
  data_treatedeligible_matchstatus %>%
  group_by(vax3_type, vax3_date) %>%
  summarise(
    n_eligible = n(),
    n_matched = sum(matched, na.rm=TRUE),
    .groups = "keep"
  ) %>%
  ungroup() %>%
  mutate(
    n_unmatched = n_eligible - n_matched,
  ) %>%
  pivot_longer(
    cols = c(n_unmatched, n_matched),
    names_to = "status",
    names_prefix = "n_",
    values_to = "n"
  ) %>%
  arrange(vax3_type, vax3_date, status) %>%
  group_by(vax3_type, vax3_date, status) %>%
  summarise(
    n = sum(n),
    .groups = "keep"
  ) %>%
  group_by(vax3_type, status) %>%
  complete(
    vax3_date = full_seq(c(dates$start_date, dates$end_date), 1), # go X days before to
    fill = list(n=0)
  ) %>%
  mutate(
    cumuln = cumsum(n)
  ) %>%
  ungroup() %>%
  mutate(
    status = factor(status, levels=c("unmatched", "matched")),
    status_descr = fct_recoderelevel(status, recoder$status)
  ) %>%
  arrange(status_descr, vax3_date)

# plot coverage for checking (not release)
plot_coverage <- data_coverage %>%
  pivot_longer(cols = c(n, cumuln)) %>%
  mutate(across(name, factor, levels = c("n", "cumuln"))) %>%
  ggplot(aes(x=vax3_date, y=value, fill=status_descr)) +
  geom_bar(stat="identity") +
  facet_grid(name~vax3_type, scales = "free_y") +
  theme(legend.position = "bottom")
ggsave(
  filename=here("output", cohort, "table1", "coverage.png"),
  plot_coverage,
  width=15, height=20, units="cm"
)  

## round to nearest 6 for disclosure control
threshold <- 6

data_coverage_rounded <-
  data_coverage %>%
  group_by(status) %>%
  mutate(
    cumuln = roundmid_any(cumuln, to = threshold),
    n = diff(c(0,cumuln)),
  )

write_csv(data_coverage_rounded, fs::path(output_dir, "coverage.csv"))
