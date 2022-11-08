
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
  subgroup <- "all"
} else {
  cohort <- args[[1]]
  subgroup <- args[[2]]
} 

output_dir <- ghere("output", cohort, "covidtests", "summary")
fs::dir_create(output_dir)

## import data ---
data_anytest_long <- read_rds(here("output", cohort, "covidtests", "process", "data_anytest_long.rds")) %>%
  mutate(across(treated, as.factor)) 

subgroup_sym <- sym(subgroup)

# calculate rates ----

data_counts <- data_anytest_long %>%
  mutate(all="all") %>%
  group_by(treated, anytest_cut, !!subgroup_sym) %>%
  summarise(
    n = roundmid_any(n(), threshold),
    total_persondays = sum(persondays),
    anytest_rate = sum(sum_anytest) / total_persondays,
    symptomatic_rate = sum(sum_symptomatic) / total_persondays,
    postest_rate = sum(sum_postest) / total_persondays,
    firstpostest_rate = sum(sum_firstpostest) / total_persondays,
    lftonly_rate = sum(sum_lftonly) / total_persondays,
    pcronly_rate = sum(sum_pcronly) / total_persondays,
    both_rate = sum(sum_both) / total_persondays,
    .groups = "keep"
  )

write_rds(data_counts, fs::path(output_dir, "rates.rds"))

rates <- c( 
  "Any SARS-CoV-2 test" = "anytest", 
  "SARS-CoV-2 test for symptomatic case" = "symptomatic",
  "Positive SARS-CoV-2 test" = "postest", 
  "First positive SARS-CoV-2 test" = "firstpostest",
  "PCR only" = "pcronly", 
  "LFT only" = "lftonly", 
  "PCR and LFT" = "both"
  )


# plot rates ----
data_counts %>%
  pivot_longer(
    cols = ends_with("rate")
  ) %>%
  mutate(across(name, factor, levels = str_c(rates, "_rate"), labels = str_wrap(names(rates), 20))) %>%
  ggplot(aes(x = anytest_cut, y = value, group = treated, colour = treated)) +
  geom_point() +
  geom_line() +
  facet_wrap(~name, nrow=2) +
  labs(
    x = "time period (days relative to trial_date)",
    y = "rate per person-day of follow-up"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90),
    legend.position = c(0.9, 0.15)
    )
ggsave(
  filename = file.path(output_dir, "rates.png"),
  width = 15, height = 20, units = "cm"
)