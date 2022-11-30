
library(tidyverse)
library(here)
library(glue)
library(cowplot)
library(ggh4x)

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

rates <- c( 
  "First positive SARS-CoV-2 test" = "firstpostest",
  "PCR only" = "pcronly", 
  "LFT only" = "lftonly", 
  "PCR and LFT" = "both",
  "Any SARS-CoV-2 test" = "anytest", 
  "SARS-CoV-2 test for symptomatic case" = "symptomatic",
  "Positive SARS-CoV-2 test" = "postest"
)

if(Sys.getenv("OPENSAFELY_BACKEND") %in% c("")) {
  
  ## Import released data ----
  release_dir <- ""
  
  output_dir <- here("output", release_dir, "figures")
  fs::dir_create(output_dir)
  
  data_rates_rounded <- read_csv(fs::path(output_dir, "covidtest_rates_rounded.csv"))
  
} else {
  
  output_dir <- ghere("output", cohort, "covidtests", "summary", subgroup)
  fs::dir_create(output_dir)
  
  ## import data ---
  data_anytest_sum <- read_rds(here("output", cohort, "covidtests", "process", "data_anytest_sum.rds")) %>%
    mutate(across(treated, as.factor)) 
  
  subgroup_sym <- sym(subgroup)
  
  # calculate rates ----
  
  calc_rates <- function(.data) {
    .data %>%
      mutate(across(matches(glue("{unname(rates)}_n")), ~.x/persondays_n)) %>%
      rename_with(.f = ~str_replace(.x, "_n", "_rate"), .cols = matches(glue("{unname(rates)}_n")))
  }
  
  data_n <- data_anytest_sum %>%
    mutate(all="all") %>%
    group_by(treated, anytest_cut, !!subgroup_sym) %>%
    summarise(
      total_n = roundmid_any(n(), threshold),
      persondays_n = sum(persondays),
      anytest_n = sum(sum_anytest),
      symptomatic_n = sum(sum_symptomatic),
      postest_n = sum(sum_postest),
      firstpostest_n = sum(sum_firstpostest),
      lftonly_n = sum(sum_lftonly),
      pcronly_n = sum(sum_pcronly),
      both_n = sum(sum_both),
      .groups = "keep"
    ) %>%
    ungroup() 
  
  data_rates_unrounded <- data_n %>%
    calc_rates()
  write_csv(data_rates_unrounded, fs::path(output_dir, "covidtest_rates_unrounded.csv"))
  
  data_rates_rounded <- data_n %>%
    mutate(across(ends_with("_n"), ~roundmid_any(.x, to = threshold))) %>%
    calc_rates()
  write_csv(data_rates_rounded, fs::path(output_dir, "covidtest_rates_rounded.csv"))
  
}

# plot covidtest_rates ----

plot_rates <- function(.data, filename) {
  
  plot_data <- data_rates_rounded %>%
    pivot_longer(
      cols = ends_with("rate")
    ) %>%
    mutate(longname=name) %>%
    mutate(across(longname, factor, levels = str_c(rates, "_rate"), labels = str_wrap(names(rates[1:7]), 20)))
  
  plot_function <- function(.data, legend.position = "none") {
    .data %>%
      ggplot(aes(x = anytest_cut, y = value, group = treated, colour = treated)) +
      geom_point() +
      geom_line() +
      facet_grid(~longname,  scales = "free") +
      labs(
        x = "time period (days relative to trial_date)",
        y = "rate per person-day of follow-up"
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 90),
        legend.position = legend.position
      )
  }
  
  p1 <- plot_data %>% 
    filter(name %in% glue("{rates[5:7]}_rate")) %>%
    droplevels() %>%
    plot_function() 
  
  p2 <- plot_data %>% 
    filter(name %in% glue("{rates[1:4]}_rate")) %>%
    droplevels() %>%
    plot_function(legend.position = "right")
  
  plot_legend <- get_legend(p2)
  p2 <- p2 + theme(legend.position = "none")  
  
  plot_grid(
    plot_grid(p1, plot_legend, nrow=1, rel_widths = c(3,1)),
    p2,
    nrow=2, rel_widths = c(10,1)
    # align = "v", axis = "l", nrow=2
    )
  
  plot_grid(p1, plot_legend, p2, nrow=2, rel_widths = c(3,1,4))
  
  ggsave(
    filename = file.path(output_dir, glue("rates_{filename}.png")),
    width = 15, height = 20, units = "cm"
  )
  
}


