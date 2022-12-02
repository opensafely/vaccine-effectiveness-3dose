
library(tidyverse)
library(here)
library(glue)
library(cowplot)
# library(ggh4x)

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

plot_rates <- function(.data, filename, legend.position = "bottom") {
  
  key <- .data %>%
    distinct(anytest_cut) %>%
    mutate(
      lower = fup_params$covidtestcuts[-length(fup_params$covidtestcuts)],
      upper = fup_params$covidtestcuts[-1]
    ) %>%
    pivot_longer(
      cols = -anytest_cut
    ) %>%
    select(anytest_cut, time=value)
  
  plot_data <- .data %>%
    pivot_longer(
      cols = ends_with("rate")
    ) %>%
    mutate(longname=name) %>%
    mutate(across(longname, 
                  factor, 
                  levels = str_c(rates, "_rate"), 
                  labels = str_wrap(names(rates), 100)
                  )) %>%
    mutate(across(treated, 
                  factor, 
                  levels = c(0,1), 
                  labels = c("two doses", "three doses")
                  )) %>%
    left_join(key, by = "anytest_cut")
  
  plot_function <- function(.data) {
    .data %>%
      ggplot(aes(x = time, y = value, group = treated, colour = treated)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
      geom_step() +
      facet_wrap("longname", ncol=1) +
      labs(
        x = "time (days relative to trial date)",
        y = "rate per person-day of follow-up"
      ) +
      scale_x_continuous(
        breaks = fup_params$covidtestcuts
        ) +
      scale_y_continuous(
        limits = c(0, NA)
      ) +
      scale_color_discrete(
        name = NULL
      ) +
      theme_bw() +
      theme(
        axis.title.x = element_text(margin = margin(t=10)),
        axis.title.y = element_text(margin = margin(r=10)),
        legend.box.background = element_rect(colour = "black"),
        legend.position = legend.position
      )
  }
  
  p1 <- plot_data %>% 
    filter(name %in% glue("{rates[5:7]}_rate")) %>%
    droplevels() %>%
    plot_function() 
  
  ggsave(
    plot = p1,
    filename = file.path(output_dir, glue("rates_{filename}.png")),
    width = 15, height = 20, units = "cm"
  )
  
  return(p1)
  
}

data_rates_unrounded %>% plot_rates(filename = "unrounded", legend.position = c(0.8, 0.225))
data_rates_rounded %>% plot_rates(filename = "rounded", legend.position = c(0.8, 0.225))
