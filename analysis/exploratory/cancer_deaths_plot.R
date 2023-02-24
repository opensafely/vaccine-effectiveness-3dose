library(tidyverse)
library(here)
library(glue)

source(here("analysis", "design.R"))
source(here("lib", "functions", "utility.R"))
source(here("manuscript", "functions.R"))

if(Sys.getenv("OPENSAFELY_BACKEND") == ""){
  
  outdir <- here("manuscript")
  
  km_contrasts_rounded <- bind_rows(
    read_csv(here("release20230207", "km_contrasts_rounded.csv")),
    read_csv(here("release20230217", "km_contrasts_noncancer_rounded.csv")) %>%
      mutate(across(subgroup_level, as.character))
  ) 
  
  cox_contrasts_rounded <- bind_rows(
    read_csv(here("release20230207", "cox_contrasts_rounded.csv")),
    read_csv(here("release20230217", "cox_contrasts_noncancer_rounded.csv")) %>%
      mutate(across(subgroup_level, as.character))
  ) 
  
} else {
  
  km_contrasts_rounded <- read_csv(
    here("output", "mrna", "models", "combined", "km_contrasts_rounded.csv")
  ) 
  
  cox_contrasts_rounded <- read_csv(
    here("output", "mrna", "models", "combined", "cox_contrasts_rounded.csv")
  )
  
}

km_contrasts_rounded <- km_contrasts_rounded %>%
  filter(subgroup %in% c("all", "noncancer")) %>%
  mutate(
    subgroup_level = factor(subgroup, levels = c("all", "noncancer")),
    subgroup = factor("cancer_deaths")
  )

cox_contrasts_rounded <- cox_contrasts_rounded %>%
  filter(subgroup %in% c("all", "noncancer")) %>%
  mutate(
    subgroup_level = factor(subgroup, levels = c("all", "noncancer")),
    subgroup = factor("cancer_deaths")
  )

# redefine combine_plot
combine_plot <- function(
  subgroup_select = "all", 
  variant_option_select = "ignore"
) {
  
  if (subgroup_select == "all") {
    colour_var <- "variant_descr"
  } else {
    colour_var <- "subgroup_level"
  }
  
  print(cat("data A:\n"))
  # km plot
  A <- km_contrasts_rounded %>%
    filter(
      filename == "daily",
      subgroup %in% subgroup_select,
      variant_option %in% variant_option_select
    ) %>%
    select(
      outcome, subgroup, subgroup_level,
      variant_option, variant, 
      time = period_end,
      starts_with("risk")
    ) %>%
    pivot_longer(
      cols = starts_with("risk"),
      names_pattern = "(.*)_(\\d)",
      names_to = c(".value", "treated")
    ) %>%
    mutate(across(treated, as.integer)) %>%
    add_descr(subgroup = FALSE) %>%
    droplevels() 
  
  print(cat("plot A:\n"))
  A <- A %>%
    km_plot(colour_var)
  
  print(cat("cox data:\n"))
  # cox plot
  cox_data <- cox_contrasts_rounded %>%
    filter(filename == "cuts", model == "cox_adj")
  
  model <- "adjusted"
  
  if (variant_option_select == "split") {
    
    cox_data <- bind_rows(
      cox_data,
      # temporarily use unadjusted estimates for fracture, as the adjusted models
      # failed due to memory issues
      cox_contrasts_rounded %>%
        filter(filename == "cuts", model == "cox_unadj", outcome == "fracture") 
    )
    
  }
  
  # derive variant info
  cox_data <- cox_data %>%
    mutate(variant = str_extract(term, "delta|transition|omicron")) %>%
    mutate(across(term, ~str_trim(str_remove(.x, ";\\sdelta|;\\stransition|;\\somicron"), side="right"))) %>%
    mutate(across(variant, ~if_else(is.na(.x), "ignore", .x)))
  
  print(cat("data B:\n"))
  B <- cox_data %>%
    filter(
      subgroup %in% subgroup_select,
      variant_option %in% variant_option_select,
      str_detect(term, "^treated")
    ) %>%
    add_descr(subgroup=FALSE) %>%
    mutate(
      period_start = as.integer(str_extract(term, "\\d+")),
      period_end = as.integer(str_extract(term, "\\d+$")),
      midpoint = (period_start + period_end)/2,
      model = factor(model)
    ) %>%
    droplevels() 
  
  print(cat("plot B:\n"))
  B <- B %>%
    hr_plot(colour_var) 
  
  if ((variant_option_select != "ignore") | (subgroup_select != "all")) {
    
    plot_legend <- cowplot::get_legend(B)
    
    A <- A + theme(legend.position = "none")
    B <- B + theme(legend.position = "none")
    
  }
  
  print(cat("combine plots:\n"))
  p <- cowplot::plot_grid(
    # the NULL is reduce the space between the plots
    A, NULL, B,
    rel_widths = c(1, -0.15, 1),
    nrow = 1,
    labels = c("A", "", "B"), align = "v"
  )
  
  if ((variant_option_select != "ignore") | (subgroup_select != "all")) {
    
    p <-  cowplot::plot_grid(
      p, plot_legend,
      nrow=2, rel_heights = c(14,1)
    )
    
  }
  
  ggsave(
    filename = "combined_plot.png",
    path = outdir,
    plot = p,
    width = 17, height = 22, units = "cm"
  )
  
}

combine_plot(subgroup_select = "cancer_deaths", variant_option_select = "ignore")
