library(tidyverse)
library(readr)
library(here)
library(glue)
library(flextable)

source(here("analysis", "design.R"))
source(here("lib", "functions", "utility.R"))

# import data

table1_rounded <- read_delim(
  here("release20230105", "matching", "table1_rounded.csv"),
  delim = "\t", escape_double = FALSE, trim_ws = TRUE
)

km_contrasts_rounded <- read_csv(
  here("release20230207", "km_contrasts_rounded.csv")
  )

cox_contrasts_rounded <- read_csv(
  here("release20230207", "cox_contrasts_rounded.csv")
  )

km_cinc4dose_rounded <- read_csv(
  here("release20230207", "km_cinc4dose_rounded.csv")
)

flowchart_final_rounded <- read_csv(
  here("release20230207", "flowchart_final_rounded.csv")
)

# STUDY POPULATION
table1_rounded %>%
  filter(variable=="N", by == "Three doses") %>%
  pull(n) %>%
  scales::comma(accuracy = 1)

flowchart_final_rounded %>% 
  mutate(across(starts_with("n"), scales::comma, accuary=1)) %>%
  mutate(across(starts_with("pct"), ~round(100*.x, 1))) %>%
  transmute(
    crit,
    n = if_else(!is.na(n_exclude), paste0(n, " (", pct_all, "%)"), n),
    n_exclude = if_else(!is.na(n_exclude), paste0(n_exclude, " (", pct_exclude, "%)"), n_exclude)
  ) %>%
  print(n=Inf)

# create table and output to word doc

doc <- officer::read_docx() 

table1_out <- table1_rounded %>%
  # something weird happened with var_label so the following line required to fix
  mutate(across(var_label, ~ if_else(variable_levels %in% "under 18", "JCVI ageband", .x))) %>%
  filter(!(var_label %in% c("Age (per year)"))) %>%
  mutate(across(c(n, N), scales::comma, accuracy=1)) %>%
  rowwise() %>%
  transmute(
    Variable = var_label,
    Level = variable_levels,
    by, 
    value = as.character(glue(stat_display))
  ) %>%
  pivot_wider(
    names_from = by,
    values_from = value
  ) %>%
  flextable() %>%
  merge_v(j=~Variable) %>%
  width(j=1:4, width=15/4, unit="cm") 

doc <- flextable::body_add_flextable(doc, value = table1_out, split = FALSE)  

doc <- print(doc, target = here("manuscript", "table1.docx"))


# MAIN ANALYSIS

# event counts and persontime
table_fu <- km_contrasts_rounded %>%
  filter(
    filename == "overall",
    subgroup=="all", 
    variant_option == "ignore"
    ) %>%
  select(outcome, starts_with(c("persontime", "n.event"))) %>%
  mutate(across(starts_with("persontime"), ~.x/365.25)) %>%
  pivot_longer(
    cols = starts_with(c("persontime", "n.event")),
    names_pattern = "(.*)_(.)",
    names_to = c("name", "group")
  ) %>%
  group_by(outcome,name) %>%
  summarise(value=scales::comma(sum(value), accuracy = 1), .groups = "keep") %>%
  ungroup() %>%
  pivot_wider(
    names_from = name,
    values_from = value
  )
  

# risk and risk differences
table_risk <- km_contrasts_rounded %>%
  filter(
    filename == "overall",
    subgroup=="all", 
    variant_option == "ignore"
  ) %>%
  select(outcome, starts_with(c("risk", "rd"))) %>%
  select(-ends_with(".se")) %>%
  pivot_longer(
    cols = starts_with("risk"),
    names_pattern = "(.*)_(.*)",
    names_to = c("estimate", "group")
  ) %>%
  # per 1000
  mutate(across(where(is.numeric), ~format(round(x=.x*1000, digits=2), nsmall=2))) %>%
  pivot_wider(
    names_from = estimate,
    values_from = value
  ) %>%
  transmute(
    outcome, 
    group = factor(group, levels = c("0", "1"), labels = c("risk_unboosted", "risk_boosted")),
    rd = paste0(rd, " (", rd.ll, ", ", rd.ul, ")"),
    risk = paste0(risk, " (", risk.ll, ", ", risk.ul, ")")
  ) %>%
  pivot_wider(
    names_from = group,
    values_from = risk
  )

table_hr <- cox_contrasts_rounded %>%
  filter(
    filename == "overall",
    subgroup=="all", 
    variant_option == "ignore",
    model == "cox_adj",
    term == "treated"
    ) %>%
  mutate(across(starts_with("coxhr"), ~format(round(.x, 2), nsmall = 2))) %>%
  transmute(
    outcome, 
    hr = paste0(coxhr, " (", coxhr.ll, ", ", coxhr.ul, ")")
    )

table_6month_out <- table_fu %>%
  left_join(table_risk, by = "outcome") %>%
  left_join(table_hr, by = "outcome") %>%
  mutate(across(
    outcome, 
    factor, 
    levels = events_lookup$event, 
    labels = events_lookup$event_descr
    )) %>%
  arrange(outcome) %>%
  select(
    Outcome = outcome, 
    `Number of events` = n.event,
    `Person-time (years)` = persontime,
    `Risk in the unboosted` = risk_unboosted,
    `Risk in the boosted` = risk_boosted,
    `Risk difference (boosted - unboosted)` = rd,
    `Hazard ratio` = hr
  ) 

doc <- officer::read_docx() 
doc <- flextable::body_add_flextable(
  doc, 
  value = flextable(table_6month_out) %>% width(j=1:ncol(table_6month_out), width=15/ncol(table_6month_out), unit="cm"), 
  split = FALSE
  )  
doc <- print(doc, target = here("manuscript", "table_6month.docx"))


# create figure ans tables of estimates
source(here("manuscript", "functions.R"))
# plots of hazard ratios
combine_plot(subgroup_select = "all", variant_option_select = "ignore")
combine_plot(subgroup_select = "all", variant_option_select = "split")
combine_plot(subgroup_select = "agegroup", variant_option_select = "ignore")
combine_plot(subgroup_select = "prior_covid_infection", variant_option_select = "ignore")
combine_plot(subgroup_select = "cev_cv", variant_option_select = "ignore")
combine_plot(subgroup_select = "vax12_type", variant_option_select = "ignore")

# tables of hazard ratios
hr_table(subgroup_select = "all", variant_option_select = "ignore")
hr_table(subgroup_select = "all", variant_option_select = "split")
hr_table(subgroup_select = "agegroup", variant_option_select = "ignore")
hr_table(subgroup_select = "prior_covid_infection", variant_option_select = "ignore")
hr_table(subgroup_select = "cev_cv", variant_option_select = "ignore")
hr_table(subgroup_select = "vax12_type", variant_option_select = "ignore")
