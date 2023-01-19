library(tidyverse)
library(readr)
library(here)
library(glue)
library(flextable)

source(here("analysis", "design.R"))
source(here("lib", "functions", "utility.R"))

# STUDY POPULATION
table1_rounded <- read_delim(
  here("release20230105", "matching", "table1_rounded.csv"),
  delim = "\t", escape_double = FALSE, trim_ws = TRUE
  )

table1_rounded %>%
  filter(variable=="N", by == "Three doses") %>%
  pull(n) %>%
  scales::comma(accuracy = 1)

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
km_contrasts_overall_rounded <- read_csv(here("release20230105", "estimates", "km_contrasts_overall_rounded.csv"))

# event counts and persontime
table_fu <- km_contrasts_overall_rounded %>%
  filter(subgroup=="all", variant_option == "ignore") %>%
  select(outcome, starts_with(c("persontime", "n.event"))) %>%
  mutate(across(starts_with("persontime"), ~.x/365.25)) %>%
  pivot_longer(
    cols = starts_with(c("persontime", "n.event")),
    names_pattern = "(.*)_(.)",
    names_to = c("name", "group")
  ) %>%
  group_by(outcome,name) %>%
  summarise(value=scales::comma(sum(value), accuracy = 1)) %>%
  ungroup() %>%
  pivot_wider(
    names_from = name,
    values_from = value
  )
  

# risk and risk differences
table_risk <- km_contrasts_overall_rounded %>%
  filter(subgroup=="all", variant_option == "ignore") %>%
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

doc <- officer::read_docx() 

table_6month_out <- table_fu %>%
  left_join(table_risk, by = "outcome") %>%
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
    `Risk difference (boosted - unboosted)` = rd
  ) %>%
  mutate(`Vaccine effectiveness (100x(1 - hazard ratio))` = NA_character_) %>%
  flextable() %>%
  width(j=1:6, width=15/6, unit="cm") 

doc <- flextable::body_add_flextable(doc, value = table_6month_out, split = FALSE)  

doc <- print(doc, target = here("manuscript", "table_6month.docx"))


# cox model estimates
cox_adj_contrasts_cuts_rounded <- read_csv(here("release20230105", "estimates", "cox_adj_contrasts_cuts_rounded.csv"))

cox_adj_contrasts_cuts_rounded %>%
  filter(subgroup=="all", variant_option == "ignore") %>%
  filter(str_detect(term, "^treated")) %>%
  mutate(across(term, ~str_remove(.x, "treated:strata\\(period\\_id\\)"))) %>%
  # transform to ve
  mutate(across(starts_with("cox"), ~format(round(100*(1-.x), 1), nsmall=1, trim=TRUE))) %>%
  transmute(
    outcome, term,
    value = paste0(coxhr, "% (", coxhr.ul, ", ", coxhr.ll, ")")
  ) %>%
  pivot_wider(
    names_from = outcome,
    values_from = value
  )

doc <- officer::read_docx() 

table_hrs_all <- cox_adj_contrasts_cuts_rounded %>%
  filter(subgroup=="all", variant_option == "ignore") %>%
  filter(str_detect(term, "^treated")) %>%
  add_descr() %>%
  mutate(across(term, ~str_remove(.x, "treated:strata\\(period\\_id\\)"))) %>%
  mutate(across(starts_with("cox"), ~round(.x, 2))) %>%
  mutate(across(starts_with("cox"), ~format(.x, nsmall=2))) %>%
  transmute(
    outcome_descr, 
    `Days since booster` = term,
    value = paste0(coxhr, " (", coxhr.ll, ", ", coxhr.ul, ")")
  ) %>%
  pivot_wider(
    names_from = outcome_descr,
    values_from = value
  ) %>%
  flextable() %>%
  width(j=1:6, width=15/6, unit="cm") 

doc <- flextable::body_add_flextable(doc, value = table_hrs_all, split = FALSE)  

doc <- print(doc, target = here("manuscript", "table_hrs_all.docx"))

# figure
# A
km_estimates_rounded <- read_csv(here("release20230105", "estimates", "km_estimates_rounded.csv"))

A <- km_estimates_rounded %>%
  filter(subgroup == "all", variant_option == "ignore") %>%
  add_descr() %>%
  mutate(
    variant_descr = factor(
      variant,
      levels = c("ignore", "delta", "transition", "omicron"),
      labels = c("all variants", "delta variant", "delta-omicron transition", "omicron variant")
    )) %>%
  km_plot("variant_descr") +
  theme(legend.position = "none")

B <- cox_adj_contrasts_cuts_rounded %>%
  filter(
    subgroup == "all",
    variant_option == "ignore",
    str_detect(term, "^treated")
    ) %>%
  add_descr() %>%
  mutate(
    period_start = as.integer(str_extract(term, "\\d+")),
    period_end = as.integer(str_extract(term, "\\d+$")),
    midpoint = (period_start + period_end)/2,
    variant_descr = "all variants",
    model = "adjusted"
    ) %>%
  hr_plot() +
  theme(legend.position = "none")
  

cowplot::plot_grid(
  A, B,
  labels = c("A", "B"), align = "v"
  )

ggsave(
  here("manuscript", "main_analysis_plot.png"),
  width = 24, height = 15, units = "cm"
)


# VARIANT ANALYSIS
cox_unadj_contrasts_cuts_rounded <- read_csv(here("release20230105", "estimates", "cox_unadj_contrasts_cuts_rounded.csv"))

cox_unadj_contrasts_cuts_rounded %>%
  filter(subgroup=="all", variant_option == "split") %>%
  filter(str_detect(term, "^treated")) %>%
  mutate(variant = str_extract(term, "\\w+$")) %>%
  mutate(across(term, ~str_extract(.x, "\\d+-\\d+"))) %>%
  # transform to ve
  mutate(across(starts_with("cox"), ~signif(100*(1-.x), 3))) %>%
  transmute(
    outcome, variant, term,
    value = paste0(coxhr, "% (", coxhr.ul, ", ", coxhr.ll, ")")
  ) %>%
  pivot_wider(
    names_from = variant,
    values_from = value
  ) %>% print(n=Inf)


km_estimates_rounded <- read_csv(here("release20230105", "estimates", "km_estimates_rounded.csv"))

km_estimates_rounded %>%
  filter(subgroup=="all", variant_option=="split", outcome == "postest") %>%
  group_by(variant) %>%
  summarise(max(time))



