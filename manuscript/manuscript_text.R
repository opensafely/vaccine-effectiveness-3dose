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
    names_from = outcome,
    values_from = value
  ) %>%
  mutate(
    across(
      name, 
      factor, 
      levels = c("n.event", "persontime"), 
      labels = c("Event count", "Person-time (years)")
      )
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
  mutate(across(where(is.numeric), ~format(round(x=.x*1000, digits=2), nsmall=2, trim = TRUE))) %>%
  pivot_wider(
    names_from = estimate,
    values_from = value
  ) %>%
  transmute(
    outcome, 
    group,
    rd = paste0(rd, " (", rd.ll, ", ", rd.ul, ")"),
    risk = paste0(risk, " (", risk.ll, ", ", risk.ul, ")")
  ) %>%
  pivot_longer(
    cols = c(rd, risk)
  ) %>%
  pivot_wider(
    names_from = outcome,
    values_from = value
  ) %>%
  mutate(
    across(
      group, 
      ~ factor(
        if_else(name == "rd", "boosted_unboosted", as.character(.x)),
        levels = c("0", "1", "boosted_unboosted"),
        labels = c("Risk in the unboosted", "Risk in the boosted", "Risk difference (boosted - unboosted)")
        )
      )
    ) %>%
  distinct() %>%
  select(-name) %>%
  rename(name = group) %>%
  arrange(name)

table_hr <- cox_contrasts_rounded %>%
  filter(
    filename == "overall",
    subgroup=="all", 
    variant_option == "ignore",
    model == "cox_adj",
    (term == "treated" | str_detect(term, "time_since_infection"))
    ) %>%
  mutate(across(starts_with("coxhr"), ~format(round(.x, 2), nsmall = 2))) %>%
  transmute(
    outcome, 
    name = term,
    hr = paste0(coxhr, " (", coxhr.ll, ", ", coxhr.ul, ")")
    ) %>%
  pivot_wider(
    names_from = outcome,
    values_from = hr
  ) %>%
  mutate(
    across(
      name,
      factor,
      levels = c("treated", "time_since_infection31-90 days", "time_since_infection91+ days"),
      labels = c("HR for boosted vs unboosted", "HR for 31-90 days since prior infection vs no prior infection", "HR for 91+ days since prior infection vs no prior infection")
    )
  )
  

outcome_descr <- tibble(event=outcomes) %>%
  left_join(events_lookup) %>%
  pull(event_descr)
outcome_named <- outcomes
names(outcome_named) <- outcome_descr

table_6month_out <- bind_rows(table_fu, table_risk, table_hr) %>%
  select(name, all_of(outcome_named))

doc <- officer::read_docx() 
doc <- flextable::body_add_flextable(
  doc, 
  value = flextable(table_6month_out) %>% width(j=1:ncol(table_6month_out), width=26/ncol(table_6month_out), unit="cm"), 
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
for (s in subgroups) {
  combine_plot(subgroup_select = s, variant_option_select = "ignore")
  if (s =="all") {
    combine_plot(subgroup_select = s, variant_option_select = "split")
  }
  for (mt in c("cox_unadj", "cox_adj")) {
    cat(glue("{s}; ignore; {mt};"), "\n")
    hr_table(subgroup_select = s, variant_option_select = "ignore", model_type = mt)
    if (s == "all") {
      cat(glue("{s}; split; {mt};"), "\n")
      hr_table(subgroup_select = s, variant_option_select = "split", model_type = mt)
    }
  }
}

