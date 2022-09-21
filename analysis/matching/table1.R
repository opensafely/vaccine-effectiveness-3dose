# # # # # # # # # # # # # # # # # # # # #
# Purpose: describe matching results
# imports matching data
# reports on matching coverage, matching flowcharts, creates a "table 1", etc
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
  cohort <- "over12"
} else {
  #FIXME replace with actual eventual action variables
  removeobjects <- TRUE
  cohort <- args[[1]]
}


## get cohort-specific parameters study dates and parameters ----

dates <- map(study_dates[[cohort]], as.Date)
params <- study_params[[cohort]]



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
  group_by(vax1_date) %>%
  summarise(
    n_eligible = n(),
    n_matched = sum(matched, na.rm=TRUE),
  ) %>%
  mutate(
    n_unmatched = n_eligible - n_matched,
  ) %>%
  pivot_longer(
    cols = c(n_unmatched, n_matched),
    names_to = "status",
    names_prefix = "n_",
    values_to = "n"
  ) %>%
  arrange(vax1_date, status) %>%
  group_by(vax1_date, status) %>%
  summarise(
    n = sum(n),
  ) %>%
  group_by(status) %>%
  complete(
    vax1_date = full_seq(c(dates$start_date, dates$end_date), 1), # go X days before to
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
  arrange(status_descr, vax1_date)



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



## plot matching coverage ----

plot_coverage_n <-
  data_coverage %>%
  ggplot()+
  geom_col(
    aes(
      x=vax1_date+0.5,
      y=n,
      group=status,
      fill=status_descr,
      colour=NULL
    ),
    position=position_stack(reverse=TRUE),
    #alpha=0.8,
    width=1
  )+
  #geom_rect(xmin=dates$start_date, xmax= dates$end_date+1, ymin=-6, ymax=6, fill="grey", colour="transparent")+
  geom_hline(yintercept = 0, colour="black")+
  scale_x_date(
    breaks = unique(lubridate::ceiling_date(data_coverage$vax1_date, "1 month")),
    limits = c(dates$start_date-1, NA),
    labels = scales::label_date("%d/%m"),
    expand = expansion(add=1),
  )+
  scale_y_continuous(
    #labels = ~scales::label_number(accuracy = 1, big.mark=",")(abs(.x)),
    expand = expansion(c(0, NA))
  )+
  scale_fill_brewer(type="qual", palette="Set2")+
  scale_colour_brewer(type="qual", palette="Set2")+
  labs(
    x="Date",
    y="Booster vaccines per day",
    colour=NULL,
    fill=NULL,
    alpha=NULL
  ) +
  theme_minimal()+
  theme(
    axis.line.x.bottom = element_line(),
    axis.text.x.top=element_text(hjust=0),
    strip.text.y.right = element_text(angle = 0),
    axis.ticks.x=element_line(),
    legend.position = "bottom"
  )+
  NULL

plot_coverage_n

ggsave(plot_coverage_n, filename="coverage_count.png", path=output_dir)

plot_coverage_cumuln <-
  data_coverage %>%
  ggplot()+
  geom_col(
    aes(
      x=vax1_date+0.5,
      y=cumuln,
      group=status,
      fill=status_descr,
      colour=NULL
    ),
    position=position_stack(reverse=TRUE),
    width=1
  )+
  geom_rect(xmin=dates$start_date, xmax= dates$end_date+1, ymin=-6, ymax=6, fill="grey", colour="transparent")+
  scale_x_date(
    breaks = unique(lubridate::ceiling_date(data_coverage$vax1_date, "1 month")),
    limits = c(dates$start_date-1, NA),
    labels = scales::label_date("%d/%m"),
    expand = expansion(add=1),
  )+
  scale_y_continuous(
    #labels = ~scales::label_number(accuracy = 1, big.mark=",")(abs(.)),
    expand = expansion(c(0, NA))
  )+
  scale_fill_brewer(type="qual", palette="Set2")+
  scale_colour_brewer(type="qual", palette="Set2")+
  scale_alpha_discrete(range= c(0.8,0.4))+
  labs(
    x="Date",
    y="Cumulative booster vaccines",
    colour=NULL,
    fill=NULL,
    alpha=NULL
  ) +
  theme_minimal()+
  theme(
    axis.line.x.bottom = element_line(),
    axis.text.x.top=element_text(hjust=0),
    strip.text.y.right = element_text(angle = 0),
    axis.ticks.x=element_line(),
    legend.position = "bottom"
  )+
  NULL

plot_coverage_cumuln

ggsave(plot_coverage_cumuln, filename="coverage_stack.png", path=output_dir)



# table 1 style baseline characteristics ----

library('gt')
library('gtsummary')

var_labels <- list(
  N  ~ "Total N",
  treated ~ "Status",
  age ~ "Age",
  sex ~ "Sex",
  #ethnicity_combined ~ "Ethnicity",
  imd_Q5 ~ "Deprivation",
  region ~ "Region",
  
  #prior_tests_cat ~ "Number of SARS-CoV-2 tests",
  prior_covid_infection ~ "Prior documented SARS-CoV-2 infection"
) %>%
set_names(., map_chr(., all.vars))

map_chr(var_labels[-c(1,2)], ~last(as.character(.)))


# use gtsummary to obtain stnadardised table 1 data
tab_summary_baseline <-
  data_matched %>%
  mutate(
    N = 1L,
    #treated_descr = fct_recoderelevel(as.character(treated), recoder$treated),
    age = factor(age, levels=sort(unique(age)))
  ) %>%
  select(
    treated,
    all_of(names(var_labels)),
  ) %>%
  tbl_summary(
    by = treated,
    label = unname(var_labels[names(.)]),
    statistic = list(N = "{N}")
  ) 

raw_stats <- tab_summary_baseline$meta_data %>%
  select(var_label, df_stats) %>%
  unnest(df_stats)


raw_stats_redacted <- raw_stats %>%
  mutate(
    n=roundmid_any(n, 6),
    N=roundmid_any(N, 6),
    p=n/N,
    var_label = factor(var_label, levels=map_chr(var_labels[-c(1,2)], ~last(as.character(.)))),
    variable_levels = replace_na(as.character(variable_levels), "")
  ) 

write_csv(raw_stats_redacted, fs::path(output_dir, "table1.csv"))





## COPY THIS SECTION TO REPORTING REPO WHEN READY

# love / smd plot ----

data_smd <- 
  raw_stats_redacted %>%
  filter(
    variable != "N"
  ) %>%
  group_by(var_label, variable, variable_levels) %>%
  summarise(
    diff = diff(p),
    sd = sqrt(sum(p*(1-p))),
    smd = diff/sd
  ) %>%
  ungroup() %>%
  mutate(
    var_label = factor(var_label, levels=map_chr(var_labels[-c(1,2)], ~last(as.character(.)))),
    variable_card = as.numeric(var_label)%%2,
    variable_levels = replace_na(as.character(variable_levels), ""),
  ) %>%
  arrange(var_label) %>%
  mutate(
    level = fct_rev(fct_inorder(str_replace(paste(var_label, variable_levels, sep=": "),  "\\:\\s$", ""))),
    cardn = row_number()
  )

plot_smd <-
  ggplot(data_smd)+
  geom_point(aes(x=smd, y=level))+
  geom_rect(aes(alpha = variable_card, ymin = rev(cardn)-0.5, ymax =rev(cardn+0.5)), xmin = -Inf, xmax = Inf, fill='grey', colour="transparent") +
  scale_alpha_continuous(range=c(0,0.3), guide="none")+
  labs(
    x="Standardised mean difference",
    y=NULL,
    alpha=NULL
  )+
  theme_minimal() +
  theme(
    strip.placement = "outside",
    strip.background = element_rect(fill="transparent", colour="transparent"),
    strip.text.y.left = element_text(angle = 0, hjust=1),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.spacing = unit(0, "lines")
  )



raw_stats_redacted %>%
  filter(
    variable !="N"
  ) %>%
  # one column per treatment group
  pivot_wider(
    id_cols = c(var_label, variable, variable_levels),
    names_from = by,
    values_from = c(N, n, p, stat_display)
  ) %>%
  # add SMD stat
  left_join(
    data_smd,
    by=c("var_label", "variable", "variable_levels")
  ) %>%
  # format columns
  mutate(
    stat_0 = glue("{n_0} ({scales::label_number(0.1, 100)(p_0)})"),
    stat_1 = glue("{n_1} ({scales::label_number(0.1, 100)(p_1)})"),
    smd = scales::label_number(0.01)(smd),
    
  ) %>%
  select(var_label, variable_levels, stat_0, stat_1, smd) %>%
  gt(
    groupname_col = "var_label"
  ) %>%
  cols_label(
    var_label = "Variable",
    variable_levels = "",
    `stat_0` = glue("Unvaccinated (N={filter(raw_stats_redacted, variable=='N', by==0) %>% pull(n)})"),
    `stat_1` = glue("Vaccinated (N={filter(raw_stats_redacted, variable=='N', by==1) %>% pull(n)})"),
    `smd` = "Standardised mean difference"
  ) %>%
  cols_align(
    align = c("right"),
    columns =  c("stat_0", "stat_1")
  ) %>%
  cols_align(
    align = c("left"),
    columns =  c("var_label", "variable_levels")
  )# %>%
  #tab_options(row_group.as_column = TRUE)



# flowchart ----

# data_flowchart_match <-
#   read_rds(here("output", "data", "data_inclusioncriteria.rds")) %>%
#   left_join(
#     data_matchstatus %>% select(patient_id, matched),
#     by="patient_id"
#   ) %>%
#   mutate(
#     c7 = c6 & matched,
#   ) %>%
#   select(-patient_id, -matched) %>%
#   group_by(vax3_type) %>%
#   summarise(
#     across(.fns=sum)
#   ) %>%
#   pivot_longer(
#     cols=-vax3_type,
#     names_to="criteria",
#     values_to="n"
#   ) %>%
#   group_by(vax3_type) %>%
#   mutate(
#     n_exclude = lag(n) - n,
#     pct_exclude = n_exclude/lag(n),
#     pct_all = n / first(n),
#     pct_step = n / lag(n),
#     crit = str_extract(criteria, "^c\\d+"),
#     criteria = fct_case_when(
#       crit == "c0" ~ "Aged 18+ and received booster dose of BNT162b2 or mRNA-1273 between 29 October 2021 and 31 January 2022", # paste0("Aged 18+\n with 2 doses on or before ", format(study_dates$lastvax2_date, "%d %b %Y")),
#       crit == "c1" ~ "  with no missing demographic information",
#       crit == "c2" ~ "  with homologous primary vaccination course of BNT162b2 or ChAdOx1",
#       crit == "c3" ~ "  and not a health and social care worker",
#       crit == "c4" ~ "  and not a care/nursing home resident, end-of-life or housebound",
#       crit == "c5" ~ "  and no COVID-19-related events within 90 days",
#       crit == "c6" ~ "  and not admitted in hospital at time of booster",
#       crit == "c7" ~ "  and successfully matched",
#       TRUE ~ NA_character_
#     )
#   )

