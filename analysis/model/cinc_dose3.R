library(tidyverse)
library(here)


## import command-line arguments ----

args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  cohort <- "mrna"
} else {
  cohort <- args[[1]]
}


## import local functions and parameters ---

source(here("analysis", "design.R"))
source(here("lib", "functions", "utility.R"))
source(here("lib", "functions", "survival.R"))

# create output directories ----

output_dir <- ghere("output", cohort, "models", "cinc_dose3")
fs::dir_create(output_dir)

## read and preprocess data ----

# the people who received their third dose on the first day of recruitment
data_treated <- bind_rows(
  readr::read_rds(here("output", "pfizer", "treated", "data_treatedeligible.rds")),
  readr::read_rds(here("output", "moderna", "treated", "data_treatedeligible.rds"))
) %>% 
  select(patient_id, jcvi_group, starts_with("vax3"), death_date, dereg_date) %>%
  filter(vax3_date == study_dates$pfizer$start_date)

# everyone who is eligible and not vaccinated on the first day of recruitment
data_potentialcontrol1 <- readr::read_rds(here("output", cohort, "matchround1", "process", "data_controlpotential.rds")) %>% 
  select(patient_id, jcvi_group, starts_with("vax3"), death_date, dereg_date)

# bind datasets and create ageband
data_all <- bind_rows(
  data_treated, 
  data_potentialcontrol1
) 

# check for duplicate IDs
data_all %>% group_by(patient_id) %>% summarise(n=n()) %>% summarise(`duplicate ids` = sum(n>1)) %>% print()

## prepare time to event data ----
data_tte <- data_all %>%
  mutate(
    
    # follow-up time is up to and including censor date
    censor_date = pmin(
      dereg_date,
      death_date,
      study_dates$recruitmentend_date,
      na.rm=TRUE
    ),
    
    # "-1" term so that day 1 is the first day on which there are third doses
    tte_outcome = tte(study_dates$pfizer$start_date - 1, vax3_date, censor_date, na.censor=FALSE), 
    ind_outcome = censor_indicator(vax3_date, censor_date)
    
  )

recruitment_length <- as.integer(study_dates$recruitmentend_date - study_dates$pfizer$start_date) + 1L


## process data_surv ----
data_surv <-
  data_tte %>%
  group_by(jcvi_group) %>%
  nest() %>%
  mutate(
    surv_obj = map(data, ~{
      survfit(Surv(tte_outcome, ind_outcome) ~ 1, data = .x)
    }),
    surv_obj_tidy = map(surv_obj, ~{
      broom::tidy(.x) %>%
        complete(
          time = seq_len(recruitment_length), # fill in 1 row for each day of follow up
          fill = list(n.event = 0, n.censor = 0) # fill in zero events on those days
        ) %>%
        fill(n.risk, .direction = c("up")) # fill in n.risk on each zero-event day
    }), # return survival table for each day of follow up
  ) %>%
  select(jcvi_group, surv_obj_tidy) %>%
  unnest(surv_obj_tidy)



## define km_process ----
km_process <- function(.data, round_by){
  
  .data %>% mutate(
    
    lagtime = lag(time, 1, 0),
    leadtime = lead(time, 1, max(time)+1),
    interval = time - lagtime,
    
    N = max(n.risk, na.rm=TRUE),
    
    # rounded to `round_by - (round_by/2)`
    cml.eventcensor = roundmid_any(cumsum(n.event+n.censor), round_by),
    cml.event = roundmid_any(cumsum(n.event), round_by),
    cml.censor = cml.eventcensor - cml.event,
    
    n.event = diff(c(0, cml.event)),
    n.censor = diff(c(0, cml.censor)),
    n.risk = roundmid_any(N, round_by) - lag(cml.eventcensor, 1, 0),
    
    # KM estimate for event of interest, combining censored and competing events as censored
    summand = (1/(n.risk-n.event)) - (1/n.risk), # = n.event / ((n.risk - n.event) * n.risk) but re-written to prevent integer overflow
    surv = cumprod(1 - n.event / n.risk),
    surv.se = surv * sqrt(cumsum(summand)), #greenwood's formula
    surv.ln.se = surv.se/surv,
    
    ## standard errors on log scale
    #surv.ll = exp(log(surv) + qnorm(0.025)*surv.ln.se),
    #surv.ul = exp(log(surv) + qnorm(0.975)*surv.ln.se),
    
    llsurv = log(-log(surv)),
    llsurv.se = sqrt((1 / log(surv)^2) * cumsum(summand)),
    
    ## standard errors on complementary log-log scale
    surv.ll = exp(-exp(llsurv + qnorm(0.975)*llsurv.se)),
    surv.ul = exp(-exp(llsurv + qnorm(0.025)*llsurv.se)),
    
    risk = 1 - surv,
    risk.se = surv.se,
    risk.ln.se = surv.ln.se,
    risk.ll = 1 - surv.ul,
    risk.ul = 1 - surv.ll
  ) %>% select(
    jcvi_group, time, lagtime, leadtime, interval,
    cml.event, cml.censor,
    n.risk, n.event, n.censor,
    surv, surv.se, surv.ll, surv.ul,
    risk, risk.se, risk.ll, risk.ul
  ) %>%
    group_modify(
      ~add_row(
        .x,
        time=0,
        lagtime=0,
        leadtime=1,
        #interval=1,
        surv=1,
        surv.ll=1,
        surv.ul=1,
        risk=0,
        risk.ll=0,
        risk.ul=0,
        .before=0
      )
    ) %>%
    # transform to dates
    mutate(across(c(time, lagtime), ~ study_dates$pfizer$start_date - 1 + .x)) %>%
    # cutoff time where cumulative incidence first exceeds 0.9
    mutate(cutoff_date = min(time[risk>=0.9], na.rm = TRUE)) %>%
    mutate(across(cutoff_date, ~if_else(is.na(as.integer(.x)), study_dates$recruitmentend_date, .x)))
}

## apply km_process ----
data_surv_unrounded <- km_process(data_surv, 1)
data_surv_rounded <- km_process(data_surv, threshold)



## plot cumulative incidence of 3rd dose for checking ----
km_plot <- function(.data) {
  .data %>%
    ggplot(aes(group=jcvi_group, colour=jcvi_group, fill=jcvi_group)) +
    # don't plot confidence intervals
    # geom_rect(aes(xmin=lagtime, xmax=time, ymin=risk.ll, ymax=risk.ul), alpha=0.1, colour="transparent")+
    # add vertical line at cutoff_date
    geom_vline(aes(xintercept = cutoff_date, colour = jcvi_group), linetype = "dotted") +
    geom_step(aes(x=time, y=risk), direction="vh")+
    geom_step(aes(x=time, y=risk), direction="vh", linetype="dashed", alpha=0.5)+
    geom_label(aes(x = cutoff_date, y = 1, label = jcvi_group), fill="white") +
    scale_color_discrete(name=NULL) +
    scale_fill_discrete(guide = "none") +
    scale_x_date(date_labels = "%b %Y") +
    scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
    labs(
      x="Date",
      y="Cumulative incidence of 3rd dose"
    ) +
    guides(colour = guide_legend(ncol = 8)) +
    theme_minimal() +
    theme(
      axis.line.x = element_line(colour = "black"),
      panel.grid.minor.x = element_blank(),
      legend.position="none",
      legend.justification = c(0,1),
    )
}

km_plot_unrounded <- km_plot(data_surv_unrounded)
km_plot_rounded <- km_plot(data_surv_rounded)

ggsave(filename=fs::path(output_dir, "km_plot_unrounded.png"), km_plot_unrounded, width=20, height=15, units="cm")
ggsave(filename=fs::path(output_dir, "km_plot_rounded.png"), km_plot_rounded, width=20, height=15, units="cm")

## create dataset of cut-off dates
data_recruitmentcutoff <- data_surv_unrounded %>%
  ungroup() %>%
  filter(time==cutoff_date) %>%
  select(jcvi_group, cutoff_date) %>%
  arrange(jcvi_group)

readr::write_rds(data_recruitmentcutoff, fs::path(output_dir, "data_recruitmentcutoff.rds"))
readr::write_csv(data_recruitmentcutoff, fs::path(output_dir, "data_recruitmentcutoff.csv"))

