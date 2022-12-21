######################################

# This script:
# calculates and plots the cumulative incidence of 4th dose in the treated
######################################

# Preliminaries ----

## Import libraries ----
library('tidyverse')
# library('lubridate')
# library('arrow')
library('here')
# library('glue')

## import local functions and parameters ---

source(here("analysis", "design.R"))
source(here("lib", "functions", "utility.R"))
source(here("lib", "functions", "survival.R"))

## import command-line arguments ----

args <- commandArgs(trailingOnly=TRUE)


if(length(args)==0){
  # use for interactive testing
  cohort <- "mrna"
  
} else {
  cohort <- args[[1]]
}

# create output directories ----

output_dir <- ghere("output", cohort, "models", "cinc_dose4")
fs::dir_create(output_dir)

## read and process data ----
data_dose4 <- read_rds(here("output", cohort, "match", "data_matched.rds")) %>%
  filter(treated==1) %>%
  select(patient_id, jcvi_group, vax4_date, death_date, dereg_date, controlistreated_date, trial_date) %>%
  mutate(across(
    jcvi_group, 
    ~ factor(case_when(
      as.character(.x) %in% c("08", "09", "10", "11", "12") ~ "08+",
      TRUE ~ as.character(.x)
    ),
    levels = c("01", "02", "03", "04a", "04b", "05", "06", "07", "08+"))
      )) %>%
  mutate(
    
    # follow-up time is up to and including censor date
    censor_date = pmin(
      dereg_date,
      death_date,
      study_dates$studyend_date,
      trial_date -1 + maxfup, # I think the "-1" term is needed here, because time 0 is trial_date-1, and if the -1 is not included time goes up to maxfup+1, which results in an NA row in the contrasts output
      na.rm=TRUE
    ),
    
    matchcensor_date = pmin(censor_date, controlistreated_date -1, na.rm=TRUE), # new censor date based on whether control gets treated or not
    
    tte_outcome = tte(trial_date - 1, vax4_date, matchcensor_date, na.censor=FALSE), # -1 because we assume vax occurs at the start of the day, and so outcomes occurring on the same day as treatment are assumed "1 day" long
    ind_outcome = censor_indicator(vax4_date, matchcensor_date),
    
  )

## process data_surv ----
data_surv <-
  data_dose4 %>%
  group_by(jcvi_group) %>%
  nest() %>%
  mutate(
    surv_obj = map(data, ~{
      survfit(Surv(tte_outcome, ind_outcome) ~ 1, data = .x)
    }),
    surv_obj_tidy = map(surv_obj, ~{
      broom::tidy(.x) %>%
        complete(
          time = seq_len(maxfup), # fill in 1 row for each day of follow up
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
  ) 
}

## apply km_process ----
data_surv_unrounded <- km_process(data_surv, 1)
data_surv_rounded <- km_process(data_surv, threshold)

write_csv(data_surv_rounded, fs::path(output_dir, "km_estimates_rounded.csv"))

## plot cumulative incidence of 4th dose for checking ----
km_plot <- function(.data) {
  .data %>%
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
    ggplot(aes(group=jcvi_group, colour=jcvi_group, fill=jcvi_group)) +
    geom_step(aes(x=time, y=risk), direction="vh")+
    geom_step(aes(x=time, y=risk), direction="vh", linetype="dashed", alpha=0.5)+
    geom_rect(aes(xmin=lagtime, xmax=time, ymin=risk.ll, ymax=risk.ul), alpha=0.1, colour="transparent")+
    scale_color_brewer(type="qual", palette="Set1", na.value="grey") +
    scale_fill_brewer(type="qual", palette="Set1", guide="none", na.value="grey") +
    scale_x_continuous(breaks = seq(0,600,14))+
    scale_y_continuous(expand = expansion(mult=c(0,0.01)))+
    coord_cartesian(xlim=c(0, NA))+
    labs(
      x="Days since 3rd dose",
      y="Cumulative incidence of 4th dose",
      colour="JCVI group",
      title=NULL
    )+
    theme_minimal()+
    theme(
      axis.line.x = element_line(colour = "black"),
      panel.grid.minor.x = element_blank(),
      legend.position=c(.05,.95),
      legend.justification = c(0,1),
    )
}

km_plot_unrounded <- km_plot(data_surv_unrounded)
km_plot_rounded <- km_plot(data_surv_rounded)

ggsave(filename=fs::path(output_dir, "km_plot_unrounded.png"), km_plot_unrounded, width=20, height=15, units="cm")
ggsave(filename=fs::path(output_dir, "km_plot_rounded.png"), km_plot_rounded, width=20, height=15, units="cm")

