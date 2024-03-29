
# # # # # # # # # # # # # # # # # # # # #
# Purpose: 
#  - import matched data
#  - adds outcome variable and restricts follow-up
#  - gets KM estimates, with covid and non covid death as competing risks
#  - The script must be accompanied by three arguments:
#    `cohort` - pfizer or moderna
#    `subgroup` - prior_covid_infection, vax12_type, cev, age65plus
#    `outcome` - the dependent variable
#    `variant_option` - ignore (ignore variant era), 
#                     - split (split fup according to variant era), 
#                     - restrict (restrict fup to each variant era)

# # # # # # # # # # # # # # # # # # # # #

# Preliminaries ----


## Import libraries ----
library('tidyverse')
library('here')
library('glue')
library('survival')


## import local functions and parameters ---

source(here("analysis", "design.R"))
source(here("lib", "functions", "utility.R"))
source(here("lib", "functions", "survival.R"))


# import command-line arguments ----

args <- commandArgs(trailingOnly=TRUE)


if(length(args)==0){
  # use for interactive testing
  cohort <- "mrna"
  subgroup <- "all"
  variant_option <- "ignore" # ignore, split, restrict (delta, transition, omicron)
  outcome <- "postest"
  
} else {
  cohort <- args[[1]]
  subgroup <- args[[2]]
  variant_option <- args[[3]]
  outcome <- args[[4]]
}

if (subgroup!="all" & variant_option != "ignore") 
  stop("Must set `variant`=\"ignore\" for subgroup analyses.")

# derive symbolic arguments for programming with

cohort_sym <- sym(cohort)
subgroup_sym <- sym(subgroup)

# create output directories ----

output_dir <- ghere("output", cohort, "models", "km", subgroup, variant_option, outcome)
fs::dir_create(output_dir)


# read and process data_matched ----
source(here("analysis", "model", "process_data_model.R"))


## competing risks cumulative risk differences ----

# derive data_surv ----
cat("---- start data_surv ----\n")
data_surv <- data_surv %>%
  # this grouping is kept for passing into the km_* functions
  group_by(treated, !!subgroup_sym, variant) %>%
  nest() %>%
  mutate(
    surv_obj = map(data, ~{
      survfit(surv_formula, data = .x)
    }),
    surv_obj_tidy = map(surv_obj, ~{
      broom::tidy(.x) %>%
        complete(
          time = seq_len(maxfup), # fill in 1 row for each day of follow up
          fill = list(n.event = 0, n.censor = 0) # fill in zero events on those days
        ) %>%
        fill(n.risk, .direction = c("up")) %>% # fill in n.risk on each zero-event day
        mutate(
          # calculate new entries to risk set (people with delayed entry)
          n.entry = n.risk - lag(n.risk-n.event-n.censor,1,0)
        )
    }), # return survival table for each day of follow up
  ) %>%
  select(!!subgroup_sym, variant, treated, surv_obj_tidy) %>%
  unnest(surv_obj_tidy) 

cat("---- end data_surv ----\n")


# define km_process function ----
km_process <- function(.data, round_by){
   
  .data %>% 
    # group_by(treated, !!subgroup_sym, variant) %>%
    mutate(
    
    lagtime = lag(time, 1, 0),
    leadtime = lead(time, 1, max(time)+1),
    interval = time - lagtime,
     
    N = max(n.risk, na.rm=TRUE),
    
    # rounded to `round_by - (round_by/2)`
    cml.entry = roundmid_any(cumsum(n.entry), round_by),
    cml.eventcensor = roundmid_any(cumsum(n.event+n.censor), round_by),
    cml.event = roundmid_any(cumsum(n.event), round_by),
    cml.censor = cml.eventcensor - cml.event,

    n.event = diff(c(0, cml.event)),
    n.censor = diff(c(0, cml.censor)),
    # n.risk = roundmid_any(N, round_by) - lag(cml.eventcensor, 1, 0), # this won't work when variant_option="split", need to incorporate delayed entries
    n.risk = cml.entry - lag(cml.eventcensor, 1, 0), # (sum of new entries) minus (sum of people no longer under follow-up)
    n.entry = diff(c(0, cml.entry)),

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
    !!subgroup_sym, variant, treated, time, lagtime, leadtime, interval,
    cml.event, cml.censor,
    n.entry, n.risk, n.event, n.censor,
    surv, surv.se, surv.ll, surv.ul,
    risk, risk.se, risk.ll, risk.ul
  ) %>%
    mutate(time_max = max(if_else(is.na(surv), NA_real_, time), na.rm = TRUE)) %>%
    filter(time <= time_max) %>%
    select(-time_max)
  
 }
 
# apply function
cat("---- start data_surv_unrounded ----\n")
data_surv_unrounded <- km_process(data_surv, 1)
write_csv(data_surv_unrounded, fs::path(output_dir, "km_estimates_unrounded.csv"))
cat("---- end data_surv_unrounded ----\n")

cat("---- start data_surv_rounded ----\n")
data_surv_rounded <- km_process(data_surv, threshold)
write_csv(data_surv_rounded, fs::path(output_dir, "km_estimates_rounded.csv"))
cat("---- end data_surv_rounded ----\n")

# define km_plot function ----
km_plot <- function(.data) {
  
  # define variable to facet by
  if (subgroup == "all") facet_sym <- sym("variant") else facet_sym <- subgroup_sym 
  
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
    mutate(
      treated_descr = fct_recoderelevel(treated, recoder$treated),
    ) %>%
    ggplot(aes(group=treated_descr, colour=treated_descr, fill=treated_descr)) +
    geom_step(aes(x=time, y=risk), direction="vh")+
    geom_step(aes(x=time, y=risk), direction="vh", linetype="dashed", alpha=0.5)+
    geom_rect(aes(xmin=lagtime, xmax=time, ymin=risk.ll, ymax=risk.ul), alpha=0.1, colour="transparent")+
    facet_grid(rows=vars(!!facet_sym))+
    scale_color_brewer(type="qual", palette="Set1", na.value="grey") +
    scale_fill_brewer(type="qual", palette="Set1", guide="none", na.value="grey") +
    scale_x_continuous(breaks = seq(0,600,14))+
    scale_y_continuous(expand = expansion(mult=c(0,0.01)))+
    coord_cartesian(xlim=c(0, NA))+
    labs(
      x="Days",
      y="Cumulative incidence",
      colour=NULL,
      title=NULL
    )+
    theme_bw()+
    theme(
      axis.line.x = element_line(colour = "black"),
      panel.grid.minor.x = element_blank(),
      legend.position=c(.05,.95),
      legend.justification = c(0,1),
    )
}

# apply function
cat("---- start km_plot_unrounded ----\n")
km_plot_unrounded <- km_plot(data_surv_unrounded)
ggsave(filename=fs::path(output_dir, "km_plot_unrounded.png"), km_plot_unrounded, width=20, height=15, units="cm")
cat("---- end km_plot_unrounded ----\n")

cat("---- start km_plot_rounded ----\n")
km_plot_rounded <- km_plot(data_surv_rounded)
ggsave(filename=fs::path(output_dir, "km_plot_rounded.png"), km_plot_rounded, width=20, height=15, units="cm")
cat("---- end km_plot_rounded ----\n")


# define km contrast function ----
# calculate quantities relating to cumulative incidence curve and their ratio / difference / etc
kmcontrasts <- function(data, cuts=NULL){

  # if cuts=NULL then function provides daily estimates
  # if eg c(0,14,28,42,...) then follow up is split on these days
  # c(0, 140)
  
  if(is.null(cuts)){cuts <- unique(c(0,data$time))}

  data %>%
    filter(time!=0) %>%
    transmute(
      !!subgroup_sym, variant,
      treated,

      time, lagtime, interval,
      period_start = as.integer(as.character(cut(time, cuts, right=TRUE, label=cuts[-length(cuts)]))),
      period_end = as.integer(as.character(cut(time, cuts, right=TRUE, label=cuts[-1]))),
      period = cut(time, cuts, right=TRUE, label=paste0(cuts[-length(cuts)]+1, " - ", cuts[-1])),

      n.atrisk = n.risk,
      n.event, n.censor,

      cml.persontime = cumsum(n.atrisk*interval),
      cml.event = cumsum(replace_na(n.event, 0)),
      cml.censor = cumsum(replace_na(n.censor, 0)),

      rate = n.event / n.atrisk,
      cml.rate = cml.event / cml.persontime,

      surv, surv.se, surv.ll, surv.ul,
      risk, risk.se, risk.ll, risk.ul,

      inc = -(surv-lag(surv,1,1))/lag(surv,1,1),

      inc2 = diff(c(0,-log(surv)))

    ) %>%
    group_by(!!subgroup_sym, variant, treated, period_start, period_end, period) %>%
    summarise(

      ## time-period-specific quantities

      persontime = sum(n.atrisk*interval), # total person-time at risk within time period

      inc = weighted.mean(inc, n.atrisk*interval),
      inc2 = weighted.mean(inc2, n.atrisk*interval),

      n.atrisk = first(n.atrisk), # number at risk at start of time period
      n.event = sum(n.event, na.rm=TRUE), # number of events within time period
      n.censor = sum(n.censor, na.rm=TRUE), # number censored within time period

      inc = n.event/persontime, # = weighted.mean(kmhaz, n.atrisk*interval), incidence rate. this is equivalent to a weighted average of the hazard ratio, with time-exposed as the weights

      interval = sum(interval), # width of time period

      ## quantities calculated from time zero until end of time period
      # these should be the same as the daily values as at the end of the time period


      surv = last(surv),
      surv.se = last(surv.se),
      surv.ll = last(surv.ll),
      surv.ul = last(surv.ul),

      risk = last(risk),
      risk.se = last(risk.se),
      risk.ll = last(risk.ll),
      risk.ul = last(risk.ul),

      
      #cml.haz = last(cml.haz),  # cumulative hazard from time zero to end of time period

      cml.rate = last(cml.rate), # event rate from time zero to end of time period

      # cml.persontime = last(cml.persontime), # total person-time at risk from time zero to end of time period
       cml.event = last(cml.event), # number of events from time zero to end of time period
      # cml.censor = last(cml.censor), # number censored from time zero to end of time period

      # cml.summand = last(cml.summand), # summand used for estimation of SE of survival

      .groups="drop"
    ) %>%
    ungroup() %>%
    pivot_wider(
      id_cols= all_of(c(subgroup, "variant", "period_start", "period_end", "period",  "interval")),
      names_from=treated,
      names_glue="{.value}_{treated}",
      values_from=c(

        persontime, n.atrisk, n.event, n.censor,
        inc, inc2,

        surv, surv.se, surv.ll, surv.ul,
        risk, risk.se, risk.ll, risk.ul,


        cml.event, cml.rate
        )
    ) %>%
    mutate(
      n.nonevent_0 = n.atrisk_0 - n.event_0,
      n.nonevent_1 = n.atrisk_1 - n.event_1,

      ## time-period-specific quantities

      # incidence rate ratio
      irr = inc_1 / inc_0,
      irr.ln.se = sqrt((1/n.event_0) + (1/n.event_1)),
      irr.ll = exp(log(irr) + qnorm(0.025)*irr.ln.se),
      irr.ul = exp(log(irr) + qnorm(0.975)*irr.ln.se),


    # incidence rate ratio, v2
      irr2 = inc2_1 / inc2_0,
      irr2.ln.se = sqrt((1/n.event_0) + (1/n.event_1)),
      irr2.ll = exp(log(irr2) + qnorm(0.025)*irr2.ln.se),
      irr2.ul = exp(log(irr2) + qnorm(0.975)*irr2.ln.se),

      # incidence rate difference
      #ird = rate_1 - rate_0,

      ## quantities calculated from time zero until end of time period
      # these should be the same as values calculated on each day of follow up


      # cumulative incidence rate ratio
      cmlirr = cml.rate_1 / cml.rate_0,
      cmlirr.ln.se = sqrt((1/cml.event_0) + (1/cml.event_1)),
      cmlirr.ll = exp(log(cmlirr) + qnorm(0.025)*cmlirr.ln.se),
      cmlirr.ul = exp(log(cmlirr) + qnorm(0.975)*cmlirr.ln.se),

      # survival ratio, standard error, and confidence limits
      sr = surv_1 / surv_0,
      #cisr.ln = log(cisr),
      sr.ln.se = (surv.se_0/surv_0) + (surv.se_1/surv_1), #because cmlhaz = -log(surv) and cmlhaz.se = surv.se/surv
      sr.ll = exp(log(sr) + qnorm(0.025)*sr.ln.se),
      sr.ul = exp(log(sr) + qnorm(0.975)*sr.ln.se),

      # risk ratio, standard error, and confidence limits, using delta method
      rr = risk_1 / risk_0,
      #cirr.ln = log(cirr),
      rr.ln.se = sqrt((risk.se_1/risk_1)^2 + (risk.se_0/risk_0)^2),
      rr.ll = exp(log(rr) + qnorm(0.025)*rr.ln.se),
      rr.ul = exp(log(rr) + qnorm(0.975)*rr.ln.se),

      # risk difference, standard error and confidence limits, using delta method
      rd = risk_1 - risk_0,
      rd.se = sqrt( (risk.se_0^2) + (risk.se_1^2) ),
      rd.ll = rd + qnorm(0.025)*rd.se,
      rd.ul = rd + qnorm(0.975)*rd.se,



      # cumulative incidence rate difference
      #cmlird = cml.rate_1 - cml.rate_0
    )
}


# apply function
cat("---- start km_contrasts_rounded_daily ----\n")
km_contrasts_rounded_daily <- kmcontrasts(data_surv_rounded)
write_csv(km_contrasts_rounded_daily, fs::path(output_dir, "km_contrasts_daily_rounded.csv"))
cat("---- end km_contrasts_rounded_daily ----\n")

cat("---- start km_contrasts_rounded_cuts ----\n")
km_contrasts_rounded_cuts <- kmcontrasts(data_surv_rounded, postbaselinecuts)
write_csv(km_contrasts_rounded_cuts, fs::path(output_dir, "km_contrasts_cuts_rounded.csv"))
cat("---- end km_contrasts_rounded_cuts ----\n")

cat("---- start km_contrasts_rounded_overall ----\n")
km_contrasts_rounded_overall <- kmcontrasts(data_surv_rounded, c(0,maxfup))
write_csv(km_contrasts_rounded_overall, fs::path(output_dir, "km_contrasts_overall_rounded.csv"))
cat("---- end km_contrasts_rounded_overall ----\n")

cat("script complete")
