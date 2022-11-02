
# # # # # # # # # # # # # # # # # # # # #
# Purpose: 
#  - import matched data
#  - adds outcome variable and restricts follow-up
#  - gets KM estimates, with covid and non covid death as competing risks
#  - COX MODELS 
#  - The script must be accompanied by three arguments:
#    `cohort` - pfizer or moderna
#    `subgroup` - prior_covid_infection, vax12_type, cev, age65plus
#    `outcome` - the dependent variable
#    `variant_option` - ignore (ignore variant era), split (split fup according to variant era), restrict (restrict fup to each variant era)

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
  variant_option <- "split" # ignore, split, restrict (delta, transition, omicron)
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


data_matched <- read_rds(ghere("output", cohort, "match", "data_matched.rds"))

if (variant_option == "restrict") {
  
  # create a duplicate row for each variant era (will filter later)
  data_matched <- data_matched %>%
    uncount(weights = nrow(variant_dates), .id="variant_id") %>%
    mutate(variant = factor(variant_id, labels = variant_dates$variant))
  
} else {
  
  data_matched <- data_matched %>%
    mutate(
      variant_id = 1L, # here variant_id is just constant
      variant = factor(variant_option)
      )
  
}

if (subgroup == "vax3_type") {
  
  # when subgroup is vax3_type, matched pairs take the vax3_type of the treated person
  data_matched <- data_matched %>%
    group_by(trial_date, match_id, matching_round) %>%
    mutate(uniquematch_id = cur_group_id()) %>% 
    ungroup() %>%
    group_by(uniquematch_id) %>%
    mutate(across(vax3_type, ~.x[treated==1])) %>%
    ungroup() %>%
    select(-uniquematch_id)
  
}

## import baseline data for matched individuals, and derive time-to-event variables
data_matched <- data_matched %>%
  # derive start and end dates for variant eras
  group_by(variant_id) %>%
  mutate(
    variantstart_date = if_else(
      variant %in% c("ignore", "split"),
      study_dates[[cohort]]$start_date,
      max(study_dates[[cohort]]$start_date, variant_dates$start_date[variant_id])
    ),
    variantend_date = if_else(
      variant %in% c("ignore", "split"),
      study_dates$studyend_date,
      min(study_dates$studyend_date, variant_dates$end_date[variant_id])
    )
  ) %>%
  ungroup() %>%
  # filter to keep only individuals with trial date during the variant era (only neede when variant_option="restrict")
  # this should take it back to 1 row per patient when variant_option="restrict", but double check this
  filter(
    trial_date >= variantstart_date,
    trial_date < variantend_date
  ) %>%
  # create a new id to account for the fact that some controls become treated (this is only needed for cox models)
  group_by(patient_id, match_id, matching_round, treated) %>% 
  mutate(new_id = cur_group_id()) %>% 
  ungroup() %>%
  mutate(all="all") %>%
  select(
    # select only variables needed for models to save space
    new_id, treated, trial_date, variantend_date, variant,
    controlistreated_date,
    vax3_date,
    death_date, dereg_date, coviddeath_date, noncoviddeath_date, vax4_date,
    all_of(covariates_model),
    all_of(c(glue("{outcome}_date"), subgroup))
  ) %>%
  mutate(

    #trial_date,
    outcome_date = .data[[glue("{outcome}_date")]],
    
    # follow-up time is up to and including censor date
    censor_date = pmin(
      dereg_date,
      # vax4_date-1, # -1 because we assume vax occurs at the start of the day
      death_date,
      variantend_date,
      trial_date + maxfup,
      na.rm=TRUE
    ),
    
    matchcensor_date = pmin(censor_date, controlistreated_date -1, na.rm=TRUE), # new censor date based on whether control gets treated or not

    tte_outcome = tte(trial_date - 1, outcome_date, matchcensor_date, na.censor=FALSE), # -1 because we assume vax occurs at the start of the day, and so outcomes occurring on the same day as treatment are assumed "1 day" long
    ind_outcome = censor_indicator(outcome_date, matchcensor_date),
    
  )


# check one row per new_id
cat("check for duplicate new_id:\n")
data_matched %>% group_by(new_id) %>% count() %>% filter(n>1) %>% nrow() %>% print()
# should always be 0

# outcome frequency
outcomes_per_treated <- table(outcome=data_matched$ind_outcome, treated=data_matched$treated)

table(
  cut(data_matched$tte_outcome, c(-Inf, 0, 1, Inf), right=FALSE, labels= c("<0", "0", ">0"))
)
# should be c(0, 0, nrow(data_matched))



## competing risks cumulative risk differences ----

# no applicable method for 'complete' applied to an object of class "c('integer', 'numeric')"

# preprocessing for data_surv ----
if (variant_option == "split") {
  
  # split follow-up time according to variant era
  # (`fup_split_variant` is also read into the `coxcontrast` function)
  fup_split_variant <-
    data_matched %>%
    select(new_id, trial_date) %>%
    uncount(weights = nrow(variant_dates), .id="variant_id") %>%
    group_by(variant_id) %>%
    mutate(
      variant = variant_dates$variant[variant_id],
      variantstart_date = variant_dates$start_date[variant_id],
      variantend_date = variant_dates$end_date[variant_id],# - 1,
      variantstart_day = as.integer(variantstart_date - trial_date),
      variantend_day = as.integer(variantend_date - trial_date)
    ) %>%
    ungroup() %>%
    # earliest variantstart_day is zero 
    mutate(across(variantstart_day, ~pmax(0, .x))) %>%
    # latest variantend_day is trial_date + maxfup
    mutate(across(variantend_day, ~pmin(maxfup, .x))) %>% 
    # filter nonsense rows after cleaning
    filter(
      # remove rows where variantend_date < trial_date
      variantend_day >= 0,
      # remove rows where cleaned variantend_day < variantstart_day
      variantend_day > variantstart_day
    ) %>%
    select(
      new_id, variant_id, variantstart_day, variantend_day
    )
  
  data_split_variant <-
    tmerge(
      data1 = data_matched,
      data2 = data_matched,
      id = new_id,
      tstart = 0,
      tstop = tte_outcome,
      ind_outcome = event(if_else(ind_outcome, tte_outcome, NA_real_))
    ) %>%
    # add post-treatment periods
    tmerge(
      data1 = .,
      data2 = fup_split_variant,
      id = new_id,
      variant_id = tdc(variantstart_day, variant_id)
    ) 
  
  data_surv <- data_split_variant %>%
    # update the variant variable from "split" to the variant using variant_id
    mutate(variant = factor(variant_dates$variant[variant_id], levels = variant_dates$variant))
  
  surv_formula <- formula(Surv(tstart, tstop, ind_outcome) ~ 1)
  
} else {
  
  data_surv <- data_matched
  
  surv_formula <- formula(Surv(tte_outcome, ind_outcome) ~ 1)
  
}

# derive data_surv ----
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
        fill(n.risk, .direction = c("up")) # fill in n.risk on each zero-event day
    }), # return survival table for each day of follow up
  ) %>%
  select(!!subgroup_sym, variant, treated, surv_obj_tidy) %>%
  unnest(surv_obj_tidy) 




km_process <- function(.data, round_by){
   
  .data %>% 
    # group_by(treated, !!subgroup_sym, variant) %>%
    mutate(
    
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
    # n.risk = roundmid_any(N, round_by) - lag(cml.eventcensor, 1, 0), # this won't work when variant_option="split", 
    # ok to use the following?
    n.risk = roundmid_any(n.risk, round_by),

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
    n.risk, n.event, n.censor,
    surv, surv.se, surv.ll, surv.ul,
    risk, risk.se, risk.ll, risk.ul
  ) %>%
    mutate(time_max = max(if_else(is.na(surv), NA_real_, time), na.rm = TRUE)) %>%
    filter(time <= time_max) %>%
    select(-time_max)
  
 }
 
 
data_surv_unrounded <- km_process(data_surv, 1)
data_surv_rounded <- km_process(data_surv, threshold)

write_rds(data_surv_unrounded, fs::path(output_dir, "km_estimates_unrounded.rds"))
write_rds(data_surv_rounded, fs::path(output_dir, "km_estimates_rounded.rds"))


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

km_plot_unrounded <- km_plot(data_surv_unrounded)
km_plot_rounded <- km_plot(data_surv_rounded)

ggsave(filename=fs::path(output_dir, "km_plot_unrounded.png"), km_plot_unrounded, width=20, height=15, units="cm")
ggsave(filename=fs::path(output_dir, "km_plot_rounded.png"), km_plot_rounded, width=20, height=15, units="cm")

# define contrast functions ----
# km
## calculate quantities relating to cumulative incidence curve and their ratio / difference / etc

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

# cox

# variant = "ignore": ignore the effect of variant
# variant = "split": split follow-up time by variants
# variant = "delta"/"transition"/"omicron": restrict to a variant era
coxcontrast <- function(data, adj = FALSE, cuts=NULL){
  
  cox_formula <- formula(Surv(tstart, tstop, ind_outcome) ~ treated)
  
  if (is.null(cuts)) {
    stop("Specify `cuts`.")
  } else if (length(cuts) < 2) {
    stop("`cuts` must specify a start and an end date")
  } 
  
  if (length(cuts) > 2 | variant_option != "ignore") {
    # stratify by fup_period if more than one follow-up period
    cox_formula <- cox_formula %>% update(as.formula(". ~ .:strata(period_id)"))
  } 
  
  # add covariates if fitting adjusted model
  if (adj) {
    cox_formula <- cox_formula %>%
      update(as.formula(str_c(". ~ . +", str_c(covariates_model, collapse = " + "))))
  }
  
  fup_period_labels <- str_c(cuts[-length(cuts)]+1, "-", lead(cuts)[-length(cuts)])
  
  data <- data %>% 
    # create variable for cuts[1] for tstart in tmerge
    mutate(time0 = cuts[1])
  
  # derive fup_split (extra processing required when variant_option %in% c("split", "restrict"))
  fup_split <-
    data %>%
    select(new_id, treated) %>%
    uncount(weights = length(cuts)-1, .id="period_id") %>%
    mutate(
      fupstart_time = cuts[period_id],
      fupend_time = cuts[period_id+1]-1,
    ) %>%
    droplevels() %>%
    select(
      new_id, period_id, fupstart_time, fupend_time# fup_time
    ) %>%
    mutate(across(period_id, factor, labels = fup_period_labels))
  
  
  # extra processing for variant_option %in% c("split", "restrict")
  if (variant_option == "split") {
    
    fup_split <- fup_split %>%
      left_join(
        fup_split_variant, by = "new_id"
      ) %>%
      # update fupstart_time and fupend_time to be within variant periods
      mutate(across(fupstart_time, ~pmax(.x, variantstart_day))) %>%
      mutate(across(fupend_time, ~pmin(.x, variantend_day))) %>%
      filter(
        # only keep rows where fup start and end times make sense
        fupstart_time <= fupend_time
      ) %>%
      mutate(variant = variant_dates$variant[variant_id]) %>%
      select(new_id, variant, period_id, fupstart_time, fupend_time)
    
  } else if (variant_option == "restrict") {
    
    fup_split <- fup_split %>%
      left_join(
        data %>% distinct(new_id, variant), by = "new_id"
      )
    
  }
  
  # add variant label for follow-up periods
  if (variant_option %in% c("split", "restrict")) {
    
    fup_period_labels <- str_c(
      rep(fup_period_labels, each = nrow(variant_dates)), 
      variant_dates$variant, 
      sep = "; "
      )
    
    fup_split <- fup_split %>%
      mutate(across(period_id, 
                    ~factor(
                      str_c(as.character(.x), variant, sep = "; "),
                      levels = fup_period_labels
                    ))) 
      
  }
  
  data_split <-
    tmerge(
      data1 = data,
      data2 = data,
      id = new_id,
      tstart = time0,
      tstop = tte_outcome,
      ind_outcome = event(if_else(ind_outcome, tte_outcome, NA_real_))
    ) %>%
    # add post-treatment periods
    tmerge(
      data1 = .,
      data2 = fup_split,
      id = new_id,
      period_id = tdc(fupstart_time, period_id)
    ) 
  
  data_cox <-
    data_split %>%
    group_by(!!subgroup_sym) %>%
    nest() %>%
    mutate(
      cox_obj = map(data, ~{
        coxph(cox_formula, data = .x, y=FALSE, robust=TRUE, id=new_id, na.action="na.fail")
      }),
      cox_obj_tidy = map(cox_obj, ~broom::tidy(.x)),
    ) %>%
    select(!!subgroup_sym, cox_obj_tidy) %>%
    unnest(cox_obj_tidy) %>%
    transmute(
      !!subgroup_sym,
      term,
      coxhazr = exp(estimate),
      coxhr.se = robust.se,
      coxhr.ll = exp(estimate + qnorm(0.025)*robust.se),
      coxhr.ul = exp(estimate + qnorm(0.975)*robust.se),
    )
  data_cox
  
}

# apply contrast functions ----

# km
km_contrasts_rounded_daily <- kmcontrasts(data_surv_rounded)
write_rds(km_contrasts_rounded_daily, fs::path(output_dir, "km_contrasts_daily_rounded.rds"))

km_contrasts_rounded_cuts <- kmcontrasts(data_surv_rounded, postbaselinecuts)
write_rds(km_contrasts_rounded_cuts, fs::path(output_dir, "km_contrasts_cuts_rounded.rds"))

km_contrasts_rounded_overall <- kmcontrasts(data_surv_rounded, c(0,maxfup))
write_rds(km_contrasts_rounded_overall, fs::path(output_dir, "km_contrasts_overall_rounded.rds"))

# cox unadjusted
cox_unadj_contrasts_cuts <- coxcontrast(data_matched, cuts = postbaselinecuts)
write_rds(cox_unadj_contrasts_cuts, fs::path(output_dir, "cox_unadj_contrasts_cuts_rounded.rds"))

cox_unadj_contrasts_overall <- coxcontrast(data_matched, cuts = c(0,maxfup))
write_rds(cox_unadj_contrasts_overall, fs::path(output_dir, "cox_unadj_contrasts_overall_rounded.rds"))

# cox adjusted
cox_adj_contrasts_cuts <- coxcontrast(data_matched, adj = TRUE, cuts = postbaselinecuts)
write_rds(cox_adj_contrasts_cuts, fs::path(output_dir, "cox_adj_contrasts_cuts_rounded.rds"))

cox_adj_contrasts_overall <- coxcontrast(data_matched, adj = TRUE, cuts = c(0,maxfup))
write_rds(cox_adj_contrasts_overall, fs::path(output_dir, "cox_adj_contrasts_overall_rounded.rds"))
