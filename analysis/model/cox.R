
# # # # # # # # # # # # # # # # # # # # #
# Purpose: 
#  - import matched data
#  - adds outcome variable and restricts follow-up
#  - Fit Cox models
#  - The script must be accompanied by three arguments:
#    `cohort` - pfizer or moderna
#    `type` - unadj, adj (Cox model adjustment)
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
  type <- "unadj"
  subgroup <- "all"
  variant_option <- "ignore" # ignore, split, restrict (delta, transition, omicron)
  outcome <- "postest"
  
} else {
  cohort <- args[[1]]
  type <- args[[2]]
  subgroup <- args[[3]]
  variant_option <- args[[4]]
  outcome <- args[[5]]
}

if (subgroup!="all" & variant_option != "ignore") 
  stop("Must set `variant`=\"ignore\" for subgroup analyses.")

# create output directories ----
output_dir <- ghere("output", cohort, "models", "cox_{type}", subgroup, variant_option, outcome)
fs::dir_create(output_dir)

# derive symbolic arguments for programming with

cohort_sym <- sym(cohort)
subgroup_sym <- sym(subgroup)

# read and process data_matched ----
source(here("analysis", "model", "process_data_model.R"))

# cox models ----

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
        coxph(cox_formula, data = .x, y=FALSE, robust=TRUE, id=patient_id, na.action="na.fail")
      }),
      cox_obj_tidy = map(cox_obj, ~broom::tidy(.x)),
    ) %>%
    select(!!subgroup_sym, cox_obj_tidy) %>%
    unnest(cox_obj_tidy) %>%
    transmute(
      !!subgroup_sym,
      term,
      coxhr = exp(estimate),
      coxhr.se = robust.se,
      coxhr.ll = exp(estimate + qnorm(0.025)*robust.se),
      coxhr.ul = exp(estimate + qnorm(0.975)*robust.se),
    )
  data_cox
  
}

# apply contrast functions ----

# cox unadjusted
if (type == "unadj") {
  
  cat("---- start cox_unadj_contrasts_cuts ----\n")
  cox_unadj_contrasts_cuts <- coxcontrast(data_matched, cuts = postbaselinecuts)
  write_rds(cox_unadj_contrasts_cuts, fs::path(output_dir, "cox_unadj_contrasts_cuts_rounded.rds"))
  cat("---- end cox_unadj_contrasts_cuts ----\n")
  
  cat("---- start cox_unadj_contrasts_overall ----\n")
  cox_unadj_contrasts_overall <- coxcontrast(data_matched, cuts = c(0,maxfup))
  write_rds(cox_unadj_contrasts_overall, fs::path(output_dir, "cox_unadj_contrasts_overall_rounded.rds"))
  cat("---- end cox_unadj_contrasts_overall ----\n")
  
}

# cox adjusted
if (type == "adj") {
  
  cat("---- start cox_adj_contrasts_cuts ----\n")
  cox_adj_contrasts_cuts <- coxcontrast(data_matched, adj = TRUE, cuts = postbaselinecuts)
  write_rds(cox_adj_contrasts_cuts, fs::path(output_dir, "cox_adj_contrasts_cuts_rounded.rds"))
  cat("---- end cox_adj_contrasts_cuts ----\n")
  
  cat("---- start cox_adj_contrasts_overall ----\n")
  cox_adj_contrasts_overall <- coxcontrast(data_matched, adj = TRUE, cuts = c(0,maxfup))
  write_rds(cox_adj_contrasts_overall, fs::path(output_dir, "cox_adj_contrasts_overall_rounded.rds"))
  cat("---- end cox_adj_contrasts_overall ----\n")
  
}

cat("script complete")
