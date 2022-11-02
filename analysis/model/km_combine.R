
# # # # # # # # # # # # # # # # # # # # #
# Purpose: Combine km and cox estimates from different outcomes
#  - The script must be accompanied by three arguments:
#    `cohort` - the cohort used
# # # # # # # # # # # # # # # # # # # # #

# Preliminaries ----

## Import libraries ----
library('tidyverse')
library('here')
library('glue')
library('survival')

## Import custom user functions from lib
source(here("lib", "functions", "utility.R"))

## Import design elements
source(here("analysis", "design.R"))

## import command-line arguments ----

args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  cohort <- "mrna"
} else {
  cohort <- args[[1]]
}



output_dir <- ghere("output", cohort, "models", "km", "combined")
fs::dir_create(output_dir)

# metaparams for all models that have been run
metaparams <-
  km_args %>%
  mutate(
    outcome_descr = fct_recoderelevel(outcome,  recoder$outcome),
    subgroup_descr = fct_recoderelevel(subgroup,  recoder$subgroups),
  )

# combine and save outputs ----
combine_and_save <- function(filename) {
  
  metaparams %>%
    mutate(
      data = pmap(
        list(cohort, subgroup, variant_option, outcome), 
        function(cohort, subgroup, variant_option, outcome) {
          subgroup <- as.character(subgroup)
          dat <- read_rds(here("output", cohort, "models", "km", subgroup, variant_option, outcome, glue("{filename}.rds")))
          dat %>%
            add_column(
              subgroup_level = as.character(.[[subgroup]]),
              subgroup_level_descr = fct_recoderelevel(.[[subgroup]], recoder[[subgroup]]),
              .before=1
            ) %>%
            select(-all_of(subgroup))
        })
    ) %>%
    unnest(data) %>%
    write_csv(fs::path(output_dir, glue("{filename}.csv")))
}


# km outputs
for (i in c("estimates", "contrasts_daily", "contrasts_cuts", "contrasts_overall")) {
  combine_and_save(glue("km_{i}_rounded")) 
}
# cox outputs
for (i in c("contrasts_cuts", "contrasts_overall")) {
  combine_and_save(glue("cox_unadj_{i}_rounded")) 
  combine_and_save(glue("cox_adj_{i}_rounded")) 
}


## move km plots to single folder ----
fs::dir_create(here("output", cohort, "models", "km", "combined"))

for (i in c("rounded", "unrounded")) {
  metaparams %>%
    mutate(
      plotdir = glue("output", cohort, "models", "km", "{subgroup}", "{variant_option}", "{outcome}", "km_plot_{i}.png", .sep="/"),
      plotnewdir = glue("output", cohort, "models", "km", "combined", "km_plot_{i}_{subgroup}_{outcome}.png", .sep="/")
    ) %>%
    {walk2(.$plotdir, .$plotnewdir, ~fs::file_copy(.x, .y, overwrite = TRUE))}
}


## plot overall estimates for inspection ----

plot_estimates <- function(.data, estimate, estimate.ll, estimate.ul, name){
  
  colour_labs <- c("ignore", variant_dates$variant)
  colour_palette <- c("#636363",  RColorBrewer::brewer.pal(n=length(variant_dates$variant), name="Dark2"))
  names(colour_palette) <- colour_labs

  plot_temp <-
    .data %>%
    mutate(across(subgroup_level_descr,
                  ~if_else(
                    subgroup=="all",glue("{subgroup_level_descr} ({variant_option})"),
                    .x))) %>%
    mutate(across(variant, factor, levels = colour_labs)) %>%
    group_by(outcome_descr) %>%
    mutate(
      outcome_descr = fct_relabel(outcome_descr, str_wrap, width=10),
      subgroup_level_descr = fct_rev(subgroup_level_descr)
    ) %>%
    ggplot(aes(y=subgroup_level_descr, colour = variant)) +
    geom_vline(aes(xintercept=0), linetype="dotted", colour="darkgrey")+
    geom_point(aes(x={{estimate}}), position=position_dodge(width=-0.3), alpha=0.7)+
    geom_linerange(aes(xmin={{estimate.ll}}, xmax={{estimate.ul}}), position=position_dodge(width=-0.3))+
    facet_grid(rows=vars(subgroup_descr), cols=vars(outcome_descr), scales="free", space="free_y", switch="y")+
    scale_x_continuous(expand = expansion(mult=c(0,0.01)))+
    scale_color_manual(values=colour_palette)+
    labs(y=NULL)+
    theme_minimal()+
    theme(
      legend.position="bottom",
      axis.text.x.top=element_text(hjust=0),

      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.background = element_blank(),
      strip.placement="outside",
      #strip.text.y.left = element_text(angle=0),
      strip.text.y.left = element_blank(),

      panel.border = element_blank(),
      panel.spacing = unit(0.3, "lines"),
    )


  ggsave(
    filename=glue("output", cohort, "models", "km", "combined", "overall_plot_rounded_{name}.png", .sep="/"),
    plot_temp,
    width=20, height=15, units="cm"
  )

  plot_temp
}

# read data and create plots ----
km_contrasts_overall <- read_csv(fs::path(output_dir, glue("km_contrasts_overall_rounded.csv")))
km_contrasts_overall %>% plot_estimates(rd, rd.ll, rd.ul, "km_rd")
km_contrasts_overall %>% plot_estimates(rr, rr.ll, rr.ul, "km_rr")


cox_unadj_contrasts_overall <- read_csv(fs::path(output_dir, glue("cox_unadj_contrasts_overall_rounded.csv"))) %>%
  mutate(variant=str_extract(term, variant_dates$variant[1]))
cox_unadj_contrasts_overall %>% plot_estimates(coxhazr, coxhr.ll, coxhr.ul, "cox_unadj")

cox_adj_contrasts_overall <- read_csv(fs::path(output_dir, glue("cox_adj_contrasts_overall_rounded.csv"))) %>%
  filter(term=="treated")
cox_adj_contrasts_overall %>% plot_estimates(coxhazr, coxhr.ll, coxhr.ul, "cox_adj")

