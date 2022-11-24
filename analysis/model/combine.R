
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

output_dir <- ghere("output", cohort, "models", "combined")
fs::dir_create(output_dir)

# metaparams for all models that have been run
metaparams <-
  km_args %>%
  select(-model) %>%
  distinct() 

# get all subgroups
subgroups <- metaparams %>% distinct(subgroup) %>% unlist() %>% unname()

# combine and save outputs ----
combine_and_save <- function(model, filename) {
  
  filename_full <- glue("{model}_{filename}_rounded")
  
  metaparams %>%
    mutate(
      data = pmap(
        list(cohort, subgroup, variant_option, outcome), 
        function(cohort, subgroup, variant_option, outcome)  {
          dat <- try(read_rds(here("output", cohort, "models", model, subgroup, variant_option, outcome, glue("{filename_full}.rds")))) 
          if (inherits(dat, "try-error")) {
            dat <- tibble()
          } else {
            dat <- dat %>%
              add_column(
                subgroup_level = as.character(.[[subgroup]]),
                .before=1
              ) 
          }
          return(dat)
        }
          
        )
    ) %>%
    unnest(data) %>%
    select(-any_of(subgroups)) %>%
    write_csv(fs::path(output_dir, glue("{filename_full}.csv")))
}


# km outputs
for (i in c("estimates", "contrasts_daily", "contrasts_cuts", "contrasts_overall")) {
  combine_and_save(model="km", filename = i) 
}
# cox outputs
for (i in c("contrasts_cuts", "contrasts_overall")) {
  for (m in c("unadj", "adj")) {
    combine_and_save(model=glue("cox_{m}"), filename = i) 
    combine_and_save(model=glue("cox_{m}"), filename = i) 
  }
}


## move km plots to single folder ----
fs::dir_create(here("output", cohort, "models", "combined"))

for (i in c("rounded", "unrounded")) {
  metaparams %>%
    mutate(
      plotdir = glue("output", cohort, "models", "km", "{subgroup}", "{variant_option}", "{outcome}", "km_plot_{i}.png", .sep="/"),
      plotnewdir = glue("output", cohort, "models", "combined", "km_plot_{i}_{subgroup}_{outcome}.png", .sep="/")
    ) %>%
    {walk2(.$plotdir, .$plotnewdir, ~fs::file_copy(.x, .y, overwrite = TRUE))}
}


## plot overall estimates for inspection ----

plot_estimates <- function(.data, estimate, estimate.ll, estimate.ul, name){
  
  colour_labs <- c("ignore", variant_dates$variant)
  colour_palette <- c("#636363",  RColorBrewer::brewer.pal(n=length(variant_dates$variant), name="Dark2"))
  names(colour_palette) <- colour_labs
  
  subgroup_levels_labels <- unlist(unname(recoder[subgroups]))
  
  plot_data <- .data %>%
    mutate(
      subgroup = fct_recoderelevel(subgroup, recoder$subgroups),
      outcome = fct_recoderelevel(outcome,  recoder$outcome),
      subgroup_level = fct_recoderelevel(subgroup_level, subgroup_levels_labels),
      yvar = as.character(if_else(
        subgroup %in% "Main",
        glue("{subgroup}\n({str_to_sentence(variant_option)})"),
        glue("{subgroup}\n({subgroup_level})")
      )),
      .before = 1
    ) %>%
    mutate(outcome = fct_relabel(outcome, str_wrap, width=10)) 
  
  ylevels <- plot_data %>%
    distinct(yvar, subgroup, subgroup_level, variant_option) %>%
    arrange(subgroup, subgroup_level, variant_option) %>%
    select(yvar) %>% unlist() %>% unname()
  
  plot_temp <- plot_data %>%
    mutate(outcome = fct_relabel(outcome, str_wrap, width=10))  %>%
    mutate(across(variant, factor, levels = colour_labs)) %>%
    mutate(across(yvar, factor, levels = ylevels)) %>%
    ggplot(aes(y=yvar, colour = variant)) +
    geom_vline(aes(xintercept=0), linetype="dotted", colour="darkgrey")+
    geom_point(aes(x={{estimate}}), position=position_dodge(width=-0.3), alpha=0.7)+
    geom_linerange(aes(xmin={{estimate.ll}}, xmax={{estimate.ul}}), position=position_dodge(width=-0.3))+
    facet_grid(rows=vars(yvar), cols=vars(outcome), scales="free", space="free_y", switch="y")+
    scale_x_continuous(expand = expansion(mult=c(0,0.01)))+
    scale_color_manual(values=colour_palette)+
    labs(y=NULL)+
    theme_bw()+
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
    filename=glue("output", cohort, "models", "combined", "overall_plot_rounded_{name}.png", .sep="/"),
    plot_temp,
    width=30, height=20, units="cm"
  )

  plot_temp
}

# read data and create plots ----
km_contrasts_overall <- read_csv(fs::path(output_dir, glue("km_contrasts_overall_rounded.csv"))) 
km_contrasts_overall %>% plot_estimates(rd, rd.ll, rd.ul, "km_rd")
km_contrasts_overall %>% plot_estimates(rr, rr.ll, rr.ul, "km_rr")


cox_unadj_contrasts_overall <- read_csv(fs::path(output_dir, glue("cox_unadj_contrasts_overall_rounded.csv"))) %>%
  mutate(variant=str_extract(term, str_c(variant_dates$variant, collapse = "|"))) %>%
  mutate(across(variant, ~if_else(is.na(.x), "ignore", .x)))
cox_unadj_contrasts_overall %>% plot_estimates(coxhr, coxhr.ll, coxhr.ul, "cox_unadj")

cox_adj_contrasts_overall <- read_csv(fs::path(output_dir, glue("cox_adj_contrasts_overall_rounded.csv"))) %>%
  filter(str_detect(term, "^treated")) %>%
  mutate(variant=str_extract(term, str_c(variant_dates$variant, collapse = "|"))) %>%
  mutate(across(variant, ~if_else(is.na(.x), "ignore", .x)))
cox_adj_contrasts_overall %>% plot_estimates(coxhr, coxhr.ll, coxhr.ul, "cox_adj")

