
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
  distinct() 

# get all subgroups
subgroups <- metaparams %>% distinct(subgroup) %>% unlist() %>% unname()

# combine and save outputs ----
combine_and_save_contrasts <- function(model_type, filenames) {
  
  metaparams <- metaparams %>%
    filter(str_detect(model, model_type)) %>%
    uncount(length(filenames), .id="filename") %>% 
    mutate(across(filename, ~filenames[.x]))
  
  metaparams %>%
    mutate(
      data = pmap(
        list(cohort, model, subgroup, variant_option, outcome, filename), 
        function(cohort, model, subgroup, variant_option, outcome, filename)  {
          dat <- try(read_csv(here("output", cohort, "models", model, subgroup, variant_option, outcome, glue("{model}_contrasts_{filename}_rounded.csv")))) 
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
    mutate(across(
      starts_with(c("surv", "risk", "inc", "cml.rate", "irr", "cmlirr", "sr", "rd", "rr", "cox")),
      round, digits=5
    )) %>%
    write_csv(fs::path(output_dir, glue("{model_type}_contrasts_rounded.csv")))
}


# km outputs
combine_and_save_contrasts(model_type="km", filenames = c("daily", "cuts", "overall")) 
# cox outputs
combine_and_save_contrasts(model_type="cox", filenames = c("cuts", "overall"))

# km for noncancer
read_csv(fs::path(output_dir, "km_contrasts_rounded.csv")) %>%
  filter(subgroup == "noncancer") %>%
  write_csv(fs::path(output_dir, "km_contrasts_noncancer_rounded.csv"))
# cox for noncancer
read_csv(fs::path(output_dir, "cox_contrasts_rounded.csv")) %>%
  filter(subgroup == "noncancer") %>%
  write_csv(fs::path(output_dir, "cox_contrasts_noncancer_rounded.csv"))

## move km plots to single folder ----
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
km_contrasts_overall <- 
  read_csv(fs::path(output_dir, glue("km_contrasts_rounded.csv"))) %>%
  filter(filename == "overall")
km_contrasts_overall %>% plot_estimates(rd, rd.ll, rd.ul, "km_rd")
km_contrasts_overall %>% plot_estimates(rr, rr.ll, rr.ul, "km_rr")


cox_contrasts_rounded <- read_csv(fs::path(output_dir, glue("cox_contrasts_rounded.csv"))) %>%
  filter(filename == "overall") %>%
  filter(str_detect(term, "^treated")) %>%
  mutate(variant=str_extract(term, str_c(variant_dates$variant, collapse = "|"))) %>%
  mutate(across(variant, ~if_else(is.na(.x), "ignore", .x)))

cox_contrasts_rounded %>% 
  filter(model == "cox_unadj") %>%
  plot_estimates(coxhr, coxhr.ll, coxhr.ul, "cox_unadj")

cox_contrasts_rounded %>% 
  filter(model == "cox_adj") %>%
  plot_estimates(coxhr, coxhr.ll, coxhr.ul, "cox_adj")
