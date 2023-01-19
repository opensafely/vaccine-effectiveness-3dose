plot_params <- function(.data, colour_var) {
  
  colour_palette <- RColorBrewer::brewer.pal(n=4, name="Dark2")
  linetype_palette <- c("Unboosted" = "dashed", "Boosted" = "solid")
  variants <- NULL
  subgroup <- NULL
  
  if (colour_var == "variant_descr") {
    
    variants <- levels(.data$variant_descr)
    names(colour_palette) <- c("all variants", "delta variant", "delta-omicron transition", "omicron variant")
    colour_palette <- colour_palette[variants]
    
  } else if (colour_var == "subgroup_level") {
    
    subgroup <- levels(.data$subgroup)
    colour_palette = colour_palette[seq_along(recoder[[subgroup]])]
    names(colour_palette) <- names(recoder[[subgroup]])
    
    .data <- .data %>%
      select(-subgroup) %>%
      mutate(
        across(subgroup_level, 
               factor,
               levels = unname(recoder[[subgroup]]), 
               labels =  names(recoder[[subgroup]])
        )
      )
    
  }
  
  .data <- .data %>%
    mutate(across(
      outcome_descr,
      factor,
      levels = levels(.data$outcome_descr),
      labels = str_replace(levels(.data$outcome_descr), "\\s", "\\\n")
    ))
  
  return(
    list(
      .data = .data,
      colour_palette = colour_palette,
      linetype_palette = linetype_palette,
      variants = variants,
      subgroup = subgroup
    )
  )
  
}


################################################################################
plot_theme <- function(...) {
  
  theme_bw() +
    theme(
      panel.border = element_blank(),
      axis.line.y = element_line(colour = "black"),
      
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.spacing = unit(0.8, "lines"),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 90),
      
      plot.title = element_text(hjust = 0),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.caption = element_text(hjust = 0, face= "italic"),
      
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y.left = element_text(margin = margin(r=10)),
      axis.title.y.right = element_text(margin = margin(l=10)),
      legend.box = "vertical"
      
    )
  
  
}

################################################################################

km_plot <- function(.data, colour_var="variant_descr"#, 
                    # end_day = max(postbaselinecuts), 
                    # padding=7, breaks = 14
                    ) {
  
  list2env(.data %>% plot_params(colour_var), envir = environment())
  
  n_levs <- .data %>% distinct(!! sym(colour_var)) %>% nrow()

  if (n_levs == 1) leg.pos <- "none" else leg.pos <- "bottom"

  p <- .data %>%
    # filter(time <= end_day) %>%
    group_by(treated_descr, outcome_descr, !! sym(colour_var)) %>%
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
    ungroup() %>%
    mutate(across(starts_with("risk"), ~1000*.x)) %>%
    ggplot(aes(
      group=paste0(treated_descr, !! sym(colour_var)),
      colour = !! sym(colour_var),
      fill = !! sym(colour_var),
      linetype = treated_descr
    )) +
    geom_step(
      aes(x=time, y=risk),
      direction="vh"
    ) +
    geom_rect(
      aes(xmin=lagtime, xmax=time, ymin=risk.ll, ymax=risk.ul),
      alpha=0.1, colour="transparent"
    ) +
    facet_grid(
      rows = vars(outcome_descr),
      switch = "y",
      scales = "free_y"
    ) +
    scale_color_manual(name = NULL, values = colour_palette) +
    scale_fill_manual(name = NULL, guide="none", values = colour_palette) +
    scale_linetype_manual(name = NULL, values = linetype_palette, guide = "none") +
    scale_x_continuous(
      breaks = postbaselinecuts,
      expand = expansion(
        add = c(0,0)
        # add=padding
        )
    ) +
    scale_y_continuous(
      # labels = ~scales::label_number(accuracy = 0.1)(1000*.x),
      expand = expansion(mult=c(0,0.01))
    ) +
    coord_cartesian(xlim=c(0, NA)) +
    labs(
      x="Days since third dose",
      y="Cumulative incidence (x 1,000)",
    ) +
    guides(nrow = 2) +
    plot_theme() +
    theme(
      legend.position=leg.pos,
      axis.line.x = element_line(colour = "black")
    )
  
  return(p)
  
}

#####################################################################################
hr_plot <- function(.data, colour_var="variant_descr") {
  
  list2env(.data %>% plot_params(colour_var), envir = environment())
  
  n_levs <- .data %>% distinct(!! sym(colour_var)) %>% nrow()
  
  if (n_levs == 1) leg.pos <- "none" else leg.pos <- "bottom"
  
  position_dodge_val <- 12
  
  primary_vax_y1 <- list(breaks = c(0.02, 0.05, 0.2, 0.5, 1, 2), limits = c(0.02, 2))
  primary_vax_y2 <- list(breaks = c(0,0.5,0.8, 0.95, 0.98))
  
  formatpercent100 <- function(x,accuracy) {
    formatx <- scales::label_percent(accuracy)(x)
    
    if_else(
      formatx==scales::label_percent(accuracy)(1),
      paste0(">",scales::label_percent(1)((100-accuracy)/100)),
      formatx
    )
  }
  
  model_levels <- levels(.data$model)
  if (length(model_levels) == 1) {
    alpha_palette <- c("unadjusted" = 0.8, "adjusted" = 0.8)[model_levels]
    scale_alpha_custom <- scale_alpha_manual(guide = "none", values = alpha_palette)
  } else {
    alpha_palette <- c("unadjusted" = 0.4, "adjusted" = 0.8)[model_levels]
    scale_alpha_custom <- scale_alpha_manual(name = NULL, values = alpha_palette)
  }
  
  p <- .data %>%
    ggplot(
      aes(x = midpoint, colour = !! sym(colour_var), alpha = model)
    ) +
    geom_hline(aes(yintercept=1), colour='grey') +
    geom_linerange(
      aes(ymin = coxhr.ll, ymax = coxhr.ul),
      position = position_dodge(width = position_dodge_val)
    ) +
    geom_point(
      aes(y = coxhr),
      position = position_dodge(width = position_dodge_val),
      size = 2
    ) +
    facet_grid(rows = vars(outcome_descr), switch = "y", scales = "free", space = "free_x") +
    scale_x_continuous(
      name = "Days since third dose",
      breaks = postbaselinecuts,
      limits = c(min(postbaselinecuts), max(postbaselinecuts)),
      expand = expansion(
        add = c(0,0)
        # add=c(
        #   min(.data$midpoint) + 7,
        #   postbaselinecuts[length(postbaselinecuts)] - max(.data$midpoint) + 7
        # )
      )
    ) +
    scale_y_log10(
      name = "Hazard ratio",
      breaks = primary_vax_y1[["breaks"]],
      limits = primary_vax_y1[["limits"]],
      oob = scales::oob_keep,
      sec.axis = sec_axis(
        ~(1-.),
        name="Vaccine effectiveness = 100 x (1 - hazard ratio)",
        breaks = primary_vax_y2[["breaks"]],
        labels = function(x){formatpercent100(x, 1)}
      )
    ) +
    scale_color_manual(name = NULL, values = colour_palette) +
    scale_alpha_custom +
    # scale_shape_manual(name = NULL, values = shape_palette) +
    # theme_bw(base_size = 12) +
    plot_theme() +
    theme(
      legend.position=leg.pos
    )
  
  # ggsave(
  #   file.path(outdir, glue::glue("cox_{filename}.png")),
  #   plot = p,
  #   height = 20, width = 15, units = "cm"
  # )
  
  return(p)
  
}

#####################################################################################

combine_plot <- function(
  subgroup_select = "all", 
  variant_option_select = "ignore"
) {
  
  if (subgroup_select == "all") {
    colour_var <- "variant_descr"
  } else {
    colour_var <- "subgroup_level"
  }
  
  A <- km_estimates_rounded %>%
    filter(
      subgroup %in% subgroup_select,
      variant_option %in% variant_option_select
      ) %>%
    add_descr() %>%
    droplevels() %>%
    km_plot(colour_var)
  
  if (variant_option_select == "split") {
    
    cox_data <- cox_unadj_contrasts_cuts_rounded %>%
      mutate(variant = str_extract(term, "delta|transition|omicron")) %>%
      mutate(across(term, ~str_trim(str_remove(.x, ";\\sdelta|;\\stransition|;\\somicron"), side="right"))) %>%
      mutate(across(variant, ~if_else(is.na(.x), "ignore", .x)))
      
    model <- "unadjusted"
    
  } else {
    
    cox_data <- cox_adj_contrasts_cuts_rounded %>%
      mutate(variant = "ignore")
    
    model <- "adjusted"
    
  }
  
  B <- cox_data %>%
    filter(
      subgroup %in% subgroup_select,
      variant_option %in% variant_option_select,
      str_detect(term, "^treated")
    ) %>%
    add_descr() %>%
    mutate(
      period_start = as.integer(str_extract(term, "\\d+")),
      period_end = as.integer(str_extract(term, "\\d+$")),
      midpoint = (period_start + period_end)/2,
      model = factor(model)
    ) %>%
    droplevels() %>%
    hr_plot(colour_var) 
  
  if ((variant_option_select != "ignore") | (subgroup_select != "all")) {
    
    plot_legend <- cowplot::get_legend(B)
    
    A <- A + theme(legend.position = "none")
    B <- B + theme(legend.position = "none")
    
  }
  
  p <- cowplot::plot_grid(
    A, B,
    labels = c("A", "B"), align = "v"
  )
  
  if ((variant_option_select != "ignore") | (subgroup_select != "all")) {
    
    p <-  cowplot::plot_grid(
     p, plot_legend,
     nrow=2, rel_heights = c(14,1)
    )
    
  }

  ggsave(
    here("manuscript", glue("plot_{subgroup_select}_{variant_option_select}.png")),
    plot = p,
    width = 24, height = 16, units = "cm"
  )

  return(p)
  
}


#####################################################################################
# table of period-specific HRs
hr_table <- function(subgroup_select, variant_option_select) {
  
  if (variant_option_select == "split") {
    grouping_var <- c("Variant era" = "variant_descr")
    cox_data <- cox_unadj_contrasts_cuts_rounded %>%
      mutate(variant = str_extract(term, "delta|transition|omicron")) %>%
      mutate(across(term, ~str_trim(str_remove(.x, ";\\sdelta|;\\stransition|;\\somicron"), side="right"))) %>%
      mutate(across(variant, ~if_else(is.na(.x), "ignore", .x)))
  } else {
    grouping_var <- c("Subgroup" = "subgroup_level")
    cox_data <- cox_adj_contrasts_cuts_rounded %>%
      mutate(across(
        subgroup_level,
        factor,
        levels = unname(recoder[[subgroup_select]]),
        labels = names(recoder[[subgroup_select]])
        ))
  }
  
  # cox model estimates
  table_data <- cox_data %>%
    filter(
      subgroup == subgroup_select, 
      variant_option == variant_option_select,
      !is.na(coxhr)
      ) %>%
    filter(str_detect(term, "^treated")) %>%
    add_descr() %>%
    mutate(across(term, ~str_remove(.x, "treated:strata\\(period\\_id\\)"))) %>%
    mutate(
      `Days since booster` = factor(
        term, 
        levels = paste0(postbaselinecuts[-length(postbaselinecuts)]+1, "-", postbaselinecuts[-1])
        )
    ) %>%
    rename(all_of(grouping_var))
  
  table_vars <- c(names(grouping_var), "outcome_descr", "Days since booster", "value")
  
  # print VE estimates
  table_data %>%
    mutate(across(starts_with("cox"), ~format(round(100*(1-.x), 1), nsmall=1, trim=TRUE))) %>%
    mutate(value = paste0(coxhr, "% (", coxhr.ul, ", ", coxhr.ll, ")")) %>%
    select(all_of(table_vars)) %>%
    pivot_wider(
      names_from = outcome_descr,
      values_from = value
    ) %>%
    arrange(!!sym(names(grouping_var)), `Days since booster`) %>%
    print(n=Inf)
  
  # save word doc table of HRs
  
  doc <- officer::read_docx() 
  
  table_hrs_all <- table_data%>%
    mutate(across(starts_with("cox"), ~round(.x, 2))) %>%
    mutate(across(starts_with("cox"), ~format(.x, nsmall=2))) %>%
    mutate(value = paste0(coxhr, " (", coxhr.ll, ", ", coxhr.ul, ")")) %>%
    select(all_of(table_vars)) %>%
    pivot_wider(
      names_from = outcome_descr,
      values_from = value
    ) %>%
    arrange(!!sym(names(grouping_var)), `Days since booster`) %>%
    flextable() %>%
    merge_v(j=1) %>%
    width(j=1:6, width=15/6, unit="cm") 
  
  doc <- flextable::body_add_flextable(doc, value = table_hrs_all, split = FALSE)  
  
  doc <- print(doc, target = here("manuscript", glue("table_hrs_{subgroup_select}_{variant_option_select}.docx")))
  
}






