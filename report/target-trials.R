library(tidyverse)
library(here)

source(here("analysis", "design.R"))


## cumulative incidence curves
km_contrasts_rounded <- read_csv(
  here("release20230207", "km_contrasts_rounded.csv")
)

km_plot_data <-
  km_contrasts_rounded %>%
  filter(
    outcome == "covidadmitted",
    variant_option == "ignore",
    subgroup == "all",
    filename == "daily"
    ) %>%
  select(
    time = period_end,
    matches(c("^risk_\\d", "^risk.ll_\\d", "^risk.ul_\\d"))
  ) %>%
  pivot_longer(
    cols = starts_with("risk"),
    names_pattern = "(.*)_(.)",
    names_to = c(".value", "group")
    ) %>%
    group_by(group) %>%
    group_modify(
      ~add_row(
        .x,
        time=0,
        risk=0,
        risk.ll=0,
        risk.ul=0,
        .before=0
      )
    ) %>%
    # use lagtime for confident intervals, use default=0 as confidence interval = c(0,0) at time 0
    # avoids error about removing x rows with missing values
    mutate(lagtime = lag(time, default=0)) %>%
    ungroup() %>%
  mutate(across(group, ~factor(.x, levels = as.character(0:1), labels = c("Unboosted", "Boosted")))) %>%
  # risk per 1000
  mutate(across(starts_with("risk"), ~ 1000*.x))

leg_title <- guide_legend(NULL)

p1 <- km_plot_data %>%
  ggplot(aes(
    group=group,
    colour=group,
    linetype=group,
    fill=group
  )) +
  geom_step(
    aes(x=time, y=risk),
    direction="vh"
  ) +
  geom_rect(
    aes(xmin=lagtime, xmax=time, ymin=risk.ll, ymax=risk.ul),
    alpha=0.1, colour="transparent"
  ) +
  scale_x_continuous(breaks = postbaselinecuts) +
  scale_linetype_manual(values = c("Unboosted" = "dashed", "Boosted" = "solid")) +
  guides(colour = leg_title, linetype = leg_title, fill = leg_title) +
  labs(
    x = "Days since trial start",
    y = "Cumulative incidence per 1,000",
    subtitle = "COVID-19 hospitalisation"
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(margin = margin(t=10)),
    axis.title.y = element_text(margin = margin(r=10)),
    legend.position = c(0.2,0.8)
  )

ggsave(
  plot = p1,
  filename = here("report", "km_full.png"),
  height = 14, width = 14, units = "cm"
)

p2 <- p1 + 
  scale_x_continuous(breaks = seq(0,14,2)) +
  coord_cartesian(xlim = c(0,14), ylim = c(0,0.5)) 
  
ggsave(
  plot = p2,
  filename = here("report", "km_28.png"),
  height = 14, width = 14, units = "cm"
)

## hazard ratios
cox_contrasts_rounded <- read_csv(
  here("release20230207", "cox_contrasts_rounded.csv")
)

cox_plot_data <-
  cox_contrasts_rounded %>%
  filter(
    outcome == "covidadmitted",
    variant_option == "ignore",
    subgroup == "all",
    model == "cox_adj",
    str_detect(term, "^treated"),
    filename == "cuts"
  ) %>%
  mutate(
    period = str_extract(term, "\\d+-\\d+"),
    period_start = as.integer(str_extract(period, "^\\d+")),
    period_end = as.integer(str_extract(period, "\\d+$")),
    midpoint = period_start + (period_end-period_start)/2
    )

position_dodge_val <- 0.8

primary_vax_y1 <- list(breaks = c(0.15, 0.25, 0.5, 0.75, 1), limits = c(0.15, 1))
primary_vax_y2 <- list(breaks = c(0, 0.25, 0.5, 0.75, 0.85))

formatpercent100 <- function(x,accuracy) {
  formatx <- scales::label_percent(accuracy)(x)
  
  if_else(
    formatx==scales::label_percent(accuracy)(1),
    paste0(">",scales::label_percent(1)((100-accuracy)/100)),
    formatx
  )
}

x_labels <- cox_plot_data %>%
  distinct(midpoint, period) %>%
  arrange(midpoint)

p3 <- cox_plot_data %>%
  ggplot(aes(x = midpoint)) +
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
  scale_x_continuous(
    name = "Days since third dose",
    breaks = x_labels$midpoint,
    labels = x_labels$period,
    limits = c(min(postbaselinecuts), max(postbaselinecuts)),
    expand = expansion(add = c(0,0))
  ) +
  scale_y_log10(
    name = "Adjusted hazard ratio (aHR)",
    breaks = primary_vax_y1[["breaks"]],
    limits = primary_vax_y1[["limits"]],
    oob = scales::oob_keep,
    sec.axis = sec_axis(
      ~(1-.),
      name="Vaccine effectiveness = 100 x (1 - aHR)",
      breaks = primary_vax_y2[["breaks"]],
      labels = function(x){formatpercent100(x, 1)}
    )
  ) +
  labs(
    subtitle = "COVID-19 hospitalisation"
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(margin = margin(t=10)),
    axis.title.y.left = element_text(margin = margin(r=10)),
    axis.title.y.right = element_text(margin = margin(l=10))
  )

ggsave(
  plot = p3,
  filename = here("report", "cox_hr.png"),
  height = 14, width = 14, units = "cm"
)
