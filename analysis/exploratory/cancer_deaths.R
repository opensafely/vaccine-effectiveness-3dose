# knit the rmarkdown file

outdir <- here::here("output", "mrna", "exploratory")
fs::dir_create(outdir)

rmarkdown::render(
  input = here::here("analysis", "exploratory", "cancer_deaths.Rmd"),
  output_file = "cancer_deaths.html",
  output_dir = outdir
  )
