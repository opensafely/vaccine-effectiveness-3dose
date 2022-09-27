######################################

# This script:
# imports data extracted by the cohort extractor (or dummy data)
# fills in unknown ethnicity from GP records with ethnicity from SUS (secondary care)
# tidies missing values
# standardises some variables (eg convert to factor) and derives some new ones
# organises vaccination date data to "vax X type", "vax X date" (rather than "pfizer X date", "az X date", ...)
######################################



# Preliminaries ----

## Import libraries ----
library("tidyverse")
library("lubridate")
library("arrow")
library("here")
library("glue")

## import local functions and parameters ---

source(here("analysis", "design.R"))

source(here("lib", "functions", "utility.R"))


## import command-line arguments ----

args <- commandArgs(trailingOnly = TRUE)


if (length(args) == 0) {
  # use for interactive testing
  cohort <- "over12"
  matching_round <- as.integer("2")
} else {
  # FIXME replace with actual eventual action variables
  cohort <- args[[1]]
  matching_round <- as.integer(args[[2]])
}

## get cohort-specific parameters study dates and parameters ----

dates <- map(study_dates[[cohort]], as.Date)
params <- study_params[[cohort]]

matching_round_date <- dates$control_extract_dates[matching_round]


## create output directory ----
fs::dir_create(ghere("output", cohort, "matchround{matching_round}", "process"))



# process ----

# use externally created dummy data if not running in the server
# check variables are as they should be
if (Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")) {

  # ideally in future this will check column existence and types from metadata,
  # rather than from a cohort-extractor-generated dummy data
  data_studydef_dummy <- read_feather(ghere("output", cohort, "matchround{matching_round}", "extract", "input_controlpotential.feather")) %>%
    # because date types are not returned consistently by cohort extractor
    mutate(across(ends_with("_date"), as.Date))

  data_custom_dummy <- read_feather(ghere("lib", "dummydata", "dummy_control_potential1_{cohort}.feather")) %>%
    mutate(
      msoa = sample(factor(c("1", "2")), size = n(), replace = TRUE) # override msoa so matching success more likely
    )


  not_in_studydef <- names(data_custom_dummy)[!(names(data_custom_dummy) %in% names(data_studydef_dummy))]
  not_in_custom <- names(data_studydef_dummy)[!(names(data_studydef_dummy) %in% names(data_custom_dummy))]


  if (length(not_in_custom) != 0) {
    stop(
      paste(
        "These variables are in studydef but not in custom: ",
        paste(not_in_custom, collapse = ", ")
      )
    )
  }

  if (length(not_in_studydef) != 0) {
    stop(
      paste(
        "These variables are in custom but not in studydef: ",
        paste(not_in_studydef, collapse = ", ")
      )
    )
  }

  # reorder columns
  data_studydef_dummy <- data_studydef_dummy[, names(data_custom_dummy)]

  unmatched_types <- cbind(
    map_chr(data_studydef_dummy, ~ paste(class(.), collapse = ", ")),
    map_chr(data_custom_dummy, ~ paste(class(.), collapse = ", "))
  )[(map_chr(data_studydef_dummy, ~ paste(class(.), collapse = ", ")) != map_chr(data_custom_dummy, ~ paste(class(.), collapse = ", "))), ] %>%
    as.data.frame() %>%
    rownames_to_column()


  if (nrow(unmatched_types) > 0) {
    stop(
      # unmatched_types
      "inconsistent typing in studydef : dummy dataset\n",
      apply(unmatched_types, 1, function(row) paste(paste(row, collapse = " : "), "\n"))
    )
  }

  data_extract <- data_custom_dummy
} else {
  data_extract <- read_feather(ghere("output", cohort, "matchround{matching_round}", "extract", "input_controlpotential.feather")) %>%
    # because date types are not returned consistently by cohort extractor
    mutate(across(ends_with("_date"), as.Date))
}


data_processed <- data_extract %>%
  mutate(
    sex = fct_case_when(
      sex == "F" ~ "Female",
      sex == "M" ~ "Male",
      # sex == "I" ~ "Inter-sex",
      # sex == "U" ~ "Unknown",
      TRUE ~ NA_character_
    ),

    # ethnicity_combined = if_else(is.na(ethnicity), ethnicity_6_sus, ethnicity),
    #
    # ethnicity_combined = fct_case_when(
    #   ethnicity_combined == "1" ~ "White",
    #   ethnicity_combined == "4" ~ "Black",
    #   ethnicity_combined == "3" ~ "South Asian",
    #   ethnicity_combined == "2" ~ "Mixed",
    #   ethnicity_combined == "5" ~ "Other",
    #   #TRUE ~ "Unknown",
    #   TRUE ~ NA_character_
    #
    # ),

    region = fct_collapse(
      region,
      `East of England` = "East",
      `London` = "London",
      `Midlands` = c("West Midlands", "East Midlands"),
      `North East and Yorkshire` = c("Yorkshire and The Humber", "North East"),
      `North West` = "North West",
      `South East` = "South East",
      `South West` = "South West"
    ),
    prior_covid_infection = (!is.na(postest_0_date)) | (!is.na(covidadmitted_0_date)) | (!is.na(primary_care_covid_case_0_date)),

    # latest covid event before study start
    anycovid_0_date = pmax(postest_0_date, covidemergency_0_date, covidadmitted_0_date, na.rm = TRUE),
    vax1_date = covid_vax_any_1_date,
  )



## select eligible patients and create flowchart ----


# Define selection criteria ----

data_criteria <- data_processed %>%
  transmute(
    patient_id,
    has_age = !is.na(age),
    has_sex = !is.na(sex) & !(sex %in% c("I", "U")),
    has_imd = imd_Q5 != "Unknown",
    vaccinated = vax1_date < matching_round_date,
    # has_ethnicity = !is.na(ethnicity_combined),
    has_region = !is.na(region),
    no_recentcovid30 = is.na(anycovid_0_date) | ((matching_round_date - anycovid_0_date) > 30),
    include = (
      has_age & has_sex & has_imd & # has_ethnicity &
        has_region &
        no_recentcovid30 &
        !vaccinated
    ),
  )


data_controlpotential <- data_criteria %>%
  filter(include) %>%
  select(patient_id) %>%
  left_join(data_processed, by = "patient_id") %>%
  droplevels()

write_rds(data_controlpotential, ghere("output", cohort, "matchround{matching_round}", "process", glue("data_controlpotential.rds")), compress = "gz")
