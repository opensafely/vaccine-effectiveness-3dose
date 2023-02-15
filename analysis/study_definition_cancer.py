# Import codelists from codelists.py
import codelists

# import json module
import json

from cohortextractor import (
  StudyDefinition,
  patients,
  codelist_from_csv,
  codelist,
  filter_codes_by_category,
  combine_codelists,
  params
)

cohort = params["cohort"]

#study_parameters
with open("./lib/design/study-dates.json") as f:
  study_dates = json.load(f)

# define variables explicitly
start_date=study_dates["pfizer"]["start_date"] # start of recruitment

# Specify study defeinition
study = StudyDefinition(
  
  # Configure the expectations framework
  default_expectations={
    "date": {"earliest": "2020-01-01", "latest": "today"},
    "incidence": 0.3,
  },
  
  # This line defines the study population 
  population = patients.which_exist_in_file(f_path=f"output/{cohort}/match/data_matched_unique.csv.gz"),

  index_date = start_date,

  # The earliest/latest approach below is so we can derive "any cancer code in 5 years before trial start date"
  # without having trial start date (to have trial start date we would have to run two separate study definitions,
  # one for the treated and one for the controls)

  # Date of latest cancer code before start of recruitment
  # Hospital admission
  cancer_hospitalisation_before_date=patients.admitted_to_hospital(
    returning= "date_admitted",
    with_these_diagnoses=codelists.cancer,
    on_or_before="index_date - 1 day",
    find_last_match_in_period=True,
    date_format="YYYY-MM-DD",
    ),
  # Primary care
  cancer_primarycare_before_date=patients.with_these_clinical_events(
    combine_codelists(
        codelists.cancer_haem_snomed, 
        codelists.cancer_nonhaem_nonlung_snomed, 
        codelists.cancer_lung_snomed
        ),
    returning="date",
    on_or_before="index_date - 1 day",
    find_last_match_in_period=True,
    date_format="YYYY-MM-DD",
    ),

  # Date of earliest cancer code after start of recruitment
  # Hospital admission
  cancer_hospitalisation_after_date=patients.admitted_to_hospital(
    returning= "date_admitted",
    with_these_diagnoses=codelists.cancer,
    on_or_after="index_date",
    find_first_match_in_period=True,
    date_format="YYYY-MM-DD",
    ),
  # Primary care
  cancer_primarycare_after_date=patients.with_these_clinical_events(
    combine_codelists(
        codelists.cancer_haem_snomed, 
        codelists.cancer_nonhaem_nonlung_snomed, 
        codelists.cancer_lung_snomed
        ),
    returning="date",
    on_or_after="index_date",
    find_first_match_in_period=True,
    date_format="YYYY-MM-DD",
    ),
  
)
