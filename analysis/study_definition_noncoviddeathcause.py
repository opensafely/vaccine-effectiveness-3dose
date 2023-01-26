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

# Specify study defeinition
study = StudyDefinition(
  
  # Configure the expectations framework
  default_expectations={
    "date": {"earliest": "2020-01-01", "latest": "today"},
    "rate": "uniform",
    "incidence": 0.2,
    "int": {"distribution": "normal", "mean": 1000, "stddev": 100},
    "float": {"distribution": "normal", "mean": 25, "stddev": 5},
  },
  
  # This line defines the study population 
  population = patients.which_exist_in_file(f_path=f"output/{cohort}/match/data_matched_unique.csv.gz"),

  # deaths with a cardiovascular icd10 code
  cvddeath_date=patients.with_these_codes_on_death_certificate(
    codelists.cvd_combined,
    returning="date_of_death",
    date_format="YYYY-MM-DD",
  ),

  # death with a cancer icd10 code
  cancerdeath_date=patients.with_these_codes_on_death_certificate(
    codelists.cancer,
    returning="date_of_death",
    date_format="YYYY-MM-DD",
  ),
  
)
