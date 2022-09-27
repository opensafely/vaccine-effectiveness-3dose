from ast import And
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

from variables_outcome import vaccination_date_X 

cohort = params["cohort"]
matching_round = params["matching_round"]
previousmatching_round = int(matching_round)-1
index_date = params["index_date"]

# import study dates defined in "./analysis/design.R" script
with open("./lib/design/study-dates.json") as f:
  study_dates = json.load(f)

# change these in design.R if necessary
start_date = study_dates[cohort]["start_date"]
end_date = study_dates[cohort]["end_date"]


# import study parameters defined in "./analysis/design.R" script  
with open("./lib/design/study-params.json") as f:
  study_params = json.load(f)

minage = study_params[cohort]["minage"]
maxage = study_params[cohort]["maxage"]
treatment = study_params[cohort]["treatment"]



############################################################
## inclusion variables
from variables_inclusion import generate_inclusion_variables 
inclusion_variables = generate_inclusion_variables(index_date="index_date")
############################################################
## matching variables
from variables_matching import generate_matching_variables 
matching_variables = generate_matching_variables(index_date="index_date")
############################################################


# Specify study defeinition
study = StudyDefinition(
  
  # Configure the expectations framework
  default_expectations={
    "date": {"earliest": "2020-01-01", "latest": end_date},
    "rate": "uniform",
    "incidence": 0.2,
    "int": {"distribution": "normal", "mean": 1000, "stddev": 100},
    "float": {"distribution": "normal", "mean": 25, "stddev": 5},
  },
  
  index_date = index_date,
  
  # This line defines the study population
  population=patients.satisfying(
    f"""
      registered
      AND
      age_aug21 >= {minage}
      AND
      age_aug21 <= {maxage}
      AND
      (NOT has_died)
      AND
      (NOT child_atrisk)
    """,
    #NOT (covid_vax_any_1_date <= index_date) # doesn't work for some reason `unknown colunm : index_date`
    #previouslymatched = patients.which_exist_in_file(f_path="output/match/cumulative_matchedcontrols{matching_round}.csv.gz"),
  ),
  
  **vaccination_date_X(
    name = "covid_vax_any",
    index_date = "1900-01-01",
    n = 1,
    target_disease_matches="SARS-2 CORONAVIRUS"
  ),
  **inclusion_variables,    
  **matching_variables,      
)
