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
inclusion_variables = generate_inclusion_variables(index_date="trial_date")
############################################################
## matching variables
from variables_matching import generate_matching_variables 
matching_variables = generate_matching_variables(index_date="trial_date")
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
  
  index_date = "2020-01-01", # this shouldn't be used anywhere!
  
  # This line defines the study population
  # FIXME this line needs to be matching_round specific -- currently it's only using data from matching_round=1
  # might be necessary to have round-specific study definitions which is a pain, but metaprogrammable.
  
  
  population = patients.satisfying(
    
    f"""
      registered
      AND
      age_aug21 >= {minage}
      AND
      age_aug21 <= {maxage}
      AND
      (NOT has_died)
      AND
      (NOT wchild)
      AND
      prematched
    """,
    
    prematched = patients.which_exist_in_file(f_path=f"output/{cohort}/matchround{matching_round}/potential/potential_matchedcontrols.csv.gz"),
    
  ),
  trial_date = patients.with_value_from_file(f_path=f"output/{cohort}/matchround{matching_round}/potential/potential_matchedcontrols.csv.gz", returning="trial_date", returning_type="date", date_format='YYYY-MM-DD'),
  
  match_id = patients.with_value_from_file(f_path=f"output/{cohort}/matchround{matching_round}/potential/potential_matchedcontrols.csv.gz", returning="match_id", returning_type="int"),
  
  **vaccination_date_X(
    name = "covid_vax_any",
    index_date = "1900-01-01",
    n = 1,
    target_disease_matches="SARS-2 CORONAVIRUS"
  ),
  **inclusion_variables,
  **matching_variables,
)
