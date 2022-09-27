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

# import study dates defined in "./analysis/design.R" script
with open("./lib/design/study-dates.json") as f:
  study_dates = json.load(f)

# change these in design.R if necessary
firstdose3_date = study_dates["pfizer"]["start_date"]
firstpossiblevax_date = study_dates["firstpossiblevax_date"]
studyend_date = study_dates["studyend_date"]
firstpfizer_date = study_dates["firstpfizer_date"]
firstaz_date = study_dates["firstaz_date"]
firstmoderna_date = study_dates["firstmoderna_date"]

############################################################
## inclusion variables
from variables_vax import generate_vax_variables 
vax_variables = generate_vax_variables(index_date="1900-01-01", n=4)
############################################################
# vax variables
from variables_inclusion import generate_inclusion_variables 
inclusion_variables = generate_inclusion_variables(index_date="covid_vax_disease_3_date")
############################################################
## jcvi variables
from variables_jcvi import generate_jcvi_variables 
jcvi_variables = generate_jcvi_variables(index_date="covid_vax_disease_3_date")
############################################################
## demographic variables
from variables_demographic import generate_demographic_variables 
demographic_variables = generate_demographic_variables(index_date="covid_vax_disease_3_date")
############################################################
## pre variables
from variables_pre import generate_pre_variables 
pre_variables = generate_pre_variables(index_date="covid_vax_disease_3_date")
############################################################
## post variables
from variables_post import generate_post_variables 
post_variables = generate_post_variables(index_date="covid_vax_disease_3_date")
############################################################

# Specify study definition
study = StudyDefinition(
  
  # Configure the expectations framework
  default_expectations={
    "date": {"earliest": "2020-01-01", "latest": studyend_date},
    "rate": "uniform",
    "incidence": 0.2,
    "int": {"distribution": "normal", "mean": 1000, "stddev": 100},
    "float": {"distribution": "normal", "mean": 25, "stddev": 5},
  },
  
  # This line defines the study population
  population=patients.satisfying(
    """
    registered
    AND
    age >= 18
    AND
    NOT has_died
    AND 
    covid_vax_disease_2_date
    """,
    
    **inclusion_variables,    

  ),
  
  #################################################################
  ## Covid vaccine dates
  #################################################################
  **vax_variables,
    
  ###############################################################################
  # jcvi variables
  ##############################################################################
  **jcvi_variables, 
  
  ###############################################################################
  # demographic variables
  ##############################################################################
  **demographic_variables,   

  ###############################################################################
  # pre variables
  ##############################################################################
  **pre_variables,      
  
  ###############################################################################
  # posts
  ##############################################################################
  **post_variables,      
  
)
