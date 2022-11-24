from cohortextractor import patients, combine_codelists
from codelists import *
import codelists

############################################################
## functions
from variables_functions import *
############################################################

# import json module
import json
# import study dates defined in "./analysis/design.R" script
with open("./lib/design/fup-params.json") as f:
  fup_params = json.load(f)

covidtestcuts = fup_params["covidtestcuts"]
n_any = int(fup_params["n_any"])
n_pos = int(fup_params["n_pos"])

def generate_covidtests_variables(index_date):
  covidtests_variables = dict(

    # number of tests
    ## number of tests in each of the covidtestcuts periods (closed on the right)
    **covidtest_n_X(
      "anytest", 
      index_date, 
      cuts=covidtestcuts,
      test_result="any"
       ),
    **covidtest_n_X(
      "postest", 
      index_date, 
      cuts=covidtestcuts,
      test_result="positive"
      ),

    # dates of tests
    ## dates of tests (to match to symoptomatic vars)
    **covidtest_returning_X(
        name="anytest",
        index_date=index_date,
        shift=int(covidtestcuts[0]),
        n=n_any,
        test_result="any",
        returning="date",
    ),
    ## whether tests were symptomatic
    **covidtest_returning_X(
        name="anytest",
        index_date=index_date,
        shift=int(covidtestcuts[0]),
        n=n_any,
        test_result="any",
        returning="symptomatic",
        return_expectations = {
            "incidence" : 1,
            # not using study def dummy data, but returns error without stating expectations
            "category": {"ratios": {"": 0.5, "Y": 0.3, "N": 0.2}},
             }
    ),
    # dates of positive tests
    **covidtest_returning_X(
        name="postest",
        index_date=index_date,
        shift=int(covidtestcuts[0]),
        n=n_pos,
        test_result="positive",
        returning="date",
    ),

    # date of first positive test (to match to case category vars, after index date only)
    firstpostest_date=patients.with_test_result_in_sgss(
        pathogen="SARS-CoV-2",
        on_or_after=index_date,
        test_result="positive",
        restrict_to_earliest_specimen_date=True,
        returning="date",
    ),
    # case-category of first positive test (after index date only)
    firstpostest_category=patients.with_test_result_in_sgss(
        pathogen="SARS-CoV-2",
        on_or_after=index_date,
        test_result="positive",
        restrict_to_earliest_specimen_date=True,
        returning="case_category",
        return_expectations = {
            "incidence" : 1,
            # not using study def dummy data, but returns error without stating expectations
            "category": {"ratios": {"": 0.3, "LFT_Only": 0.4, "PCR_Only": 0.2, "LFT_WithPCR": 0.1}},
             }
    ),
    
  )
  
  return covidtests_variables