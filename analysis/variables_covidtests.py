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

baselinedays = int(fup_params["baselinedays"])
postbaselinedays = int(fup_params["postbaselinedays"])
prebaselineperiods = int(fup_params["prebaselineperiods"])
postbaselineperiods = int(fup_params["postbaselineperiods"])
n_any = int(fup_params["n_any"])
n_pos = int(fup_params["n_pos"])


def generate_covidtests_variables(index_date):
  covidtests_variables = dict(

    # number of tests
    ## number of tests in 3 28-day pre-baseline periods
    **covidtest_n_X(
      "anytestpre", 
      index_date, 
      shift=-prebaselineperiods*postbaselinedays, 
      n=prebaselineperiods, 
      length=postbaselinedays,
       test_result="any"
       ),
    **covidtest_n_X(
      "postestpre", 
      index_date, 
      shift=-prebaselineperiods*postbaselinedays,
      n=prebaselineperiods,
      length=postbaselinedays, 
      test_result="positive"
      ),

    ## number of tests in 14 days after index date
    anytestpost_0_n=patients.with_test_result_in_sgss(
        pathogen="SARS-CoV-2",
        test_result="any",
        between=[f"{index_date}", f"{index_date} + {baselinedays-1} days"],
        find_first_match_in_period=True,
        restrict_to_earliest_specimen_date=False,
        returning="number_of_matches_in_period",
        date_format="YYYY-MM-DD",
      ),
      postestpost_0_n=patients.with_test_result_in_sgss(
        pathogen="SARS-CoV-2",
        test_result="positive",
        between=[f"{index_date}", f"{index_date} + {baselinedays-1} days"],
        find_first_match_in_period=True,
        restrict_to_earliest_specimen_date=False,
        returning="number_of_matches_in_period",
        date_format="YYYY-MM-DD",
      ),

    ## number of tests in 6 28-day post-baseline periods
    **covidtest_n_X(
      "anytestpost", 
      index_date, 
      shift=baselinedays, 
      n=postbaselineperiods, 
      length=postbaselinedays, 
      test_result="any"
      ),
    **covidtest_n_X(
      "postestpost", 
      index_date, 
      shift=baselinedays, 
      n=postbaselineperiods, 
      length=postbaselinedays, 
      test_result="positive"
      ),

    # dates of tests
    ## dates of tests (to match to symoptomatic vars)
    **covidtest_returning_X(
        name="anytest",
        date_name="anytest",
        index_date=f"{index_date} - {prebaselineperiods*postbaselinedays} days",
        n=n_any,
        test_result="any",
        restrict_to_earliest_specimen_date=False,
        returning="date",
    ),
    ## whether tests were symptomatic
    **covidtest_returning_X(
        name="anytest",
        date_name="anytest",
        index_date=f"{index_date} - {prebaselineperiods*postbaselinedays} days",
        n=n_any,
        test_result="any",
        restrict_to_earliest_specimen_date=False,
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
        date_name="postest",
        index_date=f"{index_date} - {prebaselineperiods*postbaselinedays} days",
        n=n_pos,
        test_result="any",
        restrict_to_earliest_specimen_date=False,
        returning="date",
    ),

    # date of first positive test (to match to case category vars)
    firstpostest_date=patients.with_test_result_in_sgss(
        pathogen="SARS-CoV-2",
        on_or_after=f"{index_date} - {prebaselineperiods*postbaselinedays} days",
        test_result="positive",
        restrict_to_earliest_specimen_date=True,
        returning="date",
    ),
    # case-category of first positive test
    firstpostest_category=patients.with_test_result_in_sgss(
        pathogen="SARS-CoV-2",
        on_or_after=f"{index_date} - {prebaselineperiods*postbaselinedays} days",
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