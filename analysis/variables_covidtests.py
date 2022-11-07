from cohortextractor import patients, combine_codelists
from codelists import *
import codelists

############################################################
## functions
from variables_functions import *
############################################################

n=10 # may to increase
prebaselineshift=3*28

def generate_covidtests_variables(index_date):
  covidtests_variables = dict(

    # number of tests
    ## number of tests in 3 28-day pre-baseline periods
    **covidtest_n_X("anytestpre", index_date, shift=prebaselineshift, n=3, length=28, test_result="any"),

    ## number of tests in 14 days after index date
    anytestpost_0_n=patients.with_test_result_in_sgss(
        pathogen="SARS-CoV-2",
        test_result="any",
        between=[f"{index_date}", f"{index_date} + 13 days"],
        find_first_match_in_period=True,
        restrict_to_earliest_specimen_date=False,
        returning="number_of_matches_in_period",
        date_format="YYYY-MM-DD",
      ),

    ## number of tests in 6 28-day post-baseline periods
    **covidtest_n_X("anytestpost", index_date, shift=14, n=6, length=28, test_result="any"),

    # dates of tests
    ## dates of tests (to match to symoptomatic vars)
    **covidtest_date_X(
        name="anytest",
        date_name="anytest",
        index_date=f"{index_date} - {prebaselineshift} days",
        n=n,
        test_result="any",
        restrict_to_earliest_specimen_date=False,
        returning="date",
    ),
    ## whether tests were symptomatic
    **covidtest_date_X(
        name="anytest",
        date_name="anytest",
        index_date=f"{index_date} - {prebaselineshift} days",
        n=n,
        test_result="any",
        restrict_to_earliest_specimen_date=False,
        returning="symptomatic",
        return_expectations = {
            "incidence" : 1,
            "category": {"ratios": {"": 0.5, "Y": 0.3, "N": 0.2}},
             }
    ),

    # date of first positive test (to match to case category vars)
    firstpostest_date=patients.with_test_result_in_sgss(
        pathogen="SARS-CoV-2",
        on_or_after=f"{index_date} - {prebaselineshift} days",
        test_result="positive",
        restrict_to_earliest_specimen_date=True,
        returning="date",
    ),
    # case-category of first positive test
    firstpostest_category=patients.with_test_result_in_sgss(
        pathogen="SARS-CoV-2",
        on_or_after=f"{index_date} - {prebaselineshift} days",
        test_result="positive",
        restrict_to_earliest_specimen_date=True,
        returning="case_category",
        return_expectations = {
            "incidence" : 1,
            "category": {"ratios": {"": 0.3, "LFT_Only": 0.4, "PCR_Only": 0.2, "LFT_WithPCR": 0.1}},
             }
    ),
    
  )
  
  return covidtests_variables