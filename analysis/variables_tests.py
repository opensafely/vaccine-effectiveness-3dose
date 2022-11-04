from cohortextractor import patients, combine_codelists
from codelists import *
import codelists

############################################################
## functions
from variables_functions import *
############################################################

n=10

def generate_tests_variables(index_date):
  tests_variables = dict(

    # number of tests in each post-baseline period
    **covid_test_n_X("anytest", index_date, 1, 28, "any"),

    # dates of tests (to match to symoptomatic vars)
    **covid_test_date_X(
        name="anytest",
        date_name="anytest",
        index_date=index_date,
        n=n,
        test_result="any",
        restrict_to_earliest_specimen_date=False,
        returning="date",
    ),
    # whether tests were symptomatic
    **covid_test_date_X(
        name="anytest",
        date_name="anytest",
        index_date=index_date,
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
        on_or_after=index_date,
        test_result="positive",
        restrict_to_earliest_specimen_date=True,
        returning="date",
    ),
    # case-category of first positive test
    firstpostest_category=patients.with_test_result_in_sgss(
        pathogen="SARS-CoV-2",
        on_or_after=index_date,
        test_result="positive",
        restrict_to_earliest_specimen_date=True,
        returning="case_category",
        return_expectations = {
            "incidence" : 1,
            "category": {"ratios": {"": 0.3, "LFT_Only": 0.4, "PCR_Only": 0.2, "LFT_WithPCR": 0.1}},
             }
    ),
    
  )
  
  return tests_variables