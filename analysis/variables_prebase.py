from cohortextractor import patients, combine_codelists
from codelists import *
import codelists

def generate_prebase_variables(index_date):
    prebase_variables = dict(
      ################################################################################################
  ## Pre-baseline events where event date is of interest
  ################################################################################################


  # Positive case identification prior to study start date
  primary_care_covid_case_0_date=patients.with_these_clinical_events(
    combine_codelists(
      codelists.covid_primary_care_code,
      codelists.covid_primary_care_positive_test,
      codelists.covid_primary_care_sequelae,
    ),
    returning="date",
    date_format="YYYY-MM-DD",
    on_or_before="index_date - 1 day",
    find_last_match_in_period=True,
  ),
  
  # covid PCR test dates from SGSS
  covid_test_0_date=patients.with_test_result_in_sgss(
    pathogen="SARS-CoV-2",
    test_result="any",
    on_or_before="index_date - 1 day",
    returning="date",
    date_format="YYYY-MM-DD",
    find_last_match_in_period=True,
    restrict_to_earliest_specimen_date=False,
  ),

  
  # positive covid test
  postest_0_date=patients.with_test_result_in_sgss(
      pathogen="SARS-CoV-2",
      test_result="positive",
      returning="date",
      date_format="YYYY-MM-DD",
      on_or_before="index_date - 1 day",
      find_last_match_in_period=True,
      restrict_to_earliest_specimen_date=False,
  ),
  
  # emergency attendance for covid
  covidemergency_0_date=patients.attended_emergency_care(
    returning="date_arrived",
    on_or_before="index_date - 1 day",
    with_these_diagnoses = codelists.covid_emergency,
    date_format="YYYY-MM-DD",
    find_last_match_in_period=True,
  ),
  
    # Positive covid admission prior to study start date
  covidadmitted_0_date=patients.admitted_to_hospital(
    returning="date_admitted",
    with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
    with_these_diagnoses=codelists.covid_icd10,
    on_or_before="index_date - 1 day",
    date_format="YYYY-MM-DD",
    find_last_match_in_period=True,
  ),  

            )
    return prebase_variables