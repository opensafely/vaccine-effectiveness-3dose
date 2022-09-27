from cohortextractor import patients, combine_codelists
from codelists import *
import codelists

############################################################
## functions
from variables_functions import *
############################################################

def generate_post_variables(index_date):
  post_variables = dict(
  
    # deregistration date
    dereg_date=patients.date_deregistered_from_all_supported_practices(
      on_or_after=index_date,
      date_format="YYYY-MM-DD",
    ),
  
    # positive covid test
    **covid_test_date_X(
      name = "positive_test",
      index_date = index_date,
      n = 6,
      test_result="positive",
      ),
  
    # emergency attendance
    **emergency_attendance_date_X(
      name = "emergency",
      n = 6,
      index_date = index_date,
    ),
  
    # any emergency attendance for covid
    **emergency_attendance_date_X(
      name = "covidemergency",
      n = 4,
      index_date = index_date,
      with_these_diagnoses = codelists.covid_emergency
    ),
  
    **emergency_attendance_date_X(
      name = "emergencyhosp",
      n = 4,
      index_date = index_date,
      discharged_to = codelists.discharged_to_hospital
    ),
  
    **emergency_attendance_date_X(
      name = "covidemergencyhosp",
      n = 4,
      index_date = index_date,
      with_these_diagnoses = codelists.covid_emergency,
      discharged_to = codelists.discharged_to_hospital
    ),
    
    
    # unplanned hospital admission
    **admitted_date_X(
      name = "unplanned",
      n = 6,
      index_date = index_date,
      # see https://github.com/opensafely-core/cohort-extractor/pull/497 for codes
      # see https://docs.opensafely.org/study-def-variables/#sus for more info
      with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
      with_patient_classification = ["1"], # ordinary admissions only
    ),
  
    # planned hospital admission
    **admitted_date_X(
      name = "planned",
      n = 6,
      index_date = index_date,
      # see https://github.com/opensafely-core/cohort-extractor/pull/497 for codes
      # see https://docs.opensafely.org/study-def-variables/#sus for more info
      with_admission_method=["11", "12", "13", "81"],
      with_patient_classification = ["1"], # ordinary and day-case admissions only
    ),
  
    ## Covid-related unplanned ICU hospital admissions 
    # we only need first admission for covid-related hospitalisation post,
    # but to identify first ICU / critical care admission date, we need sequential admissions
    # this assumes that a spell that is subsequent and contiguous to a covid-related admission is also coded with a code in codelists.covid_icd10
  
    # Positive covid admission prior to study start date 
    **admitted_date_X(
      name = "covid",
      n = 4,
      index_date = index_date,
      # see https://github.com/opensafely-core/cohort-extractor/pull/497 for codes
      # see https://docs.opensafely.org/study-def-variables/#sus for more info
      with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
      with_these_diagnoses=codelists.covid_icd10
    ),
  
    ## Covid-related unplanned ICU hospital admissions -- number of days in critical care for each covid-related admission
    **admitted_daysincritcare_X(
      name = "covid",
      n = 4,
      index_name = "covid",
      index_date = index_date,
      # see https://github.com/opensafely-core/cohort-extractor/pull/497 for codes
      # see https://docs.opensafely.org/study-def-variables/#sus for more info
      with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
      with_these_diagnoses=codelists.covid_icd10,
      # not filtering on patient classification as we're interested in anyone who is "really sick due to COVID"
      # most likely these are ordinary admissions but we'd want to know about other (potentially misclassified) admissions too
    ),
  
  
    # first test after index date
    covid_test_1_date=patients.with_test_result_in_sgss(
      pathogen="SARS-CoV-2",
      test_result="any",
      on_or_after=index_date,
      find_first_match_in_period=True,
      restrict_to_earliest_specimen_date=False,
      returning="date",
      date_format="YYYY-MM-DD",
    ),

    # Covid-related death
    coviddeath_date=patients.with_these_codes_on_death_certificate(
      codelists.covid_icd10,
      returning="date_of_death",
      date_format="YYYY-MM-DD",
    ),
  
    # All-cause death
    death_date=patients.died_from_any_cause(
      returning="date_of_death",
      date_format="YYYY-MM-DD",
    ),

  )
  return post_variables
