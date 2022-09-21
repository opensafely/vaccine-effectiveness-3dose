from cohortextractor import patients
from codelists import *
import codelists

def generate_demo_variables(index_date):
  demo_variables = dict(

    has_follow_up_previous_6weeks=patients.registered_with_one_practice_between(
      start_date="index_date - 42 days",
      end_date="index_date",
    ),
    
    age=patients.age_as_of( 
      "index_date - 1 day",
    ),
    
    ageband=patients.categorised_as(
              {
                  "0": "DEFAULT",
                  "5-11": """ age >= 5 AND age < 12""",
                  "12-15": """ age >= 12 AND age < 16"""
              },
              return_expectations={
                  "rate": "universal",
                  "category": {
                      "ratios": {
                          "5-11": 0.65,
                          "12-15": 0.35,
                      }
                  },
              },
      ),
    
    sex=patients.sex(
      return_expectations={
        "rate": "universal",
        "category": {"ratios": {"M": 0.49, "F": 0.51}},
        "incidence": 1,
      }
    ),
    
    # https://github.com/opensafely/risk-factors-research/issues/51
    bmi=patients.categorised_as(
      {
        "Not obese": "DEFAULT",
        "Obese I (30-34.9)": """ bmi_value >= 30 AND bmi_value < 35""",
        "Obese II (35-39.9)": """ bmi_value >= 35 AND bmi_value < 40""",
        "Obese III (40+)": """ bmi_value >= 40 AND bmi_value < 100""",
        # set maximum to avoid any impossibly extreme values being classified as obese
      },
      bmi_value=patients.most_recent_bmi(
        on_or_after="index_date - 5 years",
        minimum_age_at_measurement=16
      ),
      return_expectations={
        "rate": "universal",
        "category": {
          "ratios": {
            "Not obese": 0.7,
            "Obese I (30-34.9)": 0.1,
            "Obese II (35-39.9)": 0.1,
            "Obese III (40+)": 0.1,
          }
        },
      },
    ),
    
  
    ################################################################################################
    ## Practice and patient ID variables
    ################################################################################################
    # practice pseudo id
  
    # practice_id=patients.registered_practice_as_of(
    #   "index_date - 1 day",
    #   returning="pseudo_id",
    #   return_expectations={
    #     "int": {"distribution": "normal", "mean": 1000, "stddev": 100},
    #     "incidence": 1,
    #   },
    # ),
    
    # msoa
    
    msoa=patients.address_as_of(
      "index_date - 1 day",
      returning="msoa",
      return_expectations={
        "rate": "universal",
        "category": {"ratios": {"E02000001": 0.0625, "E02000002": 0.0625, "E02000003": 0.0625, "E02000004": 0.0625,
          "E02000005": 0.0625, "E02000007": 0.0625, "E02000008": 0.0625, "E02000009": 0.0625, 
          "E02000010": 0.0625, "E02000011": 0.0625, "E02000012": 0.0625, "E02000013": 0.0625, 
          "E02000014": 0.0625, "E02000015": 0.0625, "E02000016": 0.0625, "E02000017": 0.0625}},
      },
    ),    
    # stp is an NHS administration region based on geography
  
    stp=patients.registered_practice_as_of(
      "index_date - 1 day",
      returning="stp_code",
      return_expectations={
        "rate": "universal",
        "category": {
          "ratios": {
            "STP1": 0.1,
            "STP2": 0.1,
            "STP3": 0.1,
            "STP4": 0.1,
            "STP5": 0.1,
            "STP6": 0.1,
            "STP7": 0.1,
            "STP8": 0.1,
            "STP9": 0.1,
            "STP10": 0.1,
          }
        },
      },
    ),
    # NHS administrative region
  
    region=patients.registered_practice_as_of(
      "index_date - 1 day",
      returning="nuts1_region_name",
      return_expectations={
        "rate": "universal",
        "category": {
          "ratios": {
            "North East": 0.1,
            "North West": 0.1,
            "Yorkshire and The Humber": 0.2,
            "East Midlands": 0.1,
            "West Midlands": 0.1,
            "East": 0.1,
            "London": 0.1,
            "South East": 0.1,
            "South West": 0.1
            #"" : 0.01
          },
        },
      },
    ),
    ## IMD - quintile
    
    imd_Q5=patients.categorised_as(
      {
        "Unknown": "DEFAULT",
        "1 (most deprived)": "imd >= 0 AND imd < 32844*1/5",
        "2": "imd >= 32844*1/5 AND imd < 32844*2/5",
        "3": "imd >= 32844*2/5 AND imd < 32844*3/5",
        "4": "imd >= 32844*3/5 AND imd < 32844*4/5",
        "5 (least deprived)": "imd >= 32844*4/5 AND imd <= 32844",
      },
      return_expectations={
        "rate": "universal",
        "category": {"ratios": {"Unknown": 0.02, "1 (most deprived)": 0.18, "2": 0.2, "3": 0.2, "4": 0.2, "5 (least deprived)": 0.2}},
      },
      imd=patients.address_as_of(
      "index_date - 1 day",
      returning="index_of_multiple_deprivation",
      round_to_nearest=100,
      return_expectations={
        "category": {"ratios": {c: 1/320 for c in range(100, 32100, 100)}}
      }
      ),
  
    ),
    
    #rurality
    
    rural_urban=patients.address_as_of(
      "index_date - 1 day",
      returning="rural_urban_classification",
      return_expectations={
        "rate": "universal",
        "category": {"ratios": {1: 0.125, 2: 0.125, 3: 0.125, 4: 0.125, 5: 0.125, 6: 0.125, 7: 0.125, 8: 0.125}},
      },
    ),
  
  
    ################################################################################################
    ## occupation / residency
    ################################################################################################
  
  
    # health or social care worker  
    hscworker = patients.with_healthcare_worker_flag_on_covid_vaccine_record(returning="binary_flag"),
      
    care_home_tpp=patients.satisfying(
      "care_home='1'",
      
      care_home = patients.care_home_status_as_of(
        "index_date - 1 day",
        categorised_as={
            "1": "IsPotentialCareHome",
            "": "DEFAULT",  # use empty string
        },
      ),
    ),
    
    # Patients in long-stay nursing and residential care
    care_home_code=patients.with_these_clinical_events(
        codelists.carehome,
        on_or_before="index_date - 1 day",
        returning="binary_flag",
        return_expectations={"incidence": 0.01},
    ),

  )
  return demo_variables
