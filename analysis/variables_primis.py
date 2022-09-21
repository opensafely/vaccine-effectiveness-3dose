from cohortextractor import patients, combine_codelists
from codelists import *
import json
import codelists

def generate_primis_variables(index_date,prefix):
    primis_variables = dict(
     # From PRIMIS

  
    # asthma = patients.satisfying(
    # """
    #   astadm OR
    #   (ast AND astrxm1 AND astrxm2 AND astrxm3)
    #   """,
    # # Asthma Admission codes
    # # astadm=patients.with_these_clinical_events(
    # #   codelists.astadm,
    # #   returning="binary_flag",
    # #   on_or_before="index_date - 1 day",
    # # ),
    # # # Asthma Diagnosis code
    # # ast = patients.with_these_clinical_events(
    # #   codelists.ast,
    # #   returning="binary_flag",
    # #   on_or_before="index_date - 1 day",
    # # ),
    # # # Asthma systemic steroid prescription code in month 1
    # # astrxm1=patients.with_these_medications(
    # #   codelists.astrx,
    # #   returning="binary_flag",
    # #   between=["index_date - 30 days", "index_date - 1 day"],
    # # ),
    # # # Asthma systemic steroid prescription code in month 2
    # # astrxm2=patients.with_these_medications(
    # #   codelists.astrx,
    # #   returning="binary_flag",
    # #   between=["index_date - 60 days", "index_date - 31 days"],
    # # ),
    # # Asthma systemic steroid prescription code in month 3
    # astrxm3=patients.with_these_medications(
    #   codelists.astrx,
    #   returning="binary_flag",
    #   between= ["index_date - 90 days", "index_date - 61 days"],
    # ),

    # ),

    # Chronic Neurological Disease including Significant Learning Disorder
    chronic_neuro_disease=patients.with_these_clinical_events(
    codelists.cns_cov,
    returning="binary_flag",
    on_or_before="index_date - 1 day",
    ),

    # # Chronic Respiratory Disease
    # chronic_resp_disease = patients.satisfying(
    # "asthma OR resp_cov",
    # resp_cov=patients.with_these_clinical_events(
    #   codelists.resp_cov,
    #   returning="binary_flag",
    #   on_or_before="index_date - 1 day",
    # ),
    # ),

    sev_obesity = patients.satisfying(
    """
      sev_obesity_date > bmi_date OR
      bmi_value1 >= 40
      """,

    bmi_stage_date=patients.with_these_clinical_events(
      codelists.bmi_stage,
      returning="date",
      find_last_match_in_period=True,
      on_or_before="index_date - 1 day",
      date_format="YYYY-MM-DD",
    ),

    sev_obesity_date=patients.with_these_clinical_events(
      codelists.sev_obesity,
      returning="date",
      find_last_match_in_period=True,
      ignore_missing_values=True,
      between= ["bmi_stage_date", "index_date - 1 day"],
      date_format="YYYY-MM-DD",
    ),

    bmi_date=patients.with_these_clinical_events(
      codelists.bmi,
      returning="date",
      ignore_missing_values=True,
      find_last_match_in_period=True,
      on_or_before="index_date - 1 day",
      date_format="YYYY-MM-DD",
    ),

    bmi_value1=patients.with_these_clinical_events(
      codelists.bmi,
      returning="numeric_value",
      ignore_missing_values=True,
      find_last_match_in_period=True,
      on_or_before="index_date - 1 day",
    ),

    ),

    diabetes = patients.satisfying(
    "(dmres_date < diab_date) OR (diab_date AND (NOT dmres_date))",
    
    diab_date=patients.with_these_clinical_events(
      codelists.diab,
      returning="date",
      find_last_match_in_period=True,
      on_or_before="index_date - 1 day",
      date_format="YYYY-MM-DD",
    ),

    dmres_date=patients.with_these_clinical_events(
      codelists.dmres,
      returning="date",
      find_last_match_in_period=True,
      on_or_before="index_date - 1 day",
      date_format="YYYY-MM-DD",
    ),
    ),

    sev_mental=patients.satisfying(
    "(smhres_date < sev_mental_date) OR (sev_mental_date AND (NOT smhres_date))",

    # Severe Mental Illness codes
    sev_mental_date=patients.with_these_clinical_events(
      codelists.sev_mental,
      returning="date",
      find_last_match_in_period=True,
      on_or_before="index_date - 1 day",
      date_format="YYYY-MM-DD",
    ),
    # Remission codes relating to Severe Mental Illness
    smhres_date=patients.with_these_clinical_events(
      codelists.smhres,
      returning="date",
      find_last_match_in_period=True,
      on_or_before="index_date - 1 day",
      date_format="YYYY-MM-DD",
    ),
    ),


    # Chronic heart disease codes
    chronic_heart_disease=patients.with_these_clinical_events(
      codelists.chd_cov,
      returning="binary_flag",
      on_or_before="index_date - 1 day",
    ),

    chronic_kidney_disease=patients.satisfying(
    """
      ckd OR
      (ckd15_date AND ckd35_date >= ckd15_date)
      """,

    # Chronic kidney disease codes - all stages
    ckd15_date=patients.with_these_clinical_events(
      codelists.ckd15,
      returning="date",
      find_last_match_in_period=True,
      on_or_before="index_date - 1 day",
      date_format="YYYY-MM-DD",
    ),

    # Chronic kidney disease codes-stages 3 - 5
    ckd35_date=patients.with_these_clinical_events(
      codelists.ckd35,
      returning="date",
      find_last_match_in_period=True,
      on_or_before="index_date - 1 day",
      date_format="YYYY-MM-DD",
    ),

    # Chronic kidney disease diagnostic codes
    ckd=patients.with_these_clinical_events(
      codelists.ckd_cov,
      returning="binary_flag",
      on_or_before="index_date - 1 day",
    ),
    ),


    # Chronic Liver disease codes
    chronic_liver_disease=patients.with_these_clinical_events(
    codelists.cld,
    returning="binary_flag",
    on_or_before="index_date - 1 day",
    ),


    # immunosuppressed=patients.satisfying(
    # "immrx OR immdx",

    # # Immunosuppression diagnosis codes
    # immdx=patients.with_these_clinical_events(
    #   codelists.immdx_cov,
    #   returning="binary_flag",
    #   on_or_before="index_date - 1 day",
    # ),
    # # Immunosuppression medication codes
    # immrx=patients.with_these_medications(
    #   codelists.immrx,
    #   returning="binary_flag",
    #   between=["index_date - 182 days", "index_date - 1 day"]
    # ),
    # ),

    # Asplenia or Dysfunction of the Spleen codes
    asplenia=patients.with_these_clinical_events(
    codelists.spln_cov,
    returning="binary_flag",
    on_or_before="index_date - 1 day",
    ),

    # Wider Learning Disability
    # learndis=patients.with_these_clinical_events(
    # codelists.learndis,
    # returning="binary_flag",
    # on_or_before="index_date - 1 day",
    # ),


    # to represent household contact of shielding individual
    # hhld_imdef_dat=patients.with_these_clinical_events(
    #   codelists.hhld_imdef,
    #   returning="date",
    #   find_last_match_in_period=True,
    #   on_or_before="index_date - 1 day",
    #   date_format="YYYY-MM-DD",
    # ),
    #
  

    cancer = patients.satisfying(
    
    "cancer_primary_care",
    # cancer_hosp=patients.admitted_to_hospital(
    #   with_these_diagnoses=combine_codelists(
    #     codelists.cancer_nonhaem_icd10,
    #     codelists.cancer_haem_icd10,
    #     codelists.cancer_unspec_icd10,
    #   ),
    #   between=["index_date - 3 years", "index_date - 1 day"],
    #   returning="binary_flag",
    # ),
    cancer_primary_care=patients.with_these_clinical_events( 
      combine_codelists(
        codelists.cancer_nonhaem_snomed,
        codelists.cancer_haem_snomed
      ),
      between=["index_date - 3 years", "index_date - 1 day"],
      returning="binary_flag",
    ), 
    ),
    )
    primis_variables = {f"{prefix}{key}": val for key, val in primis_variables.items()}
    return primis_variables

