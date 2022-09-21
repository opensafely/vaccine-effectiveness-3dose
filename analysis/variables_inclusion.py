from cohortextractor import patients, combine_codelists
from codelists import *
import json
import codelists


def generate_inclusion_variables(index_date):
  inclusion_variables = dict(
    
    registered = patients.registered_as_of(
        index_date,
    ), 

    has_died = patients.died_from_any_cause(
      on_or_before=index_date,
      returning="binary_flag",
    ),

    age_aug21 = patients.age_as_of( 
        "2021-08-31",
    ),
    
    age = patients.age_as_of( 
        f"{index_date} - 1 day",
    ),
    
    wchild = patients.satisfying(
      """
      atrisk_group
      OR
      hhld_imdef
      """,
        hhld_imdef = patients.with_these_clinical_events(
            codelists.hhld_imdef_cod,
            on_or_before=f"{index_date} - 1 day",
        ),
        ########## At risk
        atrisk_group = patients.satisfying(
          """
          immunogroup
          OR
          ckd_group
          OR
          resp_group
          OR
          diab_group 
          OR
          cld
          OR
          cns_group
          OR
          chd_cov
          OR
          spln_cov
          OR
          learndis_1
          OR
          sevment_group
          OR
          preg1_group
          """,
          ##### patients with immunosuppression
          immunogroup = patients.satisfying(
          """
          immdx_1
          OR
          immrx_1
          OR
          dxt_chemo
          """,
          ###  any immunosuppressant read code is recorded
          immdx_1=patients.with_these_clinical_events(
            codelists.immdx_cov_cod,
            on_or_before=f"{index_date} - 1 day",
          ),
          ### any immunosuppression medication codes is recorded
          # The medication code should be in the last 6 months, however to 
          # prevent patients from dropping out of the audit as the vaccination 
          # campaign progresses, we look back for medications from 01/07/2020.  
          # This ensures that patients whose latest immunosuppressant medication 
          # issue was originally within the 6 month timescale but then subsequently 
          # exceeds it are still included in the at-risk group.
          immrx_1=patients.with_these_clinical_events(
            codelists.immrx_cod,
            between=["2020-07-01", f"{index_date} - 1 day"],
          ),
          ### receiving chemotherapy or radiotherapy
          dxt_chemo = patients.with_these_clinical_events(
            codelists.dxt_chemo_cod ,
            between=["2020-07-01", f"{index_date} - 1 day"],
          ),
          ),
          # patients with chronic kidney disease
          ckd_group= patients.satisfying(
          """
          ckd_cov
          OR
          (
          ckd15
          AND
          ckd35_dat>=ckd15_dat 
          )
          """,
          ### chronic kidney disease diagnostic codes
          ckd_cov = patients.with_these_clinical_events(
            codelists.ckd_cov_cod,
            on_or_before=f"{index_date} - 1 day",
          ),
          ### chronic kidney disease codes - all stages
          ckd15 = patients.with_these_clinical_events(
            codelists.ckd15_cod,
            on_or_before=f"{index_date} - 1 day",
          ),
          ### date of chronic kidney disease codes-stages 3 – 5  
          ckd35_dat=patients.with_these_clinical_events(
            codelists.ckd35_cod,
          returning="date",
          date_format="YYYY-MM-DD",
          on_or_before=f"{index_date} - 1 day",
          find_last_match_in_period=True,
          ),
          ### date of chronic kidney disease codes - all stages
          ckd15_dat=patients.with_these_clinical_events(
            codelists.ckd15_cod,
          returning="date",
          date_format="YYYY-MM-DD",
          on_or_before=f"{index_date} - 1 day",
          find_last_match_in_period=True,
          ),
          ),
          ### patients who have chronic respiratory disease
          resp_group = patients.satisfying(
          """
          ast_group
          OR
          resp_cov_1 
          """,
            ### patients with asthma 
            ast_group= patients.satisfying(
            """
            astadm_1
            OR
            (
            ast_1
            AND
            astrxm1_1
            AND
            astrxm2_1 >= 2
            )
            """,
              ### asthma admission codes
              astadm_1 = patients.with_these_clinical_events(
                  codelists.astadm_cod,
                  between=[f"{index_date} - 730 days",f"{index_date} - 1 day"],
              ),  
              ### asthma diagnosis code
              ast_1 = patients.with_these_clinical_events(
                  codelists.ast_cod,
                  on_or_before=f"{index_date} - 1 day",
              ),  
              ### asthma - inhalers in last 12 months
              astrxm1_1=patients.with_these_medications(
                  codelists.astrx_cod,
                  returning="binary_flag",
                  between=[f"{index_date} - 12 months", f"{index_date} - 1 day"],
                ),
              ### asthma - systemic oral steroid prescription codes in last 24 months
              astrxm2_1=patients.with_these_medications(
                  codelists.astrx_cod,
                  returning="number_of_matches_in_period",
                  between=[f"{index_date} - 24 months", f"{index_date} - 1 day"],
                ),
            ),
          ### chronic respiratory disease
          resp_cov_1 =patients.with_these_clinical_events(
              codelists.resp_cov_cod,
              returning="binary_flag",
              on_or_before=f"{index_date} - 1 day",
          ),
          ),
          ### patients with diabetes
          diab_group = patients.satisfying(
          """
          (
          diab
          OR
          diab_dat>dmres_dat
          )
          OR
          addis
          OR
          gdiab_group
          """,
            diab = patients.with_these_clinical_events(
              codelists.diab_cod,
              returning="binary_flag",
              on_or_before=f"{index_date} - 1 day",
            ),
            ### date any diabetes diagnosis read code is recorded
            diab_dat=patients.with_these_clinical_events(
              codelists.diab_cod,
              returning="date",
              find_last_match_in_period=True,
              on_or_before=f"{index_date} - 1 day",
              date_format="YYYY-MM-DD",
            ),
            ### date of diabetes resolved codes
            dmres_dat=patients.with_these_clinical_events(
              codelists.dmres_cod,
              returning="date",
              find_last_match_in_period=True,
              on_or_before=f"{index_date} - 1 day",
              date_format="YYYY-MM-DD",
            ),
            ### addison’s disease & pan-hypopituitary diagnosis codes
            addis = patients.with_these_clinical_events(
              codelists.addis_cod,
              returning="binary_flag",
              on_or_before=f"{index_date} - 1 day",
            ),
            ### patients who are currently pregnant with gestational diabetes 
            gdiab_group = patients.satisfying(
            """
            gdiab
            AND
            preg1_group
            """,
              ### gestational diabetes diagnosis codes
              gdiab =  patients.with_these_clinical_events(
                codelists.gdiab_cod,
                returning="binary_flag",
                between=[f"{index_date} - 254 days",f"{index_date} - 1 day"],
              ),
              ### patients who are currently pregnant 
              preg1_group = patients.satisfying(
                """
                preg
                AND
                pregdel_dat < preg_dat
                """,
                  ### pregnancy codes recorded
                  preg =  patients.with_these_clinical_events(
                  codelists.preg_cod,
                  returning="binary_flag",
                    ### Pregnancy codes recorded in the 8.5 months before the audit run date
                  between=[f"{index_date} - 254 days", f"{index_date} - 1 days"],
                  ),
                  ### pregnancy or delivery codes 
                  pregdel_dat=patients.with_these_clinical_events(
                  codelists.pregdel_cod,
                  returning="date",
                  find_last_match_in_period=True,
                  between=[f"{index_date} - 254 days", f"{index_date} - 1 days"],
                  date_format="YYYY-MM-DD",
                  ),
                  ### date of pregnancy codes recorded
                  # (Pregnancy or Delivery codes recorded in the 8.5 months before audit run date)
                  preg_dat=patients.with_these_clinical_events(
                  codelists.preg_cod,
                  returning="date",
                  find_last_match_in_period=True,
                  between=[f"{index_date} - 254 days", f"{index_date} - 1 days"],
                  date_format="YYYY-MM-DD",
                  ),
              ),
            ),
          ),  
          ### chronic liver disease codes
          cld = patients.with_these_clinical_events(
            codelists.cld_cod,
            returning="binary_flag",
            on_or_before=f"{index_date} - 1 day",
          ),
          ### patients with cns disease (including stroke/tia)
          cns_group= patients.with_these_clinical_events(
            codelists.cns_cov_cod,
            returning="binary_flag",
            on_or_before=f"{index_date} - 1 day",
          ),  
          ### chronic heart disease codes
          chd_cov = patients.with_these_clinical_events(
            codelists.chd_cov_cod,
            returning="binary_flag",
            on_or_before=f"{index_date} - 1 day",
          ),  
          ### asplenia or dysfunction of the spleen codes
          spln_cov= patients.with_these_clinical_events(
              codelists.spln_cov_cod,
              returning="binary_flag",
              on_or_before=f"{index_date} - 1 day",
          ),  
          ### wider learning disability
          learndis_1 = patients.with_these_clinical_events(
              codelists.learndis_cod,
              returning="binary_flag",
              on_or_before=f"{index_date} - 1 day",
          ),
          ### patients with severe mental health
          sevment_group = patients.satisfying(
          """
          sev_mental_dat > smhres_dat
          """,
          ### date of severe mental illness codes
          sev_mental_dat=   patients.with_these_clinical_events(
              codelists.sev_mental_cod,
              returning="date",
              find_last_match_in_period=True,
              on_or_before=f"{index_date} - 1 day",
              date_format="YYYY-MM-DD",
          ),
          ### date of remission codes relating to severe mental illness
          smhres_dat=   patients.with_these_clinical_events(
            codelists.smhres_cod,
            returning="date",
            find_last_match_in_period=True,
            on_or_before=f"{index_date} - 1 day",
            date_format="YYYY-MM-DD",
          ),
        ),
      ),
    ),
    ########## end of at risk criteria
    
          
  )
  return inclusion_variables

