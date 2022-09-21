from cohortextractor import patients
from codelists import *
import codelists


def generate_jcvi_variables(index_date):
    jcvi_variables = dict(
    cev_ever = patients.with_these_clinical_events(
      codelists.shield,
      returning="binary_flag",
      on_or_before = "index_date - 1 day",
      find_last_match_in_period = True,
    ),
    
    endoflife = patients.satisfying(
      """
      midazolam OR
      endoflife_coding
      """,
    
      midazolam = patients.with_these_medications(
        codelists.midazolam,
        returning="binary_flag",
        on_or_before = "index_date - 1 day",
      ),
      
      endoflife_coding = patients.with_these_clinical_events(
        codelists.eol,
        returning="binary_flag",
        on_or_before = "index_date - 1 day",
        find_last_match_in_period = True,
      ),
          
    ),
      
    housebound = patients.satisfying(
      """housebound_date
      AND NOT no_longer_housebound
      AND NOT moved_into_care_home
      """,
          
      housebound_date=patients.with_these_clinical_events( 
        codelists.housebound, 
        on_or_before="index_date - 1 day",
        find_last_match_in_period = True,
        returning="date",
        date_format="YYYY-MM-DD",
      ),   
      no_longer_housebound=patients.with_these_clinical_events( 
        codelists.no_longer_housebound, 
        between=["housebound_date", "index_date - 1 day"],
      ),
      moved_into_care_home=patients.with_these_clinical_events(
        codelists.carehome,
        between=["housebound_date", "index_date - 1 day"],
      ),
    ),
  
    prior_covid_test_frequency=patients.with_test_result_in_sgss(
      pathogen="SARS-CoV-2",
      test_result="any",
      between=["index_date - 182 days", "index_date - 1 day"], # 182 days = 26 weeks
      returning="number_of_matches_in_period", 
      date_format="YYYY-MM-DD",
      restrict_to_earliest_specimen_date=False,
    ),
  
  # # overnight hospital admission at time of 3rd / booster dose
  # inhospital = patients.satisfying(
  
  #   "discharged_0_date >= index_date",
    
  #   discharged_0_date=patients.admitted_to_hospital(
  #     returning="date_discharged",
  #     on_or_before="index_date", # this is the admission date
  #     # see https://github.com/opensafely-core/cohort-extractor/pull/497 for codes
  #     # see https://docs.opensafely.org/study-def-variables/#sus for more info
  #     with_admission_method = ['11', '12', '13', '21', '2A', '22', '23', '24', '25', '2D', '28', '2B', '81'],
  #     with_patient_classification = ["1"], # ordinary admissions only
  #     date_format="YYYY-MM-DD",
  #     find_last_match_in_period=True,
  #   ), 
  # ),
  
 
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
      on_or_before="index_date - 1 day",
    ),
    ### any immunosuppression medication codes is recorded
    immrx_1=patients.with_these_clinical_events(
      codelists.immrx_cod,
      on_or_before="index_date - 1 day",
    ),
    ### receiving chemotherapy or radiotherapy
    dxt_chemo = patients.with_these_clinical_events(
      codelists.dxt_chemo_cod ,
      on_or_before="index_date - 1 day",
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
      on_or_before="index_date - 1 day",
    ),
    ### chronic kidney disease codes - all stages
    ckd15 = patients.with_these_clinical_events(
      codelists.ckd15_cod,
      on_or_before="index_date - 1 day",
    ),
    ### date of chronic kidney disease codes-stages 3 – 5  
    ckd35_dat=patients.with_these_clinical_events(
      codelists.ckd35_cod,
    returning="date",
    date_format="YYYY-MM-DD",
    on_or_before="index_date - 1 day",
    find_first_match_in_period=True,
    ),
    ### date of chronic kidney disease codes - all stages
    ckd15_dat=patients.with_these_clinical_events(
      codelists.ckd15_cod,
    returning="date",
    date_format="YYYY-MM-DD",
    on_or_before="index_date - 1 day",
    find_first_match_in_period=True,
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
      astrxm2_1
      )
      """,
        ### asthma admission codes
        astadm_1 = patients.with_these_clinical_events(
            codelists.astadm_cod,
            on_or_before="index_date - 1 day",
        ),  
        ### asthma diagnosis code
        ast_1 = patients.with_these_clinical_events(
            codelists.ast_cod,
            on_or_before="index_date - 1 day",
        ),  
        ### asthma - inhalers in last 12 months
        astrxm1_1=patients.with_these_medications(
            codelists.astrx_cod,
            returning="binary_flag",
            between=["index_date - 30 days", "index_date - 1 day"],
          ),
        ### asthma - systemic oral steroid prescription codes in last 24 months
        astrxm2_1=patients.with_these_medications(
            codelists.astrx_cod,
            returning="binary_flag",
            between=["index_date - 60 days", "index_date - 31 days"],
          ),
      ),
    ### chronic respiratory disease
    resp_cov_1 =patients.with_these_clinical_events(
        codelists.resp_cov_cod,
        returning="binary_flag",
        on_or_before="index_date - 1 day",
    ),
    ),
    ### patients with diabetes
    diab_group = patients.satisfying(
    """
    diab_dat>dmres_dat
    OR
    addis
    OR
    gdiab_group
    """,
      ### date any diabetes diagnosis read code is recorded
      diab_dat=patients.with_these_clinical_events(
        codelists.diab_cod,
        returning="date",
        find_last_match_in_period=True,
        on_or_before="index_date - 1 day",
        date_format="YYYY-MM-DD",
      ),
      ### date of diabetes resolved codes
      dmres_dat=patients.with_these_clinical_events(
        codelists.dmres_cod,
        returning="date",
        find_last_match_in_period=True,
        on_or_before="index_date - 1 day",
        date_format="YYYY-MM-DD",
      ),
      ### addison’s disease & pan-hypopituitary diagnosis codes
      addis = patients.with_these_clinical_events(
        codelists.addis_cod,
        returning="binary_flag",
        on_or_before="index_date - 1 day",
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
          on_or_before="index_date - 1 day",
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
              ### pregnancy in the previous 44 weeks 
            between=["index_date - 308 days", "index_date - 1 days"],
            ),
            ### pregnancy or delivery codes 
            pregdel_dat=patients.with_these_clinical_events(
            codelists.pregdel_cod,
            returning="date",
            find_last_match_in_period=True,
            on_or_before="index_date - 1 day",
            date_format="YYYY-MM-DD",
            ),
            ### date of pregnancy codes recorded
            preg_dat=patients.with_these_clinical_events(
            codelists.preg_cod,
            returning="date",
            find_last_match_in_period=True,
            on_or_before="index_date - 1 day",
            date_format="YYYY-MM-DD",
            ),
        ),
      ),
    ),  
    ### chronic liver disease codes
    cld = patients.with_these_clinical_events(
      codelists.cld_cod,
      returning="binary_flag",
      on_or_before="index_date - 1 day",
    ),
    ### patients with cns disease (including stroke/tia)
    cns_group= patients.with_these_clinical_events(
      codelists.cns_cov_cod,
      returning="binary_flag",
      on_or_before="index_date - 1 day",
    ),  
    ### chronic heart disease codes
    chd_cov = patients.with_these_clinical_events(
      codelists.chd_cov_cod,
      returning="binary_flag",
      on_or_before="index_date - 1 day",
    ),  
    ### asplenia or dysfunction of the spleen codes
    spln_cov= patients.with_these_clinical_events(
        codelists.spln_cov_cod,
        returning="binary_flag",
        on_or_before="index_date - 1 day",
    ),  
    ### wider learning disability
    learndis_1 = patients.with_these_clinical_events(
        codelists.learndis_cod,
        returning="binary_flag",
        on_or_before="index_date - 1 day",
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
        on_or_before="index_date - 1 day",
        date_format="YYYY-MM-DD",
    ),
    ### date of remission codes relating to severe mental illness
    smhres_dat=   patients.with_these_clinical_events(
      codelists.smhres_cod,
      returning="date",
      find_last_match_in_period=True,
      on_or_before="index_date - 1 day",
      date_format="YYYY-MM-DD",
    ),
    ),
    ),
     cev = patients.satisfying(
    """severely_clinically_vulnerable AND NOT less_vulnerable""",
    ##### The shielded patient list was retired in March/April 2021 when shielding ended
    ##### so it might be worth using that as the end date instead of index_date, as we're not sure
    ##### what has happened to these codes since then, e.g. have doctors still been adding new
    ##### shielding flags or low-risk flags? Depends what you're looking for really. Could investigate separately.

    ### SHIELDED GROUP - first flag all patients with "high risk" codes
    severely_clinically_vulnerable = patients.with_these_clinical_events(
    codelists.shield,
    returning="binary_flag",
    on_or_before = "index_date - 1 day",
    find_last_match_in_period = True,
    ),

    # find date at which the high risk code was added
    date_severely_clinically_vulnerable = patients.date_of(
    "severely_clinically_vulnerable",
    date_format="YYYY-MM-DD",
    ),

    ### NOT SHIELDED GROUP (medium and low risk) - only flag if later than 'shielded'
    less_vulnerable = patients.with_these_clinical_events(
    codelists.nonshield,
    between=["date_severely_clinically_vulnerable + 1 day", "index_date - 1 day",],
    ),
    ),
    )
    return jcvi_variables