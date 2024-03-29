################################################################################
# functions for processing each of the variable groups in analysis/process_data.R

################################################################################
process_jcvi <- function(.data) {
  .data %>%
    mutate(
      
      multimorb =
        (sev_obesity) +
        (chronic_heart_disease) +
        (chronic_kidney_disease) +
        (diabetes) +
        (chronic_liver_disease) +
        (chronic_resp_disease | asthma) +
        (chronic_neuro_disease),
      multimorb = cut(multimorb, breaks = c(0, 1, 2, Inf), labels=c("0", "1", "2+"), right=FALSE),
      immunosuppressed = immunosuppressed | asplenia,
      
      # any carehome flag
      care_home_combined = care_home_tpp | care_home_code, 
      
      # clinically at-risk group
      cv = immunosuppressed | chronic_kidney_disease | chronic_resp_disease | diabetes | chronic_liver_disease |
        chronic_neuro_disease | chronic_heart_disease | asplenia | learndis | sev_mental,
      
      cev_cv = fct_case_when(
        cev ~ "Clinically extremely vulnerable",
        cv ~ "Clinically at-risk",
        TRUE ~ "Not clinically at-risk"
      ) %>% fct_rev(),
      
      # original priority groups https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/1007737/Greenbook_chapter_14a_30July2021.pdf#page=15
      # new priority groups https://www.england.nhs.uk/coronavirus/wp-content/uploads/sites/52/2021/07/C1327-covid-19-vaccination-autumn-winter-phase-3-planning.pdf
      # group 10 split into 16-39 and 40-49 because of earlier roll-out in 40+ from 15 Nov https://www.gov.uk/government/news/jcvi-issues-advice-on-covid-19-booster-vaccines-for-those-aged-40-to-49-and-second-doses-for-16-to-17-year-olds
      
      jcvi_ageband = cut(
        age_aug2021,
        breaks=c(-Inf, 18, 40, 50, 55, 60, 65, 70, 75, 80, Inf),
        labels=c("under 18", "18-39", "40-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80+"),
        right=FALSE
      ),
      
      
      jcvi_group = fct_case_when(
        care_home_combined | hscworker  ~ "01",
        age_aug2021>=80 ~ "02",
        age_aug2021>=75 ~ "03",
        age_aug2021>=70 ~ "04a",
        age_aug2021>=16 & cev ~ "04b",
        age_aug2021>=65 ~ "05",
        age_aug2021>=16 & cv ~ "06",
        age_aug2021>=60 ~ "07",
        age_aug2021>=55 ~ "08",
        age_aug2021>=50 ~ "09",
        age_aug2021>=40 ~ "10",
        age_aug2021>=30 ~ "11",
        TRUE ~ "12"
      ),
      
      jcvi_group_descr = fct_recode(
        jcvi_group,
        "Care home residents and health and social care workers"="01",
        "80+ years"="02",
        "75-79 years"="03",
        "70-74 years"="04a",
        "16-69 years and clinically extremely vulnerable"="04b",
        "65-69 years"="05",
        "16-64 years and clinically at-risk"="06",
        "60-64 years"="07",
        "55-59 years"="08",
        "50-54 years"="09",
        "40-49 years"="10",
        "30-39 years"="11",
        "18-29 years"="12"
      ),
      
    ) %>%
    select(-care_home_type, -care_home_tpp, -care_home_code)
}

################################################################################
process_covs <- function(.data) {
  .data %>%
    mutate(
      
      bmi = factor(bmi, levels = c("Not obese", "Obese I (30-34.9)", "Obese II (35-39.9)", "Obese III (40+)")),
      
      pregnancy = pregnancy & (sex == "Female") & (age < 50),
      
      prior_test_cat = cut(
        prior_test_frequency, 
        breaks=c(0, 1, 2, 3, Inf), 
        labels=c("0", "1", "2", "3+"), 
        right=FALSE
        )
      
    )  
}



################################################################################
process_demo <- function(.data) {
  .data %>%
    mutate(
      
      age65plus=age>=65,
      
      agegroup = cut(
        age, 
        breaks=c(-Inf, 18, 50, 65, 80, Inf),
        labels=c("under 18", "18-49", "50-64", "65-79", "80+"),
        right=FALSE
      ),
      
      ageband = cut(
        age,
        breaks=c(-Inf, 18, 40, 50, 60, 70, 80, 90, Inf),
        labels=c("under 18", "18-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90+"),
        right=FALSE
      ),
      
      sex = fct_case_when(
        sex == "F" ~ "Female",
        sex == "M" ~ "Male",
        #sex == "I" ~ "Inter-sex",
        #sex == "U" ~ "Unknown",
        TRUE ~ NA_character_
      ),
      
      ethnicity = if_else(is.na(ethnicity), ethnicity_6_sus, ethnicity),
      
      ethnicity = fct_case_when(
        ethnicity == "1" ~ "White",
        ethnicity == "4" ~ "Black",
        ethnicity == "3" ~ "South Asian",
        ethnicity == "2" ~ "Mixed",
        ethnicity == "5" ~ "Other",
        TRUE ~ NA_character_
        
      ),
      
      region = fct_collapse(
        region,
        `East of England` = "East",
        `London` = "London",
        `Midlands` = c("West Midlands", "East Midlands"),
        `North East and Yorkshire` = c("Yorkshire and The Humber", "North East"),
        `North West` = "North West",
        `South East` = "South East",
        `South West` = "South West"
      ),
      
      imd_Q5 = factor(imd_Q5, levels = c("1 (most deprived)", "2", "3", "4", "5 (least deprived)", "Unknown"))
      
    ) %>%
    select(-ethnicity_6_sus)
}

################################################################################
process_pre <- function(.data) {
  
  .data %>%
    mutate(
      # any covid event before study start
      prior_covid_infection = (!is.na(positive_test_0_date)) | (!is.na(admitted_covid_0_date)) | (!is.na(covidemergency_0_date)) | (!is.na(primary_care_covid_case_0_date)),
      # date of latest covid event before study start
      prior_covid_infection_date = pmax(positive_test_0_date, covidemergency_0_date, admitted_covid_0_date, na.rm=TRUE),
      time_since_infection = fct_case_when(
        is.na(prior_covid_infection_date) ~ "never",
        as.integer(index_date - prior_covid_infection_date) <= 30 ~ "1-30 days",
        as.integer(index_date - prior_covid_infection_date) <= 90 ~ "31-90 days",
        TRUE ~ "91+ days"
      )
      # note the slight discrepancy between definitions of `prior_covid_infection` (matching variable) and `anycovid_0_date` (used in exclusion criteria):
      # - `primary_care_covid_case_0_date` used to define `prior_covid_infection` but not `anycovid_0_date` 
      #    because it could refer to "history of" rather than "current"
    ) 
  
}

################################################################################
process_vax <- function(.data, stage) {
  
  data_vax <- local({
    
    data_vax_pfizer <- .data %>%
      select(patient_id, matches("covid\\_vax\\_pfizer\\_\\d+\\_date")) %>%
      pivot_longer(
        cols = -patient_id,
        names_to = c(NA, "vax_pfizer_index"),
        names_pattern = "^(.*)_(\\d+)_date",
        values_to = "date",
        values_drop_na = TRUE
      ) %>%
      arrange(patient_id, date)
    
    data_vax_az <- .data %>%
      select(patient_id, matches("covid\\_vax\\_az\\_\\d+\\_date")) %>%
      pivot_longer(
        cols = -patient_id,
        names_to = c(NA, "vax_az_index"),
        names_pattern = "^(.*)_(\\d+)_date",
        values_to = "date",
        values_drop_na = TRUE
      ) %>%
      arrange(patient_id, date)
    
    data_vax_moderna <- .data %>%
      select(patient_id, matches("covid\\_vax\\_moderna\\_\\d+\\_date")) %>%
      pivot_longer(
        cols = -patient_id,
        names_to = c(NA, "vax_moderna_index"),
        names_pattern = "^(.*)_(\\d+)_date",
        values_to = "date",
        values_drop_na = TRUE
      ) %>%
      arrange(patient_id, date)
    
    
    data_vax <-
      data_vax_pfizer %>%
      full_join(data_vax_az, by=c("patient_id", "date")) %>%
      full_join(data_vax_moderna, by=c("patient_id", "date")) %>%
      mutate(
        type = fct_case_when(
          (!is.na(vax_az_index)) & is.na(vax_pfizer_index) & is.na(vax_moderna_index) ~ "az",
          is.na(vax_az_index) & (!is.na(vax_pfizer_index)) & is.na(vax_moderna_index) ~ "pfizer",
          is.na(vax_az_index) & is.na(vax_pfizer_index) & (!is.na(vax_moderna_index)) ~ "moderna",
          (!is.na(vax_az_index)) + (!is.na(vax_pfizer_index)) + (!is.na(vax_moderna_index)) > 1 ~ "duplicate",
          TRUE ~ NA_character_
        )
      ) %>%
      arrange(patient_id, date) %>%
      group_by(patient_id) %>%
      mutate(
        vax_index=row_number()
      ) %>%
      ungroup()
    
    data_vax
    
  })
  
  data_vax_wide = data_vax %>%
    pivot_wider(
      id_cols= patient_id,
      names_from = c("vax_index"),
      values_from = c("date", "type"),
      names_glue = "covid_vax_{vax_index}_{.value}"
    )
  
  # only add variables corresponding to 4th dose if stage = treated
  if (stage == "treated") {
    vax4_vars <- rlang::quos(
      
      vax4_type = covid_vax_4_type,
      
      vax4_type_descr = fct_case_when(
        vax4_type == "pfizer" ~ "BNT162b2",
        vax4_type == "az" ~ "ChAdOx1",
        vax4_type == "moderna" ~ "Moderna",
        TRUE ~ NA_character_
      ),
      
      vax4_date = covid_vax_4_date,
      
    )
  } else if (stage == "potential") {
    vax4_vars <- rlang::quos(
      NULL
    )
  }
  
  .data %>%
    left_join(data_vax_wide, by ="patient_id") %>%
    mutate(
      vax1_type = covid_vax_1_type,
      vax2_type = covid_vax_2_type,
      vax3_type = covid_vax_3_type,
      
      vax12_type = paste0(vax1_type, "-", vax2_type),
      
      vax1_type_descr = fct_case_when(
        vax1_type == "pfizer" ~ "BNT162b2",
        vax1_type == "az" ~ "ChAdOx1",
        vax1_type == "moderna" ~ "Moderna",
        TRUE ~ NA_character_
      ),
      vax2_type_descr = fct_case_when(
        vax2_type == "pfizer" ~ "BNT162b2",
        vax2_type == "az" ~ "ChAdOx1",
        vax2_type == "moderna" ~ "Moderna",
        TRUE ~ NA_character_
      ),
      vax3_type_descr = fct_case_when(
        vax3_type == "pfizer" ~ "BNT162b2",
        vax3_type == "az" ~ "ChAdOx1",
        vax3_type == "moderna" ~ "Moderna",
        TRUE ~ NA_character_
      ),
      
      vax12_type_descr = paste0(vax1_type_descr, "-", vax2_type_descr),
      
      vax1_date = covid_vax_1_date,
      vax2_date = covid_vax_2_date,
      vax3_date = covid_vax_3_date,
      
      vax12_gap = as.integer(vax2_date - vax1_date), 
      # vax12_gap is a covariate in the model - do we want to keep it continuous, add a quadratic term, categorise it?
      
      # day of second dose relative to start of vaccination rollout (used in matching)
      vax2_day = as.integer(floor((vax2_date - as.Date("2020-12-08")))),
      
      !!! vax4_vars
      
    ) %>%
    select(
      -starts_with("covid_vax_"),
    ) 
  
}

################################################################################
process_outcome <- function(.data) {
  
  .data %>%
    mutate(
      
      # earliest covid event after study start
      anycovid_date = pmin(postest_date, covidemergency_date, covidadmitted_date, covidcritcare_date, coviddeath_date, na.rm=TRUE),
      
      noncoviddeath_date = if_else(!is.na(death_date) & is.na(coviddeath_date), death_date, as.Date(NA_character_)),
      # cvd or cancer deaths must be non-covid
      # cvddeath_date = if_else(!is.na(cvddeath_date) & is.na(coviddeath_date), death_date, as.Date(NA_character_)),
      # cancerdeath_date = if_else(!is.na(cancerdeath_date) & is.na(coviddeath_date), death_date, as.Date(NA_character_)),
      
      covidcritcareordeath_date = pmin(covidcritcare_date, coviddeath_date, na.rm=TRUE),
      
      fracture_date = pmin(fractureemergency_date, fractureadmitted_date, fracturedeath_date, na.rm=TRUE)
      
    )
}
