################################################################################
# functions for processing each of the variable groups in analysis/process_data.R

################################################################################
process_jcvi <- function(.data) {
  .data %>%
    mutate(
      
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
      
      multimorb =
        (sev_obesity) +
        (chronic_heart_disease) +
        (chronic_kidney_disease)+
        (diabetes) +
        (chronic_liver_disease)+
        (chronic_resp_disease | asthma)+
        (chronic_neuro_disease)#+
      #(learndis)+
      #(sev_mental),
      ,
      multimorb = cut(multimorb, breaks = c(0, 1, 2, Inf), labels=c("0", "1", "2+"), right=FALSE),
      immuno = immunosuppressed | asplenia,
      
      
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
        care_home_combined | hscworker  ~ "1",
        age_aug2021>=80 ~ "2",
        age_aug2021>=75 ~ "3",
        age_aug2021>=70 | (cev & (age_aug2021>=16)) ~ "4",
        age_aug2021>=65 ~ "5",
        between(age_aug2021, 16, 64.999) & cv ~ "6",
        age_aug2021>=60 ~ "7",
        age_aug2021>=55 ~ "8",
        age_aug2021>=50 ~ "9",
        age_aug2021>=40 ~ "10a",
        TRUE ~ "10b"
      ),
      
      jcvi_group_descr = fct_recode(
        jcvi_group,
        "Care home residents and health and social care workers"="1",
        "80+ years"="2",
        "75-79 years"="3",
        "70-74 years or clinically extremely vulnerable"="4",
        "65-69 years"="5",
        "16-64 years or clinically at-risk"="6",
        "60-64 years"="7",
        "55-59 years"="8",
        "50-54 years"="9",
        "40-49 years"="10a",
        "16-39 years"="10b"
      ),
    )
}


################################################################################
process_demographic <- function(.data) {
  .data %>%
    mutate(
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
      
      ethnicity_combined = if_else(is.na(ethnicity), ethnicity_6_sus, ethnicity),
      
      ethnicity_combined = fct_case_when(
        ethnicity_combined == "1" ~ "White",
        ethnicity_combined == "4" ~ "Black",
        ethnicity_combined == "3" ~ "South Asian",
        ethnicity_combined == "2" ~ "Mixed",
        ethnicity_combined == "5" ~ "Other",
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
      
      imd_Q5 = factor(imd_Q5, levels = c("1 (most deprived)", "2", "3", "4", "5 (least deprived)", "Unknown")),
      
      rural_urban_group = fct_case_when(
        rural_urban %in% c(1,2) ~ "Urban conurbation",
        rural_urban %in% c(3,4) ~ "Urban city or town",
        rural_urban %in% c(5,6,7,8) ~ "Rural town or village",
        TRUE ~ NA_character_
      ),
    ) %>%
    select(-ethnicity, -ethnicity_6_sus)
}

################################################################################
process_pre <- function(.data) {
  .data %>%
    mutate(
      prior_tests_cat = cut(prior_covid_test_frequency, breaks=c(0, 1, 2, 3, Inf), labels=c("0", "1", "2", "3+"), right=FALSE),
      # any covid event before study start
      prior_covid_infection = (!is.na(positive_test_0_date)) | (!is.na(admitted_covid_0_date)) | (!is.na(primary_care_covid_case_0_date)),
      # date of latest covid event before study start
      anycovid_0_date = pmax(positive_test_0_date, covidemergency_0_date, admitted_covid_0_date, na.rm=TRUE),
      # any reason for the discrepancy between events used to define prior_covid_infection and anycovid_0_date?
    )
}

################################################################################
process_post <- function(.data) {
  .data %>%
    mutate(
      #covidemergency_1_date = pmin(covidemergency_1_date, covidadmitted_1_date, na.rm=TRUE),
      
      # because this value is returned as a factor by the study definition
      admitted_covid_ccdays_1 = as.numeric(as.character(admitted_covid_ccdays_1)),
      admitted_covid_ccdays_2 = as.numeric(as.character(admitted_covid_ccdays_2)),
      admitted_covid_ccdays_3 = as.numeric(as.character(admitted_covid_ccdays_3)),
      admitted_covid_ccdays_4 = as.numeric(as.character(admitted_covid_ccdays_4)),
      
      
      covidcc_1_date = case_when(
        admitted_covid_ccdays_1 > 0 ~ admitted_covid_1_date,
        admitted_covid_ccdays_2 > 0 ~ admitted_covid_2_date,
        admitted_covid_ccdays_3 > 0 ~ admitted_covid_3_date,
        admitted_covid_ccdays_4 > 0 ~ admitted_covid_4_date,
        TRUE ~ as.Date(NA_character_)
      ),
      
      covidcc_2_date = case_when(
        (admitted_covid_ccdays_2 > 0) & (admitted_covid_2_date > covidcc_1_date) ~ admitted_covid_2_date,
        (admitted_covid_ccdays_3 > 0) & (admitted_covid_3_date > covidcc_1_date) ~ admitted_covid_3_date,
        (admitted_covid_ccdays_4 > 0) & (admitted_covid_4_date > covidcc_1_date) ~ admitted_covid_4_date,
        TRUE ~ as.Date(NA_character_)
      ),
      
      covidcc_3_date = case_when(
        (admitted_covid_ccdays_3 > 0) & (admitted_covid_3_date > covidcc_2_date) ~ admitted_covid_3_date,
        (admitted_covid_ccdays_4 > 0) & (admitted_covid_4_date > covidcc_2_date) ~ admitted_covid_4_date,
        TRUE ~ as.Date(NA_character_)
      ),
      
      covidcc_4_date = case_when(
        (admitted_covid_ccdays_4 > 0) & (admitted_covid_4_date > covidcc_3_date) ~ admitted_covid_4_date,
        TRUE ~ as.Date(NA_character_)
      ),
      
      # earliest covid event after study start
      anycovid_1_date = pmin(positive_test_1_date, covidemergency_1_date, admitted_covid_1_date, covidcc_1_date, coviddeath_date, na.rm=TRUE),
      
      noncoviddeath_date = if_else(!is.na(death_date) & is.na(coviddeath_date), death_date, as.Date(NA_character_)),
      
      cause_of_death = fct_case_when(
        !is.na(coviddeath_date) ~ "covid-related",
        !is.na(death_date) ~ "not covid-related",
        TRUE ~ NA_character_
      ),
    )
}
