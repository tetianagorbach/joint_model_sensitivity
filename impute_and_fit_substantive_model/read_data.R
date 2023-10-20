library(tidyverse)
library(lubridate)

# read data and change ids to the simulated -------------------------------
dat_everything <- readxl::read_excel(
  path = "/Users/gote0002/Desktop/B_TG_allaT_AgeD_kl.xlsx",
  sheet = 1,
  na = "999"
)

ids <- readr::read_csv(
  file = "/Users/gote0002/Desktop/unikt_nr_to_id.csv",
  col_types = c("d", "d")
)

# check that each initial id corresponds to only one simulated id
base::setdiff(dat_everything$unikt_nr, ids$unikt_nr) # initial ids are the same in initial data and the code file
ids$unikt_nr[duplicated(ids$unikt_nr)] # only unique values of initial ids in the code file
ids$id[duplicated(ids$id)] # only unique values of simulated ids in the code file

# substitute initial ids by the simulated ids in the data
dat_everything <- ids %>% # add data to the Betula and new ids
  left_join(dat_everything, by = "unikt_nr") %>%
  dplyr::select(-c(unikt_nr)) %>% # delete Betula id
  arrange(id, test_wave)
rm(ids)

# delete 1374 rows with no data for any of the five memory tests (9733 rows left)
dat_long <- dat_everything %>%
  rowwise() %>%
  mutate(na_count = sum(is.na(c_across(all_of(c(
    "SPT_VTCategoryCuedRecall::sptcrc", "SPT_VTCategoryCuedRecall::vtcrc",
    "SPT_VTFreeRecall::sptb", "SPT_VTFreeRecall::vtb", "WorkingMemory::wrk00"
  )))))) %>%
  filter(!na_count == 5)


# select needed variables -------------------------------------------------
dat_long <- dat_long %>%
  dplyr::select(
    id, `Participant::sex`,
    sample = `Participant::sample`,
    sex = `Participant::sex`,
    educ = `educ_T`,
    allele4 = `Genes::ApoE_fromGWAS_E4binary_FmCALC`,
    test_wave, age_cohort_T,
    rounded_age_HT_1decimal,
    rounded_age_MT_1decimal,
    MMSE = `MMT::v518`,
    EMC_FM, `SPT_VTCategoryCuedRecall::sptcrc`, `SPT_VTCategoryCuedRecall::vtcrc`,
    `SPT_VTFreeRecall::sptb`, `SPT_VTFreeRecall::vtb`, `WorkingMemory::wrk00`,
    age_dementia = `Demens::ageAtOnset`,
    `Deceased::ageAtEventOfDeath_y`,
    `Deceased::ageAtEventOfDeath_y`,
    `Deceased::vitalStatusLastUpdate_YYYY`,
    `Deceased::vitalStatusLastUpdate_MM`,
    `Deceased::vitalStatusLastUpdate_DD`,
    contains("Demens::"),
    `Demens::ageAtOnset`,
    ht_testdate_year, ht_testdate_month, HalfOfMonth_HT,
    mt_testdate_year, mt_testdate_month, HalfOfMonth_MT, `Deceased::deceasedDate`
  ) %>%
  mutate(
    sex = sex - 1, # initially 1 for women, 2 for men. Specify sex to be 0 for females and 1 for males
    EM = `SPT_VTCategoryCuedRecall::sptcrc` + `SPT_VTCategoryCuedRecall::vtcrc` +
      `SPT_VTFreeRecall::sptb` + `SPT_VTFreeRecall::vtb` + `WorkingMemory::wrk00`, # define memory score
    age = ifelse(!is.na(rounded_age_MT_1decimal), rounded_age_MT_1decimal, rounded_age_HT_1decimal),
    date_rounded = if_else(!is.na(rounded_age_MT_1decimal),
      as.Date(make_date(
        mt_testdate_year,
        mt_testdate_month,
        ifelse(HalfOfMonth_MT == "first", 1, 15)
      )),
      as.Date(make_date(
        ht_testdate_year,
        ht_testdate_month,
        ifelse(HalfOfMonth_HT == "first", 1, 15)
      ))
    ),
    date_birth = date_rounded - dyears(age), # date of birth by subtracting age from a date of testing.
    age_death_init = age + time_length( # a sum of age at testing and the difference between the date of death and the age at testing.
      difftime(
        as.Date(`Deceased::deceasedDate`),
        date_rounded
      ),
      "years"
    ),
    age_death_eval_init = age + time_length(
      difftime(
        as.Date(make_date(
          `Deceased::vitalStatusLastUpdate_YYYY`,
          `Deceased::vitalStatusLastUpdate_MM`,
          `Deceased::vitalStatusLastUpdate_DD`
        )),
        date_rounded
      ), "years"
    ),
    age_dementia_eval_init = age + time_length( # sum of age and the difference between the date of the latest dementia diagnostic evaluation and the date of testing
      difftime(
        `Demens::dateOfDiagnosticEvaluation`,
        date_rounded
      ), "years"
    )
  )


# delete people with early onset of dementia
id_early_onset <- unique(dat_long[dat_long$age_dementia < 65, ]$id)

dat_long <- dat_long %>%
  filter(!(id %in% id_early_onset))

# d <- dat_long%>%
#   filter(is.na(EM) & !is.na(EMC_FM))%>%
#   dplyr::select(id, `SPT_VTCategoryCuedRecall::sptcrc`, `SPT_VTCategoryCuedRecall::vtcrc`,
#                 `SPT_VTFreeRecall::sptb`, `SPT_VTFreeRecall::vtb`, `WorkingMemory::wrk00`)%>%
#   rowwise() %>%
#   mutate(na_count = sum(is.na(c_across(all_of(c("SPT_VTCategoryCuedRecall::sptcrc", "SPT_VTCategoryCuedRecall::vtcrc",
#                                                 "SPT_VTFreeRecall::sptb", "SPT_VTFreeRecall::vtb", "WorkingMemory::wrk00"))))))
# table(d$na_count)

# Delete 197 rows with missing memory scores
dat_long <- dat_long %>%
  filter(!is.na(EM))

# dat_long <- dat_long %>%
#   filter(!is.na(EM)) %>%
#   filter(!(id %in% id_early_onset)) %>%
#   group_by(id)%>%
#   mutate(age_dementia_range = ifelse(sum(!is.na(age_dementia))>0, max(age_dementia, na.rm=T) - min(age_dementia, na.rm=T), NA ),
#          age_death_range = ifelse(sum(!is.na(age_death_init))>0, max(age_death_init, na.rm=T) - min(age_death_init, na.rm=T), NA ),
#          age_dementia_eval_range = ifelse(sum(!is.na(age_dementia_eval_init))>0, max(age_dementia_eval_init, na.rm=T) - min(age_dementia_eval_init, na.rm=T), NA ),
#          age_death_eval_range = ifelse(sum(!is.na(age_death_eval_init))>0, max(age_death_eval_init, na.rm=T) - min(age_death_eval_init, na.rm=T), NA ))%>%
#   ungroup()
# max(dat_long$age_dementia_range, na.rm = T) # age at dementia onset is the same for all observations of a person
# max(dat_long$age_death_range, na.rm = T) # age of death  differs  by max 0.13 years between observations of a person
# max(dat_long$age_dementia_eval_range, na.rm = T) # age  at dementia evaluation differs  by max 0.14 years between observations of a person
# max(dat_long$age_death_eval_range, na.rm = T) # age  at death evaluation differs  by max 0.14 years between observations of a person

dat_surv_everything <- dat_long %>%
  group_by(id) %>%
  summarise(
    sex = unique(sex[!is.na(sex)]),
    sample = mean(sample),
    age_enrol = age[age == min(age)],
    age_dementia = ifelse(sum(!is.na(age_dementia)) > 0, max(age_dementia, na.rm = T), NA), # if no age at dementia onset recorded for any waves - NA, otherwise max age at dementia onset
    age_death = ifelse(sum(!is.na(age_death_init)) > 0, max(age_death_init, na.rm = T), NA),
    age_dementia_eval = ifelse(sum(!is.na(age_dementia_eval_init)) > 0, max(age_dementia_eval_init, na.rm = T), NA),
    age_death_eval = ifelse(sum(!is.na(age_death_eval_init)) > 0, max(age_death_eval_init, na.rm = T), NA),
    age_last_obs = max(age),
    waves_obs = paste0(test_wave, collapse = ""),
    max_wave_obs = max(test_wave),
    waves_obs_before_dementia = paste0(test_wave[is.na(age_dementia) | age < age_dementia], collapse = "")
  ) %>% # age at the last available longitudinal observation
  ungroup() %>%
  mutate(
    status = !is.na(age_dementia), # status is 1 if demented
    age_censor = ifelse(!is.na(age_dementia), # age censor is NA if demented
      NA,
      ifelse(!is.na(pmin(age_dementia_eval, age_death, na.rm = T)), # if not demented, min of age at dementia evaluation or death
        pmin(age_dementia_eval, age_death, na.rm = T),
        age_last_obs + 0.003 # if no age at dementia onset, death or dementia evaluation are available, time of the last observation + approx 1 day
      )
    ),
    age_event = ifelse(!is.na(age_dementia), age_dementia, age_censor),
  )

# define missing data pattern -------------------------------------
mis_pattern <- unique(dat_surv_everything %>% dplyr::select(sample, waves_obs, max_wave_obs)) %>%
  arrange(sample, waves_obs) %>%
  mutate(
    missing_data_pattern = ifelse(((sample == 1 | sample == 3) & max_wave_obs == 6),
      "0.a",
      ifelse((sample == 2 & waves_obs %in% c(23, 235, 25, 3)) |
        sample == 4 | sample == 5 | (sample == 6 & max_wave_obs == 6),
      "0.b",
      ifelse(
        (sample == 1 & (max_wave_obs == 2 | max_wave_obs == 1)) |
          sample == 2 |
          (sample == 3 & (max_wave_obs == 3 | max_wave_obs == 2)) |
          (sample == 6 & max_wave_obs == 5), 2, NA
      )
      )
    ),
    missing_data_pattern = ifelse(!is.na(missing_data_pattern),
      missing_data_pattern,
      ifelse(sample == 1, max_wave_obs,
        ifelse(sample == 3, max_wave_obs - 1, NA)
      )
    ),
    missing_data_pattern = factor(missing_data_pattern,
      levels = c("0.a", "5", "4", "3", "2", "0.b")
    ),
    waves_unobs = ifelse(missing_data_pattern %in% c("0.a", "0.b"), NA,
      ifelse(sample == 1 | sample == 3 | sample == 6, substring("123456", max_wave_obs + 1),
        ifelse(sample == 2 & waves_obs %in% c(2, 25), 3, NA)
      )
    )
  )

# add info about unobserved waves and missing data pattern to the survival data
dat_surv_everything <- dat_surv_everything %>%
  left_join(mis_pattern %>% dplyr::select(sample, waves_obs, waves_unobs, missing_data_pattern),
    by = c("sample", "waves_obs")
  )


## Add rows with NA for EM for "missing" observations
dat_long_to_impute <- dat_long %>%
  dplyr::select(id, sample, test_wave, age) %>%
  left_join(dat_surv_everything %>% dplyr::select(id, waves_unobs, age_dementia, age_death, age_event),
    by = "id"
  ) %>%
  filter(!is.na(waves_unobs)) %>% # do not impute for those people that have all observations
  filter(is.na(age_dementia) | age < age_dementia) %>%
  filter(age < age_death | is.na(age_death)) %>%
  group_by(id) %>%
  top_n(1, test_wave) %>%
  mutate(
    age_dropout = age
  ) %>% # calculate  age at dropout
  ungroup() # 2179 out of 4394 participants

# d <- dat_long_to_impute%>%
#   mutate(age_dementia_minus_age = age_dementia - age ,
#          age_death_minus_age = age_death - age)
# sum(d$age_dementia_minus_age<5 | d$age_death_minus_age<5 , na.rm=T) #713 people got demented or died no later than 5 years after the last observed EM score - nothing to impute for them

dat_long_to_impute <- dat_long_to_impute %>%
  slice(rep(row_number(), nchar(waves_unobs))) %>%
  add_column(test_wave_unobs = unlist(sapply(dat_long_to_impute$waves_unobs, function(x) as.numeric(unlist(strsplit(as.character(x), split = "")))))) %>%
  mutate(age = age + 5 * (test_wave_unobs - test_wave)) %>%
  filter(age < age_dementia | is.na(age_dementia)) %>%
  filter(age < age_death | is.na(age_death)) # do not impute after death or dementia, but impute after censoring time

dat_long_with_missing <- dat_long %>%
  full_join(dat_long_to_impute %>% dplyr::select(id, age, age_dropout), by = c("id", "age")) %>% # Add to longitudinal data rows with NA for EM for "missing" observations
  # specify education
  group_by(id) %>%
  mutate(
    sex = unique(sex[!is.na(sex)]),
    n_unique_educ = length(unique(educ)),
    educ_imp = ifelse(!is.na(educ[test_wave == min(test_wave)]), educ[test_wave == min(test_wave)],
      ifelse(length(educ[!is.na(educ)]) == 0, NA, min(educ, na.rm = T))
    ), # education is baseline education,
    # for 3 persons baseline education is NA, but at Betula wave 6 is 0, use 0.
    # for one person baseline education is NA, but at Betula wave 3,6 is 5, use 5.
    ind_enrol = (age == min(age)),
    age_enrol = age[age == min(age)],
    age_mci = ifelse(any(MMSE < 24) == T,
      min(age[MMSE < 24], na.rm = T),
      NA
    ),
    age_dropout = ifelse(length(age_dropout[!is.na(age_dropout)]) == 0, NA, min(age_dropout, na.rm = T)),
    age_minus_age_dropout_pos = (age - age_dropout) * ((age - age_dropout) > 0),
    cohort = mean(as.numeric(format(date_birth, "%Y")), na.rm = T)
  ) %>%
  left_join(dat_surv_everything %>% dplyr::select(id, missing_data_pattern, age_death, age_dementia_eval, age_death_eval, age_censor, age_event),
    by = "id"
  ) %>%
  ungroup()
dat_long_with_missing[dat_long_with_missing$id == 3857, "educ_imp"] <- NA # this subject has 70 as their number of years of education
# dat_long_with_missing[dat_long_with_missing$id == 3458 & dat_long_with_missing$test_wave == 3, "age_cohort_T"] <- 45 # For a subject 3458 with missing age_cohort_T put the value 45. but since the subject does not have EM measure, this does not matter
dat_long_with_missing <- dat_long_with_missing %>%
  filter(age < age_dementia | is.na(age_dementia)) %>% # do not use EM observations after dementia onset
  filter(age < age_death | is.na(age_death)) # do not use EM observations after death

rm(id_early_onset, dat_long_to_impute, mis_pattern, dat_long)


dat_long_with_missing <- dat_long_with_missing %>%
  mutate(across(c("age", "age_dropout"), ~ (.x - 20) / 100, .names = "{.col}_sc20"),
    EM_cc_mean = mean(EM[!is.na(EM) & !is.na(educ_imp)]),
    EM_cc_sd = sd(EM[!is.na(EM) & !is.na(educ_imp)]),
    EM_sc = (EM - EM_cc_mean) / EM_cc_sd, # scale memory,
    educ_cc_mean = mean(educ_imp[!is.na(EM) & !is.na(educ_imp)]),
    educ_cc_sd = sd(educ_imp[!is.na(EM) & !is.na(educ_imp)]),
    educ_imp_sc = (educ_imp - educ_cc_mean) / educ_cc_sd, # scale education
    mean_age_sc20 = (mean(age[!is.na(EM) & !is.na(educ_imp)]) - 20) / 100,
    cohort_sc = (cohort - 1937) / 100
  ) %>%
  arrange(id) # sort by id


dat_surv_everything <- dat_surv_everything %>% # add scaled education to the survival data
  left_join(dat_long_with_missing %>%
    distinct(id, educ_imp, educ_imp_sc), by = "id") %>%
  mutate(
    across(c("age_enrol", "age_event"),
      ~ (.x - 25) / 100,
      .names = "{.col}_sc25"
    ) # scaled age 0 corresponds to true age 25
  )
dat_surv <- dat_surv_everything %>%
  filter(!is.na(educ_imp)) %>%
  arrange(id)

dat_long <- dat_long_with_missing %>%
  filter(!is.na(educ_imp)) %>%
  dplyr::select(
    id, age, age_sc20, mean_age_sc20, age_cohort_T, EM, EM_sc, EM_cc_mean, EM_cc_sd,
    sex, educ, educ_imp, educ_imp_sc, educ_cc_mean, educ_cc_sd, ind_enrol, cohort_sc,
    missing_data_pattern
  )

# save(dat_long, dat_surv, file = "code/dat.Rdata")
