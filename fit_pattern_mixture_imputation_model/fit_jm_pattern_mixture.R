library(survival) # for Surv
library(tidyverse)
# read data ---------------------------------------------------------------
load("code/dat.Rdata")

# fit a pattern mixture model ---------------------------------------------------------
set.seed(0)
fit_lme_sc20 <- nlme::lme(
  EM_sc ~ poly(age_sc20 - mean_age_sc20, 3, raw = TRUE) * missing_data_pattern +
    sex +
    educ_imp_sc  +
    ind_enrol +
    I(ind_enrol * (age_sc20 - mean_age_sc20))  +
    cohort_sc ,
  random = ~ I(age_sc20 - mean_age_sc20) | id,
  data = dat_long %>% filter(!is.na(EM_sc)) %>% arrange(id), method = "ML"
)



fit_cox <- survival::coxph(
  Surv(
    time = age_enrol_sc25,
    time2 = age_event_sc25,
    event = status, origin = 0
  ) ~ sex + educ_imp_sc,
  data = dat_surv %>% arrange(id),
  x = TRUE, model = TRUE
)

fit_jm <- JMbayes2::jm(
  Surv_object = fit_cox,
  Mixed_objects = fit_lme_sc20,
  time_var = "age_sc20",
  data_Surv = dat_surv %>% arrange(id),
  functional_forms = ~ value(EM_sc) + slope(EM_sc),
  n_burnin = 100000, n_iter = 200000, n_thin = 10,
  control = list(cores = 3)
)
#save(fit_jm, fit_cox, fit_lme_sc20, file = "code/results/fit_jm_pattern_mixture100_200_10.Rdata")
