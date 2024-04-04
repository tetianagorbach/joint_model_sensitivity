library(survival) # for Surv
library(tidyverse)
library(splines)
# read data ---------------------------------------------------------------
load("dat.Rdata")


# fit a pattern mixture model ---------------------------------------------------------
set.seed(0)
fit_lme <- nlme::lme(
  EM ~ ns(age - mean_age,
    df = 3
  ) *
    missing_data_pattern +
    sex +
    educ_imp +
    ind_enrol +
    I(ind_enrol * age) +
    cohort,
  random = ~ I(age - mean_age) | id,
  data = dat_long %>% filter(!is.na(EM)) %>% arrange(id)
)

fit_cox <- survival::coxph(
  Surv(
    time = age_enrol_5,
    time2 = age_event_5,
    event = status, origin = 0
  ) ~ sex + educ_imp,
  data = dat_surv %>% arrange(id),
  x = TRUE, model = TRUE
)

fit_jm <- JMbayes2::jm(
  Surv_object = fit_cox,
  Mixed_objects = fit_lme,
  time_var = "age",
  data_Surv = dat_surv %>% arrange(id),
  functional_forms = ~ value(EM) + slope(EM),
  n_burnin = 50000, n_iter = 100000, n_thin = 50,
  control = list(cores = 3, save_random_effects = T, seed = 546271)
)
save(fit_jm, fit_cox, fit_lme, file = "results/fit_jm_pattern_mixture_saved_random_effects1.Rdata")
