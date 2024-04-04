library(survival) # for Surv
library(JMbayes2)
library(tidyverse)

load("dat.Rdata")
n_burnin <- 20000
n_iter <- 50000
n_thin <- 30
seed <- 1765482

set.seed(seed)

# run estimation ----------------------------------------------------------
fit_lme <- nlme::lme(
  EM ~ ns(age -  mean_age, df = 3) +
    sex + educ_imp + ind_enrol +
    I(ind_enrol * age)+
    cohort,
  random = ~ I(age -  mean_age)  | id,
  data = dat_long %>%
    filter(!is.na(EM)) %>% arrange(id)
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
  functional_forms = ~ value(EM) + slope(EM),
  data_Surv = dat_surv %>% arrange(id),
  n_burnin = n_burnin, n_iter = n_iter, n_thin = n_thin,
  control = list(cores = 3, seed = seed, save_random_effects = T)
)

dat_long_surv <- list(dat_long %>% arrange(id), dat_surv %>% arrange(id))
save(fit_jm, seed, dat_long_surv, dat_long, file = "results/fit_jm_observed.Rdata")


