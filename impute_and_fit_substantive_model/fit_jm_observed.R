library(survival) # for Surv
library(JMbayes2)
library(tidyverse)

load("code/dat.Rdata")
n_burnin <- 100000
n_iter <- 200000
n_thin <- 100
seed <- 1765482

set.seed(seed)

# run estimation ----------------------------------------------------------
dat_long_k <- dat_long %>%
  filter(!is.na(EM_sc))

fit_lme_sc20 <- nlme::lme(
  EM_sc ~ poly(age_sc20 - mean_age_sc20, 3, raw = TRUE) +
    sex + educ_imp_sc + ind_enrol +
    I(ind_enrol * (age_sc20 - mean_age_sc20)) +
    cohort_sc,
  random = ~ I(age_sc20 - mean_age_sc20) | id,
  data = dat_long_k %>% arrange(id), method = "ML"
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
  functional_forms = ~ value(EM_sc) + slope(EM_sc),
  data_Surv = dat_surv %>% arrange(id),
  n_burnin = n_burnin, n_iter = n_iter, n_thin = n_thin,
  control = list(cores = 3, seed = seed)
)

dat_long_surv <- list(dat_long_k %>% arrange(id), dat_surv %>% arrange(id))
save(fit_jm, seed, dat_long_surv, dat_long_k, file = "code/results/fit_jm_observed.Rdata")

