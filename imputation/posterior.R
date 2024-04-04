library(tidyverse)
library(splines)

# read data ---------------------------------------------------------------
load("dat.Rdata")
load("results/fit_jm_pattern_mixture_saved_random_effects.Rdata")

# impute ------------------------------------------------------------------
set.seed(0)
basis_in_the_estimation <- ns(dat_long %>% filter(!is.na(EM)) %>% arrange(id) %>% pull(age) -
  dat_long %>% filter(!is.na(EM)) %>% arrange(id) %>% pull(mean_age), df = 3)

model_matrix_lme <- stats::model.matrix(
  ~ ns(age - mean_age,
    df = 3,
    knots = attributes(basis_in_the_estimation)$knots,
    Boundary.knots = attributes(basis_in_the_estimation)$Boundary.knots
  ) * missing_data_pattern +
    sex +
    educ_imp +
    ind_enrol +
    I(ind_enrol * age) +
    cohort,
  random = ~ I(age - mean_age) | id,
  data = dat_long %>% arrange(id), method = "ML"
)

beta_posterior <- do.call(rbind, fit_jm[["mcmc"]][["betas1"]])
# Extract random intercepts' posterior
# For this, access the MCMC samples for parameter "b", extract the first column of each matrix
random_intercept_posterior <- do.call(
  cbind, # Combine matrices column-wise
  lapply(fit_jm[["mcmc"]][["b"]], function(x) {
    x[, 1, ]
  })
)
# Extract random slopes' posterior
# For this, access the MCMC samples for parameter "b", extract the second column of each matrix
random_slope_posterior <- do.call(
  cbind,
  lapply(fit_jm[["mcmc"]][["b"]], function(x) {
    x[, 2, ]
  })
)
dat_long_with_posterior <- dat_long
for (i in 1:nrow(beta_posterior)) {
  dat_long_with_posterior[, paste0("post", i)] <-
    model_matrix_lme %*% beta_posterior[i, ] + # fixed effects
    random_intercept_posterior[match(dat_long$id, unique(dat_long$id, )), i] +
    (dat_long$age - dat_long$mean_age) * random_slope_posterior[match(dat_long$id, unique(dat_long$id, )), i]
}
save(dat_long_with_posterior, file = "results/dat_long_with_posterior.Rdata")
