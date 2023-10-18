library(survival) # for Surv
library(JMbayes2)
library(tidyverse)

source("code/impute_and_fit_substantive_model/fit_jm_parameters.R")
# source("fit_jm_parameters.R")
# read data and lme fit with imputations ----------------------------------
load(file_data)
load(file_predictions)

set.seed(seed) 
r <- sample(1:1000,
            size = number_of_multiple_imputations,
            replace = T) # select posterior iterations to use in imputations
seeds <- sample(1:1000000, size =  number_of_multiple_imputations )


dat_long_surv <-  list()

# run estimation ----------------------------------------------------------
for (k in 1:number_of_multiple_imputations) {
  print(k)
  dat_long_k <- dat_long %>%
    left_join(dat_predicted%>%dplyr::select(c("id", "age_sc20", !!as.name(r[k]))), by = c("id", "age_sc20"))%>%
    group_by(id) %>%          
    mutate(age_dropout_sc20 = max(age_sc20[!is.na(EM_sc)]))%>% # calculate  scaled age at dropout
    ungroup()%>%
    mutate(t_minus_t_dropout_pos = (age_sc20 - age_dropout_sc20)*(age_sc20 - age_dropout_sc20>0),
           predicted_EM_delta = !!as.name(r[k]) + delta*t_minus_t_dropout_pos, # predicted by a pattern-mixture + delta*(t-t_ifi)+
           predicted_EM_lower0 = predicted_EM_delta * EM_cc_sd + EM_cc_mean < 0,  # check if predicted EM_sc corresponds to EM lower than 0
           imputed_sc = ifelse(predicted_EM_lower0, - EM_cc_mean/EM_cc_sd, predicted_EM_delta), # substitute those lower then 0 by scaled measures that corresponds to 0.
           imputed = imputed_sc * EM_cc_sd + EM_cc_mean,
           EM_sc_imp = ifelse(!is.na(EM_sc), EM_sc,  imputed_sc))%>% # imputations are named with numbers. 
    # !! tells R to treat the enclosed expression as an actual variable named k, rather than as a symbol k.
    filter(!is.na(EM_sc_imp))
  
  fit_lme_sc20 <- nlme::lme(EM_sc_imp ~ poly(age_sc20 - mean_age_sc20, 3, raw = TRUE) +
                              sex + educ_imp_sc + ind_enrol +
                              I(ind_enrol * (age_sc20 - mean_age_sc20)) +
                              cohort_sc,
                            random = ~ I(age_sc20 - mean_age_sc20) | id,
                            data = dat_long_k %>% arrange(id), method = "ML"
  )
  
  fit_cox <- survival::coxph(Surv(
    time = age_enrol_sc25,
    time2 = age_event_sc25,
    event = status, origin = 0
  ) ~ sex + educ_imp_sc,
  data = dat_surv %>% arrange(id),
  x = TRUE, model = TRUE)
  
  fit_jm <- JMbayes2::jm(
    Surv_object = fit_cox,
    Mixed_objects = fit_lme_sc20,
    time_var = "age_sc20",
    functional_forms = ~ value(EM_sc_imp) + slope(EM_sc_imp),
    data_Surv = dat_surv %>% arrange(id),
    n_burnin = n_burnin, n_iter = n_iter, n_thin = n_thin,  
    control = list(cores = 3, seed = seeds[k])
  )
  
  dat_long_surv = list(dat_long_k %>% arrange(id), dat_surv  %>% arrange(id))
  save(fit_jm, delta, seed, seeds, r, dat_long_surv, dat_long_k, file_data, file_predictions, file = paste0(output_file_name, "_", k, ".Rdata"))
  rm(fit_jm)
}


