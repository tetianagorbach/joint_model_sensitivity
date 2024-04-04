library(survival) # for Surv
library(JMbayes2)
library(tidyverse)

source("substantive_model/fit_jm_parameters.R")
# source("fit_jm_parameters.R")
# read data and lme fit with imputations ----------------------------------
load(file_data)
load(file_predictions)
rm(delta)
rm(output_file_name)
output_file_name <- gsub(":| |-", "_", paste0("results/fit_jm_delta_changing_", Sys.time()))

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
    left_join(dat_long_with_posterior%>%dplyr::select(c("id", "age", paste0("post", r[k]))), by = c("id", "age"))%>%
    group_by(id) %>%          
    mutate(age_dropout = max(age[!is.na(EM)]))%>% # calculate  scaled age at dropout
    ungroup()%>%
    mutate(t_minus_t_dropout_pos = (age - age_dropout)*(age - age_dropout>0),
           delta = -1*((age-25)/75)^3,
           predicted_EM_delta = !!as.name(paste0("post", r[k])) + delta*t_minus_t_dropout_pos, # predicted by a pattern-mixture + delta*(t-t_ifi)+
           predicted_EM_lower0 = predicted_EM_delta  < 0,  # check if predicted EM_sc corresponds to EM lower than 0
           imputed = ifelse(predicted_EM_lower0, 0, predicted_EM_delta), # substitute those lower then 0 by scaled measures that corresponds to 0.
           EM_imp = ifelse(!is.na(EM), EM,  imputed))%>% # imputations are named with numbers. 
    # !! tells R to treat the enclosed expression as an actual variable named k, rather than as a symbol k.
    filter(!is.na(EM_imp))
  
  fit_lme <- nlme::lme(
    EM_imp ~ ns(age -  mean_age, df = 3) +
      sex + educ_imp + ind_enrol +
      I(ind_enrol * age)+
      cohort,
    random = ~ I(age -  mean_age)  | id,
    data = dat_long_k %>%
      filter(!is.na(EM_imp)) %>% arrange(id)
  )
  
  fit_cox <- survival::coxph(Surv(
    time = age_enrol_5,
    time2 = age_event_5,
    event = status, origin = 0
  ) ~ sex + educ_imp,
  data = dat_surv %>% arrange(id),
  x = TRUE, model = TRUE)
  
  fit_jm <- JMbayes2::jm(
    Surv_object = fit_cox,
    Mixed_objects = fit_lme,
    time_var = "age",
    functional_forms = ~ value(EM_imp) + slope(EM_imp),
    data_Surv = dat_surv %>% arrange(id),
    n_burnin = n_burnin, n_iter = n_iter, n_thin = n_thin,  
    control = list(cores = 3, seed = seeds[k])
  )
  
  dat_long_surv = list(dat_long_k %>% arrange(id), dat_surv  %>% arrange(id))
  save(fit_jm,  seed, seeds, r, dat_long_surv, dat_long_k, file_data, file_predictions, file = paste0(output_file_name, "_", k, ".Rdata"))
  rm(fit_jm)
}


