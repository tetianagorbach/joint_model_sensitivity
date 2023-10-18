library(survival) # for Surv
library(tidyverse)

# read data ---------------------------------------------------------------
load("code/dat.Rdata")
load("code/results/fit_jm_pattern_mixture100_200_100.Rdata")

# impute ------------------------------------------------------------------
id_miss <- dat_long%>%filter(is.na(EM_sc))%>%distinct(id)%>%pull(id)
set.seed(34567)
seeds <- round(runif(length(id_miss))*100000)
dat_predicted_raw <- lapply(1:length(id_miss), function(x){
  print(x)
  dat_id <- dat_long%>%
    filter(id %in% id_miss[x])
  mcmc <-  predict(object = fit_jm,
          newdata = dat_id%>%filter(!is.na(EM_sc)),
          times = dat_id%>%filter(is.na(EM_sc))%>%pull(age_sc20),
          process = "longitudinal", 
          type = "subject_specific", 
          return_mcmc = T, 
          n_samples = 1000,
          seed = seeds[x]
  )
  other_data <-  predict(object = fit_jm,
                           newdata = dat_id%>%filter(!is.na(EM_sc)),
                           times = dat_id%>%filter(is.na(EM_sc))%>%pull(age_sc20),
                           process = "longitudinal",
                           type = "subject_specific", 
                           return_newdata = T, 
                           n_samples = 1000,
                           seed = seeds[x]
  )
  
  list(mcmc = mcmc[[2]], other_data = other_data[[2]], seed = seeds[x])
})

dat_predicted <- sapply(1:length(dat_predicted_raw), function(i){
  dat_predicted_raw_i <- dat_predicted_raw[[i]][["other_data"]][-1,]%>% # id and ages
    dplyr::select(id, age_sc20)
  
  cbind(dat_predicted_raw_i, # add imputed EM, 
        matrix(dat_predicted_raw[[i]][["mcmc"]][["mcmc"]][["EM_sc"]][-1,], 
               nrow = nrow(dat_predicted_raw_i)
               ) # "matrix" is needed for correct cbind when only one measure. is imputed (in such a case dat_predicted_raw[[i]][["mcmc"]][["mcmc"]][["EM_sc"]][-1,] is numeric)
        )
}
, simplify = F
)
dat_predicted <- do.call("rbind", dat_predicted)
save(seeds, dat_predicted_raw, dat_predicted,  file = "code/results/dat_predicted_raw.Rdata")




