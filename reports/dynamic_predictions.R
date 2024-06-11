library(tidyverse)
library(survival)
library(ggpubr)
library(splines)
library(gridExtra)
library(JMbayes2)
set.seed(11)
number_of_new_points_for_prediction <- 100
n_samples <-  1000 # 	the number of samples to use  in the prediction in predict.jm
n_mcmc <-  500 # 	 the number of Metropolis-Hastings iterations for sampling the random effects per iteration of n_samples in predict.jm

# Help functions ---------------------------------------------------------------
pool_pred <- function(vector_of_file_names, dat_id, method) { # pool predictions from the model fitted to completed after imputations data sets
  
  mcmc_pred_long <- lapply(vector_of_file_names,
                           FUN = function(x) {
                             load(x)
                             predsLong <- predict(fit_jm, # predict longitudinal data
                                                  newdata = dat_id,
                                                  times = seq(35, 100, length.out = number_of_new_points_for_prediction),
                                                  process = "longitudinal", return_mcmc = T,
                                                  type = "subject_specific",
                                                  n_samples = n_samples, n_mcmc = n_mcmc
                             )
                             
                             predsEvent <- predict(fit_jm, # predict CIF
                                                   newdata = dat_id,
                                                   times = seq(35, 100, length.out = number_of_new_points_for_prediction),
                                                   process = "event", return_mcmc = T,
                                                   n_samples = n_samples, n_mcmc = n_mcmc
                             )
                             rm(fit_jm)
                             list(  
                               predsLong[["newdata"]][["mcmc"]][["EM_imp"]],
                               predsLong[["newdata2"]][["mcmc"]][["EM_imp"]],
                               predsEvent[["mcmc"]],
                               predsEvent[["times"]],
                               dat_id
                             )
                           }
  )
  # mix all draws together
  draws_long_avail <- do.call(cbind, lapply(mcmc_pred_long, function(x) { # combine predicted memory for all imputations
    x[[1]]
  }))
  draws_long_new <- do.call(cbind, lapply(mcmc_pred_long, function(x) {
    x[[2]]
  }))
  draws_surv <- do.call(cbind, lapply(mcmc_pred_long, function(x) { # combine survival from all imputations
    x[[3]]
  }))
  times_new <- do.call(cbind, lapply(mcmc_pred_long, function(x) {
    x[[4]]
  }))
  level <- 0.95
  
  pooled_long <- data.frame( # combine predictions of the observed and the new longitudinal data
    method = method,
    times = c(mcmc_pred_long[[1]][[5]]%>%pull(age), rowMeans(times_new)),
    pred_EM = rowMeans(rbind(draws_long_avail, draws_long_new)),
    low = apply(rbind(draws_long_avail, draws_long_new),
                1,
                FUN = function(x) {
                  quantile(x, probs = (1 - level) / 2)
                }
    ),
    upp = apply(rbind(draws_long_avail, draws_long_new),
                1,
                FUN = function(x) {
                  quantile(x, probs = (1 + level) / 2)
                }
    )
  )
  pooled_surv <- data.frame(# combined predictions of the observed and the new survival data
    method = method,
    times = rowMeans(times_new),
    pred_CIF = rowMeans(draws_surv),
    low = apply(draws_surv,
                1,
                FUN = function(x) {
                  quantile(x, probs = (1 - level) / 2)
                }
    ),
    upp = apply(draws_surv,
                1,
                FUN = function(x) {
                  quantile(x, probs = (1 + level) / 2)
                }
    )
  )
  list(predicted_longitudinal = pooled_long, predicted_survival = pooled_surv)
}
pool_pred_obs <- function(vector_of_file_names, dat_id, method) { # pool predictions from the model fitted to the observed data
  mcmc_pred_long <- lapply(vector_of_file_names,
                           FUN = function(x) {
                             load(x)
                             predsLong <- predict(fit_jm, #predict longitudinal data
                                                  newdata = dat_id,
                                                  times = seq(35, 100, length.out = number_of_new_points_for_prediction),
                                                  process = "longitudinal", return_mcmc = T,
                                                  type = "subject_specific", 
                                                  n_samples = n_samples, n_mcmc = n_mcmc
                             )
                             rm(fit_jm)
                             list(
                               predsLong[["newdata"]][["mcmc"]][["EM"]],
                               predsLong[["newdata2"]][["mcmc"]][["EM"]]
                             )
                           }
  )
  mcmc_pred_surv <- lapply(vector_of_file_names, # predict CIF
                           FUN = function(x) {
                             load(x)
                             predsEvent <- predict(fit_jm,
                                                   newdata = dat_id,
                                                   times = seq(35, 100, length.out = number_of_new_points_for_prediction),
                                                   process = "event", return_mcmc = T,
                                                   n_samples = n_samples, n_mcmc = n_mcmc
                             )
                             rm(fit_jm)
                             list(predsEvent[["mcmc"]], predsEvent[["times"]])
                           }
  )
  # mix all draws together
  draws_long_avail <- do.call(cbind, lapply(mcmc_pred_long, function(x) { # choose predicted memory 
    x[[1]]
  }))
  draws_long_new <- do.call(cbind, lapply(mcmc_pred_long, function(x) {
    x[[2]]
  }))
  draws_surv <- do.call(cbind, lapply(mcmc_pred_surv, function(x) { # choose predicted CIF
    x[[1]]
  }))
  times_new <- do.call(cbind, lapply(mcmc_pred_surv, function(x) {
    x[[2]]
  }))
  level <- 0.95
  
  pooled_long <- data.frame( # combine predictions of the observed and the new longitudinal data
    method = method,
    times = c(dat_id[!is.na(dat_id$EM), ]$age, rowMeans(times_new)),
    pred_EM = rowMeans(rbind(draws_long_avail, draws_long_new)),
    low = apply(rbind(draws_long_avail, draws_long_new),
                1,
                FUN = function(x) {
                  quantile(x, probs = (1 - level) / 2)
                }
    ),
    upp = apply(rbind(draws_long_avail, draws_long_new),
                1,
                FUN = function(x) {
                  quantile(x, probs = (1 + level) / 2)
                }
    )
  )
  pooled_surv <- data.frame( # combine predictions of the observed and the new CIF data
    method = method,
    times = rowMeans(times_new),
    pred_CIF = rowMeans(draws_surv),
    low = apply(draws_surv,
                1,
                FUN = function(x) {
                  quantile(x, probs = (1 - level) / 2)
                }
    ),
    upp = apply(draws_surv,
                1,
                FUN = function(x) {
                  quantile(x, probs = (1 + level) / 2)
                }
    )
  )
  list(pooled_long, pooled_surv)
}
plot_predictions <- function(list_of_pooled_predictions, dat_id, xlim_min) {
  pred_long <- do.call(rbind, lapply(list_of_pooled_predictions[-length(list_of_pooled_predictions)], function(x) { x[[1]]}))
  pred_surv <- do.call(rbind, lapply(list_of_pooled_predictions[-length(list_of_pooled_predictions)], function(x) { x[[2]]}))

  obs_and_pred_long <- pred_long %>%
    full_join(
      data.frame(
        method = "observed",
        times = list_of_pooled_predictions[["dat_id"]]$age,
        EM = list_of_pooled_predictions[["dat_id"]]$EM
      )
    ) %>%
    mutate(
      Age = times,
      method = gsub("delta = ", "", method),
      method = factor(method, levels = c("MAR", "0", "changing", "-1", "observed"),
                      labels=c("MAR", "0", "age-varying", "-1", "observed"))
    )
  
  plot1 <- obs_and_pred_long %>%
    ggplot(aes(group = method)) +
    geom_point(
      data = subset(obs_and_pred_long, method == "observed"),
      aes(x = Age, y = EM, group = method),
    ) +
    geom_ribbon(
      data = subset(obs_and_pred_long, !method == "observed"),
      aes(x = Age, ymin = low, ymax = upp, linetype = method),
      alpha = 0.05, show.legend = F
    ) +
    geom_line(
      data = subset(obs_and_pred_long, !method == "observed"),
      aes(x = Age, y = low,  group = method, linetype = method), 
      linewidth = 0.1, col = "grey", show.legend = F
    ) +
    geom_line(
      data = subset(obs_and_pred_long, !method == "observed"),
      aes(x = Age, y = upp, group = method, linetype = method),
      linewidth = 0.1, col = "grey", show.legend = F
    ) +
    geom_line(
      data = subset(obs_and_pred_long, !method == "observed"),
      aes(x = Age, y = pred_EM, group = method, linetype = method),
      linewidth = 0.3
    ) +
    ylab("Memory") +
    ylim(-30, 60) +
    xlim(xlim_min, 100) + 
    theme_classic() +
    theme(
      legend.position = "bottom",
      plot.margin = margin(0, 0.1, 0,0 , "cm") ,
      # legend.title = element_text(size = 20),
      # legend.text = element_text(size = 20) #,
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 5),
      axis.title.y = element_text(margin = margin(r = 0) ))  + 
    guides(linetype=guide_legend(nrow=1, byrow = TRUE, keywidth = 2)) +
    scale_linetype_manual(
      values = c("MAR" = "solid" ,"0" = "longdash", "age-varying" = "dashed", "-1" = "dotted"),
      labels = c("MAR", expression(paste(Delta, "=", 0)), expression(paste(Delta, "=", -((age - 25) / 75)^3)), expression(paste(Delta, "=-1")))
    )
  
  
  plot2 <- pred_surv %>%
    mutate(Age = times) %>%
    ggplot() +
    geom_line(aes(x = Age, y = pred_CIF, group = method, linetype = method), linewidth = 0.3)+
    geom_line(aes(x = Age, y = low, group = method, linetype = method), linewidth = 0.1, col= "grey", show.legend = F) +
    geom_line(aes(x = Age, y = upp, group = method, linetype = method), linewidth = 0.1, col = "grey", show.legend = F) +
    geom_ribbon(aes(x = Age, ymin = low, ymax = upp, linetype = method), alpha = 0.05, show.legend = F) +
    ylab("CIF") +
    xlab("Age, years") +
    xlim(xlim_min, 100) + 
    theme_classic() +
    theme(
      legend.position = "bottom",
      plot.margin = margin(0, 0.1, 0, 0 , "cm"), 
      # legend.title = element_text(size = 20),
      # legend.text = element_text(size = 20) #,
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 5),
      axis.title.y = element_text(margin = margin(r = 0) ))  + 
    #guides(linetype=guide_legend(nrow=1, byrow = TRUE, keywidth = 2)) +
    scale_linetype_manual(
      values = c("MAR" = "solid" ,"0" = "longdash", "changing" = "dashed", "-1" = "dotted"),
      labels = c("MAR", expression(paste(Delta, "=", 0)), expression(paste(Delta, "=", -((age - 25) / 75)^3)), expression(paste(Delta, "=-1")))
    ) +
    scale_y_continuous(
      breaks = c(0, 0.5, 1),
      limits = c(0,1)
    ) 
  
  list(plot1, plot2)
}

# Predictions for a female with 10 y. of education born 1937 ------------------
# read data
load("dat.Rdata")
mean_dat_long <- dat_long %>%
  group_by(age_cohort_T) %>%
  summarise(EM = mean(EM)) %>%
  mutate(
    educ_imp = 10.5,# mean(dat_surv$educ_imp) = 10.48
    sex = 0,
    ind_enrol = 0,
    age = age_cohort_T,
    id = 6000,
    cohort = 0,
    mean_age = unique(dat_long$mean_age)
  ) %>%
  ungroup() 


#participants_with_two_measures <- unique(participants_with_two_measures )
plot_dynamic_prediction_memory <- list()
plot_dynamic_prediction_survival <- list()
list_of_pooled_predictions <- list()

for(i in 1:4){
  print(i)
  dat_id <- mean_dat_long %>%
    filter(age <= 40+15*(i-1) & age >= 40+15*(i-1)-15) %>%
    mutate(EM_imp = EM)
  names_all_files <- list.files("results/")
    
  # pool predictions from  the imputations
  pool_pred_observed <- pool_pred_obs(
    vector_of_file_names = "results/fit_jm_observed.Rdata",
    dat_id, method = "MAR"
  )
  pool_pred_delta_0 <- pool_pred(
    vector_of_file_names = paste0(
      "results/",
      names_all_files[grep("fit_jm_delta_0", names_all_files)]
    ),
    dat_id, method = "0"
  )
  pool_pred_delta_changing <- pool_pred(
    vector_of_file_names = paste0(
      "results/",
      names_all_files[grep("fit_jm_delta_changing", names_all_files)]
    ),
    dat_id, method = "changing"
  )
  pool_pred_delta_1 <- pool_pred(
    vector_of_file_names = paste0(
      "results/",
      names_all_files[grep("fit_jm_delta__1", names_all_files)]
    ),
    dat_id, method = "-1"
  )

  list_of_pooled_predictions[[i]] = 
    list(observed = pool_pred_observed,
         delta_0 = pool_pred_delta_0, 
         delta_changing = pool_pred_delta_changing,
         delta__1= pool_pred_delta_1,
         dat_id = dat_id)
}  

# plot for mean
for(i in 1:4){
  print(i)
  # plot pooled predictions
  plots_for_participant <-  
    plot_predictions(
      list_of_pooled_predictions =  list_of_pooled_predictions[[i]],
      xlim_min= 25
    )
  
  plot_dynamic_prediction_memory[[i]] <- plots_for_participant [[1]]
  plot_dynamic_prediction_survival[[i]] <- plots_for_participant [[2]]
}

pdf(
  file = "reports/Gorbach_Figure_3.pdf",
  width = 6.5, # The width of the plot in inches
  height = 4, onefile = FALSE
)
ggarrange(plotlist = c(plot_dynamic_prediction_memory, plot_dynamic_prediction_survival),
          ncol = 4, nrow = 2,
          common.legend = T, legend = "bottom"
) 
dev.off()

# Predictions for selected ids ------------------
# read data
load("dat.Rdata")
dat_surv_selected_ids <- dat_surv%>%filter(sex == 0 & educ_imp>=9 & educ_imp <=11 & status == F & age_last_obs > 65 & age_last_obs< 70)
# ids <- c(385, 1064, 2865, 3480, 4300) # dat_surv%>%filter(sex == 0 & educ_imp>=9 & educ_imp <=11 & status == F & age_last_obs > 65 & age_last_obs< 70)
# ids <-  c(1155, 2972, 3050, 2960, 2989, 2591) # dat_surv%>%filter(sex == 0 & educ_imp>=9 & educ_imp <=11 & status == F & age_last_obs > 83 & age_last_obs< 85)
# # ids = c(4041, 4369) - very young
# # selected ids from dat_surv%>%filter(sex == 0 & educ_imp>=9 & educ_imp <=11 & status == F & age_last_obs > 65 & age_last_obs< 70).
# # last obs between 65 and 70. 2865: better memory but only 2 obs; 4300- worse memory 2 obs; 1064: better memory, 5 obs. 
# ids <- c(2865, 4300, 1064) #  for figure 3
# # add older people: selected ids from dat_surv%>%filter(sex == 0 & educ_imp>=9 & educ_imp <=11 & status == F & age_last_obs > 83 & age_last_obs< 85)
# # 1155: only 1 obs better memory, 3050: only 1 obs, worse memory; 2591: better memory, 5 mobservations, 2989 - decreasing memory 5 obs; 2960: low memory, 5 obs.
# ids <- append(ids, c(1155, 3050, 2591, 2989, 2960)) # for figure 3
ids <- c(2865, 4300, 1064, 2591, 2960)
dat_long_selected_ids <- dat_long%>%filter(id %in% ids)
# dat_long_selected_ids <- dat_long%>%filter(id %in% a$id)

dat_long_selected_ids%>% # plot observed memory for the selected ids
  filter(!is.na(EM))%>%
  ggplot()+
  geom_line(aes(x = age, y = EM)) +
  facet_wrap(id~.)

#participants_with_two_measures <- unique(participants_with_two_measures )

list_of_pooled_predictions <- list()

for(i in 1:length(ids)){
  print(i)
  dat_id <- dat_long%>%filter(id == ids[i] & !is.na(EM))%>%
    mutate(EM_imp = EM)
  names_all_files <- list.files("results/")
  
  # pool predictions from  the imputations
  pool_pred_observed <- pool_pred_obs(
    vector_of_file_names = "results/fit_jm_observed.Rdata",
    dat_id, method = "MAR"
  )
  pool_pred_delta_0 <- pool_pred(
    vector_of_file_names = paste0(
      "results/",
      names_all_files[grep("fit_jm_delta_0", names_all_files)]
    ),
    dat_id, method = "0"
  )
  pool_pred_delta_changing <- pool_pred(
    vector_of_file_names = paste0(
      "results/",
      names_all_files[grep("fit_jm_delta_changing", names_all_files)]
    ),
    dat_id, method = "changing"
  )
  pool_pred_delta_1 <- pool_pred(
    vector_of_file_names = paste0(
      "results/",
      names_all_files[grep("fit_jm_delta__1", names_all_files)]
    ),
    dat_id, method = "-1"
  )
  list_of_pooled_predictions[[i]] = 
    list(observed = pool_pred_observed,
         delta_0 = pool_pred_delta_0, 
         delta_changing = pool_pred_delta_changing,
         delta__1= pool_pred_delta_1,
         dat_id = dat_id)
}  
# plot
i <- 0
plot_dynamic_prediction_memory <- list()
plot_dynamic_prediction_survival <- list()
for(k in c(1:5)){ # to many subjects, select a few
  i <- i+1
  # plot pooled predictions
  plots_for_participant <-  
    plot_predictions(
      list_of_pooled_predictions =  list_of_pooled_predictions[[k]],
      xlim_min = 40
    )
  
  plot_dynamic_prediction_memory[[i]] <- plots_for_participant [[1]]
  plot_dynamic_prediction_survival[[i]] <- plots_for_participant [[2]]
}

pdf(
  file = "reports/Gorbach_Figure_C3.pdf",
  width = 7, # The width of the plot in inches
  height = 4, onefile=FALSE
)
ggarrange(plotlist = c(plot_dynamic_prediction_memory, plot_dynamic_prediction_survival),
          ncol = 5, nrow = 2,
          common.legend = T, legend = "bottom"
)
dev.off()
