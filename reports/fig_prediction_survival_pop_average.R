library(tidyverse)
library(survival)
library(ggpubr)



# functions ---------------------------------------------------------------
pool_pred <- function(vector_of_file_names, dat_id, method) {
  mcmc_pred_long <- lapply(vector_of_file_names,
    FUN = function(x) {
      load(x)
      predsLong <- predict(fit_jm,
        newdata = dat_id,
        times = seq(0.4, 0.8, length.out = n),
        process = "longitudinal", return_mcmc = T,
        type = "subject_specific"
      )
      rm(fit_jm)
      list(
        predsLong[["newdata"]][["mcmc"]][["EM_sc_imp"]],
        predsLong[["newdata2"]][["mcmc"]][["EM_sc_imp"]]
      )
    }
  )
  mcmc_pred_surv <- lapply(vector_of_file_names,
    FUN = function(x) {
      load(x)
      predsEvent <- predict(fit_jm,
        newdata = dat_id,
        times = seq(0.4, 0.8, length.out = n),
        process = "event", return_mcmc = T
      )
      rm(fit_jm)
      list(predsEvent[["mcmc"]], predsEvent[["times"]])
    }
  )
  # mix all draws together
  draws_long_avail <- do.call(cbind, lapply(mcmc_pred_long, function(x) { # combine predicted memory for all imputations
    x[[1]]
  }))
  draws_long_new <- do.call(cbind, lapply(mcmc_pred_long, function(x) {
    x[[2]]
  }))
  draws_surv <- do.call(cbind, lapply(mcmc_pred_surv, function(x) { # combine survival from all imputations
    x[[1]]
  }))
  times_new <- do.call(cbind, lapply(mcmc_pred_surv, function(x) {
    x[[2]]
  }))
  level <- 0.95

  pooled_long <- data.frame(
    method = method,
    times = c(dat_id[!is.na(dat_id$EM_sc_imp), ]$age_sc20, rowMeans(times_new)),
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
  pooled_surv <- data.frame(
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
plot_pred <- function(pool_pred_1, pool_pred_2) {
  pred_long <- rbind(pool_pred_1[[1]], pool_pred_2[[1]])
  pred_surv <- rbind(pool_pred_1[[2]], pool_pred_2[[2]])
  # }
  obs_and_pred_long <- pred_long %>%
    full_join(
      data.frame(
        method = "observed",
        times = dat_id$age_sc20,
        EM = dat_id$EM_sc_imp
      )
    ) %>%
    mutate(
      Age = times * 100 + 20,
      method = gsub("delta = ", "", method)
    )
  require(gridExtra)
  EM_cc_mean <- unique(dat_long$EM_cc_mean)
  EM_cc_sd <- unique(dat_long$EM_cc_sd)

  plot1 <- obs_and_pred_long %>%
    ggplot(aes(group = method)) +
    geom_point(
      data = subset(obs_and_pred_long, method == "observed"),
      aes(x = Age, y = EM * EM_cc_sd + EM_cc_mean, group = method)
    ) +
    geom_line(
      data = subset(obs_and_pred_long, !method == "observed"),
      aes(x = Age, y = pred_EM * EM_cc_sd + EM_cc_mean, group = method, linetype = method)
    ) +
    geom_ribbon(
      data = subset(obs_and_pred_long, !method == "observed"),
      aes(x = Age, ymin = low * EM_cc_sd + EM_cc_mean, ymax = upp * EM_cc_sd + EM_cc_mean, linetype = method), alpha = 0.1, lty = 1
    ) +
    geom_line(
      data = subset(obs_and_pred_long, !method == "observed"),
      aes(x = Age, y = low * EM_cc_sd + EM_cc_mean, group = method, linetype = method), linewidth = 0.1
    ) +
    geom_line(
      data = subset(obs_and_pred_long, !method == "observed"),
      aes(x = Age, y = upp * EM_cc_sd + EM_cc_mean, group = method, linetype = method), linewidth = 0.1
    ) +
    # scale_linetype_manual(
    #   values = c("0" = 1, "-10" = 2),
    #   name = bquote(Delta == .(value))
    # ) +
    ylab("Memory") +
    theme_classic() +
    theme(
      legend.position = "bottom" # ,
      # legend.title = element_text(size = 20),
      # legend.text = element_text(size = 20) #,
      # axis.title = element_text(size = 20),
      # axis.text = element_text(size = 20)
    ) # +
  # ggplot(obs_and_pred_long%>%filter(method == "observed")%>%
  #          mutate(Age = times*100 + 20),
  #        aes(x = Age, y = pred_EM)) +
  # geom_point()

  plot2 <- pred_surv %>%
    mutate(Age = times * 100 + 20) %>%
    ggplot(aes(x = Age, y = pred_CIF, group = method, linetype = method)) +
    geom_line() +
    geom_line(aes(x = Age, y = low, group = method, linetype = method), linewidth = 0.1) +
    geom_line(aes(x = Age, y = upp, group = method, linetype = method), linewidth = 0.1) +
    geom_ribbon(aes(ymin = low, ymax = upp, linetype = method), alpha = 0.1, lty = 1) +
    theme_classic() +
    # scale_linetype_manual(
    #   values = c("0" = 1, "-10" = 2),
    #   name = bquote(Delta == .(value))
    # ) +
    ylab("CIF") +
    xlab("Age") +
    ylim(0, 1) +
    theme(
       legend.position = "bottom" # ,
      # legend.title = element_text(size = 20),
      # legend.text = element_text(size = 20) #,
      # axis.title = element_text(size = 20),
      # axis.text = element_text(size = 20)
    )

  list(plot1, plot2)
}

# analysis ------------------
# read data
load("code/dat.Rdata")

# calculate mean profile
dat_id <- dat_long %>%
  group_by(age_cohort_T) %>%
  summarise(EM_sc_imp = mean(EM_sc)) %>%
  mutate(
    educ_imp_sc = 0,
    sex = 0,
    ind_enrol = 0,
    age = age_cohort_T,
    id = 6000,
    cohort_sc = 0,
    mean_age_sc20 = unique(dat_long$mean_age_sc20),
    EM_sc_imp = 0
  ) %>%
  ungroup() %>%
  filter(age <= 65) %>%
  mutate(age_sc20 = (age - 20) / 100)

names_all_files <- list.files("code/results")
set.seed(11)
n <- 50

# pool predictions from  the imputations
pool_pred_1 <- pool_pred(
  vector_of_file_names = paste0(
    "code/results/",
    names_all_files[grep("fit_jm_delta0_", names_all_files)]
  ),
  dat_id, method = "0"
)

pool_pred_2 <- pool_pred(
  vector_of_file_names = paste0(
    "code/results/",
    names_all_files[grep("fit_jm_delta_2_", names_all_files)]
  ),
  dat_id, method = "-10"
)

# plot pooled predictions
plots <- plot_pred(
  pool_pred_1,
  pool_pred_2
)

# pdf(
#   file = "code/reports/Figure_survival_pop_average_prediction.pdf", # The directory you want to save the file in
#   width = 6, # The width of the plot in inches
#   height = 3,
#   onefile = F
# )
ggarrange(plots[[1]], plots[[2]],
  ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom"
)

dev.off()
