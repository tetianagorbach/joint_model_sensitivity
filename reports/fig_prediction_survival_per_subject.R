library(tidyverse)
library(survival)
library(ggpubr)


# custom functions --------------------------------------------------------
calc_mean_assoc_at_65 <- function(dat_for_em_prediction, vector_of_file_names) {
  em_pred <- data.frame(id = factor(dat_for_em_prediction$id))
  K <- length(vector_of_file_names)

  for (k in 1:length(vector_of_file_names)) {
    load(vector_of_file_names[k])
    col_list <- paste0("EM_sc_pred", k)
    par_est <- fit_jm[["statistics"]][["Mean"]][["betas1"]]
    est_re <- data.frame(
      id = fit_jm[["model_data"]][["idT"]],
      b = fit_jm[["statistics"]][["Mean"]][["b"]]
    )
    dat <- dat_for_em_prediction %>%
      left_join(est_re, by = "id") %>%
      mutate(
        EM_sc_pred = par_est[["(Intercept)"]] +
          par_est[["poly(age_sc20 - mean_age_sc20, 3, raw = TRUE)1"]] * (age_sc20 - mean_age_sc20) +
          par_est[["poly(age_sc20 - mean_age_sc20, 3, raw = TRUE)2"]] * (age_sc20 - mean_age_sc20)^2 +
          par_est[["poly(age_sc20 - mean_age_sc20, 3, raw = TRUE)3"]] * (age_sc20 - mean_age_sc20)^3 +
          par_est[["sex"]] * sex +
          par_est[["educ_imp_sc"]] * educ_imp_sc +
          par_est[["ind_enrolTRUE"]] * ind_enrol +
          par_est[["I(ind_enrol * (age_sc20 - mean_age_sc20))"]] * (age_sc20 - mean_age_sc20) * ind_enrol +
          par_est[["cohort_sc"]] * cohort_sc +
          b.1 +
          b.2 * (age_sc20 - mean_age_sc20),
        EM_sc_pred_der =
          par_est[["poly(age_sc20 - mean_age_sc20, 3, raw = TRUE)1"]] +
            par_est[["poly(age_sc20 - mean_age_sc20, 3, raw = TRUE)2"]] * 2 * (age_sc20 - mean_age_sc20) +
            par_est[["poly(age_sc20 - mean_age_sc20, 3, raw = TRUE)3"]] * 3 * (age_sc20 - mean_age_sc20)^2 +
            b.2,
        assoc = fit_jm[["statistics"]][["Mean"]][["alphas"]][["value(EM_sc_imp)"]] * EM_sc_pred +
          fit_jm[["statistics"]][["Mean"]][["alphas"]][["slope(EM_sc_imp)"]] * EM_sc_pred_der
      )
    em_pred[paste0("assoc", k)] <- dat$assoc
    rm(est_re)
  }
  em_pred <- em_pred %>%
    mutate(assoc_pred = rowMeans(select(em_pred, starts_with("assoc"))))
}
pool_pred_long <- function(dat_for_em_prediction, jm_fit_mi) {
  em_pred <- data.frame(id = dat_for_em_prediction$id)
  K <- length(jm_fit_mi)
  for (k in 1:K) {
    col_list <- paste0("EM_sc_pred", k)
    par_est <- fit_jm[["statistics"]][["Mean"]][["betas1"]]
    est_re <- data.frame(
      id = fit_jm[["model_data"]][["idT"]],
      b = fit_jm[["statistics"]][["Mean"]][["b"]]
    ) %>%
      mutate(id = as.numeric(id))
    dat <- dat_for_em_prediction %>%
      left_join(est_re, by = "id") %>%
      mutate(EM_sc_pred = par_est[["(Intercept)"]] +
        par_est[["sex"]] * sex +
        par_est[["educ_imp_sc"]] * educ_imp_sc +
        par_est[["age_enroll_sc20"]] * age_enroll_sc20 +
        par_est[["ind_enrollmentTRUE"]] * ind_enrollment +
        par_est[["poly(age_sc20 - 0.42, 3, raw = TRUE)1"]] * (age_sc20 - mean_age_sc20) +
        par_est[["poly(age_sc20 - 0.42, 3, raw = TRUE)2"]] * (age_sc20 - mean_age_sc20)^2 +
        par_est[["poly(age_sc20 - 0.42, 3, raw = TRUE)3"]] * (age_sc20 - mean_age_sc20)^3 +
        par_est[["age_enroll_sc20:ind_enrollmentTRUE"]] * age_enroll_sc20 * ind_enrollment +
        b.1 +
        b.2 * (age_sc20 - mean_age_sc20))
    em_pred[paste0("EM_sc_pred", k)] <- dat$EM_sc_pred
    rm(est_re)
  }
  em_pred <- em_pred %>%
    mutate(EM_sc_pred = rowMeans(select(em_pred, starts_with("EM_sc_pred"))))
  em_pred$EM_sc_pred
}

# read data and  JM estimation fit for multiple imputed data sets ------------------
load("code/dat.Rdata")


dat_surv <- dat_surv %>% # add scaled education to the survival data
  left_join(dat_long %>%
    distinct(id, cohort_sc, mean_age_sc20), by = "id")
dat_for_em_prediction <- dat_surv %>%
  dplyr::select(id, sex, educ_imp_sc, cohort_sc, mean_age_sc20) %>%
  filter(!is.na(educ_imp_sc)) %>%
  mutate(
    id = factor(id, levels = unique(id)),
    ind_enrol = 0,
    age = 65,
    across(c("age"), ~ (.x - 20) / 100, .names = "{.col}_sc20")
  ) # scaled age 0 corresponds to true age 20)

names_all_files <- list.files("code/results")

# predict the link between hazards and survival at age 65 for each imputation when delta = 0
em_pred_delta0 <- calc_mean_assoc_at_65(dat_for_em_prediction,
  vector_of_file_names = paste0("code/results/", names_all_files[grep("fit_jm_delta0_", names_all_files)])
)
# predict the link between hazards and survival at age 65 for each imputation when delta = -5
em_pred_delta5 <- calc_mean_assoc_at_65(dat_for_em_prediction,
  vector_of_file_names = paste0("code/results/", names_all_files[grep("fit_jm_delta_10_", names_all_files)])
)

# choose  people in the lower 30%  in both method (they should have better survival)
id_bottom <- intersect(
  em_pred_delta0 %>% # select people  when delta = 0
    filter(assoc_pred < quantile(assoc_pred, probs = 1 / 3)) %>%
    pull(id),
  em_pred_delta5 %>% # select people  when delta = -1.6
    filter(assoc_pred < quantile(assoc_pred, probs = 1 / 3)) %>%
    pull(id)
)

# choose  people in the middle 30%  in both
id_middle <- intersect(
  em_pred_delta0 %>% # select people  when delta = 0
    filter(assoc_pred > quantile(assoc_pred, probs = 1 / 3) & assoc_pred < quantile(assoc_pred, probs = 2 / 3)) %>%
    pull(id),
  em_pred_delta5 %>% # select people  when delta = -1.6
    filter(assoc_pred > quantile(assoc_pred, probs = 1 / 3) & assoc_pred < quantile(assoc_pred, probs = 2 / 3)) %>%
    pull(id)
)
# choose people in the higher 30% of association measure in both
id_top <- intersect(
  em_pred_delta0 %>% # select people  when delta = 0
    filter(assoc_pred > quantile(assoc_pred, probs = 2 / 3)) %>%
    pull(id),
  em_pred_delta5 %>% # select people  when delta = -1.6
    filter(assoc_pred > quantile(assoc_pred, probs = 2 / 3)) %>%
    pull(id)
)
ids_more_than_1_obs <- dat_long %>%
  filter(id %in% dat_long[is.na(EM_sc), ]$id & sex == 0 & abs(educ_imp_sc) < 0.5) %>%
  group_by(id) %>%
  summarise(
    n = sum(!is.na(EM_sc)),
    max_age = max(age[!is.na(EM_sc)])
  ) %>%
  ungroup() %>%
  filter(max_age > 70 & n > 1 & max_age < 80) %>%
  pull(id)

id_bottom <- intersect(
  intersect(id_bottom, dat_surv %>% filter(status == 0) %>% pull(id)),
  ids_more_than_1_obs
) # in bottom 30% but not demented
id_middle <- intersect(
  intersect(id_middle, dat_surv %>% filter(status == 0) %>% pull(id)),
  ids_more_than_1_obs
)
id_top <- intersect(
  intersect(id_top, dat_surv %>% filter(status == 0) %>% pull(id)),
  ids_more_than_1_obs
)

# plot pooled predictions and CI after 65
# plot pooled predictions and CI after 65
set.seed(11)
ids <- c(
  sample(id_bottom, 1),
  sample(id_middle, 1),
  sample(id_top, 1)
)

# rm(ids_more_than_3_obs, calc_mean_assoc_at_65, pool_pred_long, em_pred_delta0, em_pred_delta5)
n <- 50
plots <- list()
for (i in 1:length(ids)) {
  dat_id <- dat_long %>%
    filter(id == ids[i]) %>%
    mutate(across(c("age"), ~ (.x - 20) / 100, .names = "{.col}_sc20"),
      EM_sc_imp = EM_sc
    )
  time_death_dementia <- min((dat_surv %>% filter(id == ids[i]) %>% pull(age_dementia) - 20) / 100,
    (dat_surv %>% filter(id == ids[i]) %>% pull(age_death) - 20) / 100,
    na.rm = T
  )
  pool_pred <- function(vector_of_file_names, dat_id, method) {
    mcmc_pred_long <- lapply(vector_of_file_names,
      FUN = function(x) {
        load(x)
        predsLong <- predict(fit_jm,
          newdata = dat_id,
          times = seq(max(dat_id$age_sc20), min(0.8, time_death_dementia, na.rm = T), length.out = n),
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
          times = seq(max(dat_id$age_sc20), min(0.8, time_death_dementia, na.rm = T), length.out = n),
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
      times = c(dat_id[!is.na(dat_id$EM_sc), ]$age_sc20, rowMeans(times_new)),
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

  plot_pred <- function(pool_pred_0, pool_pred_m5) {
    pred_long <- rbind(pool_pred_0[[1]], pool_pred_m5[[1]])
    pred_surv <- rbind(pool_pred_0[[2]], pool_pred_m5[[2]])
    # }
    obs_and_pred_long <- pred_long %>%
      full_join(
        data.frame(
          method = "observed",
          times = dat_id$age_sc20,
          EM = dat_id$EM_sc
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
      scale_linetype_manual(
        values = c("0" = 1, "-10" = 2),
        name = bquote(Delta == .(value))
      ) +
      ylab("Memory") +
      xlim(50, 100) +
      ylim(-20, 60) +
      theme_classic() +
      theme(
        legend.position = "bottom",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20)
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
      scale_linetype_manual(
        values = c("0" = 1, "-10" = 2),
        name = bquote(Delta == .(value))
      ) +
      xlim(50, 100) +
      ylab("CIF") +
      xlab("Age") +
      ylim(0, 1) +
      theme(
        legend.position = "bottom",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20)
      )

    list(plot1, plot2)
  }

  plots <- append(
    plots,
    plot_pred(
      pool_pred_0 = pool_pred(
        vector_of_file_names = paste0("code/results/", names_all_files[grep("fit_jm_delta0_", names_all_files)]),
        dat_id, method = "0"
      ),
      pool_pred_m5 = pool_pred(
        vector_of_file_names = paste0("code/results/", names_all_files[grep("fit_jm_delta_10_", names_all_files)]),
        dat_id, method = "-10"
      )
    )
  )
}



pdf(
  file = "code/reports/Figure_survival_per_subject.pdf", # The directory you want to save the file in
  width = 12, # The width of the plot in inches
  height = 7,
  onefile = F
)
ggarrange(plots[[1]], plots[[3]],
  plots[[5]], plots[[2]],
  plots[[4]], plots[[6]],
  ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom"
)

dev.off()
