# read the observed and predicted data ------------------------------------
load("code/dat.Rdata") # load observed data
load("code/results/dat_predicted_raw.Rdata") # load predictions for the missing data


# impute ------------------------------------------------------------------
dat_posterior_mean <- dat_predicted %>%
  mutate(posterior_mean = rowMeans(dplyr::select(dat_predicted, -c(id, age_sc20)))) %>%
  dplyr::select(id, age_sc20, posterior_mean)


dat_long_with_imputed <- dat_long %>%
  left_join(dat_posterior_mean, by = c("id", "age_sc20")) %>%
  group_by(id) %>%
  mutate(age_dropout_sc20 = max(age_sc20[!is.na(EM_sc)])) %>% # calculate  scaled age at dropout
  ungroup() %>%
  mutate(
    t_minus_t_dropout_pos = (age_sc20 - age_dropout_sc20) * (age_sc20 - age_dropout_sc20 > 0),
    imputed0 = posterior_mean - 0 * t_minus_t_dropout_pos,
    imputed0.5 = posterior_mean - 0.5 * t_minus_t_dropout_pos,
    imputed1 = posterior_mean - 1 * t_minus_t_dropout_pos,
    imputed2 = posterior_mean - 2 * t_minus_t_dropout_pos,
    imputed5 = posterior_mean - 5 * t_minus_t_dropout_pos,
    imputed10 = posterior_mean - 10 * t_minus_t_dropout_pos
  ) # %>% # imputations are named with numbers.
# !! tells R to treat the enclosed expression as an actual variable named k, rather than as a symbol k.
# filter(!is.na(imputed0)) # 7 obs are not imputed since they are after the biggest last_obs_time.


# Figure: observed and imputed by missing data pattern --------------------
ids <- intersect(
  unique(dat_long[is.na(dat_long$EM_sc), ] %>% pull(id)), # select ids of people with more than one one observed memory score
  unique(dat_long[!is.na(dat_long$EM_sc), ] %>%
    group_by(id) %>%
    filter(n() > 1) %>%
    ungroup() %>%
    pull(id))
)

# Figure C2 ---------------------------------------------------------------
dat_long_with_imputed <- dat_long_with_imputed %>%
  pivot_longer(
    cols = c(
      imputed0, imputed0.5,
      imputed1, imputed2, imputed5, imputed10
    ),
    names_to = "delta", names_prefix = "imputed",
    values_to = "imputed"
  ) %>%
  mutate(
    delta = factor(-as.numeric(delta), levels = c("0", "-0.5", "-1", "-2", "-5", "-10")),
    imputed = (ifelse(imputed * EM_cc_sd + EM_cc_mean < 0, -EM_cc_mean/EM_cc_sd, imputed))
  )
dat_to_plot <- dat_long_with_imputed %>%
  filter(id %in% ids[1:6])
pdf(
  file = "code/reports/Figure_obs_imputed_per_subject.pdf",
  width = 6, # The width of the plot in inches
  height = 3
)

dat_to_plot %>%
  ggplot() +
  geom_point(
    data = subset(dat_to_plot, is.na(imputed)), mapping =
      aes(x = age_sc20 * 100 + 20, y = EM_sc * 11.23567 + 35.0212), shape = 16, cex = 1
  ) +
  geom_point(
    data = subset(dat_to_plot, !is.na(imputed)), mapping =
      aes(x = age_sc20 * 100 + 20, y = imputed * 11.23567 + 35.0212, shape = delta), cex = 1
  ) +
  facet_wrap(facets = vars(id), scales = "free_y") +
  xlab("Age") +
  ylab("Memory") +
  scale_shape_manual(
    values = c("0" = 1, "-0.5" = 2, "-1" = 3, "-2" = 4, "-5" = 5, "-10" = 6),
    name = expression(paste(Delta, " ="))
  ) +
  theme_classic()
dev.off()


# Figure 2 ----------------------------------------------------------------
# dat_long_with_imputed %>%
#   ggplot() +
#   geom_point(
#     data = subset(dat_long_with_imputed, delta %in% c(0, -10) & is.na(imputed)), mapping =
#       aes(x = age_sc20 * 100 + 20, y = EM_sc * EM_cc_sd + EM_cc_mean), shape = 16, cex = 1
#   ) +
#   geom_point(
#     data = subset(dat_long_with_imputed, delta %in% c(0, -10) & !is.na(imputed)), mapping =
#       aes(x = age_sc20 * 100 + 20, y = imputed * EM_cc_sd + EM_cc_mean, col = delta), cex = 1
#   ) +
#   facet_wrap(facets = vars(missing_data_pattern), scales = "free_y") +
#   xlab("Age") +
#   ylab("Memory") +
#   scale_color_manual(
#     values = c("0" = 10, "-0.5" = 2, "-1" = 3, "-2" = 4, "-5" = 5, "-10" = 6),
#     name = expression(paste(Delta, "="))
#   ) +
#   theme_classic()





pdf(
  file = "code/reports/Figure_predictions_per_pattern.pdf",
  width = 6, # The width of the plot in inches
  height = 3
)
dat_to_plot_imp <- dat_long_with_imputed %>%
  filter(delta %in% c(0, -10)) %>%
  mutate(
    EM_original_scale = ifelse(!is.na(EM_sc), EM_sc * EM_cc_sd + EM_cc_mean, imputed * EM_cc_sd + EM_cc_mean),
    age = age_sc20 * 100 + 20
  )
dat_to_plot_imp %>%
  ggplot() +
  geom_point(
    data = subset(dat_to_plot_imp, !is.na(EM_sc)),
    mapping = aes(x = age, y = EM_original_scale), shape = 16, cex = 0.1
  ) +
  geom_point(
    data = subset(dat_to_plot_imp, !is.na(imputed)),
    mapping = aes(x = age, y = EM_original_scale, col = delta), cex = 0.01
  ) +
  geom_smooth(
    mapping = aes(x = age, y = EM_original_scale, linetype = delta),
    method = "loess", col = "black", alpha = 0.1, linewidth = 0.5
  ) +
  facet_wrap(facets = vars(missing_data_pattern)) +
  xlab("Age") +
  ylab("Memory") +
  scale_linetype_manual(
    values = c("0" = 1, "-10" = 2),
    name = expression(paste(Delta, "="))
  ) +
  scale_color_manual(
    values = c("0" = "#0C7BDC", "-10" = "#FFC20A"),
    name = expression(paste(Delta, "="))
  ) +
  guides(colour = guide_legend(override.aes = list(size = 1))) +
  theme_classic()
dev.off()
