library(tidyverse)
# read the observed and predicted data ------------------------------------
load("dat.Rdata") # load observed data
load("results/dat_long_with_posterior.Rdata") # load predictions for the missing data

# impute ------------------------------------------------------------------
dat_long_with_posterior <- dat_long_with_posterior %>%
  mutate(EM_posterior = rowMeans(dplyr::select(dat_long_with_posterior, starts_with("post")))) %>%
  dplyr::select(-starts_with("post"))


dat_long_with_imputed <- dat_long_with_posterior %>%
  group_by(id) %>%
  mutate(age_dropout = max(age[!is.na(EM)])) %>% # calculate  scaled age at dropout
  ungroup() %>%
  mutate(
    t_minus_t_dropout_pos = (age - age_dropout) * (age - age_dropout > 0) * (age_dropout > 60),
    dropout_after_60 = age_dropout > 60,
    imputed_0 = EM_posterior,
    imputed_1 = EM_posterior + (-1) * t_minus_t_dropout_pos,
    imputed_changing = EM_posterior - 1 * ((age - 25) / 75)^3 * t_minus_t_dropout_pos
  ) # %>% # imputations are named with numbers.
# !! tells R to treat the enclosed expression as an actual variable named k, rather than as a symbol k.
# filter(!is.na(imputed0)) # 7 obs are not imputed since they are after the biggest last_obs_time.
# d <- dat_long_with_imputed%>%filter(id == 1570)
dat_long_with_imputed <- dat_long_with_imputed %>%
  pivot_longer(
    cols = c(
      imputed_0, imputed_changing, imputed_1
    ),
    names_to = "Delta", names_prefix = "a",
    values_to = "imputed"
  ) %>%
  mutate(
    # delta = factor(-as.numeric(delta), levels = c("0", "-0.5", "-1", "-2", "-5", "-10")),
    imputed = (ifelse(imputed < 0, 0, imputed))
  ) %>%
  mutate(
    Delta = gsub("imputed_", "", Delta),
    Delta = ifelse(Delta == 1, -1, Delta),
    Delta = factor(Delta, levels = c("0", "changing", "-1"), labels = c("0", "changing", "-1"))
  ) %>%
  mutate(imputed = ifelse(!is.na(EM), EM, imputed))



# Figure: observed and imputed by missing data pattern --------------------
# plot(x = 25:100, y = -1*((25:100-25)/75)^3, type = "l",
#      xlab = "Age, years", ylab = expression(Delta))


# Figure 2 imputations and mean trajectories per pattern ------------------

imputations_per_pattern <- dat_long_with_imputed %>%
  ggplot() +
  geom_point(
    data = subset(dat_long_with_imputed, !is.na(EM)),
    mapping = aes(x = age, y = EM), shape = 16, cex = 0.5
  ) +
  geom_point(
    data = subset(dat_long_with_imputed, is.na(EM)),
    mapping = aes(x = age, y = imputed, shape = Delta, col = Delta), cex = 0.1
  ) +
  geom_smooth(
    mapping = aes(x = age, y = imputed, linetype = Delta),
    method = "loess", col = "black", alpha = 0.2, linewidth = 0.2
  ) +
  facet_wrap(facets = vars(missing_data_pattern)) +
  xlab("Age, years") +
  ylab("Memory") +
  scale_linetype_manual(
    values = c("0" = "longdash", "changing" = "dashed", "-1" = "dotted"),
    name = expression(paste(Delta, "="))
  ) +
  scale_color_manual(
    values = c("0" = "yellow", "changing" = "blue", "-1" = "#D55E00"),
    name = expression(paste(Delta, "="))
  ) +
  scale_shape_manual(
    values = c("0" = 1, "changing" = 3, "-1" = 2),
    name = expression(paste(Delta, "="))
  ) +
  guides(colour = guide_legend(override.aes = list(size = 1))) +
  theme_classic()
pdf(
  file = "reports/Gorbach_Figure_2.pdf",
  width = 6, # The width of the plot in inches
  height = 3
)
imputations_per_pattern
dev.off()


# Figure C1:imputations for specific people -------------------------------
plot_individuals <-
  ggplot() +
  geom_point(
    data = subset(dat_long_with_imputed, !is.na(EM) & id %in% c(6, 28, 36, 41, 42, 43)), mapping =
      aes(x = age, y = EM), shape = 16, cex = 1
  ) +
  geom_point(
    data = subset(dat_long_with_imputed, is.na(EM) & id %in% c(6, 28, 36, 41, 42, 43)), mapping =
      aes(x = age, y = imputed, shape = Delta), cex = 1
  ) +
  facet_wrap(facets = vars(id)) +
  xlab("Age, years") +
  ylab("Memory") +
  scale_shape_manual(
    values = c("0" = 1, "changing" = 3, "-1" = 2),
    name = expression(paste(Delta, "="))
  ) +
  theme_classic()
pdf(
  file = "reports/Gorbach_Figure_C1.pdf",
  width = 6, # The width of the plot in inches
  height = 3
)
plot_individuals
dev.off()


