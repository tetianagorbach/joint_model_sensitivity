delta <- -1
seed <- 374986
file_data <- "dat.Rdata"
file_predictions <-  "results/dat_long_with_posterior.Rdata"
output_file_name <-  gsub(":| |-", "_", paste0("results/fit_jm_delta_",  delta, "_", Sys.time()))
n_burnin <-  20000
n_iter <-  50000
n_thin <-  30
number_of_multiple_imputations <- 1
