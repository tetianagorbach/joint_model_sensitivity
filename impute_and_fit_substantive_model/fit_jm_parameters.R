delta <- -0.5
seed <-  1166277
file_data <- "code/dat.Rdata"
file_predictions <-  "code/results/dat_predicted_raw.Rdata"
output_file_name <-  gsub(":| |-", "_", paste0("code/results/fit_jm_delta",  delta, "_", Sys.time()))
n_burnin <-  100000
n_iter <-  200000
n_thin <-  100
number_of_multiple_imputations <- 5
