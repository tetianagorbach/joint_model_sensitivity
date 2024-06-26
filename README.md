read_data.r. reads and preprocesses the Betula data


Folder "imputation" contains 
	fit_jm_pattern_mixture.r  with the code to fit the imputation model to data from dat.Rdata  (not supplied).
	posterior. r with the code to calculate the posterior for  the missing longitudinal observations in dat.Rdata from the imputation models and save the predicted values in dat_long_with_posterior.Rdata
	
Folder "substantive_model" contains 
	fit_jm_observed.r - the code to fit the substantive model to the observed data.
	fit_jm_parameters.r initialises the values of the parameters to fit the substantive model. 
	seeds.rtf provides the seeds used to fit the substantive models parametrised by the sensitivity parameter delta.
	fit_jm_delta_when_dropout_after_60 - the code for multiple imputations: impute missing longitudinal observations  and fit the substantive model to the completed data sets.
	fit_jm_delta_changes_with_time  -  the code for multiple imputations: impute missing longitudinal observations under changing delta  and fit the substantive model to the completed data sets.
	
reports contains R-files to construct figures and tables in the paper. 
	descriptives_lme.Rmd -  the code for  descriptives and fits models for longitudinal data.
	iluutsrations_imputations.R - the code for the figures to describe imputed values.
	pooled_estimates.Rmd  - the code to analyse substantive models  after multiple imputations.
	dynamic_predictions.R - provides the code for the figures with dynamic predictions.