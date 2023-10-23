read_data.r. read the initial data


fit_pattern_mixture_imputation_model contains 
	fit_jm_pattern_mixture.r  with the code to fit the imputation model to data from dat.Rdata  (not supplied) 
	predict. r with the code to predict the missing longitudinal observations in dat.Rdata from the imputation models and save the predicted values in dat_predicted_raw.Rdata
	
impute_and_fit_substantive_model contains 
	fit_jm_observed.r - the code to fit the substantive model to the observed data
	fit_jm_parameters.r initialises the values of the parameters to fit the substantive model. 
	seeds. provides the seeds used to fit the substantive models parametrised by the sensiitivity parameter delta.
	fit_jm - the code for multiple imputations: impute missing longitudinal observations using  dat.Rdata and dat_predicted_raw.Rdata and fit the substantive model to the completed data sets.
	
reports contains R-files to construct figures and tables in the paper. 
	descriptives_lme.Rmd -  the code for  descriptives and fits models for longitudinal data.
	imputations_report.R - the code for the figures to describe imputed values
	substative_model_after_imputation.Rmd  - the code to analyse substantive models  after multiple imputations
	fig_prediction_survival_pop_average.R - provides the code for the figure with dynamic predictions.