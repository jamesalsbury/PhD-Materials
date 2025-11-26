library(DTEAssurance)

# Trial A ----------------------------

n_c <- 276
n_t <- 275

censoring_model <- list(method = "Events", events = 291)

recruitment_model <- list(method = "power",
                          period = 31.5,
                          power = 1)

analysis_model <- list(method = "LRT",
                       alpha = 0.025,
                       alternative_hypothesis = "one.sided")


## Nominal Power ----------------------

target_HR <- 0.75

control_model <- list(dist = "Exponential",
                      parameter_model = "Fixed",
                      fixed_type = "Landmark",
                      t1 = 14,
                      surv_t1 = 0.5)

effect_model <- list(delay_SHELF = SHELF::fitdist(c(3, 4, 5),
                                                  probs = c(0.25, 0.5, 0.75), lower = 0, upper = 10),
                     delay_dist = "gamma",
                     HR_SHELF = SHELF::fitdist(c(target_HR - 0.001, target_HR, target_HR + 0.001),
                                               probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1.5),
                     HR_dist = "gamma",
                     P_S = 1, P_DTE = 0)




calc_dte_assurance(n_c =  n_c,
                   n_t = n_t,
                   control_model = control_model,
                   effect_model = effect_model,
                   censoring_model = censoring_model,
                   recruitment_model = recruitment_model,
                   analysis_model = analysis_model,
                   n_sims = 10000)


## Realised Power - Best Fitting ----------------------

observed_delay <- 5.3
observed_lambda <- 0.0420
observed_shape <- 0.923
observed_HR <- 0.794

control_model <- list(dist = "Weibull",
                      parameter_model = "Fixed",
                      fixed_type = "Parameters",
                      lambda = observed_lambda,
                      gamma = observed_shape)

effect_model <- list(delay_SHELF = SHELF::fitdist(c(observed_delay - 0.001, observed_delay, observed_delay + 0.001),
                                                  probs = c(0.25, 0.5, 0.75), lower = 0, upper = 10),
                     delay_dist = "gamma",
                     HR_SHELF = SHELF::fitdist(c(observed_HR - 0.001, observed_HR, observed_HR + 0.001), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1.5),
                     HR_dist = "gamma",
                     P_S = 1, P_DTE = 1)

calc_dte_assurance(n_c =  n_c,
                   n_t = n_t,
                   control_model = control_model,
                   effect_model = effect_model,
                   censoring_model = censoring_model,
                   recruitment_model = recruitment_model,
                   analysis_model = analysis_model,
                   n_sims = 10000)

## Realised Power - Design with observed delay ----------------------

target_HR <- 0.75
observed_delay <- 5.3


control_model <- list(dist = "Exponential",
                      parameter_model = "Fixed",
                      fixed_type = "Landmark",
                      t1 = 14,
                      surv_t1 = 0.5)

effect_model <- list(delay_SHELF = SHELF::fitdist(c(observed_delay - 0.001, observed_delay, observed_delay + 0.001),
                                                  probs = c(0.25, 0.5, 0.75), lower = 0, upper = 10),
                     delay_dist = "gamma",
                     HR_SHELF = SHELF::fitdist(c(observed_HR - 0.001, observed_HR, observed_HR + 0.001), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1.5),
                     HR_dist = "gamma",
                     P_S = 1, P_DTE = 1)

calc_dte_assurance(n_c =  n_c,
                   n_t = n_t,
                   control_model = control_model,
                   effect_model = effect_model,
                   censoring_model = censoring_model,
                   recruitment_model = recruitment_model,
                   analysis_model = analysis_model,
                   n_sims = 10000)


# Trial B ----------------------------
n_c <- 411
n_t <- 410

censoring_model <- list(method = "Events", events = 569)

recruitment_model <- list(method = "power",
                          period = 17.4,
                          power = 1)

analysis_model <- list(method = "LRT",
                       alpha = 0.025,
                       alternative_hypothesis = "one.sided")

## Nominal Power ----------------------


target_HR <- 0.76

control_model <- list(dist = "Exponential",
                      parameter_model = "Fixed",
                      fixed_type = "Landmark",
                      t1 = 14.8,
                      surv_t1 = 0.5)
effect_model <- list(delay_SHELF = SHELF::fitdist(c(3, 4, 5),
                                                  probs = c(0.25, 0.5, 0.75), lower = 0, upper = 10),
                     delay_dist = "gamma",
                     HR_SHELF = SHELF::fitdist(c(target_HR - 0.001, target_HR, target_HR + 0.001),
                                               probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1.5),
                     HR_dist = "gamma",
                     P_S = 1, P_DTE = 0)


calc_dte_assurance(n_c =  n_c,
                   n_t = n_t,
                   control_model = control_model,
                   effect_model = effect_model,
                   censoring_model = censoring_model,
                   recruitment_model = recruitment_model,
                   analysis_model = analysis_model,
                   n_sims = 10000)


## Realised Power - Best Fitting ----------------------

observed_delay <- 0.1
observed_lambda <- 0.0303
observed_HR <- 0.743

control_model <- list(dist = "Exponential",
                      parameter_model = "Fixed",
                      fixed_type = "Parameters",
                      lambda = observed_lambda)

effect_model <- list(delay_SHELF = SHELF::fitdist(c(observed_delay - 0.001, observed_delay, observed_delay + 0.001),
                                                  probs = c(0.25, 0.5, 0.75), lower = 0, upper = 10),
                     delay_dist = "gamma",
                     HR_SHELF = SHELF::fitdist(c(observed_HR - 0.001, observed_HR, observed_HR + 0.001), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1.5),
                     HR_dist = "gamma",
                     P_S = 1, P_DTE = 0)


calc_dte_assurance(n_c =  n_c,
                   n_t = n_t,
                   control_model = control_model,
                   effect_model = effect_model,
                   censoring_model = censoring_model,
                   recruitment_model = recruitment_model,
                   analysis_model = analysis_model,
                   n_sims = 10000)

# Trial C ----------------------------

n_c <- 410/2
n_t <- 410/2

censoring_model <- list(method = "Events", events = 312)

recruitment_model <- list(method = "power",
                          period = 13.2,
                          power = 1)

analysis_model <- list(method = "LRT",
                       alpha = 0.025,
                       alternative_hypothesis = "one.sided")

## Nominal Power ----------------------

target_HR <- 0.69

control_model <- list(dist = "Exponential",
                      parameter_model = "Fixed",
                      fixed_type = "Landmark",
                      t1 = 10,
                      surv_t1 = 0.5)

effect_model <- list(delay_SHELF = SHELF::fitdist(c(3, 4, 5),
                                                  probs = c(0.25, 0.5, 0.75), lower = 0, upper = 10),
                     delay_dist = "gamma",
                     HR_SHELF = SHELF::fitdist(c(target_HR - 0.001, target_HR, target_HR + 0.001),
                                               probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1.5),
                     HR_dist = "gamma",
                     P_S = 1, P_DTE = 0)


calc_dte_assurance(n_c =  n_c,
                   n_t = n_t,
                   control_model = control_model,
                   effect_model = effect_model,
                   censoring_model = censoring_model,
                   recruitment_model = recruitment_model,
                   analysis_model = analysis_model,
                   n_sims = 10000)

## Realised Power - Best Fitting ----------------------

observed_delay <- 0.1
observed_lambda <- 0.0419
observed_shape <- 0.784
observed_HR <- 0.475


control_model <- list(dist = "Weibull",
                      parameter_model = "Fixed",
                      fixed_type = "Parameters",
                      lambda = observed_lambda,
                      gamma = observed_shape)

effect_model <- list(delay_SHELF = SHELF::fitdist(c(observed_delay - 0.001, observed_delay, observed_delay + 0.001),
                                                  probs = c(0.25, 0.5, 0.75), lower = 0, upper = 10),
                     delay_dist = "gamma",
                     HR_SHELF = SHELF::fitdist(c(observed_HR - 0.001, observed_HR, observed_HR + 0.001), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1.5),
                     HR_dist = "gamma",
                     P_S = 1, P_DTE = 0)


calc_dte_assurance(n_c =  n_c,
                   n_t = n_t,
                   control_model = control_model,
                   effect_model = effect_model,
                   censoring_model = censoring_model,
                   recruitment_model = recruitment_model,
                   analysis_model = analysis_model,
                   n_sims = 10000)


# Trial D ----------------------------

n_c <- 915/3
n_t <- 915/3

censoring_model <- list(method = "Events", events = 460)

recruitment_model <- list(method = "power",
                          period = 8.1,
                          power = 1)

analysis_model <- list(method = "LRT",
                       alpha = 0.025,
                       alternative_hypothesis = "one.sided")

## Nominal Power ----------------------

target_HR <- 0.72

control_model <- list(dist = "Exponential",
                      parameter_model = "Fixed",
                      fixed_type = "Landmark",
                      t1 = 14,
                      surv_t1 = 0.5)
effect_model <- list(delay_SHELF = SHELF::fitdist(c(3, 4, 5),
                                                  probs = c(0.25, 0.5, 0.75), lower = 0, upper = 10),
                     delay_dist = "gamma",
                     HR_SHELF = SHELF::fitdist(c(target_HR - 0.001, target_HR, target_HR + 0.001),
                                               probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1.5),
                     HR_dist = "gamma",
                     P_S = 1, P_DTE = 0)


calc_dte_assurance(n_c =  n_c,
                   n_t = n_t,
                   control_model = control_model,
                   effect_model = effect_model,
                   censoring_model = censoring_model,
                   recruitment_model = recruitment_model,
                   analysis_model = analysis_model,
                   n_sims = 10000)


#Trial D1 ------------------------


## Realised Power - Best Fitting ----------------------

observed_delay <- 12.3
observed_lambda <- 0.0248
observed_shape <- 0.895
observed_HR <- 0.314


control_model <- list(dist = "Weibull",
                      parameter_model = "Fixed",
                      fixed_type = "Parameters",
                      lambda = observed_lambda,
                      gamma = observed_shape)

effect_model <- list(delay_SHELF = SHELF::fitdist(c(observed_delay - 0.001, observed_delay, observed_delay + 0.001),
                                                  probs = c(0.25, 0.5, 0.75), lower = 0, upper = 15),
                     delay_dist = "gamma",
                     HR_SHELF = SHELF::fitdist(c(observed_HR - 0.001, observed_HR, observed_HR + 0.001), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1.5),
                     HR_dist = "gamma",
                     P_S = 1, P_DTE = 1)


calc_dte_assurance(n_c =  n_c,
                   n_t = n_t,
                   control_model = control_model,
                   effect_model = effect_model,
                   censoring_model = censoring_model,
                   recruitment_model = recruitment_model,
                   analysis_model = analysis_model,
                   n_sims = 10000)

## Realised Power - Design with observed delay ----------------------

target_HR <- 0.72

control_model <- list(dist = "Exponential",
                      parameter_model = "Fixed",
                      fixed_type = "Landmark",
                      t1 = 14,
                      surv_t1 = 0.5)

observed_delay <- 12.3


effect_model <- list(delay_SHELF = SHELF::fitdist(c(observed_delay - 0.001, observed_delay, observed_delay + 0.001),
                                                  probs = c(0.25, 0.5, 0.75), lower = 0, upper = 15),
                     delay_dist = "gamma",
                     HR_SHELF = SHELF::fitdist(c(target_HR - 0.001, target_HR, target_HR + 0.001), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1.5),
                     HR_dist = "gamma",
                     P_S = 1, P_DTE = 1)


calc_dte_assurance(n_c =  n_c,
                   n_t = n_t,
                   control_model = control_model,
                   effect_model = effect_model,
                   censoring_model = censoring_model,
                   recruitment_model = recruitment_model,
                   analysis_model = analysis_model,
                   n_sims = 10000)


#Trial D2 ----------------------


## Realised Power - Best Fitting ----------------------

observed_delay <- 8.2
observed_lambda <- 0.0252
observed_shape <- 0.887
observed_HR <- 0.507


control_model <- list(dist = "Weibull",
                      parameter_model = "Fixed",
                      fixed_type = "Parameters",
                      lambda = observed_lambda,
                      gamma = observed_shape)

effect_model <- list(delay_SHELF = SHELF::fitdist(c(observed_delay - 0.001, observed_delay, observed_delay + 0.001),
                                                  probs = c(0.25, 0.5, 0.75), lower = 0, upper = 10),
                     delay_dist = "gamma",
                     HR_SHELF = SHELF::fitdist(c(observed_HR - 0.001, observed_HR, observed_HR + 0.001), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1.5),
                     HR_dist = "gamma",
                     P_S = 1, P_DTE = 1)


calc_dte_assurance(n_c =  n_c,
                   n_t = n_t,
                   control_model = control_model,
                   effect_model = effect_model,
                   censoring_model = censoring_model,
                   recruitment_model = recruitment_model,
                   analysis_model = analysis_model,
                   n_sims = 10000)

## Realised Power - Design with observed delay ----------------------

target_HR <- 0.72

control_model <- list(dist = "Exponential",
                      parameter_model = "Fixed",
                      fixed_type = "Landmark",
                      t1 = 14,
                      surv_t1 = 0.5)

observed_delay <- 8.2


effect_model <- list(delay_SHELF = SHELF::fitdist(c(observed_delay - 0.001, observed_delay, observed_delay + 0.001),
                                                  probs = c(0.25, 0.5, 0.75), lower = 0, upper = 15),
                     delay_dist = "gamma",
                     HR_SHELF = SHELF::fitdist(c(target_HR - 0.001, target_HR, target_HR + 0.001), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1.5),
                     HR_dist = "gamma",
                     P_S = 1, P_DTE = 1)


calc_dte_assurance(n_c =  n_c,
                   n_t = n_t,
                   control_model = control_model,
                   effect_model = effect_model,
                   censoring_model = censoring_model,
                   recruitment_model = recruitment_model,
                   analysis_model = analysis_model,
                   n_sims = 10000)


# Trial E ----------------------------

n_c <- 264/2
n_t <- 264/2

censoring_model <- list(method = "Events", events = 189)

recruitment_model <- list(method = "power",
                          period = 14.2,
                          power = 1)

analysis_model <- list(method = "LRT",
                       alpha = 0.025,
                       alternative_hypothesis = "one.sided")


## Nominal Power ----------------------


target_HR <- 0.61

control_model <- list(dist = "Exponential",
                      parameter_model = "Fixed",
                      fixed_type = "Landmark",
                      t1 = 7,
                      surv_t1 = 0.5)
effect_model <- list(delay_SHELF = SHELF::fitdist(c(3, 4, 5),
                                                  probs = c(0.25, 0.5, 0.75), lower = 0, upper = 10),
                     delay_dist = "gamma",
                     HR_SHELF = SHELF::fitdist(c(target_HR - 0.001, target_HR, target_HR + 0.001),
                                               probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1.5),
                     HR_dist = "gamma",
                     P_S = 1, P_DTE = 0)




calc_dte_assurance(n_c =  n_c,
                   n_t = n_t,
                   control_model = control_model,
                   effect_model = effect_model,
                   censoring_model = censoring_model,
                   recruitment_model = recruitment_model,
                   analysis_model = analysis_model,
                   n_sims = 10000)


## Realised Power - Best Fitting ----------------------

observed_delay <- 0.1
observed_lambda <- 0.0994
observed_shape <- 0.866 
observed_HR <- 0.572

control_model <- list(dist = "Weibull",
                      parameter_model = "Fixed",
                      fixed_type = "Parameters",
                      lambda = observed_lambda,
                      gamma = observed_shape)

effect_model <- list(delay_SHELF = SHELF::fitdist(c(observed_delay - 0.001, observed_delay, observed_delay + 0.001),
                                                  probs = c(0.25, 0.5, 0.75), lower = 0, upper = 10),
                     delay_dist = "gamma",
                     HR_SHELF = SHELF::fitdist(c(observed_HR - 0.001, observed_HR, observed_HR + 0.001), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1.5),
                     HR_dist = "gamma",
                     P_S = 1, P_DTE = 0)


calc_dte_assurance(n_c =  n_c,
                   n_t = n_t,
                   control_model = control_model,
                   effect_model = effect_model,
                   censoring_model = censoring_model,
                   recruitment_model = recruitment_model,
                   analysis_model = analysis_model,
                   n_sims = 10000)



# Trial G ----------------------------

n_c <- 120
n_t <- 240

censoring_model <- list(method = "Events", events = 281)

recruitment_model <- list(method = "power",
                          period = 14.2,
                          power = 1)

analysis_model <- list(method = "LRT",
                       alpha = 0.025,
                       alternative_hypothesis = "one.sided")

## Nominal Power ----------------------


target_HR <- 2/3

control_model <- list(dist = "Exponential",
                      parameter_model = "Fixed",
                      fixed_type = "Landmark",
                      t1 = 6,
                      surv_t1 = 0.5)
effect_model <- list(delay_SHELF = SHELF::fitdist(c(3, 4, 5),
                                                  probs = c(0.25, 0.5, 0.75), lower = 0, upper = 10),
                     delay_dist = "gamma",
                     HR_SHELF = SHELF::fitdist(c(target_HR - 0.001, target_HR, target_HR + 0.001),
                                               probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1.5),
                     HR_dist = "gamma",
                     P_S = 1, P_DTE = 0)



calc_dte_assurance(n_c =  n_c,
                   n_t = n_t,
                   control_model = control_model,
                   effect_model = effect_model,
                   censoring_model = censoring_model,
                   recruitment_model = recruitment_model,
                   analysis_model = analysis_model,
                   n_sims = 10000)


## Realised Power - Best Fitting ----------------------

observed_delay <- 3.3
observed_lambda <- 0.1182
observed_HR <- 0.539

control_model <- list(dist = "Exponential",
                      parameter_model = "Fixed",
                      fixed_type = "Parameters",
                      lambda = observed_lambda)

effect_model <- list(delay_SHELF = SHELF::fitdist(c(observed_delay - 0.001, observed_delay, observed_delay + 0.001),
                                                  probs = c(0.25, 0.5, 0.75), lower = 0, upper = 10),
                     delay_dist = "gamma",
                     HR_SHELF = SHELF::fitdist(c(observed_HR - 0.001, observed_HR, observed_HR + 0.001), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1.5),
                     HR_dist = "gamma",
                     P_S = 1, P_DTE = 1)


calc_dte_assurance(n_c =  n_c,
                   n_t = n_t,
                   control_model = control_model,
                   effect_model = effect_model,
                   censoring_model = censoring_model,
                   recruitment_model = recruitment_model,
                   analysis_model = analysis_model,
                   n_sims = 10000)


## Realised Power - Design with observed delay ----------------------

observed_delay <- 3.3
observed_lambda <- log(2)/6
observed_HR <- 2/3

control_model <- list(dist = "Exponential",
                      parameter_model = "Fixed",
                      fixed_type = "Parameters",
                      lambda = observed_lambda)

effect_model <- list(delay_SHELF = SHELF::fitdist(c(observed_delay - 0.001, observed_delay, observed_delay + 0.001),
                                                  probs = c(0.25, 0.5, 0.75), lower = 0, upper = 10),
                     delay_dist = "gamma",
                     HR_SHELF = SHELF::fitdist(c(observed_HR - 0.001, observed_HR, observed_HR + 0.001), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1.5),
                     HR_dist = "gamma",
                     P_S = 1, P_DTE = 1)


calc_dte_assurance(n_c =  n_c,
                   n_t = n_t,
                   control_model = control_model,
                   effect_model = effect_model,
                   censoring_model = censoring_model,
                   recruitment_model = recruitment_model,
                   analysis_model = analysis_model,
                   n_sims = 10000)


# Trial H ----------------------------

n_c <- 574/2
n_t <- 574/2

censoring_model <- list(method = "Events", events = 403)

recruitment_model <- list(method = "power",
                          period = 13.1,
                          power = 1)

analysis_model <- list(method = "LRT",
                       alpha = 0.025,
                       alternative_hypothesis = "one.sided")

## Nominal Power ----------------------

target_HR <- 0.72

control_model <- list(dist = "Exponential",
                      parameter_model = "Fixed",
                      fixed_type = "Landmark",
                      t1 = 8,
                      surv_t1 = 0.5)
effect_model <- list(delay_SHELF = SHELF::fitdist(c(3, 4, 5),
                                                  probs = c(0.25, 0.5, 0.75), lower = 0, upper = 10),
                     delay_dist = "gamma",
                     HR_SHELF = SHELF::fitdist(c(target_HR - 0.001, target_HR, target_HR + 0.001),
                                               probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1.5),
                     HR_dist = "gamma",
                     P_S = 1, P_DTE = 0)


calc_dte_assurance(n_c =  n_c,
                   n_t = n_t,
                   control_model = control_model,
                   effect_model = effect_model,
                   censoring_model = censoring_model,
                   recruitment_model = recruitment_model,
                   analysis_model = analysis_model,
                   n_sims = 10000)


## Realised Power - Best Fitting ----------------------


observed_delay <- 6.6
observed_lambda <- 0.0683
observed_HR <- 0.483

control_model <- list(dist = "Exponential",
                      parameter_model = "Fixed",
                      fixed_type = "Parameters",
                      lambda = observed_lambda)

effect_model <- list(delay_SHELF = SHELF::fitdist(c(observed_delay - 0.001, observed_delay, observed_delay + 0.001),
                                                  probs = c(0.25, 0.5, 0.75), lower = 0, upper = 10),
                     delay_dist = "gamma",
                     HR_SHELF = SHELF::fitdist(c(observed_HR - 0.001, observed_HR, observed_HR + 0.001), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1.5),
                     HR_dist = "gamma",
                     P_S = 1, P_DTE = 1)


calc_dte_assurance(n_c =  n_c,
                   n_t = n_t,
                   control_model = control_model,
                   effect_model = effect_model,
                   censoring_model = censoring_model,
                   recruitment_model = recruitment_model,
                   analysis_model = analysis_model,
                   n_sims = 10000)

## Realised Power - Design with observed delay ----------------------

observed_delay <- 3.3
observed_lambda <- log(2)/8
observed_HR <- 0.72

control_model <- list(dist = "Exponential",
                      parameter_model = "Fixed",
                      fixed_type = "Parameters",
                      lambda = observed_lambda)

effect_model <- list(delay_SHELF = SHELF::fitdist(c(observed_delay - 0.001, observed_delay, observed_delay + 0.001),
                                                  probs = c(0.25, 0.5, 0.75), lower = 0, upper = 10),
                     delay_dist = "gamma",
                     HR_SHELF = SHELF::fitdist(c(observed_HR - 0.001, observed_HR, observed_HR + 0.001), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1.5),
                     HR_dist = "gamma",
                     P_S = 1, P_DTE = 1)


calc_dte_assurance(n_c =  n_c,
                   n_t = n_t,
                   control_model = control_model,
                   effect_model = effect_model,
                   censoring_model = censoring_model,
                   recruitment_model = recruitment_model,
                   analysis_model = analysis_model,
                   n_sims = 10000)

