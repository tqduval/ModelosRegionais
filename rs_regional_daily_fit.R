
# PACKAGES ----------------------------------------------------------------

if(!require("pacman")) install.packages("pacman"); pacman::p_load(pacman, tidyverse)


# READ DATA ---------------------------------------------------------------

# Clear environment
rm(list = ls())

# Data
df_regional_kappa <- readRDS("Dados Gerados/df_regional_kappa.rds") # regional kappa models
df_wls_result <- as_tibble(readRDS("Dados Gerados/wls_result.rds")) # wls results
df_rs_daily <- readRDS("Dados Gerados/df_rs_daily.rds")             # daily rainfall data

# Functions
# R File 0_Funcoes.r has all auxiliary functions for the main GEV fitting function
source("C:/Users/tomas/OneDrive/1 - Acadêmico/Mestrado/Tese/3-R/Precipitacão Máxima Anual/PrecipitacaoMaximaAnual/0_Funcoes.r")

# Function for fitting the GEV with LMOM, MLE and GMLE with geophysical prior and regional prior
# This function is based on the ones in the "2_Ajuste.r" file named fun.process() and fun.ajuste.pmax(),
# but it doesnt calculate any performance metrics like LRT, AIC or BIC
source("C:/Users/Tomas/OneDrive/1 - Acadêmico/Mestrado/Tese/3-R/Precipitacão Máxima Anual/PrecipitacaoMaximaAnual/Funções/fun_gev_fit_regional.r")

# Function for computing CI for parameters
source("C:/Users/Tomas/OneDrive/1 - Acadêmico/Mestrado/Tese/3-R/Precipitacão Máxima Anual/PrecipitacaoMaximaAnual/Funções/fun_gev_ci_parameters.r")

# Function for computing GEV quantiles
# This function is based on the one in "3_Quantis.r", but with regional quantiles added
source("C:/Users/Tomas/OneDrive/1 - Acadêmico/Mestrado/Tese/3-R/Precipitacão Máxima Anual/PrecipitacaoMaximaAnual/Funções/fun_gev_quantiles_regional.r")



# FIT GEV -----------------------------------------------------------------

# Fit GEV via L-MOM, MLE and GMLE with geophysical prior and regional prior
df.fit <- fun.gev.fit.regional(df = df_rs_daily,
                               regional.models = df_regional_kappa,
                               data.names = c("gauge_code", "peak_24h"),
                               model.names = c("gauges", "mu_kappa_r", "sd_kappa_r"))

# There still needs to be a step here to compute kappa variances using the hessian
df.kappa.ci <- fun.gev.ci.par(df = df.fit,
                              gauge.names = "gauge.code",
                              models = c("l.gev", "gl.gev", "gl.gev.r"),
                              pars = "kappa")


# QUANTILES ---------------------------------------------------------------

# Tempos de retorno
tempo.retorno <- c(10.0, 20.0, 50.0, 100.0, 200.0, 500.0)
df.quantiles <- fun.gev.quantile.regional(df = df.fit,
                                          gauge.names = "gauge.code",
                                          tempo.retorno = tempo.retorno,
                                          signf.level = 0.05,
                                          model = c("l.gev", "gl.gev", "gl.gev.r"))


# SAVE RESULTS ------------------------------------------------------------

saveRDS(df.fit, "Dados Gerados/df.fit.rds")
saveRDS(df.quantiles, "Dados Gerados/df.quantiles.rds")
saveRDS(df.kappa.ci, "Dados Gerados/df.kappa.ci.rds")