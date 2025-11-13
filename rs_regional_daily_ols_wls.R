
# PACKAGES ----------------------------------------------------------------

pacman::p_load(pacman, tidyverse, psych, readxl, janitor)


# FUNCTIONS ---------------------------------------------------------------

source("Funções/fun_wls_scovmatrix.r") # build sampling covariance matrix for WLS
source("0_Funcoes.r")                  # fun.ex.var and fun.gls.var.delta

# READ AND PREPARE DATA ---------------------------------------------------

# Daily files RS
# path_rs_daily <- file.choose()
path_rs_daily <- "C:\\Users\\tomas\\OneDrive\\1 - Acadêmico\\Mestrado\\Tese\\3-R\\Precipitacão Máxima Anual\\PrecipitacaoMaximaAnual\\Fonte de Dados\\Dados_RS\\Dados máximos diários consolidados\\Máximos anuais para IDF.xlsx"
df_rs_daily <- read_excel(path_rs_daily)
col_names <- names(df_rs_daily)[c(1,4,6,7,12,13,17)]
df_rs_daily <- as_tibble(clean_names(df_rs_daily[, col_names]))

# Build predictor variables matrix
df_exp_variables <- df_rs_daily[!duplicated(df_rs_daily[1]),]                             # remove duplicates
gauges <- df_exp_variables[[1]]                                                           # get station names vector
df_exp_variables <-  fun_exp_var(df = df_exp_variables,                                   # apply transformations
                                 which_var = c("gauge_code", "elevation", "long", "lat"),
                                 which_log = c("elevation"),
                                 which_scale = c("elevation", "long", "lat"))

# Build WLS sampling covariance matrix and get lmom kappa estimates
scov <- fun_wls_scovmatrix(df = df_rs_daily, col = c("gauge_code", "peak_24h"))
df_lmom_fit <- scov[[1]]    # extract lmom estimates for kappa other things
wls_scovmatrix <- scov[[2]] # extract wls sampling covariance matrix
kappa <- df_lmom_fit$kappa  # vector with kappa estimates

# Build models
model_names <- c("ones", "elevation", "long", "lat")                      # variable names
p <- length(model_names) - 1                                              # number of variables
models <- cbind(rep(1, 2^p), do.call(expand.grid, rep(list(c(0, 1)), p))) # build models matrix
colnames(models) <- model_names                                           # add names
n_predictors <- ncol(models)
n_models <- nrow(models)

# Set column names for results
res_col_names <- c(paste0("b", 0:(n_predictors - 1)),   # model parameters
                   paste0("se", 0:(n_predictors - 1)),  # model parameters standard error
                   "var_delta",                         # model error variance
                   "asve", "avp1", "avp2")              # evaluation metrics

# Create empty result matrices
ols_result <- matrix(0, nrow = n_models, ncol = length(res_col_names), dimnames = list(NULL, res_col_names)) # OLS
wls_result <- matrix(0, nrow = n_models, ncol = length(res_col_names), dimnames = list(NULL, res_col_names)) # WLS


# BUILD MODELS ------------------------------------------------------------

for(kmodel in 1:n_models){
  
  # Numeric vector indicating which model to use
  model_index <- which(models[kmodel,] > 0)
  
  # Extract predictor variables matrix X for the current model
  X <- df_exp_variables[, model_index, drop = FALSE]
  n <- nrow(X)
  p <- ncol(X)

  # # Pre-define dimensions for output parameters (b)
  ols_res_par <- matrix(0, nrow = ncol(models), ncol = 1)    # parâmetros do modelo OLS (beta.ols)
  ols_res_se_par <- matrix(0, nrow = ncol(models), ncol = 1) # standard error nos parâmetros do modelo OLS (sigma.delta.ols)
  wls_res_par <- matrix(0, nrow = ncol(models), ncol = 1)    # parâmetros do modelo WLS (beta.wls)
  wls_res_se_par <- matrix(0, nrow = ncol(models), ncol = 1) # standard error nos parâmetros do modelo WLS (sigma.delta.wls)
  
  Xt <- t(X)
  

  # OLS --------------------------------------------------------------------
  
  # Compute parameter, variance, covariance and standard error
  ols_par <- solve(Xt %*% X) %*% Xt %*% kappa                             # matriz com parametros do modelo
  ols_residual <- kappa - X %*% ols_par                                   # aux residuals matrix estimates - model
  ols_var_delta <- as.numeric((t(ols_residual) %*% ols_residual)/(n - p)) # model error variance
  ols_cov_par <- ols_var_delta * solve(Xt %*% X)                          # parameter covariance matrix
  ols_se_par <- sqrt(diag(ols_cov_par))                                   # standard error
  
  # Performance metrics
  ols_diag_var <- ols_var_delta*diag(n)                      # model error variance * identity matrix
  ols_inv_lambda <- solve(ols_diag_var)                      # model error variance matrix Λ-1
  ols_sve <- X %*% solve(Xt %*% ols_inv_lambda %*% X) %*% Xt # sampling variance
  ols_asve <- tr(ols_sve)/n                                  # average sampling variance
  ols_avp1 <- ols_var_delta + ols_asve                       # average variance of prediction
  term <- 2 * ols_var_delta * ols_sve %*% ols_inv_lambda
  ols_avp2 <- ols_avp1 - tr(term)/n
  
  # R^2
  

  # WLS --------------------------------------------------------------------
  
  # Extract sampling covariance matrix for WLS
  # wls_scov <- scov_matrix*diag(n)
  wls_scov <- wls_scovmatrix
  
  # Check if var_delta = 0 is enough, if not, optimize for check_delta > 0
  sd_delta <- 0
  wls_est <- fun_gls_var_delta(sd_delta, n, wls_scov, kappa, X, p)
  check_delta <- wls_est$check.delta
  sd_delta <- wls_est$std.delta
  wls_par <- wls_est$b
  wls_lambda <- wls_est$lambda
  wls_inv_lambda <- wls_est$inv.lambda
  
  if(check_delta < 0){
    
    sd_delta <- 0
    
    } else{
    
      # Define check_delta as a function of just sd_delta (g)
      f_root <- function(g){fun_gls_var_delta(g, n, wls_scov, kappa, X, p)$check.delta}
      wls_est_root <- uniroot(f_root, c(0, 20)) # find root of objective function f_root
      sd_delta <- wls_est_root$root             # extract root -> optimized sd_delta
      
      # Run parameter function again for new sd_delta
      wls_est <- fun_gls_var_delta(sd_delta, n, wls_scov, kappa, X, p)
      check_delta <- wls_est$check.delta
      sd_delta <- wls_est$std.delta
      wls_par <- wls_est$b
      wls_lambda <- wls_est$lambda
      wls_inv_lambda <- wls_est$inv.lambda
      
    }
  
  # Covariance and standard error for WLS parameters
  wls_cov_par <- solve(Xt %*% wls_inv_lambda %*% X)
  wls_se_par <- sqrt(diag(wls_cov_par))
  wls_var_delta <- sd_delta^2
  
  # Performance metrics
  wls_sve <- X %*% wls_cov_par %*% Xt
  wls_asve <- tr(wls_sve)/n
  wls_avp1 <- wls_var_delta + wls_asve
  term <- 2 * wls_var_delta * wls_sve %*% wls_inv_lambda
  wls_avp2 <- wls_avp1 - tr(term)/n


  # RESULTS ----------------------------------------------------------------
  
  # OLS
  ols_res_par[model_index] <- ols_par
  ols_res_se_par[model_index] <- ols_se_par
  ols_result[kmodel,] <- c(t(ols_res_par), t(ols_res_se_par), ols_var_delta, ols_asve, ols_avp1, ols_avp2)
  
  # WLS
  wls_res_par[model_index] <- wls_par
  wls_res_se_par[model_index] <- wls_se_par
  wls_result[kmodel,] <- c(t(wls_res_par), t(wls_res_se_par), wls_var_delta, wls_asve, wls_avp1, wls_avp2)
  
}

rm(ols_asve, ols_avp1, ols_avp2, ols_cov_par, ols_diag_var, ols_inv_lambda, ols_par, ols_res_par, ols_res_se_par, ols_residual, ols_sve, ols_se_par, ols_var_delta)
rm(wls_asve, wls_avp1, wls_avp2, wls_cov_par, wls_inv_lambda, wls_par, wls_res_par, wls_res_se_par, wls_sve, wls_scov, wls_lambda, wls_se_par, wls_var_delta, wls_scovmatrix, wls_est, wls_est_root)
rm(X, Xt, check_delta, col_names, kmodel, sd_delta, model_index, term)
gc()
