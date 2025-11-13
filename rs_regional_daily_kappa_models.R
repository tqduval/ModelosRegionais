
# RUN OLS/WLS SCRIPT ------------------------------------------------------

source("rs_regional_daily_ols_wls.r")
# saveRDS(df_exp_variables, "Dados Gerados/df_exp_variables.rds")
# saveRDS(wls_result, "Dados Gerados/wls_result.rds")
# saveRDS(ols_result, "Dados Gerados/ols_result.rds"_


# ASSESS RESULTS ----------------------------------------------------------

# Choose WLS model with smallest AVP1
df_wls_result <- as_tibble(wls_result)                        # turn into tbl_df
df_wls_model <- df_wls_result[which.min(df_wls_result$avp1),] # extract row with smallest avp1

# Get dynamic vectors for b and se
b <- t(df_wls_model)[seq(1, p)]       # from one (b0) to the last parameter (bp)
se <- t(df_wls_model)[seq(p + 1,2*p)] # from the next column to bp with same length as p

# Test if bp (excluding b0) are statistically different from zero
sigf_level <- 0.01    # set significance level of 1%
t <- rep(0, p)        # pre-define t statistic vector
pvalue <- rep(0, p)   # pre-define pvalue vector
signf <- rep(TRUE, p) # pre-define significance logical vector

for(k in 1:p){
  
  t[k] <- b[k]/se[k] # t-statistic
  pvalue[k] <- 2*pt(abs(t[k]),df = n - p - 1, lower.tail = FALSE)
  signf[k] <- !is.na(pvalue[k]) & (pvalue[k] < sigf_level)
  
  if(signf[k] == FALSE) b[k] <- 0
  
}
# passar p/ ols_wls e escolher o segundo melhor caso esse não seja significativo

# BUILD MODELS FOR KAPPA --------------------------------------------------

# Append gauge vector to df_exp_variables
df_regional_kappa <- cbind(gauges = gauges, as_tibble(df_exp_variables), mu_kappa_r = rep(0, n), sd_kappa_r = rep(0, n))
kappa_r_names <- c("mu_kappa_r", "sd_kappa_r")
wls_var_delta <- df_wls_model$var_delta

# Compute wls_cov_par again with the choosen model
wls_inv_lambda <- fun_gls_var_delta(std.delta = wls_var_delta^0.5, n, scov[[2]], kappa, df_exp_variables, p)$inv.lambda
wls_cov_par <- solve(t(df_exp_variables) %*% wls_inv_lambda %*% df_exp_variables)

for(i in seq_along(gauges)){
  
  gg <- gauges[i]
  xi <- as.matrix(df_regional_kappa[i, model_names])
  
  df_regional_kappa[i, kappa_r_names[1]] <- xi %*% b
  df_regional_kappa[i, kappa_r_names[2]] <- (wls_var_delta + xi %*% wls_cov_par %*% t(xi))^0.5
  
}

# Save
saveRDS(df_regional_kappa, "Dados Gerados/df_regional_kappa.rds")

hist(df_regional_kappa$mu_kappa_r, breaks = 15)