
# RUN OLS/WLS SCRIPT ------------------------------------------------------

source("rs_regional_daily_ols_wls.r")


# ASSESS RESULTS ----------------------------------------------------------

# Choose WLS model with smallest AVP1
df_wls_result <- as_tibble(wls_result)
df_wls_model <- df_wls_result[which.min(df_wls_result$avp1),] # extract only row with lowest APV1

# The model with lowest average prediction variance represents the one with smallest model and sampling variance combined.
# To check how we should approach building a model for kappa, we can also check if the parameters are significantlly 
# different from zero by running the t-test

# Get the variables and their standard error
b_cols <- grep("^b[1-9]", names(df_wls_model), value = TRUE)     # get parameter columns above b0
se_cols <- grep("^se[1-9]", names(df_wls_model), value = TRUE)   # get standard error columns above se0
df_wls_par <- tibble(var = sub("b", "", b_cols),                 # table with variables and their respective standard errors
                     b = as.numeric(df_wls_model[1, b_cols]),
                     se = as.numeric(df_wls_model[1, se_cols]))
df_wls_par <- df_wls_par[df_wls_par[[2]] != 0, ]                 # select row where b != 0
# não gostei dessa forma de calcular

# Determine t-statistic and p-values
sig_level <- 0.05
df_wls_par$t_stat <- df_wls_par$b/df_wls_par$se
df_wls_par$p_value <- 2*pt(abs(df_wls_par$t_stat), df = n - p - 1, lower.tail = FALSE)
df_wls_par$signf <- df_wls_par$p_value < sig_level

# For all models with signf = TRUE, extract explanatory variables