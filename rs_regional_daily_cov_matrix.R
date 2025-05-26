
# PACKAGES ----------------------------------------------------------------

pacman::p_load(pacman, tidyverse, lubridate, readxl, janitor, lmom, gsl)


# READ DATA ---------------------------------------------------------------

# Daily files RS
path_rs_daily <- file.choose()
df_rs_daily <- read_excel(path_rs_daily)
col_names <- names(df_rs_daily)[c(1,4,6,7,12,13,17)]
df_rs_daily <- as_tibble(clean_names(df_rs_daily[, col_names]))


# GET L-MOMENTS AND KAPPA -------------------------------------------------

gauges <- unique(df_rs_daily$gauge_code)
m <- length(gauges)
kappas <- rep(0, m)
tau3 <- rep(0, m)
n <- rep(0, m)

# Get L-moments for all stations and fit kappa
for(gg in seq_along(gauges)){
  
  # Get am series for current gauge
  cod <- gauges[[gg]]
  curr_gg <- df_rs_daily[df_rs_daily[col_names[1]] == cod, col_names[5], drop = TRUE] # vector
  n[gg] <- length(curr_gg)
  
  # Calculate sample lmoments and ratios
  lmom <- lmom::samlmu(x = curr_gg, sort.data = TRUE) # sample lmoments
  tau3[gg] <- lmom[[3]]                               # l-skewness (3rd lmoment ratio)
  kappas[gg] <- lmom::pelgev(lmom)[[3]]               # fit GEV and extract kappa
  
  
}

# Clear memory
rm(gg, lmom, cod, curr_gg)
gc()

# Average tau3 and kappa
tau3_a <- mean(tau3)                          # average tau3
c = 2/(3 + tau3_a) - log(2)/log(3)
k <- 7.859*c + 2.9554*c^2                     # average kappa
n_a <- mean(n)                                # average number of years

# Tests Table 1 – Chowdhury et al. (1991)
# k <- -0.2                                   # average kappa
# tau3_a <- (2*(1 - 3^(-k)))/(1 - 2^(-k)) - 3 # average tau3 for testing
# c = 2/(3 + tau3_a) - log(2)/log(3)
# n_a <- 500                                  # average number of years for testing


# L-SKEWNESS AND KAPPA VARIANCES ------------------------------------------

# Compute var(tau3) and var(kappa) – Chowdhury et al. (1991) and Fill [PhD Thesis Cornell] (1994)
r <- gamma(1 + 2*k)/(gamma(1 + k)^2) # R(k) function
aux <- (1 - 2^(-k))^2                # auxiliary denominator that shows up a lot

# Using the gsl::hyperg_2F1() function as MatLabs 'mhygfx'
hy05 <- gsl::hyperg_2F1(k, 2*k, 1 + k, -0.5)
hy13 <- gsl::hyperg_2F1(k, 2*k, 1 + k, -1/3)
hy23 <- gsl::hyperg_2F1(k, 2*k, 1 + k, -2/3)

# Elements fij from the symetric variance matrix for l-moment ratios T
f11 <- (r - 1)/aux
f22 <- (r*hy05 - 1)/(2^(2*k)*aux)
f33 <- (r*hy23 - 1)/(3^(2*k)*aux)
f12 <- 0.5*(r - 2^(1 + k) + 2^(2*k))/(2^(2*k)*aux)
f13 <- 0.5*((r - 2*3^k)/(3^(2*k)*aux) - (r*hy05 - 2^(1 + k))/(2^(2*k)*aux))
f23 <- 0.5*((r*hy13 - 2*(3/2)^k)/(3^(2*k)*aux) + 1/(2^(2*k)*aux))

# Variance of tau3 and kappa
var_tau3_a <- 1/n_a*((f11 - 4*f12 + 4*f22)*tau3_a^2 + 2*(f11 - 8*f12 + 12*f22 + 6*f13 - 12*f23)*tau3_a + f11 - 12*f12 + 36*f22 + 12*f13 - 72*f23 + 36*f33)
dk_dtau3_a <- -(2.872/(tau3_a + 3))^2*(2.872/(tau3_a + 3) + 1) # partial kappa in relation to tau3
var_tau3 <- (n_a/n)*var_tau3_a                                 # vector of tau3 variances
var_k_a <- dk_dtau3_a^2*var_tau3_a                             # average variance of kappa
var_k <- (n_a/n)*var_k_a                                       # variance of kappa for each station (vector)
na_var_tau3 <- n_a*var_tau3_a


# BUILD WLS SCOV MATRIX ---------------------------------------------------

# Create empty matrix and fill diagonals with kappa variance
wls_scovmatrix <- matrix(0, ncol = m, nrow = m, dimnames = list(NULL, gauges))
for(gg in seq_along(gauges)){
  wls_scovmatrix[gg,gg] <- var_k[gg]
}


# TEST FUNCTION -----------------------------------------------------------

# Get fun_wls_scovmatrix()
source("Funções/fun_wls_scovmatrix.r")

# Test function
results <- fun_wls_scovmatrix(df = df_rs_daily, col = c("gauge_code", "peak_24h"))
wls_scovmatrix <- results[["scovmatrix"]]