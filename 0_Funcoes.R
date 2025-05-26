# Função p/ transformar um data.frame e variáveis preditoras p/ cada estação
# Essa função transforma a primeira coluna (coluna geralmente com o nome
# das estações) em uma coluna com 1s e permite escolher um conjunto de colunas
# seguintes para aplicar log e outro conjunto para centralizar
fun_exp_var <- function(df, which_var = NULL, which_log, which_scale){
  
  if(!is.character(c(which_var, which_log, which_scale))){
    stop("\nArgumentos 'which_' devem conter vetor com nome das colunas")
  }
  
  df[,1] <- 1                                                                 # primeira coluna com 1s
  df <- df[, which_var]                                                       # caso haja outras informações na função
  colnames(df)[1] <- "ones"                                                   # renomear primeira coluna
  df[, which_log] <- log(df[, which_log])                                     # aplicar log nas colunas especificadas
  df[, which_scale] <- scale(df[, which_scale], center = TRUE, scale = FALSE) # centralizar colunas especificadas
  
  return(as.matrix(df))
  
}

# Função p/ estimar beta e a variância do erro do modelo – Stedinger e Tasker (1985)
fun_gls_var_delta <- function(std.delta,  # initial model error variance
                              n,          # number of stations
                              scovmatrix, # sample covariance matrix of theta
                              theta,      # hydrologic variable of interest
                              X,          # explanatory variables matrix 
                              p           # number of GLS parameters (or explanatory variables)
                              ){
  # Calcular matriz diagonal somente com variância dos erros do modelo (sigma_delta^2*In)
  diag.var.delta <- matrix(0, ncol = n, nrow = n)
  for(i in 1:n) diag.var.delta[i,i] <- std.delta^2
  # diag.var.delta <- std.delta^2*diag()
  
  # Matriz de covariância do erro total
  lambda <- diag.var.delta + scovmatrix
  
  # Calcular estimador de beta e variancia do erro do modelo sigma_delta^2
  inv.lambda <- solve(lambda)                                             # inversa da matriz de covariancia dos erros
  b <- solve(t(X) %*% inv.lambda %*% X) %*% t(X) %*% inv.lambda %*% theta # estimador de beta – Stedinger e Tasker (1985)
  fo <- t(theta - X %*% b) %*% inv.lambda %*% (theta - X %*% b) - (n - p) # função objetivo fo = 0
  fo <- as.numeric(fo)
  
  return(list(check.delta = fo,
              std.delta = std.delta,
              b = b,
              lambda = lambda,
              inv.lambda = inv.lambda))
  
}