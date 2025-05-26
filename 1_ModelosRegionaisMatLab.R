

# PACOTES ----------------------------------------------------------------

pacman::p_load(pacman,
               tidyverse,
               beepr,
               R.matlab,  # importar arquivos do MatLab
               psych)     # função tr() para somar diagonal principal de uma matriz


# CARREGAR TABELAS MATLAB -------------------------------------------------

# Diretório c/ arquivos MatLab
dir.matlab <- "C:/Users/tomas/OneDrive/1 - Acadêmico/Mestrado/Tese/0-BASE/MatLab_Modelos Regionais"

# Listar arquivos que terminam em .mat
files.matlab <- list.files(path = dir.matlab, pattern = "\\.mat$")

# Importar arquivos
dfs.matlab <- list()
for(file in files.matlab){
  dados <- R.matlab::readMat(paste0(dir.matlab, "/", file))
  nome.matlab <- gsub(".mat", "", file) # extrair nome dos arquivos
  dfs.matlab[[nome.matlab]] <- dados
  rm(dados)
}

# Carregar individualmente as tabelas
models <- dfs.matlab$models9var$models # matriz 512x10 com todas as combinações possíveis de 9 variáveis preenchida com 1s e 0s
# A matriz models é criada usando números de 1 a 512 no formato binário
scovmatrix <- dfs.matlab$ScCov0$ScCov0 # sampling covariance matrix 89x89 (são 89 estações)
skew <- dfs.matlab$thetaSc$theta       # variável hidrológica de interesse em cada estação (assimetria)

# Carregar matriz c/ 9 variáveis preditoras de 89 estações
df.variaveis <- read.table(file = paste0(dir.matlab, "/Scwholenew.txt"),
                           header = TRUE,
                           fileEncoding = "UTF-8")

# Matriz de variáveis preditoras (explanatórias)
df.expvar <- fun.exp.var(df = df.variaveis,
                         varexp = seq(1, ncol(df.variaveis)))


# CALCULAR MODELOS DE REGRESSÃO -------------------------------------------

# Pré-definir dimensões da tabela de resultados
n.predictors <- ncol(models) # nº de variáveis de preditoras
n.models <- nrow(models)     # nº de modelos

# Nome das colunas do dataframe final
# col.names <- c(paste0("par.", 1:n.predictors),
#                paste0("se.", 1:n.predictors),
#                "lambda", "asve", "avp1", "avp2", "pseudo.r2", "sst", "sse", "ssr")

col.names <- c(paste0("par.", 1:n.predictors),
               paste0("se.", 1:n.predictors),
               "lambda", "asve", "avp1", "avp2")

# Inicializar as matrizes de resultados
ols.result <- matrix(0,
                     nrow = nrow(models),
                     # ncol = 2*n.predictors + 8,
                     ncol = 2*n.predictors + 4, # tirei os relacionados ao R² pq dava dando problema
                     dimnames = list(NULL, col.names))

wls.result <- matrix(0,
                     nrow = nrow(models),
                     ncol = 2*n.predictors + 4,
                     dimnames = list(NULL, col.names))

for(kmodel in 1:nrow(models)){ # kmodel" varia de 1 a 512
  
  cat("Modelo", kmodel, "de", nrow(models),"avaliado...\n")
  
  # Vetor numérico que indica o índice das linhas (quais colunas) que serão usadas (iguais a 1)
  model.index <- which(models[kmodel,] > 0)
  
  # Matriz 89x10 com os valores das variáveis preditoras
  X <- matrix(data = 0,
              nrow = length(skew),
              ncol = length(model.index))
  
  X <- df.expvar[, model.index, drop = FALSE] %>% as.matrix # no fim essa vai ser a matriz X
  
  n <- nrow(X)
  p <- ncol(X)
  
  # Pré-definir as dimensões de dos parâmetros de saída (beta)
  # Essas matrizes terão o número de linhas igual ao número de colunas de 'models' e 1 coluna
  res.param.ols <- matrix(0, nrow = ncol(models), ncol = 1)    # parâmetros do modelo OLS (beta.ols)
  res.se.param.ols <- matrix(0, nrow = ncol(models), ncol = 1) # standard error nos parâmetros do modelo OLS (sigma.delta.ols)
  res.param.wls <- matrix(0, nrow = ncol(models), ncol = 1)    # parâmetros do modelo WLS (beta.wls)
  res.se.param.wls <- matrix(0, nrow = ncol(models), ncol = 1) # standard error nos parâmetros do modelo WLS (sigma.delta.ols)
  res.param.gls <- matrix(0, nrow = ncol(models), ncol = 1)    # parâmetros do modelo GLS (beta.gls)
  res.se.param.gls <- matrix(0, nrow = ncol(models), ncol = 1) # standard error nos parâmetros do modelo GLS (sigma.delta.gls)
  
  # Matriz auxiliar → [(X^T*X)^(-1)]*X^T
  beta.aux <- solve(t(X) %*% X) %*% t(X)
  
  
  ## Parâmetros do OLS -------------------------------------------------------
  
  # Calcular matrizes de parâmetros, variância, covariância e erro padrão
  ols.param <- beta.aux %*% skew                                          # produto entre beta.aux e a estimativa inicial de sk
  ols.residuals <- skew - X %*% ols.param                                 # resíduos
  ols.gamma <- as.numeric((t(ols.residuals) %*% (ols.residuals))/(n - p)) # variância nos resíduos
  ols.error.var <- ols.gamma                                              # variância no erro do modelo (sigma.delta^2)
  ols.cov.param <- as.numeric(ols.gamma) * solve(t(X) %*% X)              # matriz de covariância dos parâmetros
  ols.se.param <- sqrt(diag(ols.cov.param))                               # erros padrão dos parâmetros
  
  # Calcular métricas de avaliação
  ols.diag.var <- ols.gamma*diag(n) # matriz nxn c/ elementos da diagonal principal iguais a variância dos resíduos
  ols.inv.lambda <- solve(ols.diag.var)
  ols.sve <- X %*% solve(t(X) %*% ols.inv.lambda %*% X) %*% t(X)
  ols.asve <- tr(ols.sve)/n         # tr() calcula a soma dos elementos na diagonal principal
  ols.avp1 <- ols.gamma + ols.asve  # average variance of prediction
  ols.term <- 2 * ols.gamma * ols.sve %*% ols.inv.lambda
  ols.avp2 <- ols.avp1 - tr(ols.term)/n
  
  # # Calcular pseudo-R²
  # ols.sse <- t(ols.residuals) %*% ols.residuals
  # 
  # if(length(model.index) == 1){
  #   ols.mean.y <- mean(X)*ols.param %>% as.numeric
  #   ols.ssr <- t(skew - ols.mean.y) %*% (skew - ols.mean.y)
  #   ols.denominador.r2 <- ols.sse + ols.ssr
  # }
  # else{
  #   ols.mean.y <- mean(skew) * ols.param %>% as.numeric
  #   ols.predicted <- X %*% ols.param
  #   ols.ssr <- t(ols.mean.y - ols.predicted) %*% (ols.mean.y - ols.predicted)
  #   ols.denominador.r2 <- ols.sse + ols.ssr
  # }
  # ols.pseudo.r2 <- 1 - ols.sse/ols.denominador.r2
  
  ## Parâmetros do WLS -------------------------------------------------------
  
  # Retém somente os valores na diagonal principal da matriz de covariância
  wls.cov.matrix <- scovmatrix*diag(n)
  
  # Confere se gamma = 0 é suficiente (pode ser que o otimizador retorne < 0)
  std.delta <- 0
  
  wls.est <- fun.gls.var.delta(std.delta, n, wls.cov.matrix, skew, X, p) # resultadod c/ chute
  
  check.delta <- wls.est$check.delta
  check.delta0 <- check.delta
  std.delta <- wls.est$std.delta
  wls.par <- wls.est$b
  lambda <- wls.est$lambda
  inv.lambda <- wls.est$inv.lambda
  
  # Otimizar
  if(check.delta < 0){ std.delta <- 0 } else{
    
    # Definir outra função que retorna somente o check.delta p/ usar no uniroot (equivalente ao fzero)
    f.root <- function(g){fun.gls.var.delta(g, n, wls.cov.matrix, skew, X, p)$check.delta}
    
    wls.result.root <- uniroot(f = f.root, interval = c(0,20)) # minimizar fo
    std.delta <- wls.est$root                              # extrair 
    
    # Recalcular no novo delta
    wls.res.opt <- fun.gls.var.delta(std.delta, n, wls.cov.matrix, skew, X, p)
    
    check.gamma <- wls.est$check.delta
    std.delta <- wls.est$std.delta
    wls.par <- wls.est$b
    lambda <- wls.est$lambda
    inv.lambda <- wls.est$inv.lambda
    
  }
  
  # Covariância e erro padrão dos parâmetros do WLS
  wls.cov.par <- solve(t(X) %*% inv.lambda %*% X)
  wls.std.error.par <- (diag(wls.cov.par))^0.5
  wls.var.delta <- std.delta^2
  
  # Métricas de avaliação
  wls.sve <- X %*% wls.cov.par %*% t(X)
  wls.asve <- tr(wls.sve)/n
  wls.avp1 <- wls.var.delta + wls.asve
  wls.term <- 2*wls.var.delta*wls.sve*inv.lambda
  wls.avp2 <- wls.avp1 - tr(wls.term)/n
  
  # # Calcular R^2
  # wls.sse <- n*wls.var.delta
  # if(length(model.index) == 1){
  #   wls.mean.y <- mean(X) * wls.par %>% as.numeric
  #   wls.ssto <- t(skew - wls.mean.y) %*% (skew - wls.mean.y)
  # } else{
  #   wls.mean.y <- mean(X) %*% wls.par
  #   wls.ssr <- t(wls.mean.y - X %*% wls.par) %*% (wls.mean.y - X %*% wls.par)
  #   wls.denominador.r2 <- wls.sse + wls.ssr
  #   wls.r2 <- 1 - wls.sse/wls.denominador.r2
  # }
  
  
  ## Parâmetros do GLS -------------------------------------------------------
  
  
  
  

  ## Resultados --------------------------------------------------------------
  
  res.param.ols[model.index] <- ols.param
  res.se.param.ols[model.index] <- ols.se.param
  ols.result[kmodel,] <- c(t(res.param.ols), t(res.se.param.ols), ols.gamma, ols.asve, ols.avp1, ols.avp2)
  
  res.param.wls[model.index] <- wls.par
  res.se.param.wls[model.index] <- wls.std.error.par
  wls.result[kmodel,] <- c(t(res.param.wls), t(res.se.param.wls), wls.var.delta, wls.asve, wls.avp1, wls.avp2)
  
}

beepr::beep(sound = 10)

gc()
