

# PACOTES ----------------------------------------------------------------

pacman::p_load(pacman, tidyverse, R.matlab)


# CARREGAR TABELAS MATLAB -------------------------------------------------

# Diretório c/ arquivos MatLab
dir.matlab <- "C:/Users/tomas/OneDrive/1 - Acadêmico/Mestrado/Tese/MatLab_Modelos Regionais"

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

# Carregar individualmente as tabelas (somente os carregados no código MatLab)
models <- dfs.matlab$models9var$models # matriz 512x10 com todas as combinações possíveis de 9 variáveis (modelos)
scovmatrix <- dfs.matlab$ScCov0$ScCov0 # sampling covariance matrix 89x89
skew <- dfs.matlab$thetaSc$theta       # variável hidrológica de intereesse em cada estação (assimetria)


