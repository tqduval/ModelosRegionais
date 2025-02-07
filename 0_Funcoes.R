# Função p/ transformar as variáveis preditoras
fun.exp.var <- function(df, varexp){
  
  # Transformações na matriz de variáveis preditoras
  nr <- nrow(df)
  nc <- ncol(df)
  
  df[,1] <- 1
  df[,5:nc] <- log(df[,5:nc])
  
  # Centrralizar a matriz utilizando scale() p/ subtrair a média de cada coluna
  df[,5:nc] <- scale(df[,5:nc], center = TRUE, scale = FALSE)
  
  # Selecionar colunas especificadas em varexp
  x <- df[, varexp, drop = FALSE]
  
  return(x)
  
}

varexp.teste <- fun.exp.var(df = df.variaveis, varexp = c(1,2,3,4,5,6,7,8,9,10))
