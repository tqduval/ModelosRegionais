
# PACKAGES ----------------------------------------------------------------

if(!require("pacman")) install.packages("pacman")
pacman::p_load(pacman, tidyverse, ggplot2, extrafont, patchwork, sf)


# READ DATA ---------------------------------------------------------------

# Clear environment
rm(list = ls())

df_regional_kappa <- readRDS("Dados Gerados/df_regional_kappa.rds") # regional kappa models
df_wls_result <- as_tibble(readRDS("Dados Gerados/wls_result.rds")) # wls results
df_rs_daily <- readRDS("Dados Gerados/df_rs_daily.rds")             # daily rainfall data
df.fit <- readRDS("Dados Gerados/df.fit.rds")                       # GEV parameters
df.kappa.ci <- readRDS("Dados Gerados/df.kappa.ci.rds")             # GEV kappa CI
df.quantiles <- readRDS("Dados Gerados/df.quantiles.rds")           # GEV quantiles

# Read script that contains fun.priori.beta() function
source("C:/Users/tomas/OneDrive/1_Academico/Mestrado/Tese/3-R/Precipitacão Máxima Anual/PrecipitacaoMaximaAnual/0_Funcoes.r")

# Function for computing GEV quantiles
# This function is based on the one in "3_Quantis.r", but with regional quantiles added
source("C:/Users/Tomas/OneDrive/1_Academico/Mestrado/Tese/3-R/Precipitacão Máxima Anual/PrecipitacaoMaximaAnual/Funções/fun_gev_quantiles_regional.r")


# PRE-VIS OPERATIONS ------------------------------------------------------

# Empty matrix w/ output size
n.row <- nrow(df_wls_result)
n.par <- 3
n.col <- (n.par + 1)*2
wls.print <- list()

# Set up a matrix w/ parameters on the top row and se right below
for(i in 1:n.row){
  
  row <- unlist(df_wls_result[i, 1:n.col])         # get rows
  row.fold <- matrix(row, ncol = n.col/2, byrow = TRUE) # fold into 2 lines
  row.fold <- as.data.frame(row.fold)
  wls.print[[i]] <- row.fold
  
}; rm(row, row.fold)

# Join rows into matrix
wls.print.names <- names(df_wls_result)[1:(n.par + 1)]
wls.print <- bind_rows(wls.print)
colnames(wls.print) <- wls.print.names
is.list(wls.print)

# Pseudo-R²
var.delta_0 <- df_wls_result$var_delta[1] # extract \sigma_\delta^2(0)
df_wls_result$p_r2 <- 1 - df_wls_result$var_delta/var.delta_0

# Performance metrics w/ empty rows between
wls.metrics <- list()
for(i in 1:((n.par + 1)*2)){
  
  row <- df_wls_result[i, (n.col + 1):ncol(df_wls_result)]
  empty <- rep(NA, length(row))
  wls.metrics[[i]] <- rbind(row, empty)
  
}; rm(row, empty)

wls.metrics <- bind_rows(wls.metrics)

# Bind columns
df_wls_print <- bind_cols(wls.print, wls.metrics)
is_list <- sapply(df_wls_print, is.list)

# Print out df_wls_result
write.table(x = df_wls_print,
            file = "Dados Gerados/wls_par_results.csv",
            row.names = FALSE,
            sep = ";",
            dec = ",",
            fileEncoding = "UTF-8")

# Compute a relative difference between the variance of kappa estimates between with MLE and RGMLE and GMLE and RGMLE
# This difference will always have the original estimator (or basis of comparison) like
# RD = (MLE - RGMLR)/MLE and RD = (GMLE - RGMLE)/GMLE so that gains or differences can be proporly assessed

# The main focus needs to be in terms of uncertainty both in quantiles and in kappa estimates.

# Compute the uncertainty reduction on kappa and quantile estimates
# Preguiça de fazer sem ser usando pipes
# Kappa
df.kappa.sd.redct <- ({
  df.kappa.ci %>% 
    pivot_wider(id_cols = gauge.code,
                names_from = model,
                values_from = sd,
                names_prefix = "sd.") %>% 
    mutate(redct.gl.gev = (sd.l.gev - sd.gl.gev)/sd.l.gev,         # get sd reduction between RGMLE and MLE
           redct.gl.gev.r = (sd.l.gev - sd.gl.gev.r)/sd.l.gev) %>% # get sd reduction between RGMLE and GMLE
    pivot_longer(cols = starts_with("redct."),                 # pivot longer again
                 names_to = "base.model",
                 values_to = "sd.reduction") %>% 
    select(gauge.code, base.model, sd.reduction) %>% 
    mutate(base.model = sub("redct.", "", base.model))
})

# Investigate how kappa estimates vary for within the same station
df.kappa.rank <- ({
  df.kappa.ci %>% 
    pivot_wider(id_cols = gauge.code,
                names_from = model,
                values_from = par,
                names_prefix = "kappa.") %>% 
    mutate(kappa.sd.gg = pmap_dbl(.l = list(kappa.l.gev, kappa.gl.gev, kappa.gl.gev.r), # list of columns
                                  .f = ~ sd(c(...), na.rm = TRUE))) %>%                 # calculate sd between kappa values
    mutate(rank = rank(kappa.sd.gg, ties.method = "min")) %>%                           # new column for ranking smallest sd
    arrange(rank) %>% 
    select(gauge.code, rank, kappa.l.gev, kappa.gl.gev, kappa.gl.gev.r) %>% 
    left_join(df.fit %>% select(gauge.code, n.serie), by = "gauge.code") %>%            # bring series length information
    mutate(diff.kappa = kappa.gl.gev.r - kappa.l.gev,
           group = ifelse(diff.kappa < 0, "A", "B"),
           abs.diff = abs(diff.kappa))
})

# Quantiles
df.q.sd.redct <- ({
  df.quantiles %>% 
    pivot_wider(id_cols = c(gauge.code, tempo.retorno),
                names_from = model,
                values_from = sd.q,
                names_prefix = "sd.") %>% 
    mutate(redct.gl.gev = (sd.l.gev - sd.gl.gev)/sd.l.gev,         # get sd reduction between RGMLE and MLE
           redct.gl.gev.r = (sd.l.gev - sd.gl.gev.r)/sd.l.gev) %>% # get sd reduction between RGMLE and GMLE
    pivot_longer(cols = starts_with("redct."),
                 names_to = "base.model",
                 values_to = "sd.reduction") %>% 
    mutate(base.model = sub("redct.", "", base.model)) %>% 
    select(gauge.code, tempo.retorno, base.model, sd.reduction) %>% 
    left_join(df.kappa.rank %>% select(gauge.code, diff.kappa, group), by = "gauge.code")
})


# RECORD LENGTH DISTRIBUTION ----------------------------------------------

# Record length distribution
df.record.length <- 
  df_rs_daily %>% 
  group_by(gauge_code) %>% 
  summarise(n_years = n())

max.n <- max(df.record.length$n_years)
breaks <- c(seq(30, 100, 10), max.n)

# Plot
plot.n <- {(
  
  df.record.length %>% 
    ggplot(aes(x = n_years)) +
    geom_histogram(binwidth = 10,
                   breaks = breaks,
                   color = "steelblue",
                   linewidth = 0.5,
                   fill = "lightblue") +
    geom_vline(xintercept = max.n,
               color = "#D6A249",
               linewidth = 0.8) +
    labs(y = "N° estações", x = "N [anos]") +
    geom_text(stat = "bin",
              aes(label = after_stat(count)),
              vjust = -0.8,
              breaks = breaks,
              family = "serif") +
    annotate(geom = "text",
             x = 102, y = 10
             , angle = 90,
             label = "Estação mais longa do conjunto (N = 104)",
             hjust = "left", family = "serif", size = 4) +
    scale_x_continuous(breaks = seq(30, 110, 10)) +
    ylim(c(0, 55)) +
    theme_minimal() +
    theme(panel.background = NULL,
          legend.position =  "none",
          plot.background = element_rect(color = "white"),
          panel.border = element_rect(color = "black", fill = NA),
          axis.text = element_text(color = "black"),
          text = element_text(family = "serif", size = 12,))
  
)}; plot.n

# ggsave(plot = plot.n, filename = "Plotagens/plot.n.png",
       # width = 160, height = 110, units = "mm", dpi = 300)


# FREQUENCY PLOT SETUP ----------------------------------------------------

# Get station with the smallest difference between MLE and GMLE-R kappa estimates
gg.fq <- df.kappa.rank[df.kappa.rank[["abs.diff"]] == sort(df.kappa.rank[["abs.diff"]])[1], "gauge.code", drop = TRUE] # or use the dplyr::pull function or " %>% .[[]]"

# Make observed frequency table
df.fq.observed <-  df_rs_daily %>% 
  filter(gauge_code == gg.fq) %>% 
  select(gauge_code, peak_24h) %>% 
  arrange(desc(peak_24h)) %>%                                               # arrange series decreasing
  mutate(n = length(gauge_code), i = row_number(), p = i/(n + 1), tr = 1/p) # weibull plotting position

# Return period vector
tr <- c(1.1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 10, 20, 50, 100, 200, 500)

# Make estimated parameters table for station gg.fq
df.fq.par <- df.fit %>% 
  filter(gauge.code == gg.fq) %>% 
  select(gauge.code, rain.mm, n.serie, par.l.gev, par.gl.gev.r, hess.l.gev, hess.gl.gev.r)

# Use fun.gev.quantile.regional to estimate quantiles
df.fq.estimated <- fun.gev.quantile.regional(df = df.fq.par,
                                             tempo.retorno = tr,
                                             gauge.names = "gauge.code",
                                             model = c("l.gev", "gl.gev.r"))

# Contabilizar estações
df.q.sd.redct %>% filter(tempo.retorno == 50 & base.model == "gl.gev.r") %>%
  filter(diff.kappa < 0 & sd.reduction < 0) %>% nrow # 34
df.q.sd.redct %>% filter(tempo.retorno == 50 & base.model == "gl.gev.r") %>%
  filter(diff.kappa < 0 & sd.reduction > 0) %>% nrow # 51 
df.q.sd.redct %>% filter(tempo.retorno == 50 & base.model == "gl.gev.r") %>%
  filter(diff.kappa > 0 & sd.reduction > 0) %>% nrow # 100 

# Get percentual reduction of CI's for return periods of 30 and 100 yrs quantiles
df.ci.compare <- tibble(tr = c(20, 100), prct.ci.variation = rep(0, 2))
for(i in 1:nrow(df.ci.compare)){
  
  tr <- df.ci.compare[[i,1]]
  
  ic.l <- df.fq.estimated[["ic.l"]][df.fq.estimated[["tempo.retorno"]] == tr]
  ic.u <- df.fq.estimated[["ic.u"]][df.fq.estimated[["tempo.retorno"]] == tr]
  q <- df.fq.estimated[["q"]][df.fq.estimated[["tempo.retorno"]] == tr]
  
  width <- ic.u - ic.l
  df.ci.compare[i,2] <- (width[1] - width[2]) / width[1] * 100
  
}

# Choose one station to illustrate a prior distribution
gg <- "83942"
mu.kappa.r <- df_regional_kappa[df_regional_kappa[["gauges"]] == gg, "mu_kappa_r"] 
sd.kappa.r <- df_regional_kappa[df_regional_kappa[["gauges"]] == gg, "sd_kappa_r"] 


# ANALYSES ON WEIBULL GAUGES ----------------------------------------------

# If kappa > 0, GEV takes the form of the Weibull distribution with a fixed right limit
# Analyse how close the observed maximum of each station with kappa > 0 comes to this limit
df.weibull <- ({
  df.fit %>% 
    select(gauge.code, rain.mm, n.serie, par.gl.gev.r) %>% 
    mutate(csi = sapply(par.gl.gev.r, function(par) par[[1]]),   # extract csi
           alpha = sapply(par.gl.gev.r, function(par) par[[2]]), # extract alpha
           kappa = sapply(par.gl.gev.r, function(par) par[[3]]), # extract kappa
           ymax = sapply(rain.mm, function(y) max(y[[1]]))) %>%  # extract local ymax
    filter(kappa > 0) %>%                                        # filter out kappa < 0
    mutate(limit = csi + alpha/kappa,                            # theoretical Weibull pdf limit
           prct.limit = ymax/limit)                              # how close ymax gets to the limit
})

# Simple plot
hist(x = df.weibull$prct.limit)
boxplot(x = df.weibull$prct.limit)

font.size <- 9

# I. MAP OF THE GAUGES ----------------------------------------------------

sf_daily <- st_as_sf(df_rs_daily, coords = c("long", "lat"), crs = st_crs(4674))[, c("gauge_code", "geometry")]
sf_uf_path <- "C:/Users/tomas/OneDrive/7 - Engenharia/Programas/GIS/0-BASE/IBGE/BR_UF_2022/BR_UF_2022.shp"
sf_rs <- sf::st_read(dsn = sf_uf_path) %>% filter(SIGLA_UF == "RS")
# st_write(obj = sf_daily, dsn = "Dados Gerados/sf_daily.shp")

plot.gg.map <- ({
  ggplot() +
    geom_sf(data = sf_rs, fill = "white", color = "black", linewidth = 0.6) +
    geom_sf(data = sf_daily, fill = "steelblue2", color = "steelblue4", alpha = 0.8, shape = 21) +
    theme_minimal() +
    labs(x = "Longitude", y = "Latitude", fill = "", color = "") +
    theme(panel.background = NULL,
          legend.position =  "none",
          plot.background = element_rect(color = "white"),
          aspect.ratio = 1,
          panel.border = element_rect(color = "black", fill = NA),
          text = element_text(family = "serif", size = font.size))
})

# ggsave(filename = "Plotagens/plot.gg.map.png",
#        dpi = 300, width = 8, height = 8, unit = "cm")


# II. KAPPA ----------------------------------------------------

# Get some fonts
# extrafont::font_import()
# loadfonts(device = "win") # this takes some time

colors <- c("MLE" = "#cdc720", "GMLE-G" = "#12bce5", "GMLE-R" = "darkorange") # assign colors
model.names <- factor(names(colors), levels = c("MLE", "GMLE-G", "GMLE-R"))   # factor to stablish an order


plot.kappa <- ({
  
  # Boxplot of standard deviation of kappa estimates
  p1 <- df.kappa.ci %>%  
    mutate(model = factor(model, levels = c("l.gev", "gl.gev", "gl.gev.r")),
           model = recode(model,  "l.gev" = "MLE", "gl.gev" = "GMLE-G", "gl.gev.r" = "GMLE-R")) %>% # recode changes individual factors
    ggplot() +
    # geom_jitter(aes(x = model, y = sd), alpha = 0.2, size = 1.5, color = "grey70") +
    stat_boxplot(aes(x = model, y = sd), geom = "errorbar", width = 0.3) +
    geom_boxplot(aes(x = model, y = sd, fill = model), linewidth = 0.5, outlier.shape = 4, outlier.size = 2, show.legend = FALSE) +
    scale_fill_manual(values = colors) +
    scale_x_discrete(labels = model.names) +
    labs(x = "", y = expression("Desvio-padrão de" ~ hat(kappa))) +
    theme_minimal() +
    theme(legend.position =  "none",
          plot.background = element_rect(color = "white"),
          panel.border = element_rect(color = "black", fill = NA),
          text = element_text(family = "Times New Roman", color = "black", size = font.size)); p1
  
  # ggsave(filename = "Plotagens/separado/kappa.boxplot.png", plot = p1,
  #        width = 160, height = 110, dpi = 200, units = "mm")
  
  # Plot of reduction in standard deviation from MLE to the other two GML estimators
  p2 <- df.kappa.sd.redct %>% 
    mutate(base.model = factor(base.model, levels = c("gl.gev", "gl.gev.r")),
           base.model = recode(base.model,  "gl.gev" = "GMLE-G", "gl.gev.r" = "GMLE-R")) %>% # recode changes individual factors
    ggplot(aes(x = sd.reduction, color = base.model)) +
    stat_ecdf(geom = "smooth", linewidth = 0.7, show.legend = FALSE) +
    labs(x = "Redução na incerteza [%]", y = "Estações [%]", color = "") +
    coord_flip() +
    scale_color_manual(values = c("GMLE-G" = "#12bce5", "GMLE-R" = "darkorange")) +
    scale_x_continuous(labels = scales::percent) +
    scale_y_continuous(labels = scales::percent) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.background = element_rect(color = "white"),
          panel.border = element_rect(color = "black", fill = NA),
          text = element_text(family = "Times New Roman", color = "black", size = font.size)); p2
  
  # ggsave(filename = "Plotagens/separado/sigma.reduc.kappa.png", plot = p2,
  #        width = 160, height = 110, dpi = 200, units = "mm")
  
  # Prior distributions (geophysical and regional for one specific station)
  p3 <- ggplot() +
    # geom_blank(aes(color = "MLE"), inherit.aes = TRUE) +
    geom_hline(aes(yintercept = 0, color = "MLE"), linewidth = 0.7, alpha = 0) +
    stat_function(fun = fun.priori.beta, n = 1000, args = list(a = 0.5, mu = -0.10, sd = 0.122), linewidth = 0.7, aes(color = "GMLE-G")) +
    stat_function(fun = fun.priori.beta, n = 1000, args = list(a = 0.5, mu = mu.kappa.r, sd = sd.kappa.r), linewidth = 0.7, aes(color = "GMLE-R")) +
    scale_color_manual(values = colors,                                                 # use mapping from the colors object
                       breaks = model.names,                                            # breaks to order the elements
                       guide = guide_legend(override.aes = list(alpha = c(1, 1, 1)))) + # override.aes to make sure MLE appears in the legend
    scale_x_continuous(limits = c(-0.5, 0.5)) +
    labs(x = "κ", y = "π(κ)", color = "") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.background = element_rect(color = "white"),
          legend.margin = margin(0, 0, 0, 0),
          legend.justification = "center",
          panel.border = element_rect(color = "black", fill = NA),
          text = element_text(family = "Times New Roman", color = "black", size = font.size)); p3
  
  # ggsave(filename = "Plotagens/separado/prioris.kappa.png", plot = p3,
  #        width = 160, height = 110, dpi = 200, units = "mm")
  
  # Create final plot with patchwork
  (p1 + p2)/p3 +
    plot_layout(guides = "collect", heights = c(-1, -1)) &
    plot_annotation(tag_levels = "a", tag_suffix = ")") &
    theme(legend.position = "bottom")
  
}); plot.kappa

# ggsave(plot = plot.kappa, filename = "Plotagens/plot.kappa.png",
#        width = 160, height = 110, units = "mm", dpi = 300)


# III. QUANTILES ----------------------------------------------------------

colors.tr <- c("A" = "#bf4f10", "B" = "#ffc727") # set colors for plot
kappa.gg.fq <- c("l.gev" = sapply(df.fq.par$par.l.gev, function(p) p[[3]]), "gl.gev.r" = sapply(df.fq.par$par.gl.gev.r, function(p) p[[3]]))
annotate.font.size <- 3

plot.quantiles <- ({
  
  p4 <- df.q.sd.redct %>% 
    mutate(base.model = factor(base.model, levels = c("gl.gev", "gl.gev.r")),
           base.model = recode(base.model,  "gl.gev" = "GMLE-G", "gl.gev.r" = "GMLE-R")) %>%          # recode changes individual factors
    filter(tempo.retorno == 50 & base.model == "GMLE-R") %>%                                          # keep only T = 50 yrs
    arrange(sd.reduction) %>%                                                                         # arrange increasing
    filter(!row_number() %in% 1:3) %>%                                                                # remove first three observations that are just insane
    ggplot(aes(x = diff.kappa, y = sd.reduction, color = group)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = 0, linewidth = 0.2, linetype = "dashed", alpha = 0.8, color = "grey30") + # add axis centered on (0,0)
    geom_vline(xintercept = 0, linewidth = 0.2, linetype = "dashed", alpha = 0.8, color = "grey30") + # x axis
    annotate("curve", xend = -0.33, x = -0.2, yend = -0.8, y = -0.71,                                 # arrow indicating points outside plot limits
             arrow = arrow(type = "open", length = unit(1.8, "mm")),
             linewidth = 0.35, curvature = 0.1) +
    annotate("text", x = -0.21, y = -0.7, label = "02955006\n02753016\n03055005",                     # stations outside plot limits
             hjust = -0.3, family = "Times New Roman", size = unit(annotate.font.size, "mm")) +
    scale_color_manual(values = colors.tr,
                       labels = c("A" = expression(kappa[GMLE-R] < kappa[MLE]),                       # change the label values
                                  "B" = expression(kappa[GMLE-R] > kappa[MLE]))) +
    # scale_color_manual(values = colors.tr, labels = c("A" = "GMLE-R < MLE", "B" = "GMLE-R > MLE")) +
    scale_y_continuous(labels = scales::percent, limits = c(-0.8, 0.7)) +
    labs(x = "Diferença entre estimativas de κ", y = "Redução na incerteza [%]", color = "") +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.background = element_rect(color = "white"),
          panel.border = element_rect(color = "black", fill = NA),
          # aspect.ratio = 1,                                                                           # ensure square plot shape when plotting
          text = element_text(family = "Times New Roman", color = "black", size = font.size)); p4
  
  # ggsave(filename = "Plotagens/separado/sigma.q.reduc.png", plot = p4,
  #        width = 160, height = 110, dpi = 200, units = "mm")
  # 
  # Gráfico de frequência p/ estações cujo |diff.kappa| é próximo de 0
  p5 <- ggplot() +
    geom_ribbon(data = df.fq.estimated %>% filter(model == "l.gev"), alpha = 0.1, linewidth = 0.4,
                aes(x = tempo.retorno, ymin = ic.l, ymax = ic.u,
                    linetype = "IC 95% MLE", fill = "IC 95% MLE"), color = "black") +
    geom_ribbon(data = df.fq.estimated %>% filter(model == "gl.gev.r"), alpha = 0.3, linewidth = 0.4,
                aes(x = tempo.retorno, ymin = ic.l, ymax = ic.u,
                    linetype = "IC 95% GMLE-R", fill = "IC 95% GMLE-R"), color = "black") +
    geom_line(data = df.fq.estimated %>% filter(model == "gl.gev.r"),
              aes(x = tempo.retorno, y = q, linewidth = "Estimado"), color = "black") +
    geom_point(data = df.fq.observed,
               aes(x = tr, y = peak_24h, color = "Observado", shape = "Observado")) +
    scale_x_log10(breaks = c(1, 10, 20, 30, 50, 100, 200, 300, 500),
                  labels = c("1", "10", "20", "30", "50", "100", "200", "300", "500")) +
    geom_text(aes(x = 1, y = 265, label = "kappa[MLE]==0.0528"), parse = TRUE, # parse interprets expressions in labels
              hjust = -0.05, family = "Times New Roman", size = unit(annotate.font.size, "mm")) +
    geom_text(aes(x = 1, y = 254, label = "kappa[GMLE-R]==0.0538"), parse = TRUE,
              hjust = -0.04, family = "Times New Roman", size = unit(annotate.font.size, "mm")) +
    geom_text(aes(x = 1, y = 243, label = paste('Tamanho da série:', nrow(df.fq.observed), 'anos')),
              hjust = -0.03, family = "Times New Roman", size = unit(annotate.font.size, "mm")) +
    scale_linetype_manual(values = c("Estimado" = "solid", "IC 95% GMLE-R" = "solid", "IC 95% MLE" = "dashed")) +
    scale_linewidth_manual(values = c("Estimado" = 1)) +
    scale_color_manual(values = c("Observado" = "red")) +
    scale_shape_manual(values = c("Observado" = 8)) +
    scale_fill_manual(values = c("IC 95% GMLE-R" = "grey20", "IC 95% MLE" = "grey50")) +
    labs(x = "Tempo de retorno [anos]", y = "Precipitação diária máxima anual [mm]",
         color = "", shape = "", fill = "", linetype = "", linewidth = "") +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.background = element_rect(color = "white"),
          panel.border = element_rect(color = "black", fill = NA),
          text = element_text(family = "Times New Roman", color = "black", size = font.size)); p5
  
  ggsave(filename = "Plotagens/separado/freq.plot.png", plot = p5,
         width = 160, height = 110, dpi = 200, units = "mm")
  

  p4 + p5 +
    plot_annotation(tag_levels = "a", tag_suffix = ")") +
    plot_layout(ncol = 2, guides = "collect",) &           # change spacing between plot
    theme(text = element_text(family = "Times New Roman", color = "black", size = font.size),
          legend.position = "bottom",
          legend.box.spacing = unit(5, "mm"),
          legend.spacing.x = unit(2, "mm"),
          legend.key.size = unit(4, "mm"),
          legend.margin = margin(0, 0, 0, 0))                # remove plot margins
    
}); plot.quantiles

ggsave(plot = plot.quantiles, filename = "Plotagens/plot.quantiles.png",
       width = 160, height = 90, units = "mm", dpi = 300)
