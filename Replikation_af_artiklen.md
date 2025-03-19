Pakker
```{r}
library(esreg)
library(stats)
library(ggplot2)
library(dplyr)
```

VaR og ES
```{r}
VaR <- function(alpha = 0.975, sigma2 = 1, mu = 0, shift = 0, scale = 1, df = NULL , Norm = FALSE) {
  if (!Norm) {
      VaR <- qt(1 - alpha, df = df) * scale+shift
  } else {
      VaR <- sqrt((df - 2) / df) * qt(1 - alpha, df = df)*scale+shift
  }
  return(-VaR)
}

ES <- function(alpha = 0.975, sigma2 = 1, mu = 0, shift = 0, scale = 1, df = NULL , Norm = FALSE) {
  if (!Norm) {
      x <- qt(1 - alpha, df = df)
      ES <- -(dt(x, df) / (1 - alpha) * (df + x^2) / (df - 1)) * scale + shift
  } else {
      x <- qt(1 - alpha, df = df)
      ES <- -sqrt((df - 2) / df) * (dt(x, df) / (1 - alpha) * (df + x^2) / (df - 1))*scale+shift
  }
  return(-ES)
}
```

Shift og scale
```{r}
shift <- function(alpha, df, dfnull, Norm=FALSE) {
  VaR1 <- VaR(1-alpha, df = dfnull, Norm = Norm )
  VaR2 <- VaR(1-alpha, df = df, Norm = Norm )
  res <- (VaR1 - VaR2)
  return(res)
}

scale <- function(signi = 0.05, df = 100, Norm=FALSE) {
  res <- ES(alpha = 0.975, sigma2 = 1, mu = 0, df = df , Norm=Norm) / 
    ES(alpha = 1 - signi, sigma2 = 1, mu = 0, df = df , Norm=Norm)
  return(res)
}
```

Test 1
```{r}

Z_1 <- function(quan, sigma2, mu, T, M, scale=1, shift = 0, df = NULL, dfnull = 100, Norm=TRUE ) {
  Z_1_list <- c()
  var <- VaR(alpha = 0.975, sigma2 = sigma2, mu = mu, df = dfnull )
  es <- ES(alpha = 0.975, sigma2 = sigma2, mu = mu, df = dfnull )

  for (i in 1:M) {
    if(!Norm){
      X_t <- rt(T, df = df)+shift
    } else {
      X_t<-rt(T,df=df)*sqrt((df-2)/df)+shift
    }
    q <- 0
    P <- var + X_t
    N_t <- sum(P < 0)

    for (j in 1:T) {
      if (X_t[j] + var < 0) {
        q <- q + X_t[j] / es
      }
    }
    Z_1 <- q / N_t + 1
    Z_1_list <- c(Z_1_list, Z_1)
  }

  Z_1_list <- Z_1_list[!is.na(Z_1_list)]

  return(list(Z_1_list = Z_1_list))
}
```

Test 2
```{r}
Z_2 <- function(quan, sigma2, mu, T, M, shift = 0, scale = 1, df = NULL, dfnull = 100, Norm = FALSE ) {
  Z_2_list <- c()
  var <- VaR(alpha = 0.975, sigma2 = sigma2,mu = mu, df = dfnull , Norm = Norm)
  es <- ES(alpha = 0.975, sigma2 = sigma2,mu = mu, df = dfnull , Norm = Norm)

  for (i in 1:M) {
    if (!Norm) {
      X_t <- rt(T, df = df) * scale+shift
    } else {
      X_t <- rt(T, df = df) * sqrt((df - 2) / df)+shift
    }

    q <- 0
    Talpha <- T * 0.025

    for (j in 1:T) {
      if (X_t[j] + var < 0) {
        q <- q + X_t[j] / es
      }
    }
    Z_2 <- q / Talpha + 1
    Z_2_list <- c(Z_2_list, Z_2)
  }

  Z_2_list <- Z_2_list[!is.na(Z_2_list)]
  return(list(Z_2_list = Z_2_list))
}
```

Test 3
```{r}
Z_3 <- function(quan, T, M, df, dfnull, shift = 0, scale = 1, Norm = FALSE) {
  Z_3_list <- numeric(M)
  TAlpha <- round(T * 0.025)
  
  if (!Norm) {
    integral <- -T / TAlpha * integrate(
      function(p) pbeta(1 - p, shape1 = T - TAlpha, shape2 = TAlpha) * qt(p, df = dfnull),
      lower = 0, upper = 1
    )$value
  } else {
    integral <- -T / TAlpha * integrate(
      function(p) pbeta(1 - p, shape1 = T - TAlpha, shape2 = TAlpha) * qt(p, df = dfnull, ncp = 0) * sqrt((dfnull - 2) / dfnull),
      lower = 0, upper = 1
    )$value
  }
  
  for (j in 1:M) {
    U <- runif(T)
    PU <- if (!Norm) {
      qt(U, df = df) * scale+shift
    } else {
      qt(U, df = df) * sqrt((df - 2) / df)+shift
    }
    
    PU <- sort(PU)
    q <- sum(PU[1:TAlpha])
    ES_hat <- -q / TAlpha
    
    Z_3_list[j] <- -ES_hat / integral + 1
  }
  
  Z_3_list <- Z_3_list[!is.na(Z_3_list)]
  
  return(list(Z_3_list = Z_3_list))
}
```

VaR backtest
```{r}
VaR_backtest <- function(sigma2, mu, T, M, shift = 0, scale = 1, df = NULL, dfnull = 100 , Norm = FALSE, alpha=0.99) {
  var <- VaR(alpha = alpha, sigma2 = sigma2, mu = mu, df = dfnull , Norm = Norm)
  var_back <- numeric(M)
  
  if (!Norm) {
    for (j in 1:M) {
      X_t <- rt(T, df = df) * scale+shift
      exceedances <- sum(X_t + var < 0)
      var_back[j] <- exceedances
    }
  } else {
    for (j in 1:M) {
      X_t <- rt(T, df = df) * sqrt((df - 2) / df)+shift
      exceedances <- sum(X_t + var < 0)
      var_back[j] <- exceedances
    }
  }
  
  var_back <- var_back[!is.na(var_back)]
  return(list(var_back = var_back, mean_exceedances = mean(var_back)))
}

```

Power
```{r}
Power <- function(data_null, data_alternatives, significance_level, alternative_names = NULL) {
  sorted_null <- sort(data_null)
  crit_value <- quantile(sorted_null, probs = significance_level)

  powers <- list(Significance_Level = significance_level)
  for (i in seq_along(data_alternatives)) {
    power <- mean(data_alternatives[[i]] <= crit_value)
    powers[[alternative_names[i]]] <- power
  }

  return(as.data.frame(powers))
}
```

Til print af tabeller
```{r}
compute_table <- function(data_null, data_alt, significance_levels, df_value, source_name, power_labels) {
  lapply(significance_levels, function(sl) {
    Power(data_null, data_alt, sl, power_labels)
  }) %>%
    bind_rows() %>%
    mutate(H0_df = df_value, Source = source_name)
}
```

Til plots
```{r}
smoothed_cdf <- function(data, M = 100000) {
  density_data <- density(data, M = M, adjust = 3, from = min(data), to = max(data))
  x <- density_data$x
  y <- cumsum(density_data$y) / sum(density_data$y)
  return(data.frame(value = x, cdf = y))
}

plot_cdfs <- function(null_data, alternatives, critical_value1, critical_value2, pct_5 = NULL, pct_10 = NULL, alternative_names, plot_title = "Comparison of CDF (Null) and 1 - CDF (Alternatives)", xlim=NULL, ylim=NULL, colors = c("purple", "blue", "red", "green4"),h0_name="Null Hypo (CDF)") {
  cdf_null <- smoothed_cdf(null_data)
  cdf_null$group <- h0_name
  
  cdf_alternatives <- do.call(rbind, lapply(1:length(alternatives), function(i) {
    cdf_data <- smoothed_cdf(alternatives[[i]])
    cdf_data$cdf <- 1 - cdf_data$cdf
    cdf_data$group <- alternative_names[i]
    return(cdf_data)
  }))
  
  combined_data <- rbind(cdf_null, cdf_alternatives)
  
  p <- ggplot(combined_data, aes(x = value, y = cdf, color = group)) +
    geom_line(size = 1) +
    geom_vline(xintercept = critical_value1, linetype = "dashed", color = "black") +
    geom_vline(xintercept = critical_value2, linetype = "dashed", color = "black")+
    geom_vline(xintercept = pct_5, linetype = "dashed", color = "green4")+
    geom_vline(xintercept = pct_10, linetype = "dashed", color = "green4")+
    labs(
    title = plot_title,
    x = "Value",
    y = "Probability",
    color = "Group"
  ) +
    theme_minimal() +
    theme(legend.position = "bottom")+
    scale_color_manual(values = colors)
  
  if (!is.null(xlim)) {
    p <- p + scale_x_continuous(limits = xlim)
  }
  
  if (!is.null(ylim)) {
    p <- p + scale_y_continuous(limits = ylim)
  }
  
  print(p)
}

```

Standard antagelser
```{r}
T <- 250
M <- 100000
alpha <- 0.975
beta <- 0.99
```

Tabel 1 del 1:
```{r}
nu <- c(5, 100)
alpha_prime <- c(0.975, 0.95, 0.9)

results <- data.frame()

for (df in nu) {
  
  for (ap in alpha_prime) {
    
    gamma <- ES(alpha, df = df , Norm = FALSE)/ES(ap, df = df , Norm = FALSE)  
    
    ES_new <- ES(alpha, df = df , Norm = FALSE) * gamma
    
    VaR_new1pct <- VaR(alpha=beta, df = df , Norm = FALSE) * gamma
    
    results <- rbind(results, data.frame(
      nu = df,
      alpha_prime = ap,
      gamma = gamma,
      VaR_new = VaR_new1pct,
      ES_new = ES_new
    ))
  }
}

print(results)
```

Tabel 1 del 2:
```{r}
#scaled tests for df = 100 og df = 5 - Z2
Z2_tab1_100_25 <- Z_2(quan = 0.05, sigma2 = 1, mu = 0, T = 250, M = 100000, shift = 0, scale = scale(0.025, 100), df = 100, dfnull = 100, Norm = FALSE )
Z2_tab1_100_5 <- Z_2(quan = 0.05, sigma2 = 1, mu = 0, T = 250, M = 100000, shift = 0, scale = scale(0.05, 100), df = 100, dfnull = 100, Norm = FALSE )
Z2_tab1_100_10 <- Z_2(quan = 0.05, sigma2 = 1, mu = 0, T = 250, M = 100000, shift = 0, scale = scale(0.1, 100 ), df = 100, dfnull = 100, Norm = FALSE )

Z2_tab1_5_25 <- Z_2(quan = 0.05, sigma2 = 1, mu = 0, T = 250, M = 100000, shift = 0, scale = scale(0.025, 5 ), df = 5, dfnull = 5, Norm = FALSE )
Z2_tab1_5_5 <- Z_2(quan = 0.05, sigma2 = 1, mu = 0, T = 250, M = 100000, shift = 0, scale = scale(0.05, 5 ), df = 5, dfnull = 5, Norm = FALSE )
Z2_tab1_5_10 <- Z_2(quan = 0.05, sigma2 = 1, mu = 0, T = 250, M = 100000, shift = 0, scale = scale(0.1, 5 ), df = 5, dfnull = 5, Norm = FALSE )

#Powers for df = 5,10
names <- c("Power 5%", "Power 10%")
significance_levels5 <- c(0.041, 0.106)
significance_levels100 <- c(0.04, 0.108)

#Z2: H0 = 5 df
tabel1_5_Z2 <- compute_table(
  Z2_tab1_5_25$Z_2_list,
  list(Z2_tab1_5_5$Z_2_list, Z2_tab1_5_10$Z_2_list),
  significance_levels5, 5, "Z2", names
)

#Z2: H0 = 100 df
tabel1_100_Z2 <- compute_table(
  Z2_tab1_100_25$Z_2_list,
  list(Z2_tab1_100_5$Z_2_list, Z2_tab1_100_10$Z_2_list),
  significance_levels100, 100, "Z2", names
)

#scaled tests for df = 100 og df = 5 - Z3
Z3_tab1_100_25 <- Z_3(quan = 0.05, T = 250, M = 100000, shift = 0, scale = scale(0.025, 100 ), df = 100, dfnull = 100, Norm = FALSE)
Z3_tab1_100_5 <- Z_3(quan = 0.05, T = 250, M = 100000, shift = 0, scale = scale(0.05, 100 ), df = 100, dfnull = 100, Norm = FALSE)
Z3_tab1_100_10 <- Z_3(quan = 0.05, T = 250, M = 100000, shift = 0, scale = scale(0.1, 100 ), df = 100, dfnull = 100, Norm = FALSE)

Z3_tab1_5_25 <- Z_3(quan = 0.05, T = 250, M = 100000, shift = 0, scale = scale(0.025, 5 ), df = 5, dfnull = 5, Norm = FALSE)
Z3_tab1_5_5 <- Z_3(quan = 0.05, T = 250, M = 100000, shift = 0, scale = scale(0.05, 5 ), df = 5, dfnull = 5, Norm = FALSE)
Z3_tab1_5_10 <- Z_3(quan = 0.05, T = 250, M = 100000, shift = 0, scale = scale(0.1, 5 ), df = 5, dfnull = 5, Norm = FALSE)

#Z3: H0 = 5 df
tabel1_5_Z3 <- compute_table(
  Z3_tab1_5_25$Z_3_list,
  list(Z3_tab1_5_5$Z_3_list, Z3_tab1_5_10$Z_3_list),
  significance_levels5, 5, "Z3", names
)

#Z3: H0 = 100 df
tabel1_100_Z3 <- compute_table(
  Z3_tab1_100_25$Z_3_list,
  list(Z3_tab1_100_5$Z_3_list, Z3_tab1_100_10$Z_3_list),
  significance_levels100, 100, "Z3", names
)

#scaled tests for df = 100 og df = 5 - VaR
VaR_tab1_100_25 <- VaR_backtest(sigma2 = 1, mu = 0, T = 250, M = 100000, shift = 0, scale = scale(0.025, 100), df = 100, dfnull = 100, Norm = FALSE )
VaR_tab1_100_5 <- VaR_backtest(sigma2 = 1, mu = 0, T = 250, M = 100000, shift = 0, scale = scale(0.05, 100), df = 100, dfnull = 100, Norm = FALSE )
VaR_tab1_100_10 <- VaR_backtest(sigma2 = 1, mu = 0, T = 250, M = 100000, shift = 0, scale = scale(0.1, 100), df = 100, dfnull = 100, Norm = FALSE )

VaR_tab1_5_25 <- VaR_backtest(sigma2 = 1, mu = 0, T = 250, M = 100000, shift = 0, scale = scale(0.025, 5), df = 5, dfnull = 5, Norm = FALSE )
VaR_tab1_5_5 <- VaR_backtest(sigma2 = 1, mu = 0, T = 250, M = 100000, shift = 0, scale = scale(0.05, 5), df = 5, dfnull = 5, Norm = FALSE )
VaR_tab1_5_10 <- VaR_backtest(sigma2 = 1, mu = 0, T = 250, M = 100000, shift = 0, scale = scale(0.1, 5), df = 5, dfnull = 5, Norm = FALSE )

#VaR: H0 = 5 df
tabel1_5_VaR <- compute_table(
  -VaR_tab1_5_25$var_back,
  list(-VaR_tab1_5_5$var_back, -VaR_tab1_5_10$var_back),
  significance_levels5, 5, "VaR", names
)

#VaR: H0 = 100 df
tabel1_100_VaR <- compute_table(
  -VaR_tab1_100_25$var_back,
  list(-VaR_tab1_100_5$var_back, -VaR_tab1_100_10$var_back),
  significance_levels100, 100, "VaR", names
)

#Flet tabeller
tabel1 <- bind_rows(
  tabel1_5_Z2, tabel1_100_Z2, tabel1_5_Z3,
  tabel1_100_Z3, tabel1_5_VaR, tabel1_100_VaR
) %>%
  mutate(
    Significance_Level = Significance_Level * 100,
    Source = factor(Source, levels = c("Z2", "Z3", "VaR"))
  ) %>%
  arrange(H0_df, Significance_Level, Source)

print(tabel1)
```

Figur 2
```{r}
scaled_t_distribution <- function(df, gamma_values, x_range) {
  data <- data.frame()
  
  for (gamma in gamma_values) {
    x <- seq(x_range[1], x_range[2], length.out = 500)
    y <- dt(x / gamma, df = df) / gamma
    temp <- data.frame(x = x, y = y, gamma = factor(gamma))
    data <- rbind(data, temp)
  }
  
  return(data)
}

normal_distribution <- function(x_range) {
  x <- seq(x_range[1], x_range[2], length.out = 500)
  y <- dnorm(x) 
  return(data.frame(x = x, y = y, gamma = "Normal"))
}


df <- 100 
gamma_values <- c(1, 1.1, 1.2, 1.5)
x_range <- c(-6, 6) 

data_t <- scaled_t_distribution(df, gamma_values, x_range)
data_normal <- normal_distribution(x_range)
data_combined <- rbind(data_t, data_normal)

ggplot(data_combined, aes(x = x, y = y, color = gamma)) +
  geom_line(size = 1.2, aes(linetype = gamma)) +
  labs(
    title = expression(paste("Scaled ", italic("t"), "-distributions with Normal Distribution, df = 100")),
    x = NULL, y = NULL, color = expression(gamma), linetype = expression(gamma)
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  scale_color_manual(
    values = c("red", "purple", "blue", "green4", "black"),
    labels = c("γ = 1", "γ = 1.1", "γ = 1.2", "γ = 1.5", "Normal")
  ) +
  scale_linetype_manual(
    values = c("solid", "solid", "solid", "solid", "dotted"),
    labels = c("γ = 1", "γ = 1.1", "γ = 1.2", "γ = 1.5", "Normal")
  )
```

Figur 3 og 4
```{r}
# Figur 3 for Z2 (df = 100)
plot_cdfs(
  null_data = Z2_tab1_100_25$Z_2_list,
  alternatives = list(Z2_tab1_100_5$Z_2_list, Z2_tab1_100_10$Z_2_list),
  critical_value1 = quantile(Z2_tab1_100_25$Z_2_list, probs = 0.040),
  critical_value2 = quantile(Z2_tab1_100_25$Z_2_list, probs = 0.108),
  pct_5 = quantile(Z2_tab1_100_25$Z_2_list, probs = 0.05),
  pct_10 = quantile(Z2_tab1_100_25$Z_2_list, probs = 0.1),
  alternative_names = c("H1 alpha'=5%","H1 alpha'=10%"),
  plot_title = "Figure 3: CDF Plot for Z2 (df = 100)", 
  xlim=c(-4,1),h0_name="H0 alpha=2,5%"
)

# Figur 3 for Z3 (df = 100)
plot_cdfs(
  null_data = Z3_tab1_100_25$Z_3_list,
  alternatives = list(Z3_tab1_100_5$Z_3_list, Z3_tab1_100_10$Z_3_list),
  critical_value1 = quantile(Z3_tab1_100_25$Z_3_list, probs = 0.040),
  critical_value2 = quantile(Z3_tab1_100_25$Z_3_list, probs = 0.108),
  pct_5 = quantile(Z3_tab1_100_25$Z_3_list, probs = 0.05),
  pct_10 = quantile(Z3_tab1_100_25$Z_3_list, probs = 0.1),
  alternative_names = c("H1 alpha'=5%","H1 alpha'=10%"),
  plot_title = "Figure 3: CDF Plot for Z3 (df = 100)",
  xlim=c(-0.65, 0.35),h0_name="H0 alpha=2,5%"
)

# Figur 3 for VaR (df = 100)
plot_cdfs(
  null_data = -VaR_tab1_100_25$var_back,
  alternatives = list(-VaR_tab1_100_5$var_back, -VaR_tab1_100_10$var_back),
  critical_value1 = quantile(-VaR_tab1_100_25$var_back, probs = 0.041),
  critical_value2 = quantile(-VaR_tab1_100_25$var_back, probs = 0.106),
  pct_5 = NULL,
  pct_10 = NULL,
  alternative_names = c("H1 alpha'=5%","H1 alpha'=10%"),
  plot_title = "Figure 3: CDF Plot for VaR (df = 100)",
  xlim=c(-18,0),h0_name="H0 alpha=2,5%"
)

# Figur 4 for Z2 (df = 5)
plot_cdfs(
  null_data = Z2_tab1_5_25$Z_2_list,
  alternatives = list(Z2_tab1_5_5$Z_2_list, Z2_tab1_5_10$Z_2_list),
  critical_value1 = quantile(Z2_tab1_5_25$Z_2_list, probs = 0.041),
  critical_value2 = quantile(Z2_tab1_5_25$Z_2_list, probs = 0.106),
  pct_5 = quantile(Z2_tab1_5_25$Z_2_list, probs = 0.05),
  pct_10 = quantile(Z2_tab1_5_25$Z_2_list, probs = 0.1),
  alternative_names = c("H1 alpha'=5%","H1 alpha'=10%"),
  plot_title = "Figure 4: CDF Plot for Z2 (df = 5)",
  xlim=c(-4.25,1),h0_name="H0 alpha=2,5%"
)

# Figur 4 for Z3 (df = 5)
plot_cdfs(
  null_data = Z3_tab1_5_25$Z_3_list,
  alternatives = list(Z3_tab1_5_5$Z_3_list, Z3_tab1_5_10$Z_3_list),
  critical_value1 = quantile(Z3_tab1_5_25$Z_3_list, probs = 0.041),
  critical_value2 = quantile(Z3_tab1_5_25$Z_3_list, probs = 0.106),
  pct_5 = quantile(Z3_tab1_5_25$Z_3_list, probs = 0.05),
  pct_10 = quantile(Z3_tab1_5_25$Z_3_list, probs = 0.1),
  alternative_names = c("H1 alpha'=5%","H1 alpha'=10%"),
  plot_title = "Figure 4: CDF Plot for Z3 (df = 5)",
  xlim=c(-1.25, 0.3),h0_name="H0 alpha=2,5%"
)

# Figur 4 for VaR (df = 5)
plot_cdfs(
  null_data = -VaR_tab1_5_25$var_back,
  alternatives = list(-VaR_tab1_5_5$var_back, -VaR_tab1_5_10$var_back),
  critical_value1 = quantile(-VaR_tab1_5_25$var_back, probs = 0.041),
  critical_value2 = quantile(-VaR_tab1_5_25$var_back, probs = 0.106),
  pct_5 = NULL,
  pct_10 = NULL,
  alternative_names = c("H1 alpha'=5%","H1 alpha'=10%"),
  plot_title = "Figure 4: CDF Plot for VaR (df = 5)",
  xlim=c(-17,0),h0_name="H0 alpha=2,5%"
)
```

Tabel 2 del 1:
```{r}
new_nu <- c(100,10,5,3)

results <- data.frame()

for (df in new_nu) {
  
    ES_2.5pct <- ES(alpha, df = df , Norm = FALSE)
    
    VaR_1pct <- VaR(alpha=beta, df = df , Norm = FALSE)
    
    ES_2.5pct_n <- ES(alpha, df = df , Norm = TRUE)
  
    VaR_1pct_n <- VaR(alpha=beta, df = df , Norm = TRUE)
    
    results <- rbind(results, data.frame(
      nu = df,
      VaR_1pct = VaR_1pct,
      ES_2.5pct = ES_2.5pct,
      VaR_1pct_n=VaR_1pct_n,
      ES_2.5pct_n=ES_2.5pct_n
    ))
  }

print(results)
```

Tabel 2 del 2 - std student t
```{r}
#Z2 df 100:
Z2_tab2_100_100 <-Z_2(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = 0 , scale  = 1, df = 100, dfnull = 100, Norm = FALSE )
Z2_tab2_100_10<-Z_2(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = 0 , scale  = 1, df = 10,  dfnull = 100, Norm = FALSE )
Z2_tab2_100_5<-Z_2(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = 0 , scale  = 1, df = 5,  dfnull = 100, Norm = FALSE )
Z2_tab2_100_3<-Z_2(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = 0 , scale  = 1, df = 3,   dfnull = 100, Norm = FALSE )

#Z2 df 10:
Z2_tab2_10_10 <-Z_2(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = 0 , scale  = 1, df = 10, dfnull = 10, Norm = FALSE )
Z2_tab2_10_5 <-Z_2(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = 0 , scale  = 1, df = 5,  dfnull = 10, Norm = FALSE )
Z2_tab2_10_3 <-Z_2(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = 0 , scale  = 1, df = 3,  dfnull = 10, Norm = FALSE )


#Z3 df 100:
Z3_tab2_100_100 <-Z_3(quan = 0.05, T=250, M = 100000, df=100, dfnull = 100, shift = 0, scale = 1)
Z3_tab2_100_10<-Z_3(quan = 0.05, T=250, M = 100000, df=10,  dfnull = 100, shift = 0, scale = 1)
Z3_tab2_100_5<-Z_3(quan = 0.05, T=250, M = 100000, df=5,  dfnull = 100, shift = 0, scale = 1)
Z3_tab2_100_3 <-Z_3(quan = 0.05, T=250, M = 100000, df=3,   dfnull = 100, shift = 0, scale = 1)

#Z3 df 10:
Z3_tab2_10_10<-Z_3(quan = 0.05, T=250, M = 100000, df=10,  dfnull = 10, shift = 0, scale = 1)
Z3_tab2_10_5 <-Z_3(quan = 0.05, T=250, M = 100000, df=5,   dfnull = 10, shift = 0, scale = 1)
Z3_tab2_10_3 <-Z_3(quan = 0.05, T=250, M = 100000, df=3,   dfnull = 10, shift = 0, scale = 1)

#VaR df 100:
VaR_tab2_100_100 <-VaR_backtest( sigma2= 1 , mu = 0, T = 250, M = 100000, shift = 0 , scale  = 1, df = 100, dfnull = 100, Norm = FALSE )
VaR_tab2_100_10 <-VaR_backtest( sigma2= 1 , mu = 0, T = 250, M = 100000, shift = 0 , scale  = 1, df = 10,  dfnull = 100, Norm = FALSE )
VaR_tab2_100_5 <-VaR_backtest( sigma2= 1 , mu = 0, T = 250, M = 100000, shift = 0 , scale  = 1, df = 5,  dfnull = 100, Norm = FALSE )
VaR_tab2_100_3 <-VaR_backtest( sigma2= 1 , mu = 0, T = 250, M = 100000, shift = 0 , scale  = 1, df = 3,   dfnull = 100, Norm = FALSE )

#VaR df 10:
VaR_tab2_10_10 <-VaR_backtest(sigma2= 1 , mu = 0, T = 250, M = 100000, shift = 0 , scale  = 1, df = 10, dfnull = 10, Norm = FALSE )
VaR_tab2_10_5 <-VaR_backtest(sigma2= 1 , mu = 0, T = 250, M = 100000, shift = 0 , scale  = 1, df = 5,  dfnull = 10, Norm = FALSE )
VaR_tab2_10_3<-VaR_backtest(sigma2= 1 , mu = 0, T = 250, M = 100000, shift = 0 , scale  = 1, df = 3,  dfnull = 10, Norm = FALSE )
```

Printer tabel
```{r}
names100 <- c("Powerdf10", "Powerdf3")
names10 <- c("Powerdf5", "Powerdf3")
significance_levels100_std <- c(0.041, 0.104)
significance_levels10_std <- c(0.04, 0.106)

#H0 = df 100
tabel2_100_Z2 <- compute_table(
  Z2_tab2_100_100$Z_2_list,
  list(Z2_tab2_100_10$Z_2_list, Z2_tab2_100_3$Z_2_list),
  significance_levels100_std, 100, "Z2", names100
)

tabel2_100_Z3 <- compute_table(
  Z3_tab2_100_100$Z_3_list,
  list(Z3_tab2_100_10$Z_3_list, Z3_tab2_100_3$Z_3_list),
  significance_levels100_std, 100, "Z3", names100
)

tabel2_100_VaR <- compute_table(
  -VaR_tab2_100_100$var_back,
  list(-VaR_tab2_100_10$var_back, -VaR_tab2_100_3$var_back),
  significance_levels100_std, 100, "VaR", names100
)

#H0 = df 10
tabel2_10_Z2 <- compute_table(
  Z2_tab2_10_10$Z_2_list,
  list(Z2_tab2_10_5$Z_2_list, Z2_tab2_10_3$Z_2_list),
  significance_levels10_std, 10, "Z2", names10
)

tabel2_10_Z3 <- compute_table(
  Z3_tab2_10_10$Z_3_list,
  list(Z3_tab2_10_5$Z_3_list, Z3_tab2_10_3$Z_3_list),
  significance_levels10_std, 10, "Z3", names10
)

tabel2_10_VaR <- compute_table(
  -VaR_tab2_10_10$var_back,
  list(-VaR_tab2_10_5$var_back, -VaR_tab2_10_3$var_back),
  significance_levels10_std, 10, "VaR", names10
)

tabel2_std <- bind_rows(
  tabel2_100_Z2, tabel2_100_Z3, tabel2_100_VaR,
  tabel2_10_Z2, tabel2_10_Z3, tabel2_10_VaR
) %>%
  mutate(
    Significance_Level = Significance_Level * 100,
    Source = factor(Source, levels = c("Z2", "Z3", "VaR"))
  ) %>%
  arrange(H0_df, Significance_Level, Source)

tabel2_std <- tabel2_std[, c("Significance_Level", "Powerdf3", "Powerdf5", "Powerdf10", "H0_df", "Source")]

print(tabel2_std)
```

Tabel 2 del 2 - normalized student t
```{r}
#Z2: df = 100
Z2_tab2_100_100_n <-Z_2(quan = 0.05, sigma2= 1 , mu = 0, T = 250, n = 100000, shift = 0 , scale  = 1, df = 100, dfnull = 100, Norm = TRUE )
Z2_tab2_100_10_n<-Z_2(quan = 0.05, sigma2= 1 , mu = 0, T = 250, n = 100000, shift = 0 , scale  = 1, df = 10,  dfnull = 100, Norm = TRUE )
Z2_tab2_100_5_n <-Z_2(quan = 0.05, sigma2= 1 , mu = 0, T = 250, n = 100000, shift = 0 , scale  = 1, df = 5, dfnull = 100, Norm = TRUE )
Z2_tab2_100_3_n<-Z_2(quan = 0.05, sigma2= 1 , mu = 0, T = 250, n = 100000, shift = 0 , scale  = 1, df = 3,   dfnull = 100, Norm = TRUE )

#Z3: df = 100
Z3_tab2_100_100_n<-Z_3(quan = 0.05, T=250, n = 100000, df=100,  dfnull = 100, shift = 0, scale = 1, Norm=TRUE)
Z3_tab2_100_10_n <-Z_3(quan = 0.05, T=250, n = 100000, df=10,   dfnull = 100, shift = 0, scale = 1, Norm=TRUE)
Z3_tab2_100_5_n <-Z_3(quan = 0.05, T=250, n = 100000, df=5,   dfnull = 100, shift = 0, scale = 1, Norm=TRUE)
Z3_tab2_100_3_n <-Z_3(quan = 0.05, T=250, n = 100000, df=3,   dfnull = 100, shift = 0, scale = 1, Norm=TRUE)

#VaR: df = 100
VaR_tab2_100_100_n <-VaR_backtest( sigma2 = 1 , mu = 0, T = 250, n = 100000, shift = 0 , scale  = 1, df = 100, dfnull = 100, Norm = TRUE )
VaR_tab2_100_10_n <-VaR_backtest( sigma2 = 1 , mu = 0, T = 250, n = 100000, shift = 0 , scale  = 1, df = 10,  dfnull = 100, Norm = TRUE )
VaR_tab2_100_5_n <-VaR_backtest( sigma2 = 1 , mu = 0, T = 250, n = 100000, shift = 0 , scale  = 1, df = 5,  dfnull = 100, Norm = TRUE )
VaR_tab2_100_3_n <-VaR_backtest( sigma2 = 1 , mu = 0, T = 250, n = 100000, shift = 0 , scale  = 1, df = 3,   dfnull = 100, Norm = TRUE )

#Z2: df = 10
Z2_tab2_10_10_n <-Z_2(quan = 0.05, sigma2= 1 , mu = 0, T = 250, n = 100000, shift = 0 , scale  = 1, df = 10, dfnull = 10, Norm = TRUE )
Z2_tab2_10_5_n <-Z_2(quan = 0.05, sigma2= 1 , mu = 0, T = 250, n = 100000, shift = 0 , scale  = 1, df = 5,  dfnull = 10, Norm = TRUE )
Z2_tab2_10_3_n <-Z_2(quan = 0.05, sigma2= 1 , mu = 0, T = 250, n = 100000, shift = 0 , scale  = 1, df = 3,  dfnull = 10, Norm = TRUE )

#Z3: df = 10
Z3_tab2_10_10_n<-Z_3(quan = 0.05, T=250, n = 100000, df=10,  dfnull = 10, shift = 0, scale = 1, Norm=TRUE)
Z3_tab2_10_5_n <-Z_3(quan = 0.05, T=250, n = 100000, df=5,   dfnull = 10, shift = 0, scale = 1, Norm=TRUE)
Z3_tab2_10_3_n <-Z_3(quan = 0.05, T=250, n = 100000, df=3,   dfnull = 10, shift = 0, scale = 1, Norm=TRUE)

#VaR: df = 10
VaR_tab2_10_10_n <-VaR_backtest( sigma2= 1 , mu = 0, T = 250, n = 100000, shift = 0 , scale  = 1, df = 10, dfnull = 10, Norm = TRUE )
VaR_tab2_10_5_n <-VaR_backtest( sigma2= 1 , mu = 0, T = 250, n = 100000, shift = 0 , scale  = 1, df = 5, dfnull = 10, Norm = TRUE )
VaR_tab2_10_3_n <-VaR_backtest( sigma2= 1 , mu = 0, T = 250, n = 100000, shift = 0 , scale  = 1, df = 3, dfnull = 10, Norm = TRUE )

significance_levels100_norm <- c(0.044, 0.11)
significance_levels10_norm <- c(0.044, 0.112)
significance_levels100_VaR_norm <- c(0.041, 0.104)
significance_levels10_VaR_norm <- c(0.04, 0.106)

#H_0: df = 100
tabel2_100_Z2_norm <- compute_table(
  Z2_tab2_100_100_n$Z_2_list,
  list(Z2_tab2_100_10_n$Z_2_list, Z2_tab2_100_3_n$Z_2_list),
  significance_levels100_norm, 100, "Z2", names100
)

tabel2_100_Z3_norm <- compute_table(
  Z3_tab2_100_100_n$Z_3_list,
  list(Z3_tab2_100_10_n$Z_3_list, Z3_tab2_100_3_n$Z_3_list),
  significance_levels100_norm, 100, "Z3", names100
)

tabel2_100_VaR_norm <- compute_table(
  (-VaR_tab2_100_100_n$var_back),
  list((-VaR_tab2_100_10_n$var_back), (-VaR_tab2_100_3_n$var_back)),
  significance_levels100_VaR_norm, 100, "VaR", names100
)

#H0: df = 10
tabel2_10_Z2_norm <- compute_table(
  Z2_tab2_10_10_n$Z_2_list,
  list(Z2_tab2_10_5_n$Z_2_list, Z2_tab2_10_3_n$Z_2_list),
  significance_levels10_norm, 10, "Z2", names10
)

tabel2_10_Z3_norm <- compute_table(
  Z3_tab2_10_10_n$Z_3_list,
  list(Z3_tab2_10_5_n$Z_3_list, Z3_tab2_10_3_n$Z_3_list),
  significance_levels10_norm, 10, "Z3", names10
)

tabel2_10_VaR_norm <- compute_table(
  (-VaR_tab2_10_10_n$var_back),
  list((-VaR_tab2_10_5_n$var_back), (-VaR_tab2_10_3_n$var_back)),
  significance_levels10_VaR_norm, 10, "VaR", names10
)


tabel2_norm <- bind_rows(
  tabel2_100_Z2_norm, tabel2_100_Z3_norm, tabel2_100_VaR_norm,
  tabel2_10_Z2_norm, tabel2_10_Z3_norm, tabel2_10_VaR_norm
) %>%
  mutate(
    Significance_Level = Significance_Level * 100,
    Source = factor(Source, levels = c("Z2", "Z3", "VaR"))
  ) %>%
  arrange(H0_df, Significance_Level, Source)

tabel2_norm <- tabel2_norm[, c("Significance_Level", "Powerdf3", "Powerdf5", "Powerdf10", "H0_df", "Source")]

print(tabel2_norm)
```

Figur 5
```{r}
generate_t_distributions <- function(df_values, x_range) {
  data <- data.frame()
  
  for (df in df_values) {
    x <- seq(x_range[1], x_range[2], length.out = 500)
    y <- dt(x, df = df) 
    temp <- data.frame(x = x, y = y, df = factor(df))
    data <- rbind(data, temp)
  }
  
  return(data)
}
normal_distribution <- function(x_range) {
  x <- seq(x_range[1], x_range[2], length.out = 500)
  y <- dnorm(x) 
  return(data.frame(x = x, y = y, df = "Normal"))
}

df_values <- c(100, 10, 5, 3)
x_range <- c(-6, 6)

data_t <- generate_t_distributions(df_values, x_range)
data_normal <- normal_distribution(x_range)
data_combined <- rbind(data_t, data_normal)

ggplot(data_combined, aes(x = x, y = y, color = df,linetype = df)) +
  geom_line(size = 0.8) +
  labs(
    title = "",
    x = NULL, y = NULL,
    color = expression(df), linetype = expression(df)
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  scale_color_manual(
    values = c("red", "purple", "blue", "green4", "black"),
    labels = c("100", "10", "5", "3", "Normal")
  ) +
  scale_linetype_manual(
    values = c("solid", "solid", "solid", "solid", "dotted"),
    labels = c("100", "10", "5", "3", "Normal")
  )
```

Figur 6
```{r}
# Figur 6 for Z2 (df = 100)
plot_cdfs(
  null_data = Z2_tab2_100_100$Z_2_list,
  alternatives = list(Z2_tab2_100_10$Z_2_list,Z2_tab2_100_5$Z_2_list, Z2_tab2_100_3$Z_2_list),
  critical_value1 = quantile(Z2_tab2_100_100$Z_2_list, probs = 0.041),
  critical_value2 = quantile(Z2_tab2_100_100$Z_2_list, probs = 0.104),
  pct_5 = quantile(Z2_tab2_100_100$Z_2_list, probs = 0.05),
  pct_10 = quantile(Z2_tab2_100_100$Z_2_list, probs = 0.1),
  alternative_names = c("H1 df=10", "H1 df=5", "H1 df=3"),
  plot_title = "Figure 6: Z2, df = 100",
  xlim = c(-6.2, 1),h0_name="H0 df=100"
)

# Figur 6 for Z3 (df = 100)
plot_cdfs(
  null_data = Z3_tab2_100_100$Z_3_list,
  alternatives = list(Z3_tab2_100_10$Z_3_list,Z3_tab2_100_5$Z_3_list, Z3_tab2_100_3$Z_3_list),
  critical_value1 = quantile(Z3_tab2_100_100$Z_3_list, probs = 0.041),
  critical_value2 = quantile(Z3_tab2_100_100$Z_3_list, probs = 0.104),
  pct_5 = quantile(Z3_tab2_100_100$Z_3_list, probs = 0.05),
  pct_10 = quantile(Z3_tab2_100_100$Z_3_list, probs = 0.1),
  alternative_names = c("H1 df=10", "H1 df=5", "H1 df=3"),
  plot_title = "Figure 6: Z3, df = 100", 
  xlim = c(-3.5, 0.5),h0_name="H0 df=100"
)

# Figur 6 for VaR (df = 100)
plot_cdfs(
  null_data = -VaR_tab2_100_100$var_back,
  alternatives = list(-VaR_tab2_100_10$var_back, -VaR_tab2_100_5$var_back, -VaR_tab2_100_3$var_back),
  critical_value1 = quantile(-VaR_tab2_100_100$var_back, probs = 0.041),
  critical_value2 = quantile(-VaR_tab2_100_100$var_back, probs = 0.104),
  pct_5 = NULL,
  pct_10 =NULL,
  alternative_names = c("H1 df=10", "H1 df=5", "H1 df=3"),
  plot_title = "Figure 6: VaR, df = 100",h0_name="H0 df=100"
)
```

Figur 7
```{r}
# Figur 7 for Z2 (df = 10)
plot_cdfs(
  null_data = Z2_tab2_10_10$Z_2_list,
  alternatives = list(Z2_tab2_10_5$Z_2_list, Z2_tab2_10_3$Z_2_list),
  critical_value1 = quantile(Z2_tab2_10_10$Z_2_list, probs = 0.04),
  critical_value2 = quantile(Z2_tab2_10_10$Z_2_list, probs = 0.106),
  pct_5 = quantile(Z2_tab2_10_10$Z_2_list, probs = 0.05),
  pct_10 = quantile(Z2_tab2_10_10$Z_2_list, probs = 0.1),
  alternative_names = c("H1 df=5", "H1 df=3"),
  plot_title = "Figure 7: Z2, df = 10",
  xlim=c(-5,1),h0_name="H0 df=10"
)

# Figur 7 for Z3 (df = 10)
plot_cdfs(
  null_data = Z3_tab2_10_10$Z_3_list,
  alternatives = list(Z3_tab2_10_5$Z_3_list, Z3_tab2_10_3$Z_3_list),
  critical_value1 = quantile(Z3_tab2_10_10$Z_3_list, probs = 0.04),
  critical_value2 = quantile(Z3_tab2_10_10$Z_3_list, probs = 0.106),
  pct_5 = quantile(Z3_tab2_10_10$Z_3_list, probs = 0.05),
  pct_10 = quantile(Z3_tab2_10_10$Z_3_list, probs = 0.1),
  alternative_names = c("H1 df=5", "H1 df=3"),
  plot_title = "Figure 7: Z3, df = 10",
  xlim=c(-2.5, 0.25),h0_name="H0 df=10"
)

# Figur 7 for VaR (df = 10)
plot_cdfs(
  null_data = -VaR_tab2_10_10$var_back,
  alternatives = list(-VaR_tab2_10_5$var_back, -VaR_tab2_10_3$var_back),
  critical_value1 = quantile(-VaR_tab2_10_10$var_back, probs = 0.040),
  critical_value2 = quantile(-VaR_tab2_10_10$var_back, probs = 0.106),
  pct_5 = NULL,
  pct_10 =NULL,
  alternative_names = c("H1 df=5", "H1 df=3"),
  plot_title = "Figure 7: VaR, df = 10",h0_name="H0 df=10"
)
```

Figur 8
```{r}
generate_normalized_t_distributions <- function(df_values, x_range) {
  data <- data.frame()
  
  for (df in df_values) {
    x <- seq(x_range[1], x_range[2], length.out = 500)
    scale_factor <- sqrt((df - 2) / df) 
    y <- dt(x / scale_factor, df = df) / scale_factor
    temp <- data.frame(x = x, y = y, df = factor(df))
    data <- rbind(data, temp)
  }
  
  return(data)
}

df_values <- c(100, 10, 5, 3)
x_range <- c(-6, 6)

data_normalized_t <- generate_normalized_t_distributions(df_values, x_range)

ggplot(data_normalized_t, aes(x = x, y = y, color = df)) +
  geom_line(size = 0.8) +
  labs(
    title = "Normalized t-distributions",
    x = NULL, y = NULL,
    color = "df"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  scale_color_manual(
    values = c("red", "purple", "blue", "green4")
  )

```

Figur 9
```{r}
# Figur 9 for Z2 (df = 100)
plot_cdfs(
  null_data = Z2_tab2_100_100_n$Z_2_list,
  alternatives = list(Z2_tab2_100_10_n$Z_2_list,Z2_tab2_100_5_n$Z_2_list, Z2_tab2_100_3_n$Z_2_list),
  critical_value1 = quantile(Z2_tab2_100_100_n$Z_2_list, probs = 0.044),
  critical_value2 = quantile(Z2_tab2_100_100_n$Z_2_list, probs = 0.11),
  pct_5 = quantile(Z2_tab2_100_100_n$Z_2_list, probs = 0.05),
  pct_10 = quantile(Z2_tab2_100_100_n$Z_2_list, probs = 0.1),
  alternative_names = c("H1 df=10", "H1 df=5", "H1 df=3"),
  plot_title = "Figure 9: Norm Z2, df = 100",
  xlim=c(-2, 1),h0_name="H0 df=100"
)

# Figur 9 for Z3 (df = 100)
plot_cdfs(
  null_data = Z3_tab2_100_100_n$Z_3_list,
  alternatives = list(Z3_tab2_100_10_n$Z_3_list,Z3_tab2_100_5_n$Z_3_list, Z3_tab2_100_3_n$Z_3_list),
  critical_value1 = quantile(Z3_tab2_100_100_n$Z_3_list, probs = 0.044),
  critical_value2 = quantile(Z3_tab2_100_100_n$Z_3_list, probs = 0.11),
  pct_5 = quantile(Z3_tab2_100_100_n$Z_3_list, probs = 0.05),
  pct_10 = quantile(Z3_tab2_100_100_n$Z_3_list, probs = 0.1),
  alternative_names = c("H1 df=10", "H1 df=5", "H1 df=3"),
  plot_title = "Figure 9: Norm Z3, df = 100", 
  xlim=c(-1.4, 0.3),h0_name="H0 df=100"
)

# Figur 9 for VaR (df = 100)
plot_cdfs(
  null_data = -VaR_tab2_100_100_n$var_back,
  alternatives = list(-VaR_tab2_100_10_n$var_back,-VaR_tab2_100_5_n$var_back, -VaR_tab2_100_3_n$var_back),
  critical_value1 = quantile(-VaR_tab2_100_100_n$var_back, probs = 0.041),
  critical_value2 = quantile(-VaR_tab2_100_100_n$var_back, probs = 0.106),
  pct_5 = NULL,
  pct_10 = NULL,
  alternative_names = c("H1 df=10", "H1 df=5", "H1 df=3"),
  plot_title = "Figure 9: Norm VaR, df = 100",h0_name="H0 df=100"
)
```

Figur 10
```{r}
# Figur 10 for Z2 (df = 10)
plot_cdfs(
  null_data = Z2_tab2_10_10_n$Z_2_list,
  alternatives = list(Z2_tab2_10_5_n$Z_2_list, Z2_tab2_10_3_n$Z_2_list),
  critical_value1 = quantile(Z2_tab2_10_10_n$Z_2_list, probs = 0.044),
  critical_value2 = quantile(Z2_tab2_10_10_n$Z_2_list, probs = 0.112),
  pct_5 = quantile(Z2_tab2_10_10_n$Z_2_list, probs = 0.05),
  pct_10 = quantile(Z2_tab2_10_10_n$Z_2_list, probs = 0.1),
  alternative_names = c("H1 df=5", "H1 df=3"),
  plot_title = "Figure 10: Norm Z2, df = 10", 
  xlim=c(-1.5, 0.75),h0_name="H0 df=10"
)

# Figur 10 for Z3 (df = 10)
plot_cdfs(
  null_data = Z3_tab2_10_10_n$Z_3_list,
  alternatives = list(Z3_tab2_10_5_n$Z_3_list, Z3_tab2_10_3_n$Z_3_list),
  critical_value1 = quantile(Z3_tab2_10_10_n$Z_3_list, probs = 0.044),
  critical_value2 = quantile(Z3_tab2_10_10_n$Z_3_list, probs = 0.112),
  pct_5 = quantile(Z3_tab2_10_10_n$Z_3_list, probs = 0.05),
  pct_10 = quantile(Z3_tab2_10_10_n$Z_3_list, probs = 0.1),
  alternative_names = c("H1 df=5", "H1 df=3"),
  plot_title = "Figure 10: Norm Z3, df = 10", 
  xlim=c(-1.3, 0.3),h0_name="H0 df=10"
)

# Figur 10 for VaR (df = 10)
plot_cdfs(
  null_data = -VaR_tab2_10_10_n$var_back,
  alternatives = list(-VaR_tab2_10_5_n$var_back, -VaR_tab2_10_3_n$var_back),
  critical_value1 = quantile(-VaR_tab2_10_10_n$var_back, probs = 0.04),
  critical_value2 = quantile(-VaR_tab2_10_10_n$var_back, probs = 0.106),
  pct_5 = NULL,
  pct_10 = NULL,
  alternative_names = c("H1 df=5", "H1 df=3"),
  plot_title = "Figure 10: Norm VaR, df = 10", 
  xlim=c(-8,0),h0_name="H0 df=10"
)
```

Tabel 3 - del 1
```{r}
new_nu <- c(100,10,5,3)
results <- data.frame()

for (df in new_nu) {
    shift_value <- shift(alpha, df = df, dfnull = 100, Norm=FALSE)
    shift_value_n <- shift(alpha, df = df, dfnull = 100, Norm=TRUE)
    ES_2.5pct_3 <- ES(alpha, df = df, shift=shift_value , Norm = FALSE)
    VaR_1pct_3 <- VaR(alpha=beta, shift=shift_value, df = df , Norm = FALSE)
    ES_2.5pct_n_3 <- ES(alpha, df = df, shift=shift_value_n , Norm = TRUE)
    VaR_1pct_n_3 <- VaR(alpha=beta, shift=shift_value_n, df = df , Norm =TRUE)

    results <- rbind(results, data.frame(
      nu = df,
      VaR_1pct_3 = VaR_1pct_3,
      ES_2.5pct_3 = ES_2.5pct_3,
      VaR_1pct_n_3=VaR_1pct_n_3,
      ES_2.5pct_n_3=ES_2.5pct_n_3
    ))
}

print(results)
```

Tabel 3 - del 2. Std student-t
```{r}
#Z1: df=100
Z1_tab3_100_100_VaR<-Z_1(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=100, dfnull = 100), df = 100, dfnull = 100,Norm = FALSE )
Z1_tab3_100_10_VaR <-Z_1(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=10, dfnull = 100), df = 10, dfnull = 100, Norm = FALSE )
Z1_tab3_100_5_VaR <-Z_1(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=5, dfnull = 100), df = 5, dfnull = 100,Norm = FALSE )
Z1_tab3_100_3_VaR<-Z_1(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=3, dfnull = 100), df = 3, dfnull = 100,Norm = FALSE )

#Z2: df=100
Z2_tab3_100_100_VaR <-Z_2(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=100, dfnull = 100) , scale  = 1, df = 100, dfnull = 100, Norm = FALSE )
Z2_tab3_100_10_VaR <-Z_2(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=10, dfnull = 100) , scale  = 1, df = 10, dfnull = 100, Norm = FALSE )
Z2_tab3_100_5_VaR <-Z_2(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=5, dfnull = 100) , scale  = 1, df = 5, dfnull = 100, Norm = FALSE )
Z2_tab3_100_3_VaR <-Z_2(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=3, dfnull = 100) , scale  = 1, df = 3, dfnull = 100, Norm = FALSE )

#Z3: df=100
Z3_tab3_100_100_VaR <-Z_3(quan = 0.05, T=250, M = 100000, df=100, dfnull = 100, shift = shift(0.975, df=100, dfnull = 100), scale = 1, Norm = FALSE)
Z3_tab3_100_10_VaR <-Z_3(quan = 0.05, T=250, M = 100000, df=10, dfnull = 100, shift = shift(0.975, df=10, dfnull = 100), scale = 1, Norm = FALSE)
Z3_tab3_100_5_VaR <-Z_3(quan = 0.05, T=250, M = 100000, df=5, dfnull = 100, shift = shift(0.975, df=5, dfnull = 100), scale = 1, Norm = FALSE)
Z3_tab3_100_3_VaR <-Z_3(quan = 0.05, T=250, M = 100000, df=3, dfnull = 100, shift = shift(0.975, df=3, dfnull = 100), scale = 1, Norm = FALSE)

#VaR: df=100
VaR_tab3_100_100_VaR <- VaR_backtest( sigma2 = 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=100, dfnull = 100) , scale  = 1, df = 100, dfnull = 100, Norm = FALSE )
VaR_tab3_100_10_VaR <-VaR_backtest( sigma2 = 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=10, dfnull = 100) , scale  = 1, df = 10,  dfnull = 100, Norm = FALSE )
VaR_tab3_100_5_VaR <-VaR_backtest( sigma2 = 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=5, dfnull = 100) , scale  = 1, df = 5,  dfnull = 100, Norm = FALSE )
VaR_tab3_100_3_VaR <-VaR_backtest( sigma2 = 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=3, dfnull = 100) , scale  = 1, df = 3,   dfnull = 100, Norm = FALSE )


#Z1: df = 10
Z1_tab3_10_10_VaR <-Z_1(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=10, dfnull = 10), df = 10, dfnull = 10, Norm = FALSE )
Z1_tab3_10_5_VaR <- Z_1(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=5, dfnull = 10), df = 5, dfnull = 10, Norm = FALSE )
Z1_tab3_10_3_VaR<-Z_1(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=3, dfnull = 10), df = 3, dfnull = 10,Norm = FALSE )

#Z2: df = 10
Z2_tab3_10_10_VaR <-Z_2(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=10, dfnull = 10) , scale  = 1, df = 10, dfnull = 10, Norm = FALSE )
Z2_tab3_10_5_VaR <-Z_2(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=5, dfnull = 10) , scale  = 1, df = 5, dfnull = 10, Norm = FALSE )
Z2_tab3_10_3_VaR <-Z_2(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=3, dfnull = 10) , scale  = 1, df = 3, dfnull = 10, Norm = FALSE )

#Z3: df = 10
Z3_tab3_10_10_VaR <-Z_3(quan = 0.05, T=250, M = 100000, df=10, dfnull = 10, shift = shift(0.975, df=10, dfnull = 10), scale = 1, Norm = FALSE)
Z3_tab3_10_5_VaR <-Z_3(quan = 0.05, T=250, M = 100000, df=5, dfnull = 10, shift = shift(0.975, df=5, dfnull = 10), scale = 1, Norm = FALSE)
Z3_tab3_10_3_VaR <-Z_3(quan = 0.05, T=250, M = 100000, df=3, dfnull = 10, shift = shift(0.975, df=3, dfnull = 10), scale = 1, Norm = FALSE)

#VaR: df = 10
VaR_tab3_10_10_VaR <-VaR_backtest( sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=10, dfnull = 10) , scale  = 1, df = 10, dfnull = 10, Norm = FALSE )
VaR_tab3_10_5_VaR <-VaR_backtest( sigma2= 1 , mu = 0, T = 250, M = 100000, shift =shift(0.975, df=5, dfnull = 10) , scale  = 1, df = 5, dfnull = 10, Norm = FALSE )
VaR_tab3_10_3_VaR <-VaR_backtest( sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=3, dfnull = 10) , scale  = 1, df = 3, dfnull = 10, Norm = FALSE )
```

Normalized t
```{r}
#Z1: df=100
Z1_tab3_100_100_VaR_n<-Z_1(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=100, dfnull = 100, Norm=TRUE), df = 100, dfnull = 100, Norm=TRUE )
Z1_tab3_100_10_VaR_n <-Z_1(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=10, dfnull = 100, Norm=TRUE), df = 10, dfnull = 100,  Norm=TRUE )
Z1_tab3_100_5_VaR_n <-Z_1(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=5, dfnull = 100, Norm=TRUE), df = 5, dfnull = 100,  Norm=TRUE )
Z1_tab3_100_3_VaR_n<-Z_1(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=3, dfnull = 100, Norm=TRUE), df = 3, dfnull = 100,  Norm=TRUE )

#Z2: df=100
Z2_tab3_100_100_VaR_n <-Z_2(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=100, dfnull = 100, Norm=TRUE) , scale  = 1, df = 100, dfnull = 100, Norm = TRUE )
Z2_tab3_100_10_VaR_n <-Z_2(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=10, dfnull = 100, Norm=TRUE) , scale  = 1, df = 10, dfnull = 100, Norm = TRUE )
Z2_tab3_100_5_VaR_n <-Z_2(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=5, dfnull = 100, Norm=TRUE) , scale  = 1, df = 5, dfnull = 100, Norm = TRUE )
Z2_tab3_100_3_VaR_n <-Z_2(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=3, dfnull = 100, Norm=TRUE) , scale  = 1, df = 3, dfnull = 100, Norm = TRUE )

#Z3: df=100
Z3_tab3_100_100_VaR_n <-Z_3(quan = 0.05, T=250, M = 100000, df=100, dfnull = 100, shift = shift(0.975, df=100, dfnull = 100, Norm=TRUE), scale = 1, Norm = TRUE)
Z3_tab3_100_10_VaR_n <-Z_3(quan = 0.05, T=250, M = 100000, df=10, dfnull = 100, shift = shift(0.975, df=10, dfnull = 100, Norm=TRUE), scale = 1, Norm = TRUE)
Z3_tab3_100_5_VaR_n <-Z_3(quan = 0.05, T=250, M = 100000, df=5, dfnull = 100, shift = shift(0.975, df=5, dfnull = 100, Norm=TRUE), scale = 1, Norm = TRUE)
Z3_tab3_100_3_VaR_n <-Z_3(quan = 0.05, T=250, M = 100000, df=3, dfnull = 100, shift = shift(0.975, df=3, dfnull = 100, Norm=TRUE), scale = 1, Norm = TRUE)

#VaR: df=100
VaR_tab3_100_100_VaR_n <-VaR_backtest( sigma2 = 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=100, dfnull = 100, Norm=TRUE) , scale  = 1, df = 100, dfnull = 100, Norm = TRUE )
VaR_tab3_100_10_VaR_n <-VaR_backtest( sigma2 = 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=10, dfnull = 100, Norm=TRUE) , scale  = 1, df = 10,  dfnull = 100, Norm = TRUE )
VaR_tab3_100_5_VaR_n <-VaR_backtest( sigma2 = 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=5, dfnull = 100, Norm=TRUE) , scale  = 1, df = 5,  dfnull = 100, Norm = TRUE )
VaR_tab3_100_3_VaR_n <-VaR_backtest( sigma2 = 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=3, dfnull = 100, Norm=TRUE) , scale  = 1, df = 3,   dfnull = 100, Norm = TRUE )

#Z1: df=10
Z1_tab3_10_10_VaR_n <-Z_1(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=10, dfnull = 10, Norm=TRUE), df = 10, dfnull = 10,  Norm=TRUE )
Z1_tab3_10_5_VaR_n <-Z_1(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=5, dfnull = 10, Norm=TRUE), df = 5, dfnull = 10,  Norm=TRUE )
Z1_tab3_10_3_VaR_n<-Z_1(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=3, dfnull = 10, Norm=TRUE), df = 3, dfnull = 10,  Norm=TRUE )

#Z2: df=10
Z2_tab3_10_10_VaR_n <-Z_2(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=10, dfnull = 10, Norm=TRUE) , scale  = 1, df = 10, dfnull = 10, Norm = TRUE )
Z2_tab3_10_5_VaR_n <-Z_2(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=5, dfnull = 10, Norm=TRUE) , scale  = 1, df = 5, dfnull = 10, Norm = TRUE )
Z2_tab3_10_3_VaR_n <-Z_2(quan = 0.05, sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=3, dfnull = 10, Norm=TRUE) , scale  = 1, df = 3, dfnull = 10, Norm = TRUE )

#Z3: df=10
Z3_tab3_10_10_VaR_n <-Z_3(quan = 0.05, T=250, M = 100000, df=10, dfnull = 10, shift = shift(0.975, df=10, dfnull = 10, Norm=TRUE), scale = 1, Norm = TRUE)
Z3_tab3_10_5_VaR_n <-Z_3(quan = 0.05, T=250, M = 100000, df=5, dfnull = 10, shift = shift(0.975, df=5, dfnull = 10, Norm=TRUE), scale = 1, Norm = TRUE)
Z3_tab3_10_3_VaR_n <-Z_3(quan = 0.05, T=250, M = 100000, df=3, dfnull = 10, shift = shift(0.975, df=3, dfnull = 10, Norm=TRUE), scale = 1, Norm = TRUE)

#VaR: df=10
VaR_tab3_10_10_VaR_n <-VaR_backtest(sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=10, dfnull = 10, Norm=TRUE), scale  = 1, df = 10, dfnull = 10 , Norm = TRUE)
VaR_tab3_10_5_VaR_n <-VaR_backtest( sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=5, dfnull = 10, Norm=TRUE) , scale  = 1, df = 5, dfnull = 10 , Norm = TRUE)
VaR_tab3_10_3_VaR_n <-VaR_backtest( sigma2= 1 , mu = 0, T = 250, M = 100000, shift = shift(0.975, df=3, dfnull = 10, Norm=TRUE) , scale  = 1, df = 3, dfnull = 10 , Norm = TRUE)
```

Printer std. student-t tabel
```{r}
names100 <- c("Powerdf10", "Powerdf3")
names10 <- c("Powerdf5", "Powerdf3")
significance_levels100_tab3 <- c(0.043, 0.109)
significance_levels10_tab3 <- c(0.041,0.107)
significance_levels100_tab3_std <- c(0.041, 0.104)
significance_levels10_tab3_std <- c(0.04, 0.106)

# H_0: df = 100
tabel3_100_Z1 <- compute_table(
  Z1_tab3_100_100_VaR$Z_1_list,
  list(Z1_tab3_100_10_VaR$Z_1_list, Z1_tab3_100_3_VaR$Z_1_list),
  significance_levels100_tab3, 100, "Z1", names100
)

tabel3_100_Z2 <- compute_table(
  Z2_tab3_100_100_VaR$Z_2_list,
  list(Z2_tab3_100_10_VaR$Z_2_list, Z2_tab3_100_3_VaR$Z_2_list),
  significance_levels100_tab3, 100, "Z2", names100
)

tabel3_100_Z3 <- compute_table(
  Z3_tab3_100_100_VaR$Z_3_list,
  list(Z3_tab3_100_10_VaR$Z_3_list, Z3_tab3_100_3_VaR$Z_3_list),
  significance_levels100_tab3, 100, "Z3", names100
)

tabel3_100_VaR <- compute_table(
  (-VaR_tab3_100_100_VaR$var_back),
  list((-VaR_tab3_100_10_VaR$var_back), (-VaR_tab3_100_3_VaR$var_back)),
  significance_levels100_tab3_std, 100, "VaR", names100
)

# H0: df = 10
tabel3_10_Z1 <- compute_table(
  Z1_tab3_10_10_VaR$Z_1_list,
  list(Z1_tab3_10_5_VaR$Z_1_list, Z1_tab3_10_3_VaR$Z_1_list),
  significance_levels10_tab3, 10, "Z1", names10
)

tabel3_10_Z2 <- compute_table(
  Z2_tab3_10_10_VaR$Z_2_list,
  list(Z2_tab3_10_5_VaR$Z_2_list, Z2_tab3_10_3_VaR$Z_2_list),
  significance_levels10_tab3, 10, "Z2", names10
)

tabel3_10_Z3 <- compute_table(
  Z3_tab3_10_10_VaR$Z_3_list,
  list(Z3_tab3_10_5_VaR$Z_3_list, Z3_tab3_10_3_VaR$Z_3_list),
  significance_levels10_tab3, 10, "Z3", names10
)

tabel3_10_VaR <- compute_table(
  (-VaR_tab3_10_10_VaR$var_back),
  list((-VaR_tab3_10_5_VaR$var_back), (-VaR_tab3_10_3_VaR$var_back)),
  significance_levels10_tab3_std, 10, "VaR", names10
)

tabel3_std <- bind_rows(
  tabel3_100_Z1, tabel3_100_Z2, tabel3_100_Z3, tabel3_100_VaR,
  tabel3_10_Z1, tabel3_10_Z2, tabel3_10_Z3, tabel3_10_VaR
) %>%
  mutate(
    Significance_Level = Significance_Level * 100,
    Source = factor(Source, levels = c("Z1","Z2", "Z3", "VaR"))
  ) %>%
  arrange(H0_df, Significance_Level, Source)

tabel3_std <- tabel3_std[, c("Significance_Level", "Powerdf3", "Powerdf5", "Powerdf10", "H0_df", "Source")]

print(tabel3_std)
```

Printer normalized student-t tabel
```{r}
significance_levels100_tab3_norm <- c(0.041, 0.111)
significance_levels10_tab3_norm <- c(0.042, 0.114)
significance_levels100_tab3_std <- c(0.041, 0.104)
significance_levels10_tab3_std <- c(0.04, 0.106)
names100 <- c("Powerdf10", "Powerdf3")
names10 <- c("Powerdf5", "Powerdf3")

# H_0: df = 100
tabel3_100_Z1_norm <- compute_table(
  Z1_tab3_100_100_VaR_n$Z_1_list,
  list(Z1_tab3_100_10_VaR_n$Z_1_list, Z1_tab3_100_3_VaR_n$Z_1_list),
  significance_levels100_tab3_norm, 100, "Z1", names100
)

tabel3_100_Z2_norm <- compute_table(
  Z2_tab3_100_100_VaR_n$Z_2_list,
  list(Z2_tab3_100_10_VaR_n$Z_2_list, Z2_tab3_100_3_VaR_n$Z_2_list),
  significance_levels100_tab3_norm, 100, "Z2", names100
)

tabel3_100_Z3_norm <- compute_table(
  Z3_tab3_100_100_VaR_n$Z_3_list,
  list(Z3_tab3_100_10_VaR_n$Z_3_list, Z3_tab3_100_3_VaR_n$Z_3_list),
  significance_levels100_tab3_norm, 100, "Z3", names100
)

tabel3_100_VaR_norm <- compute_table(
  (-VaR_tab3_100_100_VaR_n$var_back),
  list(-VaR_tab3_100_10_VaR_n$var_back, -VaR_tab3_100_3_VaR_n$var_back),
  significance_levels100_tab3_std, 100, "VaR", names100
)

# H0: df = 10
tabel3_10_Z1_norm <- compute_table(
  Z1_tab3_10_10_VaR_n$Z_1_list,
  list(Z1_tab3_10_5_VaR_n$Z_1_list, Z1_tab3_10_3_VaR_n$Z_1_list),
  significance_levels10_tab3_norm, 10, "Z1", names10
)

tabel3_10_Z2_norm <- compute_table(
  Z2_tab3_10_10_VaR_n$Z_2_list,
  list(Z2_tab3_10_5_VaR_n$Z_2_list, Z2_tab3_10_3_VaR_n$Z_2_list),
  significance_levels10_tab3_norm, 10, "Z2", names10
)

tabel3_10_Z3_norm <- compute_table(
  Z3_tab3_10_10_VaR_n$Z_3_list,
  list(Z3_tab3_10_5_VaR_n$Z_3_list, Z3_tab3_10_3_VaR_n$Z_3_list),
  significance_levels10_tab3_norm, 10, "Z3", names10
)

tabel3_10_VaR_norm <- compute_table(
  (-VaR_tab3_10_10_VaR_n$var_back),
  list((-VaR_tab3_10_5_VaR_n$var_back), (-VaR_tab3_10_3_VaR_n$var_back)),
  significance_levels10_tab3_std, 10, "VaR", names10
)

tabel3_norm <- bind_rows(
  tabel3_100_Z1_norm, tabel3_100_Z2_norm, tabel3_100_Z3_norm, tabel3_100_VaR_norm,
  tabel3_10_Z1_norm, tabel3_10_Z2_norm, tabel3_10_Z3_norm, tabel3_10_VaR_norm
) %>%
  mutate(
    Significance_Level = Significance_Level * 100,
    Source = factor(Source, levels = c("Z1","Z2", "Z3", "VaR"))
  ) %>%
  arrange(H0_df, Significance_Level, Source)

tabel3_norm <- tabel3_norm[, c("Significance_Level", "Powerdf3", "Powerdf5", "Powerdf10", "H0_df", "Source")]

print(tabel3_norm)
```

Figur 11
```{r}
generate_t_distributions_VaR <- function(df_values, x_range) {
  data <- data.frame()
  
  for (df in df_values) {
    shift_value <- shift(0.025, df = df, dfnull = 100)
    x <- seq(x_range[1], x_range[2], length.out = 500)
    shifted_x <- (x + shift_value)
    y <- dt(shifted_x, df = df) 
    temp <- data.frame(x = x, y = y, df = factor(df))
    data <- rbind(data, temp)
  }
  
  return(data)
}

df_values <- c(100, 10, 5, 3)
x_range <- c(-6, 6)
fixed_var <- -VaR(alpha = 0.975, sigma2 = 1, mu = 0, shift = 0, scale = 1, df = 100 , Norm = FALSE)

data_t <- generate_t_distributions_VaR(df_values, x_range)

ggplot(data_t, aes(x = x, y = y, color = df)) +
  geom_line(size = 0.8) + 
  geom_vline(xintercept = fixed_var, color = "black", size = 1, linetype = "solid") + 
  labs(
    title = "t-distributions with equal VaR2.5%",
    x = NULL, y = NULL,
    color = "df"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  scale_color_manual(
    values = c("red", "purple", "blue", "green4")
  )
```

Figur 12
```{r}
# Figur 12 for Z1 (df = 100)
plot_cdfs(
  null_data = Z1_tab3_100_100_VaR$Z_1_list,
  alternatives = list(Z1_tab3_100_10_VaR$Z_1_list, Z1_tab3_100_5_VaR$Z_1_list, Z1_tab3_100_3_VaR$Z_1_list),
  critical_value1 = quantile(Z1_tab3_100_100_VaR$Z_1_list, probs = 0.043),
  critical_value2 = quantile(Z1_tab3_100_100_VaR$Z_1_list, probs = 0.109),
  pct_5 = quantile(Z1_tab3_100_100_VaR$Z_1_list, probs = 0.05),
  pct_10 = quantile(Z1_tab3_100_100_VaR$Z_1_list, probs = 0.1),
  alternative_names = c("H1 df=10", "H1 df=5","H1 df=3"),
  plot_title = "Figure 12: Norm Z1, df = 100",
  xlim=c(-2.5, 0.25),h0_name="H0 df=100"
)

# Figur 12 for Z2 (df = 100)
plot_cdfs(
  null_data = Z2_tab3_100_100_VaR$Z_2_list,
  alternatives = list(Z2_tab3_100_10_VaR$Z_2_list,Z2_tab3_100_5_VaR$Z_2_list, Z2_tab3_100_3_VaR$Z_2_list),
  critical_value1 = quantile(Z2_tab3_100_100_VaR$Z_2_list, probs = 0.043),
  critical_value2 = quantile(Z2_tab3_100_100_VaR$Z_2_list, probs = 0.109),
  pct_5 = quantile(Z2_tab3_100_100_VaR$Z_2_list, probs = 0.05),
  pct_10 = quantile(Z2_tab3_100_100_VaR$Z_2_list, probs = 0.1),
  alternative_names = c("H1 df=10", "H1 df=5","H1 df=3"),
  plot_title = "Figure 12: Norm Z2, df = 100",
  xlim=c(-2.5, 0.75),h0_name="H0 df=100"
)

# Figur 12 for Z3 (df = 100)
plot_cdfs(
  null_data = Z3_tab3_100_100_VaR$Z_3_list,
  alternatives = list(Z3_tab3_100_10_VaR$Z_3_list,Z3_tab3_100_5_VaR$Z_3_list, Z3_tab3_100_3_VaR$Z_3_list),
  critical_value1 = quantile(Z3_tab3_100_100_VaR$Z_3_list, probs = 0.043),
  critical_value2 = quantile(Z3_tab3_100_100_VaR$Z_3_list, probs = 0.109),
  pct_5 = quantile(Z3_tab3_100_100_VaR$Z_3_list, probs = 0.05),
  pct_10 = quantile(Z3_tab3_100_100_VaR$Z_3_list, probs = 0.1),
  alternative_names = c("H1 df=10", "H1 df=5","H1 df=3"),
  plot_title = "Figure 12: Norm Z3, df = 100",
  xlim=c(-2.5, 0.25),h0_name="H0 df=100"
)

# Figur 12 for VaR (df = 100)
plot_cdfs(
  null_data = -VaR_tab3_100_100_VaR$var_back,
  alternatives = list(-VaR_tab3_100_10_VaR$var_back,-VaR_tab3_100_5_VaR$var_back, -VaR_tab3_100_3_VaR$var_back),
  critical_value1 = quantile(-VaR_tab3_100_100_VaR$var_back, probs = 0.043),
  critical_value2 = quantile(-VaR_tab3_100_100_VaR$var_back, probs = 0.109),
  pct_5 = NULL,
  pct_10 = NULL,
  alternative_names = c("H1 df=10", "H1 df=5","H1 df=3"),
  plot_title = "Figure 12: Norm VaR, df = 100",h0_name="H0 df=100"
)
```

Figur 13
```{r}
# Figur 13 for Z1 (df = 10)
plot_cdfs(
  null_data = Z1_tab3_10_10_VaR$Z_1_list,
  alternatives = list(Z1_tab3_10_5_VaR$Z_1_list, Z1_tab3_10_3_VaR$Z_1_list),
  critical_value1 = quantile(Z1_tab3_10_10_VaR$Z_1_list, probs = 0.041),
  critical_value2 = quantile(Z1_tab3_10_10_VaR$Z_1_list, probs = 0.107),
  pct_5 = quantile(Z1_tab3_10_10_VaR$Z_1_list, probs = 0.05),
  pct_10 = quantile(Z1_tab3_10_10_VaR$Z_1_list, probs = 0.1),
  alternative_names = c("H1 df=5", "H1 df=3"),
  plot_title = "Figure 13: Norm Z1, df = 10",
  xlim=c(-2.1, 0.25),h0_name="H0 df=10"
)

# Figur 13 for Z2 (df = 10)
plot_cdfs(
  null_data = Z2_tab3_10_10_VaR$Z_2_list,
  alternatives = list(Z2_tab3_10_5_VaR$Z_2_list, Z2_tab3_10_3_VaR$Z_2_list),
  critical_value1 = quantile(Z2_tab3_10_10_VaR$Z_2_list, probs = 0.041),
  critical_value2 = quantile(Z2_tab3_10_10_VaR$Z_2_list, probs = 0.107),
  pct_5 = quantile(Z2_tab3_10_10_VaR$Z_2_list, probs = 0.05),
  pct_10 = quantile(Z2_tab3_10_10_VaR$Z_2_list, probs = 0.1),
  alternative_names = c("H1 df=5", "H1 df=3"),
  plot_title = "Figure 13: Norm Z2, df = 10", 
  xlim=c(-2.5,0.75),h0_name="H0 df=10"
)

# Figur 13 for Z3 (df = 10)
plot_cdfs(
  null_data = Z3_tab3_10_10_VaR$Z_3_list,
  alternatives = list(Z3_tab3_10_5_VaR$Z_3_list, Z3_tab3_10_3_VaR$Z_3_list),
  critical_value1 = quantile(Z3_tab3_10_10_VaR$Z_3_list, probs = 0.041),
  critical_value2 = quantile(Z3_tab3_10_10_VaR$Z_3_list, probs = 0.107),
  pct_5 = quantile(Z3_tab3_10_10_VaR$Z_3_list, probs = 0.05),
  pct_10 = quantile(Z3_tab3_10_10_VaR$Z_3_list, probs = 0.1),
  alternative_names = c("H1 df=5", "H1 df=3"),
  plot_title = "Figure 13: Norm Z3, df = 10",
  xlim=c(-2.1, 0.25),h0_name="H0 df=10"
)

# Figur 13 for VaR (df = 10)
plot_cdfs(
  null_data = -VaR_tab3_10_10_VaR$var_back,
  alternatives = list(-VaR_tab3_10_5_VaR$var_back, -VaR_tab3_10_3_VaR$var_back),
  critical_value1 = quantile(-VaR_tab3_10_10_VaR$var_back, probs = 0.041),
  critical_value2 = quantile(-VaR_tab3_10_10_VaR$var_back, probs = 0.107),
  pct_5 = NULL,
  pct_10 = NULL,
  alternative_names = c("H1 df=5", "H1 df=3"),
  plot_title = "Figure 13: Norm VaR, df = 10",
  xlim=c(-10,0),h0_name="H0 df=10"
)
```


