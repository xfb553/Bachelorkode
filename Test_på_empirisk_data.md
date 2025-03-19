Pakker: 
```{r, warning=FALSE}
library(ggplot2)
library(gridExtra)
library(rugarch)
library(esback)
library(readr)
library(dplyr)
```

VaR og ES
```{r, warning=FALSE}
VaR <- function(alpha = 0.975, sigma2 = 1, mu = 0, df = NULL , Norm = FALSE) {
  if (!Norm) {
      VaR <- qt(1 - alpha, df = df)
  } else {
      VaR <- sqrt((df - 2) / df) * qt(1 - alpha, df = df)
  }
  return(-VaR)
}

ES <- function(alpha = 0.975, sigma2 = 1, mu = 0, df = NULL , Norm = FALSE) {
  if (!Norm) {
      x <- qt(1 - alpha, df = df)
      ES <- -(dt(x, df) / (1 - alpha) * (df + x^2) / (df - 1))
  } else {
      x <- qt(1 - alpha, df = df)
      ES <- -sqrt((df - 2) / df) * (dt(x, df) / (1 - alpha) * (df + x^2) / (df - 1))
  }
  return(-ES)
}
```

Test 1
```{r,warning=FALSE}
Z_1 <- function(R_t, var, es) {
  T<-250
    q <- 0
    P <- var + R_t
    N_t <- sum(P < 0)

    for (j in 1:T) {
      if (R_t[j] + var[j] < 0) {
        q <- q + (R_t[j] / es[j])
      }
    }
    if (N_t == 0) {
      Z_1 <- NA
      } else {
      Z_1 <- q / N_t + 1
      }


  return(list(Z_1=Z_1))
}
```

Test 2
```{r,warning=FALSE}
Z_2 <- function(R_t,var,es) {
  T<-250

    q <- 0
    Talpha <- T * 0.025

    for (j in 1:T) {
      if (R_t[j] + var[j] < 0) {
        q <- q + R_t[j] / es[j]
      }
    }
    Z_2 <- q / Talpha + 1
  return(list(Z_2 = Z_2))
}
```

CC-backtest
```{r,warning=FALSE}
cc_backtest <- function(r, q, e, s=NULL, alpha) {
  
  n <- length(r)

  V <- cbind(1-alpha - (r > q),
             q - e - (r > q) * (q - r) / (1-alpha))
    H3 <- array(NA, dim = c(n, 4, 2))  # q = 4
    H3[,, 1] <- cbind(1, abs(q), 0, 0)
    H3[,, 2] <- cbind(0, 0, 1, 1 / s)

  hV3 <- apply(H3,2, function(x) rowSums(x*V))

    omega3 <- crossprod(hV3) / n
    
    t4 <- sqrt(n) * diag(omega3)^(-1/2) * colMeans(hV3)
    p4 <- min(4 * sum(1 / (1:4)) * min(sort(1 - stats::pnorm(t4)) / (1:4)), 1)
    min_index<-which.min(sort(1 - stats::pnorm(t4)) / (1:4))
    min_t_value <- t4[min_index] 

  # Return results
  ret <- list(t2_list=t4,
    pvalue_onesided_general = p4,
    min_t2=min_t_value
  )
  ret
}
```

ESR-backtest:
```{r,warning=FALSE}
esr_backtest <- function(r, e,  alpha,
                         cov_config=list(sparsity='nid', sigma_est='scl_sp', misspec=TRUE)) {

    data <- data.frame(r = r, e = e)
    model <- I(r - e) ~ e | 1
    h0 <- c(NA, NA, 0)
    one_sided <- TRUE

  fit0 <- esreg::esreg(model, data = data, alpha = 1-alpha, g1 = 2, g2 = 1)
  cov0 <- esreg::vcovA(fit0,
                       sparsity = cov_config$sparsity,
                       sigma_est = cov_config$sigma_est,
                       misspec = cov_config$misspec)
  s0 <- fit0$coefficients - h0
  mask <- !is.na(h0)

    t0 <- as.numeric(s0[mask] / sqrt(cov0[mask, mask]))
    pv0_1s <- (1-stats::pnorm(t0))

  ret <- list(
    pvalue_onesided_asymptotic = pv0_1s,
    t0=t0
  )

  ret
}

```
Indlæsning af data til estination af parametre under AR(1)-GARCH(1,1):
```{r, warning=FALSE}

SP500_Data <- read_delim("C:/Users/Maja/OneDrive - University of Copenhagen/Desktop/Bachelor/S&P 500 Historical Data.csv",
    delim = ";", escape_double = FALSE, col_types = cols(Price = col_number(), 
        Open = col_number(), High = col_number(), 
        Low = col_number(), Vol. = col_character()), 
    trim_ws = TRUE)

SP500_Data$Date <- gsub("-", "/", SP500_Data$Date)
SP500_Data$Date <- as.Date(SP500_Data$Date, format = "%m/%d/%Y")
SP500_Data <- SP500_Data %>%
  arrange(Date) %>%
  mutate(Days = as.numeric(Date - min(Date)),
         Change = as.numeric(gsub("%", "", `Change %`)) / 100)

test_days <- 250

SP500_Train <- SP500_Data[1:(nrow(SP500_Data) - test_days), ]
SP500_Train_returns <- SP500_Train$Change*100

SP500_Test <- SP500_Data[(nrow(SP500_Data) - test_days + 1):nrow(SP500_Data), ]
SP500_Test_returns <- SP500_Test$Change*100

ARGARCH <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(1, 0)),
  distribution.model = "std"
)

SP500_ARGARCH_fit <- ugarchfit(spec = ARGARCH, data = SP500_Train_returns)
fitted_values_train <- coef(SP500_ARGARCH_fit)
```

Tester for autokorrelation i AR(1)-GARCH(1,1)
```
fit<-ugarchfit(spec = ARGARCH, data = SP500_Test$Change*100)
std_resid <- residuals(fit, standardize = TRUE)
Box.test(std_resid, lag = 21, type = "Ljung-Box")
```

Test for autokorrelation i GARCH(1,1)
```
ARGARCH <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0)),  # AR(1) del
  distribution.model = "std")

fit_train<-ugarchfit(spec = ARGARCH, data = SP500_Train$Change*100)
std_resid_train <- residuals(fit_train, standardize = TRUE)
Box.test(std_resid_train, lag = 21, type = "Ljung-Box")
```

Ligger testdata-returns inden for modellens VaR og ES?
```{r, warning=FALSE}
T<-250

omega_H0 <- fitted_values_train["omega"]
alpha_H0 <- fitted_values_train["alpha1"]
beta_H0 <- fitted_values_train["beta1"]
phi_H0 <- fitted_values_train["ar1"]
df_H0 <- fitted_values_train["shape"]
mu_H0 <- fitted_values_train["mu"]
phi_H0 <- fitted_values_train["ar1"]

sigma2_H0 <- numeric(length(SP500_Test_returns))
sigma2_H0[1] <- omega_H0 / (1 - alpha_H0- beta_H0)

R_t <- SP500_Test_returns  

for (i in 2:length(R_t)) {
  sigma2_H0[i] <- omega_H0 + alpha_H0 * (R_t[i - 1])^2 + beta_H0 * sigma2_H0[i - 1]
}

z <-  rt(T,df_H0)*sqrt((df_H0 - 2) / df_H0)

VaR_values_H0 <- VaR(alpha = 0.975, df = df_H0) * sqrt(sigma2_H0)
ES_values_H0 <- ES(alpha = 0.975, df = df_H0) * sqrt(sigma2_H0)


plot_data <- data.frame(
  Time = 1:length(R_t),
  Returns = R_t,
  VaR = -VaR_values_H0,
  ES = -ES_values_H0
)

plot_data$Exceedance <- ifelse(plot_data$Returns < plot_data$VaR, "Exceedance", "No Exceedance")

ggplot(plot_data, aes(x = Time)) +
  geom_line(aes(y = Returns, color = "SP500 Returns")) +
  geom_line(aes(y = VaR, color = "VaR (97.5%) - H0")) +
  geom_line(aes(y = ES, color = "ES (97.5%) - H0")) +
  geom_point(data = plot_data[plot_data$Exceedance == "Exceedance", ], aes(y = Returns), color = "black", size = 1) +
  labs(
    title = "SP500 Returns vs. VaR/ES calculated from H0",
    x = "Time in days",
    y = "Returns in %",
    color = "Legend"
  ) +
  scale_color_manual(values = c("SP500 Returns" = "blue", "VaR (97.5%) - H0" = "red", "ES (97.5%) - H0" = "green")) +
  theme_minimal()

```

Varrierer alpha, for at se hvornår VaR og ES over- eller undervurderes
Lav alpha
```{r,warning=FALSE}
omega_H0 <- fitted_values_train["omega"]
alpha_H1 <- 0.03
beta_H1 <- fitted_values_train["beta1"]+ fitted_values_train["alpha1"]-alpha_H1
phi_H0 <- fitted_values_train["ar1"]
df_H0 <- fitted_values_train["shape"]
phi_H0 <- fitted_values_train["ar1"]

sigma2_H1 <- numeric(length(SP500_Test$Change))
sigma2_H1[1] <- omega_H0 / (1 - alpha_H1- beta_H1)

R_t <- SP500_Test_returns  

for (i in 2:length(R_t)) {
  sigma2_H1[i] <- omega_H0 + alpha_H1 * (R_t[i - 1])^2 + beta_H1 * sigma2_H1[i - 1]
}

z <- rt(T,df_H0)*sqrt((df_H0 - 2) / df_H0)

VaR_values_H1 <- VaR(alpha = 0.975, df = df_H0) * sqrt(sigma2_H1)
ES_values_H1 <- ES(alpha = 0.975, df = df_H0) * sqrt(sigma2_H1)

plot_data_low <- data.frame(
  Time = 1:length(R_t),
  Returns = R_t,
  VaR = -VaR_values_H1, 
  ES = -ES_values_H1  
)

plot_data_low$Exceedance <- ifelse(plot_data_low$Returns < plot_data_low$VaR, "Exceedance", "No Exceedance")

ggplot(plot_data_low, aes(x = Time)) +
  geom_line(aes(y = Returns, color = "Returns")) +
  geom_line(aes(y = VaR, color = "VaR (97.5%)")) +
  geom_line(aes(y = ES, color = "ES (97.5%)")) +
  geom_point(data = plot_data_low[plot_data_low$Exceedance == "Exceedance", ], aes(y = Returns), color = "black", size = 1) +
  labs(
    x = "Time in days",
    y = "Returns in %",
    color = "Legend"
  ) +
  scale_color_manual(values = c("Returns" = "blue", "VaR (97.5%)" = "red", "ES (97.5%)" = "green")) +
  theme_minimal()
```

Høj alpha
```{r,warning=FALSE}
omega_H0 <- fitted_values_train["omega"]
alpha_H1 <- 0.25
beta_H1 <- fitted_values_train["beta1"]+fitted_values_train["alpha1"]-alpha_H1
phi_H0 <- fitted_values_train["ar1"]
df_H0 <- fitted_values_train["shape"]
phi_H0 <- fitted_values_train["ar1"]

sigma2_H1 <- numeric(length(SP500_Test$Change))
sigma2_H1[1] <- omega_H0 / (1 - alpha_H1- beta_H1)
R_t <- SP500_Test_returns

for (i in 2:length(R_t)) {
  sigma2_H1[i] <- omega_H0 + alpha_H1 * (R_t[i - 1])^2 + beta_H1 * sigma2_H1[i - 1]
}

z <- rt(T,df_H0)*sqrt((df_H0 - 2) / df_H0)

VaR_values_H1 <- VaR(alpha = 0.975, df = df_H0) * sqrt(sigma2_H1)
ES_values_H1 <- ES(alpha = 0.975, df = df_H0) * sqrt(sigma2_H1)

plot_data_high <- data.frame(
  Time = 1:length(R_t),
  Returns = R_t ,
  VaR = -VaR_values_H1, 
  ES = -ES_values_H1
)

plot_data_high$Exceedance <- ifelse(plot_data_high$Returns < plot_data_high$VaR, "Exceedance", "No Exceedance")

ggplot(plot_data_high, aes(x = Time)) +
  geom_line(aes(y = Returns, color = "Returns")) +
  geom_line(aes(y = VaR, color = "VaR (97.5%)")) +
  geom_line(aes(y = ES, color = "ES (97.5%)")) +
  geom_point(data = plot_data_high[plot_data_high$Exceedance == "Exceedance", ], aes(y = Returns), color = "black", size = 1) +
  labs(
    x = "Time in days",
    y = "Returns in %",
    color = "Legend"
  ) +
  scale_color_manual(values = c("Returns" = "blue", "VaR (97.5%)" = "red", "ES (97.5%)" = "green")) +
  theme_minimal()
```

Samlet alpha plot
```{r, warning=FALSE}
plot_data$Scenario <- "Alpha=0.166"
plot_data_high$Scenario <- "Alpha=0.25"
plot_data_low$Scenario <- "Alpha=0.03"

combined_plot_data <- rbind(plot_data, plot_data_low, plot_data_high)

ggplot(combined_plot_data, aes(x = Time)) +
  geom_line(aes(y = Returns, color = "Returns", linetype = "Returns"), size = 0.52) +
  geom_line(aes(y = ES, color = Scenario, linetype = Scenario), size = 0.52) +
  labs(
    title = "SP500 Returns vs. Expected Shortfall (ES) for Different Alpha Levels",
    x = "Time in days",
    y = "Returns and ES in %",
    color = "ARCH parameter",
    linetype = "ARCH parameter"
  ) +
  scale_color_manual(values = c("Returns" = "black",        
                                "Alpha=0.25" = "blue",   
                                "Alpha=0.166" = "green",  
                                "Alpha=0.03" = "red")) +
  scale_linetype_manual(values = c("Returns" = "solid", 
                                   "Alpha=0.25" = "solid",  
                                   "Alpha=0.166" = "solid",  
                                   "Alpha=0.03" = "solid")) +  
  theme_minimal(base_size = 10) +  
  theme(
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid.major = element_line(color = "gray80", linetype = "dotted"),  
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)  
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 2)), 
    linetype = guide_legend(override.aes = list(size = 1.2))
  )
```


Lav omega
```{r,warning=FALSE}
omega_H1 <- 0.01
alpha_H0 <- fitted_values_train["alpha1"]
beta_H0 <- fitted_values_train["beta1"]
phi_H0 <- fitted_values_train["ar1"]
mu_H0 <- fitted_values_train["mu"]
df_H0 <- fitted_values_train["shape"]
phi_H0 <- fitted_values_train["ar1"]

sigma2_H1 <- numeric(length(SP500_Test$Change))
sigma2_H1[1] <- omega_H1 / (1 - alpha_H0- beta_H0)

R_t <- SP500_Test_returns  

for (i in 2:length(R_t)) {
  sigma2_H1[i] <- omega_H1 + alpha_H0 * (R_t[i - 1])^2 + beta_H0 * sigma2_H1[i - 1]
}

z <- rt(T,df_H0)*sqrt((df_H0 - 2) / df_H0)

VaR_values_H1 <- VaR(alpha = 0.975, df = df_H0) * sqrt(sigma2_H1)
ES_values_H1 <- ES(alpha = 0.975, df = df_H0) * sqrt(sigma2_H1)

plot_data_omega_low <- data.frame(
  Time = 1:length(R_t),
  Returns = R_t ,
  VaR = -VaR_values_H1 ,
  ES = -ES_values_H1
)

plot_data_omega_low$Exceedance <- ifelse(plot_data_omega_low$Returns < plot_data_omega_low$VaR, "Exceedance", "No Exceedance")

ggplot(plot_data_omega_low, aes(x = Time)) +
  geom_line(aes(y = Returns, color = "Returns")) +
  geom_line(aes(y = VaR, color = "VaR (97.5%)")) +
  geom_line(aes(y = ES, color = "ES (97.5%)")) +
  geom_point(data = plot_data_omega_low[plot_data_omega_low$Exceedance == "Exceedance", ], aes(y = Returns), color = "black", size = 1) +
  labs(
    x = "Time in days",
    y = "Returns in %",
    color = "Legend"
  ) +
  scale_color_manual(values = c("Returns" = "blue", "VaR (97.5%)" = "red", "ES (97.5%)" = "green")) +
  theme_minimal()

```

Høj omega
```{r,warning=FALSE}
omega_H1 <- 0.08  
alpha_H0 <- fitted_values_train["alpha1"]
beta_H0 <- fitted_values_train["beta1"]
phi_H0 <- fitted_values_train["ar1"]
mu_H0 <- fitted_values_train["mu"]
df_H0 <- fitted_values_train["shape"]
phi_H0 <- fitted_values_train["ar1"]

sigma2_H1 <- numeric(length(SP500_Test$Change))
sigma2_H1[1] <- omega_H1 / (1 - alpha_H0- beta_H0)

R_t <- SP500_Test_returns  

for (i in 2:length(R_t)) {
  sigma2_H1[i] <- omega_H1 + alpha_H0 * (R_t[i - 1])^2 + beta_H0 * sigma2_H1[i - 1]
}

z <- rt(T,df_H0)*sqrt((df_H0 - 2) / df_H0)

VaR_values_H1 <- VaR(alpha = 0.975, df = df_H0) * sqrt(sigma2_H1)
ES_values_H1 <- ES(alpha = 0.975, df = df_H0) * sqrt(sigma2_H1)

plot_data_omega_high <- data.frame(
  Time = 1:length(R_t),
  Returns = R_t,
  VaR = -VaR_values_H1,
  ES = -ES_values_H1
)

plot_data_omega_high$Exceedance <- ifelse(plot_data_omega_high$Returns < plot_data_omega_high$VaR, "Exceedance", "No Exceedance")

ggplot(plot_data_omega_high, aes(x = Time)) +
  geom_line(aes(y = Returns, color = "Returns")) +
  geom_line(aes(y = VaR, color = "VaR (97.5%)")) +
  geom_line(aes(y = ES, color = "ES (97.5%)")) +
  geom_point(data = plot_data_omega_high[plot_data_omega_high$Exceedance == "Exceedance", ], aes(y = Returns), color = "black", size = 1) +
  labs(
    x = "Time in days",
    y = "Returns in %",
    color = "Legend"
  ) +
  scale_color_manual(values = c("Returns" = "blue", "VaR (97.5%)" = "red", "ES (97.5%)" = "green")) +
  theme_minimal()

```

Samlet omega plot
```{r, warning=FALSE}
plot_data$Scenario <- "Omega=0.042"
plot_data_omega_low$Scenario <- "Omega=0.004"
plot_data_omega_high$Scenario <- "Omega=0.1"

combined_plot_data <- rbind(plot_data, plot_data_omega_low, plot_data_omega_high)

ggplot(combined_plot_data, aes(x = Time)) +
  geom_line(aes(y = Returns, color = "Returns", linetype = "Returns"), size = 0.52) +
  geom_line(aes(y = ES, color = Scenario, linetype = Scenario), size = 0.52) +
  labs(
    title = "SP500 Returns vs. Expected Shortfall (ES) for Different Omega Levels",
    x = "Time in days",
    y = "Returns and ES in %",
    color = "Omega parameter",
    linetype = "Omega parameter"
  ) +
  scale_color_manual(values = c("Returns" = "black",        
                                "Omega=0.1" = "blue",   
                                "Omega=0.042" = "green",  
                                "Omega=0.004" = "red")) +
  scale_linetype_manual(values = c("Returns" = "solid", 
                                   "Omega=0.1" = "solid",  
                                   "Omega=0.042" = "solid",  
                                   "Omega=0.004" = "solid")) +  
  theme_minimal(base_size = 10) +  
  theme(
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid.major = element_line(color = "gray80", linetype = "dotted"),  
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)  
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 2)), 
    linetype = guide_legend(override.aes = list(size = 1.2))
  )
```


Laver 250x10.000 simulationer under AR-GARCH sande parametre "H0"
```{r, warning=FALSE}
set.seed(110694)

fitted_values_train <- coef(SP500_ARGARCH_fit)
omega <- fitted_values_train["omega"]
alpha <- fitted_values_train["alpha1"]
beta <- fitted_values_train["beta1"]
df <- fitted_values_train["shape"]
mu <- fitted_values_train["mu"]
phi <- fitted_values_train["ar1"]

M <- 10000
T <- 250  

R_t_matrix <- matrix(NA, nrow = T, ncol = M)
sigma2_matrix <- matrix(NA, nrow = T, ncol = M)
sigma_matrix <- matrix(NA, nrow = T, ncol = M)
var_matrix <- matrix(NA, nrow = T, ncol = M)
es_matrix <- matrix(NA, nrow = T, ncol = M)

for (sim in 1:M) {
  z <- sqrt((df - 2) / df) * rt(T, df)
  
  sigma2_sim <- rep(NA, T)
  sigma_sim <- rep(NA, T)
  x <- rep(NA, T)

  sigma2_sim[1] <- omega /(1 - alpha - beta)
  x[1] <- sqrt(sigma2_sim[1]) * z[1]
  sigma_sim[1]<-sqrt(sigma2_sim[1])

  for (i in 2:T) {
    sigma2_sim[i] <- omega + alpha * x[i - 1]^2 + beta * sigma2_sim[i - 1]
    x[i] <- mu+phi * x[i - 1] + sqrt(sigma2_sim[i]) * z[i]
    sigma_sim[i]<-sqrt(sigma2_sim[i])
  }
  
  R_t_matrix[, sim] <- x
  sigma2_matrix[, sim] <- sigma2_sim
  sigma_matrix[, sim] <- sigma_sim
}

for (sim in 1:M) {
  for (t in 1:T) {
    var_matrix[t, sim] <- VaR(alpha = 0.975, sigma2 = sigma_matrix[t, sim], mu = mu, df = df )*sqrt(sigma2_matrix[t, sim])
    es_matrix[t, sim]  <- ES(alpha = 0.975, sigma2 = sigma2_matrix[t, sim], mu = mu, df = df )*sqrt(sigma2_matrix[t, sim])
  }
}

```

Lav matricer for R_t, VaR og ES om til lister og finder de 10.000 værdier for teststørrelserne
```{r,warning=FALSE}
R_t_list <- split(R_t_matrix, col(R_t_matrix))
sigma_list<-split(sigma_matrix, col(sigma_matrix))
var_list <- split(var_matrix, col(var_matrix))
es_list <- split(es_matrix, col(es_matrix))

Z_1_results <- lapply(seq_along(R_t_list), function(i) {
  Z_1(R_t_list[[i]],var_list[[i]],es_list[[i]])$Z_1
  })
Z_1_results<-unlist(Z_1_results[!is.na(Z_1_results)])

Z_2_results <- lapply(seq_along(R_t_list), function(i) {
  Z_2(R_t_list[[i]],var_list[[i]],es_list[[i]])$Z_2
  })
Z_2_results<-unlist(Z_2_results)

cc_results <- lapply(seq_along(R_t_list), function(i) {
  cc_backtest(R_t_list[[i]],var_list[[i]],es_list[[i]],sigma_list[[i]], alpha=0.975)$min_t2
  })
cc_results<-unlist(cc_results)

esr_results <- lapply(seq_along(R_t_list), function(i) {
  esr_backtest(R_t_list[[i]],-es_list[[i]],alpha=0.975)$t0
  })
esr_results<-unlist(esr_results)

```
Visualisering af testsetup
```{r,warning=FALSE}
kombi_data <- tail(c(R_t_matrix[,1], SP500_Test_returns), 500)

df <- data.frame(
  Observationer = 1:500,
  Afkast = kombi_data
)

split_index <- 250
split_end <- 500

# Lav ggplot
ggplot(df, aes(x = Observationer, y = Afkast)) +
  geom_line(color = "blue", size = 0.5) +  # Tegn linjen for data
  geom_vline(xintercept = split_index, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = split_end, color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = split_index, y = min(df$Afkast), label = "02-02-2024", color = "red", hjust = 1.2, vjust = 1,size=3) + 
  annotate("text", x = split_end, y = min(df$Afkast), label = "31-01-2025", color = "red", hjust = 1.2, vjust = 1,size=3) + 
   annotate("text", x = split_index, y = max(df$Afkast), label = "250 dages AR(1)-GARCH(1,1) simulationer", color = "black", hjust = 1.01, vjust = 1,size=3.7) +
  annotate("text", x = split_end, y = max(df$Afkast), label = "S&P500 Test data", color = "black", hjust = 1.5, vjust = 1,size=3.7) + 
  labs(title = "GARCH(1,1) Simuleringer og S&P500",
       x = "Observationer",
       y = "Afkast") +
  theme_minimal()

```

Testværdier på empirisk data
```{r,warning=FALSE}
kombi_data<-tail(c(R_t_matrix[,1],SP500_Test_returns),500)

split_index <- 250
N<-500
var_hs<-rep(NA,T)
es_hs<-rep(NA,T)
varians<-rep(NA,T)
sd<-rep(NA,T)

  for (i in 251:N) {
    retwindow <- kombi_data[(i-T):(i-1)]
    
    var_hs[i] <- -quantile(retwindow, probs = 0.025)
    varians[i]<- var(retwindow)
    sd[i]<-sqrt(varians[i])
    
    es_hs[i] <- -mean(retwindow[retwindow < -var_hs[i]])
  }
sd_rec<-tail(sd,T)
var_hs_rec<-tail(var_hs,T)
es_hs_rec<-tail(es_hs,T)

cc_testværdi<-cc_backtest(r=SP500_Test_returns, q=var_hs_rec, e=es_hs_rec, s=sd_rec, alpha=0.975)$min_t2
esr_testværdi<- esr_backtest(r=SP500_Test_returns,e=-es_hs_rec,alpha=0.975)$t0
z2_testværdi<- Z_2(SP500_Test_returns,var_hs_rec,es_hs_rec )$Z_2
z1_testværdi<- Z_1(SP500_Test_returns,var_hs_rec,es_hs_rec  )$Z_1

df <- data.frame(Returns = Z_2_results)
ggplot(df, aes(x = Z_2_results)) +
  geom_density(alpha = 1, size=0.8, color="black", bw = 0.1) + 
  geom_vline(aes(xintercept = quantile(Z_2_results, 0.05), color = "Kritisk værdi 5% for Z2 ved AR(1)-GARCH(1,1)"), linetype = "dashed", size=0.8) +
  geom_vline(aes(xintercept = z2_testværdi, color = "Z2 Testværdi med empirisk data"), size=0.8) +
  scale_color_manual(values = c("Kritisk værdi 5% for Z2 ved AR(1)-GARCH(1,1)" = "black","Z2 Testværdi med empirisk data" = "purple")) +
  labs(title = "Sammenligning af kritisk værdi og empirisk testværdi for Z2",
       x = "Afkast",
       y = "Tæthed",
       color = " ") +
  theme_minimal() +
  theme(legend.position = "bottom")

df <- data.frame(Returns = Z_1_results)
ggplot(df, aes(x = Z_1_results)) +
  geom_density(alpha = 1, size=0.8, color="black") + 
  geom_vline(aes(xintercept = quantile(Z_1_results, 0.05), color = "Kritisk værdi 5% for Z1 ved AR(1)-GARCH(1,1)"), linetype = "dashed", size=0.8) +
  geom_vline(aes(xintercept = z1_testværdi, color = "Z1 Testværdi med empirisk data"), size=0.8) +
  scale_color_manual(values = c("Kritisk værdi 5% for Z1 ved AR(1)-GARCH(1,1)" = "black", "Z1 Testværdi med empirisk data" = "deeppink")) +
  labs(title = "Sammenligning af kritisk værdi og empirisk testværdi for Z1",
       x = "Afkast",
       y = "Tæthed",
       color = " ") +
  theme_minimal() +
  theme(legend.position = "bottom")

df <- data.frame(Returns = cc_results)
ggplot(df, aes(x = cc_results)) +
  geom_density(alpha = 1, size=0.8, color="black", bw=0.5) + 
  geom_vline(aes(xintercept = quantile(cc_results, 0.05), color = "Kritisk værdi 5% for CC ved AR(1)-GARCH(1,1)"), linetype = "dashed", size=0.8) +
  geom_vline(aes(xintercept = cc_testværdi, color = "CC Testværdi med empirisk data"), size=1) +
  scale_color_manual(values = c("Kritisk værdi 5% for CC ved AR(1)-GARCH(1,1)" = "black", "CC Testværdi med empirisk data" = "blue")) +
  labs(title = "Sammenligning af kritisk værdi og empirisk testværdi for CC",
       x = "Afkast",
       y = "Tæthed",
       color = " ") +
  theme_minimal() +
  theme(legend.position = "bottom")

df <- data.frame(Returns = esr_results)
ggplot(df, aes(x = esr_results)) +
  geom_density(alpha = 1, size=0.8, color="black") + 
  geom_vline(aes(xintercept = quantile(esr_results, 0.05), color = "Kritisk værdi 5% for ESR ved AR(1)-GARCH(1,1)"), linetype = "dashed", size=0.8) +
  geom_vline(aes(xintercept = esr_testværdi, color = "ESR Testværdi med empirisk data"), size=1) +
  scale_color_manual(values = c("Kritisk værdi 5% for ESR ved AR(1)-GARCH(1,1)" = "black","ESR Testværdi med empirisk data" = "red")) +
  labs(title = "Sammenligning af kritisk værdi og empirisk testværdi for Intercept ESR",
       x = "Afkast",
       y = "Tæthed",
       color = " ") +
  theme_minimal() +
  theme(legend.position = "bottom")

```




