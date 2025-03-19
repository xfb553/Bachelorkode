Pakker:
```{r}
library(ggplot2)
library(gridExtra)
library(rugarch)
library(esback)
library(readr)
library(dplyr)
library(quantmod)
library(fGarch)
```

Indlæsning af data til estination af parametre under AR(1)-GARCH(1,1):
```{r}
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

VaR og ES
```{r}
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
```{r}
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
```{r}
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

cc-backtest
```{r}
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

Esr-backtest:
```{r}
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

Power:
```{r}
Power <- function(data_null, data_alternatives, significance_level) {
  sorted_null <- sort(data_null)
  crit_value <- quantile(sorted_null, probs = significance_level)

  powers <- list(Significance_Level = significance_level)
    power <- mean(data_alternatives <= crit_value)

  return(as.data.frame(power))
}
```



Historisk simulation på AR(1)-GARCH(1,1):
```{r}
phi <- fitted_values_train["ar1"]
omega <- fitted_values_train["omega"]
alpha <- fitted_values_train["alpha1"]
beta <- fitted_values_train["beta1"]
mu <- fitted_values_train["mu"]
df <- fitted_values_train["shape"]
T <- 250
N <- 500
M <- 10000

sigma_matrix_hs<-matrix(NA, nrow = N, ncol = M)
x_matrix_hs <- matrix(NA, nrow = N, ncol = M)
var_forcast_matrix_hs<- matrix(NA, nrow = N, ncol = M)
es_forcast_matrix_hs<- matrix(NA, nrow = N, ncol = M)
var_matrix_hs <- matrix(NA, nrow = N - 250, ncol = M)
es_matrix_hs <- matrix(NA, nrow = N - 250, ncol = M)

for (sim in 1:M) {
  
  z <- rt(N,df=df)/sqrt(df/(df-2))
  
  x_hs <- rep(NA, N)
  sigma2_hs <- rep(NA, N)
  H_0_sigma_hs <- rep(NA, N)
  var_forcast_hs<-rep(NA,N)
  es_forcast_hs<-rep(NA,N)
  
  sigma2_hs[1] <- omega / (1 - alpha - beta)
  x_hs[1] <- sqrt(sigma2_hs[1]) * z[1]
  H_0_sigma_hs[1] <- sqrt(sigma2_hs[1])
  var_forcast_hs[1] <- VaR(alpha=0.975, df=df , Norm=FALSE) * H_0_sigma_hs[1]
  es_forcast_hs[1] <- ES(alpha=0.975, df=df , Norm=FALSE) * H_0_sigma_hs[1]

  for (i in 2:N) {
    sigma2_hs[i] <- omega + alpha * x_hs[i - 1]^2 + beta * sigma2_hs[i - 1]
    x_hs[i] <- mu+phi * x_hs[i - 1] + sqrt(sigma2_hs[i]) * z[i]
    H_0_sigma_hs[i] <- sqrt(sigma2_hs[i])
    var_forcast_hs[i] <- VaR(alpha=0.975, df=df , Norm=FALSE) * H_0_sigma_hs[i]
    es_forcast_hs[i] <- ES(alpha=0.975, df=df , Norm=FALSE) * H_0_sigma_hs[i]
  }
  
  x_matrix_hs[, sim] <- x_hs
  sigma_matrix_hs[,sim] <- H_0_sigma_hs
  var_forcast_matrix_hs[,sim]<-var_forcast_hs
  es_forcast_matrix_hs[,sim]<-es_forcast_hs

  for (i in 251:N) {
    retwindow <- x_hs[(i-T):(i-1)]
    
    var_matrix_hs[i - T,sim] <- -quantile(retwindow, probs = 0.025)
    
    es_matrix_hs[i - T,sim] <- -mean(retwindow[retwindow < -var_matrix_hs[i - T, sim]])
  }
}
x_recent <- tail(x_matrix_hs, T)
var0_recent<-tail(var_forcast_matrix_hs,T)
es0_recent<-tail(es_forcast_matrix_hs,T)
es1_hs_recent <- tail(es_matrix_hs, T)
var1_hs_recent <- tail(var_matrix_hs, T)
sigma_matrix_hs_recent<-tail(sigma_matrix_hs,T)
```

Power med HS:
```{r}
signi_values <- seq(0, 1, by = 0.001)

power_esr <- rep(NA,length(signi_values))
power_cc<- rep(NA,length(signi_values))
power_z1 <- rep(NA,length(signi_values))
power_z2 <- rep(NA,length(signi_values))
crit_esr <- rep(NA,length(signi_values))
crit_cc <- rep(NA,length(signi_values))
crit_z1 <- rep(NA,length(signi_values))
crit_z2 <- rep(NA,length(signi_values))

  results_H1_esr <- rep(NA, M)
  results_H0_esr <- rep(NA, M)
  results_H1_cc <- rep(NA, M)
  results_H0_cc <- rep(NA, M)
  results_H1_z1 <- rep(NA, M)
  results_H0_z1 <- rep(NA, M)
  results_H1_z2 <- rep(NA, M)
  results_H0_z2 <- rep(NA, M)
  

  for (i in 1:M) {
    results_H0_esr[[i]] <- esr_backtest(r=x_recent[,i], e=es0_recent[,i],  alpha=0.975 )$t0
    
    results_H1_esr[[i]] <- esr_backtest(r=x_recent[,i], e=es1_hs_recent[,i],  alpha=0.975 )$t0
    
    results_H0_cc[[i]] <- cc_backtest(r = x_recent[,i], q = var0_recent[,i], e = es0_recent[,i], alpha = 0.975, s = sigma_matrix_hs_recent[,i])$min_t2
    
    results_H1_cc[[i]] <- cc_backtest(r = x_recent[,i], q = var1_hs_recent[,i], e = es1_hs_recent[,i], alpha = 0.975, s = sigma_matrix_hs_recent[,i])$min_t2
    
    results_H0_z1[[i]] <- Z_1(x_recent[,i],var0_recent[,i],es0_recent[,i])$Z_1
    results_H0_z1<-results_H0_z1[!is.na(results_H0_z1)]
    
    results_H1_z1[[i]] <- Z_1(x_recent[,i],var1_hs_recent[,i],es1_hs_recent[,i])$Z_1
    results_H1_z1<-results_H1_z1[!is.na(results_H1_z1)]
    
    results_H0_z2[[i]] <- Z_2(x_recent[,i],var0_recent[,i],es0_recent[,i])$Z_2
    
    results_H1_z2[[i]] <- Z_2(x_recent[,i],var1_hs_recent[,i],es1_hs_recent[,i])$Z_2
    
  }

  
for (s in seq_along(signi_values)) {
  
  signi <- signi_values[s] 

  power_esr[s] <- Power(results_H0_esr, results_H1_esr, signi)
  power_esr<-unlist(power_esr)
  
  crit_esr[s] <- Power(results_H0_esr, results_H0_esr, signi)
  crit_esr<-unlist(crit_esr)
  
  power_cc[s] <- Power(results_H0_cc, results_H1_cc, signi)
  power_cc<-unlist(power_cc)
  
  crit_cc[s] <- Power(results_H0_cc, results_H0_cc, signi)
  crit_cc<-unlist(crit_cc)
  
  power_z1[s] <- Power(results_H0_z1,results_H1_z1, signi)
  power_z1<-unlist(power_z1)
  
  crit_z1[s] <- Power(results_H0_z1, results_H0_z1, signi)
  crit_z1<-unlist(crit_z1)
  
  power_z2[s] <- Power(results_H0_z2, results_H1_z2, signi)
  power_z2<-unlist(power_z2)
  
  crit_z2[s] <- Power(results_H0_z2, results_H0_z2, signi)
  crit_z2<-unlist(crit_z2)

}
resultater<-data.frame(signi=signi_values, power_esr=power_esr,crit_esr=crit_esr,power_cc=power_cc,crit_cc=crit_cc, power_z1=power_z1,crit_z1=crit_z1,power_z2=power_z2,crit_z2=crit_z2)  
```

Plot for HS:
```{r}
power_hs <- data.frame(
  crit_values_hs = c(resultater$crit_esr, resultater$crit_cc, resultater$crit_z1, resultater$crit_z2),
  power_hs = c(resultater$power_esr, resultater$power_cc, resultater$power_z1, resultater$power_z2),
  Test = rep(c("ESR Test", "CC Test", "Z1 Test", "Z2 Test"), c(length(crit_esr), length(crit_cc), length(crit_z1), length(crit_z2)))
)

ggplot(power_hs, aes(x = crit_values_hs, y = power_hs, color = Test))+
  geom_line(size = 1) +
  labs(title = "Power vs Empiriske størrelser med HS",
       x = "Empirisk størrelse",
       y = "Power",
       color = "Test Type") +
  theme_minimal()+
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("ESR Test" = "red", "CC Test" = "blue", "Z1 Test"="deeppink", "Z2 Test"="purple"))+geom_abline(slope = 1, intercept = 0, color = "black")+
  coord_cartesian(ylim = c(0, 1))

```


Filtreret historisk simulation på AR(1)-GARCH(1,1):
```{r}
phi <- fitted_values_train["ar1"]
omega <- fitted_values_train["omega"]
alpha <- fitted_values_train["alpha1"]
beta <- fitted_values_train["beta1"]
mu <- fitted_values_train["mu"]
df <- fitted_values_train["shape"]
T <- 250
N <- 500
M <- 1000

sigma_fhs_matrix<-matrix(NA, nrow = N, ncol = M)
x_fhs_matrix <- matrix(NA, nrow = N, ncol = M)
var1_fhs_matrix <- matrix(NA, nrow = N - 250, ncol = M)
es1_fhs_matrix <- matrix(NA, nrow = N - 250, ncol = M)
var0_fhs_matrix <-  matrix(NA, nrow = N, ncol = M)
es0_fhs_matrix <-  matrix(NA, nrow = N, ncol = M)

for (sim in 1:M) {

  z <- rt(N,df=df)/sqrt(df/(df-2))  
  
  x_fhs <- rep(NA, N)
  sigma2_fhs <- rep(NA, N)
  H_0_sigma_fhs <- rep(NA, N)
  var0_fhs_forcast<-rep(NA,N)
  es0_fhs_forcast<-rep(NA,N)
  
  sigma2_fhs[1] <- omega / (1 - alpha - beta)
  x_fhs[1] <- sqrt(sigma2_fhs[1]) * z[1]
  H_0_sigma_fhs[1] <- sqrt(sigma2_fhs[1])
  var0_fhs_forcast[1] <- VaR(alpha=0.975, df=df , Norm=FALSE) * H_0_sigma_fhs[1]
  es0_fhs_forcast[1] <- ES(alpha=0.975, df=df , Norm=FALSE) * H_0_sigma_fhs[1]

  for (i in 2:N) {
    sigma2_fhs[i] <- omega + alpha * x_fhs[i - 1]^2 + beta * sigma2_fhs[i - 1]
    x_fhs[i] <- mu+phi * x_fhs[i - 1] + sqrt(sigma2_fhs[i]) * z[i]
    H_0_sigma_fhs[i] <- sqrt(sigma2_fhs[i])
    var0_fhs_forcast[i] <- VaR(alpha=0.975, df=df , Norm=FALSE) * H_0_sigma_fhs[i]
    es0_fhs_forcast[i] <- ES(alpha=0.975, df=df , Norm=FALSE) * H_0_sigma_fhs[i]
  }

  x_fhs_matrix[, sim] <- x_fhs
  sigma_fhs_matrix[,sim] <- H_0_sigma_fhs
  var0_fhs_matrix[,sim]<-var0_fhs_forcast
  es0_fhs_matrix[,sim]<-es0_fhs_forcast

  for (i in 251:N) {
    retwindow <- x_fhs[(i-T):(i-1)]
    
  argarch_spec <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
    mean.model = list(armaOrder = c(1,0)),
    distribution.model = "std",  
    fixed.pars = list(omega = omega)
  )

  fit_argarch <- tryCatch(
    ugarchfit(spec = argarch_spec, data = retwindow, solver = "hybrid"),
    error = function(e) return(NA)
  )

  if (inherits(fit_argarch, "try-error")) next

  sigma_forecast_argarch <- as.numeric(ugarchforecast(fit_argarch, n.ahead=1)@forecast$sigmaFor[1])

  sigma_t_argarch <- tail(fit_argarch@fit$sigma, length(retwindow))
  ret_standardized_argarch <- retwindow / sigma_t_argarch

  var1_fhs_matrix[i - T, sim] <- - ret_standardized_argarch*quantile(sigma_forecast_argarch, probs=0.025)

  threshold_std <- quantile(ret_standardized_argarch, probs=0.025)  

  es_obs <- ret_standardized_argarch[ret_standardized_argarch < threshold_std]
  
if (length(es_obs) > 2) {
  es1_fhs_matrix[i - T, sim] <- -sigma_forecast_argarch * mean(es_obs)
} else {
  es1_fhs_matrix[i - T, sim] <- NA
}
    
  }
}
x_fhs_recent <- tail(x_fhs_matrix, T)
var0_fhs_recent<-tail(var0_fhs_matrix,T)
es0_fhs_recent<-tail(es0_fhs_matrix,T)
es1_fhs_recent <- tail(es1_fhs_matrix, T)
var1_fhs_recent <- tail(var1_fhs_matrix, T)
sigma_fhs_matrix_recent<-tail(sigma_fhs_matrix,T)
```


Power med FHS:
```{r}
signi_values <- seq(0, 1, by = 0.0001)

power_esr_fhs <- numeric(length(signi_values))
power_cc_fhs <- numeric(length(signi_values))
power_z1_fhs <- numeric(length(signi_values))
power_z2_fhs <- numeric(length(signi_values))
crit_esr_fhs <- numeric(length(signi_values))
crit_cc_fhs <- numeric(length(signi_values))
crit_z2_fhs <- numeric(length(signi_values))
crit_z1_fhs <- numeric(length(signi_values))

  H1_esr_fhs <- rep(NA, M)
  H0_esr_fhs <- rep(NA, M)
  H1_cc_fhs <- rep(NA, M)
  H0_cc_fhs <- rep(NA, M)
  H1_z1_fhs <- rep(NA, M)
  H0_z1_fhs <- rep(NA, M)
  H1_z2_fhs <- rep(NA, M)
  H0_z2_fhs <- rep(NA, M)

  for (i in 1:M) {
    H0_esr_fhs[[i]] <- esr_backtest(r=x_fhs_recent[,i], e=-es0_fhs_recent[,i],  alpha=0.975 )$t0
    
    H1_esr_fhs[[i]] <- esr_backtest(r=x_fhs_recent[,i], e=-es1_fhs_recent[,i],  alpha=0.975 )$t0
    
    H0_cc_fhs[[i]] <- cc_backtest(r = x_fhs_recent[,i], q = var0_fhs_recent[,i], e = es0_fhs_recent[,i], alpha = 0.975, s = sigma_fhs_matrix_recent[,i])$min_t2
    
    H1_cc_fhs[[i]] <- cc_backtest(r = x_fhs_recent[,i], q = var1_fhs_recent[,i], e = es1_fhs_recent[,i], alpha = 0.975, s = sigma_fhs_matrix_recent[,i])$min_t2
    
    H0_z1_fhs[[i]] <- Z_1(x_fhs_recent[,i],var0_fhs_recent[,i],es0_fhs_recent[,i])$Z_1
    H0_z1_fhs<-H0_z1_fhs[!is.na(H0_z1_fhs)]
    
    H1_z1_fhs[[i]] <- Z_1(x_fhs_recent[,i],var1_fhs_recent[,i],es1_fhs_recent[,i])$Z_1
    H1_z1_fhs<-H1_z1_fhs[!is.na(H1_z1_fhs)]
    
    H0_z2_fhs[[i]] <- Z_2(x_fhs_recent[,i],var0_fhs_recent[,i],es0_fhs_recent[,i])$Z_2
    
    H1_z2_fhs[[i]] <- Z_2(x_fhs_recent[,i],var1_fhs_recent[,i],es1_fhs_recent[,i])$Z_2
  }

for (s in seq_along(signi_values)) {
  
  signi <- signi_values[s]
 
  power_esr_fhs[s] <- Power(H0_esr_fhs, H1_esr_fhs, signi)
  power_esr_fhs<-unlist(power_esr_fhs)
  
  crit_esr_fhs[s] <- Power(H0_esr_fhs, H0_esr_fhs, signi)
  crit_esr_fhs<-unlist(crit_esr_fhs)
  
  power_cc_fhs[s] <- Power(H0_cc_fhs, H1_cc_fhs, signi)
  power_cc_fhs<-unlist(power_cc_fhs)

  crit_cc_fhs[s] <- Power(H0_cc_fhs, H0_cc_fhs, signi)
  crit_cc_fhs<-unlist(crit_cc_fhs)
  
  power_z1_fhs[s] <- Power(H0_z1_fhs, H1_z1_fhs, signi)
  power_z1_fhs<-unlist(power_z1_fhs)

  crit_z1_fhs[s] <- Power(H0_z1_fhs, H0_z1_fhs, signi)
  crit_z1_fhs<-unlist(crit_z1_fhs)
  
  power_z2_fhs[s] <- Power(H0_z2_fhs, H1_z2_fhs, signi)
  power_z2_fhs<-unlist(power_z2_fhs)

  crit_z2_fhs[s] <- Power(H0_z2_fhs, H0_z2_fhs, signi)
  crit_z2_fhs<-unlist(crit_z2_fhs)
  
}
df_power_fhs<-data.frame(signi=signi_values, power_cc_fhs=power_cc_fhs, crit_cc_fhs=crit_cc_fhs, power_esr_fhs=power_esr_fhs, crit_esr_fhs=crit_esr_fhs, power_z1_fhs=power_z1_fhs, crit_z1_fhs=crit_z1_fhs, power_z2_fhs=power_z2_fhs, crit_z2_fhs=crit_z2_fhs)
```


Plot for FHS:
```{r}
power_fhs_ <- data.frame(
  crit_fhs = c(df_power_fhs$crit_cc_fhs, df_power_fhs$crit_esr_fhs, df_power_fhs$crit_z1_fhs, df_power_fhs$crit_z2_fhs),
  power_fhs = c(df_power_fhs$power_cc_fhs, df_power_fhs$power_esr_fhs,df_power_fhs$power_z1_fhs, df_power_fhs$power_z2_fhs),
  Test = rep(c("CC Test","ESR Test","Z1 Test", "Z2 Test"), c(length(crit_esr_fhs), length(crit_cc_fhs), length(crit_z1_fhs), length(crit_z2_fhs))
))

ggplot(power_fhs_, aes(x = crit_fhs, y = power_fhs, color = Test)) +
  geom_line(size = 1) +
  labs(title = "Power vs empiriske størrelser for ESR og CC Tests med FHS",
       x = "Empirisk størrelse",
       y = "Power",
       color = "Test Type") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("ESR Test" = "red", "CC Test" = "blue", "Z1 Test"="deeppink", "Z2 Test"="purple"))+geom_abline(slope = 1, intercept = 0, color = "black")+
  coord_cartesian(ylim = c(0, 1)) 
```


Plotter tætheder under H0 og H1 ved HS:
```{r}
M <- 10000
cc0 <- rep(NA,M)
cc1 <- rep(NA,M)
esr0<-rep(NA,M)
esr1<-rep(NA,M)
z20<-rep(NA,M)
z21<-rep(NA,M)
z10<-rep(NA,M)
z11<-rep(NA,M)

for (i in 1:M){
cc0[[i]]<-cc_backtest(x_recent[,i], var0_recent[,i], es0_recent[,i], s=sigma_matrix_hs_recent[,i], alpha=0.975)$min_t2

cc1[[i]]<-cc_backtest(x_recent[,i], var1_hs_recent[,i], es1_hs_recent[,i], s=sigma_matrix_hs_recent[,i], alpha=0.975)$min_t2

esr0[[i]]<-esr_backtest(r=x_recent[,i], e=-es0_recent[,i],  alpha=0.975)$t0

esr1[[i]]<-esr_backtest(r=x_recent[,i], e=-es1_hs_recent[,i],  alpha=0.975)$t0

z20[[i]]<-Z_2(x_recent[,i],var0_recent[,i],es0_recent[,i])$Z_2

z21[[i]]<-Z_2(x_recent[,i],var1_hs_recent[,i],es1_hs_recent[,i])$Z_2

z10[[i]]<-Z_1(x_recent[,i],var0_recent[,i],es0_recent[,i])$Z_1

z11[[i]]<-Z_1(x_recent[,i],var1_hs_recent[,i],es1_hs_recent[,i])$Z_1


}

data_cc <- data.frame(
  Value_cc = c(cc0, cc1),
  Type = rep(c("CC H0", "CC H1"), c(length(cc0), length(cc1)))
)
data_esr <- data.frame(
  Value_esr = c(esr0, esr1),
  Type = rep(c("Intercept ESR under H0", "Intercept ESR under H1"), c(length(esr0), length(esr1)))
)
data_z2 <- data.frame(
  Value_z2 = c(z20, z21),
  Type= rep(c("Z2 H0", "Z2 H1"), c(length(z20), length(z21)))
)
data_z1 <- data.frame(
  Value_z1 = c(z10, z11),
  Type= rep(c("Z1 H0", "Z1 H1"), c(length(z10), length(z11)))
)

ggplot(data_cc, aes(x = Value_cc, color = Type)) +
  geom_density(size = 1, bw=0.8) +
  labs(title = "Tætheden for CC backtest under H0 og H1 ved HS",
       x = "Value", y = "Density") +
  scale_color_manual(values = c("CC H0" = "black", "CC H1" = "blue")) +
    scale_x_continuous(limits = c(-10, 15)) +
  theme_minimal()+
  theme(legend.position = "bottom")

ggplot(data_esr, aes(x = Value_esr, color = Type)) +
  geom_density(size = 1) +
  labs(title = "Tætheden for Intercept ESR backtest under H0 og H1 ved HS",
       x = "Value", y = "Density") +
  scale_color_manual(values = c("Intercept ESR under H0" = "black", "Intercept ESR under H1" = "red")) +
  scale_x_continuous(limits = c(-3, 6)) +
  theme_minimal()+
  theme(legend.position = "bottom")

ggplot(data_z2, aes(x = Value_z2, color = Type, bw = 0.1)) +
  geom_density(size = 1) +
  labs(title = "Tætheden for Z2 backtest under H0 og H1 ved HS",
       x = "Value", y = "Density") +
  scale_color_manual(values = c("Z2 H0" = "black", "Z2 H1" = "purple")) +
  scale_x_continuous(limits = c(-6, 1)) +
  theme_minimal()+
  theme(legend.position = "bottom")

  ggplot(data_z1, aes(x = Value_z1, color = Type)) +
  geom_density(size = 1) +
  labs(title = "Tætheden for Z1 backtest under H0 og H1 ved HS",
       x = "Value", y = "Density") +
  scale_color_manual(values = c("Z1 H0" = "black", "Z1 H1" = "deeppink")) +
      scale_x_continuous(limits = c(-0.75, 0.75))+
  theme_minimal()+
  theme(legend.position = "bottom")

```


Plotter tætheder under H0 og H1 ved FHS:
```{r}
M <- M <- 1000
cc0_f <- rep(NA,M)
cc1_f <- rep(NA,M)
esr0_f<-rep(NA,M)
esr1_f<-rep(NA,M)
z20_f<-rep(NA,M)
z21_f<-rep(NA,M)
z10_f<-rep(NA,M)
z11_f<-rep(NA,M)

for (i in 1:M){
cc0_f[[i]]<-cc_backtest(x_fhs_recent[,i], var0_fhs_recent[,i], es0_fhs_recent[,i], s=sigma_fhs_matrix_recent[,i], alpha=0.975)$min_t2
cc1_f[[i]]<-cc_backtest(x_fhs_recent[,i], var1_fhs_recent[,i], es1_fhs_recent[,i], s=sigma_fhs_matrix_recent[,i], alpha=0.975)$min_t2
esr0_f[[i]]<-esr_backtest(r=x_fhs_recent[,i], e=-es0_fhs_recent[,i],  alpha=0.975 )$t0
esr1_f[[i]]<-esr_backtest(r=x_fhs_recent[,i], e=-es1_fhs_recent[,i],  alpha=0.975 )$t0
z20_f[[i]]<-Z_2(x_fhs_recent[,i],var0_fhs_recent[,i],es0_fhs_recent[,i])$Z_2
z21_f[[i]]<-Z_2(x_fhs_recent[,i],var1_fhs_recent[,i],es1_fhs_recent[,i])$Z_2
z10_f[[i]]<-Z_1(x_fhs_recent[,i],var0_fhs_recent[,i],es0_fhs_recent[,i])$Z_1
z11_f[[i]]<-Z_1(x_fhs_recent[,i],var1_fhs_recent[,i],es1_fhs_recent[,i])$Z_1
}

data_cc_f <- data.frame(
  Value_cc_f = c(cc0_f, cc1_f),
  Type = rep(c("CC H0", "CC H1"), c(length(cc0_f), length(cc1_f)))
)
data_esr_f <- data.frame(
  Value_esr_f = c(esr0_f, esr1_f),
  Type = rep(c("Intercept ESR under H0", "Intercept ESR under H1"), c(length(esr0_f), length(esr1_f)))
)
data_z2_f <- data.frame(
  Value_z2_f = c(z20_f, z21_f),
  Type = rep(c("Z2 under H0", "Z2 under H1"), c(length(z20_f), length(z21_f)))
)
data_z1_f <- data.frame(
  Value_z1_f = c(z10_f, z11_f),
  Type = rep(c("Z1 under H0", "Z1 under H1"), c(length(z10_f), length(z11_f)))
)

ggplot(data_cc_f, aes(x = Value_cc_f, color = Type)) +
  geom_density(size = 1.2, bw = 0.5) +
  labs(title = "Tætheden for CC backtest under H0 og H1 ved FHS",
       x = "Value", y = "Density") +
  scale_color_manual(values = c("CC H0" = "black", "CC H1" = "blue"))+
  scale_x_continuous(limits = c(-10, 15)) +
  theme_minimal()+
  theme(legend.position = "bottom")
  
ggplot(data_esr_f, aes(x = Value_esr_f, color = Type)) +
  geom_density(size = 1.2) +
  labs(title = "Tætheden for Intercept ESR backtest under H0 og H1 ved FHS",
       x = "Value", y = "Density") +
  scale_color_manual(values = c("Intercept ESR under H0" = "black", "Intercept ESR under H1" = "red"))+
  scale_x_continuous(limits = c(-3, 6)) +
  theme_minimal()+
  theme(legend.position = "bottom")

ggplot(data_z2_f, aes(x = Value_z2_f, color = Type)) +
  geom_density(size = 1.2) +
  labs(title = "Tætheden for Z2 backtest under H0 og H1 ved FHS",
       x = "Value", y = "Density") +
  scale_color_manual(values = c("Z2 under H0" = "black", "Z2 under H1" = "purple"))+
  scale_x_continuous(limits = c(-6, 1), breaks = seq(-5, 5, by = 1)) +
  theme_minimal()+
  theme(legend.position = "bottom")
  
ggplot(data_z1_f, aes(x = Value_z1_f, color = Type)) +
  geom_density(size = 1.2) +
  labs(title = "Tætheden for Z1 backtest under H0 og H1 ved FHS",
       x = "Value", y = "Density") +
  scale_color_manual(values = c("Z1 under H0" = "black", "Z1 under H1" = "deeppink"))+
      scale_x_continuous(limits = c(-0.75, 0.75))+
  theme_minimal()+
  theme(legend.position = "bottom")
```
