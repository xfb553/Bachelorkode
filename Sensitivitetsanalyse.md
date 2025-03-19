Pakker
```{r}
library(esreg)
library(stats)
library(ggplot2)
library(dplyr)
library(rugarch)
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

Test 3
```{r,warning=FALSE}
Z_3 <- function(R_t) {

  T<-250
  TAlpha <- round(T * 0.025)

    integral <- -T / TAlpha * integrate(
      function(p) pbeta(1 - p, shape1 = T - TAlpha, shape2 = TAlpha) * quantile(R_t,p),
      lower = 0, upper = 1
    )$value
  
  
    U <- runif(T)
      PU<-quantile(R_t,U)
  
    
    PU <- sort(PU)
    q <- sum(PU[1:TAlpha])
    ES_hat <- -q / TAlpha
    
  Z_3 <- -ES_hat / integral + 1
  
  return(list(Z_3 = Z_3))
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

Rejection_rate
```{r,warning=FALSE}
rejection_rate <- function(data_null, data_alternatives, significance_level,sim) {
  data_alternatives<-data_alternatives[!is.na(data_alternatives)]
  sorted_null <- sort(data_null)
  crit_value <- quantile(sorted_null, probs = significance_level,na.rm = TRUE)

  powers <- list(Significance_Level = significance_level)
    power <- sum(data_alternatives <= crit_value)/sim

  return(as.data.frame(power))
}

#Standard antagelser
T <- 250
M <- 10000
```


Indlæsning af data
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

Ændringer af alpha
```{r,warning=FALSE}
alpha_values <- seq(0.03,0.25, 0.02) 

k_alpha <- function(alpha_values, days = T, simulations = M) {

fitted_values_train <- coef(SP500_ARGARCH_fit)
omega <- fitted_values_train["omega"]
alpha <- fitted_values_train["alpha1"]
beta <- fitted_values_train["beta1"]
alpha_beta_sum<-alpha+beta
mu <- fitted_values_train["mu"]
df <- fitted_values_train["shape"]
phi <- fitted_values_train["ar1"]
  
  cc_p_values_matrix <- matrix(NA, nrow = simulations, ncol = length(alpha_values))
  cc_rejection_rates <- numeric(length(alpha_values))
  esr_p_values_matrix <- matrix(NA, nrow = simulations, ncol = length(alpha_values))
  esr_rejection_rates <- numeric(length(alpha_values))
  Z1_A_matrix<-matrix(NA, nrow = simulations, ncol = length(alpha_values))
  Z1_rejection_rates<-numeric(length(alpha_values))
  Z2_A_matrix<-matrix(NA, nrow = simulations, ncol = length(alpha_values))
  Z2_rejection_rates<-numeric(length(alpha_values))
  Z3_A_matrix<-matrix(NA, nrow = simulations, ncol = length(alpha_values))
  Z3_rejection_rates<-numeric(length(alpha_values))
  
  for (sim in 1:M) {
    
    z <- sqrt((df - 2) / df) * rt(T, df)
    
    x <- rep(NA, T)
    sigma2 <- rep(NA, T)
    H_0_sigma<-rep(NA,T)
    var_forcast<-rep(NA,T)
    es_forcast<-rep(NA,T)
    sigma2[1] <- omega / (1 - alpha - beta)
    x[1] <- sqrt(sigma2[1]) * z[1]
    H_0_sigma[1] <- sqrt(sigma2[1])
    var_forcast[1] <- VaR(alpha=0.975, df=df, Norm=FALSE) * H_0_sigma[1]
    es_forcast[1] <- ES(alpha=0.975, df=df, Norm=FALSE) * H_0_sigma[1]
    
    for (i in 2:T) {
      sigma2[i] <- omega + alpha * x[i - 1]^2 + beta * sigma2[i - 1]
      x[i] <- mu+phi * x[i - 1] + sqrt(sigma2[i]) * z[i]
      H_0_sigma[i] <- sqrt(sigma2[i])
      var_forcast[i] <- VaR(alpha=0.975,df=df, Norm=FALSE) * H_0_sigma[i]
      es_forcast[i] <- ES(alpha=0.975, df=df, Norm=FALSE) * H_0_sigma[i]
    }
    
    for (j in seq_along(alpha_values)) {
      alpha_A <- alpha_values[j]
      beta_A <- alpha_beta_sum - alpha_A
      
      x_A <- rep(NA, T)
      sigma2_A <- rep(NA, T)
      sigma2_A[1] <- omega / (1 - alpha_A - beta_A)
      x_A[1] <- sqrt(sigma2_A[1]) * z[1]
      
      for (i in 2:T) {
        sigma2_A[i] <- omega + alpha_A * x_A[i - 1]^2 + beta_A * sigma2_A[i - 1]
        x_A[i] <- mu+phi * x_A[i - 1] + sqrt(sigma2_A[i]) * z[i]
      }
      
      cc_p_value <- cc_backtest(r=x_A, q=var_forcast, e=es_forcast, s=H_0_sigma, alpha=0.975)$pvalue_onesided_general
      cc_p_values_matrix[sim, j] <- cc_p_value
      esr_p_value <- esr_backtest(r=x_A, e=-es_forcast, alpha=0.975)$pvalue_onesided_asymptotic
      esr_p_values_matrix[sim, j] <- esr_p_value
      
      Z1_A<-Z_1(x_A,var_forcast,es_forcast, quan = 0.05)$Z_1
      Z1_A_matrix[sim, j]<-Z1_A
      Z1<-Z_1(x,var_forcast,es_forcast, quan = 0.05)$Z_1
      Z1 <- Z1[!is.na(Z1)]
      Z2_A<-Z_2(x_A,var_forcast,es_forcast, quan = 0.05)$Z_2
      Z2_A_matrix[sim, j]<-Z2_A
      Z2<-Z_2(x,var_forcast,es_forcast, quan = 0.05)$Z_2
    }
  }

  for (j in seq_along(alpha_values)) {
    cc_rejection_rates[j] <- sum(cc_p_values_matrix[, j] < 0.05, na.rm = TRUE) / simulations
    esr_rejection_rates[j] <- sum(esr_p_values_matrix[, j] < 0.05, na.rm = TRUE) / simulations
    Z1_rejection_rates[j]<-rejection_rate(Z1, Z1_A_matrix[,j], 0.05, simulations)
    Z2_rejection_rates[j]<-rejection_rate(Z2, Z2_A_matrix[,j], 0.05, simulations)
  }

plot(alpha_values, cc_rejection_rates, type = "o", pch = 16, col = "blue", 
       xlab = expression(alpha), ylab = "Rejection Rate", 
       main = "Rejection Rate vs. ARCH Parameter (Alpha)", ylim = c(0, 1))
lines(alpha_values, esr_rejection_rates, type = "o", pch = 16, col = "red")
lines(alpha_values, Z1_rejection_rates, type = "o", pch = 16, col = "deeppink")
lines(alpha_values, Z2_rejection_rates, type = "o", pch = 16, col = "purple")
abline(v = fitted_values_train["alpha1"], col = "darkgray", lty = 2, lwd = 2)

legend("topleft", legend = c("CC Backtest", "ESR Backtest","Z1 backtest","Z2 backtest" ,"Alpha Reference"), 
       col = c("blue", "red","deeppink", "purple", "darkgray"), pch = c(16, 16, 16, 16,NA), lty = c(1, 1,1, 1,2))

}
k_alpha(alpha_values, 250)
```
Ændringer i phi
```{r, warning=FALSE}
phi_values <- seq(-0.3,0.5, 0.1)  

k_phi <- function(phi_values, days = T, simulations = M) {
  
fitted_values_train <- coef(SP500_ARGARCH_fit)
omega <- fitted_values_train["omega"]
alpha <- fitted_values_train["alpha1"]
beta <- fitted_values_train["beta1"]
df <- fitted_values_train["shape"]
mu <- fitted_values_train["mu"]
phi <- fitted_values_train["ar1"]
  
  cc_p_values_matrix <- matrix(NA, nrow = simulations, ncol = length(phi_values))
  cc_rejection_rates <- numeric(length(phi_values))
  esr_p_values_matrix <- matrix(NA, nrow = simulations, ncol = length(phi_values))
  esr_rejection_rates <- numeric(length(phi_values))
  Z2_A_matrix<-matrix(NA, nrow = simulations, ncol = length(phi_values))
  Z2_rejection_rates<-numeric(length(phi_values))
  Z1_A_matrix<-matrix(NA, nrow = simulations, ncol = length(phi_values))
  Z1_rejection_rates<-numeric(length(phi_values))
  
  for (sim in 1:M) {
    
    z <- sqrt((df - 2) / df) * rt(T, df)
    
    x <- rep(NA, T)
    sigma2 <- rep(NA, T)
    H_0_sigma<-rep(NA,T)
    var_forcast<-rep(NA,T)
    es_forcast<-rep(NA,T)
    sigma2[1] <- omega / (1 - alpha - beta)
    x[1] <- sqrt(sigma2[1]) * z[1]
    H_0_sigma[1] <- sqrt(sigma2[1])
    var_forcast[1] <- VaR(alpha=0.975,df=df, Norm=FALSE) * H_0_sigma[1]
    es_forcast[1] <- ES(alpha=0.975, df=df, Norm=FALSE) * H_0_sigma[1]
    
    for (i in 2:T) {
      sigma2[i] <- omega + alpha * x[i - 1]^2 + beta * sigma2[i - 1]
      x[i] <- mu+phi * x[i - 1] + sqrt(sigma2[i]) * z[i]
      H_0_sigma[i] <- sqrt(sigma2[i])
      var_forcast[i] <- VaR(alpha=0.975, df=df, Norm=FALSE) * H_0_sigma[i]
      es_forcast[i] <- ES(alpha=0.975, df=df, Norm=FALSE) * H_0_sigma[i]
    }
    
    for (j in seq_along(phi_values)) {
      phi_A <- phi_values[j]
      
      x_A <- rep(NA, T)
      sigma2_A <- rep(NA, T)
      sigma2_A[1] <- omega / (1 - alpha - beta)
      x_A[1] <- sqrt(sigma2_A[1]) * z[1]
      
      for (i in 2:T) {
        sigma2_A[i] <- omega + alpha * x_A[i - 1]^2 + beta * sigma2_A[i - 1]
        x_A[i] <- mu+phi_A * x_A[i - 1] + sqrt(sigma2_A[i]) * z[i]
      }
      
      cc_p_value <- cc_backtest(r=x_A, q=var_forcast, e=es_forcast, s=H_0_sigma, alpha=0.975)$pvalue_onesided_general
      cc_p_values_matrix[sim, j] <- cc_p_value
      esr_p_value <- esr_backtest(r=x_A, e=-es_forcast, alpha=0.975)$pvalue_onesided_asymptotic
      esr_p_values_matrix[sim, j] <- esr_p_value
      
      Z2_A<-Z_2(x_A,var_forcast,es_forcast, quan = 0.05)$Z_2
      Z2_A_matrix[sim, j]<-Z2_A
      Z2<-Z_2(x,var_forcast,es_forcast, quan = 0.05)$Z_2
      Z1_A<-Z_1(x_A,var_forcast,es_forcast, quan = 0.05)$Z_1
      Z1_A_matrix[sim, j]<-Z1_A
      Z1<-Z_1(x,var_forcast,es_forcast, quan = 0.05)$Z_1
    }
  }
  
  for (j in seq_along(phi_values)) {
    cc_rejection_rates[j] <- sum(cc_p_values_matrix[, j] < 0.05, na.rm = TRUE) / simulations
    esr_rejection_rates[j] <- sum(esr_p_values_matrix[, j] < 0.05, na.rm = TRUE) / simulations
    Z2_rejection_rates[j]<-rejection_rate(Z2, Z2_A_matrix[,j], 0.05, simulations)
    Z1_rejection_rates[j]<-rejection_rate(Z1, Z1_A_matrix[,j], 0.05, simulations)
  }

plot(phi_values, cc_rejection_rates, type = "o", pch = 16, col = "blue", 
       xlab = expression(phi), ylab = "Rejection Rate", 
       main = "Rejection Rate vs. AR Parameter (phi)", ylim = c(0, 1))
lines(phi_values, esr_rejection_rates, type = "o", pch = 16, col = "red")
lines(phi_values, Z2_rejection_rates, type = "o", pch = 16, col = "purple")
lines(phi_values, Z1_rejection_rates, type = "o", pch = 16, col = "deeppink")
abline(v = fitted_values_train["ar1"], col = "darkgray", lty = 2, lwd = 2)

# Tilføj en forklarende legend
legend("topleft", legend = c("CC Backtest", "ESR Backtest", "Z1 backtest","Z2 backtest" ,"Phi Reference"), 
       col = c("blue", "red", "deeppink", "purple", "darkgray"), pch = c(16, 16, 16, 16,NA), lty = c(1, 1,1, 1,2))

}
k_phi(phi_values, 250)
```

Ændringer i omega
```{r, warning=FALSE}
omega_values <- seq(0.004,0.1, length.out = 12)
k_omega <- function(omega_values, days=T, simulations = M) {

fitted_values_train <- coef(SP500_ARGARCH_fit)
omega <- fitted_values_train["omega"]
alpha <- fitted_values_train["alpha1"]
beta <- fitted_values_train["beta1"]
df <- fitted_values_train["shape"]
mu <- fitted_values_train["mu"]
phi <- fitted_values_train["ar1"]
 
  cc_p_values_matrix <- matrix(NA, nrow = M, ncol = length(omega_values))
  cc_rejection_rates <- numeric(length(omega_values))
  esr_p_values_matrix <- matrix(NA, nrow = M, ncol = length(omega_values))
  esr_rejection_rates <- numeric(length(omega_values))
  Z2_A_matrix<-matrix(NA, nrow = M, ncol = length(omega_values))
  Z2_rejection_rates<-numeric(length(omega_values))
  Z1_A_matrix<-matrix(NA, nrow = M, ncol = length(omega_values))
  Z1_rejection_rates<-numeric(length(omega_values))
  
  for (sim in 1:M) {

    z <- sqrt((df - 2) / df) * rt(T, df)
    
    x <- rep(NA, T)
    sigma2 <- rep(NA, T)
    H_0_sigma<-rep(NA,T)
    var_forcast<-rep(NA,T)
    es_forcast<-rep(NA,T)
    sigma2[1] <- omega / (1 - alpha - beta)
    x[1] <- sqrt(sigma2[1]) * z[1]
    H_0_sigma[1] <- sqrt(sigma2[1])
    var_forcast[1] <- VaR(alpha=0.975, df=df, Norm=FALSE) * H_0_sigma[1]
    es_forcast[1] <- ES(alpha=0.975, df=df, Norm=FALSE) * H_0_sigma[1]
    
    for (i in 2:T) {
      sigma2[i] <- omega + alpha * x[i - 1]^2 + beta * sigma2[i - 1]
      x[i] <- mu+phi * x[i - 1] + sqrt(sigma2[i]) * z[i]
      H_0_sigma[i] <- sqrt(sigma2[i])
      var_forcast[i] <- VaR(alpha=0.975,df=df, Norm=FALSE) * H_0_sigma[i]
      es_forcast[i] <- ES(alpha=0.975, df=df, Norm=FALSE) * H_0_sigma[i]
    }
    
    for (j in seq_along(omega_values)) {
      omega_A <- omega_values[j]
      
      x_A <- rep(NA, T)
      sigma2_A <- rep(NA, T)
      sigma2_A[1] <- omega_A / (1 - alpha - beta)
      x_A[1] <- sqrt(sigma2_A[1]) * z[1]
      
      for (i in 2:N) {
        sigma2_A[i] <- omega_A + alpha * x_A[i - 1]^2 + beta * sigma2_A[i - 1]
        x_A[i] <- mu+phi * x_A[i - 1] + sqrt(sigma2_A[i]) * z[i]
      }
     
      cc_p_value <- cc_backtest(r=x_A, q=var_forcast, e=es_forcast, s=H_0_sigma, alpha=0.975)$pvalue_onesided_general
      cc_p_values_matrix[sim, j] <- cc_p_value
      esr_p_value <- esr_backtest(r=x_A, e=-es_forcast, alpha=0.975)$pvalue_onesided_asymptotic
      esr_p_values_matrix[sim, j] <- esr_p_value
      
      Z2_A<-Z_2(x_A,var_forcast,es_forcast, quan = 0.05)$Z_2
      Z2_A_matrix[sim, j]<-Z2_A
      Z2<-Z_2(x,var_forcast,es_forcast, quan = 0.05)$Z_2
      Z1_A<-Z_1(x_A,var_forcast,es_forcast, quan = 0.05)$Z_1
      Z1_A_matrix[sim, j]<-Z1_A
      Z1<-Z_1(x,var_forcast,es_forcast, quan = 0.05)$Z_1
      Z1 <- Z1[!is.na(Z1)]
    }
  }
  
  for (j in seq_along(omega_values)) {
    cc_rejection_rates[j] <- sum(cc_p_values_matrix[, j] < 0.05, na.rm = TRUE) / simulations
    esr_rejection_rates[j] <- sum(esr_p_values_matrix[, j] < 0.05, na.rm = TRUE) / simulations
    Z2_rejection_rates[j]<-rejection_rate(Z2, Z2_A_matrix[,j], 0.05, simulations)
    Z1_rejection_rates[j]<-rejection_rate(Z1, Z1_A_matrix[,j], 0.05, simulations)
  }
  
plot(omega_values, cc_rejection_rates, type = "o", pch = 16, col = "blue", 
       xlab = expression(omega), ylab = "Rejection Rate", 
       main = "Rejection Rate vs. Basisvariansen (omega)", ylim = c(0, 1))
lines(omega_values, esr_rejection_rates, type = "o", pch = 16, col = "red")
lines(omega_values, Z2_rejection_rates, type = "o", pch = 16, col = "purple")
lines(omega_values, Z1_rejection_rates, type = "o", pch = 16, col = "deeppink")
abline(v = fitted_values_train["omega"], col = "darkgray", lty = 2, lwd = 2)

# Tilføj en forklarende legend
legend("topright", legend = c("CC Backtest", "ESR Backtest",  "Z1 backtest", "Z2 backtest" ,"Omega Reference"), 
       col = c("blue", "red", "deeppink", "purple", "darkgray"), pch = c(16, 16, 16, 16, NA), lty = c(1, 1,1, 1,2))

}
k_omega(omega_values, 250)
```

Tæheder med forskellige værdier af alpha
```{r,warning=FALSE}
T<-250
M<-10000

cc0_s <- rep(NA,M)
cc1_s <- rep(NA,M)
cc2_s <- rep(NA,M)
cc3_s <- rep(NA,M)
cc4_s <- rep(NA,M)
esr0_s <- rep(NA,M)
esr1_s <- rep(NA,M)
esr2_s <- rep(NA,M)
esr3_s <- rep(NA,M)
esr4_s <- rep(NA,M)
Z10_s <- rep(NA,M)
Z11_s <- rep(NA,M)
Z12_s <- rep(NA,M)
Z13_s <- rep(NA,M)
Z14_s <- rep(NA,M)
Z20_s <- rep(NA,M)
Z21_s <- rep(NA,M)
Z22_s <- rep(NA,M)
Z23_s <- rep(NA,M)
Z24_s <- rep(NA,M)

  phi <- fitted_values_train["ar1"]
  omega <- fitted_values_train["omega"]
  alpha <- fitted_values_train["alpha1"]
  beta <- fitted_values_train["beta1"]
  alpha_beta_sum <- alpha + beta
  df <- fitted_values_train["shape"]
  mu <- fitted_values_train["mu"]

sigma_matrix<-matrix(NA, nrow = T, ncol = M)
x_matrix <- matrix(NA, nrow = T, ncol = M)
var_forcast_matrix<- matrix(NA, nrow = T, ncol = M)
es_forcast_matrix<- matrix(NA, nrow = T, ncol = M)
xA1_matrix<-matrix(NA,nrow=T,ncol=M)
xA2_matrix<-matrix(NA,nrow=T,ncol=M)
xA3_matrix<-matrix(NA,nrow=T,ncol=M)
xA4_matrix<-matrix(NA,nrow=T,ncol=M)
  
  for (sim in 1:M) {

    z <- rt(T, df = df) / sqrt(df / (df - 2)) 
    
    # Simulate original process
    x <- rep(NA, T)
    sigma2 <- rep(NA, T)
    H_0_sigma<-rep(NA,T)
    var_forcast<-rep(NA,T)
    es_forcast<-rep(NA,T)
    sigma2[1] <- omega / (1 - alpha - beta)
    x[1] <- sqrt(sigma2[1]) * z[1]
    H_0_sigma[1] <- sqrt(sigma2[1])
    var_forcast[1] <- VaR(alpha=0.975, df=df, Norm=FALSE) * H_0_sigma[1]
    es_forcast[1] <- ES(alpha=0.975, df=df, Norm=FALSE) * H_0_sigma[1]
    
    for (i in 2:T) {
      sigma2[i] <- omega + alpha * x[i - 1]^2 + beta * sigma2[i - 1]
      x[i] <- mu+phi * x[i - 1] + sqrt(sigma2[i]) * z[i]
      H_0_sigma[i] <- sqrt(sigma2[i])
      var_forcast[i] <- VaR(alpha=0.975, df=df, Norm=FALSE) * H_0_sigma[i]
      es_forcast[i] <- ES(alpha=0.975, df=df, Norm=FALSE) * H_0_sigma[i]
    }
    
  x_matrix[, sim] <- x
  sigma_matrix[,sim] <- H_0_sigma
  var_forcast_matrix[,sim]<-var_forcast
  es_forcast_matrix[,sim]<-es_forcast
    
      alpha_1 <- 0.03
      alpha_2 <- 0.1
      alpha_3 <- 0.2
      alpha_4 <- 0.25
      beta_1 <- alpha_beta_sum - alpha_1
      beta_2 <- alpha_beta_sum - alpha_2
      beta_3 <- alpha_beta_sum - alpha_3
      beta_4 <- alpha_beta_sum - alpha_4
      
      x_A1 <- rep(NA, T)
      x_A2 <- rep(NA, T)
      x_A3 <- rep(NA, T)
      x_A4 <- rep(NA, T)
      sigma2_A1 <- rep(NA, T)
      sigma2_A2 <- rep(NA, T)
      sigma2_A3 <- rep(NA, T)
      sigma2_A4 <- rep(NA, T)
      sigma2_A1[1] <- omega / (1 - alpha_1 - beta_1)
      sigma2_A2[1] <- omega / (1 - alpha_1 - beta_1)
      sigma2_A3[1] <- omega / (1 - alpha_3 - beta_3)
      sigma2_A4[1] <- omega / (1 - alpha_4 - beta_4)
      x_A1[1] <- sqrt(sigma2_A1[1]) * z[1]
      x_A2[1] <- sqrt(sigma2_A2[1]) * z[1]
      x_A3[1] <- sqrt(sigma2_A3[1]) * z[1]
      x_A4[1] <- sqrt(sigma2_A4[1]) * z[1]
      
      for (i in 2:T) {
        sigma2_A1[i] <- omega + alpha_1 * x_A1[i - 1]^2 + beta_1 * sigma2_A1[i - 1]
        x_A1[i] <- mu + phi * x_A1[i - 1] + sqrt(sigma2_A1[i]) * z[i]
        sigma2_A2[i] <- omega + alpha_2 * x_A2[i - 1]^2 + beta_2 * sigma2_A2[i - 1]
        x_A2[i] <- mu + phi * x_A2[i - 1] + sqrt(sigma2_A2[i]) * z[i]
        sigma2_A3[i] <- omega + alpha_3 * x_A3[i - 1]^2 + beta_3 * sigma2_A3[i - 1]
        x_A3[i] <- mu + phi * x_A3[i - 1] + sqrt(sigma2_A3[i]) * z[i]
        sigma2_A4[i] <- omega + alpha_4 * x_A4[i - 1]^2 + beta_4 * sigma2_A4[i - 1]
        x_A4[i] <- mu + phi * x_A4[i - 1] + sqrt(sigma2_A4[i]) * z[i]
      }
      xA1_matrix[, sim] <- x_A1
      xA2_matrix[, sim] <- x_A2
      xA3_matrix[, sim] <- x_A3
      xA4_matrix[, sim] <- x_A4
  }

 for (i in 1:M){
   cc0_s[[i]] <- cc_backtest(r=x_matrix[,i], q=var_forcast_matrix[,i], e=es_forcast_matrix[,i], s=sigma_matrix[,i], alpha=0.975)$min_t2
   cc1_s[[i]] <- cc_backtest(r=xA1_matrix[,i], q=var_forcast_matrix[,i], e=es_forcast_matrix[,i], s=sigma_matrix[,i], alpha=0.975)$min_t2
   cc2_s[[i]] <- cc_backtest(r=xA2_matrix[,i], q=var_forcast_matrix[,i], e=es_forcast_matrix[,i], s=sigma_matrix[,i], alpha=0.975)$min_t2
   cc3_s[[i]] <- cc_backtest(r=xA3_matrix[,i], q=var_forcast_matrix[,i], e=es_forcast_matrix[,i], s=sigma_matrix[,i], alpha=0.975)$min_t2
   cc4_s[[i]] <- cc_backtest(r=xA4_matrix[,i], q=var_forcast_matrix[,i], e=es_forcast_matrix[,i], s=sigma_matrix[,i], alpha=0.975)$min_t2

 }
for (i in 1:M){
  esr0_s[[i]] <- esr_backtest(r=x_matrix[,i], e=-es_forcast_matrix[,i], alpha=0.975)$t0
  esr1_s[[i]] <- esr_backtest(r=xA1_matrix[,i], e=-es_forcast_matrix[,i], alpha=0.975)$t0
  esr2_s[[i]] <- esr_backtest(r=xA2_matrix[,i], e=-es_forcast_matrix[,i], alpha=0.975)$t0
  esr3_s[[i]] <- esr_backtest(r=xA3_matrix[,i], e=-es_forcast_matrix[,i], alpha=0.975)$t0
 esr4_s[[i]] <- esr_backtest(r=xA4_matrix[,i], e=-es_forcast_matrix[,i], alpha=0.975)$t0
}

for (i in 1:M){
Z10_s[[i]]<-Z_1(x_matrix[,i],var_forcast_matrix[,i],es_forcast_matrix[,i])$Z_1
Z20_s[[i]]<-Z_2(x_matrix[,i],var_forcast_matrix[,i],es_forcast_matrix[,i])$Z_2
Z11_s[[i]]<-Z_1(xA1_matrix[,i],var_forcast_matrix[,i],es_forcast_matrix[,i])$Z_1
Z21_s[[i]]<-Z_2(xA1_matrix[,i],var_forcast_matrix[,i],es_forcast_matrix[,i])$Z_2
Z12_s[[i]]<-Z_1(xA2_matrix[,i],var_forcast_matrix[,i],es_forcast_matrix[,i])$Z_1
Z22_s[[i]]<-Z_2(xA2_matrix[,i],var_forcast_matrix[,i],es_forcast_matrix[,i])$Z_2
Z13_s[[i]]<-Z_1(xA3_matrix[,i],var_forcast_matrix[,i],es_forcast_matrix[,i])$Z_1
Z23_s[[i]]<-Z_2(xA3_matrix[,i],var_forcast_matrix[,i],es_forcast_matrix[,i])$Z_2
Z14_s[[i]]<-Z_1(xA4_matrix[,i],var_forcast_matrix[,i],es_forcast_matrix[,i])$Z_1
Z24_s[[i]]<-Z_2(xA4_matrix[,i],var_forcast_matrix[,i],es_forcast_matrix[,i])$Z_2
}
       
       
  data_cc_s <- data.frame(
  Value_cc_s = c(cc0_s, cc1_s,cc2_s,cc3_s,cc4_s),
  Type_cc_s = rep(c("CC H0", "CC 0.03 H1", "CC 0.1 H1" , "CC 0.2 H1",  "CC 0.25 H1"), c(length(cc0_s), length(cc1_s), length(cc2_s), length(cc3_s), length(cc4_s)))
)
data_esr_s <- data.frame(
  Value_esr_s = c(esr0_s, esr1_s, esr2_s, esr3_s, esr4_s),
  Type_esr_s = rep(c("Intercept ESR H0", "Intercept ESR 0.03 H1", "Intercept ESR 0.1 H1", "Intercept ESR 0.2 H1", "Intercept ESR 0.25 H1"), c(length(esr0_s), length(esr1_s), length(esr2_s), length(esr3_s), length(esr4_s)))
)
data_Z1_s <- data.frame(
  Value_Z1_s = c(Z10_s, Z11_s, Z12_s, Z13_s, Z14_s),
  Type_Z1_s = rep(c("Z1 H0", "Z1 0.03 H1", "Z1 0.1 H1", "Z1 0.2 H1", "Z1 0.25 H1"), c(length(Z10_s), length(Z11_s), length(Z12_s), length(Z13_s), length(Z14_s)))
)
data_Z2_s <- data.frame(
  Value_Z2_s = c(Z20_s, Z21_s, Z22_s, Z23_s, Z24_s),
  Type_Z2_s = rep(c("Z2 H0", "Z2 0.03 H1", "Z2 0.1 H1", "Z2 0.2 H1", "Z2 0.25 H1"), c(length(Z20_s), length(Z21_s), length(Z22_s), length(Z23_s), length(Z24_s)))
)

ggplot(data_cc_s, aes(x = Value_cc_s, color = Type_cc_s)) +
  geom_density(size = 0.8, adjust = 2, bw=0.5) +
  labs(title = "Tætheden for CC backtest under H0 og H1 ved forskellig alpha",
       x = "Value", y = "Density") +
  scale_color_manual(values = c("CC H0"="black", "CC 0.03 H1"="blue", "CC 0.1 H1"="red", "CC 0.2 H1"="green4",  "CC 0.25 H1"="purple"))+
  xlim(-8, 1)+
  theme_minimal() +theme(legend.title = element_blank(),legend.position = "bottom")

ggplot(data_esr_s, aes(x = Value_esr_s, color = Type_esr_s)) +
  geom_density(size = 0.8, adjust = 2) +
  labs(title = "Tætheden for Intercept ESR backtest under H0 og H1 ved forskellige alpha",
       x = "Value", y = "Density") +
  scale_color_manual(values = c("Intercept ESR H0"="black", "Intercept ESR 0.03 H1"="blue", "Intercept ESR 0.1 H1"="red", "Intercept ESR 0.2 H1"="green4", "Intercept ESR 0.25 H1"="purple")) +
  theme_minimal() +theme(legend.title = element_blank(),legend.position = "bottom", legend.text = element_text(size = 6))

ggplot(data_Z1_s, aes(x = Value_Z1_s, color = Type_Z1_s)) +
  geom_density(size = 0.8, adjust = 2) +
  labs(title = "Tætheden for Z1 backtest under H0 og H1 ved forskellig alpha",
       x = "Value", y = "Density") +
  scale_color_manual(values = c("Z1 H0"="black", "Z1 0.03 H1"="blue", "Z1 0.1 H1"="red", "Z1 0.2 H1"="green4",  "Z1 0.25 H1"="purple"))+
  #xlim(-0.5, 0.75)+
  theme_minimal() +theme(legend.title = element_blank(),legend.position = "bottom")

ggplot(data_Z2_s, aes(x = Value_Z2_s, color = Type_Z2_s)) +
  geom_density(size = 0.8, adjust = 2) +
  labs(title = "Tætheden for Z2 backtest under H0 og H1 ved forskellig alpha",
       x = "Value", y = "Density") +
  scale_color_manual(values = c("Z2 H0"="black", "Z2 0.03 H1"="blue", "Z2 0.1 H1"="red", "Z2 0.2 H1"="green4",  "Z2 0.25 H1"="purple"))+
  #xlim(-2, 3)+
  theme_minimal() +theme(legend.title = element_blank(),legend.position = "bottom")
```

Tæheder med forskellige værdier af phi
```{r,warning=FALSE}
T<-250
M<-10000

cc0_s_phi <- rep(NA,M)
cc1_s_phi <- rep(NA,M)
cc2_s_phi <- rep(NA,M)
cc3_s_phi <- rep(NA,M)
cc4_s_phi <- rep(NA,M)
esr0_s_phi <- rep(NA,M)
esr1_s_phi <- rep(NA,M)
esr2_s_phi <- rep(NA,M)
esr3_s_phi <- rep(NA,M)
esr4_s_phi <- rep(NA,M)
Z10_s_phi <- rep(NA,M)
Z11_s_phi <- rep(NA,M)
Z12_s_phi <- rep(NA,M)
Z13_s_phi <- rep(NA,M)
Z14_s_phi <- rep(NA,M)
Z20_s_phi <- rep(NA,M)
Z21_s_phi <- rep(NA,M)
Z22_s_phi <- rep(NA,M)
Z23_s_phi <- rep(NA,M)
Z24_s_phi <- rep(NA,M)

  phi <- fitted_values_train["ar1"]
  omega <- fitted_values_train["omega"]
  alpha <- fitted_values_train["alpha1"]
  beta <- fitted_values_train["beta1"]
  df <- fitted_values_train["shape"]
  mu <- fitted_values_train["mu"]

sigma_matrix<-matrix(NA, nrow = T, ncol = M)
x_matrix <- matrix(NA, nrow = T, ncol = M)
var_forcast_matrix<- matrix(NA, nrow = T, ncol = M)
es_forcast_matrix<- matrix(NA, nrow = T, ncol = M)
xA1_matrix<-matrix(NA,nrow=T,ncol=M)
xA2_matrix<-matrix(NA,nrow=T,ncol=M)
xA3_matrix<-matrix(NA,nrow=T,ncol=M)
xA4_matrix<-matrix(NA,nrow=T,ncol=M)
  
  for (sim in 1:M) {

    z <-rt(T, df = df) / sqrt(df / (df - 2)) 
    
    x <- rep(NA, T)
    sigma2 <- rep(NA, T)
    H_0_sigma<-rep(NA,T)
    var_forcast<-rep(NA,T)
    es_forcast<-rep(NA,T)
    sigma2[1] <- omega / (1 - alpha - beta)
    x[1] <- sqrt(sigma2[1]) * z[1]
    H_0_sigma[1] <- sqrt(sigma2[1])
    var_forcast[1] <- VaR(alpha=0.975, df=df, Norm=FALSE) * H_0_sigma[1]
    es_forcast[1] <- ES(alpha=0.975, df=df, Norm=FALSE) * H_0_sigma[1]
    
    for (i in 2:T) {
      sigma2[i] <- omega + alpha * x[i - 1]^2 + beta * sigma2[i - 1]
      x[i] <- mu+phi * x[i - 1] + sqrt(sigma2[i]) * z[i]
      H_0_sigma[i] <- sqrt(sigma2[i])
      var_forcast[i] <- VaR(alpha=0.975, df=df, Norm=FALSE) * H_0_sigma[i]
      es_forcast[i] <- ES(alpha=0.975, df=df, Norm=FALSE) * H_0_sigma[i]
    }
    
  x_matrix[, sim] <- x
  sigma_matrix[,sim] <- H_0_sigma
  var_forcast_matrix[,sim]<-var_forcast
  es_forcast_matrix[,sim]<-es_forcast
    
      phi_1 <- -0.3
      phi_2 <- 0
      phi_3 <- 0.3
      phi_4 <- 0.5
      
      x_A1 <- rep(NA, T)
      x_A2 <- rep(NA, T)
      x_A3 <- rep(NA, T)
      x_A4 <- rep(NA, T)
      sigma2_A <- rep(NA, T)
      sigma2_A[1] <- omega / (1 - alpha - beta)
      x_A1[1] <- sqrt(sigma2_A1[1]) * z[1]
      x_A2[1] <- sqrt(sigma2_A2[1]) * z[1]
      x_A3[1] <- sqrt(sigma2_A3[1]) * z[1]
      x_A4[1] <- sqrt(sigma2_A4[1]) * z[1]
      
      for (i in 2:T) {
        sigma2_A[i] <- omega + alpha * x_A1[i - 1]^2 + beta * sigma2_A1[i - 1]
        x_A1[i] <- mu + phi_1 * x_A1[i - 1] + sqrt(sigma2_A[i]) * z[i]
        x_A2[i] <- mu + phi_2 * x_A2[i - 1] + sqrt(sigma2_A[i]) * z[i]
        x_A3[i] <- mu + phi_3 * x_A3[i - 1] + sqrt(sigma2_A[i]) * z[i]
        x_A4[i] <- mu + phi_4 * x_A4[i - 1] + sqrt(sigma2_A[i]) * z[i]
      }
      xA1_matrix[, sim] <- x_A1
      xA2_matrix[, sim] <- x_A2
      xA3_matrix[, sim] <- x_A3
      xA4_matrix[, sim] <- x_A4
  }

 for (i in 1:M){
   cc0_s_phi[[i]] <- cc_backtest(r=x_matrix[,i], q=var_forcast_matrix[,i], e=es_forcast_matrix[,i], s=sigma_matrix[,i], alpha=0.975)$min_t2
      cc1_s_phi[[i]] <- cc_backtest(r=xA1_matrix[,i], q=var_forcast_matrix[,i], e=es_forcast_matrix[,i], s=sigma_matrix[,i], alpha=0.975)$min_t2
      cc2_s_phi[[i]] <- cc_backtest(r=xA2_matrix[,i], q=var_forcast_matrix[,i], e=es_forcast_matrix[,i], s=sigma_matrix[,i], alpha=0.975)$min_t2
      cc3_s_phi[[i]] <- cc_backtest(r=xA3_matrix[,i], q=var_forcast_matrix[,i], e=es_forcast_matrix[,i], s=sigma_matrix[,i], alpha=0.975)$min_t2
      cc4_s_phi[[i]] <- cc_backtest(r=xA4_matrix[,i], q=var_forcast_matrix[,i], e=es_forcast_matrix[,i], s=sigma_matrix[,i], alpha=0.975)$min_t2

 }
for (i in 1:M){
  esr0_s_phi[[i]] <- esr_backtest(r=x_matrix[,i], e=-es_forcast_matrix[,i], alpha=0.975)$t0
  esr1_s_phi[[i]] <- esr_backtest(r=xA1_matrix[,i], e=-es_forcast_matrix[,i], alpha=0.975)$t0
  esr2_s_phi[[i]] <- esr_backtest(r=xA2_matrix[,i], e=-es_forcast_matrix[,i], alpha=0.975)$t0
  esr3_s_phi[[i]] <- esr_backtest(r=xA3_matrix[,i], e=-es_forcast_matrix[,i], alpha=0.975)$t0
 esr4_s_phi[[i]] <- esr_backtest(r=xA4_matrix[,i], e=-es_forcast_matrix[,i], alpha=0.975)$t0
 }

for (i in 1:M){
Z10_s_phi[[i]]<-Z_1(x_matrix[,i],var_forcast_matrix[,i],es_forcast_matrix[,i])$Z_1
Z20_s_phi[[i]]<-Z_2(x_matrix[,i],var_forcast_matrix[,i],es_forcast_matrix[,i])$Z_2
Z11_s_phi[[i]]<-Z_1(xA1_matrix[,i],var_forcast_matrix[,i],es_forcast_matrix[,i])$Z_1
Z21_s_phi[[i]]<-Z_2(xA1_matrix[,i],var_forcast_matrix[,i],es_forcast_matrix[,i])$Z_2
Z12_s_phi[[i]]<-Z_1(xA2_matrix[,i],var_forcast_matrix[,i],es_forcast_matrix[,i])$Z_1
Z22_s_phi[[i]]<-Z_2(xA2_matrix[,i],var_forcast_matrix[,i],es_forcast_matrix[,i])$Z_2
Z13_s_phi[[i]]<-Z_1(xA3_matrix[,i],var_forcast_matrix[,i],es_forcast_matrix[,i])$Z_1
Z23_s_phi[[i]]<-Z_2(xA3_matrix[,i],var_forcast_matrix[,i],es_forcast_matrix[,i])$Z_2
Z14_s_phi[[i]]<-Z_1(xA4_matrix[,i],var_forcast_matrix[,i],es_forcast_matrix[,i])$Z_1
Z24_s_phi[[i]]<-Z_2(xA4_matrix[,i],var_forcast_matrix[,i],es_forcast_matrix[,i])$Z_2
}      
       
data_cc_s_phi <- data.frame(
  Value_cc_s_phi = c(cc0_s_phi, cc1_s_phi,cc2_s_phi,cc3_s_phi,cc4_s_phi),
  Type_cc_s_phi = rep(c("CC H0", "CC -0.3 H1", "CC 0 H1" , "CC 0.3 H1",  "CC 0.5 H1"), c(length(cc0_s_phi), length(cc1_s_phi), length(cc2_s_phi), length(cc3_s_phi), length(cc4_s_phi)))
)
data_esr_s_phi <- data.frame(
  Value_esr_s_phi = c(esr0_s_phi, esr1_s_phi, esr2_s_phi, esr3_s_phi, esr4_s_phi),
  Type_esr_s_phi = rep(c("Intercept ESR H0", "Intercept ESR -0.3 H1", "Intercept ESR 0 H1", "Intercept ESR 0.3 H1", "Intercept ESR 0.5 H1"), c(length(esr0_s_phi), length(esr1_s_phi), length(esr2_s_phi), length(esr3_s_phi), length(esr4_s_phi)))
)

data_Z1_s_phi <- data.frame(
  Value_Z1_s_phi = c(Z10_s_phi, Z11_s_phi, Z12_s_phi, Z13_s_phi, Z14_s_phi),
  Type_Z1_s_phi = rep(c("Z1 H0", "Z1 -0.3 H1", "Z1 0 H1", "Z1 0.3 H1", "Z1 0.5 H1"), c(length(Z10_s_phi), length(Z11_s_phi), length(Z12_s_phi), length(Z13_s_phi), length(Z14_s_phi)))
)
data_Z2_s_phi <- data.frame(
  Value_Z2_s_phi = c(Z20_s_phi, Z21_s_phi, Z22_s_phi, Z23_s_phi, Z24_s_phi),
  Type_Z2_s_phi = rep(c("Z2 H0", "Z2 -0.3 H1", "Z2 0 H1", "Z2 0.3 H1", "Z2 0.5 H1"), c(length(Z20_s_phi), length(Z21_s_phi), length(Z22_s_phi), length(Z23_s_phi), length(Z24_s_phi)))
)

ggplot(data_cc_s_phi, aes(x = Value_cc_s_phi, color = Type_cc_s_phi)) +
  geom_density(size = 0.8, bw=0.5) +
  labs(title = "Tætheden for CC backtest under H0 og H1 ved forskellig phi",
       x = "Value", y = "Density") +
  scale_color_manual(values = c("CC H0"="black", "CC -0.3 H1"="blue", "CC 0 H1"="red", "CC 0.3 H1"="green4",  "CC 0.5 H1"="purple"))+
  xlim(-8,2)+
  theme_minimal() +theme(legend.title = element_blank(),legend.position = "bottom")

ggplot(data_esr_s_phi, aes(x = Value_esr_s_phi, color = Type_esr_s_phi)) +
  geom_density(size = 0.8) +
  labs(title = "Tætheden for Intercept ESR backtest under H0 og H1 ved forskellige phi",
       x = "Value", y = "Density") +
  scale_color_manual(values = c("Intercept ESR H0"="black", "Intercept ESR -0.3 H1"="blue", "Intercept ESR 0 H1"="red", "Intercept ESR 0.3 H1"="green4", "Intercept ESR 0.5 H1"="purple")) +
  theme_minimal() +theme(legend.title = element_blank(),legend.position = "bottom")

ggplot(data_Z1_s_phi, aes(x = Value_Z1_s_phi, color = Type_Z1_s_phi)) +
  geom_density(size = 0.8, adjust = 2) +
  labs(title = "Tætheden for Z1 backtest under H0 og H1 ved forskellig phi",
       x = "Value", y = "Density") +
  scale_color_manual(values = c("Z1 H0"="black", "Z1 -0.3 H1"="blue", "Z1 0 H1"="red", "Z1 0.3 H1"="green4",  "Z1 0.5 H1"="purple"))+
  theme_minimal() +theme(legend.title = element_blank(),legend.position = "bottom")

ggplot(data_Z2_s_phi, aes(x = Value_Z2_s_phi, color = Type_Z2_s_phi)) +
  geom_density(size = 0.8, adjust = 2) +
  labs(title = "Tætheden for Z2 backtest under H0 og H1 ved forskellig phi",
       x = "Value", y = "Density") +
  scale_color_manual(values = c("Z2 H0"="black", "Z2 -0.3 H1"="blue", "Z2 0 H1"="red", "Z2 0.3 H1"="green4",  "Z2 0.5 H1"="purple"))+
  theme_minimal() +theme(legend.title = element_blank(),legend.position = "bottom")
```

Tæheder med forskellige værdier af omega
```{r,warning=FALSE}
T<-250
M<-10000

cc0_omega <- rep(NA,M)
cc1_omega <- rep(NA,M)
cc2_omega <- rep(NA,M)
cc3_omega <- rep(NA,M)
cc4_omega <- rep(NA,M)
esr0_omega <- rep(NA,M)
esr1_omega <- rep(NA,M)
esr2_omega <- rep(NA,M)
esr3_omega <- rep(NA,M)
esr4_omega <- rep(NA,M)
Z10_omega <- rep(NA,M)
Z11_omega <- rep(NA,M)
Z12_omega <- rep(NA,M)
Z13_omega <- rep(NA,M)
Z14_omega <- rep(NA,M)
Z20_omega <- rep(NA,M)
Z21_omega <- rep(NA,M)
Z22_omega <- rep(NA,M)
Z23_omega <- rep(NA,M)
Z24_omega <- rep(NA,M)

  phi <- fitted_values_train["ar1"]
  omega <- fitted_values_train["omega"]
  alpha <- fitted_values_train["alpha1"]
  beta <- fitted_values_train["beta1"]
  alpha_beta_sum <- alpha + beta
  df <- fitted_values_train["shape"]
  mu <- fitted_values_train["mu"]
  
sigma_matrix<-matrix(NA, nrow = T, ncol = M)
x_matrix <- matrix(NA, nrow = T, ncol = M)
var_forcast_matrix<- matrix(NA, nrow = T, ncol = M)
es_forcast_matrix<- matrix(NA, nrow = T, ncol = M)
xA1_matrix<-matrix(NA,nrow=T,ncol=M)
xA2_matrix<-matrix(NA,nrow=T,ncol=M)
xA3_matrix<-matrix(NA,nrow=T,ncol=M)
xA4_matrix<-matrix(NA,nrow=T,ncol=M)
  
  for (sim in 1:M) {
    
    z <- rt(T, df = df) / sqrt(df / (df - 2))
  
    x <- rep(NA, T)
    sigma2 <- rep(NA, T)
    H_0_sigma<-rep(NA,T)
    var_forcast<-rep(NA,T)
    es_forcast<-rep(NA,T)
    sigma2[1] <- omega / (1 - alpha - beta)
    x[1] <- sqrt(sigma2[1]) * z[1]
    H_0_sigma[1] <- sqrt(sigma2[1])
    var_forcast[1] <- VaR(alpha=0.975, df=df, Norm=FALSE) * H_0_sigma[1]
    es_forcast[1] <- ES(alpha=0.975, df=df, Norm=FALSE) * H_0_sigma[1]
    
    for (i in 2:T) {
      sigma2[i] <- omega + alpha * x[i - 1]^2 + beta * sigma2[i - 1]
      x[i] <- mu+phi * x[i - 1] + sqrt(sigma2[i]) * z[i]
      H_0_sigma[i] <- sqrt(sigma2[i])
      var_forcast[i] <- VaR(alpha=0.975, df=df, Norm=FALSE) * H_0_sigma[i]
      es_forcast[i] <- ES(alpha=0.975, df=df, Norm=FALSE) * H_0_sigma[i]
    }
    
  x_matrix[, sim] <- x
  sigma_matrix[,sim] <- H_0_sigma
  var_forcast_matrix[,sim]<-var_forcast
  es_forcast_matrix[,sim]<-es_forcast
    
      omega_1 <- 0.01
      omega_2 <- 0.03 
      omega_3 <- 0.06
      omega_4 <- 0.08
      
      x_A1 <- rep(NA, T)
      x_A2 <- rep(NA, T)
      x_A3 <- rep(NA, T)
      x_A4 <- rep(NA, T)
      sigma2_A1 <- rep(NA, T)
      sigma2_A2 <- rep(NA, T)
      sigma2_A3 <- rep(NA, T)
      sigma2_A4 <- rep(NA, T)
      sigma2_A1[1] <- omega_1 / (1 - alpha - beta)
      sigma2_A2[1] <- omega_2 / (1 - alpha - beta)
      sigma2_A3[1] <- omega_3 / (1 - alpha - beta)
      sigma2_A4[1] <- omega_4 / (1 - alpha - beta)
      x_A1[1] <- sqrt(sigma2_A1[1]) * z[1]
      x_A2[1] <- sqrt(sigma2_A2[1]) * z[1]
      x_A3[1] <- sqrt(sigma2_A3[1]) * z[1]
      x_A4[1] <- sqrt(sigma2_A4[1]) * z[1]
      
      for (i in 2:T) {
        sigma2_A1[i] <- omega_1 + alpha * x_A1[i - 1]^2 + beta * sigma2_A1[i - 1]
        x_A1[i] <- mu + phi * x_A1[i - 1] + sqrt(sigma2_A1[i]) * z[i]
        
        sigma2_A2[i] <- omega_2 + alpha * x_A2[i - 1]^2 + beta * sigma2_A2[i - 1]
        x_A2[i] <- mu + phi * x_A2[i - 1] + sqrt(sigma2_A2[i]) * z[i]
        
        sigma2_A3[i] <- omega_3 + alpha * x_A3[i - 1]^2 + beta * sigma2_A3[i - 1]
        x_A3[i] <- mu + phi * x_A3[i - 1] + sqrt(sigma2_A3[i]) * z[i]
        
        sigma2_A4[i] <- omega_4 + alpha * x_A4[i - 1]^2 + beta * sigma2_A4[i - 1]
        x_A4[i] <- mu + phi * x_A4[i - 1] + sqrt(sigma2_A4[i]) * z[i]
      }
      xA1_matrix[, sim] <- x_A1
      xA2_matrix[, sim] <- x_A2
      xA3_matrix[, sim] <- x_A3
      xA4_matrix[, sim] <- x_A4
  }

 for (i in 1:M){
   cc0_omega[[i]] <- cc_backtest(r=x_matrix[,i], q=var_forcast_matrix[,i], e=es_forcast_matrix[,i], s=sigma_matrix[,i], alpha=0.975)$min_t2
      cc1_omega[[i]] <- cc_backtest(r=xA1_matrix[,i], q=var_forcast_matrix[,i], e=es_forcast_matrix[,i], s=sigma_matrix[,i], alpha=0.975)$min_t2
      cc2_omega[[i]] <- cc_backtest(r=xA2_matrix[,i], q=var_forcast_matrix[,i], e=es_forcast_matrix[,i], s=sigma_matrix[,i], alpha=0.975)$min_t2
      cc3_omega[[i]] <- cc_backtest(r=xA3_matrix[,i], q=var_forcast_matrix[,i], e=es_forcast_matrix[,i], s=sigma_matrix[,i], alpha=0.975)$min_t2
      cc4_omega[[i]] <- cc_backtest(r=xA4_matrix[,i], q=var_forcast_matrix[,i], e=es_forcast_matrix[,i], s=sigma_matrix[,i], alpha=0.975)$min_t2

 }

for (i in 1:M){
  esr0_omega[[i]] <- esr_backtest(r=x_matrix[,i], e=-es_forcast_matrix[,i], alpha=0.975)$t0
  esr1_omega[[i]] <- esr_backtest(r=xA1_matrix[,i], e=-es_forcast_matrix[,i], alpha=0.975)$t0
  esr2_omega[[i]] <- esr_backtest(r=xA2_matrix[,i], e=-es_forcast_matrix[,i], alpha=0.975)$t0
  esr3_omega[[i]] <- esr_backtest(r=xA3_matrix[,i], e=-es_forcast_matrix[,i], alpha=0.975)$t0
 esr4_omega[[i]] <- esr_backtest(r=xA4_matrix[,i], e=-es_forcast_matrix[,i], alpha=0.975)$t0
}

for (i in 1:M){
Z10_omega[[i]]<-Z_1(x_matrix[,i],var_forcast_matrix[,i],es_forcast_matrix[,i])$Z_1
Z20_omega[[i]]<-Z_2(x_matrix[,i],var_forcast_matrix[,i],es_forcast_matrix[,i])$Z_2
Z11_omega[[i]]<-Z_1(xA1_matrix[,i],var_forcast_matrix[,i],es_forcast_matrix[,i])$Z_1
Z21_omega[[i]]<-Z_2(xA1_matrix[,i],var_forcast_matrix[,i],es_forcast_matrix[,i])$Z_2
Z12_omega[[i]]<-Z_1(xA2_matrix[,i],var_forcast_matrix[,i],es_forcast_matrix[,i])$Z_1
Z22_omega[[i]]<-Z_2(xA2_matrix[,i],var_forcast_matrix[,i],es_forcast_matrix[,i])$Z_2
Z13_omega[[i]]<-Z_1(xA3_matrix[,i],var_forcast_matrix[,i],es_forcast_matrix[,i])$Z_1
Z23_omega[[i]]<-Z_2(xA3_matrix[,i],var_forcast_matrix[,i],es_forcast_matrix[,i])$Z_2
Z14_omega[[i]]<-Z_1(xA4_matrix[,i],var_forcast_matrix[,i],es_forcast_matrix[,i])$Z_1
Z24_omega[[i]]<-Z_2(xA4_matrix[,i],var_forcast_matrix[,i],es_forcast_matrix[,i])$Z_2
}   
       
data_cc_omega <- data.frame(
  Value_cc_omega = c(cc0_omega, cc1_omega,cc2_omega,cc3_omega,cc4_omega),
  Type_cc_omega = rep(c("CC H0", "CC 0.01 H1", "CC 0.03 H1" , "CC 0.06 H1",  "CC 0.08 H1"), c(length(cc0_omega), length(cc1_omega), length(cc2_omega), length(cc3_omega), length(cc4_omega)))
)
data_esr_omega <- data.frame(
  Value_esr_omega = c(esr0_omega, esr1_omega, esr2_omega, esr3_omega, esr4_omega),
  Type_esr_omega = rep(c("Intercept ESR H0", "Intercept ESR 0.01 H1", "Intercept ESR 0.03 H1", "Intercept ESR 0.06 H1", "Intercept ESR 0.08 H1"), c(length(esr0_omega), length(esr1_omega), length(esr2_omega), length(esr3_omega), length(esr4_omega)))
)
data_Z1_omega <- data.frame(
  Value_Z1_omega = c(Z10_omega, Z11_omega, Z12_omega, Z13_omega, Z14_omega),
  Type_Z1_omega = rep(c("Z1 H0", "Z1 0.01 H1", "Z1 0.03 H1", "Z1 0.06 H1", "Z1 0.08 H1"), c(length(Z10_omega), length(Z11_omega), length(Z12_omega), length(Z13_omega), length(Z14_omega)))
)
data_Z2_omega <- data.frame(
  Value_Z2_omega = c(Z20_omega, Z21_omega, Z22_omega, Z23_omega, Z24_omega),
  Type_Z2_omega = rep(c("Z2 H0", "Z2 0.01 H1", "Z2 0.03 H1", "Z2 0.06 H1", "Z2 0.08 H1"), c(length(Z20_omega), length(Z21_omega), length(Z22_omega), length(Z23_omega), length(Z24_omega)))
)

ggplot(data_cc_omega, aes(x = Value_cc_omega, color = Type_cc_omega)) +
  geom_density(size = 0.8) +
  labs(title = "Tætheden for CC 3 backtest under H0 og H1 ved forskellig omega",
       x = "Value", y = "Density") +
  scale_color_manual(values = c("CC H0"="black", "CC 0.01 H1"="blue", "CC 0.03 H1"="red", "CC 0.06 H1"="green4",  "CC 0.08 H1"="purple"))+
  #xlim(-10,1)+
  theme_minimal() +theme(legend.title = element_blank(),legend.position = "bottom")

ggplot(data_esr_omega, aes(x = Value_esr_omega, color = Type_esr_omega)) +
  geom_density(size = 0.8) +
  labs(title = "Tætheden for Intercept ESR backtest under H0 og H1 ved forskellige omega",
       x = "Value", y = "Density") +
  scale_color_manual(values = c("Intercept ESR H0"="black", "Intercept ESR 0.01 H1"="blue", "Intercept ESR 0.03 H1"="red", "Intercept ESR 0.06 H1"="green4", "Intercept ESR 0.08 H1"="purple")) +
  xlim(-10,3)+
  theme_minimal() +theme(legend.title = element_blank(),legend.position = "bottom")

ggplot(data_Z1_omega, aes(x = Value_Z1_omega, color = Type_Z1_omega)) +
  geom_density(size = 0.8, adjust = 2) +
  labs(title = "Tætheden for Z1 backtest under H0 og H1 ved forskellig omega",
       x = "Value", y = "Density") +
  scale_color_manual(values = c("Z1 H0"="black", "Z1 0.01 H1"="blue", "Z1 0.03 H1"="red", "Z1 0.06 H1"="green4",  "Z1 0.08 H1"="purple"))+
  theme_minimal() +theme(legend.title = element_blank(),legend.position = "bottom")

ggplot(data_Z2_omega, aes(x = Value_Z2_omega, color = Type_Z2_omega)) +
  geom_density(size = 0.8, adjust = 2) +
  labs(title = "Tætheden for Z2 backtest under H0 og H1 ved forskellig omega",
       x = "Value", y = "Density") +
  scale_color_manual(values = c("Z2 H0"="black", "Z2 0.01 H1"="blue", "Z2 0.03 H1"="red", "Z2 0.06 H1"="green4",  "Z2 0.08 H1"="purple"))+
  xlim(-1.5,4)+
  theme_minimal() +theme(legend.title = element_blank(),legend.position = "bottom")
  
```
