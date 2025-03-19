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

Indlæsning af data
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
```

Lukkepriser og afkast
```{r, warning=FALSE}
x_breaks <- seq(min(SP500_Data$Date), max(SP500_Data$Date), by = "6 months")

par(mar = c(6, 4, 4, 4) + 0.1)

plot(SP500_Data$Date, SP500_Data$Price, type = "l", col = "red", 
     xlab = "", ylab = "Closing Price (USD)", lwd = 2, xaxt = "n")

axis(1, at = x_breaks, labels = format(x_breaks, "%Y-%m"), las = 2, cex.axis = 0.7)

mtext("Date", side = 1, line = 4, cex = 1)

par(new = TRUE)
plot(SP500_Data$Date, SP500_Data$Change * 100, type = "l", col = "blue", 
     axes = FALSE, xlab = "", ylab = "", lwd = 1)

axis(4)
mtext("Daily Returns (%)", side = 4, line = 3)

title(main = "SP500 Closing Price and Daily Returns")

```

Deler data op i træningsdata og testdata
```{r, warning=FALSE}
test_days <- 250

#Træningsdata
SP500_Train <- SP500_Data[1:(nrow(SP500_Data) - test_days), ]
SP500_Train_returns <- SP500_Train$Change*100

#Testdata
SP500_Test <- SP500_Data[(nrow(SP500_Data) - test_days + 1):nrow(SP500_Data), ]
SP500_Test_returns <- SP500_Test$Change*100
```

AR-GARCH - Train-data
```{r, warning=FALSE}
ARGARCH <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(1, 0)),
  distribution.model = "std"
)

SP500_ARGARCH_fit <- ugarchfit(spec = ARGARCH, data = SP500_Train_returns, solver="hybrid")
coef(SP500_ARGARCH_fit)
```

Test af signifikans
```{r, warning=FALSE}
SP500_ARGARCH_fit@fit$matcoef
```

Test for autokorrelation
```
fit_train<-ugarchfit(spec = ARGARCH, data = SP500_Train$Change*100)
std_resid_train <- residuals(fit_train, standardize = TRUE)
Box.test(std_resid_train, lag = 21, type = "Ljung-Box")
```

Train-data plot
```{r, warning=FALSE}
fitted_values_train <- coef(SP500_ARGARCH_fit)
omega <- fitted_values_train["omega"]
alpha1 <- fitted_values_train["alpha1"]
beta1 <- fitted_values_train["beta1"]
mu <- fitted_values_train["mu"]
df <- fitted_values_train["shape"]
phi <- fitted_values_train["ar1"]

sigma2 <- numeric(length(SP500_Train_returns))
sigma2[1] <- omega/(1 - alpha1 - beta1)

R_t <- SP500_Train_returns

for (i in 2:length(R_t)) {
  sigma2[i] <- omega + alpha1 * (R_t[i - 1])^2 + beta1 * sigma2[i - 1]
}

VaR_values_train <- VaR(alpha = 0.975, df = df) * sqrt(sigma2)
ES_values_train <- ES(alpha = 0.975, df = df) * sqrt(sigma2)

plot_data <- data.frame(
  Time = 1:length(R_t),
  Returns = R_t,
  VaR = -VaR_values_train,
  ES = -ES_values_train
)
plot_data$Exceedance <- ifelse(plot_data$Returns < plot_data$VaR, "Exceedance", "No Exceedance")
ggplot(plot_data, aes(x = Time)) +
  geom_line(aes(y = Returns, color = "Returns")) +
  geom_line(aes(y = VaR, color = "VaR (97.5%)")) +
  geom_line(aes(y = ES, color = "ES (97.5%)")) +
  geom_point(data = plot_data[plot_data$Exceedance == "Exceedance",], aes(y = Returns), color = "black", size = 1)+
  labs(
    title = "Returns, VaR, and ES from AR(1)-GARCH(1,1)",
    x = "Time in days",
    y = "Returns in %",
    color = "Legend"
  ) +
  scale_color_manual(values = c("Returns" = "blue", "VaR (97.5%)" = "red", "ES (97.5%)" = "green")) +
  theme_minimal()
```

Historgram og AR-GARCH fordeling
```{r, warning=FALSE}
n_sim <- length(SP500_Train_returns)
simulated_returns <- ugarchsim(SP500_ARGARCH_fit, n.sim = n_sim)@simulation$seriesSim

df_hist <- data.frame(returns = SP500_Train_returns, type = "Observed")
df_sim <- data.frame(returns = simulated_returns, type = "ARGARCH Simulated")
df_combined <- rbind(df_hist, df_sim)

ggplot(df_combined, aes(x = returns)) +
  geom_histogram(data = df_hist, aes(y = ..density..), bins = 50, fill = "blue", color = "black", alpha = 0.5) +
  geom_density(data = df_sim, aes(y = ..density..), color = "red", size = 1.2) +
  labs(title = "Histogram af SP500 afkast med ARGARCH-simuleret tæthed",
       x = "Afkast", y = "Frekvens") +
  theme_minimal() +
  scale_fill_manual(values = c("Observed" = "blue", "ARGARCH Simulated" = "red"))

```

Laver en simulation under AR(1)GARCh(1,1) for H0 parametre (traindata)
```{r}
T <- 250

phi<-fitted_values_train["ar1"]
omega <- fitted_values_train["omega"]
alpha <- fitted_values_train["alpha1"]
beta <- fitted_values_train["beta1"]
mu <- fitted_values_train["mu"]
df <- fitted_values_train["shape"]

z <- rt(N,df)*sqrt((df - 2) / df)

x <- rep(NA,T)
sigma2_sim <- rep(NA,T)

sigma2_sim[1] <- omega / (1 - alpha - beta)
x[1] <- sqrt(sigma2_sim[1]) * z[1]

for (i in 2:T) {
  sigma2_sim[i] = omega + alpha * x[i-1]^2+ beta*sigma2_sim[i-1]
  x[i] = mu+phi*x[i-1]+sqrt(sigma2_sim[i]) * z[i]
}

```

Passer simulerede afkast med vores model?
Sim under AR-GARCH, VaR/ES fra traindata
```{r}
sim_ARGARCH_fit <- ugarchfit(spec = ARGARCH, data = x)
simulation_coefficients <- coef(sim_ARGARCH_fit)

omega <- fitted_values_train["omega"]
alpha1 <- fitted_values_train["alpha1"]
beta1 <- fitted_values_train["beta1"]
mu <- fitted_values_train["mu"]
df <- fitted_values_train["shape"]
phi <- fitted_values_train["ar1"]

sigma2_new <- numeric(length(x))
sigma2_new[1] <-  omega /(1 - alpha1 - beta1)

R_t <- x

for (i in 2:length(R_t)) {
  sigma2_new[i] <- omega + alpha1 * (R_t[i - 1])^2 + beta1 * sigma2_new[i - 1]
}

VaR_values_x <- VaR(alpha = 0.975, df = df) * sqrt(sigma2_new)
ES_values_x <- ES(alpha = 0.975, df = df) * sqrt(sigma2_new)

plot_data <- data.frame(
  Time = 1:length(R_t),
  Returns = R_t,
  VaR = -VaR_values_x,
  ES = -ES_values_x
)
plot_data$Exceedance <- ifelse(plot_data$Returns < plot_data$VaR, "Exceedance", "No Exceedance")

ggplot(plot_data, aes(x = Time)) +
  geom_line(aes(y = Returns, color = "Returns")) +
  geom_line(aes(y = VaR, color = "VaR (97.5%)")) +
  geom_line(aes(y = ES, color = "ES (97.5%)")) +
  geom_point(data = plot_data[plot_data$Exceedance == "Exceedance",], aes(y = Returns), color = "black", size = 1)+
  labs(
    title = "Returns, VaR, and ES from AR(1)-GARCH(1,1)",
    x = "Time in days",
    y = "Returns in %",
    color = "Legend"
  ) +
  scale_color_manual(values = c("Returns" = "blue", "VaR (97.5%)" = "red", "ES (97.5%)" = "green")) +
  theme_minimal()

```

Er parameterne for vores simulation tæt på de sande parametre?
```{r}
orig_coef <- abs(coef(SP500_ARGARCH_fit))
sim_coef <- abs(coef(sim_ARGARCH_fit))

param_names <- names(coef(SP500_ARGARCH_fit))

par(mfrow=c(3, 2), mar=c(4, 4, 2, 1))

for (i in 1:length(param_names)) {
  values <- c(orig_coef[i], sim_coef[i])
  bar_pos <- barplot(values, 
                     names.arg=c("Træning", "Simulation"),
                     col=c("blue", "red"),
                     main=param_names[i],
                     ylab="Værdi",
                     ylim=c(0, max(values) * 1.5))
  
  text(x=bar_pos, 
       y=values * 1.05,
       labels=round(values, 4), 
       pos=3)
}
par(mfrow=c(1,1))

```

```{r}
M <- 10000  
T <- 250

param_results <- matrix(NA, nrow = M, ncol = 6)
colnames(param_results) <- c("mu", "ar1", "omega", "alpha1", "beta1", "shape")

for (sim in 1:M) {
  scaling <- sqrt((df - 2) / df)
  z <- scaling * rt(T, df)

  sigma2_sim <- rep(NA, T)
  x <- rep(NA, T)

  sigma2_sim[1] <- omega / (1 - alpha - beta)
  x[1] <- sqrt(sigma2_sim[1]) * z[1]
 
  for (i in 2:T) {
    sigma2_sim[i] <- omega + alpha * x[i-1]^2 + beta * sigma2_sim[i-1]
    x[i] <- mu+phi * x[i-1] + sqrt(sigma2_sim[i]) * z[i]
    
    if (is.na(x[i]) || is.infinite(x[i])) {
      cat("Simulation", sim, "mislykkedes: NA eller Inf i x", "\n")
      next
    }
  }
  
  spec <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(1, 0)),
    distribution.model = "std"
  )
  
  fit <- tryCatch(
    ugarchfit(spec, data = x, solver = "hybrid"),
    error = function(e) {
      cat("Simulation", sim, "mislykkedes: ", e$message, "\n")
      return(NULL)
    }
  )
  
  if (!is.null(fit)) {
    sim_coefs <- coef(fit)
    
    param_names <- c("mu", "ar1", "omega", "alpha1", "beta1", "shape")
    sim_values <- sapply(param_names, function(p) {
        return(ifelse(p %in% names(sim_coefs), sim_coefs[p], NA))
    })
    
    param_results[sim, ] <- sim_values
  }
}

param_results <- param_results[complete.cases(param_results), ]

cat("Antal succesfulde simulationer:", nrow(param_results), "ud af", M, "\n")

new_params_sim_mean <- colMeans(param_results, na.rm = TRUE)

```

Histogrammer og CI for simulationer - ligger H0 parametre indenfor CI?
```{r}
plot_histogram_with_CI <- function(param_data, param_name) {
  param_data <- param_data[is.finite(param_data)]
  
  CI <- quantile(param_data, probs = c(0.025, 0.975), na.rm = TRUE)
  mean_val <- mean(param_data, na.rm = TRUE)
  
  p <- ggplot(data = data.frame(Value = param_data), aes(x = Value)) +
    geom_histogram(color = "black", fill = "blue", bins = 30) +
    geom_vline(xintercept = CI, color = "red", linetype = "dashed") +
    geom_vline(xintercept = mean_val, color = "green", linetype = "solid") +
    labs(title = paste("Histogram og 95% CI for", param_name),
         x = param_name, 
         y = "Frekvens")  

  x_right <- max(param_data, na.rm = TRUE) * 0.8 
  y_max <- max(hist(param_data, plot = FALSE)$counts) * 1.01
  
  p <- p + 
    annotate("label", x = x_right, y = y_max, 
             label = paste("Mean:", round(mean_val, 4)), 
             hjust = 1, vjust = 1, color = "green", fill = "white", size = 3, label.size = 0.3) +
    annotate("label", x = x_right, y = y_max * 0.9, 
             label = paste("2.5%:", round(CI[1], 4)), 
             hjust = 1, vjust = 1, color = "red", fill = "white", size = 3, label.size = 0.3) +
    annotate("label", x = x_right, y = y_max * 0.8, 
             label = paste("97.5%:", round(CI[2], 4)), 
             hjust = 1, vjust = 1, color = "red", fill = "white", size = 3, label.size = 0.3)

  return(p)
}

par_names <- colnames(param_results)
plot_list <- list()
for (i in 1:ncol(param_results)) {
  plot <- plot_histogram_with_CI(param_results[, i], par_names[i])
  if (!is.null(plot)) {
    plot_list[[i]] <- plot
  }
}

plot_list <- Filter(Negate(is.null), plot_list)
if (length(plot_list) > 0) {
  library(gridExtra)
  library(grid)
  num_plots <- length(plot_list)
  half <- ceiling(num_plots / 2)
  
  grid.newpage()
  grid.arrange(grobs = plot_list[1:half], ncol = 3)
  
  if (half < num_plots) {
    grid.newpage()
    grid.arrange(grobs = plot_list[(half+1):num_plots], ncol = 3)
  }
}

```
