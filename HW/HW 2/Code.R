library(quantmod)
library(tseries)
library(rugarch)

# ------------------------------------------------

# 1.1 Get four years of adjusted closing prices
# Define your start and end dates (you can adjust as needed)
start_date <- as.Date("2019-04-01")
end_date   <- as.Date("2023-04-01")

# Download data from Yahoo Finance
getSymbols("NKE", 
           src  = "yahoo", 
           from = start_date, 
           to   = end_date)

# Extract the Adjusted Close column
NKE_adj <- Ad(NKE)  # Ad() extracts the adjusted close prices

# 1.2 Calculate log returns
# We use dailyReturn() with type = "log" or manually compute diff(log(.))
NKE_returns <- dailyReturn(NKE_adj, type = "log")

# You might remove any NAs that appear (usually at the first data point)
NKE_returns <- na.omit(NKE_returns)

# 1.3 Perform the Augmented Dickey-Fuller test
adf_result <- adf.test(NKE_returns)

# Print ADF test results
print(adf_result)

# ------------------------------------------------

# 2.1. Split data into training (N-5) and test (5) sets

N <- length(NKE_returns)
train_returns <- NKE_returns[1:(N-5)]  # first N-5 for model estimation
test_returns  <- NKE_returns[(N-4):N]  # last 5 for out-of-sample checks (optional)

# 2.2 Run ARMA(p, q) + GARCH(1,1) over p,q in [0..4]
# We'll store the AIC and BIC in a data frame
results <- data.frame(
  p   = integer(),
  q   = integer(),
  AIC = numeric(),
  BIC = numeric(),
  stringsAsFactors = FALSE
)

# Loop over all p, q combinations
for (p in 0:4) {
  for (q in 0:4) {
    # Specify an sGARCH(1,1) model with ARMA(p, q) in the mean
    spec <- ugarchspec(
      variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
      mean.model     = list(armaOrder = c(p,q), include.mean = TRUE),
      distribution.model = "norm"  # Normal distribution assumption
    )
    
    # Fit the model to the training set
    fit <- ugarchfit(spec = spec, data = train_returns, solver = "hybrid")
    
    # Extract info criteria
    ic <- infocriteria(fit)
    # By default, infocriteria returns: c(Akaike, Bayes, Shibata, Hannan-Quinn)
    AIC_val <- ic[1]  # AIC
    BIC_val <- ic[2]  # BIC
    
    # Append to results data frame
    results <- rbind(results, data.frame(
      p   = p,
      q   = q,
      AIC = AIC_val,
      BIC = BIC_val
    ))
  }
}
# 2.3 Review the table of AIC/BIC and pick best models
print(results)

# Find the best model by AIC
best_AIC <- results[which.min(results$AIC), ]
cat("Best model by AIC:\n")
print(best_AIC)

# Find the best model by BIC
best_BIC <- results[which.min(results$BIC), ]
cat("Best model by BIC:\n")
print(best_BIC)

# ------------------------------------------------

# 3.1 Use the actual p,q you found in your previous loops. 
# For demonstration, let's assume:
pAIC <- best_AIC$p
qAIC <- best_AIC$q
pBIC <- best_BIC$p
qBIC <- best_BIC$q

# 3.2 Define and fit each model on the training data
# Best AIC model
specAIC <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model     = list(armaOrder = c(pAIC, qAIC), include.mean = TRUE),
  distribution.model = "norm"
)
fitAIC <- ugarchfit(spec = specAIC, data = train_returns, solver = "hybrid")

# Best BIC model
specBIC <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model     = list(armaOrder = c(pBIC, qBIC), include.mean = TRUE),
  distribution.model = "norm"
)
fitBIC <- ugarchfit(spec = specBIC, data = train_returns, solver = "hybrid")

# 3.3 Forecast for the last 5 days using ugarchforecast()
fcastAIC <- ugarchforecast(fitAIC, n.ahead = 5)
fcastBIC <- ugarchforecast(fitBIC, n.ahead = 5)

# Extract the forecasted returns (mean forecasts)
# seriesFor is a matrix with columns for each step
aic_pred <- as.numeric(fcastAIC@forecast$seriesFor)
bic_pred <- as.numeric(fcastBIC@forecast$seriesFor)

# Actual returns for the test period
test_actual <- as.numeric(test_returns)  # length should be 5

# 3.3 Define a "random walk" forecast
last_train_return <- tail(train_returns, 1)
rw_pred <- rep(last_train_return, 5)

# 3.4 Compute accuracy metrics: bias, MAD, MSE
calc_errors <- function(actual, forecast) {
  e <- forecast - actual
  bias <- mean(e)
  mad  <- mean(abs(e))
  mse  <- mean(e^2)
  return(c(bias = bias, MAD = mad, MSE = mse))
}

aic_err <- calc_errors(test_actual, aic_pred)
bic_err <- calc_errors(test_actual, bic_pred)
rw_err  <- calc_errors(test_actual, rw_pred)

cat("Accuracy metrics (Bias, MAD, MSE):\n\n")

cat("Best AIC Model:\n")
print(aic_err)
cat("\n")

cat("Best BIC Model:\n")
print(bic_err)
cat("\n")

cat("Random Walk Model:\n")
print(rw_err)
cat("\n")

