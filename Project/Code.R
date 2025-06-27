# --------------------------------------------------------------------------
# FRE 7831 Project 
# Author:Charlie Yi
# Date: Sys.Date()
# Security: NKE
# --------------------------------------------------------------------------

# --- 0. Setup ---

# Load necessary libraries
library(quantmod)
library(PerformanceAnalytics)
library(forecast)
library(tseries)
library(zoo)
library(dplyr)
library(lubridate)
library(ggplot2)

# Define parameters
ticker <- "NKE"
start_date <- "2005-03-01" # Start date for ~20 years before end date
end_date <- "2025-02-28"
lookback_tsm <- 252
lookahead_tsm <- 21
sma_pairs <- list(c(63, 15), c(126, 30), c(252, 63)) # Short/Long pairs
window_years <- 10
risk_free_rate <- 0.0 # Annual risk-free rate assumption

# Function to calculate annualized Sharpe Ratio from daily returns
calculate_sharpe <- function(daily_returns, rf = risk_free_rate) {
  # Ensure numeric type and handle NAs
  daily_returns <- na.omit(as.numeric(daily_returns))
  if (length(daily_returns) < 2 || sd(daily_returns) == 0) {
    return(NA) # Cannot calculate Sharpe if not enough data or no volatility
  }
  # Handle case where all returns are zero
  if (all(daily_returns == 0)) {
    return(0)
  }
  mean_daily_return <- mean(daily_returns)
  sd_daily_return <- sd(daily_returns)
  daily_rf <- (1 + rf)^(1/252) - 1 # Approximate daily risk-free rate
  sharpe_ratio <- (mean_daily_return - daily_rf) / sd_daily_return
  annualized_sharpe <- sharpe_ratio * sqrt(252)
  return(annualized_sharpe)
}

# Function to calculate strategy returns (handles cash periods with 0 return)
calculate_strategy_returns <- function(prices, signals) {
  # Ensure prices and signals are xts objects aligned by time
  if (!is.xts(prices)) prices <- xts(prices, order.by = index(signals))
  
  # Align signals with prices, taking next day's price for return calculation
  returns <- lag(signals, 1) * dailyReturn(prices)
  # Replace NA signals (e.g., initial period) with 0 (holding cash)
  returns[is.na(returns)] <- 0
  colnames(returns) <- "StrategyReturn"
  return(returns)
}


# --- 1. Data Retrieval (Implicit Step) ---

# Download data
cat("Fetching data for", ticker, "from", start_date, "to", end_date, "\n")
getSymbols(ticker, src = "yahoo", from = start_date, to = end_date)
prices <- Ad(get(ticker)) # Use adjusted closing prices
colnames(prices) <- ticker

# Calculate daily log returns (often preferred for financial time series)
daily_log_returns <- dailyReturn(prices, type = 'log')
daily_log_returns <- daily_log_returns[-1, ] # Remove first NA value
colnames(daily_log_returns) <- paste0(ticker, "_LogReturn")

# Calculate daily arithmetic returns (needed for strategy performance)
daily_arith_returns <- dailyReturn(prices, type = 'arithmetic')
daily_arith_returns <- daily_arith_returns[-1, ] # Remove first NA value
colnames(daily_arith_returns) <- paste0(ticker, "_ArithReturn")

cat("Data fetched. Rows:", nrow(prices), "\n")
cat("Date Range:", format(start(prices)), "to", format(end(prices)), "\n\n")


# --- 2. Statistics for the Sample --- [cite: 3]

cat("--- Step 2: Descriptive Statistics (based on Daily Log Returns) ---\n")

# Ensure returns are numeric vector for some functions
returns_vec <- as.numeric(daily_log_returns)
returns_vec <- returns_vec[!is.na(returns_vec)] # Remove NAs if any

stats_summary <- data.frame(
  Mean = mean(returns_vec) * 252, # Annualized Mean
  Sigma = sd(returns_vec) * sqrt(252), # Annualized Std Dev (Volatility)
  Skewness = skewness(returns_vec),
  Kurtosis = kurtosis(returns_vec) # Excess Kurtosis
)

print(stats_summary)

cat("\nACF Plot for Daily Log Returns:\n")
# Plot ACF
acf(returns_vec, main = paste("ACF of", ticker, "Daily Log Returns"), lag.max = 20)
# You can also get the numeric values: acf_values <- acf(returns_vec, plot = FALSE)

cat("\n--- End of Step 2 ---\n\n")


# --- 3. ARMA Estimation and Unit Root Test --- [cite: 4]

cat("--- Step 3: ARMA Estimation and Unit Root Test (based on Daily Log Returns) ---\n")

# Estimate ARMA coefficients using auto.arima
cat("Running auto.arima...\n")
# Using the time series object directly
arma_model <- auto.arima(daily_log_returns, stationary = TRUE, seasonal = FALSE, stepwise = TRUE, approximation = FALSE)
cat("Best ARMA Model found by auto.arima:\n")
print(summary(arma_model))
# Note: auto.arima often centers the data, check coefficients like 'intercept' or 'drift'

# Check for unit root using Augmented Dickey-Fuller (ADF) test
cat("\nPerforming ADF Test for Unit Root...\n")
adf_test_result <- adf.test(daily_log_returns, alternative = "stationary")
print(adf_test_result)
# Interpretation: A small p-value (e.g., < 0.05) suggests rejecting the null hypothesis (presence of a unit root), indicating stationarity.

# Check for drift: The ADF test can include drift, or check the ARMA model's intercept/drift term significance.
# If the 'drift' term in auto.arima is significant (check Pr(>|z|) or similar), it indicates drift.
# Alternatively, run adf.test with drift explicitly:
# adf_test_with_drift <- ur.df(daily_log_returns, type = "drift", lags = ...) # Using urca package is another option

cat("--- End of Step 3 ---\n\n")


# --- 4. Implement Strategies & Choose Best SMA (Full Sample) ---

cat("--- Step 4: Implementing Strategies and Choosing Best SMA (Full Sample) ---\n")

# --- 4a. TSM Strategy Signal Generation --- [cite: 5]
# Signal: 1 if past 'lookback' return > 0, else 0 (momentum)
# We apply the signal for the next 'lookahead' period.

# Calculate rolling returns over the lookback period
# ROC computes rate of change: (price_t / price_{t-n}) - 1
# We use ROC on prices directly to mimic momentum calculation over lookback
tsm_momentum <- ROC(prices, n = lookback_tsm, type = "discrete")
tsm_momentum <- na.omit(tsm_momentum)

# Generate signals: 1 if momentum is positive, 0 otherwise (invest or cash)
# Signal is generated at time t, applied for t+1 to t+lookahead_tsm
# Simple TSM: If lookback return > 0, invest for next 'lookahead' days.
# A more common TSM implementation buys if return > 0, sells/shorts if < 0.
# Project description just says "Time Series Momentum", usually implies long/short or long/flat.
# Let's assume long/flat: invest (1) if past return > 0, else cash (0).
# Note: This is a simplified application. Real TSM often ranks assets or uses volatility scaling.

# Forward fill the signal for the lookahead period
# Create signal: 1 if momentum > 0, else 0. Lagged by 1 to apply for next day.
signal_tsm <- lag(ifelse(tsm_momentum > 0, 1, 0), 1)
signal_tsm[is.na(signal_tsm)] <- 0 # Start with cash

# Align signal with price index
signal_tsm <- signal_tsm[index(daily_arith_returns)]
signal_tsm[is.na(signal_tsm)] <- 0 # Ensure no NAs at start/end

# --- 4b. SMA Strategy Signal Generation --- [cite: 6]

# Function to generate SMA signals
generate_sma_signal <- function(prices, short_n, long_n) {
  sma_short <- SMA(prices, n = short_n)
  sma_long <- SMA(prices, n = long_n)
  # Signal: 1 if short SMA > long SMA, else 0 (invest or cash)
  signal <- ifelse(sma_short > sma_long, 1, 0)
  # Lag signal by 1 day (trade on next day's price)
  signal <- lag(signal, 1)
  signal[is.na(signal)] <- 0 # Start with cash
  colnames(signal) <- paste0("SMA_", short_n, "_", long_n)
  return(signal)
}

# Generate signals for all SMA pairs
sma_signals <- list()
for (pair in sma_pairs) {
  short_period <- pair[1]
  long_period <- pair[2]
  sma_signals[[paste0("SMA_", short_period, "_", long_period)]] <-
    generate_sma_signal(prices, short_period, long_period)
}

# --- 4c. Calculate SMA Performance (Full Sample) & Choose Best ---

cat("Calculating Sharpe Ratios for SMA strategies over the full sample...\n")
sma_performance <- data.frame(
  SMA_Combo = character(),
  Sharpe_Ratio = numeric()
)

# Align price returns with the longest signal period possible
aligned_returns <- daily_arith_returns[index(sma_signals[[length(sma_signals)]])] # Use longest SMA to align

for (i in 1:length(sma_signals)) {
  combo_name <- names(sma_signals)[i]
  signal <- sma_signals[[i]]
  
  # Align signal and returns
  # Ensure signal index is within price index
  signal <- signal[index(aligned_returns)]
  signal[is.na(signal)] <- 0 # Fill any potential misalignments with 0
  
  # Calculate strategy returns (0 when signal is 0)
  strategy_returns <- signal * aligned_returns
  colnames(strategy_returns) <- combo_name
  
  # Calculate Sharpe Ratio
  sharpe <- calculate_sharpe(strategy_returns)
  
  sma_performance <- rbind(sma_performance, data.frame(
    SMA_Combo = combo_name,
    Sharpe_Ratio = sharpe
  ))
}

print(sma_performance)

# Choose the best SMA combo based on Sharpe Ratio
best_sma_combo_name <- sma_performance$SMA_Combo[which.max(sma_performance$Sharpe_Ratio)]
best_sma_signal <- sma_signals[[best_sma_combo_name]]

# Extract best periods for clarity
best_periods <- as.numeric(unlist(strsplit(gsub("SMA_", "", best_sma_combo_name), "_")))
best_short_sma <- best_periods[1]
best_long_sma <- best_periods[2]


cat("\nBest SMA Combo based on Sharpe Ratio:", best_sma_combo_name, "\n")
cat("--- End of Step 4 ---\n\n")


# --- 5. Compare Performance (Rolling Windows) --- [cite: 7]

cat("--- Step 5: Comparing Strategy Performance over Rolling 10-Year Windows ---\n")

# Ensure prices/returns are xts objects
if (!is.xts(prices)) prices <- xts(prices, order.by = as.Date(rownames(prices)))
if (!is.xts(daily_arith_returns)) daily_arith_returns <- xts(daily_arith_returns, order.by = index(prices)[-1])

# Prepare data frame to store results
results_rolling <- data.frame()

# Get dates for rolling windows
dates <- index(daily_arith_returns)
n_days <- length(dates)
window_length_days <- floor(window_years * 252) # Approximate days in 10 years

cat("Starting rolling window analysis (Window size:", window_length_days, "days)...\n")
pb <- txtProgressBar(min = 0, max = n_days - window_length_days + 1, style = 3) # Progress bar

# Loop through all possible start dates for a 10-year window
for (i in 1:(n_days - window_length_days + 1)) {
  # Define window start and end dates
  window_start_date <- dates[i]
  window_end_date <- dates[i + window_length_days - 1]
  
  # Subset data for the current window
  window_prices <- prices[paste0(window_start_date, "/", window_end_date)]
  window_returns_arith <- daily_arith_returns[paste0(window_start_date, "/", window_end_date)]
  
  if(nrow(window_returns_arith) < window_length_days * 0.9) next # Skip if window too small (e.g., due to holidays)
  
  # --- 5a. B&H Performance for the window ---
  # Note: B&H daily returns are just the arithmetic returns
  perf_bh <- tryCatch({
    list(
      Annualized_Return = Return.annualized(window_returns_arith)[1,1],
      Sharpe_Ratio = calculate_sharpe(window_returns_arith), # SharpeRatio.annualized(window_returns_arith, Rf = risk_free_rate/252, scale=252)[1,1] -> gives error sometimes
      MDD = maxDrawdown(window_returns_arith)
    )
  }, error = function(e) list(Annualized_Return=NA, Sharpe_Ratio=NA, MDD=NA))
  
  
  # --- 5b. TSM Performance for the window ---
  # Need to generate signals *based on data available up to the start of the window* or *within* the window?
  # Project implies using the same strategy parameters throughout. Apply full-sample signals to the window.
  window_signal_tsm <- signal_tsm[index(window_returns_arith)]
  window_signal_tsm[is.na(window_signal_tsm)] <- 0 # Handle potential NAs
  
  tsm_strategy_returns <- window_signal_tsm * window_returns_arith
  colnames(tsm_strategy_returns) <- "TSM_Return"
  tsm_strategy_returns[is.na(tsm_strategy_returns)] <- 0 # Critical: 0 return when signal is 0 or NA
  
  perf_tsm <- tryCatch({
    list(
      Annualized_Return = Return.annualized(tsm_strategy_returns)[1,1],
      Sharpe_Ratio = calculate_sharpe(tsm_strategy_returns),
      MDD = maxDrawdown(tsm_strategy_returns)
    )
  }, error = function(e) list(Annualized_Return=NA, Sharpe_Ratio=NA, MDD=NA))
  
  
  # --- 5c. Best SMA Performance for the window ---
  # Apply the best full-sample SMA signal to the window
  window_signal_sma <- best_sma_signal[index(window_returns_arith)]
  window_signal_sma[is.na(window_signal_sma)] <- 0 # Handle potential NAs
  
  sma_strategy_returns <- window_signal_sma * window_returns_arith
  colnames(sma_strategy_returns) <- "SMA_Return"
  sma_strategy_returns[is.na(sma_strategy_returns)] <- 0 # Critical: 0 return when signal is 0 or NA
  
  perf_sma <- tryCatch({
    list(
      Annualized_Return = Return.annualized(sma_strategy_returns)[1,1],
      Sharpe_Ratio = calculate_sharpe(sma_strategy_returns),
      MDD = maxDrawdown(sma_strategy_returns)
    )
  }, error = function(e) list(Annualized_Return=NA, Sharpe_Ratio=NA, MDD=NA))
  
  
  # --- 5d. Store results for the window ---
  results_rolling <- rbind(results_rolling, data.frame(
    Window_Start = window_start_date,
    Window_End = window_end_date,
    BH_AnnReturn = perf_bh$Annualized_Return,
    BH_Sharpe = perf_bh$Sharpe_Ratio,
    BH_MDD = perf_bh$MDD,
    TSM_AnnReturn = perf_tsm$Annualized_Return,
    TSM_Sharpe = perf_tsm$Sharpe_Ratio,
    TSM_MDD = perf_tsm$MDD,
    SMA_AnnReturn = perf_sma$Annualized_Return,
    SMA_Sharpe = perf_sma$Sharpe_Ratio,
    SMA_MDD = perf_sma$MDD
  ))
  
  setTxtProgressBar(pb, i) # Update progress bar
}
close(pb) # Close progress bar
cat("\nRolling window analysis complete.\n")

# Remove rows with NA results if any (can happen at edges or if errors occurred)
results_rolling <- na.omit(results_rolling)

# --- 5e. Analyze Distributions and Statistical Significance --- [cite: 8]

cat("\nAnalyzing distributions of performance metrics across all windows...\n")

# Summary statistics for Annualized Returns
summary_returns <- results_rolling %>%
  select(BH_AnnReturn, TSM_AnnReturn, SMA_AnnReturn) %>%
  summary()
cat("\nSummary Statistics for Annualized Returns:\n")
print(summary_returns)

# Summary statistics for Sharpe Ratios
summary_sharpe <- results_rolling %>%
  select(BH_Sharpe, TSM_Sharpe, SMA_Sharpe) %>%
  summary()
cat("\nSummary Statistics for Sharpe Ratios:\n")
print(summary_sharpe)

# Summary statistics for Max Drawdowns
summary_mdd <- results_rolling %>%
  select(BH_MDD, TSM_MDD, SMA_MDD) %>%
  summary()
cat("\nSummary Statistics for Maximum Drawdowns:\n")
print(summary_mdd)





# --- Plotting Distributions (Optional but recommended for report) ---
# Requires ggplot2
# Reshape data for easier plotting
library(tidyr)
results_long_ret <- results_rolling %>% select(Window_Start, BH_AnnReturn, TSM_AnnReturn, SMA_AnnReturn) %>% gather(key="Strategy", value="Return", -Window_Start)
results_long_sharpe <- results_rolling %>% select(Window_Start, BH_Sharpe, TSM_Sharpe, SMA_Sharpe) %>% gather(key="Strategy", value="Sharpe", -Window_Start)

# Plot Histograms
hist_returns <- ggplot(results_long_ret, aes(x=Return, fill=Strategy)) + geom_histogram(alpha=0.6, position="identity", bins=30) + facet_wrap(~Strategy) + ggtitle("Distribution of Annualized Returns (10Y Rolling Windows)") + theme_minimal()
print(hist_returns)

hist_sharpe <- ggplot(results_long_sharpe, aes(x=Sharpe, fill=Strategy)) + geom_histogram(alpha=0.6, position="identity", bins=30) + facet_wrap(~Strategy) + ggtitle("Distribution of Sharpe Ratios (10Y Rolling Windows)") + theme_minimal()
print(hist_sharpe)

# --- Statistical Significance Tests ---
# Using paired t-tests because the returns/Sharpes are calculated over the same time windows

cat("\nPerforming Paired t-tests for difference in means...\n")

# Compare Annualized Returns
t_test_ret_tsm_bh <- t.test(results_rolling$TSM_AnnReturn, results_rolling$BH_AnnReturn, paired = TRUE)
t_test_ret_sma_bh <- t.test(results_rolling$SMA_AnnReturn, results_rolling$BH_AnnReturn, paired = TRUE)
t_test_ret_tsm_sma <- t.test(results_rolling$TSM_AnnReturn, results_rolling$SMA_AnnReturn, paired = TRUE)

cat("\nPaired t-test: TSM AnnReturn vs BH AnnReturn\n")
print(t_test_ret_tsm_bh)
cat("\nPaired t-test: SMA AnnReturn vs BH AnnReturn\n")
print(t_test_ret_sma_bh)
cat("\nPaired t-test: TSM AnnReturn vs SMA AnnReturn\n")
print(t_test_ret_tsm_sma)


# Compare Sharpe Ratios
t_test_sharpe_tsm_bh <- t.test(results_rolling$TSM_Sharpe, results_rolling$BH_Sharpe, paired = TRUE)
t_test_sharpe_sma_bh <- t.test(results_rolling$SMA_Sharpe, results_rolling$BH_Sharpe, paired = TRUE)
t_test_sharpe_tsm_sma <- t.test(results_rolling$TSM_Sharpe, results_rolling$SMA_Sharpe, paired = TRUE)

cat("\nPaired t-test: TSM Sharpe vs BH Sharpe\n")
print(t_test_sharpe_tsm_bh)
cat("\nPaired t-test: SMA Sharpe vs BH Sharpe\n")
print(t_test_sharpe_sma_bh)
cat("\nPaired t-test: TSM Sharpe vs SMA Sharpe\n")
print(t_test_sharpe_tsm_sma)

# Interpretation: Look at the p-value. A small p-value (e.g., < 0.05) indicates a statistically significant difference between the means of the two compared groups.

cat("\n--- End of Step 5 ---\n\n")


# --- End of Script ---
# The data frame 'results_rolling' contains the key metrics for each window.
# Use this data frame and the printed summaries/test results for your report.