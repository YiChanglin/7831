################################################################################
# Changlin Yi
# cy2578@nyu.edu
# FRE 7831 HW-3
################################################################################

rm(list = ls())

# Core packages
library(quantmod)
library(PortfolioAnalytics)
library(ROI)
library(ROI.plugin.quadprog)
library(ROI.plugin.glpk)
library(xts)
library(lubridate)


# Top 10 SOXX ETF constituents
tkr_list <- c("AVGO", "NVDA", "TXN", "AMD", "QCOM", 
              "KLAC", "AMAT", "LRCX", "MPWR", "ASML")

# Use 5 years of data: set endDate to today's date and startDate 5 years ago.
endDate   <- Sys.Date()                    # e.g., 2025-04-15 (today)
startDate <- endDate - years(5)              # 5 years back

# Download historical daily adjusted close prices from Yahoo Finance
getSymbols(tkr_list, src = "yahoo", from = startDate, to = endDate)

# Merge the adjusted close prices into one xts object
prices <- do.call(merge, lapply(tkr_list, function(sym) Ad(get(sym))))
colnames(prices) <- tkr_list

# Compute daily returns (using the discrete return method)
returns <- ROC(prices, type = "discrete")[-1,]

# Use the column names from returns as asset names
fund.names <- colnames(returns)

# Create a portfolio specification with full investment and long-only constraints.
pspec <- portfolio.spec(assets = fund.names)
pspec <- add.constraint(portfolio = pspec, type = "full_investment")
pspec <- add.constraint(portfolio = pspec, type = "long_only")
pspec <- add.objective(portfolio = pspec, type = "return", name = "mean")
pspec <- add.objective(portfolio = pspec, type = "risk", name = "StdDev")

# Look-back: 3 years ≈ 252 trading days/year * 3 = 756 days.
# Look-ahead: 1 month ≈ 21 trading days.
lookback_days  <- 252 * 3  # 756 days for optimization
lookahead_days <- 21       # 21 days to test out-of-sample performance

# Determine available observations
n_obs <- nrow(returns)
# Create a sequence of start indices for rolling windows.
# Here we roll forward by 21 days each time.
roll_starts <- seq(1, n_obs - (lookback_days + lookahead_days) + 1, by = lookahead_days)

# Pre-allocate vectors to store performance metrics (annualized)
n_rolls <- length(roll_starts)
mvp_ann_return <- rep(NA, n_rolls)
mvp_ann_vol    <- rep(NA, n_rolls)
mvp_ann_sharpe <- rep(NA, n_rolls)
ewp_ann_return <- rep(NA, n_rolls)
ewp_ann_vol    <- rep(NA, n_rolls)
ewp_ann_sharpe <- rep(NA, n_rolls)
# Pre-allocate a vector to store the count of significant weights (greater than 5%)
sig_counts <- rep(NA, n_rolls)

for(i in seq_along(roll_starts)) {
  
  # Define indices for look-back (in-sample) and look-ahead (out-of-sample) periods
  lb_start <- roll_starts[i]
  lb_end   <- lb_start + lookback_days - 1
  la_start <- lb_end + 1
  la_end   <- lb_end + lookahead_days
  
  # Ensure we do not exceed the dataset
  if(la_end > n_obs) break
  
  # Subset the returns data for the look-back and look-ahead periods
  lb_returns <- returns[lb_start:lb_end, ]
  la_returns <- returns[la_start:la_end, ]
  
  mvp_opt <- optimize.portfolio(R = lb_returns,
                                portfolio = pspec,
                                optimize_method = "ROI",
                                maxSR = TRUE,
                                trace = FALSE)
  mvp_weights <- mvp_opt$weights  # these are the optimized weights
  
  # --- Count the number of significant weights (> 5%) ---
  sig_counts[i] <- sum(as.numeric(mvp_weights) > 0.05)
  
  # Compute portfolio daily returns for the MVP portfolio
  mvp_daily <- as.numeric(la_returns %*% mvp_weights)
  
  # Compute performance metrics for MVP (annualized return, volatility, Sharpe ratio)
  mvp_mean   <- mean(mvp_daily)
  mvp_sd     <- sd(mvp_daily)
  mvp_ann_return[i] <- mvp_mean * 252
  mvp_ann_vol[i]    <- mvp_sd * sqrt(252)
  mvp_ann_sharpe[i] <- ifelse(mvp_ann_vol[i] != 0, mvp_ann_return[i] / mvp_ann_vol[i], NA)
  
  # Compute the same for an Equally Weighted Portfolio (EWP)
  num_assets <- length(fund.names)
  ewp_weights <- rep(1/num_assets, num_assets)
  ewp_daily   <- as.numeric(la_returns %*% ewp_weights)
  ewp_mean    <- mean(ewp_daily)
  ewp_sd      <- sd(ewp_daily)
  ewp_ann_return[i] <- ewp_mean * 252
  ewp_ann_vol[i]    <- ewp_sd * sqrt(252)
  ewp_ann_sharpe[i] <- ifelse(ewp_ann_vol[i] != 0, ewp_ann_return[i] / ewp_ann_vol[i], NA)
  
  # Optionally, print progress:
  cat("Window", i, "of", n_rolls, ": MVP Sharpe =", round(mvp_ann_sharpe[i], 3),
      "; EWP Sharpe =", round(ewp_ann_sharpe[i], 3), "\n")
}

# Create a results data frame using the date that marks the beginning of the look-ahead period
result_dates <- index(returns)[roll_starts + lookback_days]
results <- data.frame(Date       = result_dates,
                      MVP_Return = mvp_ann_return,
                      MVP_Vol    = mvp_ann_vol,
                      MVP_Sharpe = mvp_ann_sharpe,
                      EWP_Return = ewp_ann_return,
                      EWP_Vol    = ewp_ann_vol,
                      EWP_Sharpe = ewp_ann_sharpe)
print(results)

# Plot the time series of the number of assets with weight > 5%
plot(result_dates, sig_counts, type = "o", col = "blue",
     xlab = "Date", ylab = "Number of Assets with Weight > 5%",
     main = "Evolution of Significant Portfolio Weights (MVP)")

# Compare Sharpe Ratio Distributions and Statistical Tests 
par(mfrow = c(1, 2))  # Side-by-side plots
hist(mvp_ann_sharpe, main = "MVP Sharpe Distribution",
     xlab = "MVP Annualized Sharpe", col = "lightblue", breaks = 10)
hist(ewp_ann_sharpe, main = "EWP Sharpe Distribution",
     xlab = "EWP Annualized Sharpe", col = "lightgreen", breaks = 10)
par(mfrow = c(1, 1))  # Reset plotting layout

# Perform paired t-test and Wilcoxon signed-rank test comparing MVP vs. EWP Sharpe ratios
t_test_result <- t.test(mvp_ann_sharpe, ewp_ann_sharpe, paired = TRUE)
wilcox_test_result <- wilcox.test(mvp_ann_sharpe, ewp_ann_sharpe, paired = TRUE)

cat("\nPaired t-test on Sharpe Ratios:\n")
print(t_test_result)

cat("\nPaired Wilcoxon Test on Sharpe Ratios:\n")
print(wilcox_test_result)
