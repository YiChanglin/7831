# Fama-French Factor Model Estimation Script
# -------------------------------------------
# 1. Clear workspace and load required libraries
rm(list = ls())
library(quantmod)
library(xts)
library(tseries)

# 2. Read the cleaned Fama–French 5-factor data (last 4 years)
ff_file <- "/Users/charlie/Desktop/NYU 2025 Spring/spring 下半/7831/HW/HW 4/F-F_Research_Data_5_Factors2021_2024.CSV"
ff <- read.csv(ff_file, header = TRUE, stringsAsFactors = FALSE)

# 3. Convert Date column: first to character, then parse as Date (YYYYMMDD)
ff$Date <- as.character(ff$Date)
ff$Date <- as.Date(ff$Date, format = "%Y%m%d")

# 4. Create an xts object and extract factor series
ff_xts <- xts(ff[, -1], order.by = ff$Date)
rmrf <- ff_xts$Mkt.RF
smb  <- ff_xts$SMB
hml  <- ff_xts$HML
rmw  <- ff_xts$RMW
cma  <- ff_xts$CMA
rf   <- ff_xts$RF

# 5. Define function to estimate factor models for a given ticker
estimate_models <- function(ticker, start_date, end_date) {
  # Download daily adjusted close prices and compute daily returns (%)
  stock_data <- tryCatch(
    getSymbols(ticker, src = "yahoo", from = start_date, to = end_date, auto.assign = FALSE),
    error = function(e) stop(paste("Error downloading data for", ticker, ":", e$message))
  )
  returns <- dailyReturn(Ad(stock_data)) * 100
  
  # Align with risk-free rate and compute excess returns
  df <- merge(returns, rf, join = "inner")
  colnames(df) <- c("Return", "RF")
  df$ExcessReturn <- df$Return - df$RF
  
  # Merge excess returns with factor series (inner joins)
  data_all <- merge(df$ExcessReturn, rmrf, join = "inner")
  data_all <- merge(data_all, smb,  join = "inner")
  data_all <- merge(data_all, hml,  join = "inner")
  data_all <- merge(data_all, rmw,  join = "inner")
  data_all <- merge(data_all, cma,  join = "inner")
  colnames(data_all) <- c("ExcessReturn", "Mkt.RF", "SMB", "HML", "RMW", "CMA")
  data_all <- na.omit(data_all)
  
  # Estimate models: CAPM, FF3, FF5
  capm <- lm(ExcessReturn ~ Mkt.RF,                     data = data_all)
  ff3  <- lm(ExcessReturn ~ Mkt.RF + SMB + HML,         data = data_all)
  ff5  <- lm(ExcessReturn ~ Mkt.RF + SMB + HML + RMW + CMA, data = data_all)
  
  # Print summaries
  cat("\n==============================\n")
  cat("Ticker:", ticker, "  Period:", start_date, "to", end_date, "\n")
  cat("-- CAPM (1-factor) Summary --\n"); print(summary(capm))
  cat("\n-- FF3 (3-factor) Summary --\n");   print(summary(ff3))
  cat("\n-- FF5 (5-factor) Summary --\n");   print(summary(ff5))
}

# 6. Specify tickers and date range
tickers    <- c("TSM", "OXY", "MRVI", "SJW")
start_date <- first(index(ff_xts))
end_date   <- last(index(ff_xts))

# 7. Run estimation for each ticker
for (tkr in tickers) {
  estimate_models(tkr, start_date, end_date)
}
