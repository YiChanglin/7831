# Load required packages
if (!require("quantmod", quietly = TRUE)) {
  install.packages("quantmod")
  library(quantmod)
}
if (!require("moments", quietly = TRUE)) {
  install.packages("moments")
  library(moments)
}

# Choose a unique security ticker (not AAPL or MSFT)
ticker <- "AMZN"  # change as needed

# Define the start date for four years of data
start_date <- Sys.Date() - 365 * 4

# Download data from Yahoo Finance (adjusted closing prices will be included)
data <- getSymbols(ticker, src = "yahoo", from = start_date, auto.assign = FALSE)

# Extract the Adjusted closing prices
prices <- Ad(data)

# Calculate daily returns (proportion returns)
returns <- dailyReturn(prices)
returns <- na.omit(returns)  # Remove missing values

# Compute statistics: mean, standard deviation, skewness, and kurtosis
mean_return  <- mean(returns)
sigma_return <- sd(returns)
skew_return  <- skewness(returns)
kurt_return  <- kurtosis(returns)

# Display the statistics
cat("Statistics for", ticker, "daily returns (last 4 years):\n")
cat("Mean         :", mean_return, "\n")
cat("Std Dev      :", sigma_return, "\n")
cat("Skewness     :", skew_return, "\n")
cat("Kurtosis     :", kurt_return, "\n\n")

# Plot the distribution of returns
hist(returns, breaks = 50, probability = TRUE,
     main = paste("Distribution of", ticker, "Daily Returns"),
     xlab = "Daily Return", col = "lightblue", border = "white")

# Generate a sequence of x-values for the Normal density curve
x_vals <- seq(min(returns), max(returns), length.out = 100)

# Calculate the Normal density with the computed mean and sigma
normal_density <- dnorm(x_vals, mean = mean_return, sd = sigma_return)

# Overlay the Normal density curve on the histogram
lines(x_vals, normal_density, col = "red", lwd = 2)

# Add a legend for clarity
legend("topright", legend = c("Empirical Returns", "Normal Distribution"),
       col = c("lightblue", "red"), lwd = c(10, 2), bty = "n")


