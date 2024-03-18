#exercise A

library(MASS)

# Define the function
predict_property_prices_manual <- function(y, X) {
  # Add a column of 1s to X for the intercept
  X <- cbind(1, X)
  
  # Calculate the coefficients (beta) using OLS formula: beta = (X'X)^(-1)X'y
  beta <- solve(t(X) %*% X) %*% t(X) %*% y
  
  # Calculate predictions
  y_hat <- X %*% beta
  
  # Calculate residuals
  residuals <- y - y_hat
  
  # Calculate error variance
  error_variance <- sum(residuals^2) / (nrow(X) - ncol(X))
  
  # Calculate standard errors of coefficients
  se_beta <- sqrt(diag(error_variance * solve(t(X) %*% X)))
  
  # Calculate t-values
  t_values <- beta / se_beta
  
  # Calculate p-values
  p_values <- 2 * pt(-abs(t_values), df = nrow(X) - ncol(X))
  
  # Calculate 95% confidence intervals
  conf_int_lower <- beta - qt(0.975, nrow(X) - ncol(X)) * se_beta
  conf_int_upper <- beta + qt(0.975, nrow(X) - ncol(X)) * se_beta
  conf_intervals <- cbind(conf_int_lower, conf_int_upper)
  
  return(list(
    Coefficients = beta,
    Error_Variance = error_variance,
    Standard_Errors = se_beta,
    t_Values = t_values,
    p_Values = p_values,
    Confidence_Intervals = conf_intervals
  ))
}


# Assuming y is your dependent variable vector and X is a matrix with each column being one covariate
# For example:
y <- Boston$medv
X <- as.matrix(Boston[, c("lstat", "rm", "age", "tax", "dis")])

result_manual <- predict_property_prices_manual(y, X)
print(result_manual)
