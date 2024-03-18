#exercise A

library(MASS)

Boston

# Load the MASS package to access the Boston dataset
library(MASS)

# Define the function
predict_property_prices <- function(dependent_var, covariates, data) {
  # Construct the formula for the linear model
  formula <- as.formula(paste(dependent_var, "~", paste(covariates, collapse = "+")))
  
  # Fit the linear model
  model <- lm(formula, data = data)
  
  # Get the summary of the model
  summary_model <- summary(model)
  
  # Extract the coefficients summary
  coefficients_summary <- summary_model$coefficients
  
  # Calculate the error variance
  error_variance <- summary_model$sigma^2
  
  # 95% Confidence Intervals
  conf_intervals <- confint(model, level = 0.95)
  
  # Construct the return list
  result <- list(
    OLS_Estimates = list(
      Intercept = coefficients_summary[1, 1],
      Slope_Parameters = coefficients_summary[-1, 1],
      Error_Variance = error_variance
    ),
    Test_Statistics = list(
      t_Values = coefficients_summary[, 3],
      p_Values = coefficients_summary[, 4]
    ),
    Confidence_Intervals = conf_intervals
  )
  
  return(result)
}

# Example of using the function with the Boston dataset
# Let's say we're predicting 'medv' (Median value of owner-occupied homes in $1000's) using 'lstat' (lower status of the population), 
# 'rm' (average number of rooms per dwelling), 'age' (proportion of owner-occupied units built prior to 1940), and 'tax' (full-value property-tax rate per $10,000)
result <- predict_property_prices("medv", c("lstat", "rm", "age", "tax"), Boston)

# Print the result
print(result)
