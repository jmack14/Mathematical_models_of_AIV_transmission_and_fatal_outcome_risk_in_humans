# ============================================================
# Sensitivity analysis of R0, Rstar, and a assumptions
# ============================================================

# Fixed random seed for reproducibility
set.seed(20260828)

sensitivity_scenarios <- list(
  
  # R0 assumptions
  `R0 rate = 15` = list(
    R0 = list(rate = 15),
    a = c(0.000012, 0.000024),
    Rstar = c(1, 2, 1.1),
    psi = list(shape = delta/theta0, scale = theta0)
  ),
  
  `R0 rate = 25` = list(
    R0 = list(rate = 25),
    a = c(0.000012, 0.000024),
    Rstar = c(1, 2, 1.1),
    psi = list(shape = delta/theta0, scale = theta0)
  ),
  
  # Rstar assumptions
  `Rstar mode = 1.05` = list(
    R0 = list(rate = 20),
    a = c(0.000012, 0.000024),
    Rstar = c(1, 2, 1.05),
    psi = list(shape = delta/theta0, scale = theta0)
  ),
  
  `Rstar mode = 1.20` = list(
    R0 = list(rate = 20),
    a = c(0.000012, 0.000024),
    Rstar = c(1, 2, 1.20),
    psi = list(shape = delta/theta0, scale = theta0)
  ),
  
  # a assumptions
  `a half` = list(
    R0 = list(rate = 20),
    a = c(0.000006, 0.000012),
    Rstar = c(1, 2, 1.1),
    psi = list(shape = delta/theta0, scale = theta0)
  ),
  
  `a double` = list(
    R0 = list(rate = 20),
    a = c(0.000024, 0.000048),
    Rstar = c(1, 2, 1.1),
    psi = list(shape = delta/theta0, scale = theta0)
  )
)


# ============================================================
# Empirical quantile function
# ============================================================

get_empirical_quantile <- function(x, p) {
  
  # Sort values from lowest to highest
  x_sorted <- sort(x)
  
  # Cumulative probabilities
  cum.prob <- seq(
    1 / length(x_sorted),
    1,
    1 / length(x_sorted)
  )
  
  # First observation with cumulative probability >= p
  idx <- min(which(cum.prob >= p))
  
  x_sorted[idx]
}


# ============================================================
# Run sensitivity analysis
# ============================================================

sensitivity_results <- bind_rows(
  
  lapply(names(sensitivity_scenarios), function(s){
    
    # Run the model
    res <- run_scenario(
      sensitivity_scenarios[[s]],
      num_samples = 1000
    )
    
    # --------------------------------------------------------
    # Infection estimates
    # --------------------------------------------------------
    
    meannh <- res$transformed_samples$meannh
    
    median_inf <- get_empirical_quantile(meannh, 0.50)
    lower_inf  <- get_empirical_quantile(meannh, 0.025)
    upper_inf  <- get_empirical_quantile(meannh, 0.975)
    
    
    # --------------------------------------------------------
    # Severity proxy
    # --------------------------------------------------------
    
    proxy <- 100 * 16.7 / meannh
    
    # Directly calculate quantiles from proxy distribution
    median_proxy <- get_empirical_quantile(proxy, 0.50)
    lower_proxy  <- get_empirical_quantile(proxy, 0.025)
    upper_proxy  <- get_empirical_quantile(proxy, 0.975)
    
    
    # --------------------------------------------------------
    # Return results
    # --------------------------------------------------------
    
    data.frame(
      Scenario = s,
      
      # Mean and uncertainty for human infections
      Mean = mean(meannh),
      Median = median_inf,
      Lower95 = lower_inf,
      Upper95 = upper_inf,
      
      # Mean and uncertainty for severity / IFR proxy
      Median_proxy = median_proxy,
      proxy_Lower95 = lower_proxy,
      proxy_Upper95 = upper_proxy
    )
    
  })
  
)


# ============================================================
# Print results
# ============================================================

print(sensitivity_results)


# ============================================================
# Save results
# ============================================================

write.csv(
  sensitivity_results,
  "Parameter_sensitivity.csv",
  row.names = FALSE
)
