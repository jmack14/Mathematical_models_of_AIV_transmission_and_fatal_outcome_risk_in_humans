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
# Run sensitivity analysis
# ============================================================

sensitivity_results <- bind_rows(
  
  lapply(names(sensitivity_scenarios), function(s){
    
    res <- run_scenario(
      sensitivity_scenarios[[s]],
      num_samples = 1000
    )
    
    median_inf <- res$quarts$n[3]
    lower_inf  <- res$quarts$n[1]
    upper_inf  <- res$quarts$n[5]
    
    data.frame(
      Scenario = s,
      Mean = mean(res$transformed_samples$meannh),
      Median = median_inf,
      Lower95 = lower_inf,
      Upper95 = upper_inf,
      Median_IFR = 100 * 20.6 / median_inf,
      IFR_Lower95 = 100 * 20.6 / upper_inf,
      IFR_Upper95 = 100 * 20.6 / lower_inf
    )
    
  })
  
)

print(sensitivity_results)

# ============================================================
# Save results
# ============================================================

write.csv(
  sensitivity_results,
  "Parameter_sensitivity.csv",
  row.names = FALSE
)