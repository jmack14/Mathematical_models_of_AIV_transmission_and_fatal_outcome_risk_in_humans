# ============================================================
# Sensitivity analysis: incomplete reporting of human AIV deaths
# ============================================================

# Load infection estimates
dat38 <- read.csv("transformed_samples_38.csv", stringsAsFactors = FALSE)

# ------------------------------------------------------------
# Input observed deaths
# ------------------------------------------------------------

reported_deaths <- 16.7


# ------------------------------------------------------------
# Calculate severity proxy under reporting scenarios
# ------------------------------------------------------------

reporting_scenarios <- c(
  "100%" = 1.0,
  "80%"  = 0.8,
  "60%"  = 0.6,
  "40%"  = 0.4
)


# ------------------------------------------------------------
# Quantile function matching the main analysis
# ------------------------------------------------------------

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


# ------------------------------------------------------------
# Calculate severity proxies and uncertainty intervals
# ------------------------------------------------------------

proxy_results <- lapply(
  reporting_scenarios,
  function(reporting_fraction) {
    
    # Correct deaths for incomplete reporting
    adjusted_deaths <- reported_deaths / reporting_fraction
    
    # Severity proxy per 10,000 infections
    proxy <- 100 * adjusted_deaths / dat38$meannh
    
    data.frame(
      Reporting = reporting_fraction,
      Median_IFR = get_empirical_quantile(proxy, 0.50),
      Lower_95   = get_empirical_quantile(proxy, 0.025),
      Upper_95   = get_empirical_quantile(proxy, 0.975)
    )
  }
)


# ------------------------------------------------------------
# Combine results
# ------------------------------------------------------------

proxy_results <- do.call(rbind, proxy_results)

rownames(proxy_results) <- names(reporting_scenarios)


# ------------------------------------------------------------
# Print results
# ------------------------------------------------------------

print(proxy_results)


# ------------------------------------------------------------
# Save results
# ------------------------------------------------------------

write.csv(
  proxy_results,
  "Reporting_sensitivity.csv",
  row.names = TRUE
)
