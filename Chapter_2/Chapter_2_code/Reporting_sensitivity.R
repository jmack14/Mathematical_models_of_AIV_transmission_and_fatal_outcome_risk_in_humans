# ============================================================
# Sensitivity analysis: incomplete reporting of human AIV deaths
# ============================================================

# Path to transformed samples relative to project root
dat38 <- read.csv(
  "C:/Users/marti/Desktop/Joshua/Code/Clean_code/Output/transformed_samples_38.csv",
  stringsAsFactors = FALSE
)

# Load infection estimates
dat38 <- read.csv(
  "C:/Users/ER/Desktop/Fall_2025/Grad_school/thesis/Code/Clean_code/Output/transformed_samples_38.csv"
)

# ------------------------------------------------------------
# Input observed deaths
# ------------------------------------------------------------

reported_deaths <- 20.6  


# ------------------------------------------------------------
# Calculate IFR under reporting scenarios
# ------------------------------------------------------------

reporting_scenarios <- c(
  "100%" = 1.0,
  "80%"  = 0.8,
  "60%"  = 0.6,
  "40%"  = 0.4
)


ifr_results <- lapply(
  reporting_scenarios,
  function(reporting_fraction) {
    
    # Correct deaths for incomplete reporting
    adjusted_deaths <- reported_deaths / reporting_fraction
    
    # IFR
    IFR <- 100 * adjusted_deaths / dat38$meannh
    
    data.frame(
      Reporting = reporting_fraction,
      Median_IFR = median(IFR),
      Lower_95 = quantile(IFR, 0.025),
      Upper_95 = quantile(IFR, 0.975)
    )
  }
)


ifr_results <- do.call(rbind, ifr_results)

rownames(ifr_results) <- names(reporting_scenarios)

print(ifr_results)


# Save
write.csv(
  ifr_results,
  "Reporting_sensitivity.csv",
  row.names = TRUE
)