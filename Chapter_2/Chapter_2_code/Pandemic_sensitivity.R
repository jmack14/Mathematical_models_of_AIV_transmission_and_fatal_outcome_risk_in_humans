# ============================================================
# Sensitivity analysis: historical pandemic dataset
# Excluding uncertain pandemics (1889 and/or 2009)
# ============================================================

set.seed(20260828)

library(bbmle)

fit_interpandemic <- function(pan_year, scenario){
  
  # ----------------------------------------------------------
  # Interpandemic intervals
  # ----------------------------------------------------------
  
  dd <- diff(pan_year)
  
  # Remove last pandemic (no following interval)
  dates <- head(pan_year, -1)
  year  <- dates - dates[1]
  
  # ==========================================================
  # Model without effect of date 
  # ==========================================================
  
  delta1_fixed <- 0
  
  Mgammanll0 <- function(alpha = 1, delta0 = 30){
    
    nll <- -sum(
      log(
        (dd^(alpha - 1)) *
          exp(-dd / ((delta0 + delta1_fixed * year) / alpha)) /
          ((((delta0 + delta1_fixed * year) / alpha)^alpha) *
             gamma(alpha))
      )
    )
    
    return(nll)
  }
  
  mod0 <- mle2(
    Mgammanll0,
    start = list(
      alpha = 2.7156292423542,
      delta0 = mean(dd)
    ),
    control = list(maxit = 1000)
  )
  
  alpha0 <- unname(coef(mod0)["alpha"])
  delta  <- unname(coef(mod0)["delta0"])
  theta0 <- delta / alpha0
  
  # ==========================================================
  # Model with effect of date 
  # ==========================================================
  
  Mgammanll <- function(alpha = 1,
                        delta0 = 30,
                        delta1 = 0){
    
    nll <- -sum(
      log(
        (dd^(alpha - 1)) *
          exp(-dd / ((delta0 + delta1 * year) / alpha)) /
          ((((delta0 + delta1 * year) / alpha)^alpha) *
             gamma(alpha))
      )
    )
    
    return(nll)
    
  }
  
  mod <- mle2(
    Mgammanll,
    start = list(
      alpha = 5.8,
      delta0 = mean(dd),
      delta1 = -0.1
    ),
    control = list(maxit = 1000)
  )
  
  alpha  <- unname(coef(mod)["alpha"])
  delta0 <- unname(coef(mod)["delta0"])
  delta1 <- unname(coef(mod)["delta1"])
  theta  <- delta0 / alpha
  
  # ----------------------------------------------------------
  # Save parameter estimates for IFR analysis
  # ----------------------------------------------------------
  
  save(
    mod0,
    mod,
    theta0,
    theta,
    delta,
    delta0,
    delta1,
    file = paste0(
      "Interpandemic_",
      gsub("[ &]", "_", scenario),
      ".RData"
    )
  )
  
  # ----------------------------------------------------------
  # Return summary table
  # ----------------------------------------------------------
  
  data.frame(
    Scenario = scenario,
    Pandemics = length(pan_year),
    Mean_Interval = round(mean(dd), 1),
    alpha0 = round(alpha0, 3),
    theta0 = round(theta0, 3),
    alpha = round(alpha, 3),
    theta = round(theta, 3),
    delta1 = round(delta1, 3),
    stringsAsFactors = FALSE
  )
  
}

# ============================================================
# Historical pandemic scenarios
# ============================================================

baseline <- c(1781, 1830, 1889, 1918, 1957, 1968, 2009)

no1889 <- c(1781, 1830, 1918, 1957, 1968, 2009)

no2009 <- c(1781, 1830, 1889, 1918, 1957, 1968)

no1889_2009 <- c(1781, 1830, 1918, 1957, 1968)

# ============================================================
# Run sensitivity analysis
# ============================================================

results <- rbind(
  fit_interpandemic(baseline, "Baseline"),
  fit_interpandemic(no1889, "Exclude1889"),
  fit_interpandemic(no2009, "Exclude2009"),
  fit_interpandemic(no1889_2009, "Exclude1889_2009")
)

print(results)

write.csv(
  results,
  "Historical_pandemic_sensitivity.csv",
  row.names = FALSE
)

# ============================================================
# Run the 38-year IFR scenario for each pandemic sensitivity analysis
# ============================================================

num_samples <- 1000

scenario_files <- c(
  Exclude1889 = "Interpandemic_Exclude1889.RData",
  Exclude2009 = "Interpandemic_Exclude2009.RData",
  Exclude1889_2009 = "Interpandemic_Exclude1889_2009.RData"
)

for (s in names(scenario_files)) {
  
  load(scenario_files[s])
  
  message("Running pandemic scenario: ", s)
  
  res <- run_scenario(
    list(
      R0 = list(rate = 20),
      a = c(0.000012, 0.000024),
      Rstar = c(1, 2, 1.1),
      psi = list(
        shape = delta / theta0,
        scale = theta0
      )
    ),
    num_samples
  )
  
  write.csv(
    res$transformed_samples,
    paste0("transformed_samples_38_", s, ".csv"),
    row.names = FALSE,
    quote = FALSE
  )
  
  write.csv(
    res$sort_output,
    paste0("mean_nh_38_", s, ".csv"),
    row.names = FALSE,
    quote = FALSE
  )
}

# ============================================================
# Summary of annual human infections across scenarios
# ============================================================

files <- list.files(pattern = "mean_nh_38_.*\\.csv$")

comparison_results <- lapply(files, function(f) {
  
  dat <- read.csv(f)
  
  data.frame(
    Scenario = sub("^mean_nh_38_(.*)\\.csv$", "\\1", f),
    Median = median(dat$meannh),
    Lower95 = quantile(dat$meannh, 0.025),
    Upper95 = quantile(dat$meannh, 0.975),
    Median_IFR = 100 * 20.6 / median(dat$meannh),
    IFR_Lower95 = 100 * 20.6 / quantile(dat$meannh, 0.975),
    IFR_Upper95 = 100 * 20.6 / quantile(dat$meannh, 0.025)
  )
  
})

comparison_results <- do.call(rbind, comparison_results)

print(comparison_results)

write.csv(
  comparison_results,
  "Pandemic_sensitivity.csv",
  row.names = FALSE
)
