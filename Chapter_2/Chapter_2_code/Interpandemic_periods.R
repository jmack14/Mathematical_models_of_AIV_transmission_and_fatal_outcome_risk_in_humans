# ============================================================
# Interpandemic period models
# ============================================================

# ============================================================
# Load libraries
# ============================================================

library(bbmle)
library(patchwork)
library(dplyr)
library(ggplot2)
library(ggtext)
library(scales)
library(lubridate)

# ============================================================
# Interpandemic period data
# ============================================================

pan_year <- c(1781, 1830, 1889, 1918, 1957, 1968, 2009)

# Interpandemic intervals
dd <- diff(pan_year)

# Remove last pandemic (no following interval)
dates <- head(pan_year, -1)
year  <- dates - dates[1]

# ============================================================
# Model with no effect of date on the interpandemic period
# ============================================================

delta1 <- 0

Mgammanll0 <- function(alpha = 1, delta0 = 30) {
  nll <- -sum(
    log(
      (dd^(alpha - 1)) *
        exp(-dd / ((delta0 + delta1 * year) / alpha)) /
        ((((delta0 + delta1 * year) / alpha)^alpha) * gamma(alpha))
    )
  )
  return(nll)
}

MgammaMLE0 <- function(...) {
  fit <- mle2(Mgammanll0, control = list(maxit = 1000))
  if (conv <- fit@details$convergence != 0) {
    stop(paste0("Convergence code", conv, "in gammaMLE"))
  }
  return(fit)
}

mmfit0 <- function() {
  dat <- dd
  MgammaMLE0(data = dat, start = list(alpha = 2.7156292423542, delta0 = mean(dat)))
}

mod0 <- mmfit0()
alpha0 <- unname(coef(mod0)["alpha"]) # shape parameter
delta  <- unname(coef(mod0)["delta0"]) # mean  
theta0 <- delta / alpha0 # scale parameter

summary(mod0)
coef(mod0)
confint(mod0)

# Determine confidence interval for theta0

# Log(theta0)
log_theta0 <- log(theta0)

# Gradient of log(theta0)
grad_log_theta0 <- c(
  d_alpha  = -1 / alpha0,
  d_delta0 =  1 / delta
)

# Variance-covariance matrix of MLEs 
V0 <- vcov(mod0)

# Delta-method variance on log scale
var_log_theta0 <- t(grad_log_theta0) %*% V0 %*% grad_log_theta0
se_log_theta0  <- sqrt(var_log_theta0)

# 95% CI on log scale
CI_log <- c(
  lower = log_theta0 - 1.96 * se_log_theta0,
  upper = log_theta0 + 1.96 * se_log_theta0
)

# Back-transform
CI_theta0 <- exp(CI_log)

theta0
CI_theta0

# ============================================================
# Model with effect of date on the interpandemic period
# ============================================================

Mgammanll <- function(alpha = 1, delta0 = 30, delta1 = 0) {
  nll <- -sum(
    log(
      (dd^(alpha - 1)) *
        exp(-dd / ((delta0 + delta1 * year) / alpha)) /
        ((((delta0 + delta1 * year) / alpha)^alpha) * gamma(alpha))
    )
  )
  return(nll)
}

MgammaMLE <- function(...) {
  fit <- mle2(Mgammanll, control = list(maxit = 1000))
  if (conv <- fit@details$convergence != 0) {
    stop(paste0("Convergence code", conv, "in gammaMLE"))
  }
  return(fit)
}

mmfit <- function() {
  dat <- dd
  MgammaMLE(data = dat, start = list(alpha = 5.8, delta0 = 53.5, delta1 = -0.1))
}

mod <- mmfit()
alpha  <- unname(coef(mod)["alpha"]) 
delta0 <- unname(coef(mod)["delta0"]) 
delta1 <- unname(coef(mod)["delta1"]) 
theta  <- delta0 / alpha # scale

summary(mod)
coef(mod)
confint(mod)
confint(mod, level = 0.75)

# Determine confidence interval for theta

# Log(theta)
log_theta <- log(theta)

# Gradient of log(theta) with respect to all parameters
grad_log_theta <- c(
  d_alpha  = -1 / alpha,
  d_delta0 =  1 / delta0,
  d_delta1 =  0
)

# Variance-covariance matrix of MLEs 
V1 <- vcov(mod)

# Delta-method variance on log scale
var_log_theta <- t(grad_log_theta) %*% V1 %*% grad_log_theta
se_log_theta  <- sqrt(var_log_theta)

# 95% CI on log scale
CI_log <- c(
  lower = log_theta - 1.96 * se_log_theta,
  upper = log_theta + 1.96 * se_log_theta
)

# Back-transform
CI_theta <- exp(CI_log)

theta
CI_theta

# ============================================================
# Likelihood ratio test
# ============================================================

LRT <- -2 * (logLik(mod0)[1] - logLik(mod)[1])
anova(mod0, mod)

# ============================================================
# Save results
# ============================================================

save(
  mod0, mod,
  theta0, theta, delta, delta0, delta1,
  file = "Interpandemic_period_results.RData"
)

# ============================================================
# Plot pandemic timeline: Fig1 
# ============================================================

Timeline <- data.frame(
  what = c(
    "1781 - Russian Catarrh", "1830 - Influenza Pandemic", "1889 - Russian flu\n(Subtype uncertain; reported as H3)",
    "1918 - Spanish flu (H1)", "1957 - Asiatic flu (H2)", "1968 - Hong Kong flu (H3)", "2009 - H1N1"
  ),
  when = ymd(c(
    "1781-12-01", "1830-12-01", "1889-12-01",
    "1918-12-01", "1957-12-01", "1968-12-01", "2009-12-01"
  )),
  event.type = ""
)

Periods <- data.frame(
  what = c(
    "1781 - Russian Catarrh/1830 - Influenza Pandemic",
    "1830 - Influenza Pandemic/1889 - Russian flu (H3)",
    "1889 - Russian flu\n(Subtype uncertain; reported as H3)/1918 - Spanish flu (H1)",
    "1918 - Spanish flu (H1)/1957 - Asiatic flu (H2)",
    "1957 - Asiatic flu (H2)/1968 - Hong Kong flu (H3)",
    "1968 - Hong Kong flu (H3)/2009 - H1N1"
  ),
  start = ymd(c(
    "1781-12-01", "1830-12-01", "1889-12-02",
    "1918-12-01", "1957-12-01", "1968-12-01"
  )),
  end = ymd(c(
    "1830-12-01", "1889-12-02", "1918-12-01",
    "1957-12-01", "1968-12-01", "2009-12-01"
  ))
)

# Midpoints and interpandemic period lengths
Periods <- Periods %>%
  mutate(
    length_years = round(as.numeric(difftime(end, start, units = "days")) / 365.25, 1),
    mid_date     = start + (end - start) / 2,
    length_label = ifelse(length_years == 41, "41 years", as.character(length_years))
  )

Fig1 <- ggplot(Timeline, aes(x = when, y = 0)) +
  geom_line(color = "black") +
  geom_segment(
    data = Periods,
    aes(x = start, xend = end, y = 0, yend = 0),
    linewidth = 2, color = "black"
  ) +
  geom_segment(
    data = Timeline,
    aes(x = when, xend = when, y = 0, yend = 0.001),
    linewidth = 0.7, color = "black"
  ) +
  geom_text(
    data = Timeline,
    aes(x = when, y = 0.002, label = what),
    angle = 45, hjust = 0, size = 4) +
  geom_text(
    data = Periods,
    aes(x = mid_date, y = -0.001, label = length_label),
    size = 4, fontface = "plain"
  ) +
  geom_point(aes(y = 0), color = "black") +
  scale_x_date(
    name = "Year",
    limits = c(ymd("1781-01-01"), ymd("2009-12-01")),
    expand = expansion(mult = c(0.05, 0.1))
  ) +
  scale_y_continuous(
    name = "", limits = c(-0.003, 0.03), breaks = NULL, expand = expansion(mult = c(0, 0))
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14, face = "plain"),
    axis.text.x = element_text(size = 14, face = "plain", margin = margin(t = 5)),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

ggsave("Fig1.tiff", Fig1, width = 10, height = 6, dpi = 300, device = "tiff", compression = "lzw")

# ============================================================
# Plot interpandemic period distributions: FigA2
# ============================================================

period <- seq(0, 100)
psi.fit <- data.frame(period, prob = dgamma(period, shape = alpha0, scale = theta0))

FigA2 <- ggplot(psi.fit, aes(period, prob)) +
  geom_ribbon(aes(ymin = 0, ymax = prob), fill = "steelblue", color = "black") +
  annotate("point", x = 11, y = 0, color = "black", size = 3) +
  annotate("point", x = 49, y = 0, color = "black", size = 3) +
  annotate("point", x = 59, y = 0, color = "black", size = 3) +
  annotate("point", x = 29, y = 0, color = "black", size = 3) +
  annotate("point", x = 39, y = 0, color = "black", size = 3) +
  annotate("point", x = 41, y = 0, color = "black", size = 3) +
  xlab("Interpandemic period, *ψ* (years)") +
  ylab("Probability") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 16, face = "plain"),
    axis.title = element_text(size = 18, face = "plain"),
    legend.position = "none"
  ) +
  theme(axis.title.x = element_markdown())

ggsave("FigA2.tiff", FigA2, width = 10, height = 4, dpi = 300, device = "tiff", compression = "lzw")

# ============================================================
# Plot the effect of date on the interpandemic period: FigA3
# ============================================================

plot_data <- data.frame(
  dates = dates,
  dd    = dd,
  mean  = delta0 + delta1 * year
)

FigA3 <- ggplot(plot_data, aes(x = dates, y = dd)) +
  geom_point(size = 3, color = "black") +
  geom_line(aes(y = mean), color = "black", linewidth = 1.2) +
  xlab("Year of pandemic start") +
  ylab("Interpandemic period, *Ψ*") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 16, face = "plain"),
    axis.title = element_text(size = 18, face = "plain"),
    axis.title.y = element_markdown()
  )

ggsave("FigA3.tiff", FigA3, width = 10, height = 6, dpi = 300, device = "tiff", compression = "lzw")
