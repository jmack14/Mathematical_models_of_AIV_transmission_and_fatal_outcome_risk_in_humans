# ============================================================
# Comparison of mean annual number of zoonotic spillovers of AIV with Day et al., 2006
# ============================================================

# ============================================================
# Load libraries
# ============================================================

library(patchwork)
library(dplyr)
library(ggplot2)
library(ggridges)
library(ggtext)
library(lhs)
library(sensitivity)
library(triangle)
library(here)

# ============================================================
# Load interpandemic period results
# ============================================================

# Path relative to project root
Interpandemic_results <- here("Code", "Clean_code", "Output", "Interpandemic_period_results.RData")

# Check if file exists before loading
if (!file.exists(Interpandemic_results)) {
  stop("Error: Interpandemic_period_results.RData not found.
       Make sure the Output folder is inside Code/Clean_code/")
}

# Load the file
load(Interpandemic_results)
message("Interpandemic period results loaded successfully!") 

# ============================================================
# Helper functions
# ============================================================

# Core spillover model
spillover_model <- function(parameters, R0) {
  psi   <- parameters[1]
  a     <- parameters[2]
  Rstar <- parameters[3]
  (1 - R0) / (a * psi * (1 - (1 / Rstar)))
}

# Wrapper to run Latin hypercube sampling
run_lhs <- function(R0) {
  
  param_ranges <- list(
    psi   = list(shape = delta/theta0, scale = theta0),
    a     = c(0.000012, 0.000024),
    Rstar = c(1, 2, 1.1)
  )
  
  lhs_samples <- randomLHS(10000, length(param_ranges))
  
  transformed_samples <- data.frame(
    psi   = qgamma(lhs_samples[,1], shape = param_ranges$psi$shape, scale = param_ranges$psi$scale),
    a     = qunif(lhs_samples[,2], min = param_ranges$a[1], max = param_ranges$a[2]),
    Rstar = qtriangle(lhs_samples[,3],
                      a = param_ranges$Rstar[1],
                      b = param_ranges$Rstar[2],
                      c = param_ranges$Rstar[3])
  )
  
  output <- apply(transformed_samples, 1, spillover_model, R0 = R0)
  cum.prob <- seq(1/length(output), 1, 1/length(output))
  
  data.frame(
    n        = sort(output),
    cum.prob = cum.prob,
    R0       = R0
  )
}

# ============================================================
# Run all R0 scenarios and combine results
# ============================================================

sort.output  <- run_lhs(0.0) %>% mutate(n1=292, n2=583,  n3=3167, n4=6333)
sort.output2 <- run_lhs(0.2) %>% mutate(n1=233, n2=467,  n3=2533, n4=5067)
sort.output4 <- run_lhs(0.4) %>% mutate(n1=175, n2=350,  n3=1900, n4=3800)
sort.output6 <- run_lhs(0.6) %>% mutate(n1=117, n2=233,  n3=1267, n4=2533)
sort.output8 <- run_lhs(0.8) %>% mutate(n1=58,  n2=117,  n3=633,  n4=1267)

ridge <- bind_rows(
  sort.output,
  sort.output2,
  sort.output4,
  sort.output6,
  sort.output8
)

# ============================================================
# Shading
# ============================================================

light.grey <- data.frame(
  R0 = c(seq(0,1,.2), seq(.8,0,-.2), 0),
  n  = c(292,233,175,117,58,0,633,1267,1900,2533,3167,292)
)

dark.grey <- data.frame(
  R0 = c(seq(0,1,.2), seq(.8,0,-.2), 0),
  n  = c(583,467,350,233,117,0,1267,2533,3800,5067,6333,583)
)

# ============================================================
# Plot minimum annual spillovers: Fig5A
# ============================================================

ridges_min <- ggplot() +
  geom_polygon(data = dark.grey, aes(x = n, y = R0),
               fill = "darkgrey", color = "black", lwd = 0.2) +
  geom_polygon(data = light.grey, aes(x = n, y = R0),
               fill = "lightgrey", color = "black", lwd = 0.2, alpha = 0.5) +
  geom_density_ridges_gradient(
    data = ridgeP,
    aes(x = n, y = R0, group = R0, fill = after_stat(x)),
    scale = 1,
    rel_min_height = 0.01
  ) +
  scale_x_continuous(expand = c(0,0), limits = c(0,7000)) +
  scale_fill_viridis_c(
    name = expression(italic(bar(n))[italic(z)]),
    option = "C"
  ) +
  xlab(expression("Mean number of spillovers, " * italic(bar(n))[italic(z)])) +
  ylab(expression(italic(R[0]))) +
  ggtitle("Minimum annual number of zoonotic spillovers") +
  labs(tag = "A") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0, size = 18, face = "plain"),
    axis.text  = element_text(size = 16, face = "plain"),
    axis.title = element_text(size = 18, face = "plain"),
    plot.tag   = element_text(size = 18, face = "bold"),
    plot.tag.position = c(0,1),
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 16),
    legend.position = "none"  
  )

# ============================================================
# Plot annual spillovers: Fig5B
# ============================================================

ridges <- ggplot() +
  geom_polygon(data = dark.grey, aes(x = n, y = R0),
               fill = "darkgrey", color = "black", lwd = 0.2) +
  geom_polygon(data = light.grey, aes(x = n, y = R0),
               fill = "lightgrey", color = "black", lwd = 0.2, alpha = 0.5) +
  geom_density_ridges_gradient(
    data = ridge,
    aes(x = n, y = R0, group = R0, fill = after_stat(x)),
    scale = 1,
    rel_min_height = 0.01
  ) +
  scale_x_continuous(expand = c(0,0), limits = c(0,7000)) +
  scale_fill_viridis_c(
    name = expression(italic(bar(n))[italic(z)]),
    option = "C"
  ) +
  xlab(expression("Mean number of spillovers, " * italic(bar(n))[italic(z)])) +
  ylab(expression(italic(R[0]))) +
  ggtitle("Annual number of zoonotic spillovers") +
  labs(tag = "B") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0, size = 18, face = "plain"),
    axis.text  = element_text(size = 16, face = "plain"),
    axis.title = element_text(size = 18, face = "plain"),
    legend.title = element_text(size = 18),
    legend.text  = element_text(size = 16),
    plot.tag = element_text(size = 18, face = "bold"),
    plot.tag.position = c(0,1)
  )

# ============================================================
# Combine Fig5A and Fig5B
# ============================================================

Fig5 <- ridges_min + ridges

ggsave(
  filename = "Fig5.tiff",
  plot = Fig5,
  width = 14,
  height = 7,
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)