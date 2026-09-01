# ============================================================
# Estimating the infection fatality ratio for humans infected with avian influenza viruses
# ============================================================

# Fixed random seed 
set.seed(20260828)

# ============================================================
# Load libraries
# ============================================================

library(patchwork)
library(dplyr)
library(ggplot2)
library(ggtext)
library(scales)
library(ggrepel)
library(grid)
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
# pi_i function
# ============================================================

model <- function(parameters) {
  a     <- parameters[1]
  Rstar <- parameters[2]
  R0    <- parameters[3]
  psi   <- parameters[4]
  
  a1 <- (R0 - 1 - a)
  a2 <- (R0 + 1 + a)^2
  P  <- (1 - (1 / Rstar))
  a3 <- 1 + (a * (1 - P))
  
  (a1 + sqrt(a2 - 4 * R0 * a3)) / (2 * R0)
}

# ============================================================
# Run a single scenario
# ============================================================

run_scenario <- function(param_ranges, num_samples) {
  
  # Latin Hypercube Sampling
  lhs_samples <- randomLHS(num_samples, length(param_ranges))
  
  transformed_samples <- data.frame(
    a = qunif(lhs_samples[,1],
              min = param_ranges$a[1],
              max = param_ranges$a[2]),
    Rstar = qtriangle(lhs_samples[,2],
                      a = param_ranges$Rstar[1],
                      b = param_ranges$Rstar[2],
                      c = param_ranges$Rstar[3]),
    R0 = qexp(lhs_samples[,3],
              rate = param_ranges$R0$rate),
    psi = qgamma(lhs_samples[,4],
                 shape = param_ranges$psi$shape,
                 scale = param_ranges$psi$scale)
  )
  
  # Calculate pi_i
  transformed_samples$pi <- apply(transformed_samples, 1, model)
  bar_pi <- mean(transformed_samples$pi)
  
  # Solve for mean number of zoonotic spillovers
  Lambda <- function(mean_nz, Lambda_target) {
    nz <- seq(0, 1e5)
    PNz_nz <- dpois(nz, mean_nz)
    1 - sum(PNz_nz * (1 - bar_pi)^nz) - Lambda_target
  }
  
  Lambda_target_vec <- 1 / transformed_samples$psi
  mean_nz_vec <- numeric(length(Lambda_target_vec))
  
  for (i in seq_along(Lambda_target_vec)) {
    mean_nz_vec[i] <- uniroot(Lambda, interval = c(0, 1e6), Lambda_target = Lambda_target_vec[i])$root
  }
  
  transformed_samples$mean_nz <- mean_nz_vec
  
  # Mean number of human infections
  bar_m <- mean(1 / (1 - transformed_samples$R0))
  transformed_samples$meannh <- bar_m * transformed_samples$mean_nz
  
  # Sort and calculate cumulative probabilities
  cum.prob <- seq(1 / num_samples, 1, 1 / num_samples)
  sort.output <- data.frame(
    meannh = sort(transformed_samples$meannh),
    cum.prob = cum.prob
  )
  
  # Quantiles
  i025 <- min(which(sort.output$cum.prob >= 0.025))
  i25  <- min(which(sort.output$cum.prob >= 0.25))
  i50  <- min(which(sort.output$cum.prob >= 0.5))
  i75  <- min(which(sort.output$cum.prob >= 0.75))
  i975 <- min(which(sort.output$cum.prob >= 0.975))
  
  quarts <- data.frame(
    cprob = c(sort.output$cum.prob[i025], sort.output$cum.prob[i25],
              sort.output$cum.prob[i50], sort.output$cum.prob[i75],
              sort.output$cum.prob[i975]),
    n = c(sort.output$meannh[i025], sort.output$meannh[i25],
          sort.output$meannh[i50], sort.output$meannh[i75],
          sort.output$meannh[i975])
  )
  
  list(
    transformed_samples = transformed_samples,
    sort_output = sort.output,
    quarts = quarts
  )
}

# ============================================================
# Run scenarios
# ============================================================

num_samples <- 1000

res_1 <- run_scenario(
  list(
    R0 = list(rate = 20),
    a = c(0.000012, 0.000024),
    Rstar = c(1, 2, 1.1),
    psi = list(shape = delta/theta0, scale = theta0)
  ),
  num_samples
)

res_2 <- run_scenario(
  list(
    R0 = list(rate = 20),
    a = c(0.000012, 0.000024),
    Rstar = c(1, 2, 1.1),
    psi = list(shape = delta0/theta, scale = theta)
  ),
  num_samples
)

res_3 <- run_scenario(
  list(
    R0 = list(rate = 20),
    a = c(0.000012, 0.000024),
    Rstar = c(1, 2, 1.1),
    psi = list(shape = delta0/theta, scale = 18.8/(delta0/theta))
  ),
  num_samples
)

# ============================================================
# Save results
# ============================================================

write.csv(res_1$transformed_samples, "transformed_samples_38.csv", row.names = FALSE, quote = FALSE)
write.csv(res_2$transformed_samples, "transformed_samples_53.5.csv", row.names = FALSE, quote = FALSE)
write.csv(res_3$transformed_samples, "transformed_samples_18.8.csv", row.names = FALSE, quote = FALSE)

write.csv(res_1$sort_output, "mean_nh_38.csv", row.names = FALSE, quote = FALSE)
write.csv(res_2$sort_output, "mean_nh_53.5.csv", row.names = FALSE, quote = FALSE)
write.csv(res_3$sort_output, "mean_nh_18.8.csv", row.names = FALSE, quote = FALSE)

# ============================================================
# Plot the mean annual number of human infections with AIV: Fig2A
# ============================================================

cols <- c(
  "53.5" = "#1f78b4",  # blue
  "38"   = "#e31a1c",  # red
  "18.8" = "#ff7f00"   # orange
)

base_theme <- theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    plot.tag = element_text(size = 12, face = "bold"),
    plot.tag.position = c(0,1)
  )

Fig2A <- ggplot() +
  geom_line(data = res_1$sort_output, aes(x = meannh, y = cum.prob, colour = "38"), size = 0.9) +
  geom_line(data = res_2$sort_output, aes(x = meannh, y = cum.prob, colour = "53.5"), size = 0.9) +
  geom_line(data = res_3$sort_output, aes(x = meannh, y = cum.prob, colour = "18.8"), size = 0.9) +
  
  geom_point(data = res_1$quarts, aes(x = n, y = cprob), size = 1.5) +
  geom_point(data = res_2$quarts, aes(x = n, y = cprob), size = 1.5) +
  geom_point(data = res_3$quarts, aes(x = n, y = cprob), size = 1.5) +
  
  geom_text_repel(data = res_1$quarts, aes(x = n, y = cprob, label = round(n,0)),
                  size = 3.5, hjust = 0, vjust = 1) +
  geom_text_repel(data = res_2$quarts, aes(x = n, y = cprob, label = round(n,0)),
                  size = 3.5, hjust = 0, vjust = 1) +
  geom_text_repel(data = res_3$quarts, aes(x = n, y = cprob, label = round(n,0)),
                  size = 3.5, hjust = 1, vjust = 0) +
  
  ylab("Cumulative probability") +
  
  xlab(expression(Mean~number~of~human~infections*","~bar(italic(n))[h])) +
  
  scale_x_log10(
    breaks = c(1000, 10000, 100000),
    labels = c("1,000","10,000","100,000"),
    limits = c(1000, 100000)
  ) +
  
  scale_color_manual(
    values = cols,
    breaks = c("18.8","38","53.5"),
    labels = c("18.8","38","53.5")
  ) +
  
  labs(color = "Mean interpandemic\nperiod (years)",
       tag = "A") +
  
  base_theme +
  theme(
    legend.position = c(0.78,0.35),
    legend.key.size = unit(0.6,"cm")
  )

# ============================================================
# Prepare IFR distributions
# ============================================================

deaths <- 20.6

make_ifr_df <- function(res){
  df <- res$sort_output %>%
    mutate(IFR = deaths / meannh * 10000)
  
  data.frame(
    IFR = sort(df$IFR),
    cum.prob = seq(1 / nrow(df), 1, 1 / nrow(df))
  )
}

ifr_1 <- make_ifr_df(res_1)
ifr_2 <- make_ifr_df(res_2)
ifr_3 <- make_ifr_df(res_3)

# ----------------------------
# Quantiles 
# ----------------------------

get_quarts <- function(ifr_df){
  probs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  idx <- sapply(probs, function(p)
    min(which(ifr_df$cum.prob >= p))
  )
  
  data.frame(
    cprob = ifr_df$cum.prob[idx],
    IFR   = ifr_df$IFR[idx]
  )
}

q1 <- get_quarts(ifr_1)
q2 <- get_quarts(ifr_2)
q3 <- get_quarts(ifr_3)

# ============================================================
# Plot the IFRs: Fig2B
# ============================================================

Fig2B <- ggplot() +
  geom_line(data = ifr_1, aes(x = IFR, y = cum.prob, colour = "38"), size = 0.9) +
  geom_line(data = ifr_2, aes(x = IFR, y = cum.prob, colour = "53.5"), size = 0.9) +
  geom_line(data = ifr_3, aes(x = IFR, y = cum.prob, colour = "18.8"), size = 0.9) +
  
  geom_point(data = q1, aes(x = IFR, y = cprob), size = 1.5) +
  geom_point(data = q2, aes(x = IFR, y = cprob), size = 1.5) +
  geom_point(data = q3, aes(x = IFR, y = cprob), size = 1.5) +
  
  geom_text_repel(data = q1, aes(x = IFR, y = cprob, label = round(IFR,1)),
                  size = 3.5, hjust = 0, vjust = 1) +
  geom_text_repel(data = q2, aes(x = IFR, y = cprob, label = round(IFR,1)),
                  size = 3.5, hjust = 0, vjust = 1) +
  geom_text_repel(data = q3, aes(x = IFR, y = cprob, label = round(IFR,1)),
                  size = 3.5, hjust = 1, vjust = 0) +
  
  ylab("Cumulative probability") +
  xlab("Infection fatality ratio (per 10,000 infections)") +
  
  scale_x_log10(
    breaks = c(1,10,100,1000)
  ) +
  
  scale_color_manual(values = cols) +
  
  labs(tag = "B") +
  
  base_theme +
  theme(legend.position = "none")

# ============================================================
# Combine Fig2A and Fig2B 
# ============================================================
Fig2 <- Fig2A + Fig2B +
  plot_layout(ncol = 1)

ggsave("Fig2.tiff",
       Fig2,
       width = 7,
       height = 8,
       dpi = 300,
       device = "tiff",
       compression = "lzw")

# ============================================================
# Plot IFR comparison: Fig3
# ============================================================

# IFR data
ifr_data <- data.frame(
  Virus = c("Seasonal influenza",
            "SARS-CoV-2",
            "AIV (53.5 years)",
            "AIV (38 years)",
            "AIV (18.8 years)"),
  estimate = c(4.7, 16, 46, 32, 16),
  lower = c(2.9, 16, 17, 9.6, 5.6),
  upper = c(6.5, 16, 97, 75, 33)
)

# Order for plotting
ifr_data$Virus <- factor(
  ifr_data$Virus,
  levels = rev(c("AIV (53.5 years)",
                 "AIV (38 years)",
                 "AIV (18.8 years)",
                 "SARS-CoV-2",
                 "Seasonal influenza"))
)

base_theme <- theme_bw() +
  theme(
    axis.text = element_text(size = 14, face = "plain", color = "black"),   # bigger tick labels
    axis.title = element_text(size = 14, face = "plain", color = "black"), # bigger axis titles
    plot.tag = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

Fig3 <- ggplot(ifr_data, aes(x = estimate, y = Virus)) +
  geom_errorbarh(aes(
    xmin = ifelse(Virus != "SARS-CoV-2", lower, NA),
    xmax = ifelse(Virus != "SARS-CoV-2", upper, NA)
  ),
  height = 0.2,
  size = 1.2) +
  geom_point(size = 2.5, color = "red") +
  geom_text_repel(aes(label = estimate),
                  size = 4.5,       # slightly larger labels for readability
                  fontface = "plain",
                  nudge_x = 0.5,
                  nudge_y = 0.20,
                  direction = "y",
                  segment.size = 0.2) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  xlab("Infection fatality ratio (per 10,000 infections)") +
  ylab(NULL) +
  base_theme

ggsave(
  filename = "Fig3.tiff",
  plot = Fig3,
  width = 8,
  height = 5,
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)

# ============================================================
# Plot parameter distributions: FigA1
# ============================================================

plot_distributions <- function(transformed_samples) {
  
  plot_a <- ggplot(transformed_samples, aes(x = a)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "black") +
    labs(x = expression(italic(a)), y = "Count") +
    ggtitle("Probability that a strain of AIV with pandemic potential emerges") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18),     
      axis.text = element_text(size = 16),                    
      axis.title = element_text(size = 18),                   
      plot.margin = margin(t = 15, r = 10, b = 15, l = 10)
    )
  
  plot_Rstar <- ggplot(transformed_samples, aes(x = Rstar)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "black") +
    labs(x = expression(italic(R)^"*"), y = "Count") +
    ggtitle("Reproduction number of a reassortant AIV in humans") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 18),
      plot.margin = margin(t = 15, r = 10, b = 15, l = 10)
    )
  
  plot_R0 <- ggplot(transformed_samples, aes(x = R0)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "black") +
    labs(x = expression(italic(R)[0]), y = "Count") +
    ggtitle("Reproduction number of an AIV in humans prior to evolutionary change") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 18),
      plot.margin = margin(t = 15, r = 10, b = 20, l = 10)
    )
  
  combined <- (plot_a / plot_Rstar / plot_R0) +
    plot_annotation(tag_levels = 'A') &       
    theme(
      plot.tag = element_text(face = "bold", size = 18)  
    )
  
  return(combined)
}

FigA1 <- plot_distributions(res_1$transformed_samples)

ggsave(
  filename = "FigA1.tiff",
  plot = FigA1,
  device = "tiff",
  height = 10,
  width = 12,
  dpi = 300,
  compression = "lzw"
)
