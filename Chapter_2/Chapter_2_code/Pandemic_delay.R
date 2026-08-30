# ============================================================
# Estimating the effect of preventing AIV spillovers on pandemic delay time
# ============================================================

# ============================================================
# Load libraries
# ============================================================

library(dplyr)
library(ggplot2)
library(here)

# ============================================================
# Read transformed samples
# ============================================================

# Path to transformed samples relative to project root
transformed_samples <- read.csv(
  "C:/Users/marti/Desktop/Joshua/Code/Clean_code/Output/transformed_samples_38.csv",
  stringsAsFactors = FALSE
)


<- here("Code", "Clean_code", "Output", "transformed_samples_38.csv")

# Check if file exists before reading
if (!file.exists(transformed_file)) {
  stop("Error: transformed_samples_38.csv not found.
       Make sure the Output folder is inside the project root.")
}

transformed_samples <- read.csv(transformed_file, stringsAsFactors = FALSE) 
message("Transformed samples loaded successfully!")

# ============================================================
# Calculate baseline mean pi
# ============================================================

bar_pi <- mean(transformed_samples$pi)

# Vector of fraction of zoonotic spillovers averted
Ovec <- seq(0, 0.5, 0.05)

# ============================================================
# Lambda function
# ============================================================

Lambda <- function(mean_nz) {
  nz <- seq(1, 1e5)
  PNz_nz <- dpois(nz, mean_nz * (1 - O))
  Lambda_val <- 1 - sum(PNz_nz * (1 - bar_pi)^nz)
  return(Lambda_val)
}

# ============================================================
# Compute Lambda and delay for all scenarios
# ============================================================

Lambda.vals <- NULL

for (i in seq_along(Ovec)) {
  O <- Ovec[i]
  for (j in seq_along(transformed_samples$mean_nz)) {
    mean_nz <- transformed_samples$mean_nz[j]
    Lambda1 <- Lambda(mean_nz)
    delay <- 1 / Lambda1
    Lambda.vals <- rbind(Lambda.vals, 
                         data.frame(O = O, mean_nz = mean_nz, Lamda = Lambda1, delay = delay))
  }
}

# ============================================================
# Compute baseline interpandemic period
# ============================================================

Lambda.vals_0 <- filter(Lambda.vals, O == 0)
mean.inter.pan <- mean(Lambda.vals_0$delay)

# ============================================================
# Adjust delays and convert O to percentage
# ============================================================

Lambda.vals <- Lambda.vals %>%
  mutate(delay0 = delay - mean.inter.pan,
         O = 100 * O)

# ============================================================
# Summarize mean and quantiles for plotting
# ============================================================

mean.Lambda <- Lambda.vals %>%
  group_by(O) %>%
  summarise(
    mean.delay0 = mean(delay0),
    max.delay0 = quantile(delay0, probs = 0.025),
    min.delay0 = quantile(delay0, probs = 0.975),
    lowerq = quantile(delay0, 0.25),
    upperq = quantile(delay0, 0.75),
    .groups = "drop"
  ) %>%
  mutate(O = round(O, 0))

# Select subset for labeling
subset <- filter(mean.Lambda, O %in% c(5, 10, 20, 30, 40, 50)) %>%
  mutate(rounded = round(mean.delay0, 1))

# ============================================================
# Plot the effect of the number of zoonotic spillovers prevented on pandemic delay time: Fig4
# ============================================================

# Prepare labels for selected points
subset <- filter(mean.Lambda, O %in% c(5, 10, 20, 30, 40, 50)) %>%
  mutate(rounded = round(mean.delay0, 1))

Fig4 <- ggplot(data = mean.Lambda, aes(x = O)) +
  geom_ribbon(aes(ymin = min.delay0, ymax = max.delay0), fill = "dodgerblue", alpha = 0.25) +
  geom_ribbon(aes(ymin = lowerq, ymax = upperq), fill = "dodgerblue", alpha = 0.75) +
  geom_line(aes(y = mean.delay0), size = 1) +
  geom_point(data = subset, aes(x = O, y = mean.delay0), size = 3, color = "black") +
  geom_text(data = subset, aes(x = O, y = mean.delay0, label = rounded),
            hjust = 0.8, vjust = -1, size = 4, fontface = "plain") +
  xlab("% Zoonotic spillovers to humans prevented") +
  ylab("Delay to next pandemic (years)") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 14),       
    axis.title = element_text(size = 14, face = "plain"),  
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  )

ggsave(
  filename = "Fig4.tiff",
  plot = Fig4,
  device = "tiff",
  width = 8,      
  height = 5,     
  dpi = 300,
  compression = "lzw"
)
