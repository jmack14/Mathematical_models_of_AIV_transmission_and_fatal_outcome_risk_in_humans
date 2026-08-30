# ============================================================
# Mechanistic boundary analysis 
# ============================================================

# Function to calculate pi
calc_pi <- function(a, Rstar, R0){
  
  a1 <- (R0 - 1 - a)
  a2 <- (R0 + 1 + a)^2
  
  P <- (1 - (1/Rstar))
  
  a3 <- 1 + (a*(1-P))
  
  (a1 + sqrt(a2 - 4*R0*a3))/(2*R0)
}


# ------------------------------------------------------------
# Extreme values of pi
# ------------------------------------------------------------

R0_min <- qexp(0.025, rate = 20)
R0_max <- qexp(0.975, rate = 20)

a_min <- 0.000012
a_max <- 0.000024

Rstar_min <- 1.1
Rstar_max <- 2


# Calculate pi across extreme combinations
pi_values <- expand.grid(
  a = c(a_min, a_max),
  Rstar = c(Rstar_min, Rstar_max),
  R0 = c(R0_min, R0_max)
)

pi_values$pi <- mapply(
  calc_pi,
  pi_values$a,
  pi_values$Rstar,
  pi_values$R0
)


pi_min <- min(pi_values$pi)
pi_max <- max(pi_values$pi)


# ------------------------------------------------------------
# Convert pi into required zoonotic spillovers
# ------------------------------------------------------------

pandemic_prob <- 1/53.5


solve_nz <- function(pi, pandemic_prob){
  
  f <- function(nz){
    1 - (1-pi)^nz - pandemic_prob
  }
  
  uniroot(
    f,
    interval = c(0, 1e8)
  )$root
}


# Lower pi requires more spillovers
nz_max <- solve_nz(pi_min, pandemic_prob)

# Higher pi requires fewer spillovers
nz_min <- solve_nz(pi_max, pandemic_prob)


# ------------------------------------------------------------
# Extreme values of m
# ------------------------------------------------------------

m_min <- 1/(1-R0_min)
m_max <- 1/(1-R0_max)


# ------------------------------------------------------------
# Extreme annual human infection estimates
# ------------------------------------------------------------

# Minimum infections:
# fewest spillovers C least transmission after spillover
nh_min <- nz_min * m_min

# Maximum infections:
# most spillovers C greatest transmission after spillover
nh_max <- nz_max * m_max


# ------------------------------------------------------------
# Extreme IFR estimates
# ------------------------------------------------------------

annual_deaths <- 20.6

# IFR highest when infections lowest
ifr_max <- annual_deaths / nh_min * 100

# IFR lowest when infections highest
ifr_min <- annual_deaths / nh_max * 100


# ------------------------------------------------------------
# Boundary results
# ------------------------------------------------------------

boundary_results <- data.frame(
  
  Parameter = c(
    "Minimum pi",
    "Maximum pi",
    "Minimum zoonotic spillovers",
    "Maximum zoonotic spillovers",
    "Minimum m",
    "Maximum m",
    "Minimum annual infections",
    "Maximum annual infections",
    "Minimum IFR (%)",
    "Maximum IFR (%)"
  ),
  
  Value = c(
    pi_min,
    pi_max,
    nz_min,
    nz_max,
    m_min,
    m_max,
    nh_min,
    nh_max,
    ifr_min,
    ifr_max
  )
)


print(boundary_results)


write.csv(
  boundary_results,
  "Mechanistic_boundary.csv",
  row.names = FALSE
)