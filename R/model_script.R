# R/model_script.R
#
# Run from the project root:
# source("R/model_script.R")
#
# Generates the symbolic model.
# Simulation is performed separately in R/dmod_model.R.


source("R/equation_generator.R")


# ------------------------------------------------------------
# 1. Reference geometry
# ------------------------------------------------------------

bounds <- list(
  min = c(1, 1),
  max = c(7, 5)
)

spacing <- c(0.5, 0.8)
origin <- c(0, 0)

# Out-of-plane thickness for conversion:
# sarcomere force -> nodal force density.
thickness <- 1


# ------------------------------------------------------------
# 2. Material, density and pacing parameters
#
# Illustrative values, not calibrated physical parameters.
# All quantities must use consistent units.
# ------------------------------------------------------------

parameters <- list(
  # Spatially varying positive mass density.
  rho = function(X, parameters) {
    1 + 0.1 * X[1]
  },
  
  # Isotropic example materials.
  # General anisotropic tensors can be supplied through E and D.
  lambda_E = 1,
  mu_E = 1,
  
  lambda_D = 0.1,
  mu_D = 0.1,
  
  kappa = 0.05,
  
  # Pacing: 1 Hz, active during one third of each cycle.
  pace_period = 1,
  pace_duty = 1 / 3,
  pace_start = 0,
  
  # Peak force per sarcomere.
  Tmax = 0.1,
  
  # Relative strength of one selected sarcomere.
  weak_fraction = 0.2
)


# ------------------------------------------------------------
# 3. Active sarcomere forces
#
# Weak bond: node (3,3) -> node (4,3).
# Its left endpoint has reference position (1.0, 1.6).
# ------------------------------------------------------------

sarcomere_tension <- function(X_left, X_right, parameters) {
  Tmax <- parameter_value(parameters, "Tmax")
  weak_fraction <- parameter_value(parameters, "weak_fraction")
  
  check_scalar(Tmax, "Tmax")
  check_scalar(weak_fraction, "weak_fraction")
  
  if (Tmax < 0 ||
      weak_fraction < 0 ||
      weak_fraction > 1) {
    stop("Require Tmax >= 0 and 0 <= weak_fraction <= 1.")
  }
  
  weak <- abs(X_left[1] - 1.0) < 1e-10 &&
    abs(X_left[2] - 1.6) < 1e-10
  
  Tmax * (if (weak) weak_fraction else 1)
}

parameters$P_amplitude <- make_sarcomere_force_profile(
  grid_bounds = bounds,
  spacing = spacing,
  tension = sarcomere_tension,
  origin = origin,
  thickness = thickness
)


# ------------------------------------------------------------
# 4. Generate the full model
#
# du/dt = v
# dv/dt = [div(E:grad(u) + D:grad(v)) - kappa*u + P] / rho
#
# Left boundary fixed; other boundaries free.
# ------------------------------------------------------------

model <- generate_equations(
  grid_bounds = bounds,
  spacing = spacing,
  fields = default_fields(2),
  parameters = parameters,
  origin = origin,
  fixed_side = "left"
)


# ------------------------------------------------------------
# 5. Basic checks
# ------------------------------------------------------------

stopifnot(
  nrow(model$system) == 120L,
  length(model$fixed_values) == 20L,
  all(model$fixed_values == 0),
  all(is.finite(model$parameter_values)),
  all(model$initial_values == 0)
)

# All nodal density values must be positive.
rho_names <- grep(
  "^rho_",
  names(model$parameter_values),
  value = TRUE
)

stopifnot(
  length(rho_names) > 0L,
  all(model$parameter_values[rho_names] > 0)
)

# Pacing peak and inactive phase.
stopifnot(
  abs(pacing_pulse(1 / 6, parameters) - 1) < 1e-12,
  pacing_pulse(0.6, parameters) == 0
)

cat(
  "Model generated:",
  nrow(model$system), "dynamic states,",
  length(model$fixed_values), "fixed components.\n"
)

print(head(model$system[, c("state", "boundary")]))