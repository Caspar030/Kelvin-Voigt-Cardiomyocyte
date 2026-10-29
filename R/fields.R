# R/fields.R
#

#         ###
#        ####
#       ## ##
#      ##  ##
#     ##   ##
#    ##    ##
#   ##########
#          ##
#          ##
#          ##


#In fields.R we define the symbolic fields of operators.R numerically.





# ============================================================
# FILE OVERVIEW — fields.R
#
# Numerical geometry, material properties, and active forces
# for the Kelvin-Voigt model.
# ============================================================
#
# 1. Helpers and validation
#    Retrieve parameters, supply defaults, and validate inputs.
#
# 2. Reference and current geometry
#    Calculate reference distances, displaced positions,
#    and distances between displaced nodes.
#
# 3. Material tensors
#    Construct isotropic fourth-order tensors for E and D.
#    Validate tensor dimensions and finite numerical values.
#
# 4. Pacing
#    Define the time-dependent activation pulse and its
#    symbolic expression.
#
# 5. Active sarcomere forces
#    Build equal-and-opposite tensile forces between
#    longitudinal neighbors and convert net nodal forces
#    into force densities using control volumes.
#
# 6. Numerical field definitions
#    Provide a common interface for material properties,
#    density, restoring coefficients, and active forces.
#    Support constant values and spatially varying functions.
#
# 7. Numerical-to-symbolic field binding
#    Evaluate spatial fields at reference-grid nodes and
#    associate their values with named ODE coefficients.
#    Keep the activation pulse symbolic.
#
# ------------------------------------------------------------
# This file defines functions only.
# No model generation, simulation, or plotting is run here.
# ============================================================



# Reference configuration; vector-valued displacement.
#
# E, D:
#   fourth-order arrays or functions(X, parameters)
#
# rho, kappa, ceta:
#   scalars or functions(X, parameters)
#
# Active force density:
#   P(X,t) = P_amplitude(X) * pulse(t)
#
# make_sarcomere_force_profile() constructs P_amplitude from
# equal and opposite tensile forces between longitudinal neighbors.
#
# No source() calls or simulations in this file.


# ------------------------------------------------------------
# 1. Helpers and validation
# ------------------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

parameter_value <- function(parameters, name, default = NULL) {
  if (!is.list(parameters)) {
    stop("parameters must be a list.")
  }
  
  parameters[[name, exact = TRUE]] %||% default
}

check_vector <- function(x, name) {
  if (!is.numeric(x) ||
      !is.null(dim(x)) ||
      length(x) == 0L ||
      any(!is.finite(x))) {
    stop(name, " must be a finite numeric vector.")
  }
  
  invisible(TRUE)
}

check_scalar <- function(x, name) {
  if (!is.numeric(x) ||
      length(x) != 1L ||
      !is.null(dim(x)) ||
      !is.finite(x)) {
    stop(name, " must be a finite numeric scalar.")
  }
  
  invisible(TRUE)
}

check_spatial_dimension <- function(dimension) {
  check_scalar(dimension, "dimension")
  
  if (!(dimension %in% 1:3)) {
    stop("dimension must be 1, 2 or 3.")
  }
  
  as.integer(dimension)
}

check_field_geometry <- function(grid_bounds, spacing, origin) {
  check_vector(spacing, "spacing")
  check_vector(origin, "origin")
  
  dimension <- check_spatial_dimension(length(spacing))
  
  if (any(spacing <= 0)) {
    stop("spacing must contain positive values.")
  }
  
  if (!is.list(grid_bounds)) {
    stop("grid_bounds must be a list.")
  }
  
  lower <- grid_bounds[["min", exact = TRUE]]
  upper <- grid_bounds[["max", exact = TRUE]]
  
  check_vector(lower, "grid_bounds$min")
  check_vector(upper, "grid_bounds$max")
  
  if (length(origin) != dimension ||
      length(lower) != dimension ||
      length(upper) != dimension ||
      any(lower != floor(lower)) ||
      any(upper != floor(upper)) ||
      any(lower >= upper)) {
    stop("Invalid geometry dimensions or grid bounds.")
  }
  
  invisible(dimension)
}


# ------------------------------------------------------------
# 2. Reference and current geometry
# ------------------------------------------------------------

# euclidean distance between two points in reference configuration.

A0 <- function(X, Y, parameters = list()) {
  check_vector(X, "X")
  check_vector(Y, "Y")
  
  if (length(X) != length(Y)) {
    stop("X and Y must have the same dimension.")
  }
  
  sqrt(sum((Y - X)^2))
}


current_position <- function(X, u) {
  check_vector(X, "X")
  check_vector(u, "u")
  
  if (length(X) != length(u)) {
    stop("X and u must have the same dimension.")
  }
  
  X + u
}

#current distance of nodes, AFTER their displacement from reference configuration
a_current <- function(X, Y, uX, uY, parameters = list()) {
  x <- current_position(X, uX)
  y <- current_position(Y, uY)
  
  if (length(x) != length(y)) {
    stop("Both nodes must have the same dimension.")
  }
  
  sqrt(sum((y - x)^2))
}


# ------------------------------------------------------------
# 3. Material tensors
# ------------------------------------------------------------


#Using Lame parameters lambda und mu um E zu konstruieren.
# They restrict your material to ISOTRPOIC LINEAR ELASTICITY.
# The array E can store more general tensors, but isotropic_tensor() generates only a particular subset of them.

# So far, this is the fallback option to when no explicit E is supplied: Happens in Section 6, inside default_fields() → material_field()


validate_material <- function(C, dimension, name) {
  expected <- rep(as.integer(dimension), 4L)
  
  if (!is.numeric(C) ||
      !identical(dim(C), expected) ||
      any(!is.finite(C))) {
    stop(
      name, " must be a finite numeric array with dimensions ",
      paste(expected, collapse = " x "), "."
    )
  }
  
  # Shape and finite values are checked.
  # Physical material symmetries and positivity are NOT checked.
  C
}

isotropic_tensor <- function(dimension, lambda, mu) {
  dimension <- check_spatial_dimension(dimension)
  check_scalar(lambda, "lambda")
  check_scalar(mu, "mu")
  
  C <- array(0, dim = rep(dimension, 4L))
  
  for (i in seq_len(dimension)) {
    for (j in seq_len(dimension)) {
      for (k in seq_len(dimension)) {
        for (l in seq_len(dimension)) {
          C[i, j, k, l] <-
            lambda * (i == j) * (k == l) +
            mu * (
              (i == k) * (j == l) +
                (i == l) * (j == k)
            )
        }
      }
    }
  }
  
  C
}


# ------------------------------------------------------------
# 4. Pacing
#
# Sinus-squared contraction pulse:
#   active during pace_duty * pace_period
#   zero during the remaining cycle
#
# Value and first time derivative are continuous.
# This is an activation pulse, not an electrical action potential.
# ------------------------------------------------------------

pacing_parameters <- function(parameters = list()) {
  period <- parameter_value(parameters, "pace_period", 1)
  duty <- parameter_value(parameters, "pace_duty", 1 / 3)
  start <- parameter_value(parameters, "pace_start", 0)
  
  check_scalar(period, "pace_period")
  check_scalar(duty, "pace_duty")
  check_scalar(start, "pace_start")
  
  if (period <= 0) {
    stop("pace_period must be positive.")
  }
  
  if (duty <= 0 || duty >= 1) {
    stop("pace_duty must lie strictly between 0 and 1.")
  }
  
  c(
    pace_period = unname(period),
    pace_duty = unname(duty),
    pace_start = unname(start)
  )
}

pacing_expression <- function() {
  cycles <- "((t - pace_start) / pace_period)"
  phase <- paste0("(", cycles, " - floor(", cycles, "))")
  
  paste0(
    "ifelse((t >= pace_start) & (", phase, " < pace_duty), ",
    "sin(pi * ", phase, " / pace_duty)^2, 0)"
  )
}

pacing_pulse <- function(t, parameters = list()) {
  check_vector(t, "t")
  
  values <- c(
    list(t = t),
    as.list(pacing_parameters(parameters))
  )
  
  eval(
    parse(text = pacing_expression()),
    envir = values,
    enclos = baseenv()
  )
}


# ------------------------------------------------------------
# 5. Active longitudinal sarcomere forces
#
# For a 2D grid, longitudinal direction is the first axis.
#
# tension(X_left, X_right, parameters):
#   returns a nonnegative peak tensile FORCE for one bond.
#   Must return the same value for identical inputs.
#
# Each bond pulls its endpoints toward each other.
# Forces act along the reference x-direction.
#
# The resulting node force is divided by the control volume
# to obtain a FORCE DENSITY, consistent with mass density rho.
#
# This does not prescribe displacement or shortening.
# There is currently no length/velocity dependence.
#
# thickness:
#   out-of-plane thickness in consistent physical units.
# ------------------------------------------------------------


#

make_sarcomere_force_profile <- function(
    grid_bounds,
    spacing,
    tension,
    origin = c(0, 0),
    thickness = 1
) {
  dimension <- check_field_geometry(
    grid_bounds, spacing, origin
  )
  
  if (dimension != 2L) {
    stop("The sarcomere force profile currently requires a 2D grid.")
  }
  
  if (!is.function(tension)) {
    stop("tension must be a function(X_left, X_right, parameters).")
  }
  
  check_scalar(thickness, "thickness")
  
  if (thickness <= 0) {
    stop("thickness must be positive.")
  }
  
  lower <- grid_bounds[["min", exact = TRUE]]
  upper <- grid_bounds[["max", exact = TRUE]]
  
  # Capture the configuration at profile creation.
  force(tension)
  force(spacing)
  force(origin)
  force(thickness)
  force(lower)
  force(upper)
  
  function(X, parameters = list()) {
    check_vector(X, "X")
    
    if (length(X) != 2L) {
      stop("X must have two components.")
    }
    
    grid_index <- lower + (X - origin) / spacing
    index <- round(grid_index)
    
    if (any(abs(grid_index - index) > 1e-8) ||
        any(index < lower) ||
        any(index > upper)) {
      stop("X must be a reference-grid node.")
    }
    
    bond_tension <- function(left_index) {
      right_index <- left_index + c(1, 0)
      
      X_left <- origin + (left_index - lower) * spacing
      X_right <- origin + (right_index - lower) * spacing
      
      value <- tension(X_left, X_right, parameters)
      check_scalar(value, "Sarcomere tension")
      
      if (value < 0) {
        stop("Active tensile force must be nonnegative.")
      }
      
      as.numeric(value)
    }
    
    # The right bond pulls the node toward +x.
    right_force <- if (index[1] < upper[1]) {
      bond_tension(index)
    } else {
      0
    }
    
    # The left bond pulls the node toward -x.
    left_force <- if (index[1] > lower[1]) {
      bond_tension(index - c(1, 0))
    } else {
      0
    }
    
    # Same control-volume convention as discretization.R.
    widths <- spacing
    on_boundary <- index == lower | index == upper
    widths[on_boundary] <- widths[on_boundary] / 2
    
    volume <- prod(widths) * thickness
    
    c(
      (right_force - left_force) / volume,
      0
    )
  }
}


# ------------------------------------------------------------
# 6. Numerical field definitions
#
# Example parameter entries:
#
# parameters$E = fourth-order array OR function(X, parameters)
# parameters$D = fourth-order array OR function(X, parameters)
# parameters$rho = scalar OR function(X, parameters)
#
# parameters$P_amplitude =
#   vector OR function(X, parameters)
#
# For sarcomere-based forces, assign the function returned by
# make_sarcomere_force_profile() to parameters$P_amplitude.
#
# Default active force is zero until a profile is supplied.
# ------------------------------------------------------------



# Here, numerical field definitions can be made.
default_fields <- function(dim) {
  dimension <- check_spatial_dimension(dim)
  
  check_position <- function(X) {
    check_vector(X, "X")
    
    if (length(X) != dimension) {
      stop("Position dimension does not match field dimension.")
    }
  }
  
  material_field <- function(name, lambda_default, mu_default) {
    force(name)
    force(lambda_default)
    force(mu_default)
    
    function(X, parameters = list()) {
      check_position(X)
      
      C <- parameter_value(parameters, name)
      
      if (is.null(C)) {
        C <- isotropic_tensor(
          dimension,
          lambda = parameter_value(
            parameters,
            paste0("lambda_", name),
            lambda_default
          ),
          mu = parameter_value(
            parameters,
            paste0("mu_", name),
            mu_default
          )
        )
      } else if (is.function(C)) {
        C <- C(X, parameters)
      }
      
      validate_material(C, dimension, name)
    }
  }
  
  scalar_field <- function(name, default,
                           positive = FALSE,
                           nonnegative = FALSE) {
    force(name)
    force(default)
    force(positive)
    force(nonnegative)
    
    function(X, parameters = list()) {
      check_position(X)
      
      value <- parameter_value(parameters, name, default)
      
      if (is.function(value)) {
        value <- value(X, parameters)
      }
      
      check_scalar(value, name)
      
      if (positive && value <= 0) {
        stop(name, " must be positive.")
      }
      
      if (nonnegative && value < 0) {
        stop(name, " must be nonnegative.")
      }
      
      as.numeric(value)
    }
  }
  
  amplitude <- function(X, parameters = list()) {
    check_position(X)
    
    # Exact lookup: P_amplitude is not mistaken for P.
    if (!is.null(parameter_value(parameters, "P"))) {
      stop(
        "Use parameters$P_amplitude for the spatial force profile. ",
        "Time dependence is specified by the pacing parameters."
      )
    }
    
    value <- parameter_value(
      parameters,
      "P_amplitude",
      rep(0, dimension)
    )
    
    if (is.function(value)) {
      value <- value(X, parameters)
    }
    
    check_vector(value, "P_amplitude")
    
    if (length(value) != dimension) {
      stop("P_amplitude needs one component per spatial dimension.")
    }
    
    as.numeric(value)
  }
  
  active_force <- function(X, t, parameters = list()) {
    check_scalar(t, "t")
    
    amplitude(X, parameters) * pacing_pulse(t, parameters)
  }
  
  list(
    A0 = A0,
    
    # Illustrative isotropic defaults, not calibrated parameters.
    E = material_field("E", 1, 1),
    D = material_field("D", 1, 1),
    
    rho = scalar_field("rho", 1, positive = TRUE),
    kappa = scalar_field("kappa", 0.05, nonnegative = TRUE),
    
    # Available as a field, but not automatically included
    # in the equation generator.
    ceta = scalar_field("ceta", 0.8, nonnegative = TRUE),
    
    P_amplitude = amplitude,
    P = active_force,
    
    pulse = pacing_pulse,
    pulse_expression = pacing_expression,
    pulse_parameters = pacing_parameters
  )
}


# ------------------------------------------------------------
# 7. Numerical fields -> symbolic coefficient names
#
# Spatial fields provide numerical defaults at grid nodes.
# The generated ODE retains named coefficient symbols.
#
# The time-dependent pulse is kept symbolic, not evaluated here.
#
# Note:
# Spatial-function parameter dependencies are sampled numerically.
# A later fitting adapter must explicitly preserve/reconstruct
# shared dependencies such as Tmax or material parameters.
# ------------------------------------------------------------


#builds the connection between “a coefficient name in the equation” and “the numerical value belonging to that name.”


make_field_binding <- function(
    fields,
    parameters,
    grid_bounds,
    spacing,
    origin = c(0, 0)
) {
  required <- c(
    "E", "D", "rho", "kappa",
    "P_amplitude", "pulse_expression", "pulse_parameters"
  )
  
  if (!is.list(fields) ||
      !all(required %in% names(fields)) ||
      !all(vapply(fields[required], is.function, logical(1)))) {
    stop("Missing field functions required for symbolic generation.")
  }
  
  if (!is.list(parameters)) {
    stop("parameters must be a list.")
  }
  
  dimension <- check_field_geometry(
    grid_bounds, spacing, origin
  )
  
  lower <- grid_bounds[["min", exact = TRUE]]
  
  pulse <- fields[["pulse_expression", exact = TRUE]]()
  
  if (!is.character(pulse) ||
      length(pulse) != 1L ||
      is.na(pulse) ||
      !nzchar(pulse)) {
    stop("pulse_expression must return one expression string.")
  }
  
  if (length(parse(text = pulse)) != 1L) {
    stop("pulse_expression must contain exactly one expression.")
  }
  
  pulse_values <- fields[["pulse_parameters", exact = TRUE]](
    parameters
  )
  
  if (!is.numeric(pulse_values) ||
      any(!is.finite(pulse_values)) ||
      is.null(names(pulse_values)) ||
      anyNA(names(pulse_values)) ||
      any(!nzchar(names(pulse_values))) ||
      anyDuplicated(names(pulse_values)) > 0L) {
    stop("pulse_parameters must return a named finite numeric vector.")
  }
  
  coefficient_values <- as.list(pulse_values)
  cache <- new.env(parent = emptyenv())
  
  coefficient <- function(name, index,
                          component = integer(0)) {
    key <- paste(c(name, index), collapse = "_")
    
    if (!exists(key, envir = cache, inherits = FALSE)) {
      X <- origin + (index - lower) * spacing
      
      field_function <- fields[[name, exact = TRUE]]
      
      if (!is.function(field_function)) {
        stop("Missing numerical field function: ", name)
      }
      
      value <- field_function(X, parameters)
      
      if (name %in% c("E", "D")) {
        value <- validate_material(value, dimension, name)
        
      } else if (name == "P_amplitude") {
        check_vector(value, name)
        
        if (length(value) != dimension) {
          stop("P_amplitude has the wrong dimension.")
        }
        
      } else {
        check_scalar(value, name)
        
        if (name == "rho" && value <= 0) {
          stop("rho must be positive.")
        }
        
        if (name %in% c("kappa", "ceta") && value < 0) {
          stop(name, " must be nonnegative.")
        }
      }
      
      assign(key, value, envir = cache)
    }
    
    value <- get(key, envir = cache, inherits = FALSE)
    
    selected <- if (length(component) == 0L) {
      value
    } else if (is.null(dim(value))) {
      value[component]
    } else {
      value[matrix(component, nrow = 1L)]
    }
    
    check_scalar(selected, paste0(name, " component"))
    
    symbol <- paste(c(name, component, index), collapse = "_")
    coefficient_values[[symbol]] <<- as.numeric(selected)
    
    symbol
  }
  
  symbol <- function(
    name,
    index,
    shape,
    component = integer(0),
    state_variables = c("u", "v"),
    fields = NULL,
    parameters = list()
  ) {
    if (name == "P") {
      amplitude_symbol <- coefficient(
        "P_amplitude",
        index,
        component
      )
      
      return(paste0(
        "(", amplitude_symbol, ") * (", pulse, ")"
      ))
    }
    
    if (name %in% c("E", "D", "rho", "kappa", "ceta")) {
      return(coefficient(name, index, component))
    }
    
    stop("Unknown coefficient field: ", name)
  }
  
  list(
    symbol = symbol,
    
    parameter_values = function() {
      unlist(coefficient_values, use.names = TRUE)
    }
  )
}