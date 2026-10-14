# R/test_analytic_solution.R



#Testing Analytical vs numerical implementation of the kelvin-voigt model. Is this a validation of the model?
#TODO: do some more tests to compare analytical with numerical solution as part of validation



# Compare the implemented Kelvin–Voigt system with an
# analytical damped longitudinal mode.
#
# Boundary conditions:
#   left fixed, other boundaries free.
#
# No active force and no external forcing required.


source("R/equation_generator.R")
library(dMod2)


# ------------------------------------------------------------
# 1. Exact continuum solution
# ------------------------------------------------------------

A <- 0.01
q <- pi / 2

alpha <- 0.1 * q^2
omega <- sqrt(q^2 + 0.3 - alpha^2)

exact_u <- function(x, t) {
  A * sin(q * x) * exp(-alpha * t) *
    (cos(omega * t) + alpha / omega * sin(omega * t))
}

exact_v <- function(x, t) {
  -A * sin(q * x) * exp(-alpha * t) *
    (q^2 + 0.3) / omega * sin(omega * t)
}


# ------------------------------------------------------------
# 2. Build and simulate one grid
# ------------------------------------------------------------

run_analytic_test <- function(Nx = 17L, Ny = 3L) {
  
  stopifnot(
    Nx >= 3, Ny >= 3,
    Nx == floor(Nx), Ny == floor(Ny)
  )
  
  bounds <- list(
    min = c(1, 1),
    max = c(Nx, Ny)
  )
  
  h <- c(1 / (Nx - 1), 1 / (Ny - 1))
  
  parameters <- list(
    rho = 1,
    kappa = 0.3,
    
    lambda_E = 0,
    mu_E = 0.5,
    
    lambda_D = 0,
    mu_D = 0.1,
    
    P_amplitude = c(0, 0)
  )
  
  model <- generate_equations(
    grid_bounds = bounds,
    spacing = h,
    parameters = parameters,
    fixed_side = "left"
  )
  
  system <- model$system
  
  indices <- t(vapply(
    strsplit(system$index, ",", fixed = TRUE),
    as.numeric,
    numeric(2)
  ))
  
  x_coordinates <- (indices[, 1] - 1) * h[1]
  
  # All velocities and transverse displacements initially zero.
  initial <- model$initial_values
  initial[] <- 0
  
  longitudinal_u <- system$variable == "u" &
    system$component == 1
  
  initial[system$state[longitudinal_u]] <-
    exact_u(x_coordinates[longitudinal_u], 0)
  
  
  # --------------------------------------------------------
  # Remove pacing expression; active force is zero.
  # --------------------------------------------------------
  
  rhs <- setNames(system$rhs, system$state)
  
  rhs[] <- gsub(
    pacing_expression(),
    "0",
    rhs,
    fixed = TRUE
  )
  
  # Insert constant coefficient values directly.
  # This also avoids compiling unused nodal fit parameters.
  constants <- as.list(model$parameter_values)
  
  replace_constants <- function(e) {
    if (is.symbol(e)) {
      name <- as.character(e)
      
      if (name %in% names(constants)) {
        return(unname(constants[[name]]))
      }
      
      return(e)
    }
    
    if (is.call(e)) {
      for (k in seq_along(e)[-1L]) {
        e[[k]] <- replace_constants(e[[k]])
      }
    }
    
    e
  }
  
  rhs <- vapply(rhs, function(s) {
    expression <- replace_constants(parse(text = s)[[1]])
    
    paste(
      deparse(expression, width.cutoff = 500L),
      collapse = " "
    )
  }, character(1))
  
  
  # --------------------------------------------------------
  # Compile and simulate
  # --------------------------------------------------------
  
  outdir <- file.path(
    tempdir(),
    paste0("kv_analytic_", Nx, "_", Ny)
  )
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  ode_model <- dMod2::odemodel(
    f = rhs,
    deriv = FALSE,
    backend = "cppDE",
    modelname = paste0("kv_analytic_", Nx, "_", Ny),
    compile = FALSE,
    outdir = outdir
  )
  
  x <- dMod2::Xs(
    ode_model,
    optionsOde = list(
      atol = 1e-11,
      rtol = 1e-9
    )
  )
  
  dMod2::compile(x)
  
  required <- dMod2::getParameters(x)
  missing <- setdiff(required, names(initial))
  
  if (length(missing)) {
    stop("Missing values: ", paste(missing, collapse = ", "))
  }
  
  times <- seq(0, 5, by = 0.01)
  
  prediction <- x(
    times,
    initial[required],
    deriv = FALSE
  )
  
  df <- as.data.frame(prediction)
  
  stopifnot(all(is.finite(df$value)))
  
  
  # --------------------------------------------------------
  # Compare every dynamic state with its exact solution
  # --------------------------------------------------------
  
  row <- match(df$name, system$state)
  
  if (anyNA(row)) {
    stop("Unexpected state names in simulation output.")
  }
  
  df$exact <- 0
  
  is_u1 <- system$variable[row] == "u" &
    system$component[row] == 1
  
  is_v1 <- system$variable[row] == "v" &
    system$component[row] == 1
  
  df$exact[is_u1] <- exact_u(
    x_coordinates[row[is_u1]],
    df$time[is_u1]
  )
  
  df$exact[is_v1] <- exact_v(
    x_coordinates[row[is_v1]],
    df$time[is_v1]
  )
  
  df$error <- df$value - df$exact
  
  u_rows <- system$variable[row] == "u"
  v_rows <- system$variable[row] == "v"
  transverse <- system$component[row] == 2
  
  summary <- data.frame(
    Nx = Nx,
    Ny = Ny,
    hx = h[1],
    max_u_error = max(abs(df$error[u_rows])),
    max_u_error_percent_of_A =
      100 * max(abs(df$error[u_rows])) / A,
    max_v_error = max(abs(df$error[v_rows])),
    max_transverse_value = max(abs(df$value[transverse]))
  )
  
  print(summary)
  
  list(
    data = df,
    summary = summary,
    Nx = Nx,
    Ny = Ny
  )
}


# ------------------------------------------------------------
# 3. Run the comparison
# ------------------------------------------------------------

result <- run_analytic_test(Nx = 17, Ny = 3)


# ------------------------------------------------------------
# 4. Plot displacement at the right boundary
# ------------------------------------------------------------

selected_state <- paste(
  "u", 1, result$Nx, 2,
  sep = "_"
)

curve <- result$data[
  result$data$name == selected_state,
]

curve <- curve[order(curve$time), ]

draw_comparison <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old))
  
  par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))
  
  plot(
    curve$time, curve$exact,
    type = "l",
    lwd = 2,
    col = "black",
    xlab = "Zeit [s]",
    ylab = "Longitudinale Verschiebung",
    main = "Analytische und numerische Lösung"
  )
  
  lines(
    curve$time, curve$value,
    col = "#D1495B",
    lwd = 2,
    lty = 2
  )
  
  legend(
    "topright",
    legend = c("Analytisch", "Numerisch"),
    col = c("black", "#D1495B"),
    lty = c(1, 2),
    lwd = 2,
    bty = "n"
  )
  
  plot(
    curve$time, curve$error,
    type = "l",
    col = "#2585B5",
    lwd = 2,
    xlab = "Zeit [s]",
    ylab = "Numerisch − analytisch",
    main = "Diskretisierungs- und Integrationsfehler"
  )
  
  abline(h = 0, lty = 3)
}

dir.create("output", showWarnings = FALSE)

pdf("output/analytic_comparison.pdf", width = 9, height = 7)
draw_comparison()
dev.off()

if (interactive()) {
  draw_comparison()
}