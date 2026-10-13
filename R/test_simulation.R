# R/test_simulation.R
#
# Requires our latest model files and dMod2.
#
# Geometry:
#   8 horizontal x 5 vertical nodes
#   left edge fixed, all other edges free
#
# Activation:
#   smooth rise -> plateau -> smooth fall -> relaxation
#
# Output:
#   output/contraction_dynamics.pdf
#   output/contraction_grid.pdf
#
# Live visualization:
#   animate_grid()


source("R/equation_generator.R")
library(dMod2)


# ------------------------------------------------------------
# 1. Settings
# ------------------------------------------------------------

bounds <- list(min = c(1, 1), max = c(8, 5))
spacing <- c(0.5, 0.8)

duration <- 4
output_step <- 0.005

period <- 1
active_duration <- period / 3
ramp_duration <- 0.025

# Start after a short initial relaxation phase.
pulse_start <- 0.3

# Target maximum longitudinal shortening at static equilibrium.
target_shortening <- 0.15

parameters <- list(
  rho = 0.0001,
  
  lambda_E = 0.8,
  mu_E = 0.8,
  
  # D/E = 0.02 s for the material terms.
  # In a nearly overdamped proportional mode, this corresponds
  # to roughly 0.08 s to reach 98% of a step response.
  # The full lattice response is checked below.
  lambda_D = 0.016,
  mu_D = 0.016,
  
  kappa = 0.1,
  
  # Unit force first; amplitude is calibrated below.
  Tmax = 1,
  
  pace_period = period,
  pace_duty = active_duration / period,
  pace_start = pulse_start
)

parameters$P_amplitude <- make_sarcomere_force_profile(
  grid_bounds = bounds,
  spacing = spacing,
  thickness = 1,
  
  tension = function(X_left, X_right, parameters) {
    parameter_value(parameters, "Tmax")
  }
)

kv_model <- generate_equations(
  grid_bounds = bounds,
  spacing = spacing,
  parameters = parameters,
  fixed_side = "left"
)


# ------------------------------------------------------------
# 2. Replace the original pulse with an external activation
# ------------------------------------------------------------

rhs <- setNames(
  kv_model$system$rhs,
  kv_model$system$state
)

pulse_text <- pacing_expression()
stopifnot(any(grepl(pulse_text, rhs, fixed = TRUE)))

rhs[] <- gsub(
  pulse_text,
  "pace_input",
  rhs,
  fixed = TRUE
)

activation <- function(t) {
  tau <- (t - pulse_start) %% period
  
  smoothstep <- function(z) {
    z <- pmax(0, pmin(1, z))
    z^2 * (3 - 2*z)
  }
  
  ifelse(
    t < pulse_start,
    0,
    smoothstep(tau / ramp_duration) *
      smoothstep((active_duration - tau) / ramp_duration)
  )
}


# ------------------------------------------------------------
# 3. Calibrate active force using static equilibrium
#
# With v = 0 and activation = 1:
#   0 = K*u + b
#
# All sarcomere forces are scaled together, preserving
# their equal-and-opposite force balance.
# ------------------------------------------------------------

system <- kv_model$system

u_rows <- which(system$variable == "u")
v_rows <- which(system$variable == "v")

u_names <- system$state[u_rows]
v_names <- system$state[v_rows]

stopifnot(
  length(u_names) == length(v_names),
  identical(sub("^u_", "v_", u_names), v_names)
)

acceleration_expressions <- lapply(
  rhs[v_names],
  function(s) parse(text = s)[[1]]
)

env <- list2env(
  c(
    as.list(kv_model$parameter_values),
    as.list(kv_model$initial_values)
  ),
  parent = baseenv()
)

env$pace_input <- 0

evaluate_acceleration <- function() {
  vapply(acceleration_expressions, function(e) {
    eval(e, envir = env)
  }, numeric(1))
}

baseline <- evaluate_acceleration()
K <- matrix(0, length(v_names), length(u_names))

for (j in seq_along(u_names)) {
  env[[u_names[j]]] <- 1
  K[, j] <- evaluate_acceleration() - baseline
  env[[u_names[j]]] <- 0
}

env$pace_input <- 1
b <- evaluate_acceleration() - baseline

u_static_unit <- setNames(
  as.numeric(solve(K, -b)),
  u_names
)

static_component <- function(values, component, i, j) {
  name <- paste("u", component, i, j, sep = "_")
  
  if (name %in% names(kv_model$fixed_values)) {
    return(0)
  }
  
  unname(values[name])
}

unit_shortening <- numeric(0)

for (j in 1:5) {
  for (i in 1:7) {
    du <- static_component(u_static_unit, 1, i + 1, j) -
      static_component(u_static_unit, 1, i, j)
    
    unit_shortening <- c(unit_shortening, -du / spacing[1])
  }
}

if (max(unit_shortening) <= 0) {
  stop("The static solution does not show longitudinal shortening.")
}

force_scale <- target_shortening / max(unit_shortening)

amplitude_names <- grep(
  "^P_amplitude_",
  names(kv_model$parameter_values),
  value = TRUE
)

stopifnot(length(amplitude_names) > 0)

kv_model$parameter_values[amplitude_names] <-
  kv_model$parameter_values[amplitude_names] * force_scale

u_static <- u_static_unit * force_scale

cat(
  "Peak force per sarcomere:", force_scale, "\n",
  "Target static longitudinal shortening:",
  100 * target_shortening, "%\n"
)


# ------------------------------------------------------------
# 4. Initial asymmetry
# ------------------------------------------------------------

set.seed(42)

kv_model$initial_values[u_names] <- runif(
  length(u_names), -0.003, 0.003
)

# Initial velocities remain zero.


# ------------------------------------------------------------
# 5. Forcing table and dMod2 simulation
# ------------------------------------------------------------

times <- seq(0, duration, by = output_step)

cycle_starts <- pulse_start +
  seq(-1, ceiling(duration / period)) * period

edges <- c(
  cycle_starts,
  cycle_starts + ramp_duration,
  cycle_starts + active_duration - ramp_duration,
  cycle_starts + active_duration
)

edges <- edges[edges >= 0 & edges <= duration]

forcing_times <- sort(unique(c(
  seq(0, duration, by = 0.0005),
  times,
  edges
)))

forcing_data <- data.frame(
  name = "pace_input",
  time = forcing_times,
  value = activation(forcing_times)
)

outdir <- file.path(tempdir(), "kv_fast_plateau")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ode_model <- dMod2::odemodel(
  f = rhs,
  forcings = "pace_input",
  deriv = FALSE,
  backend = "cppDE",
  modelname = "kv_fast_plateau",
  compile = FALSE,
  outdir = outdir
)

x <- dMod2::Xs(
  ode_model,
  forcing_data,
  condition = "pacing",
  optionsOde = list(
    atol = 1e-10,
    rtol = 1e-8
  )
)

dMod2::compile(x)

available <- c(
  kv_model$parameter_values,
  kv_model$initial_values
)

needed <- dMod2::getParameters(x)
missing <- setdiff(needed, names(available))

if (length(missing)) {
  stop("Missing parameters: ", paste(missing, collapse = ", "))
}

prediction <- x(
  times,
  available[needed],
  deriv = FALSE
)

prediction_df <- as.data.frame(prediction)

df <- prediction_df[
  prediction_df$condition == "pacing",
]

stopifnot(all(is.finite(df$value)))


# ------------------------------------------------------------
# 6. Reconstruct node positions
# ------------------------------------------------------------

nodes <- expand.grid(i = 1:8, j = 1:5)
nodes$X <- (nodes$i - 1) * spacing[1]
nodes$Y <- (nodes$j - 1) * spacing[2]

displacement <- function(component, i, j) {
  name <- paste("u", component, i, j, sep = "_")
  
  if (name %in% names(kv_model$fixed_values)) {
    return(rep(0, length(times)))
  }
  
  rows <- df[df$name == name, ]
  stopifnot(!anyDuplicated(rows$time))
  
  values <- rows$value[match(times, rows$time)]
  stopifnot(all(is.finite(values)))
  
  values
}

ux <- vapply(seq_len(nrow(nodes)), function(k) {
  displacement(1, nodes$i[k], nodes$j[k])
}, numeric(length(times)))

uy <- vapply(seq_len(nrow(nodes)), function(k) {
  displacement(2, nodes$i[k], nodes$j[k])
}, numeric(length(times)))

xx <- sweep(ux, 2, nodes$X, "+")
yy <- sweep(uy, 2, nodes$Y, "+")


# ------------------------------------------------------------
# 7. Length diagnostics
# No arbitrary 10% cutoff.
# ------------------------------------------------------------

strains <- list()
reversed_order <- FALSE

for (j in 1:5) {
  ids <- which(nodes$j == j)
  
  for (i in 1:7) {
    a <- ids[i]
    b <- ids[i + 1]
    
    dx <- xx[, b] - xx[, a]
    dy <- yy[, b] - yy[, a]
    
    strains[[length(strains) + 1L]] <-
      sqrt(dx^2 + dy^2) / spacing[1] - 1
    
    reversed_order <- reversed_order || any(dx <= 0)
  }
}

cat(
  "Observed sarcomere length changes [%]:",
  100 * range(unlist(strains)), "\n"
)

if (reversed_order) {
  warning("Longitudinal node order reverses: inspect the geometry.")
}


# ------------------------------------------------------------
# 8. Quantify transition times at the right boundary
#
# Time to enter AND remain in a 2% band during the tested phase.
# This measures one observable, not every internal mode.
# ------------------------------------------------------------

right_ids <- which(nodes$i == 8)
shortening <- -rowMeans(ux[, right_ids, drop = FALSE])

static_right_shortening <- -mean(vapply(1:5, function(j) {
  static_component(u_static, 1, 8, j)
}, numeric(1)))

if (static_right_shortening <= 0) {
  stop("Expected a positive static shortening of the right boundary.")
}

# Analyze a later cycle to reduce initial-transient influence.
cycle_start <- pulse_start + 2 * period

settling_time <- function(start, end, target) {
  ids <- which(times >= start & times <= end)
  
  if (!length(ids)) return(NA_real_)
  
  tolerance <- 0.02 * abs(static_right_shortening)
  inside <- abs(shortening[ids] - target) <= tolerance
  
  # TRUE only if all subsequent samples remain inside the band.
  remains_inside <- rev(cumprod(as.integer(rev(inside)))) == 1
  first <- which(remains_inside)[1]
  
  if (is.na(first)) return(NA_real_)
  
  times[ids[first]] - start
}

contraction_time <- settling_time(
  start = cycle_start,
  end = cycle_start + active_duration - ramp_duration,
  target = static_right_shortening
)

relaxation_start <- cycle_start + active_duration - ramp_duration

relaxation_time <- settling_time(
  start = relaxation_start,
  end = cycle_start + period - output_step,
  target = 0
)

cat(
  "Contraction settling time [s]:", contraction_time, "\n",
  "Relaxation settling time [s]:", relaxation_time, "\n",
  "NA means the 2% band was not reached and maintained.\n"
)


# ------------------------------------------------------------
# 9. Draw actual node positions
# ------------------------------------------------------------

palette <- colorRampPalette(c("#2585B5", "#E45132"))(101)

xlim <- range(c(nodes$X, xx)) + c(-0.15, 0.15)
ylim <- range(c(nodes$Y, yy)) + c(-0.15, 0.15)

draw_frame <- function(k) {
  level <- activation(times[k])
  bond_colour <- palette[1 + round(100 * level)]
  
  plot(
    NA,
    xlim = xlim, ylim = ylim, asp = 1,
    xlab = "Longitudinale Position",
    ylab = "Transversale Position",
    main = sprintf("Gitterkontraktion | t = %.2f s", times[k])
  )
  
  for (j in 1:5) {
    ids <- which(nodes$j == j)
    
    lines(
      nodes$X[ids], nodes$Y[ids],
      col = "grey80", lty = 3
    )
    
    lines(
      xx[k, ids], yy[k, ids],
      col = bond_colour, lwd = 2
    )
  }
  
  for (i in 1:8) {
    ids <- which(nodes$i == i)
    
    lines(
      nodes$X[ids], nodes$Y[ids],
      col = "grey80", lty = 3
    )
    
    lines(
      xx[k, ids], yy[k, ids],
      col = "grey55"
    )
  }
  
  points(
    xx[k, ], yy[k, ],
    pch = ifelse(nodes$i == 1, 15, 16),
    col = "#182B49",
    cex = 0.9
  )
  
  legend(
    "topright",
    legend = sprintf("Aktivierung: %.2f", level),
    pch = 21,
    pt.bg = if (level > 0.01) "#F2C94C" else "grey85",
    col = "grey35",
    pt.cex = 1.6,
    bty = "n"
  )
}

animate_grid <- function(pause = 0.03) {
  old <- par(no.readonly = TRUE)
  on.exit(par(old))
  
  par(mfrow = c(1, 1))
  
  for (k in seq(1, length(times), by = 6)) {
    draw_frame(k)
    Sys.sleep(pause)
  }
  
  invisible(NULL)
}


# ------------------------------------------------------------
# 10. Save plots
# ------------------------------------------------------------

dir.create("output", showWarnings = FALSE)

pdf("output/contraction_dynamics.pdf", width = 10, height = 6)
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

plot(
  times, activation(times),
  type = "l", lwd = 2, col = "#E45132",
  xlab = "Zeit [s]", ylab = "Aktivierung",
  main = "Aktive Kraft: Anstieg, Haltephase, Abfall"
)

plot(
  times, shortening,
  type = "l", lwd = 2, col = "#2585B5",
  xlab = "Zeit [s]", ylab = "Mittlere Verkürzung am rechten Rand",
  main = "Mechanische Antwort"
)

abline(h = 0, col = "grey60", lty = 2)
abline(
  h = static_right_shortening,
  col = "#E45132", lty = 2
)

legend(
  "topright",
  legend = c("Simulation", "Statisches kontrahiertes Gleichgewicht"),
  col = c("#2585B5", "#E45132"),
  lty = c(1, 2),
  bty = "n"
)

dev.off()

pdf("output/contraction_grid.pdf", width = 12, height = 8)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

snapshot_times <- cycle_start + c(0, 0.05, 0.15, 0.30, 0.40, 0.65)

for (t in snapshot_times) {
  draw_frame(which.min(abs(times - t)))
}

dev.off()

cat(
  "\nPlots saved in output/.\n",
  "Live animation: animate_grid()\n"
)











#creting high-res gif. #########Needs: sudo apt install ffmpeg

create_grid_gif <- function(
    from = max(min(times), max(times) - 3),
    to = max(times),
    fps = 30,
    slowdown = 4,
    width = 1000,
    height = 800,
    filename = "output/grid_contraction_highres.gif"
) {
  
  # ----------------------------------------------------------
  # 1. Check requirements
  # ----------------------------------------------------------
  
  ffmpeg <- Sys.which("ffmpeg")
  
  if (!nzchar(ffmpeg)) {
    stop(
      "FFmpeg wurde nicht gefunden.\n",
      "Unter Ubuntu/Debian im Terminal installieren:\n",
      "sudo apt install ffmpeg"
    )
  }
  
  if (!exists("draw_frame", mode = "function")) {
    stop("Bitte zuerst die Simulation mit draw_frame() ausführen.")
  }
  
  if (!is.numeric(times) ||
      length(times) < 2L ||
      any(!is.finite(times)) ||
      any(diff(times) <= 0)) {
    stop("times muss endlich und streng aufsteigend sein.")
  }
  
  settings <- list(
    from = from,
    to = to,
    fps = fps,
    slowdown = slowdown,
    width = width,
    height = height
  )
  
  valid <- vapply(settings, function(value) {
    is.numeric(value) &&
      length(value) == 1L &&
      is.finite(value)
  }, logical(1))
  
  if (!all(valid) ||
      from < min(times) ||
      to > max(times) ||
      to <= from ||
      fps <= 0 ||
      slowdown <= 0 ||
      width < 100 ||
      height < 100) {
    stop("Ungültige Export-Einstellungen.")
  }
  
  width <- as.integer(round(width))
  height <- as.integer(round(height))
  
  
  # ----------------------------------------------------------
  # 2. Prepare output and temporary files
  # ----------------------------------------------------------
  
  dir.create(
    dirname(filename),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  gif_file <- file.path(
    normalizePath(dirname(filename)),
    basename(filename)
  )
  
  frame_dir <- tempfile("grid_gif_")
  dir.create(frame_dir)
  
  on.exit(
    unlink(frame_dir, recursive = TRUE),
    add = TRUE
  )
  
  palette_file <- file.path(frame_dir, "palette.png")
  log_file <- file.path(frame_dir, "ffmpeg.log")
  
  frame_pattern <- file.path(
    frame_dir,
    "frame_%05d.png"
  )
  
  
  # ----------------------------------------------------------
  # 3. Select simulated frames
  #
  # fps refers to simulated time.
  # slowdown changes playback speed only.
  # ----------------------------------------------------------
  
  nframes <- max(2L, round((to - from) * fps))
  
  target_times <- from +
    (0:(nframes - 1L)) * (to - from) / nframes
  
  indices <- vapply(target_times, function(t) {
    which.min(abs(times - t))
  }, integer(1))
  
  if (anyDuplicated(indices)) {
    warning(
      "Einige Frames verwenden denselben Simulationszeitpunkt. ",
      "Für zusätzliche echte Zwischenbilder feinere Ausgabezeiten ",
      "in der Simulation verwenden."
    )
  }
  
  playback_fps <- nframes / ((to - from) * slowdown)
  
  cat(
    "Erzeuge", nframes, "Bilder mit",
    width, "x", height, "Pixeln.\n",
    "Abspieldauer ungefähr:",
    round((to - from) * slowdown, 2), "Sekunden.\n"
  )
  
  
  # ----------------------------------------------------------
  # 4. Save individual PNGs
  # ----------------------------------------------------------
  
  save_frame <- function(k, path) {
    png(
      filename = path,
      width = width,
      height = height,
      res = 120,
      bg = "white"
    )
    
    on.exit(dev.off(), add = TRUE)
    
    par(
      mfrow = c(1, 1),
      mar = c(4, 4, 3, 1),
      cex = 0.9
    )
    
    draw_frame(k)
  }
  
  for (m in seq_along(indices)) {
    path <- file.path(
      frame_dir,
      sprintf("frame_%05d.png", m)
    )
    
    save_frame(indices[m], path)
    
    if (m %% 20L == 0L || m == nframes) {
      cat("Bilder:", m, "/", nframes, "\n")
    }
  }
  
  
  # ----------------------------------------------------------
  # 5. FFmpeg helper
  # ----------------------------------------------------------
  
  run_ffmpeg <- function(arguments) {
    status <- system2(
      command = ffmpeg,
      args = arguments,
      stdout = log_file,
      stderr = log_file
    )
    
    if (status != 0L) {
      if (file.exists(log_file)) {
        cat(
          paste(
            tail(readLines(log_file, warn = FALSE), 30),
            collapse = "\n"
          ),
          "\n"
        )
      }
      
      stop("FFmpeg-Export fehlgeschlagen; siehe Ausgabe oben.")
    }
  }
  
  rate <- sprintf("%.8f", playback_fps)
  
  
  # ----------------------------------------------------------
  # 6. First pass: generate a shared colour palette
  #
  # Separate passes avoid buffering the entire animation
  # while waiting for the palette.
  # ----------------------------------------------------------
  
  cat("Erzeuge Farbpalette ...\n")
  
  run_ffmpeg(c(
    "-hide_banner",
    "-loglevel", "error",
    "-y",
    "-framerate", rate,
    "-start_number", "1",
    "-i", shQuote(frame_pattern),
    "-vf", shQuote("palettegen=max_colors=256:stats_mode=full"),
    "-frames:v", "1",
    "-update", "1",
    shQuote(palette_file)
  ))
  
  
  # ----------------------------------------------------------
  # 7. Second pass: encode GIF using the palette
  # ----------------------------------------------------------
  
  cat("Erzeuge GIF ...\n")
  
  run_ffmpeg(c(
    "-hide_banner",
    "-loglevel", "error",
    "-y",
    "-framerate", rate,
    "-start_number", "1",
    "-i", shQuote(frame_pattern),
    "-i", shQuote(palette_file),
    "-lavfi",
    shQuote("paletteuse=dither=sierra2_4a"),
    "-frames:v", as.character(nframes),
    "-loop", "0",
    shQuote(gif_file)
  ))
  
  if (!file.exists(gif_file) ||
      file.info(gif_file)$size <= 0) {
    stop("Es wurde keine gültige GIF-Datei erzeugt.")
  }
  
  size_mb <- file.info(gif_file)$size / 1024^2
  
  cat(
    "\nGIF gespeichert:\n",
    gif_file,
    "\nDateigröße:",
    round(size_mb, 2), "MB\n"
  )
  
  invisible(gif_file)
}


# ------------------------------------------------------------
# Export: last 3 simulated seconds, four-times slower playback
# ------------------------------------------------------------

gif_file <- create_grid_gif(
  fps = 30,
  slowdown = 4,
  width = 1000,
  height = 800
)