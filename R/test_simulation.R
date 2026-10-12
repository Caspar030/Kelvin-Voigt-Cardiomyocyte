# R/test_simulation.R
#
# 4 horizontal x 3 vertical nodes.
# Left edge fixed.
# Periodic active sarcomere forces.
# Small initial displacement of one internal node.
#
# Output:
#   output/grid_3x4.gif
#   output/grid_3x4_displacement.png
#
# Run from the project root:
# source("R/test_simulation.R")


source("R/equation_generator.R")

required_packages <- c("dMod2", "ggplot2", "gganimate", "gifski")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages)) {
  stop(
    "Please install: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(ggplot2)


# ------------------------------------------------------------
# 1. Geometry and parameters
# ------------------------------------------------------------

bounds <- list(
  min = c(1, 1),
  max = c(4, 3)
)

spacing <- c(0.5, 0.8)
origin <- c(0, 0)

parameters <- list(
  rho = 1,
  
  lambda_E = 1,
  mu_E = 1,
  
  lambda_D = 0.3,
  mu_D = 0.3,
  
  kappa = 0.05,
  
  pace_period = 1,
  pace_duty = 1 / 3,
  pace_start = 0,
  
  Tmax = 0.05
)

# All sarcomeres have the same active strength.
parameters$P_amplitude <- make_sarcomere_force_profile(
  grid_bounds = bounds,
  spacing = spacing,
  origin = origin,
  thickness = 1,
  
  tension = function(X_left, X_right, parameters) {
    parameter_value(parameters, "Tmax")
  }
)

kv_model <- generate_equations(
  grid_bounds = bounds,
  spacing = spacing,
  parameters = parameters,
  origin = origin,
  fixed_side = "left"
)

# 12 nodes - 3 fixed nodes = 9 dynamic nodes.
# Each dynamic node has 2 displacement and 2 velocity components.
stopifnot(
  nrow(kv_model$system) == 36L,
  length(kv_model$fixed_values) == 12L
)

# Small initial inhomogeneity at internal node (2,2).
# This changes initial conditions, not the ongoing forcing.
kv_model$initial_values["u_1_2_2"] <- 0.005
kv_model$initial_values["u_2_2_2"] <- 0.008


# ------------------------------------------------------------
# 2. Replace inline pacing with an external forcing
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


# ------------------------------------------------------------
# 3. Simulation time and forcing
# ------------------------------------------------------------

duration <- 8
times <- seq(0, duration, by = 0.05)

pace_parameters <- as.list(
  kv_model$parameter_values[
    c("pace_period", "pace_duty", "pace_start")
  ]
)

period <- pace_parameters$pace_period
duty <- pace_parameters$pace_duty
start <- pace_parameters$pace_start

cycles <- seq(
  floor((0 - start) / period) - 1,
  ceiling((duration - start) / period) + 1
)

edges <- c(
  start + cycles * period,
  start + (cycles + duty) * period
)

forcing_times <- sort(unique(c(
  seq(0, duration, length.out = 8001),
  times,
  edges[edges >= 0 & edges <= duration]
)))

forcing_data <- data.frame(
  name = "pace_input",
  time = forcing_times,
  value = pacing_pulse(forcing_times, pace_parameters)
)


# ------------------------------------------------------------
# 4. Compile and simulate with dMod2
# ------------------------------------------------------------

build_dir <- file.path(tempdir(), "kv_grid_3x4")
dir.create(build_dir, recursive = TRUE, showWarnings = FALSE)

ode_model <- dMod2::odemodel(
  f = rhs,
  deriv = FALSE,
  forcings = "pace_input",
  modelname = "kv_grid_3x4",
  backend = "cppDE",
  compile = FALSE,
  outdir = build_dir
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
  stop("Missing parameter values: ", paste(missing, collapse = ", "))
}

prediction <- x(
  times,
  available[needed],
  deriv = FALSE
)

prediction_df <- as.data.frame(prediction)
prediction_df <- prediction_df[
  prediction_df$condition == "pacing",
]

stopifnot(all(is.finite(prediction_df$value)))


# ------------------------------------------------------------
# 5. Reconstruct all node positions, including fixed nodes
# ------------------------------------------------------------

nodes <- expand.grid(i = 1:4, j = 1:3)
nodes$id <- seq_len(nrow(nodes))
nodes$X <- origin[1] + (nodes$i - 1) * spacing[1]
nodes$Y <- origin[2] + (nodes$j - 1) * spacing[2]
nodes$fixed <- nodes$i == 1

state_values <- function(variable, component, i, j) {
  name <- paste(variable, component, i, j, sep = "_")
  
  if (name %in% names(kv_model$fixed_values)) {
    return(rep(unname(kv_model$fixed_values[name]), length(times)))
  }
  
  rows <- prediction_df[prediction_df$name == name, ]
  
  if (anyDuplicated(rows$time)) {
    stop("Duplicate time entries for ", name)
  }
  
  values <- rows$value[match(times, rows$time)]
  
  if (anyNA(values)) {
    stop("Missing simulated values for ", name)
  }
  
  values
}

# 1 = actual displacement.
# Set, for example, to 3 only for a clearly labelled magnification.
display_gain <- 1

frames <- do.call(rbind, lapply(seq_len(nrow(nodes)), function(k) {
  node <- nodes[k, ]
  
  ux <- state_values("u", 1, node$i, node$j)
  uy <- state_values("u", 2, node$i, node$j)
  
  data.frame(
    id = node$id,
    i = node$i,
    j = node$j,
    time = times,
    frame = seq_along(times),
    X = node$X,
    Y = node$Y,
    ux = ux,
    uy = uy,
    x = node$X + display_gain * ux,
    y = node$Y + display_gain * uy,
    fixed = node$fixed,
    activation = pacing_pulse(times, pace_parameters)
  )
}))

frame_labels <- sprintf("%.2f s", times)

frames$frame_label <- factor(
  frame_labels[frames$frame],
  levels = frame_labels
)


# ------------------------------------------------------------
# 6. Connect neighboring nodes for visualization
# ------------------------------------------------------------

node_id <- function(i, j) {
  nodes$id[nodes$i == i & nodes$j == j]
}

bonds <- list()

for (j in 1:3) {
  for (i in 1:3) {
    bonds[[length(bonds) + 1L]] <- data.frame(
      from = node_id(i, j),
      to = node_id(i + 1, j),
      type = "Sarkomer"
    )
  }
}

for (j in 1:2) {
  for (i in 1:4) {
    bonds[[length(bonds) + 1L]] <- data.frame(
      from = node_id(i, j),
      to = node_id(i, j + 1),
      type = "Transversal"
    )
  }
}

bonds <- do.call(rbind, bonds)

segments <- do.call(rbind, lapply(seq_len(nrow(bonds)), function(b) {
  left <- frames[frames$id == bonds$from[b], ]
  right <- frames[frames$id == bonds$to[b], ]
  
  stopifnot(identical(left$frame, right$frame))
  
  data.frame(
    frame_label = left$frame_label,
    x = left$x,
    y = left$y,
    xend = right$x,
    yend = right$y,
    type = bonds$type[b],
    activation = left$activation
  )
}))

reference_segments <- data.frame(
  x = nodes$X[bonds$from],
  y = nodes$Y[bonds$from],
  xend = nodes$X[bonds$to],
  yend = nodes$Y[bonds$to]
)

# Fixed axis limits prevent artificial motion from rescaling.
x_limits <- range(c(nodes$X, frames$x)) + c(-0.12, 0.12)
y_limits <- range(c(nodes$Y, frames$y)) + c(-0.12, 0.12)


# ------------------------------------------------------------
# 7. Animated grid
# ------------------------------------------------------------

grid_plot <- ggplot() +
  geom_segment(
    data = reference_segments,
    aes(x = x, y = y, xend = xend, yend = yend),
    colour = "grey80",
    linetype = "dotted",
    linewidth = 0.6
  ) +
  geom_segment(
    data = segments[segments$type == "Transversal", ],
    aes(x = x, y = y, xend = xend, yend = yend),
    colour = "grey55",
    linewidth = 0.7
  ) +
  geom_segment(
    data = segments[segments$type == "Sarkomer", ],
    aes(
      x = x, y = y, xend = xend, yend = yend,
      colour = activation
    ),
    linewidth = 1.5
  ) +
  geom_point(
    data = frames,
    aes(x = x, y = y, shape = fixed),
    size = 4,
    colour = "#182B49"
  ) +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 15),
    breaks = c(FALSE, TRUE),
    labels = c("Beweglich", "Fixiert"),
    name = NULL
  ) +
  scale_colour_gradient(
    low = "#2585B5",
    high = "#E45132",
    limits = c(0, 1),
    name = "Aktivierung"
  ) +
  coord_equal(
    xlim = x_limits,
    ylim = y_limits,
    expand = FALSE
  ) +
  labs(
    title = "Periodische Kontraktion – Zeit: {current_frame}",
    subtitle = paste0(
      "4 × 3 Knoten | links fixiert | Auslenkungsfaktor: ",
      display_gain
    ),
    x = "Longitudinale Position",
    y = "Transversale Position",
    caption = "Gepunktet: Referenzgitter. Rot: aktive Phase."
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  ) +
  gganimate::transition_manual(frame_label)

output_dir <- "output"
dir.create(output_dir, showWarnings = FALSE)

animation <- gganimate::animate(
  grid_plot,
  nframes = length(times),
  fps = 20,
  width = 850,
  height = 650,
  renderer = gganimate::gifski_renderer()
)

gif_file <- file.path(output_dir, "grid_3x4.gif")

gganimate::anim_save(
  filename = gif_file,
  animation = animation
)


# ------------------------------------------------------------
# 8. Time courses: initial transient and periodic response
# ------------------------------------------------------------

trace_data <- rbind(
  transform(
    frames[frames$i == 2 & frames$j == 2, ],
    node = "Innerer Knoten (2,2)"
  ),
  transform(
    frames[frames$i == 4 & frames$j == 2, ],
    node = "Rechter Knoten (4,2)"
  )
)

trace_plot <- ggplot(
  trace_data,
  aes(x = time, y = ux, colour = node)
) +
  geom_hline(yintercept = 0, colour = "grey75") +
  geom_line(linewidth = 0.8) +
  labs(
    title = "Einschwingverhalten und periodische Bewegung",
    x = "Zeit [s]",
    y = "Longitudinale Verschiebung",
    colour = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(
  file.path(output_dir, "grid_3x4_displacement.png"),
  trace_plot,
  width = 9,
  height = 4,
  dpi = 180
)

print(trace_plot)

cat(
  "\nAnimation gespeichert:\n",
  normalizePath(gif_file),
  "\n"
)