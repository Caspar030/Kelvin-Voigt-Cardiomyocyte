# R/operators.R
# Requires expression_tree.R


#    ########
#   ##      ##
#           ##
#          ##
#        ##
#      ##
#    ##
#   ##
#   ##########
#   ##########


# Diese Datei baut mithilfe der Expressions_tree den Ausdrucksbaum auf. Hier werden die Operationen symbolisch definiert.
# So etwa Addition, Skalarmultiplikation, das Tensorprodukt oder eben die Divergenz.
# Die Funktion arbeitet mit den "math_expression" objects aus expression_tree.R

check_expr <- function(x) {
  if (!inherits(x, "math_expr")) {
    stop("Expected a math_expr object.")
  }
  
  invisible(TRUE)
}

check_dimension <- function(dimension) {
  if (!is.numeric(dimension) ||
      length(dimension) != 1L ||
      !is.finite(dimension) ||
      !(dimension %in% 1:3)) {
    stop("dimension must be 1, 2 or 3.")
  }
  
  as.integer(dimension)
}


# ------------------------------------------------------------
# Gradient
# grad(u)[i,j] = derivative of u[i] with respect to X[j]
# ------------------------------------------------------------

grad <- function(x, dimension = 2L) {
  check_expr(x)
  dimension <- check_dimension(dimension)
  
  new_expr(
    "grad",
    args = list(x),
    shape = c(x$shape, dimension)
  )
}


# ------------------------------------------------------------
# Divergence: contracts the last index with spatial derivative
# div(S)[i] = sum_j derivative of S[i,j] with respect to X[j]
# ------------------------------------------------------------

div <- function(x, dimension = 2L) {
  check_expr(x)
  dimension <- check_dimension(dimension)
  
  if (length(x$shape) == 0L ||
      tail(x$shape, 1L) != dimension) {
    stop("Last tensor dimension must match spatial dimension.")
  }
  
  new_expr(
    "div",
    args = list(x),
    shape = head(x$shape, -1L)
  )
}


# ------------------------------------------------------------
# Matrix-vector or matrix-matrix multiplication
# Not the fourth-order material contraction!
# ------------------------------------------------------------

matmul <- function(A, B) {
  check_expr(A)
  check_expr(B)
  
  if (length(A$shape) != 2L ||
      !(length(B$shape) %in% c(1L, 2L))) {
    stop("matmul requires a matrix and a vector or matrix.")
  }
  
  if (A$shape[2] != B$shape[1]) {
    stop("Incompatible dimensions in matmul.")
  }
  
  new_expr(
    "matmul",
    args = list(A, B),
    shape = c(A$shape[1], B$shape[-1L])
  )
}


# ------------------------------------------------------------
# Double contraction: fourth-order tensor with a matrix
# result[i,j] = sum_{k,l} A[i,j,k,l] * B[k,l]
# ------------------------------------------------------------

contract <- function(A, B) {
  check_expr(A)
  check_expr(B)
  
  if (length(A$shape) != 4L || length(B$shape) != 2L) {
    stop("contract requires a fourth-order tensor and a matrix.")
  }
  
  if (!identical(A$shape[3:4], B$shape)) {
    stop("Contracted dimensions do not match.")
  }
  
  new_expr(
    "contract",
    args = list(A, B),
    shape = A$shape[1:2]
  )
}


# ------------------------------------------------------------
# Addition
# ------------------------------------------------------------

add <- function(x, y) {
  check_expr(x)
  check_expr(y)
  
  if (!identical(x$shape, y$shape)) {
    stop("Addition requires identical tensor shapes.")
  }
  
  new_expr("add", args = list(x, y), shape = x$shape)
}


# ------------------------------------------------------------
# Multiplication by a scalar expression
# ------------------------------------------------------------

multiply <- function(x, y) {
  check_expr(x)
  check_expr(y)
  
  x_scalar <- length(x$shape) == 0L
  y_scalar <- length(y$shape) == 0L
  
  if (!x_scalar && !y_scalar) {
    stop("multiply requires at least one scalar.")
  }
  
  new_expr(
    "multiply",
    args = list(x, y),
    shape = if (x_scalar) y$shape else x$shape
  )
}