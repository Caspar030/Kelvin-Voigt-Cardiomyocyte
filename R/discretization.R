# R/discretization.R
#

#   ##########
#   ##########
#   ##
#   ##
#   ########
#          ##
#           ##
#           ##
#   ##      ##
#    ########


#translates symbolic expressions like div(E grad u) into specific operations for a node.

# Supported operators:
#   field, grad, div, contract, matmul, add, multiply
#
# Convention:
#   grad(u)[i,j] = d u[i] / d X[j]
#   contract(C,G)[i,j] = sum_{k,l} C[i,j,k,l] * G[k,l]
#
# Discretization:
#   - uniform Cartesian grid
#   - divergence from shared face fluxes
#   - normal derivatives from neighboring node differences
#   - tangential face derivatives from averaged node derivatives
#   - arithmetic interpolation of coefficients to faces
#   - zero exterior flux when grid_bounds is supplied
#   - half-width control volumes at boundary nodes
#
# Output:
#   scalar: character string
#   vector: character vector
#   higher-order tensor: character array
#
# Material coefficients remain symbolic unless field_symbol
# explicitly replaces them by numerical values.


# ------------------------------------------------------------
# 1. Naming and grid utilities
# ------------------------------------------------------------

state_name <- function(variable, index) {
  paste(c(variable, as.integer(index)), collapse = "_")
}

neighbor_index <- function(index, dimension, direction) {
  result <- index
  result[dimension] <- result[dimension] + direction
  result
}

default_field_symbol <- function(
    name,
    index,
    shape,
    component = integer(0),
    state_variables = c("u", "v"),
    fields = NULL,
    parameters = list()
) {
  # Examples:
  #   scalar u:       u_3_4
  #   vector u[1]:    u_1_3_4
  #   E[1,2,1,2]:     E_1_2_1_2_3_4
  paste(
    c(name, as.integer(component), as.integer(index)),
    collapse = "_"
  )
}




# ------------------------------------------------------------
# 2. First-derivative stencils at nodes
# ------------------------------------------------------------

derivative_stencil <- function(
    index,
    dimension,
    grid_bounds = NULL
) {
  if (!is.null(grid_bounds)) {
    lower <- grid_bounds$min[dimension]
    upper <- grid_bounds$max[dimension]
    position <- index[dimension]
    
    if (position < lower || position > upper) {
      stop("Index lies outside grid_bounds.")
    }
    
    if (upper <= lower) {
      stop("Each axis needs at least two nodes.")
    }
    
    if (position == lower) {
      return(list(
        points = list(
          index,
          neighbor_index(index, dimension, 1L)
        ),
        coefficients = c(-1, 1),
        denominator = 1
      ))
    }
    
    if (position == upper) {
      return(list(
        points = list(
          neighbor_index(index, dimension, -1L),
          index
        ),
        coefficients = c(-1, 1),
        denominator = 1
      ))
    }
  }
  
  list(
    points = list(
      neighbor_index(index, dimension, -1L),
      neighbor_index(index, dimension, 1L)
    ),
    coefficients = c(-1, 1),
    denominator = 2
  )
}


# ------------------------------------------------------------
# 3. Recursive symbolic discretization
# ------------------------------------------------------------

discretize <- function(
    expr,
    index,
    spacing,
    fields = NULL,
    parameters = list(),
    field_symbol = default_field_symbol,
    state_variables = c("u", "v"),
    grid_bounds = NULL
) {
  if (!inherits(expr, "math_expr")) {
    stop("expr must be a math_expr.")
  }
  
  if (!is.numeric(index) ||
      length(index) < 1L ||
      any(!is.finite(index)) ||
      any(index != floor(index))) {
    stop("index must contain finite integers.")
  }
  
  dimension <- length(index)
  
  if (!is.numeric(spacing) ||
      length(spacing) != dimension ||
      any(!is.finite(spacing)) ||
      any(spacing <= 0)) {
    stop("spacing must contain one positive finite value per axis.")
  }
  
  if (!is.function(field_symbol)) {
    stop("field_symbol must be a function.")
  }
  
  if (!is.null(grid_bounds)) {
    lower <- grid_bounds$min
    upper <- grid_bounds$max
    
    if (!is.numeric(lower) ||
        !is.numeric(upper) ||
        length(lower) != dimension ||
        length(upper) != dimension ||
        any(!is.finite(lower)) ||
        any(!is.finite(upper)) ||
        any(lower != floor(lower)) ||
        any(upper != floor(upper)) ||
        any(lower >= upper)) {
      stop("Invalid grid bounds: need at least two nodes per axis.")
    }
    
    if (any(index < lower) || any(index > upper)) {
      stop("index lies outside grid_bounds.")
    }
  }
  
  # ----------------------------------------------------------
  # Symbolic arithmetic
  # ----------------------------------------------------------
  
  number <- function(x) {
    sprintf("%.17g", x)
  }
  
  sum_strings <- function(x) {
    x <- x[x != "0"]
    
    if (length(x) == 0L) return("0")
    if (length(x) == 1L) return(x)
    
    paste0("(", paste(x, collapse = " + "), ")")
  }
  
  product <- function(a, b) {
    if (a == "0" || b == "0") return("0")
    if (a == "1") return(b)
    if (b == "1") return(a)
    
    paste0("(", a, ") * (", b, ")")
  }
  
  difference <- function(a, b) {
    if (identical(a, b)) return("0")
    
    paste0("((", a, ") - (", b, "))")
  }
  
  quotient <- function(a, denominator) {
    if (a == "0") return("0")
    
    paste0("(", a, ") / (", number(denominator), ")")
  }
  
  average <- function(a, b) {
    if (identical(a, b)) return(a)
    
    quotient(sum_strings(c(a, b)), 2)
  }
  
  # ----------------------------------------------------------
  # Evaluate one field component at a node
  # ----------------------------------------------------------
  
  field_at_node <- function(e, idx, component) {
    value <- field_symbol(
      name = e$name,
      index = idx,
      shape = e$shape,
      component = component,
      state_variables = state_variables,
      fields = fields,
      parameters = parameters
    )
    
    if (!is.character(value) ||
        length(value) != 1L ||
        is.na(value) ||
        !nzchar(value)) {
      stop("field_symbol must return one nonempty expression string.")
    }
    
    value
  }
  
  # ----------------------------------------------------------
  # Differentiate one component at a node
  # ----------------------------------------------------------
  
  node_derivative <- function(e, idx, component, axis) {
    stencil <- derivative_stencil(
      index = idx,
      dimension = axis,
      grid_bounds = grid_bounds
    )
    
    terms <- vapply(seq_along(stencil$points), function(m) {
      product(
        number(stencil$coefficients[m]),
        value_at(e, stencil$points[[m]], component)
      )
    }, character(1))
    
    quotient(
      sum_strings(terms),
      stencil$denominator * spacing[axis]
    )
  }
  
  # ----------------------------------------------------------
  # Evaluate one expression component
  #
  # face = NULL: evaluate at node idx
  # face = axis: evaluate halfway between idx and idx + e_axis
  # ----------------------------------------------------------
  
  value_at <- function(e, idx, component = integer(0),
                       face = NULL) {
    
    # Field
    if (e$op == "field") {
      left <- field_at_node(e, idx, component)
      
      if (is.null(face)) {
        return(left)
      }
      
      right <- field_at_node(
        e,
        neighbor_index(idx, face, 1L),
        component
      )
      
      return(average(left, right))
    }
    
    # Addition
    if (e$op == "add") {
      return(sum_strings(c(
        value_at(e$args[[1]], idx, component, face),
        value_at(e$args[[2]], idx, component, face)
      )))
    }
    
    # Scalar multiplication
    if (e$op == "multiply") {
      a <- e$args[[1]]
      b <- e$args[[2]]
      
      a_scalar <- length(a$shape) == 0L
      b_scalar <- length(b$shape) == 0L
      
      if (!a_scalar && !b_scalar) {
        stop("multiply requires at least one scalar.")
      }
      
      return(product(
        value_at(
          a, idx,
          if (a_scalar) integer(0) else component,
          face
        ),
        value_at(
          b, idx,
          if (b_scalar) integer(0) else component,
          face
        )
      ))
    }
    
    # Gradient
    if (e$op == "grad") {
      if (tail(e$shape, 1L) != dimension) {
        stop("Gradient dimension does not match grid dimension.")
      }
      
      child <- e$args[[1]]
      axis <- tail(component, 1L)
      child_component <- head(component, -1L)
      
      if (is.null(face)) {
        return(node_derivative(
          child, idx, child_component, axis
        ))
      }
      
      right_idx <- neighbor_index(idx, face, 1L)
      
      # Normal derivative at a face.
      if (axis == face) {
        return(quotient(
          difference(
            value_at(child, right_idx, child_component),
            value_at(child, idx, child_component)
          ),
          spacing[axis]
        ))
      }
      
      # Tangential derivative at a face.
      # Uses immediate and diagonal neighboring nodes.
      return(average(
        node_derivative(child, idx, child_component, axis),
        node_derivative(child, right_idx, child_component, axis)
      ))
    }
    
    # Fourth-order tensor contracted with a matrix.
    if (e$op == "contract") {
      A <- e$args[[1]]
      B <- e$args[[2]]
      
      if (length(A$shape) != 4L ||
          length(B$shape) != 2L ||
          any(A$shape[3:4] != B$shape)) {
        stop("Invalid shapes for contract.")
      }
      
      terms <- character(0)
      
      for (k in seq_len(B$shape[1])) {
        for (l in seq_len(B$shape[2])) {
          terms <- c(terms, product(
            value_at(A, idx, c(component, k, l), face),
            value_at(B, idx, c(k, l), face)
          ))
        }
      }
      
      return(sum_strings(terms))
    }
    
    # Matrix-vector or matrix-matrix multiplication.
    if (e$op == "matmul") {
      A <- e$args[[1]]
      B <- e$args[[2]]
      
      if (length(A$shape) != 2L ||
          !(length(B$shape) %in% c(1L, 2L)) ||
          A$shape[2] != B$shape[1]) {
        stop("Invalid shapes for matmul.")
      }
      
      terms <- vapply(seq_len(A$shape[2]), function(k) {
        b_component <- if (length(B$shape) == 1L) {
          k
        } else {
          c(k, component[2])
        }
        
        product(
          value_at(A, idx, c(component[1], k), face),
          value_at(B, idx, b_component, face)
        )
      }, character(1))
      
      return(sum_strings(terms))
    }
    
    # Conservative divergence.
    if (e$op == "div") {
      if (!is.null(face)) {
        stop("A divergence inside a face expression is unsupported.")
      }
      
      child <- e$args[[1]]
      
      if (length(child$shape) == 0L ||
          tail(child$shape, 1L) != dimension) {
        stop("Last flux index must match spatial dimension.")
      }
      
      terms <- vapply(seq_len(dimension), function(axis) {
        at_lower <- !is.null(grid_bounds) &&
          idx[axis] == grid_bounds$min[axis]
        
        at_upper <- !is.null(grid_bounds) &&
          idx[axis] == grid_bounds$max[axis]
        
        flux_component <- c(component, axis)
        
        flux_plus <- if (at_upper) {
          "0"
        } else {
          value_at(child, idx, flux_component, face = axis)
        }
        
        flux_minus <- if (at_lower) {
          "0"
        } else {
          value_at(
            child,
            neighbor_index(idx, axis, -1L),
            flux_component,
            face = axis
          )
        }
        
        width <- spacing[axis]
        
        if (at_lower || at_upper) {
          width <- width / 2
        }
        
        quotient(
          difference(flux_plus, flux_minus),
          width
        )
      }, character(1))
      
      return(sum_strings(terms))
    }
    
    stop("Unknown expression operator: ", e$op)
  }
  
  # ----------------------------------------------------------
  # Assemble all output components in R array order
  # ----------------------------------------------------------
  
  shape <- expr$shape
  
  if (length(shape) == 0L) {
    return(value_at(expr, index))
  }
  
  components <- as.matrix(expand.grid(
    lapply(shape, seq_len),
    KEEP.OUT.ATTRS = FALSE
  ))
  
  result <- vapply(seq_len(nrow(components)), function(m) {
    value_at(
      expr,
      index,
      as.integer(components[m, ])
    )
  }, character(1))
  
  if (length(shape) == 1L) {
    return(unname(result))
  }
  
  array(result, dim = shape)
}

