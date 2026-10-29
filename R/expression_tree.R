# R/expression_tree.R
#
#        ##
#      ####
#    ######
#        ##
#        ##
#        ##
#        ##
#        ##
#    ##########
#    ##########

#                           
# shape:
#   integer(0)   scalar
#   c(2)         vector
#   c(2, 2)      matrix
#   c(2,2,2,2)   fourth-order tensor

new_expr <- function(op, args = list(),
                     shape = integer(0), name = NULL) {
  
  if (!is.character(op) || length(op) != 1L ||
      is.na(op) || !nzchar(op)) {
    stop("op must be a nonempty string.")
  }
  
  if (!is.list(args) ||
      !all(vapply(args, inherits, logical(1), "math_expr"))) {
    stop("args must be a list of math_expr objects.")
  }
  
  if (!is.numeric(shape) ||
      any(!is.finite(shape)) ||
      any(shape <= 0) ||
      any(shape != floor(shape)) ||
      any(shape > .Machine$integer.max)) {
    stop("shape must contain positive integer dimensions.")
  }
  
  if (!is.null(name) &&
      (!is.character(name) || length(name) != 1L ||
       is.na(name) || !nzchar(name))) {
    stop("name must be NULL or a nonempty string.")
  }
  
  structure(
    list(
      op = op,
      args = args,
      shape = as.integer(shape),
      name = name
    ),
    class = "math_expr"
  )
}

field <- function(name, shape = integer(0)) {
  if (missing(name) || is.null(name)) {
    stop("A field needs a name.")
  }
  
  new_expr("field", shape = shape, name = name)
}

shape_name <- function(shape) {
  if (length(shape) == 0L) {
    return("scalar")
  }
  
  paste(shape, collapse = " x ")
}

show_type <- function(x) {
  stopifnot(inherits(x, "math_expr"))
  
  cat(
    "Operation:", x$op,
    "\nName:", if (is.null(x$name)) "-" else x$name,
    "\nShape:", shape_name(x$shape),
    "\n"
  )
  
  invisible(x)
}


#Die Funktion ist rekursiv. Sie zeigt verschachtelte Ausdrücke an.
print_expr <- function(x, indent = 0L) {
  stopifnot(inherits(x, "math_expr"))
  
  label <- if (x$op == "field") {
    paste0("Field: ", x$name)
  } else {
    paste0("Operator: ", x$op)
  }
  
  cat(
    strrep("  ", indent),
    label, " [", shape_name(x$shape), "]\n",
    sep = ""
  )
  
  for (arg in x$args) {
    print_expr(arg, indent + 1L)
  }
  
  invisible(x)
}