#' Rearrange Columns in a Data Frame
#'
#' This function rearranges the columns of a data frame based on a specified order.
#'
#' @param data A data frame in which columns need to be rearranged.
#' @param vars A named vector specifying the new positions of the columns. The names of the vector should be the column names, and the values should be the desired positions.
#' @return A data frame with the columns rearranged according to the specified order.
#' @details The function takes a data frame and a named vector as input. The named vector indicates the desired column names and their new positions. It performs sanity checks to ensure that the specified columns exist in the data frame, that the positions are within valid range, and that there are no duplicates. The function returns the data frame with columns rearranged as specified.
#' @export
#'

arrange.vars <- function(data, vars){
  ##stop if not a data.frame (but should work for matrices as well)
  stopifnot(is.data.frame(data))

  ##sort out inputs
  data.nms <- names(data)
  var.nr <- length(data.nms)
  var.nms <- names(vars)
  var.pos <- vars
  ##sanity checks
  stopifnot( !any(duplicated(var.nms)),
             !any(duplicated(var.pos)) )
  stopifnot( is.character(var.nms),
             is.numeric(var.pos) )
  stopifnot( all(var.nms %in% data.nms) )
  stopifnot( all(var.pos > 0),
             all(var.pos <= var.nr) )

  ##prepare output
  out.vec <- character(var.nr)
  out.vec[var.pos] <- var.nms
  out.vec[-var.pos] <- data.nms[ !(data.nms %in% var.nms) ]
  stopifnot( length(out.vec)==var.nr )

  ##re-arrange vars by position
  data <- data[ , out.vec]
  return(data)
}
