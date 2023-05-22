



#' Calculate mode of vector
#'
#' @param v A numeric vector
#'
#' @return A number
#' @export
#'
#' @examples
#' getMode(c(1,2,3,2,3,4,2))
getMode = function(v) {
  uniqv = unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
