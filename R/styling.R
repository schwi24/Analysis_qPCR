#' Change the scientific number style
#'
#' I take no credit for this function. It was taken from Stackoverflow from
#' an answer by Jack Aidley:
#' https://stackoverflow.com/questions/11610377/how-do-i-change-the-formatting-of-numbers-on-an-axis-with-ggplot
#' The function changes the style of numbers from '0e+0' to '0x10^0'.
#'
#' @return A string in "fancy" scientific notation.
#' @param x A string in scientific notation.
#' @export
fancy_scientific <- function(x) {

  ## use scientific notation
  x <- format(x, scientific = TRUE)

  ## keep all the digits  with appropriate quotes
  x <- gsub("0e\\+00","0",x)
  x <- gsub("^(.*)e", "'\\1'e", x)
  ##turn the 'e+' into plotmath format
  x <- gsub("e", "%*%10^", x)

  ## return this as an expression
  parse(text=x)
}
