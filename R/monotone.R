#' Monotone Regression Function
#'
#' \code{monotone} performs simple linear ordered monotone or isotonic regression.
#' The function follows the up-and-down-blocks implementation (Kruskal, 1964)
#' of the pool-adjacent-violators algorithm (Ayer, Brunk, Ewing, Reid, and Silverman, 1955)
#' with additional lookaheads (Busing, 2021).
#'
#' @param x a real-valued vector.
#' @param w a real-valued vector with positive weights (default a vector with ones).
#'
#' @details Error checking on \code{x} or \code{w} is not present.
#'
#' @return Returns a real-valued vector with values of \code{x} in increasing order.
#'
#' @references Ayer M., H.D. Brunk, G.M. Ewing, W.T. Reid, and E. Silverman (1955). 
#'   An empirical distribution function for sampling with incomplete information. 
#'   \emph{The Annals of Mathematical Statistics}, pp. 641-647. 
#'   URL https://www.jstor.org/stable/pdf/2236377.pdf.
#'
#'   Busing, F.M.T.A. (2022). 
#'   Monotone Regression: A Simple and Fast O(n) PAVA Implementation. 
#'   \emph{Journal of Statistical Software, Code Snippets, 102 (1)}, pp. 1-25. 
#'   (<doi:10.18637/jss.v102.c01>)
#'   
#'
#'   Kruskal, J.B. (1964). 
#'   Nonmetric multidimensional scaling: a numerical method.
#'   \emph{Psychometrika, 29(2)}, pp. 115-129. 
#'   URL http://cda.psych.uiuc.edu/psychometrika_highly_cited_articles/kruskal_1964b.pdf.
#'
#' @examples
#' y <- c( 8, 4, 8, 2, 2, 0, 8 )
#' x <- monotone( y )
#' print( x )
#' 
#' @export
#' 
#' @useDynLib monotone monotoneC
monotone <- function(x, w = rep( 1, length( x ) ) ) {
  n <- length( x )
  return( .C( "monotoneC", as.integer( n ), x = as.double( x ), as.double( w ), PACKAGE = "monotone" )$x )
}