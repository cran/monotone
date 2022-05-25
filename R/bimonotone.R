#' Bivariate Monotone Regression Function
#'
#' \code{bimonotone} performs bivariate monotone regression.
#' The function uses the up-and-down-blocks implementation (Kruskal, 1964)
#' of the pool-adjacent-violators algorithm (Ayer, Brunk, Ewing, Reid, and Silverman, 1955),
#' with additional lookaheads, repeatedly, for both rows and columns, until convergence.
#'
#' @param x a real-valued matrix.
#' @param w a real-valued matrix with positive weights (default a matrix with ones).
#' @param maxiter maximum number of iterations (default = 65536)
#' @param eps precision of estimates (default = 1.4901161193847656e-08)
#'
#' @details Error checking on \code{x}, \code{w}, \code{maxiter}, or \code{eps} is not present.
#'
#' @return Returns a real-valued matrix with values of both rows and columns of \code{x} in monotone order.
#'
#' @references Bril G, Dykstra R, Pillers C, Robertson T (1984).
#'   Algorithm AS 206: isotonic regression in two independent variables.
#'   Journal of the Royal Statistical Society. Series C (Applied Statistics), 33(3), 352-357. 
#'   URL https://www.jstor.org/stable/pdf/2347723.pdf.
#' 
#'   Busing, F.M.T.A. (2022). 
#'   Monotone Regression: A Simple and Fast O(n) PAVA Implementation. 
#'   \emph{Journal of Statistical Software, Code Snippets, 102 (1)}, pp. 1-25. 
#'   (<doi:10.18637/jss.v102.c01>)
#'
#'   Dykstra R.L., Robertson T. (1982).
#'   An algorithm for isotonic regression for two or more independent variables.
#'   The Annals of Statistics, 10(3), 708-716.
#'   URL https:  //projecteuclid.org/download/pdf_1/euclid.aos/1176345866.
#'   
#'   Turner, T.R. (2019).
#'   Iso: Functions to Perform Isotonic Regression.
#'   R package version 0.0-18.
#'   URL https://cran.r-project.org/package=Iso
#'   
#' @examples
#' G <- matrix( c( 1, 5.2, 0.1, 0.1, 5, 0, 6, 2, 3, 5.2, 5, 7, 4, 5.5, 6, 6 ), 4, 4 )
#' print( G )
#' H <- bimonotone( G )
#' print( H )
#' y <- c( 8, 4, 8, 2, 2, 0, 8 )
#' x <- bimonotone( as.matrix( y ) )
#' print( x )
#' x <- bimonotone( t( as.matrix( y ) ) )
#' print( x )
#'
#' @export
#'
#' @useDynLib monotone bimonotone
bimonotone <- function( x, w = matrix( 1, nrow( x ), ncol( x ) ), maxiter = 65536, eps = 1.4901161193847656e-08 ) {
  n <- nrow( x )
  m <- ncol( x )
  return( matrix( .C( "bimonotoneC", as.integer( n ), as.integer( m ), x = as.double( x ), as.double( w ), as.integer( maxiter ), as.double( eps ), PACKAGE = "monotone" )$x, n, m ) )
}