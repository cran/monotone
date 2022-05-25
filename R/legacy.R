#' Monotone Regression Legacy Function
#'
#' \code{legacy} provides some functions for monotone regression from the past.
#' Current implementations have been translated into C for proper comparison in Busing (2022).
#'
#' @param x a real-valued vector.
#' @param w a real-valued vector with positive weights (default a vector with ones).
#' @param number function number (specifications below).
#'
#' @details Legacy implementations by number, function, author, and year:
#' \itemize{
#'  \item{0 = default (do nothing)}  
#'  \item{1 = fitm() by Kruskal}{ (1964).}
#'  \item{2 = wmrmnh() by van Waning}{ (1976).}
#'  \item{3 = amalgm() by Cran}{ (1980).}
#'  \item{4 = pav() by Bril}{ (1984).}
#'  \item{5 = isoreg() by Gupta}{ (1995).}
#'  \item{6 = iso_pava() by Turner}{ (1997).}
#'  \item{7 = isotonic() by Kincaid}{ (2001).}
#'  \item{8 = isomean() by Strimmer}{ (2008).}
#'  \item{9 = pooled_pava() by Pedregosa}{ (2011).}
#'  \item{10 = linear_pava() by Tulloch}{ (2014).}
#'  \item{11 = inplace_pava() by Varoquaux}{ (2016).}
#'  \item{12 = md_pava() by Danish}{ (2016).}
#'  \item{13 = reg_1d_l2() by Xu (2017).}
#'  \item{14 = jbkpava() by de Leeuw}{ (2017).}
#' }
#'
#' @details Error checking on \code{w} or \code{x} is not present.
#'
#' @return Returns a real-valued vector with values of \code{x} in increasing order.
#'
#' @references 
#' 
#'   Busing, F.M.T.A. (2022). 
#'   Monotone Regression: A Simple and Fast O(n) PAVA Implementation. 
#'   \emph{Journal of Statistical Software, Code Snippets, 102 (1)}, pp. 1-25. 
#'   (<doi:10.18637/jss.v102.c01>)
#'
#' @examples
#' y <- c( 8, 4, 8, 2, 2, 0, 8 )
#' x <- legacy( y, number = 1 )
#' print( x )
#' 
#' @export
#'
#' @useDynLib monotone legacyC
legacy <- function( x, w = rep( 1, length( x ) ), number = 0 )
{
  if ( number < 0 || number > 14 ) number <- 0
  n <- length( x )
  if ( number == 0 ) return( .C( "fake", n = as.integer(n), x = as.double(x), w = as.double(w), PACKAGE = "monotone" )$x )
  else if ( number ==  1 ) return( .C( "fitm", n = as.integer( n ), x = as.double( x ), w = as.double( w ), PACKAGE = "monotone" )$x )
  else if ( number ==  2 ) return( .C( "wmrmnh", n = as.integer( n ), x = as.double( x ), w = as.double( w ), PACKAGE = "monotone" )$x )
  else if ( number ==  3 ) return( .C( "amalgm", n = as.integer( n ), x = as.double( x ), w = as.double( w ), PACKAGE = "monotone" )$x )
  else if ( number ==  4 ) return( .C( "pav", n = as.integer( n ), x = as.double( x ), w = as.double( w ), PACKAGE = "monotone" )$x )
  else if ( number ==  5 ) return( .C( "isoreg", n = as.integer( n ), x = as.double( x ), w = as.double( w ), PACKAGE = "monotone" )$x )
  else if ( number ==  6 ) return( .C( "iso_pava", n = as.integer( n ), x = as.double( x ), w = as.double( w ), PACKAGE = "monotone" )$x )
  else if ( number ==  7 ) return( .C( "isotonic", n = as.integer( n ), x = as.double( x ), w = as.double( w ), PACKAGE = "monotone" )$x )
  else if ( number ==  8 ) return( .C( "isomean", n = as.integer( n ), x = as.double( x ), w = as.double( w ), PACKAGE = "monotone" )$x )
  else if ( number ==  9 ) return( .C( "pooled_pava", n = as.integer( n ), x = as.double( x ), w = as.double( w ), PACKAGE = "monotone" )$x )
  else if ( number == 10 ) return( .C( "linear_pava", n = as.integer( n ), x = as.double( x ), w = as.double( w ), PACKAGE = "monotone" )$x )
  else if ( number == 11 ) return( .C( "inplace_pava", n = as.integer( n ), x = as.double( x ), w = as.double( w ), PACKAGE = "monotone" )$x )
  else if ( number == 12 ) return( .C( "md_pava", n = as.integer( n ), x = as.double( x ), w = as.double( w ), PACKAGE = "monotone" )$x )
  else if ( number == 13 ) return( .C( "reg_1d_l2", n = as.integer( n ), x = as.double( x ), w = as.double( w ), PACKAGE = "monotone" )$x )
  else if ( number == 14 ) return( .C( "jbkpava", n = as.integer( n ), x = as.double( x ), w = as.double( w ), PACKAGE = "monotone" )$x )
}
