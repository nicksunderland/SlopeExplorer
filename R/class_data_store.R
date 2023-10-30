# define a data stricuture for each iteration
DataStore <- function(points=NULL, pt_data=NULL, sx0=NULL, sy0=NULL, sxy0=NULL, dir0=NULL, sx1=NULL, sy1=NULL, sxy1=NULL, pi0=NULL) {
  structure(list(
    points = points,
    pt_data = pt_data,
    sx0 = sx0,
    sy0 = sy0,
    sxy0 = sxy0,
    dir0 = dir0,
    sx1 = sx1,
    sy1 = sy1,
    sxy1 = sxy1,
    pi0 = pi0
  ), class = "DataStore")
}
