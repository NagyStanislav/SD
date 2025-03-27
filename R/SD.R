SD2m = function(x){
  return(cl(x))
}

#### SD2 ####

#' \code{K}-hull depth for bivariate data
#'
#' C++ and R implementations of the algorithm for
#' the computation of the \code{K}-hull depth for
#' two-dimensional data. The C++ implementation
#' is the fastest one that should be used in practice; its
#' complexity is O(n*log(n)). For comparison,
#' also the (slower) R version of the same complexity
#' is provided. Compared to function 
#' \link[ddalpha:depth.simplicial]{depth.simplicial} from 
#' package \code{ddalpha}, this function works also for higer
#' sample sizes (the function \link[ddalpha:depth.simplicial]{depth.simplicial})
#' overflows and gives negative results with higher sample sizes).
#'
#' @param data A matrix of dimension \code{n}-times-\code{2}, where
#' \code{n} is the sample size. The depth of all the points in this matrix
#' is computed, with respect to the remaining \code{n-1} points.
#'
#' @param k Positive integer, parameter of the \code{k}-hull depth. By
#' default set to \code{k=3}, which gives the standard bivariate simplicial 
#' depth.
#'
#' @param normalize Logical indicator whether the resulting depth should
#' be given as an integer value (\code{normalize=FALSE}, the default value), or
#' whether this number should be divided by \code{choose(n-1,k)} to obtain
#' the standardized version of the depth between \code{0} and \code{1}.
#' 
#' @return A vector of length \code{n} of the resulting \code{k}-hull depths.
#'
#' @references Stanislav Nagy, Martin Wendler, and Carsten Jentsch (2025). 
#' \emph{Resampling simplicial depth.} Under review.
#' 
#' @references Erik Mendros and Stanislav Nagy (2025). 
#' \emph{k-hull depth: In between simplicial and halfspace depth.} Statistical 
#' Papers. To appear.
#'
#' @examples
#' n = 165
#' data = matrix(rnorm(n*2),ncol=2)
#' 
#' res = SD2(data,3)
#' i = 1
#' res[i]
#' ddalpha::depth.simplicial(data[i,],data[-i,],exact=TRUE)*choose(n-1,3)

SD2 = function(data, k=3, normalize=FALSE, vrs=c("C","R")){
  n = nrow(data)
  if(normalize) cst = choose(n-1,k) else cst = 1
  vrs = match.arg(vrs)
  if(vrs=="C") return(SD2C(data,n,k)/cst)
  if(vrs=="R"){
    data2 = t(data)
    return(apply(data2,2,function (x){
      data2 = data2 - x
      xt = data2[1,order(data2[2,] / data2[1,])]
      SD::cl(sign(xt[xt!=0]))
    })/cst)    
  }
}

# n = 1654
# data = matrix(rnorm(n*2),ncol=2)
# 
# i = 1
# res[i]
# SD2b(data)[i]
# ddalpha::depth.simplicial(data[i,],data[-i,],exact=TRUE)*choose(n-1,3)
# 
# (bm = bench::mark(
#   SD2C(data,n,3),
#   SD2b(data),
#   ddalpha::depth.simplicial(data,data,exact=TRUE)*choose(n-1,3),
#   check=FALSE
# ))
# 
# plot(bm)