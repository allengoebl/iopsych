# These functions use fuse() to calculate correlations between item composites.


# External functions ----------------------------------------------------------


#' Computes the correlation between two composites of items.
#' 
#' Computes the correlation between two composites of items. Composites may 
#' contain overalapping items. Items weights for each composite may be specified.
#'
#' @param r_mat A correlation matrix.
#' @param a The items used for composite A specified as a vector of column numbers.
#' @param b The items used for composite B specified as a vector of column numbers.
#' @param wt_a A vector containing the weights of each item in composite A.
#' @param wt_b A vector containing the weights of each item in composite B.
#' @return A correlation coefficient.
#'   \describe{
#'     \item{correlation}{The correlation between two composites.}
#'     \item{covariance}{The covariance between two composites.}
#'     \item{variance_a}{The variance of composite A.}
#'     \item{variance_b}{The variance of composite B.}
#'   }
#' @note This function is entended to be used for single cases. See \code{fuse2()} 
#'       for a vectorized alternative to this function. 
#' @author Allen Goebl and Jeff Jones
#' @references Lord, F.M. & Novick, M.R. (1968). \emph{Statisticl theories of 
#'             menal test scores}.
#' @examples
#' Rxx <- matrix(c(1.00, 0.25,  0.50,  0.61,
#'                 0.25, 1.00,  0.30,  0.10,
#'                 0.50, 0.30,  1.00, -0.30,
#'                 0.61, 0.10, -0.30,  1.00), 4, 4)
#' a   <- c(1, 3)
#' b   <- c(2, 4)
#' 
#' # Example using overlapping items and weights
#' Rxx  <- matrix(.3, 4, 4); diag(Rxx) <- 1
#' a    <- c(1, 2, 4)
#' b    <- c(2, 3)
#' wt_a <- c(.60, .25, .15)
#' wt_b <- c(2, 3)
#' 
#' fuse(r_mat = Rxx, a = a, b = b, wt_a = wt_a, wt_b = wt_b)
#' 
#' @export
fuse <- function(r_mat, a, b, wt_a=rep(1, length(a)), wt_b=rep(1, length(b))) {
    #Sanity check
    .isCorMat(r_mat)
    if(length(a) != length(wt_a)) {stop("wt_a must be the same length as a")}
    if(length(b) != length(wt_b)) {stop("wt_b must be the same length as b")}
    #Variances
    ab <- (wt_a %*% r_mat[a,b] %*% wt_b)
    aa <- (wt_a %*% r_mat[a,a] %*% wt_a)
    bb <- (wt_b %*% r_mat[b,b] %*% wt_b)
    #Equation
    return(ab / (sqrt(aa * bb)))
}

#' Computes the correlation between two composites of items using weights.
#' 
#' Computes the correlation between two composites of items. Composites may 
#' contain overalapping items. Items weights for each composite may be specified.
#'
#' @param r_mat A correlation matrix.
#' @param wt_a A vector containing the weights of each item in composite A. 
#'        Items which are not included in the composite should be assigned a
#'        weight of 0.
#' @param wt_b A vector containing the weights of each item in composite B. 
#'        Items which are not included in the composite should be assigned a
#'        weight of 0.
#' @return A correlation coefficient.
#' @note This is an alternative version of \code{fuse()} which uses weight vectors
#'       to specify both item selection and weight. This syntax maybe be preferable
#'       to some users. Furthermore, this function is more powerful in that it
#'       can return values for multiple sets of weights.
#' @author Allen Goebl and Jeff Jones
#' @references Lord, F.M. & Novick, M.R. (1968). \emph{Statisticl theories of 
#'             menal test scores}.
#' @examples
#' Rxx <- matrix(c(1.00, 0.25,  0.50,  0.61,
#'                 0.25, 1.00,  0.30,  0.10,
#'                 0.50, 0.30,  1.00, -0.30,
#'                 0.61, 0.10, -0.30,  1.00), 4, 4)
#' wt_a   <- c(1, 0, 1, 0)
#' wt_b   <- c(0, 1, 0, 1)
#' 
#' # Example using overlapping items and weights
#' Rxx  <- matrix(.3, 4, 4); diag(Rxx) <- 1
#' wt_a <- c(.60, .25, 0, .15)
#' wt_b <- c(0, 2, 3, 0)
#' 
#' fuse2(r_mat = Rxx, wt_a = wt_a, wt_b = wt_b)
#' 
#' @export
fuse2 <- function(r_mat, wt_a, wt_b) {
    #Sanity check
    .isCorMat(r_mat)
    if (is.vector(wt_a)) {dim(wt_a) <- c(1, length(wt_a))}
    if (is.vector(wt_b)) {dim(wt_b) <- c(1, length(wt_b))}
    #Variances
    ab <- (wt_a %*% r_mat %*% t(wt_b))
    aa <- diag(wt_a %*% r_mat %*% t(wt_a))
    bb <- diag(wt_b %*% r_mat %*% t(wt_b))
    #Equation
    return(ab / (sqrt(aa %o% bb)))
}

#' Computes the correlations between a correlation matrix and a weighted composite
#' of items from the matrix.
#' 
#' @param r_mat A correlation matrix.
#' @param wt A vector containing the weights of each item in composite A or a
#'        matrix with one row per weight vector.
#' @return A vector of correlation coefficients will be returned if \code{wt_a}
#'         is a vector. If /code{wt_b} is a matrix, a matrix of correlation 
#'         coefficients with one row for each weight vector will be returned.
#' @author Allen Goebl and Jeff Jones
#' @references Lord, F.M. & Novick, M.R. (1968). \emph{Statisticl theories of 
#'             mental test scores}.
#' @examples
#' data(dls2007)
#' dat <- dls2007
#' rxx <- dat[1:4, 2:5]
#' wt1 <- c(1, 1, 1, 1)
#' wt2 <- c(2, 0, 1, 0)
#' wt  <- rbind(wt1, wt2)
#' 
#' fuseVec(r_mat=rxx, wt=wt1)
#' @export
fuseVec <- function(r_mat, wt) {
    #Sanity check
    .isCorMat(r_mat)
    if (is.vector(wt)) {dim(wt) <- c(1, length(wt))}
    #Variances
    ab <- wt %*% r_mat
    aa <- diag(wt %*% r_mat %*% t(wt))
    #Equation
    return(ab / (sqrt(aa)))
}

#' The intercorrelation among items and composites made of these items.
#' 
#' The key matrix is used to specify any number of weighted item composites.
#' A correlation matrix of these composites and the original correlation matrix
#' is then computed and returned.
#'
#' @param r_mat A correlation matrix.
#' @param wt A matrix with one row for each composite and one column for 
#'        each item contained in \code{r_mat}. The value if each element corresponds
#'        to the weight given to an item.
#' @param type The type of output desired.
#' @return If \code{type} = "cxc" then a matrix of the intercorrelations between the specified 
#'         composites are returned. 
#'         If \code{type} = "full" then all the intercorrelations between both the
#'         original items and the specified composites are returned.
#' @author Allen Goebl and Jeff Jones
#' @examples
#' Rxx <- matrix(c(1.00, 0.25,  0.50,  0.61,
#'                 0.25, 1.00,  0.30,  0.10,
#'                 0.50, 0.30,  1.00, -0.30,
#'                 0.61, 0.10, -0.30,  1.00), 4, 4); Rxx
#' 
#' # Single composite
#' wt <- matrix(c(1, 2, 3, -1), 1, 4); wt
#' 
#' fuseMat(r_mat = Rxx, wt = wt)
#' 
#' # Three composites
#' wt  <- matrix(c(1, 2, 3, -1,
#'                 2, 1, 0, -2,
#'                 1, 1, 0,  0), 3, 4, byrow = TRUE)
#' 
#' fuseMat(Rxx, wt)
#' @export
fuseMat <- function(r_mat, wt, type="full"){
    #Check input
    .isCorMat(r_mat)
    if(nrow(r_mat) != ncol(wt)) { stop("key_mat does not match r_mat") }
    #Fuse
    if(type == "full") { wt <- rbind(diag(nrow(r_mat)), wt) }
    v <- (wt %*% r_mat %*% t(wt))
    diagv <- diag(v)
    return(v / (sqrt(diagv %o% diagv)))
}

