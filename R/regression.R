# This file contains all the functions for calculating R2 and beta weights
utils::globalVariables(c("rxy", "rxx", "ryy"))


### Internal Functions --------------------------------------------------------


#' Correlation between weighted predictor composite and criterion.
#' 
#' @param rxx A matrix of predictor intercorrelations.
#' @param rxy A vector of predictor criterion intercorrelations, or a
#'        matrix with one row per criterion.
#' @param wt A vector of predictor weights, or a matrix with one set of 
#'        predictor weights per row. 
#' @return A matrix of correlation coefficent with one row per weight vector
#'         and one column per rxy vector.
#' @author Allen Goebl Jeff Jones
#' @examples
#' library(iopsych)
#' data(dls2007)
#' dat <- dls2007[1:6, 2:7]
#' rxx <- dat[1:4, 1:4]
#' rxy <- dat[5:6, 1:4]
#' 
#' wt1 <- c(1,1,1,1)
#' wt2 <- c(1,2,3,4)
#' wt_mat <- rbind(wt1, wt2)
#'
#' solveWtPred(rxx=rxx, rxy=rxy, wt=wt_mat)
#' @export
solveWtPred <- function(rxx, rxy, wt) {
        if (is.vector(wt)) {dim(wt) <- c(1, length(wt))}
        if (is.vector(rxy)) {dim(rxy) <- c(1, length(rxy))}
        numer <- tcrossprod(wt, rxy)
        denom <- sqrt(diag(wt %*% rxx %*% t(wt)))
        return(numer / denom)
}

#' Correlation between weighted criterion composite and predictors.
#' 
#' @param ryy A matrix of criterion intercorrelations.
#' @param rxy A vector of predictor criterion intercorrelations, or a
#'        matrix with one row per criterion.
#' @param wt A vector of criterion weights, or a matrix with one set of 
#'        criterion weights per row. 
#' @return A matrix of correlation coefficent with one row per weight vector
#'         and one column per predictor.
#' @author Allen Goebl Jeff Jones
#' @examples
#' library(iopsych)
#' data(dls2007)
#' dat <- dls2007[1:6, 2:7]
#' ryy <- dat[5:6, 5:6]
#' rxy <- dat[5:6, 1:4]
#' 
#' wt1 <- c(.25, .75)
#' wt2 <- c(.75, .25)
#' wt_mat <- rbind(wt1, wt2)
#'
#' solveWtCrit(ryy=ryy, rxy=rxy, wt=wt_mat)
#' @export
solveWtCrit <- function(ryy, rxy, wt) {
        if (is.vector(wt)) {dim(wt) <- c(1, length(wt))}
        if (is.vector(rxy)) {dim(rxy) <- c(1, length(rxy))}
        numer <- (wt %*% rxy)
        denom <- sqrt(diag(wt %*% ryy %*% t(wt)))
        return(numer / denom)
}

#SolveWtPredCrit?

#' Find beta weights
#' 
#' @param rxx A matrix of predictor intercorrelations.
#' @param rxy A vector of predictor criterion intercorrelations, or a
#'        matrix with one row per criterion.
#' @return A vector or matrix of regression weights.
#' @author Allen Goebl and Jeff Jones
#' @examples
#' library(iopsych)
#' data(dls2007)
#' dat <- dls2007[1:6, 2:7]
#' rxx <- dat[1:4, 1:4]
#' rxy <- dat[1:4, 5]
#'
#' solveBeta(rxx=rxx, rxy=rxy)
#' @export
solveBeta <- function(rxx, rxy) tcrossprod(rxy, chol2inv(chol(rxx)))

#' Find R2
#' 
#' @param rxx A matrix of predictor intercorrelations.
#' @param rxy rxy A vector of predictor criterion intercorrelations, or a
#'        matrix with one row per criterion.
#' @return R2 and Regression weights
#' @author Allen Goebl and Jeff Jones
#' @examples
#' print("example needed")
#' @export
solveR2 <- function(rxx, rxy) {
    if (is.vector(rxy)) {dim(rxy) <- c(1, length(rxy))}
    beta <- solveBeta(rxx, rxy)
    return(diag(tcrossprod(beta, rxy)))
}

#' Find beta weights and R2
#' 
#' @param rxx A matrix of predictor intercorrelations.
#' @param rxy A vector of predictor criterion intercorrelations, or a
#'        matrix with one row per criterion.
#' @return R2 and Regression weights
#' @author Allen Goebl and Jeff Jones
#' @examples
#' print("example needed")
#' @export
solveReg <- function(rxx, rxy) {
    if (is.vector(rxy)) {dim(rxy) <- c(1, length(rxy))}
    beta <- solveBeta(rxx, rxy)
    R2 <- diag(tcrossprod(beta, rxy))
    out <- (list(R2, beta))
    names(out) <- c("R2", "beta")
    return(out)
}


# External functions ----------------------------------------------------------


#' Regression
#' 
#' @param r_mat A correlation matrix.
#' @param y_col The column representing the criterion variable.
#' @param x_col A vector of columns representing predictor variables.
#' @param N Number of observations
#' @param alpha alpha value for (1 - alpha)\% Confidence Interval
#' @return Regression beta weights and \eqn{R^2}. If N is supplied, the standard error of 
#' the beta weights as well as the confidence intervals are returned as well.
#' @author Allen Goebl and Jeff Jones
#' @note If N is non-null the function will compute corrected standard errors
#' for the standardized regression coefficients using the delta method. For additional 
#' details on the calculation of the corrected standard errors see Yuan and Chan (2011) 
#' and Jones and Waller (2015).
#' @references Jones, J. A. & Waller, N. G. (2015). The normal-theory and asymptotic 
#' distribution-free covariance matrix of standardized regression coefficients: 
#' Theoretical extensions and finite sample behavior. \emph{Psychometrika, 80}, 365-378.
#' @references Yuan, K. and Chan, W. (2011). Biases and standard errors of
#' standardized regression coefficients. \emph{Psychometrika, 76}, 670-690.
#' @examples
#' Rs <- matrix(c(1.0, 0.2,  0.3, 0.4, -0.4,
#'                0.2, 1.0,  0.5, 0.1,  0.1,
#'                0.3, 0.5,  1.0, 0.2, -0.3,
#'                0.4, 0.1,  0.2, 1.0,  0.4,
#'               -0.4, 0.1, -0.3, 0.4,  1.0), 5, 5)
#' ys <- 5
#' xs <- 1:4
#'
#' rmatReg(Rs, ys, xs)
#' 
#' # Example with standard errors
#' rmatReg(Rs, ys, xs, N = 100)
#' @export
#' @importFrom utils capture.output
rmatReg <- function(r_mat, y_col, x_col, N = NULL, alpha = 0.05) {
    .checkIndex(r_mat=r_mat, y_col=y_col, x_col=x_col)
    list2env(x=.indexMat(r_mat, y_col, x_col), envir=environment())
    if(!is.null(N)) {
        invisible(capture.output(seOut <- fungible::seBetaCor(R=rxx, rxy=rxy, Nobs=N)))
    }
    else seOut <- NULL    
    out <- solveReg(rxx=rxx, rxy=rxy)
    #Format Output
      #NOTE: Maybe use class here
    return(list(R2 = out$R2, beta = out$beta, se_beta = seOut$se.Beta, ci_beta = seOut$CI.beta))
}

#' Find \eqn{R^2} given arbitrary predictor weights
#' 
#' @param r_mat A correlation matrix.
#' @param y_col A vector of columns representing criterion variables.
#' @param x_col A vector of columns representing predictor variables.
#' @param wt A vector of predictor weights or a list of multiple vectors.
#' @return Regression R2.
#' @note This is a wrapper for solveWt().
#' @author Allen Goebl and Jeff Jones
#' @examples
#' library(iopsych)
#' #Get Data
#' data(dls2007)
#' r_mat <- dls2007[1:6, 2:7]
#' 
#' #Get weights
#' unit_wt <- c(1,1,1,1)
#' other_wt <- c(1,2,1,.5)
#' wt_mat <- rbind(unit_wt, other_wt)
#'
#' #Solve
#' rmatWtR2(r_mat=r_mat, y_col=6, x_col=1:4, wt=unit_wt)
#' rmatWtR2(r_mat=r_mat, y_col=6, x_col=1:4, wt=other_wt)
#' rmatWtR2(r_mat=r_mat, y_col=6, x_col=1:4, wt=wt_mat)
#' @export
rmatWtR2 <- function(r_mat, y_col, x_col, wt) {
    list2env(x=.indexMat(r_mat, y_col, x_col), envir=environment()) 
    r <- solveWtPred(rxx=rxx, rxy=rxy, wt=wt)
    r2 <- (r^2)
    return(r2)
}

