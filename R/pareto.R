# These functions calaculate pareto tradeoffs
utils::globalVariables(c("rxy", "rxx", "ryy"))


# Internal ---------------------------------------------------------------------


#' An internal function for generating criterion weights for XX pareto plots.
#' 
#' @param pts How many different pairs of criterion weights to calculate.
#' @return A matrix of of criterion weights with one column per criteiron. 
#' @author Allen Goebl and Jeff Jones
#' @examples
#' print("example needed")
#' @keywords internal
#' @rdname internal.wtPair
.wtPair <- function(pts) {
    wt_one <- seq(from=0, to=1, by=(1/(pts-1)))
    wt_two <- (1-wt_one)
    return(cbind(wt_one, wt_two))
}


# External --------------------------------------------------------------------


#' Computes data needed for a XX Pareto plot.
#' 
#' @param r_mat A correlation matrix.
#' @param y_col A vector of columns representing criterion variables.
#' @param x_col A vector of columns representing predictor variables.
#' @param pts The number of points used. Determines accuracy.
#' @return  
#'   \describe{
#'      \item{betas}{A matrix of beta weights for each criteria weight}
#'      \item{wt_one}{The weight given to the first criterion}
#'      \item{multiple_r}{The correlation between the predictor and criterion composites}
#'    }
#' @author Allen Goebl and Jeff Jones
#' @examples
#' # Setup Data
#' data(dls2007)
#' r_mat <- dls2007[1:6, 2:7]
#'
#' #Run Model
#' XX1 <- paretoXX(r_mat=r_mat, x_col=1:4, y_col=5:6)
#' # Plot Multiple correlations
#' plot(c(0,1), c(.3,.5), type="n", xlab="C1 Wt", ylab="mr") 
#' lines((XX1$wts)[,1], (XX1$multiple_r)[,1])
#' lines((XX1$wts)[,1], (XX1$multiple_r)[,2])
#' @export
paretoXX <- function(r_mat, x_col, y_col, pts=100) {
    #Check Input
    .isCorMat(r_mat)
    if(length(y_col) != 2) {stop("y_col should be a vector of length 2")}
    if(max(c(x_col)) > ncol(r_mat)) {stop("x_col is out of bounds")}
    if(max(c(y_col)) > ncol(r_mat)) {stop("y_col is out of bounds")}
    #Terms
    list2env(x=.indexMat(r_mat, y_col, x_col), envir=environment())
    #Solve
    wts <- .wtPair(pts)
    rxy_composite <- solveWtCrit(ryy=ryy, rxy=rxy, wt=wts)
    betas <- solveBeta(rxx=rxx, rxy=rxy_composite)
    mr <- solveWtPred(rxx=rxx, rxy=rxy, wt=betas)
    #Format output
    out <- list(betas, wts,  mr)
    names(out) <- c("betas", "wts", "multiple_r")
    colnames(out[["betas"]]) <- colnames(r_mat)[x_col]
    return(out)
}


#Pareto XY --------------------------------------------------------------------


#' Computes data needed for a XY Pareto plot.
#'
#' @param r_mat A correlation matrix.
#' @param y_col A vector of columns representing criterion variables.
#' @param x_col A vector of columns representing predictor variables.
#' @param d_vec A vector of d scores.
#' @param gen The number of iterations used by the algorithim.
#' @param pop The population or number of cases used by the algorithim.
#' @param pred_lower The minimum weight allowed for each predictor. 
#' @param pred_upper The maximum weight allowed for each predictor.
#' @return  
#'   \describe{
#'     \item{betas}{A matrix of beta weights for each criteria weight}
#'     \item{mr_d}{A matrix of multiple correlations or d values 
#'                 corresponding to each row of beta weights.}
#'     \item{pareto_optimal}{A vector indicating whether each value is
#'           pareto optimal}
#'   }
#' @author Allen Goebl Jeff Jones
#' @examples
#' data(dls2007)
#' dat <- dls2007
#' r_mat <- dat[1:6, 2:7]
#' x_col <- 1:4 
#' y_col <- 5:6
#' d_vec <- -dat[1:4, 1]
#'
#' paretoXY(r_mat=r_mat, x_col=1:4, y_col=5, d_vec=d_vec, pred_lower=c(0,0,0,0))
#' paretoXY(r_mat=r_mat, x_col=1:4, y_col=c(5,6))
#' @export
#Pareto function
paretoXY <- function(r_mat, x_col, y_col, d_vec=NULL, gen=100, pop=100, 
                     pred_lower=rep(-2, length(x_col)),
                     pred_upper=rep( 2, length(x_col)), ...) {
  #Terms
  rxx <- r_mat[x_col, x_col]
  rxy <- rbind(r_mat[y_col, x_col], d_vec)
  #Set Objective function
  obj <- function(a) -t(solveWtPred(rxx=rxx, rxy=rxy, wt=a))
  #Optimize Objective
  out <- mco::nsga2(fn           = obj, 
                    idim         = length(x_col), 
                    odim         = nrow(rxy), 
                    generations  = gen, 
                    popsize      = pop,
                    lower.bounds = pred_lower,
                    upper.bounds = pred_upper,
                    vectorized   = TRUE,
                    ...)
  #Rename output -- More complex output needed if gen > 1
  if(!is.null(names(out))) {
    names(out) <- c("betas", "mr_d", "pareto_optimal")
    out$mr_d <- -out$mr_d
    return(out)
  } else {
    for (i in 1:length(out)){
      names(out[[i]]) <- c("betas", "mr_d", "pareto_optimal")
      out[[i]]$mr_d <- -out[[i]]$mr_d
    }
  }
  out
}

