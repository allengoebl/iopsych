# This file contains functions for converting between r and d.

utils::globalVariables(c("uc", "sr_majority", "sr_minority", "ai", "rxx", "ryy"))


#' Convert from r to d
#' 
#' @param r A r-value or a vector of r values.
#' @return A d value or a vector of d values.
#' @author Allen Goebl and Jeff Jones
#' @examples
#' cor2d(.3)
#' cor2d(((1:9)/10))
#' @export
cor2d <- function (r) ((2 * r) / sqrt(1 - (r ^ 2)))

#' Convert from d to r
#' 
#' @param d A d-value or a vector of d values.
#' @return A r value or a vector of r values.
#' @author Allen Goebl and Jeff Jones
#' @examples
#' d2cor(.3)
#' d2cor(((1:9)))
#' @export
d2cor <- function (d) (sqrt((d ^ 2) / (4 + (d ^ 2))))

#' Estimates the d of a composite.
#'
#' @param rxx A matrix of predictor intercorrelations. 
#' @param d_vec A vector containing d's for each predictor.
#' @param wt_vec A vector containing the weights of each item in rxx.
#' @return A vector of correlation coefficients.
#' @author Jeff Jones and Allen Goebl
#' @references Sackett, P. R., & Ellingson, J. E. (1997). \emph{Personnel 
#' Psychology., 50(3)}, 707-721.
#' @note This is essentially the same function as solveWt().
#' @examples
#' Rxx <- matrix(.3, 3, 3); diag(Rxx) <- 1
#' ds  <- c(.2, .4, .3)
#' dComposite(rxx = Rxx, d_vec = ds)
#' 
#' Rxx <- matrix(c(1.0, 0.3, 0.2, 
#'                 0.3, 1.0, 0.1,
#'                 0.2, 0.1, 1.0), 3, 3)
#' ds  <- c(.1, .3, .7)
#' ws  <- c(1, .5, .5)
#' dComposite(rxx = Rxx, d_vec = ds, wt_vec = ws)
#' @export
dComposite <- function(rxx, d_vec, wt_vec=rep(1, length(d_vec))) {
    numer <- t(wt_vec) %*% d_vec
    denom <- sqrt(t(wt_vec) %*% rxx %*% wt_vec)
    return(numer/denom)
}

#' Estimate adverse impact given d and sr
#'
#' @param d Subgroup difference. 
#' @param sr The percentage of the applicant population that is selected.
#' @param pct_minority The percentage of the applicant population that is part of
#'        a given minority group.
#' @return 
#'   \describe{
#'     \item{ai}{The adverse impact ratio.}
#'     \item{overall_sr}{The overall selection ratio.}
#'     \item{sr_majority}{The selection ratio for the majority group.}
#'     \item{sr_minority}{The selection ratio for the minority group.}
#'     \item{uc}{the predictor cutoff value that corresponds to the given overall selection ratio}
#'   }
#' @author Jeff Jones and Allen Goebl
#' @references De Corte, W., Lievens, F.(2003). A practical procedure to estimate
#' the quality and the adverse impact of single-stage selection decisions.
#' \emph{International Journal of Selection and Assessment, 11(1)}, 87-95.
#' @examples
#' aiEst(d = -0.15, sr = 0.25, pct_minority = 0.30)
#' 
#' aiEst(d = -0.40, sr = 0.10, pct_minority = 0.15)
#' @export
#' @importFrom stats pnorm
#' @importFrom stats uniroot
aiEst <- function(d, sr, pct_minority){
    pct_majority <- (1 - pct_minority)
    #Function to be minimized
    stRoot <- function(x) {
        pct_majority*(1 - pnorm(x)) + pct_minority*(1 - pnorm(x - d)) - sr
    }
    #Selected
    uc <- uniroot(stRoot, interval=c(-3,3))$root
    sr_majority <- (1 - pnorm(uc))
    sr_minority <- (1 - pnorm(uc - d))
    ai <- (sr_minority / sr_majority)
    out <- list(ai, sr, sr_majority, sr_minority, uc)
    names(out) <- c("ai","overall_sr", "sr_majority", "sr_minority", "uc")
    return(out)
}

#' Estimate ai and average criterion scores for majority and minority groups.
#'
#' @param mr The correlation between the predictor and criterion composites.
#' @param dx A vector of d values for the predictors. These d values are expected
#'           to have been computed in the direction of Minority - Majority.
#' @param dy A vector of d values for the criteria These d values are expected
<<<<<<< 0064e9a28fc3f3cb2c53edbb97ff99cc5cbc1fcb
#'           to have been computed in the direction of Minority - Majority
#' @param sr The percentage of the applicant population who are selected.
#' @param pct_minority The percentage of the applicant population who are part of
=======
#'           to have been computed in the direction of Majority - Minority.
#' @param sr The percentage of the applicant population that is selected.
#' @param pct_minority The percentage of the applicant population that is part of
>>>>>>> FO changes
#'        a given minority group.
#' @return \describe{
#'             \item{AI}{Adverse Impact}
#'             \item{Overeall_sr}{The overall selection ratio set by the user}
#'             \item{Majority_sr}{Majority Selection Rate}
#'             \item{Minority_sr}{Minority Selection Rate}
#'             \item{Majority_Standardized}{Predicted composite criterion score relative to the majority population}
#'             \item{Global_Standardized}{Predicted composite criterion score relative to the overall population} 
#'          }
#' @author Jeff Jones and Allen Goebl
#' @references De Corte, W., Lievens, F.(2003). A practical procedure to estimate
#' the quality and the adverse impact of single-stage selection decisions.
#' \emph{International Journal of Selection and Assessment.}, 11(1), 87-95.
#' @examples
#' aiPux(0.6, dx=-0.8, sr=0.3, pct_minority=0.25)
#' aiPux(0.6, dx=-0.8, dy=-0.2, sr=0.3, pct_minority=0.25)
#' @export
aiPux <- function(mr, dx, dy = -1, sr, pct_minority) {
    list2env(x=aiEst(dx, sr, pct_minority), envir=environment())
    pct_majority <- 1 - pct_minority; pct_majority
    #Expected criterion score of applicant groups (relative to majority) 
    meanZi <-  dy + mr * dnorm(uc - dx)/(1 - pnorm(uc - dx))
    meanZa <-  mr * dnorm(uc) / (1 - pnorm(uc))
    meanZt <- (pct_majority*sr_majority*meanZa+pct_minority*sr_minority*meanZi)/sr
    majority_std <- matrix(c(meanZi, meanZa, meanZt), 3, 1)
    rownames(majority_std) <- c('Zi', 'Za', 'Zt')
    #Expected criterion score of applicant groups (relative to the total group)
    sd_global <- sqrt(1 + pct_minority * pct_majority * dy^2)
    meanZig <- (meanZi - pct_minority * dy) / sd_global
    meanZag <- (meanZa - pct_minority * dy) / sd_global
    meanZtg <- (meanZt - pct_minority * dy) / sd_global
    global_std <- matrix(c(meanZig, meanZag, meanZtg), 3, 1)
    #Format Output
    rownames(global_std) <- c('Zi', 'Za', 'Zt')
    out <- list(ai, sr, sr_majority, sr_minority, majority_std, global_std)
    names(out) <- c("AI","Overall_sr", "Majority_sr", "Minority_sr", 
                    "Majority_Standardized", "Global_Standardized")
    return(out)
}

#' Estimate ai and average criterion scores for majority and minority groups.
#'
#' @param r_mat Super correlation matrix between the predictors and criteria.
#'        This argument assumes that the predictors come first in the matrix. 
#' @param y_col A vector of columns representing criterion variables.
#' @param x_col A vector of columns representing predictor variables.
#' @param dX A vector of d values for the predictors. These d values are expected
#'           to have been computed in the direction of Minority - Majority.
#' @param dY A vector of d values for the criteria These d values are expected
#'           to have been computed in the direction of Minority - Majority.
#' @param wt_x Weights for the predictors to form the overall composite predictor.
#' @param wt_y Weights for the criteria to form the overall composite criterion.
#' @param sr The percentage of the applicant population that is selected.
#' @param pct_minority The percentage of the applicant population that is part of
#'        a given minority group.
#' @return \describe{
#'             \item{AI}{Adverse Impact}
#'             \item{Overeall_sr}{The overall selection ratio set by the user}
#'             \item{Majority_sr}{Majority Selection Rate}
#'             \item{Minority_sr}{Minority Selection Rate}
#'             \item{Majority_Standardized}{Predicted composite criterion score relative to the majority population}
#'             \item{Global_Standardized}{Predicted composite criterion score relative to the overall population} 
#'          }
#' @author Jeff Jones and Allen Goebl
#' @references De Corte, W., Lievens, F.(2003). A Practical procedure to estimate
#' the quality and the adverse impact of single-stage selection decisions.
#' \emph{International Journal of Selection and Assessment, 11(1)}, 87-95.
#' @references De Corte, W. (2003). Caiqs user's guide. http://allserv.rug.ac.be/~wdecorte/software.html
#' @examples
#' # Example taken from De Corte, W. (2003)
#' R <- matrix(c(1.000, 0.170, 0.000, 0.100, 0.290, 0.160, 
#'               0.170, 1.000, 0.120, 0.160, 0.300, 0.260, 
#'               0.000, 0.120, 1.000, 0.470, 0.120, 0.200, 
#'               0.100, 0.160, 0.470, 1.000, 0.240, 0.250, 
#'               0.290, 0.300, 0.120, 0.240, 1.000, 0.170, 
#'               0.160, 0.260, 0.200, 0.250, 0.170, 1.000), 6, 6)
#'
#' wt_x <- c(.244, .270, .039, .206) 
#' wt_y <- c(6, 2)
#' sr    <- 0.25
#' pct_minority <- .20
#'
#' dX   <- c(-1, -0.09, -0.09, -0.20)
#' dY   <- c(-0.450, 0.0)
#'
#' aiPuxComposite(R, 5:6, 1:4, dX, dY, wt_x, wt_y, sr, pct_minority)
#' 
#' # compare the output from predictAI with the output in the CAIQS manual on page 7 where SR = .250
#' 
#' @export
aiPuxComposite <- function(r_mat, y_col, x_col, dX, dY, wt_x, wt_y, sr, pct_minority){
    #Define rxx, ryy, rxy
    list2env(x=.indexMat(r_mat, y_col, x_col), envir=environment())
    #Compute subgroup differences in the predictor and criterion composites
    dx <- dComposite(rxx=rxx, d_vec=dX, wt_vec=wt_x)
    dy <- dComposite(rxx=ryy, d_vec=dY, wt_vec=wt_y)
    #Computer correlation between predictor and criterion composites.
    mr <- fuse(r_mat=r_mat, a=x_col, b=y_col, wt_a=wt_x, wt_b=wt_y)
    return(aiPux(mr=mr, dx=dx, dy=dy, sr=sr, pct_minority=pct_minority))
}   

#' The average score of selected applicants on a predictor composite.
#'
#' When scores on the predictor composite are assumed to be normally
#' distributed, the average score of selected applicants can be computed for 
#' an arbitrary selection ratio using the ordinate of the normal curve.
#'
#' @param sr A selection ratio or a vector of selection ratios.
#' @return ux: The average score of those selected on a predicter composite.
#' @author Allen Goebl and Jeff Jones
#' @references Naylor, J. C., & Shine, L. C. (1965). A table for determining the
#'  increase in mean criterion score obtained by using a selection device. 
#'  \emph{Journal of Industrial Psychology}, 78-109.
#' @examples
#' ux(.6)
#' @export
#' @importFrom stats dnorm
#' @importFrom stats qnorm
ux <- function(sr) (dnorm(qnorm(1-sr)) / sr)





