#' Bayes  risk optimal  testing
#' 
#' Carries  out   the  Bayes   risk  optimal  test   of  the   composite  null
#' "\eqn{\delta_x    \times     \delta_y=0}"    against     its    alternative
#' "\eqn{\delta_x \times \delta_y\neq  0}" based on the test  statistic in the
#' real plane.
#' 
#' @param t  A  \code{vector}  consisting of  two  \code{numeric}s, the  test
#'   statistic in  the real  plane, or a  'n x 2'  \code{matrix} of  such test
#'   statistics.
#'
#' @param map  The "map" of rejection  probabilities -- the 'map'  item of the
#'     output of function \code{compute_map_rejection_probs}.
#'
#' @param truncation A nonnegative  \code{numeric} used to bound the rejection
#'     region away  from the null  hypothesis space.  Defaults to 0,  in which
#'     case the  rejection region is minimax optimal.
#' 
#' @details For details, we refer  to the technical report  "Optimal Tests of
#'     the Composite Null Hypothesis Arising in Mediation Analysis", by
#'     Miles & Chambaz (2024), https://arxiv.org/abs/2107.07575 
#' 
#' @return A list, consisting  of: \describe{ \item{t:}{a \code{vector} of two
#'     \code{numeric}s, the test statistic, or a 'n x 2' \code{matrix} of such
#'     test  statistics;} \item{alpha:}{a  \code{numeric}, the  type-I error;}
#'     \item{truncation:}{a  nonnegative  \code{numeric},  used to  bound  the
#'     rejection    region   away    from   the    null   hypothesis    space}
#'     \item{decision:}{a  \code{vector} of  \code{logical}s, \code{FALSE}  if
#'     the  null hypothesis  can  be  rejected for  the  alternative at  level
#'     'alpha'  and \code{TRUE}  otherwise;}  \item{pval:}{a \code{vector}  of
#'     \code{numeric}s,  the  p-values  of  the tests,  'NA'  in  this  case;}
#'     \item{method:}{the  \code{character} "Bayes";}\item{map:}{The  "map" of
#'     rejection probabilities  -- the  'map' item of  the output  of function
#'     \code{compute_map_rejection_probs}.} }
#'
#' @examples
#' n <- 10
#' x <- MASS::mvrnorm(2 * n, mu = c(0, 0), Sigma = diag(c(1, 1)))
#' delta <- matrix(stats::runif(4 * n, min = -3, max = 3), ncol = 2)
#' epsilon <- stats::rbinom(n, 1, 1/2)
#' delta <- delta * cbind(c(epsilon, rep(1, n)),
#'                        c(1 - epsilon, rep(1, n)))
#' x <- x + delta
#' (mt_01_0.05 <- mediation_test_Bayes(x, map = map_01_0.05))
#' plot(mt_01_0.05)
#' (mt_quad_0.05_0.1 <- mediation_test_Bayes(x, map = map_quad_0.05_0.1))
#' plot(mt_quad_0.05_0.1)
#'
#' @export
mediation_test_Bayes <- function(t, map, truncation = 0) {
    t <- R.utils::Arguments$getNumerics(t)
    if (!inherits(map, "map.rejection.probabilities")) {
        R.oo::throw("Argument 'map' should be an output of function 'compute_map_rejection_probs', not ", map) 
    }
    truncation <- R.utils::Arguments$getNumeric(truncation, c(0, Inf))
    if (is.vector(t)) {
        t <- matrix(t, ncol = 2)
    }
    if (!ncol(t) == 2) {
        R.oo::throw("Argument 't' should be a vector consisting of two real numbers or a n x 2 matrix or data frame (the test statistic(s) in the real plane), not ", t) 
    }
    alpha <- attr(map, "alpha")
    ##
    make_decision <- function(tt, aalpha, mmap) {
        K <- ncol(mmap)/2
        B <- 2 * stats::qnorm(1 - aalpha/2)
        abstt <- abs(tt)
        max_abstt <- pmax(abstt[, 1], abstt[, 2])
        min_abstt <- pmin(abstt[, 1], abstt[, 2])
        ##
        decision <- (max_abstt > B) & (min_abstt > B/2)
        idx <- which(max_abstt <= B)
        if (length(idx) >= 1) {
            ij <- cbind(ceiling((tt[idx, 1] + B) * K/B),
                        ceiling((tt[idx, 2] + B) * K/B))
            decision[idx] <- as.logical(stats::rbinom(length(idx), 1, map[ij]))
        }
        names(decision) <- rep("rejection", length(decision))
        return(decision)
    }
    compute_pval <- function(tt) {
        NA
    }
    decision <- make_decision(t, alpha, map)
    if (truncation > 0) {
        decision <- decision & (apply(abs(t), 1, min) >= truncation)
    }
    pvals <- compute_pval(t)
    
    out <- list(t = t,
                alpha = alpha,
                truncation = truncation,
                decision = decision,
                pvals = pvals,
                method = "Bayes",
                map = map)
    class(out) <- "mediation.test"
    return(out)
}

