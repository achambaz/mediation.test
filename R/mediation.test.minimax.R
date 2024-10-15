#' Carries   out   the   minimax   optimal  test   of   the   composite   null
#' "\eqn{\delta_x    \times     \delta_y=0}"    against     its    alternative
#' "\eqn{\delta_x \times \delta_y\neq  0}" based on the test  statistic in the
#' real plane.
#' 
#' @param t  A  \code{vector}  consisting of  two  \code{numeric}s, the  test
#'   statistic in  the real  plane, or a  'n x 2'  \code{matrix} of  such test
#'   statistics.
#'
#' @param alpha A positive \code{numeric}, the wished type-I error.
#'
#' @param truncation A nonnegative  \code{numeric} used to bound the rejection
#'     region away  from the null  hypothesis space.  Defaults to 0,  in which
#'     case the  rejection region is minimax optimal.
#' 
#' @param sample_size An \code{integer}  (larger than  one), the size  of the
#'     sample used  to derive the  test statistic. Defaults to  'Inf', meaning
#'     that, under the  null hypothesis, the test statistic is  drawn from the
#'     \eqn{N_2(0,I_2)} law.  If  the integer is finite, then,  under the null
#'     hypothesis, the test statistic is drawn from the product of two Student
#'     laws with 'sample_size-1' degrees of freedom.
#'
#' @details For details, we refer  to the technical report  "Optimal Tests of
#'     the Composite Null Hypothesis Arising in Mediation Analysis", by
#'     Miles & Chambaz (2021), https://arxiv.org/abs/2107.07575 
#' 
#' @return A list, consisting  of: \describe{ \item{t:}{a \code{vector} of two
#'     \code{numeric}s, the test statistic, or a 'n x 2' \code{matrix} of such
#'     test  statistics;} \item{alpha:}{a  \code{numeric}, the  type-I error;}
#'     \item{truncation:}{a  nonnegative  \code{numeric},  used to  bound  the
#'     rejection    region   away    from   the    null   hypothesis    space}
#'     \item{sample_size:}{an \code{integer},  the size of the  sample used to
#'     derive   the  test   statistic}  \item{decision:}{a   \code{vector}  of
#'     \code{logical}s, \code{FALSE}  if the  null hypothesis can  be rejected
#'     for  the  alternative  at  level 'alpha'  and  \code{TRUE}  otherwise;}
#'     \item{pval:}{a \code{vector}  of \code{numeric}s,  the p-values  of the
#'     tests;} \item{method:}{the \code{character} "minimax".} }
#'
#' @seealso [BH_mediation_test_minimax()], which builds upon the present function to implement a Benjamini-Hochberg procedure to control false discovery rate.
#' @examples
#' n <- 10
#' x <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = diag(c(1, 1)))
#' delta <- matrix(stats::runif(2 * n, min = -3, max = 3), ncol = 2)
#' epsilon <- stats::rbinom(n, 1, 1/2)
#' delta <- delta * cbind(epsilon, 1 - epsilon)
#' x <- x + delta
#' (mt <- mediation_test_minimax(x, alpha = 1/20))
#' plot(mt)
#'
#' @export
mediation_test_minimax <- function(t, alpha = 0.05, truncation = 0, sample_size = Inf) {
    t <- R.utils::Arguments$getNumerics(t)
    if (is.vector(t)) {
        t <- matrix(t, ncol = 2)
    }
    if (!ncol(t) == 2) {
        R.oo::throw("Argument 't' should be a vector consisting of two real numbers or a n x 2 matrix or data frame (the test statistic(s) in the real plane), not ", t) 
    }
    alpha <- R.utils::Arguments$getNumeric(alpha, c(0, 1/2))
    truncation <- R.utils::Arguments$getNumeric(truncation, c(0, Inf))
    sample_size <- R.utils::Arguments$getNumeric(sample_size, c(2, Inf))
    if (!is.infinite(sample_size)) {
        if (!is.integer(sample_size)) {
            sample_size <- as.integer(sample_size)
            warning("Coercing argument 'sample_size' to an integer.")
        }
    }
    make_decision <- function(tt, aalpha) {
        aalpha_inverse <- floor(1/aalpha)
        qtls <- stats::qt(seq(from = 1 - aalpha_inverse * aalpha/2,
                              to = 1 - aalpha_inverse * alpha/2 + aalpha_inverse * aalpha/2,
                              by = aalpha/2),
                          df = sample_size)
        tt <- abs(tt)
        int_x <- findInterval(tt[, 1], qtls)
        int_y <- findInterval(tt[, 2], qtls)
        decision <- (int_x == int_y)
        names(decision) <- rep("rejection", length(decision))
        return(decision)
    }
    compute_pval <- function(tt) {
        tt <- abs(tt)
        phi_maxtt <- stats::pt(pmax(tt[, 1], tt[, 2]), df = sample_size)
        phi_mintt <- stats::pt(pmin(tt[, 1], tt[, 2]), df = sample_size)
        K <- floor((1 - phi_maxtt)/(phi_maxtt - phi_mintt))
        sums <- c(0, cumsum(1/(1:max(K))))
        ##
        pvals <- 2*(1 - phi_mintt)/(K+1) + 2*(phi_maxtt-phi_mintt)*sums[K+1]
        return(pvals)
    }
    decision <- make_decision(t, alpha)
    if (truncation > 0) {
        decision <- decision & (apply(abs(t), 1, min) >= truncation)
    }
    pvals <- compute_pval(t)
    
    out <- list(t = t,
                alpha = alpha,
                truncation = truncation,
                sample_size = sample_size,
                decision = decision,
                pvals = pvals,
                method = "minimax")
    class(out) <- "mediation.test"
    return(out)
}

#' Carries   out   the   minimax   optimal  test   of   the   composite   null
#' "\eqn{\delta_x    \times     \delta_y=0}"    against     its    alternative
#' "\eqn{\delta_x \times \delta_y\neq  0}" based on the test  statistic in the
#' real  plane,   aiming  for  a  control   of  false  discovery  rate   à  la
#' Benjamini-Hochberg.
#' 
#' @param t  A  \code{vector}  consisting of  two  \code{numeric}s, the  test
#'   statistic in  the real  plane, or a  'n x 2'  \code{matrix} of  such test
#'   statistics.
#'
#' @param alpha A positive \code{numeric}, the wished false discovery rate.
#'
#' @param truncation A nonnegative  \code{numeric} used to bound the rejection
#'     region away  from the null  hypothesis space.  Defaults to 0,  in which
#'     case the  rejection region is minimax optimal.
#' 
#' @param sample_size An \code{integer}  (larger than  one), the size  of the
#'     sample used  to derive the  test statistic. Defaults to  'Inf', meaning
#'     that, under the  null hypothesis, the test statistic is  drawn from the
#'     \eqn{N_2(0,I_2)} law.  If  the integer is finite, then,  under the null
#'     hypothesis, the test statistic is drawn from the product of two Student
#'     laws with 'sample_size-1' degrees of freedom.
#'
#' @details For details, we refer  to the technical report  "Optimal Tests of
#'     the Composite Null Hypothesis Arising in Mediation Analysis", by
#'     Miles & Chambaz (2021), https://arxiv.org/abs/2107.07575 
#' 
#' @return A list, consisting  of: \describe{ \item{t:}{a \code{vector} of two
#'     \code{numeric}s, the test statistic, or a 'n x 2' \code{matrix} of such
#'     test statistics;}  \item{alpha:}{a \code{numeric}, the  false discovery
#'     rate;} \item{truncation:}{a  nonnegative \code{numeric}, used  to bound
#'     the   rejection   region  away   from   the   null  hypothesis   space}
#'     \item{sample_size:}{an \code{integer},  the size of the  sample used to
#'     derive   the  test   statistic}  \item{decision:}{a   \code{vector}  of
#'     \code{logical}s, \code{FALSE}  if the  null hypothesis can  be rejected
#'     for the  alternative at  false discovery  rate 'alpha'  and \code{TRUE}
#'     otherwise;}  \item{pval:}{a   \code{vector}  of   \code{numeric}s,  the
#'     p-values   of    the   tests;}    \item{method:}{the   \code{character}
#'     "minimax+BH".}}
#' 
#' @seealso [mediation_test_minimax()], upon which this function builds.
#' @examples
#' n <- 10
#' x <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = diag(c(1, 1)))
#' delta <- matrix(stats::runif(2 * n, min = -3, max = 3), ncol = 2)
#' epsilon <- stats::rbinom(n, 1, 1/2)
#' delta <- delta * cbind(epsilon, 1 - epsilon)
#' x <- x + delta
#' (mt <- BH_mediation_test_minimax(x, alpha = 0.05))
#' plot(mt)
#'
#' @export
BH_mediation_test_minimax <- function(t, alpha = 0.05, truncation = 0, sample_size = Inf) {
    mt <- mediation_test_minimax(t = t, alpha = alpha, truncation = truncation, sample_size = sample_size)
    pval <- stats::p.adjust(mt$pval)
    mt$pval <- pval
    mt$decision <- (pval <= alpha) & (apply(abs(t), 1, min) >= truncation)
    mt$method <- "minimax+BH"
    return(mt)
}
