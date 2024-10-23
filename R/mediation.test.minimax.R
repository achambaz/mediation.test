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
#' @param compute_pvals A   \code{logical}  indicating   whether  or   not
#'     (conservative) p-values should be computed. 
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
#'     \item{sample_size:}{an \code{integer},  the size of the  sample used to
#'     derive   the  test   statistic}  \item{decision:}{a   \code{vector}  of
#'     \code{logical}s, \code{FALSE}  if the  null hypothesis can  be rejected
#'     for  the  alternative  at  level 'alpha'  and  \code{TRUE}  otherwise;}
#'     \item{pval:}{a \code{vector}  of \code{numeric}s,  the (conservative)
#'    p-values  of the tests;} \item{method:}{the \code{character} "minimax".}}
#'
#' @seealso [mediation_test_minimax_BH()], which builds upon the present function to implement a Benjamini-Hochberg procedure to control false discovery rate.
#' @examples
#' n <- 10
#' x <- MASS::mvrnorm(2 * n, mu = c(0, 0), Sigma = diag(c(1, 1)))
#' delta <- matrix(stats::runif(4 * n, min = -3, max = 3), ncol = 2)
#' epsilon <- stats::rbinom(n, size = 1, prob = 1/2)
#' delta <- delta * cbind(c(epsilon, rep(1, n)),
#'                        c(1 - epsilon, rep(1, n)))
#' x <- x + delta
#' (mt <- mediation_test_minimax(x, alpha = 1/20))
#' plot(mt)
#'
#' @export
mediation_test_minimax <- function(t, alpha = 0.05, truncation = 0, sample_size = Inf, compute_pvals = TRUE) {
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
    compute_pvals <- R.utils::Arguments$getLogical(compute_pvals)
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
    if (compute_pvals) {
        pvals <- compute_pval(t)
    } else {
        pvals <- NA
    }
    
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
#' @param BH A \code{character},  either "statistics" or  "pval", determining
#'     whether the Benjamini-Hochberg  testing procedure is based  on the test
#'     statistics or on the (conservative) p-values.
#' 
#' @param sample_size A \code{integer}  (larger than  one), the size  of the
#'     sample used  to derive the  test statistic. Defaults to  'Inf', meaning
#'     that, under the  null hypothesis, the test statistic is  drawn from the
#'     \eqn{N_2(0,I_2)} law.  If  the integer is finite, then,  under the null
#'     hypothesis, the test statistic is drawn from the product of two Student
#'     laws with 'sample_size-1' degrees of freedom.
#'
#' @param K An \code{integer} (10 by  default) used internally to speed up the
#'     Benjamini-Hochberg  testing procedure  when  it is  based  on the  test
#'     statistics (more details in the  Details section). Must be smaller than
#'     'nrow(t)'.
#' 
#'  @details  Suppose  we  are  carrying out  \eqn{J}  tests.   Starting  from
#'     \eqn{j=1},  we   increment  \eqn{j}  until  fewer   than  \eqn{j}  null
#'     hypotheses  are   rejected  at  level  \eqn{\alpha\times   j/J}  for  K
#'     consecutive iterations.   Eventually, we reject all  null hypotheses at
#'     level \eqn{\alpha \times  jj/J}, where \eqn{jj} is  the largest integer
#'     so far for which at least \eqn{jj} null hypotheses are rejected at this
#'     level.  For  further   details,  we  refer  to   the  technical  report
#'     "Optimal Tests of the Composite Null Hypothesis Arising in Mediation Analysis",
#'     by Miles & Chambaz (2024), https://arxiv.org/abs/2107.07575
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
#'     original  (conservative)  p-values  of the  tests;}  \item{method:}{the
#'     \code{character}  "minimax_BH",} \item{BH:}{a  \code{character}, either
#'     'statistics'  or  'pval',  determining whether  the  Benjamini-Hochberg
#'     testing  procedure  is   based  on  the  test  statistics   or  on  the
#'     (conservative) p-values;}\item{K:}{the \code{integer} K used internally
#'     to speed up the Benjamini-Hochberg procedure based on the test statistics.}}
#' 
#' @seealso [mediation_test_minimax()], upon which this function builds.
#' @examples
#' n <- 100
#' x <- MASS::mvrnorm(2 * n, mu = c(0, 0), Sigma = diag(c(1, 1)))
#' delta <- matrix(stats::runif(4 * n, min = -5, max = 5), ncol = 2)
#' epsilon <- stats::rbinom(n, size = 1, prob = 1/2)
#' delta <- delta * cbind(c(epsilon, rep(1, n)),
#'                        c(1 - epsilon, rep(1, n)))
#' x <- x + delta
#' mt_BH_stats <- mediation_test_minimax_BH(x, alpha = 1/20, BH = "statistics")
#' (sum(mt_BH_stats$decision))
#' mt_BH_pval <- mediation_test_minimax_BH(x, alpha = 1/20, BH = "pval")
#' (sum(mt_BH_pval$decision))
#'
#' @export
mediation_test_minimax_BH <- function(t, alpha = 0.05, truncation = 0, BH = c("statistics", "pval"),
                                      sample_size = Inf, K = 10L) {
    t <- R.utils::Arguments$getNumerics(t)
    alpha <- R.utils::Arguments$getNumeric(alpha, c(0, 1/2))
    K <- R.utils::Arguments$getInteger(K, c(1, Inf))
    if (is.vector(t)) {
        t <- matrix(t, ncol = 2)
    }
    nn <- nrow(t)
    if (nn == 1) {
        warning("One single test statistic. Switching to 'mediation_test_minimax'.")
        mt <- mediation_test_minimax(t, alpha, truncation, sample_size)
    } else {
        if (K >= nn) {
            K <- nn - 1
            warning("Replacing the user-supplied 'K' with 'nrow(tt)-1'.")
        }
        BH <- match.arg(BH)
        if (BH == "pval") {## regular BH testing procedure based on the (conservative) p-values
            mt <- mediation_test_minimax(t, alpha, truncation, sample_size)
            pval <- stats::p.adjust(mt$pval)
            decision <- (pval <= alpha) & (apply(abs(t), 1, min) >= truncation)
            names(decision) <- rep("rejection", length(decision))
            mt$decision <- decision
        } else {## BH testing procedure based on the test statistics
            fun <- function(aalpha) {
                decision <- mediation_test_minimax(t, aalpha, truncation, sample_size, FALSE)$decision
                sum(decision)
            }
            carry_on <- TRUE
            from <- -K + 1
            to <- 0
            adjusted_alpha <- NA
            nb_rejections <- c()
            valid_levels <- c()
            while (carry_on & to < nn) {
                from <- from + K
                to <- min(c(to + K, nn))
                alphas <- alpha * seq(from = from/nn, to = to/nn, by = 1/nn)
                nb_rejec <- vapply(alphas, fun, integer(1))
                nb_rejections <- c(nb_rejections, nb_rejec)
                valid_lev <- (nb_rejec >= seq.int(from = from, to = to, by = 1))
                valid_levels <- c(valid_levels, valid_lev)
                cum_valid_levels <- cumsum(valid_levels)
                test <- any(cum_valid_levels -
                            dplyr::lag(cum_valid_levels, K) <= 0, na.rm = TRUE)
                if (test) {## otherwise, we carry on
                    carry_on <- FALSE
                    if (max(cum_valid_levels, na.rm = TRUE) == 0) {
                        adjusted_alpha <- alpha
                        mt <- mediation_test_minimax(t, adjusted_alpha, truncation, sample_size)
                        decision <- rep(FALSE, nn)
                        names(decision) <- rep("rejection", length(decision))
                        mt$decision <- decision
                    } else {
                        idx <- min(which.max(cum_valid_levels))
                        adjusted_alpha <- alpha *  idx / nn
                        mt <- mediation_test_minimax(t, adjusted_alpha, truncation, sample_size)
                        mt$alpha <- alpha
                    }
                }
            }
            if (is.na(adjusted_alpha)) {
                mt <- mediation_test_minimax(t, alpha, truncation, sample_size)
            }
        }
        mt$method <- "minimax_BH"
        mt$BH <- BH
        mt$K <- K
    }
    return(mt)
}
