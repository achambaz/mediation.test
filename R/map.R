#' Computing 'map's of rejection probabilities for the Bayes risk optimal test
#' 
#' Computes the  "map" of rejection  probabilities for the Bayes  risk optimal
#' test of the  composite null "\eqn{\delta_x \times  \delta_y=0}" against its
#' alternative  "\eqn{\delta_x  \times  \delta_y\neq  0}" based  on  the  test
#' statistic in the real  plane.  The Bayes risk is induced  by either the 0-1
#' or  the  bounded  quadratic  loss  function, a  wished  type-I  error  and,
#' possibly, a truncation parameter.
#'
#' @param alpha A positive  \code{numeric}, the wished type-I  error (default
#'     value 0.05).
#'
#' @param K An  \code{integer} parametrizing the  discretization of  the null
#'     hypothesis  test and  of the  plane.  The  complexity of  the resulting
#'     optimization  problem is  \eqn{O(K^2)}, so  'K' should  neither be  too
#'     small nor too large. The default  value is 'K=16'.  We recommend 'K=64'
#'     as a sensible value and enforce \eqn{2 \le K \le 128}.
#'
#' @param loss  A  \code{character},   either  "0-1"  (default   value)  or
#'     "quadratic", to indicate which loss function to consider.
#'
#' @param truncation A nonnegative  \code{numeric}, enabling to bound away the
#'     rejection region from  the null hypothesis space. Its  default value is
#'     0, meaning that no truncation is required.
#' 
#' @param return_solver A \code{logical}, to request that the 'lpSolve' object
#'     be returned  (if 'TRUE') or  not (if  'FALSE', default). Note  that the
#'     'lpSolve' object can be quite big.
#' 
#'  @return A  '2*K  x 2*K'  \code{matrix} whose  '(i,j)'  coefficient is  the
#'    probability to  reject the null  when the  test statistic falls  in the
#'    square \eqn{[b(i-1-K)/K, b(i-K)/K] \times [b(j-1-K)/K, b(j-K)/K]}, with
#'    'b' set  to '2 *  stats::qnorm(1 - alpha  / 2)'. If  'return_solver' is
#'    'TRUE', then the  matrix has an attribute called 'solver'  which is the
#'    complete  output of  the  optimization function  'lpSolve::lp' used  to
#'    determine the probabilities.
#'
#' @examples
#' ## one of the four outputs of 'compute_map_rejection_probs' stored in the package
#' head(map_01_0.05, c(5, 5)) 
#' map <- compute_map_rejection_probs(alpha = 0.05, K = 16, loss = "0-1")
#' plot(map)
#' 
#' @export
compute_map_rejection_probs <- function(alpha = 0.05, K = 16, loss = c("0-1", "quadratic"),
                                        truncation = 0,
                                        return_solver = FALSE) {
    alpha <- R.utils::Arguments$getNumeric(alpha, c(0, 1/2))
    K <- R.utils::Arguments$getInteger(K, c(2, 128))
    truncation <- R.utils::Arguments$getNumeric(truncation, c(0, Inf))
    loss <- match.arg(loss)
    return_solver <- R.utils::Arguments$getLogical(return_solver)
    ##
    B <- 2 * stats::qnorm(1 - alpha / 2)
    sigma <- sqrt(B) # prior standard deviation
    ##
    r_xy <- seq(-B, B, length.out = 2 * K + 1)
    g_xy <- seq(-2 * B, 2 * B, length.out = 4 * K + 1)
    G_x <- c(g_xy, rep(0, 4 * K))
    G_y <- c(rep(0, 4 * K + 1), g_xy[-(2 * K + 1)])
    ##
    fun_rhs <- function(delta, z) {
        z <- R.utils::Arguments$getNumeric(z, c(0, Inf))
        1 - stats::pnorm(z - delta) + stats::pnorm(-z - delta)
    }
    rhs <- alpha - (fun_rhs(G_x, B/2) * fun_rhs(G_y, B/2) -
                    (fun_rhs(G_x, B/2) - fun_rhs(G_x, B)) * (fun_rhs(G_y, B/2) - fun_rhs(G_y, B)))
    if (truncation > 0) {
        not_truncated <- (r_xy[1:(2 * K)] > truncation) | (r_xy[2:(2 * K + 1)] < -truncation)
        not_truncated <- matrix(not_truncated, ncol = 1) %*% not_truncated
        not_truncated <- as.vector(not_truncated)
        const.rhs <- c(rhs, not_truncated)
    } else {
        const.rhs <- c(rhs, rep(1, 4 * K^2))
    }
    const.dir <- rep("<=", 4 * K^2 + 8 * K + 1)
    ##
    weights <- array(0, c(2 * K, 2 * K, 8 * K + 1))
    for (i in 1:(2 * K)) {
        for (j in 1:(2 * K)) {
            weights[i,j,] <- (stats::pnorm(r_xy[i + 1] - G_x) - stats::pnorm(r_xy[i] - G_x)) *
                (stats::pnorm(r_xy[j + 1] - G_y) - stats::pnorm(r_xy[j] - G_y))
        }
    }
    const.mat <- cbind(apply(weights, 3, as.vector), ## type-I error
                       diag(4 * K^2)) ## probabilities ('<=1' and '>=0' necessarily)
    ##
    if (loss == "0-1") {
        SD <- sqrt(1 + sigma^2)
        risk <- stats::pnorm(r_xy[2:(2 * K + 1)], sd = SD) - stats::pnorm(r_xy[1:(2 * K)], sd = SD)
        risk <- matrix(risk, ncol = 1) %*% risk
    } else if (loss == "quadratic") {
        SD <- sqrt(1 + sigma^2)
        risk <- sigma^2 * (stats::pnorm(r_xy[2:(2 * K + 1)]/SD) -
                           stats::pnorm(r_xy[1:(2 * K)]/SD)) -
            sigma^4/(1 + sigma ^ 2)^(3/2) * (r_xy[2:(2 * K + 1)] * stats::dnorm(r_xy[2:(2 * K + 1)]/SD) -
                                             r_xy[1:(2 * K)] * stats::dnorm(r_xy[1:(2 * K)]/SD))
        risk <- matrix(risk, ncol = 1) %*% risk
    } else {## should never happen
        R.oo::throw("Argument 'loss' should be either '0-1' or 'quadratic', not ", loss) 
    }

    if (!return_solver) {
        map <- lpSolve::lp(direction = "max",
                           objective.in = as.vector(risk),
                           const.mat = const.mat, 
                           const.dir = const.dir,
                           const.rhs = const.rhs,
                           transpose.constraints = FALSE)$solution
        map <- matrix(pmin(1, pmax(0, map)), nrow = 2 * K, ncol = 2 * K)
        attr(map, "solver") <- NA
    } else {
        full_map <- lpSolve::lp(direction = "max",
                                objective.in = as.vector(risk),
                                const.mat = const.mat, 
                                const.dir = const.dir,
                                const.rhs = const.rhs,
                                transpose.constraints = FALSE)
        map <- matrix(pmin(1, pmax(0, full_map$solution)),
                      nrow = 2 * K, ncol = 2 * K)
        attr(map, "solver") <- full_map
    }
    attr(map, "loss") <- loss
    attr(map, "alpha") <- alpha
    attr(map, "truncation") <- truncation
    class(map) <- "map.rejection.probabilities"
    
    return(map)
}
