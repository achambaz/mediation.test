#' Prints the output of \code{function}s \code{mediation_test_minimax} and \code{mediation_test_Bayes}.
#' 
#' @param x An output of \code{function} \code{mediation_test_minimax} or \code{mediation_test_Bayes}.
#'
#' @param ... Not used.
#' 
#' @return Nothing.
#'
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
#' @method print mediation.test
#' 
#' @export
print.mediation.test <- function(x, ...) {
    single <- nrow(x$t) == 1
    cat("Testing the composite null 'delta_x * delta_y = 0' against its alternative 'delta_x * delta_y != 0':\n")
    cat("* method:\n")
    print(x$method)
    cat("* test statictic:\n")
    if (!single) {
        print(utils::head(x$t))
    } else {
        print(x$t)
    } 
    cat("* wished type-I error:\n")
    if (x$method == "minimax") {
        print(x$alpha)
    } else if (x$method == "Bayes") {
        alpha <- stringr::str_extract(attr(x$map, "info"), "[0-9\\.]+$")
        alpha <- as.numeric(alpha)
        print(x$alpha)
        cat("* loss function:\n")
        loss <- stringr::str_extract(attr(x$map, "info"), "^[^;]+")
        print(loss)
    } else {
        R.oo::throw("Item 'method' should be either 'minimax' or 'Bayes', not ", x$method) 
    }
    cat("* user-supplied truncation parameter:\n")
    print(x$truncation)
    if (x$method == "minimax") {
        cat("* size  of the sample used to derive the  test statistic ('Inf' to use a Gaussian approximation; otherwise, use a product of Student laws):\n")
        print(x$sample_size)
    }
    cat("* decision:\n")
    DECISION <- c(sprintf("cannot reject the null for its alternative with confidence %1.3f\n", x$alpha),
                  sprintf("can reject the null for its alternative with confidence %1.3f\n", x$alpha))
    decision <- DECISION[x$decision + 1]
    cat(utils::head(decision))
    if (!single) {
        cat("...\n")
    }
    cat("* p-value:\n")
    print(utils::head(x$pval))
    if (!single) {
        cat("...\n")
    }
    invisible()
}

#' Plots the output of \code{function} \code{mediation_test_minimax} with the rejection region. 
#' 
#' @param x An output of \code{function} \code{mediation_test_minimax}.
#'
#' @param filename  Either \code{NULL} (defaults) or a file  name to create on
#'   disk.
#'
#' @param return_fig  A \code{logical}, to request that the  'ggplot2' object be
#'     returned (if 'TRUE') or not (if 'FALSE'). 
#' 
#' @param xlim,ylim  Two \code{vectors} of \code{numeric}s,  the wished x-axis
#'   and y-axis ranges (both default to 'c(-4, 4)').
#' 
#' @param ... Not used.
#' 
#' @return Nothing unless 'return_fig' is set to 'TRUE', in which case the function returns the 'ggplot2' object.
#'
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
#' @method plot mediation.test
#' 
#' @export
plot.mediation.test <- function(x, filename = NULL, return_fig = FALSE, xlim = c(-4, 4), ylim = c(-4, 4), ...) {
    if (!is.null(filename)) {
        file <- R.utils::Arguments$getFilename(filename)
    }
    return_fig <- R.utils::Arguments$getLogical(return_fig)
    xlim <- R.utils::Arguments$getNumerics(xlim)
    ylim <- R.utils::Arguments$getNumerics(ylim)
    if (length(xlim) != 2 | length(ylim) != 2) {
        R.oo::throw("Arguments 'xlim' and 'ylim' must be vectors of length 2.\n")
    }
    xlim <- sort(xlim)
    ylim <- sort(ylim)
    ## to please R CMD
    xmin <- NULL
    xmax <- NULL
    ymin <- NULL
    ymax <- NULL
    Xn <- NULL
    Yn <- NULL
    ##
    if (x$method == "minimax") {
        alpha <- x$alpha
        alpha_inverse <- floor(1/alpha)
        ## a <- stats::qt(seq(from = 1 - alpha_inverse * x$alpha/2,
        ##                    to = 1 - alpha_inverse * x$alpha/2 + alpha_inverse * x$alpha/2,
        ##                    by = x$alpha/2),
        ##                df = x$sample_size)
        a <- stats::qt(seq(from = 1 - alpha_inverse * alpha/2,
                           to = 1 - alpha_inverse * alpha/2 + alpha_inverse * alpha/2,
                           by = alpha/2),
                       df = x$sample_size)
        a <- pmax(a, x$truncation)
        df <- tibble::tibble(
                          xmin = c(rep(a, 2), rep(-a, 2)),
                          xmax = c(rep(c(a[-1], Inf), 2), rep(c(-a[-1], -Inf), 2)),
                          ymin = rep(c(a, -a), 2),
                          ymax = rep(c(a[-1], Inf, -a[-1], -Inf), 2)
                      )
        scatter_df <- tibble::tibble(Xn = x$t[, 1], Yn = x$t[, 2], decision = x$decision)
        ##
        msg <- sprintf("type-I error: %.3f; truncation: %.3f", alpha, x$truncation)
        ##
        fig <- ggplot2::ggplot() +
            ggplot2::geom_rect(data = df,
                               ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = x$alpha),
                               alpha = 1) +
            ggplot2::geom_point(data = scatter_df,
                                ggplot2::aes(x = Xn, y = Yn, color = !x$decision)) +
            ggplot2::theme(legend.position = "none") +
            ggplot2::labs(x = bquote(X[n]), y = bquote(Y[n])) +
            ggplot2::geom_hline(yintercept = 0) + 
            ggplot2::geom_vline(xintercept = 0) + 
            ggplot2::scale_x_continuous(limits = xlim, expand = c(0, 0)) +
            ggplot2::scale_y_continuous(limits = ylim, expand = c(0, 0)) +
            ggplot2::coord_fixed() +
            ggplot2::ggtitle("Minimax optimal test",
                             subtitle = msg)
    } else if (x$method == "Bayes") {
        scatter_df <- tibble::tibble(Xn = x$t[, 1], Yn = x$t[, 2], decision = x$decision)
        fig <- plot(x$map, return_fig = TRUE)
        fig <- fig +
            ggplot2::geom_point(data = scatter_df,
                                ggplot2::aes(x = Xn, y = Yn, color = !x$decision)) +
            ggplot2::theme(legend.position = "none") +
            ggplot2::labs(x = bquote(X[n]), y = bquote(Y[n])) +
            ggplot2::geom_hline(yintercept = 0) + 
            ggplot2::geom_vline(xintercept = 0) + 
            ggplot2::scale_x_continuous(limits = xlim, expand = c(0, 0)) +
            ggplot2::scale_y_continuous(limits = ylim, expand = c(0, 0)) +
            ggplot2::ggtitle("Bayes risk optimal test",
                             subtitle = paste0(attr(x$map, "info"), "; ",
                                              sprintf("truncation: %.3f", x$truncation)))
    } else {## should never happen
        R.oo::throw("Argument 'method' should be either 'minimax' or 'Bayes', not ", x$method)         
    }
    if (!is.null(filename)) {
        ggplot2::ggsave(fig, file = file)
    } else {
        outside_of_range <- any(x$t[, 1] < xlim[1] | x$t[, 1] > xlim[2]) |
            any(x$t[, 2] < ylim[1] | x$t[, 2] > ylim[2]) 
        if (outside_of_range) {
            message("Some points fall outside the range of the figure.\n")
        }
        suppressWarnings(print(fig))
    }
    if (!return_fig) {
        invisible()
    } else {
        return(fig)
    }
}


## ####################################################### #
## TEMPORARILY REMOVED DUE TO ISSUES WITH LIBRARY "plotly" #
## ####################################################### #

## #' Plots  the   output  of  \code{function}  \code{mediation_test_minimax}   with  the
## #' corresponding 3D power surface.
## #' 
## #' @param x An output of \code{function} \code{mediation_test_minimax}.
## #'
## #' @param filename  Either \code{NULL} (defaults) or a file  name to create on
## #'   disk. If \code{NULL}, then the plot is sent to the default web browser.
## #'
## #' @param nbins_xy  An \code{integer} (chosen between 50  and 1000, defaulting
## #'   to 501), the number of bins used to discretize the along each axis.
## #'
## #' @param range_delta_x,range_delta_y Two \code{vector}s  of \code{numeric}s
## #'   describing the  ranges of \eqn{\delta_x} and  \eqn{\delta_y}. By default,
## #'   they are both set to \code{c(-6, 6)}.
## #' 
## #' @param ... Not used.
## #'
## #' @details For simplicity, the wished type-I error 'alpha' is approximated by the smallest unit fraction larger than 'alpha'.
## #' 
## #' @return Nothing.
## #'
## #' @examples
## #' n <- 10
## #' x <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = diag(c(1, 1)))
## #' delta <- matrix(stats::runif(2 * n, min = -3, max = 3), ncol = 2)
## #' epsilon <- stats::rbinom(n, 1, 1/2)
## #' delta <- delta * cbind(epsilon, 1 - epsilon)
## #' x <- x + delta
## #' (mt <- mediation_test(x, alpha = 1/20))
## #' show.3Dplot <- FALSE # change to 'TRUE' to make and show the 3D plot
## #' if (show.3Dplot) {
## #'   plot3d(mt)
## #' }
## #'
## #' 
## #' # @aliases plot3d
## #' 
## #' # @export plot3d
## #'
## #' # @export
## R.methodsS3::setMethodS3(
##                  "plot3d", "mediation.test",
##                  function(x, filename = NULL, nbins_xy = 501, range_delta_x = c(-6, 6), range_delta_y = c(-6, 6), ...) {
##                      if (!is.null(filename)) {
##                          file <- R.utils::Arguments$getFilename(filename)
##                      }
##                      nbins_xy <- R.utils::Arguments$getInteger(nbins_xy, c(50, 1e3))
##                      range_delta_x <- R.utils::Arguments$getNumerics(range_delta_x)
##                      range_delta_y <- R.utils::Arguments$getNumerics(range_delta_y)
##                      if (length(range_delta_x) != 2 | length(range_delta_y) != 2) {
##                          R.oo::throw("Arguments 'range_delta_x' and 'range_delta_y' should be two vectors consisting of two real numbers, the ranges of 'delta_x' and 'delta_y', not ",
##                                      delta_x, delta_y)
##                      }
##                      ## to please R CMD
##                      xmin <- NULL
##                      xmax <- NULL
##                      ymin <- NULL
##                      ymax <- NULL
##                      ##
##                      compute_power_surface <- Vectorize(
##                          function(d_x, d_y, tib) {
##                              sum(abs((stats::pnorm(tib$xmax - d_x) -
##                                       stats::pnorm(tib$xmin - d_x)) * 
##                                      (stats::pnorm(tib$ymax - d_y) -
##                                       stats::pnorm(tib$ymin - d_y))))
##                          },
##                          c("d_x", "d_y"))
##                      ##
##                      alpha <- 1/floor(1/x$alpha)
##                      a <- stats::qnorm(seq(.5, 1-alpha/2, alpha/2))
##                      df <- tibble::tibble(
##                                        xmin = c(rep(a, 2), rep(-a, 2)),
##                                        xmax = c(rep(c(a[-1], Inf), 2), rep(c(-a[-1], -Inf), 2)),
##                                        ymin = rep(c(a, -a), 2),
##                                        ymax = rep(c(a[-1], Inf, -a[-1], -Inf), 2)
##                                    )
##                      delta_x <- seq(min(range_delta_x), max(range_delta_x), length = nbins_xy)
##                      delta_y <- seq(min(range_delta_y), max(range_delta_y), length = nbins_xy)
##                      power_surface <- outer(delta_x, delta_y, compute_power_surface, tib = df)
##                      ##
##                      fig <- plotly::plot_ly(x = ~delta_x, y = ~delta_y, z = ~power_surface)
##                      fig <- plotly::add_surface(fig,
##                                                 contours = list(
##                                                     z = list(
##                                                         show = TRUE,
##                                                         usecolormap = TRUE,
##                                                         highlightcolor = "#ff0000",
##                                                         project = list(z = TRUE)
##                                                     )
##                                                 ))
##                      fig <- plotly::add_trace(fig,
##                                               x = x$t[, 1], y = x$t[, 2],
##                                               z = compute_power_surface(x$t[, 1], x$t[, 2], tib = df),
##                                               mode = "markers", type = "scatter3d", 
##                                               marker = list(
##                                                   size = 5, color = "red", symbol = 104))
##                      fig <- plotly::layout(fig,
##                                            scene = list(
##                                                camera = list(
##                                                    eye = list(x = 0, y = -1.75, z = 1.75)
##                                                ),
##                                                xaxis = list(title = "delta*_x"),
##                                                yaxis = list(title = "delta*_y"),
##                                                zaxis = list(title = "Rej. prob.")
##                                            )
##                                            )
##                      if (!is.null(filename)) {
##                          save_to_disk <- try(plotly::orca(fig, file = file))
##                          if (inherits(save_to_disk, "try-error")) {
##                              warning("Saving image to disk failed.")
##                          }
##                      } else {
##                          print(fig)
##                      }
##                      invisible()
##                  }
##              )

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
#' @param  sample_size An \code{integer}  (larger than  one), the size  of the
#'     sample used  to derive the  test statistic. Defaults to  'Inf', meaning
#'     that, under the  null hypothesis, the test statistic is  drawn from the
#'     \eqn{N_2(0,I_2)} law.  If  the integer is finite, then,  under the null
#'     hypothesis, the test statistic is drawn from the product of two Student
#'     laws with 'sample_size-1' degrees of freedom.
#'
#' @details  For details, we refer  to the technical report  "Optimal Tests of
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

#' Plots the "map" of rejection probabilities  for the Bayes risk optimal test
#' of  the  composite  null  "\eqn{\delta_x \times  \delta_y=0}"  against  its
#' alternative  "\eqn{\delta_x  \times  \delta_y\neq  0}" based  on  the  test
#' statistic in the real  plane.  
#' 
#' @param x An output of \code{function} \code{}.
#'
#' @param filename  Either \code{NULL} (defaults) or a file  name to create on
#'   disk.
#'
#' @param return_fig  A \code{logical}, to request that the  'ggplot2' object be
#'     returned (if 'TRUE') or not (if 'FALSE'). 
#' 
#' @param ... Not used.
#' 
#' @return Nothing unless 'return_fig' is set to 'TRUE', in which case the function returns the 'ggplot2' object.
#'
#' @examples
#' 
#' map <- compute_map_rejection_probs(alpha = 0.05, K = 16, loss = "0-1")
#' plot(map)
#'
#' @method plot map.rejection.probabilities
#' 
#' @export
plot.map.rejection.probabilities <- function(x, filename = NULL, return_fig = FALSE, ...) {
    if (!is.null(filename)) {
        file <- R.utils::Arguments$getFilename(filename)
    }
    return_fig <- R.utils::Arguments$getLogical(return_fig)
    ##
    K <- nrow(x)/2
    msg <- attr(x, "info")
    alpha <- stringr::str_extract(msg, "[0-9\\.]+$")
    alpha <- as.numeric(alpha)
    B <- 2 * stats::qnorm(1 - alpha / 2)
    r_xy <- seq(-B, B, length.out = 2 * K)
    sbttl <- sprintf("(restricted to [%.3f,%.3f])", -B, B)
    ##
    ## to please R CMD
    delta_x <- NULL
    delta_y <- NULL
    prob <- NULL
    
    df <- tidyr::expand_grid(
                     delta_x = r_xy,
                     delta_y = r_xy
                 ) 
    df <- dplyr::mutate(df, 
                        probs = as.vector(x)
                        )
    ##
    fig <- ggplot2::ggplot(df) +
        ggplot2::geom_raster(ggplot2::aes(x = delta_x, y = delta_y, fill = probs)) +
        ggplot2::labs(x = bquote(delta[x]), y = bquote(delta[y])) +
        ggplot2::geom_hline(yintercept = 0) + 
        ggplot2::geom_vline(xintercept = 0) + 
        ggplot2::coord_fixed() +
        ggplot2::ggtitle("Map of rejection probabilities",
                         subtitle = paste(msg, sbttl)) +
        ggplot2::guides(fill = ggplot2::guide_legend(title = "rejection\nprob."))
    if (!is.null(filename)) {
        ggplot2::ggsave(fig, file = file)
    } else {
        suppressWarnings(print(fig))
    }
    if (!return_fig) {
        invisible()
    } else {
        return(fig)
    }
}


#' Computes the  "map" of rejection  probabilities for the Bayes  risk optimal
#' test of the  composite null "\eqn{\delta_x \times  \delta_y=0}" against its
#' alternative  "\eqn{\delta_x  \times  \delta_y\neq  0}" based  on  the  test
#' statistic in the real  plane.  The Bayes risk is induced  by either the 0-1
#' or the bounded quadratic loss function.
#'
#' @param alpha A positive  \code{numeric}, the wished type-I  error (default
#'     value 0.05).
#'
#' @param K An  \code{integer} parametrizing the  discretization of  the null
#'     hypothesis  test and  of the  plane.  The  complexity of  the resulting
#'     optimization  problem is  \eqn{O(K^2)}, so  'K' should  neither be  too
#'     small nor too large. The default  value is 'K=16'.  We recommend 'K=64'
#'     as a sensible value and enforce \eqn{2 \leq K \leq 128}.
#'
#' @param loss  A  \code{character},   either  "0-1"  (default   value)  or
#'     "quadratic", to indicate which loss function to consider.
#'
#' @param return_solver A \code{logical}, to request that the 'lpSolve' object
#'     be returned  (if 'TRUE') or  not (if  'FALSE', default). Note  that the
#'     'lpSolve' object can be quite big.
#' 
#' @return A  '2*K  x 2*K'  \code{matrix} whose  '(i,j)'  coefficient is  the
#'     probability to  reject the null  when the  test statistic falls  in the
#'     square \eqn{[b(i-1-K)/K,  b(i-K)/K] \times [b(j-1-K)/K,  b(j-K)/K]}. If
#'     'return_solver'  is 'TRUE',  then the  matrix has  an attribute  called
#'     'solver'  which is  the complete  output of  the optimization  function
#'     'lpSolve::lp' used to determine the  probabilities.
#'
#' @examples
#' map <- compute_map_rejection_probs(alpha = 0.05, K = 16, loss = "0-1")
#' plot(map)
#' 
#' @export
compute_map_rejection_probs <- function(alpha = 0.05, K = 16, loss = c("0-1", "quadratic"),
                                        return_solver = FALSE) {
    alpha <- R.utils::Arguments$getNumeric(alpha, c(0, 1/2))
    K <- R.utils::Arguments$getInteger(K, c(2, 128))
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
    ##
    weights <- array(0, c(2 * K, 2 * K, 8 * K + 1))
    for (i in 1:(2 * K)) {
        for (j in 1:(2 * K)) {
            weights[i,j,] <- (stats::pnorm(r_xy[i + 1] - G_x) - stats::pnorm(r_xy[i] - G_x)) *
                (stats::pnorm(r_xy[j + 1] - G_y) - stats::pnorm(r_xy[j] - G_y))
        }
    }
    ##
    if (loss == "0-1") {
        SD <- sqrt(1 + sigma^2)
        risk <- stats::pnorm(r_xy[1:(2 * K)], sd = SD) - stats::pnorm(r_xy[2:(2 * K + 1)], sd = SD)
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

    full_map <- lpSolve::lp(direction = "max",
                            objective.in = as.vector(risk),
                            const.mat = cbind(apply(weights, 3, as.vector),
                                              diag(4 * K^2)), 
                            const.dir = rep("<=", 4 * K^2 + 8 * K + 1),
                            const.rhs = c(rhs, rep(1, 4 * K^2)),
                            transpose.constraints = FALSE)
    
    map <- matrix(pmin(1, pmax(0, full_map$solution)),
                  nrow = 2 * K, ncol = 2 * K)
    attr(map, "info") <- paste(loss, "loss;", "type-I error:", alpha)
    class(map) <- "map.rejection.probabilities"

    if (!return_solver) {
        attr(map, "solver") <- NA
    } else {
        attr(map, "solver") <- full_map
    }
    
    return(map)
}

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
#' @details For details, we refer  to the technical report  "Optimal Tests of
#'     the Composite Null Hypothesis Arising in Mediation Analysis", by
#'     Miles & Chambaz (2021), https://arxiv.org/abs/2107.07575 
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
#' x <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = diag(c(1, 1)))
#' delta <- matrix(stats::runif(2 * n, min = -3, max = 3), ncol = 2)
#' epsilon <- stats::rbinom(n, 1, 1/2)
#' delta <- delta * cbind(epsilon, 1 - epsilon)
#' x <- x + delta
#' (mt <- mediation_test_Bayes(x, map = map_01_0.05))
#' plot(mt)
#'
#' @export
mediation_test_Bayes <- function(t, map = map_01_0.05, truncation = 0) {
    t <- R.utils::Arguments$getNumerics(t)
    truncation <- R.utils::Arguments$getNumeric(truncation, c(0, Inf))
    if (is.vector(t)) {
        t <- matrix(t, ncol = 2)
    }
    if (!ncol(t) == 2) {
        R.oo::throw("Argument 't' should be a vector consisting of two real numbers or a n x 2 matrix or data frame (the test statistic(s) in the real plane), not ", t) 
    }
    msg <- attr(map, "info")
    alpha <- stringr::str_extract(msg, "[0-9\\.]+$")
    alpha <- as.numeric(alpha)
    ##
    make_decision <- function(tt, aalpha, mmap) {
        K <- ncol(mmap)/2
        B <- 2 * qnorm(1 - aalpha/2)
        abstt <- abs(tt)
        max_abstt <- pmax(abstt[, 1], abstt[, 2])
        min_abstt <- pmin(abstt[, 1], abstt[, 2])
        ##
        decision <- (max_abstt > B) & (min_abstt > B/2)
        idx <- which(max_abstt <= B)
        if (length(idx) >= 1) {
            ij <- cbind(ceiling((tt[idx, 1] + B) * K/B),
                        ceiling((tt[idx, 2] + B) * K/B))
            decision[idx] <- as.logical(rbinom(length(idx), 1, map[ij]))
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

