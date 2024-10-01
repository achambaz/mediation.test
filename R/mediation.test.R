#' Prints the output of \code{function} \code{mediation_test}.
#' 
#' @param x An output of \code{function} \code{mediation_test}.
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
#' (mt <- mediation_test(x, alpha = 1/20))
#' plot(mt)
#' 
#' @method print mediation.test
#' 
#' @export
print.mediation.test <- function(x, ...) {
    single <- nrow(x$t) == 1
    cat("Testing the composite null 'delta_x * delta_y = 0' against its alternative 'delta_x * delta_y != 0':\n")
    cat("* test statictic:\n")
    print(utils::head(x$t))
    if (!single) {
        print(utils::head(x$t))
    } 
    cat("* wished type-I error:\n")
    print(x$alpha)
    cat("* user-supplied truncation parameter:\n")
    print(x$truncation)
    cat("* size  of the sample used to derive the  test statistic ('Inf' to use a Gaussian approximation; otherwise, use a product of Student laws):\n")
    print(x$sample_size)    
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

#' Plots the output of \code{function} \code{mediation_test} with the rejection region. 
#' 
#' @param x An output of \code{function} \code{mediation_test}.
#'
#' @param filename  Either \code{NULL} (defaults) or a file  name to create on
#'   disk.
#'
#' @param return_fig  A \code{logical}, to request that the  ggplot2 object be
#'     returned (if 'TRUE') or not (if 'FALSE'). 
#' 
#' @param xlim,ylim  Two \code{vectors} of \code{numeric}s,  the wished x-axis
#'   and y-axis ranges (both default to 'c(-4, 4)').
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
#' (mt <- mediation_test(x, alpha = 1/20))
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
    alpha_inverse <- floor(1/x$alpha)
    a <- stats::qt(seq(from = 1 - alpha_inverse * x$alpha/2,
                       to = 1 - alpha_inverse * x$alpha/2 + alpha_inverse * x$alpha/2,
                       by = x$alpha/2),
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
        ggplot2::coord_fixed()
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

## #' Plots  the   output  of  \code{function}  \code{mediation_test}   with  the
## #' corresponding 3D power surface.
## #' 
## #' @param x An output of \code{function} \code{mediation_test}.
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

#' Carries      out     the      test      of      the     composite      null
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
#'     tests.}  }
#'
#' @examples
#' n <- 10
#' x <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = diag(c(1, 1)))
#' delta <- matrix(stats::runif(2 * n, min = -3, max = 3), ncol = 2)
#' epsilon <- stats::rbinom(n, 1, 1/2)
#' delta <- delta * cbind(epsilon, 1 - epsilon)
#' x <- x + delta
#' (mt <- mediation_test(x, alpha = 1/20))
#' plot(mt)
#'
#' @export
mediation_test <- function(t, alpha = 0.05, truncation = 0, sample_size = Inf) {
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
                pvals = pvals)
    class(out) <- "mediation.test"
    return(out)
}

