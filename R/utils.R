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
    method <- stringr::str_split(x$method, "\\+")[[1]]
    if (length(method) == 1) {## "minimax" or "Bayes"
        cat(method, "\n")
    } else if (length(method) == 2){## "minimax" and "BH"
        cat("minimax & Benjamini-Hochberg procedure\n")
    } else {
        R.oo::throw("Item 'method' should be either 'minimax' or 'minimax+BH' or 'Bayes', not ", x$method) 
    }
    cat("* test statictic:\n")
    if (!single) {
        print(utils::head(x$t))
    } else {
        print(x$t)
    }
    if (length(method) == 1) {
        cat("* wished type-I error:\n")
        if (method == "minimax") {
            print(x$alpha)
        } else if (method == "Bayes") {
            alpha <- stringr::str_extract(attr(x$map, "info"), "[0-9\\.]+$")
            alpha <- as.numeric(alpha)
            cat("* wished type-I error:\n")
            print(x$alpha)
            cat("* loss function:\n")
            loss <- stringr::str_extract(attr(x$map, "info"), "^[^;]+")
            print(loss)
        }
    } else if (length(method) == 2) {
        cat("* wished false discovery rate:\n")
        print(x$alpha)
    } else {
        R.oo::throw("Item 'method' should be either 'minimax' or 'minimax+BH' or 'Bayes', not ", x$method) 
    }
    cat("* user-supplied truncation parameter:\n")
    print(x$truncation)
    if (method[1] == "minimax") {
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

#' Plots the output of \code{function}s \code{mediation_test_minimax} and \code{mediation_test_Bayes} with the rejection region. 
#' 
#' @param x An output of \code{function}s \code{mediation_test_minimax} or \code{mediation_test_Bayes}.
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
    method <- stringr::str_split(x$method, "\\+")[[1]]
    if (method[1] == "minimax") {
        alpha <- x$alpha
        alpha_inverse <- floor(1/alpha)
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
        if (length(method) == 1) {
            msg <- sprintf("type-I error: %.3f; truncation: %.3f", alpha, x$truncation)
        } else if (length(method) == 2) {
            msg <- sprintf("false discovery rate: %.3f (standard BH procedure); truncation: %.3f", alpha, x$truncation)
        } else {
            R.oo::throw("Item 'method' should be either 'minimax' or 'minimax+BH' or 'Bayes', not ", x$method) 
        }
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
    } else if (method[1] == "Bayes") {
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
        R.oo::throw("Argument 'method' should be either 'minimax' or 'minimax+BH' or 'Bayes', not ", x$method)         
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
    probs <- NULL
    
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

