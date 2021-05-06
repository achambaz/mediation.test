#' mediation.test:  More  Powerful  Tests  of the  Composite  Null  Hypothesis
#' Arising in Mediation Analysis
#'
#' The mediation.test package implements more  powerful tests of the composite
#' null hypothesis  arising in mediation  analysis. It contains  the following
#' functions: mediation_test, print.mediation.test, plot.mediation.test
#' @docType package
#' @name mediation.test
NULL
#> NULL

#' Prints the output of \code{function} \code{mediation_test}.
#' 
#' @param x An output of \code{function} \code{mediation_test}.
#'
#' @param ... Not used.
#' 
#' @return Nothing.
#'
#' @examples
#' (mt <- mediation_test(c(1.1, 2.2), alpha = 1/20))
#' plot(mt)
#' 
#' @method print mediation.test
#' 
#' @export
print.mediation.test <- function(x, ...) {
  cat("Testing the composite null 'x * y = 0' against its alternative 'x * y != 0':\n")
  cat("* test statictic:\n")
  print(x$t)
  cat("* wished type-I error:\n")
  print(x$alpha)
  cat("* decision:\n")
  decision <- ifelse(x$decision,
                     sprintf("can reject the null for its alternative with confidence %1.3f\n", x$alpha),
                     sprintf("cannot reject the null for its alternative with confidence %1.3f\n", x$alpha))
  cat(decision)
  cat("* (random) approximate p-value:\n")
  print(x$pval)
  int_pval <- sprintf("  falling deterministically in interval [%1.3f, %1.3f]\n", x$int_pval[1], x$int_pval[2])
  cat(int_pval)
  invisible()
}

#' Plots the output of \code{function} \code{mediation_test} with the rejection region. 
#' 
#' @param x An output of \code{function} \code{mediation_test}.
#'
#' @param filename  Either \code{NULL} (defaults) or a file  name to create on
#'   disk.
#' 
#' @param ... Not used.
#' 
#' @return Nothing.
#'
#' @examples 
#' (mt <- mediation_test(c(1.1, 2.2), alpha = 1/20))
#' plot(mt)

#' @method plot mediation.test
#' 
#' @export
plot.mediation.test <- function(x, filename = NULL, ...) {
  if (!is.null(filename)) {
    file <- R.utils::Arguments$getFilename(filename)
  }
  ## to please R CMD
  xmin <- NULL
  xmax <- NULL
  ymin <- NULL
  ymax <- NULL
  Xn <- NULL
  Yn <- NULL
  ##
  a <- stats::qnorm(seq(.5, 1-x$alpha/2, x$alpha/2))
  df <- tibble::tibble(
                  xmin = c(rep(a, 2), rep(-a, 2)),
                  xmax = c(rep(c(a[-1], Inf), 2), rep(c(-a[-1], -Inf), 2)),
                  ymin = rep(c(a, -a), 2),
                  ymax = rep(c(a[-1], Inf, -a[-1], -Inf), 2)
                )
  scatter_df <- tibble::tibble(Xn = x$t[1], Yn = x$t[2], decision = x$decision)
  ##
  fig <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = df,
                       ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = x$alpha),
                       alpha = 1) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(x = bquote(X[n]), y = bquote(Y[n])) +
    ggplot2::geom_hline(yintercept = 0) + 
    ggplot2::geom_vline(xintercept = 0) + 
    ggplot2::scale_x_continuous(limits = c(-4, 4), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(-4, 4), expand = c(0, 0)) +
    ggplot2::geom_point(data = scatter_df,
                        ggplot2::aes(x = Xn, y = Yn, color = !x$decision)) +
    ggplot2::coord_fixed()
  if (!is.null(filename)) {
    ggplot2::ggsave(fig, file = file)
  } else {
    print(fig)
  }
  invisible()
}

#' Plots  the   output  of  \code{function}  \code{mediation_test}   with  the
#' corresponding 3D power surface.
#' 
#' @param x An output of \code{function} \code{mediation_test}.
#'
#' @param filename  Either \code{NULL} (defaults) or a file  name to create on
#'   disk. If \code{NULL}, then the plot is sent to the default web browser.
#'
#' @param nbins_xy  An \code{integer} (chosen between 50  and 1000, defaulting
#'   to 501), the number of bins used to discretize the along each axis.
#'
#' @param  range_delta_x,range_delta_y Two \code{vector}s  of \code{numeric}s
#'   describing the  ranges of \eqn{\delta_x} and  \eqn{\delta_y}. By default,
#'   they are both set to \code{c(-6, 6)}.
#' 
#' @param ... Not used.
#' 
#' @return Nothing.
#'
#' @examples 
#' (mt <- mediation_test(c(1.1, 2.2), alpha = 1/20))
#' plot3d(mt)
#'
#' @aliases plot3d
#' 
#' @export plot3d
#'
#' @export
R.methodsS3::setMethodS3(
               "plot3d", "mediation.test",
               function(x, filename = NULL, nbins_xy = 501, range_delta_x = c(-6, 6), range_delta_y = c(-6, 6), ...) {
                 if (!is.null(filename)) {
                   file <- R.utils::Arguments$getFilename(filename)
                 }
                 nbins_xy <- R.utils::Arguments$getInteger(nbins_xy, c(50, 1e3))
                 range_delta_x <- R.utils::Arguments$getNumerics(range_delta_x)
                 range_delta_y <- R.utils::Arguments$getNumerics(range_delta_y)
                 if (length(range_delta_x) != 2 | length(range_delta_y) != 2) {
                   R.oo::throw("Arguments 'range_delta_x' and 'range_delta_y' should be two vectors consisting of two real numbers, the ranges of 'delta_x' and 'delta_y', not ",
                               delta_x, delta_y)
                 }
                 ## to please R CMD
                 xmin <- NULL
                 xmax <- NULL
                 ymin <- NULL
                 ymax <- NULL
                 ##
                 compute_power_surface <- Vectorize(
                   function(d_x, d_y, tib) {
                     sum(abs((stats::pnorm(tib$xmax - d_x) -
                              stats::pnorm(tib$xmin - d_x)) * 
                             (stats::pnorm(tib$ymax - d_y) -
                              stats::pnorm(tib$ymin - d_y))))
                   },
                   c("d_x", "d_y"))
                 ##
                 a <- stats::qnorm(seq(.5, 1-x$alpha/2, x$alpha/2))
                 df <- tibble::tibble(
                                 xmin = c(rep(a, 2), rep(-a, 2)),
                                 xmax = c(rep(c(a[-1], Inf), 2), rep(c(-a[-1], -Inf), 2)),
                                 ymin = rep(c(a, -a), 2),
                                 ymax = rep(c(a[-1], Inf, -a[-1], -Inf), 2)
                               )
                 delta_x <- seq(min(range_delta_x), max(range_delta_x), length = nbins_xy)
                 delta_y <- seq(min(range_delta_y), max(range_delta_y), length = nbins_xy)
                 power_surface <- outer(delta_x, delta_y, compute_power_surface, tib = df)
                 ##
                 fig <- plotly::plot_ly(x = ~delta_x, y = ~delta_y, z = ~power_surface)
                 fig <- plotly::add_surface(fig,
                                            contours = list(
                                              z = list(
                                                show = TRUE,
                                                usecolormap = TRUE,
                                                highlightcolor = "#ff0000",
                                                project = list(z = TRUE)
                                              )
                                            ))
                 fig <- plotly::add_trace(fig,
                                          x = x$t[1], y = x$t[2],
                                          z = compute_power_surface(x$t[1], x$t[2], tib = df),
                                          mode = "markers", type = "scatter3d", 
                                          marker = list(
                                            size = 5, color = "red", symbol = 104))
                 fig <- plotly::layout(fig,
                                       scene = list(
                                         camera = list(
                                           eye = list(x = 0, y = -1.75, z = 1.75)
                                         ),
                                         xaxis = list(title = "delta*_x"),
                                         yaxis = list(title = "delta*_y"),
                                         zaxis = list(title = "Rej. prob.")
                                       )
                                       )
                 if (!is.null(filename)) {
                   save_to_disk <- try(plotly::orca(fig, file = file))
                   if (inherits(save_to_disk, "try-error")) {
                     warning("Saving image to disk failed.")
                   }
                 } else {
                   print(fig)
                 }
                 invisible()
               }
             )

#' Carries out the test of the composite null "\eqn{x \times y=0}" against its
#' alternative "\eqn{x  \times y\neq 0}"  based on  the test statistic  in the
#' real plane.
#' 
#' @param t  A  \code{vector}  consisting of  two  \code{numeric}s, the  test
#'   statistic in the real plane.
#'
#' @param alpha A positive \code{numeric}, the wished type-I error, which must
#'   be the inverse of an integer larger  than 2 (defaults to 1/20=5%).  If it
#'   is not the  inverse of an integer, then it  is automatically rounded down
#'   to the closer inverse of an integer.
#' 
#' @return A list, consisting  of: \describe{ \item{t:}{a \code{vector} of two
#'   \code{numeric}s, the test statistic;} \item{alpha:}{a \code{numeric}, the
#'   type-I  error  (possibly  rounded  down  to  the  closer  inverse  of  an
#'   integer);}  \item{decision:}{a \code{logical},  \code{FALSE} if  the null
#'   hypothesis  can be  rejected for  the  alternative at  level 'alpha'  and
#'   \code{TRUE}   otherwise;}   \item{int_pval:}{a   \code{vector}   of   two
#'   \code{numeric}s,  a lower  bound  and  an upper  bound  on the  p-value;}
#'   \item{pval:}{a \code{numeric}, a (random) p-value drawn uniformly between
#'   the two aforementiond bounds.}}
#'
#' @examples
#' (mt <- mediation_test(c(1.1, 2.2), alpha = 1/20))
#' plot(mt)
#'
#' @export
mediation_test <- function(t, alpha = 0.05) {
  t <- R.utils::Arguments$getNumerics(t)
  if (!length(t) == 2) {
    R.oo::throw("Argument 't' should be a vector consisting of two real numbers, the test statistic in the real plane, not ", t) 
  }
  alpha <- R.utils::Arguments$getNumeric(alpha, c(0, 1/2))
  if ((1/alpha) %% 1 != 0) {
    alpha <- 1/floor(1/alpha)
    warning("Argument 'alpha' is not the inverse of an integer larger than 2. It has been rounded to the closer such number. New value of alpha: ",
            alpha, "= 1 /", 1/alpha, ").")
  }
  make_decision <- function(tt, aalpha) {
    qtls <- stats::qnorm(seq(0.5, 1 - aalpha/2, aalpha/2))
    tt1 <- max(abs(tt))
    tt2 <- min(abs(tt))
    if (tt1 %in% qtls | tt2 %in% qtls) {
      decision <- FALSE
    } else {
      decision <- tt2 > qtls[max(which(qtls <= tt1))]
    }
    names(decision) <- "rejection"
    return(decision)
  }
  compute_pval <- function(tt) {
    decision <- TRUE
    alpha_inverse <- 1
    while (decision) {
      pval_upper_bound <- 1/alpha_inverse
      alpha_inverse <- alpha_inverse + 1
      decision <- make_decision(tt, 1/alpha_inverse)
    }
    pval_lower_bound <- 1/alpha_inverse
    pval <- stats::runif(1, pval_lower_bound, pval_upper_bound)
    pvals <- list(lower_bound = pval_lower_bound,
                  pval = pval,
                  upper_bound = pval_upper_bound)
    return(pvals)
  }
  decision <- make_decision(t, alpha)
  pvals <- compute_pval(t)
  
  out <- list(t = t, alpha = alpha, decision = decision,
              pval = pvals$pval,
              int_pval = c(pvals$lower_bound, pvals$upper_bound))
  class(out) <- "mediation.test"
  return(out)
}

