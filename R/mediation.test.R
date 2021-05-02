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
                     sprintf("can reject the null for its alternative with confidence %1.3f", x$alpha),
                     sprintf("cannot reject the null for its alternative with confidence %1.3f", x$alpha))
  print(decision)
  cat("* approximate p-value:\n")
  print(x$pval)
  invisible()
}

#' Plots the output of \code{function} \code{mediation_test}.
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
#' @return  A list,  consisting  of  the test  statistic,  the type-I  error
#'   (possibly rounded down  to the closer inverse of an  integer), a logical,
#'   \code{FALSE} if the  null hypothesis can be rejected  for the alternative
#'   at level 'alpha' and \code{TRUE} otherwise, and an approximate p-value.
#'
#' @details Note that the p-value is random.  
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
  compute_pval <- function(tt, aalpha) {
    decision <- TRUE
    alpha_inverse <- 1
    while (decision) {
      alpha_inverse <- alpha_inverse + 1
      decision <- make_decision(tt, 1/alpha_inverse)
    }
    pval <- stats::runif(1, aalpha, 1/(alpha_inverse - 1))
    names(pval) <- "p-value"
    return(pval)
  }
  decision <- make_decision(t, alpha)
  pval <- compute_pval(t, alpha)
  out <- list(t = t, alpha = alpha, decision = decision, pval = pval)
  class(out) <- "mediation.test"
  return(out)
}

