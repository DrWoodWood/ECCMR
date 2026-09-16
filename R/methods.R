## S3 methods for objects of class "eccmr".

#' @export
print.eccmr <- function(x, digits = 4, ...) {
  cat("ECC-MR: error-correcting code inspired Mendelian randomization\n")
  cat(sprintf("  instruments: %d | parity checks: %d\n",
              x$n_instruments, x$n_checks))
  cat(sprintf("  causal effect theta = %.*f", digits, x$theta))
  if (!is.na(x$se)) {
    cat(sprintf("  (SE %.*f, 95%% CI [%.*f, %.*f])",
                digits, x$se, digits, x$ci[1L], digits, x$ci[2L]))
  }
  cat("\n")
  cat(sprintf("  pleiotropic instruments detected: %d\n",
              length(x$pleiotropic_index)))
  cat(sprintf("  syndrome test: T = %.3f, df = %d, p = %.3g\n",
              x$syndrome$statistic, x$syndrome$df, x$syndrome$p_value))
  cat(sprintf("  lambda = %.4g, gamma = %.4g | converged: %s (%d iterations)\n",
              x$lambda, x$gamma, x$fit$converged, x$fit$iter))
  invisible(x)
}

#' @export
summary.eccmr <- function(object, digits = 4, top_n = 10L, ...) {
  cat("ECC-MR summary\n")
  cat("==============\n")
  cat(sprintf("Causal estimate: theta = %.*f\n", digits, object$theta))
  if (!is.na(object$se)) {
    cat(sprintf("  SE = %.*f, 95%% percentile CI = [%.*f, %.*f] (%d bootstrap replicates)\n",
                digits, object$se, digits, object$ci[1L], digits, object$ci[2L],
                length(object$boot_theta)))
  }
  cat(sprintf("Instruments: %d | parity checks: %d | pleiotropic detected: %d\n",
              object$n_instruments, object$n_checks,
              length(object$pleiotropic_index)))
  cat(sprintf("Global syndrome test of pleiotropy: T = %.3f, df = %d, p = %.3g\n",
              object$syndrome$statistic, object$syndrome$df,
              object$syndrome$p_value))
  cat(sprintf("Tuning: lambda = %.4g, gamma = %.4g\n",
              object$lambda, object$gamma))
  if (length(object$pleiotropic) > 0L) {
    cat("\nTop pleiotropic instruments (by |alpha|):\n")
    ord <- order(abs(object$pleiotropic), decreasing = TRUE)
    show <- utils::head(ord, top_n)
    tab <- data.frame(
      instrument = if (is.null(names(object$pleiotropic)))
        object$pleiotropic_index[show] else names(object$pleiotropic)[show],
      alpha = signif(object$pleiotropic[show], digits)
    )
    print(tab, row.names = FALSE)
  }
  out <- list(theta = object$theta, se = object$se, ci = object$ci,
              syndrome = object$syndrome, lambda = object$lambda,
              gamma = object$gamma, pleiotropic = object$pleiotropic)
  class(out) <- "summary.eccmr"
  invisible(out)
}

#' @export
print.summary.eccmr <- function(x, ...) invisible(x)

#' @export
coef.eccmr <- function(object, ...) {
  c(theta = object$theta)
}

#' @export
plot.eccmr <- function(x, ...) {
  z <- x$syndrome$z
  if (length(z) == 0L) {
    message("No parity checks available to plot.")
    return(invisible(x))
  }
  zc <- stats::qnorm(0.975)
  col <- ifelse(abs(z) > zc, "#B2182B", "#2166AC")
  graphics::plot(z, pch = 16, cex = 0.7, col = col,
                 xlab = "Parity-check index",
                 ylab = "Standardized syndrome",
                 main = "ECC-MR syndrome diagnostics", ...)
  graphics::abline(h = c(-zc, 0, zc), lty = c(2, 1, 2),
                   col = c("grey50", "grey30", "grey50"))
  invisible(x)
}
