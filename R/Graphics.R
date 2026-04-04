library(ggplot2)

#' Plot KM and Copula PFS estimates side by side with confidence intervals
#'
#' Produces two ggplot2 panels — one for the Kaplan–Meier estimate and one for
#' the copula-based estimate — each showing the mean survival curve and its
#' pointwise confidence interval across simulation iterations.
#'
#' @param km_results     Data frame returned by \code{run_iterations()}, with
#'   columns \code{Rate}, \code{CI_Width}, and \code{Coverage}.
#' @param copula_results Data frame returned by \code{run_copula_iterations()},
#'   with the same column structure.
#' @param true_pfs       Optional numeric vector of true PFS values (length
#'   \code{n_times}). When supplied, the true curve is overlaid on both panels
#'   as a dashed reference line.
#' @param title          Optional character string used as an overall plot title.
#'
#' @return A \code{ggplot} object (patchwork-combined two-panel figure).
#'
#' @details
#' Both \code{km_results} and \code{copula_results} must have the same number
#' of rows (one per time point). The CI bands are reconstructed as
#' \code{Rate ± CI_Width / 2}.
#'
#' @examples
#' \dontrun{
#' km_res     <- run_iterations(...)
#' copula_res <- run_copula_iterations(...)
#' true_pfs   <- get_true_rates(...)$surv
#'
#' plot_pfs_estimates(km_res, copula_res, true_pfs = true_pfs)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_point
#'   scale_y_continuous scale_x_continuous labs theme_minimal theme
#'   element_text element_line element_blank
#' @importFrom patchwork plot_annotation wrap_plots
#' @export
plot_pfs_estimates <- function(km_results,
                               copula_results,
                               true_pfs = NULL,
                               title    = NULL) {

  # ---- input checks ----------------------------------------------------------
  required_cols <- c("Rate", "CI_Width", "Coverage")

  if (!all(required_cols %in% names(km_results)))
    stop("km_results must have columns: Rate, CI_Width, Coverage")

  if (!all(required_cols %in% names(copula_results)))
    stop("copula_results must have columns: Rate, CI_Width, Coverage")

  if (nrow(km_results) != nrow(copula_results))
    stop("km_results and copula_results must have the same number of rows (time points).")

  n_times <- nrow(km_results)

  if (!is.null(true_pfs) && length(true_pfs) != n_times)
    stop("true_pfs must have length equal to the number of time points.")

  df_km     <- make_df(km_results,     "Kaplan-Meier")
  df_copula <- make_df(copula_results, "Copula")

  # clamp CIs to [0, 1]
  df_km$ci_lower     <- pmax(df_km$ci_lower,     0)
  df_km$ci_upper     <- pmin(df_km$ci_upper,     1)
  df_copula$ci_lower <- pmax(df_copula$ci_lower, 0)
  df_copula$ci_upper <- pmin(df_copula$ci_upper, 1)

  # ---- shared theme ----------------------------------------------------------
  pfs_theme <- theme_minimal(base_size = 13) +
    theme(
      plot.title       = element_text(face = "bold", size = 14, margin = margin(b = 8)),
      plot.subtitle    = element_text(size = 11, colour = "grey40", margin = margin(b = 10)),
      axis.title       = element_text(size = 11),
      axis.text        = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey90"),
      legend.position  = "bottom",
      legend.title     = element_blank()
    )

  # ---- helper to build one panel ---------------------------------------------
  make_panel <- function(df, panel_title, colour_line, colour_ribbon) {

    p <- ggplot(df, aes(x = time)) +
      geom_ribbon(
        aes(ymin = ci_lower, ymax = ci_upper),
        fill    = colour_ribbon,
        alpha   = 0.25
      ) +
      geom_line(
        aes(y = surv, colour = "Estimated PFS"),
        linewidth = 0.9
      ) +
      geom_point(
        aes(y = surv, colour = "Estimated PFS"),
        size = 2
      )

    # overlay true PFS if provided
    if (!is.null(true_pfs)) {
      true_df <- data.frame(time = seq_len(n_times), surv = true_pfs)
      p <- p +
        geom_line(
          data      = true_df,
          aes(x = time, y = surv, colour = "True PFS"),
          linetype  = "dashed",
          linewidth = 0.8
        )
    }

    p +
      scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2),
        labels = scales::percent_format(accuracy = 1)
      ) +
      scale_x_continuous(
        breaks = seq_len(n_times)
      ) +
      scale_colour_manual(
        values = c(
          "Estimated PFS" = colour_line,
          "True PFS"      = "grey30"
        )
      ) +
      labs(
        title    = panel_title,
        subtitle = sprintf("Mean ± 95%% CI across simulation iterations"),
        x        = "Time",
        y        = "Survival Probability"
      ) +
      pfs_theme
  }

  # ---- build each panel ------------------------------------------------------
  p_km     <- make_panel(df_km,     "Kaplan–Meier Estimate",  "#2166AC", "#2166AC")
  p_copula <- make_panel(df_copula, "Copula-Based Estimate",  "#D6604D", "#D6604D")

  # ---- combine with patchwork ------------------------------------------------
  combined <- patchwork::wrap_plots(p_km, p_copula, ncol = 2)

  if (!is.null(title)) {
    combined <- combined +
      patchwork::plot_annotation(
        title = title,
        theme = theme(
          plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
        )
      )
  }

  return(combined)
}
