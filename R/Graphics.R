#' Plot KM and Copula PFS estimates side by side with confidence intervals
#'
#' Produces two ggplot2 panels -- one for the Kaplan-Meier estimate and one for
#' the copula-based estimate -- each showing the mean survival curve and its
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
#' \code{Rate +/- CI_Width / 2}.
#'
#' Survival curves are drawn as step functions to match the standard
#' Kaplan-Meier presentation. The CI ribbon is manually expanded into
#' staircase coordinates so it aligns correctly with \code{geom_step()}.
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
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_step
#'   scale_y_continuous scale_x_continuous scale_colour_manual
#'   labs theme_minimal theme element_text element_line element_blank
#'   margin
#' @importFrom patchwork plot_annotation wrap_plots
#' @importFrom scales percent_format
#' @importFrom rlang .data
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

  # ---- build tidy data frames ------------------------------------------------
  df_km <- data.frame(
    time     = seq_len(n_times),
    surv     = km_results$Rate,
    ci_lower = pmax(km_results$Rate - km_results$CI_Width / 2, 0),
    ci_upper = pmin(km_results$Rate + km_results$CI_Width / 2, 1)
  )

  df_copula <- data.frame(
    time     = seq_len(n_times),
    surv     = copula_results$Rate,
    ci_lower = pmax(copula_results$Rate - copula_results$CI_Width / 2, 0),
    ci_upper = pmin(copula_results$Rate + copula_results$CI_Width / 2, 1)
  )

  # ---- convert a data frame to step-ribbon coordinates ----------------------
  # geom_ribbon has no step stat, so each row is duplicated and the second
  # copy's x is advanced by 1, producing a staircase-shaped ribbon that aligns
  # with geom_step().
  to_step_ribbon <- function(df) {
    n   <- nrow(df)
    idx <- rep(seq_len(n), each = 2)
    out <- df[idx, ]
    out$time[seq(2, nrow(out), by = 2)] <- out$time[seq(2, nrow(out), by = 2)] + 1
    out
  }

  # ---- shared theme ----------------------------------------------------------
  pfs_theme <- ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = 14,
                                               margin = ggplot2::margin(b = 8)),
      plot.subtitle    = ggplot2::element_text(size = 11, colour = "grey40",
                                               margin = ggplot2::margin(b = 10)),
      axis.title       = ggplot2::element_text(size = 11),
      axis.text        = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey90"),
      legend.position  = "bottom",
      legend.title     = ggplot2::element_blank()
    )

  # ---- helper to build one panel ---------------------------------------------
  make_panel <- function(df, panel_title, colour_line, colour_ribbon) {

    step_df <- to_step_ribbon(df)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$time)) +
      ggplot2::geom_ribbon(
        data        = step_df,
        ggplot2::aes(x = .data$time, ymin = .data$ci_lower, ymax = .data$ci_upper),
        fill        = colour_ribbon,
        alpha       = 0.25,
        inherit.aes = FALSE
      ) +
      ggplot2::geom_step(
        ggplot2::aes(y = .data$surv, colour = "Estimated PFS"),
        linewidth = 0.9
      )

    # overlay true PFS if provided
    if (!is.null(true_pfs)) {
      true_df <- data.frame(time = seq_len(n_times), surv = true_pfs)
      p <- p +
        ggplot2::geom_step(
          data      = true_df,
          ggplot2::aes(x = .data$time, y = .data$surv, colour = "True PFS"),
          linetype  = "dashed",
          linewidth = 0.8
        )
    }

    p +
      ggplot2::scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2),
        labels = scales::number_format(accuracy = 0.001)
      ) +
      ggplot2::scale_x_continuous(
        breaks = seq_len(n_times)
      ) +
      ggplot2::scale_colour_manual(
        values = c(
          "Estimated PFS" = colour_line,
          "True PFS"      = "grey30"
        )
      ) +
      ggplot2::labs(
        title    = panel_title,
        x        = "Time",
        y        = "Survival Probability"
      ) +
      pfs_theme
  }

  # ---- build each panel ------------------------------------------------------
  p_km     <- make_panel(df_km,     "Kaplan-Meier Estimate", "#2166AC", "#2166AC")
  p_copula <- make_panel(df_copula, "Copula-Based Estimate", "#D6604D", "#D6604D")

  # ---- combine with patchwork ------------------------------------------------
  combined <- patchwork::wrap_plots(p_km, p_copula, ncol = 2)

  if (!is.null(title)) {
    combined <- combined +
      patchwork::plot_annotation(
        title = title,
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5)
        )
      )
  }

  return(combined)
}

#' Plot power curves over time for one or more beta values
#'
#' Produces a ggplot2 figure showing estimated power at each follow-up time
#' point, with one line per value of \code{beta}. An optional horizontal
#' reference line marks the nominal significance level (type I error).
#'
#' @param power_df   Data frame returned by stacking multiple calls to
#'   \code{power_copula_pfs()}, with columns \code{Time}, \code{Power},
#'   \code{MeanDiff}, and \code{beta}. A single-beta data frame (no \code{beta}
#'   column) is also accepted.
#' @param alpha_level Numeric. Nominal significance level drawn as a horizontal
#'   reference line. Set to \code{NULL} to suppress. Default \code{0.05}.
#' @param title      Optional character string used as the plot title.
#' @param subtitle   Optional character string used as the plot subtitle.
#'
#' @return A \code{ggplot} object.
#'
#' @details
#' If \code{power_df} contains a \code{beta} column it is converted to a
#' factor so each level receives a distinct colour from a accessible palette.
#' Points are overlaid on lines to highlight the discrete time grid.
#'
#' @examples
#' \dontrun{
#' betas <- c(-0.2, -0.5, -0.8)
#' power_curve <- lapply(betas, function(b) {
#'   res       <- power_copula_pfs(...)
#'   res$beta  <- b
#'   res
#' })
#' power_curve_df <- do.call(rbind, power_curve)
#' plot_power_curve(power_curve_df)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline
#'   scale_y_continuous scale_x_continuous scale_colour_manual
#'   labs theme_minimal theme element_text element_line element_blank
#'   margin guide_legend
#' @importFrom rlang .data
#' @export
plot_power_curve <- function(power_df,
                             alpha_level = 0.05,
                             title       = NULL,
                             subtitle    = NULL) {

  # ---- input checks -----------------------------------------------------------
  required_cols <- c("Time", "Power")
  if (!all(required_cols %in% names(power_df)))
    stop("power_df must contain at least columns: Time, Power")

  # ---- handle single-beta input -----------------------------------------------
  if (!"beta" %in% names(power_df)) {
    power_df$beta <- "beta"          # single unlabelled group
    single_beta   <- TRUE
  } else {
    single_beta <- FALSE
  }

  power_df$beta <- factor(power_df$beta)

  # ---- colour palette (colourblind-friendly) ----------------------------------
  palette <- c(
    "#2166AC", "#D6604D", "#4DAC26",
    "#8073AC", "#E08214", "#01665E"
  )
  n_betas <- nlevels(power_df$beta)

  colours <- palette[seq_len(n_betas)]
  names(colours) <- levels(power_df$beta)

  # ---- shared theme (mirrors Graphics.R) -------------------------------------
  power_theme <- ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = 14,
                                               margin = ggplot2::margin(b = 8)),
      plot.subtitle    = ggplot2::element_text(size = 11, colour = "grey40",
                                               margin = ggplot2::margin(b = 10)),
      axis.title       = ggplot2::element_text(size = 11),
      axis.text        = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey90"),
      legend.position  = if (single_beta) "none" else "bottom",
      legend.title     = ggplot2::element_text(size = 10, face = "bold")
    )

  # ---- build plot -------------------------------------------------------------
  p <- ggplot2::ggplot(
    power_df,
    ggplot2::aes(
      x      = .data$Time,
      y      = .data$Power,
      colour = .data$beta,
      group  = .data$beta
    )
  ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2.5)

  # reference line at alpha_level
  if (!is.null(alpha_level)) {
    p <- p + ggplot2::geom_hline(
      yintercept = alpha_level,
      linetype   = "dashed",
      linewidth  = 0.6,
      colour     = "grey50"
    )
  }

  # --- define labels before the plot ---
  beta_labels <- paste0("\u03b2 = ", levels(power_df$beta))
  names(beta_labels) <- levels(power_df$beta)

  # --- then the plot block ---
  p <- p +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      labels = scales::number_format(accuracy = 0.01)
    ) +
    ggplot2::scale_x_continuous(
      breaks = unique(power_df$Time)
    ) +
    ggplot2::scale_colour_manual(
      values = colours,
      labels = if (single_beta) NULL else beta_labels,
      guide  = ggplot2::guide_legend(title = "\u03b2")
    ) +
    ggplot2::labs(
      title    = title,
      subtitle = subtitle,
      x        = "Time",
      y        = "Power"
    ) +
    power_theme

  return(p)
}





