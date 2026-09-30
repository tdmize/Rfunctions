############################################################
# Marginal-effects helper functions
#
# meineq_weights()
#   Computes prevalence weights for all pairwise comparisons
#   of a nominal predictor, as used in weighted marginal-effect
#   inequality calculations from Mize & Han (2025).
#
# Primary use:
#   weights <- meineq_weights(model, race3)
#
# The model is used instead of the raw data so weights are based
# on the estimation sample. The function returns a named numeric
# vector ordered to match marginaleffects pairwise contrast output.
############################################################

meineq_weights <- function(model, variable, model_weights = NULL) {
  if (!requireNamespace("insight", quietly = TRUE)) {
    stop("Package 'insight' is required.", call. = FALSE)
  }
  if (!is.null(model_weights) &&
      (!is.logical(model_weights) || length(model_weights) != 1L || is.na(model_weights))) {
    stop("`model_weights` must be TRUE, FALSE, or NULL.", call. = FALSE)
  }

  expr <- substitute(variable)
  if (is.character(expr) && length(expr) == 1L) {
    var <- expr
  } else {
    vars <- all.vars(expr)
    if (length(vars) != 1L) {
      stop("`variable` must identify one variable in the model.", call. = FALSE)
    }
    var <- vars
  }

  data <- tryCatch(
    insight::get_data(model, source = "frame", verbose = FALSE),
    error = function(e) NULL
  )

  if (is.null(data) || !is.data.frame(data)) {
    stop("Could not recover the model estimation data.", call. = FALSE)
  }
  if (!var %in% names(data)) {
    stop(sprintf("Variable '%s' was not found in the model data.", var), call. = FALSE)
  }

  # Survey-model fits use survey-weighted category proportions by default.
  # For ordinary weighted models, fitting weights are used only when the
  # user explicitly requests model_weights = TRUE.
  is_survey <- any(grepl("^(svy|svrep)", class(model)))
  use_weights <- if (is.null(model_weights)) is_survey else model_weights

  x <- data[[var]]
  obs_weights <- NULL

  if (isTRUE(use_weights)) {
    # For survey-family models, the fitted object retains the survey design.
    # survey::weights(..., type = "sampling") is the authoritative source
    # for both ordinary and replicate-weight survey designs.
    if (isTRUE(is_survey) && !is.null(model$survey.design)) {
      obs_weights <- tryCatch(
        stats::weights(model$survey.design, type = "sampling"),
        error = function(e) NULL
      )
    }

    # General fallback for ordinary weighted models and other supported classes.
    if (is.null(obs_weights)) {
      obs_weights <- tryCatch(
        insight::get_weights(model, source = "frame"),
        error = function(e) NULL,
        warning = function(w) NULL
      )
    }
    if (is.null(obs_weights)) {
      stop("Could not recover model weights. Use `model_weights = FALSE` for unweighted category proportions.", call. = FALSE)
    }
    obs_weights <- as.numeric(obs_weights)
    if (length(obs_weights) != length(x)) {
      stop("Recovered model weights do not align with the model estimation data.", call. = FALSE)
    }
  }

  keep <- !is.na(x)
  if (!is.null(obs_weights)) keep <- keep & !is.na(obs_weights)
  x <- x[keep]
  if (!is.null(obs_weights)) obs_weights <- obs_weights[keep]

  if (!length(x)) {
    stop(sprintf("Variable '%s' has no observed values in the estimation sample.", var), call. = FALSE)
  }

  if (is.factor(x)) {
    x <- droplevels(x)
  } else {
    x <- factor(x, levels = sort(unique(x)))
  }

  if (is.null(obs_weights)) {
    p <- prop.table(table(x))
  } else {
    if (any(!is.finite(obs_weights)) || any(obs_weights < 0)) {
      stop("Model weights must be finite and nonnegative.", call. = FALSE)
    }
    if (sum(obs_weights) <= 0) {
      stop("Model weights must have a positive sum.", call. = FALSE)
    }
    p <- tapply(obs_weights, x, sum)
    p <- p / sum(p)
  }

  k <- length(p)
  if (k < 2L) {
    stop(sprintf("Variable '%s' must have at least two observed categories.", var), call. = FALSE)
  }

  pairs <- utils::combn(seq_along(p), 2)
  weights <- (as.numeric(p[pairs[1, ]]) + as.numeric(p[pairs[2, ]])) / (k - 1)
  names(weights) <- paste(names(p)[pairs[2, ]], "-", names(p)[pairs[1, ]])

  # marginaleffects uses keyed grouping for pairwise contrasts, which
  # sorts character contrast labels in C-locale order. Base R's radix
  # ordering uses the same C-locale collation.
  weights <- weights[order(names(weights), method = "radix")]
  weights / sum(weights)
}
