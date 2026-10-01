############################################################
# Marginal-effects helper functions
#
# meineq() and totalme()
#   Build the `hypothesis` function for avg_comparisons() that returns
#   the ME inequality or the Total ME from Mize & Han (2025) for every
#   variable named in `variables`. Nominal variables use "pairwise";
#   the weighted ME inequality is the default. Weights are matched to
#   the contrasts by name, so row order never matters. Works with
#   multi-outcome models, by =, and suest() objects (each model's
#   weights come from its own estimation sample). As in Stata, the
#   category shares use the model's weights (survey design or fitting
#   weights) when it has them; model_weights = FALSE ignores them.
#
#   avg_comparisons(model, variables = list(race4 = "pairwise"),
#                   hypothesis = meineq(model))
#   avg_comparisons(model, variables = list(age = "sd", race4 = "pairwise"),
#                   hypothesis = totalme(model))
#
# meineq_weights()
#   Prevalence weights for all pairwise comparisons of a nominal
#   predictor, as a named vector ordered to match marginaleffects
#   pairwise contrast output.
#
#   weights <- meineq_weights(model, race4)
#
# The model is used instead of the raw data so weights are based
# on the estimation sample.
############################################################

meineq <- function(model, weighted = TRUE, model_weights = TRUE) {
  .me_summary(model, weighted, total = FALSE, model_weights)
}

totalme <- function(model, weighted = TRUE, model_weights = TRUE) {
  .me_summary(model, weighted, total = TRUE, model_weights)
}

.me_summary <- function(model, weighted, total, model_weights) {
  is_suest <- inherits(model, "suest_model")
  models <- if (is_suest) model$models else list(model = model)
  cache <- new.env()

  # Weight of each pair of categories, keyed by both "a - b" and "b - a"
  pair_w <- function(mod, var) {
    key <- paste(mod, var, sep = "\r")
    if (is.null(cache[[key]])) {
      p <- .meineq_shares(models[[mod]], var, NULL, model_weights)
      pr <- utils::combn(names(p), 2)
      w <- if (weighted) (p[pr[1, ]] + p[pr[2, ]]) / (length(p) - 1) else rep(1 / ncol(pr), ncol(pr))
      w <- c(stats::setNames(w, paste(pr[2, ], "-", pr[1, ])), stats::setNames(w, paste(pr[1, ], "-", pr[2, ])))
      if (anyDuplicated(names(w))) stop(sprintf("Category labels of '%s' are ambiguous.", var), call. = FALSE)
      cache[[key]] <- w
    }
    cache[[key]]
  }

  function(x) {
    x <- as.data.frame(x)
    by <- setdiff(names(x)[seq_len(match("estimate", names(x)))], c("term", "group", "contrast", "estimate"))
    g <- if (is.null(x$group)) rep("", nrow(x)) else as.character(x$group)
    mod <- if (is_suest) sub("::.*", "", g) else rep("model", nrow(x))
    if (total && (if (is_suest) !any(grepl("::", g)) else all(g == ""))) {
      stop("totalme() needs a model with several outcome categories.", call. = FALSE)
    }

    w <- rep(NA_real_, nrow(x))
    mt <- paste(mod, x$term, sep = "\r")
    for (k in unique(mt)) {
      i <- mt == k
      w[i] <- if (total && length(unique(x$contrast[i])) == 1L) 1 else pair_w(mod[i][1], x$term[i][1])[x$contrast[i]]
    }
    if (anyNA(w)) stop("Some contrasts do not match the categories in the estimation sample.", call. = FALSE)

    grp <- if (total) (if (is_suest) mod else NULL) else if (is.null(x$group)) NULL else g
    key <- do.call(paste, c(list(x$term), list(grp), x[by], sep = "\r"))
    key <- factor(key, levels = unique(key))
    a <- abs(x$estimate)
    if (total) {
      kc <- factor(paste(key, x$contrast, sep = "\r"), levels = unique(paste(key, x$contrast, sep = "\r")))
      a <- as.numeric(tapply(a, kc, sum)) / 2
      first <- !duplicated(kc)
      w <- w[first]
      key <- key[first]
      x <- x[first, , drop = FALSE]
      if (!is.null(grp)) grp <- grp[first]
    }
    bad <- abs(tapply(w, key, sum) - 1) > 1e-8
    if (any(bad)) {
      stop(sprintf("Use \"pairwise\" for '%s': not every pairwise contrast is present.", x$term[match(names(bad)[bad][1], key)]), call. = FALSE)
    }

    first <- !duplicated(key)
    out <- data.frame(term = x$term[first])
    if (!is.null(grp)) out$group <- grp[first]
    out[by] <- x[first, by, drop = FALSE]
    out$estimate <- as.numeric(tapply(w * a, key, sum))
    out
  }
}

meineq_weights <- function(model, variable, model_weights = TRUE) {
  p <- .meineq_shares(model, substitute(variable), parent.frame(), model_weights)
  k <- length(p)
  pairs <- utils::combn(seq_along(p), 2)
  weights <- (as.numeric(p[pairs[1, ]]) + as.numeric(p[pairs[2, ]])) / (k - 1)
  names(weights) <- paste(names(p)[pairs[2, ]], "-", names(p)[pairs[1, ]])

  # marginaleffects uses keyed grouping for pairwise contrasts, which
  # sorts character contrast labels in C-locale order. Base R's radix
  # ordering uses the same C-locale collation.
  weights <- weights[order(names(weights), method = "radix")]
  weights / sum(weights)
}

# Category shares of `expr` in the model's estimation sample
.meineq_shares <- function(model, expr, env, model_weights = TRUE) {
  if (!requireNamespace("insight", quietly = TRUE)) {
    stop("Package 'insight' is required.", call. = FALSE)
  }
  if (!isTRUE(model_weights) && !isFALSE(model_weights)) {
    stop("`model_weights` must be TRUE or FALSE.", call. = FALSE)
  }

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
  # A name held in an object (e.g., a loop variable) is used as the variable name
  if (!var %in% names(data)) {
    val <- tryCatch(eval(expr, env), error = function(e) NULL)
    if (is.character(val) && length(val) == 1L) var <- val
  }
  if (!var %in% names(data)) {
    stop(sprintf("Variable '%s' was not found in the model data.", var), call. = FALSE)
  }
  if (is.numeric(data[[var]]) && !var %in% attr(data, "factors")) {
    stop(sprintf("Variable '%s' is numeric in the model. Fit it as a factor for pairwise contrasts.", var), call. = FALSE)
  }

  # As in Stata, a model fit with weights (survey design or fitting
  # weights) gives weighted category shares; unweighted models do not.
  is_survey <- any(grepl("^(svy|svrep)", class(model)))

  x <- data[[var]]
  obs_weights <- NULL

  if (model_weights) {
    # For survey-family models, the fitted object retains the survey design.
    # survey::weights(..., type = "sampling") is the authoritative source
    # for both ordinary and replicate-weight survey designs.
    if (isTRUE(is_survey) && !is.null(model$survey.design)) {
      obs_weights <- tryCatch(
        stats::weights(model$survey.design, type = "sampling"),
        error = function(e) NULL
      )
    }

    # Fitting weights of ordinary models (NULL for an unweighted model).
    if (is.null(obs_weights)) {
      obs_weights <- tryCatch(
        suppressWarnings(insight::get_weights(model, source = "frame")),
        error = function(e) NULL
      )
    }
    if (!is.null(obs_weights)) {
      obs_weights <- as.numeric(obs_weights)
      if (length(obs_weights) != length(x)) {
        stop("Recovered model weights do not align with the model estimation data. Use `model_weights = FALSE` for unweighted category shares.", call. = FALSE)
      }
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

  if (length(p) < 2L) {
    stop(sprintf("Variable '%s' must have at least two observed categories.", var), call. = FALSE)
  }
  p <- stats::setNames(as.numeric(p), names(p))
  attr(p, "var") <- var
  p
}
