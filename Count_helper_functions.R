############################################################
# Count model helper functions
#
# countfit(): compares Poisson (PRM), negative binomial (NBRM),
# zero-inflated Poisson (ZIP), and zero-inflated negative binomial
# (ZINB) models, like Stata's countfit. By default the zero-inflated
# models use all predictors in both the count and inflation parts;
# use inflate = ~ x1 + x2 to choose the inflation predictors, and
# models = c("PRM", "NBRM", "ZIP") to compare only some of the models.
#
#   source("https://raw.githubusercontent.com/tdmize/Rfunctions/main/Count_helper_functions.R")
#   cf <- countfit(y ~ x1 + x2, data = mydata)
#   cf$pic   # observed - predicted graph (a ggplot object)
#
# Required packages: MASS, pscl, ggplot2
############################################################

countfit <- function(formula, data, inflate = NULL, maxcount = 9,
                     models = c("PRM", "NBRM", "ZIP", "ZINB")) {
  for (p in c("MASS", "pscl", "ggplot2")) {
    if (!requireNamespace(p, quietly = TRUE))
      stop("countfit() needs the ", p, " package: install.packages(\"", p, "\")", call. = FALSE)
  }
  if (!inherits(formula, "formula") || length(formula) != 3)
    stop("countfit() needs a formula such as countfit(y ~ x1 + x2, data = mydata)", call. = FALSE)
  if ("|" %in% all.names(formula[[3]]))
    stop("Give the inflation predictors with inflate = ~ x1 + x2, not with | in the formula", call. = FALSE)
  if (!is.null(inflate) && (!inherits(inflate, "formula") || length(inflate) != 2))
    stop("inflate must be a one-sided formula such as inflate = ~ x1 + x2", call. = FALSE)
  if (!is.numeric(maxcount) || length(maxcount) != 1 || maxcount < 0 || maxcount != round(maxcount))
    stop("maxcount must be a single whole number of 0 or more", call. = FALSE)
  alias <- c(prm = "PRM", poisson = "PRM", nbrm = "NBRM", nbreg = "NBRM", negbin = "NBRM",
    zip = "ZIP", zinb = "ZINB")
  if (!is.character(models) || length(models) == 0 || anyNA(alias[tolower(models)]))
    stop("models can include \"PRM\", \"NBRM\", \"ZIP\", and \"ZINB\"", call. = FALSE)
  sel <- intersect(c("PRM", "NBRM", "ZIP", "ZINB"), alias[tolower(models)])

  # One estimation sample for all models: complete cases on every model variable
  zformula <- sformula <- formula
  if (!is.null(inflate)) {
    zformula[[3]] <- call("|", formula[[3]], inflate[[2]])
    sformula[[3]] <- call("+", formula[[3]], inflate[[2]])
  }
  data <- as.data.frame(data)
  vars <- intersect(all.vars(terms(sformula, data = data)), names(data))
  data <- data[complete.cases(data[vars]), , drop = FALSE]
  if (nrow(data) == 0) stop("No observations remain after dropping missing data", call. = FALSE)
  y <- model.response(model.frame(formula, data))
  if (!is.numeric(y) || is.factor(y))
    stop("The outcome must be a numeric count variable, not a factor or text", call. = FALSE)
  y <- as.numeric(y)
  if (any(y < 0) || any(y != round(y)))
    stop("The outcome must be a count: whole numbers of 0 or more", call. = FALSE)
  if (all(y == 0)) stop("The outcome is 0 for every observation", call. = FALSE)
  if (!any(y == 0) && any(c("ZIP", "ZINB") %in% sel))
    stop("The outcome has no zeros, so the zero-inflated models cannot be estimated; ",
      "use models = c(\"PRM\", \"NBRM\")", call. = FALSE)

  # Fit each model; stop with a clear message if any fail or do not converge
  zctl <- pscl::zeroinfl.control(reltol = 1e-12)
  fitters <- list(
    PRM  = function() glm(formula, family = poisson, data = data, na.action = na.omit),
    NBRM = function() MASS::glm.nb(formula, data = data, na.action = na.omit),
    ZIP  = function() pscl::zeroinfl(zformula, data = data, dist = "poisson", na.action = na.omit, control = zctl),
    ZINB = function() pscl::zeroinfl(zformula, data = data, dist = "negbin", na.action = na.omit, control = zctl))
  mods <- list()
  why <- character()
  for (m in sel) {
    fit <- tryCatch(fitters[[m]](), error = function(e) conditionMessage(e))
    if (is.character(fit)) why[m] <- fit
    else if (isFALSE(fit$converged) || !is.finite(logLik(fit)) || !all(is.finite(na.omit(coef(fit)))))
      why[m] <- "did not converge"
    else mods[[m]] <- fit
  }
  if (length(why)) {
    ok <- setdiff(sel, names(why))
    stop("The ", paste(names(why), collapse = " and "), " model", if (length(why) > 1) "s",
      " did not converge (", paste0(names(why), ": ", why, collapse = "; "), ").",
      if (length(ok)) paste0(" Run countfit() again without ", if (length(why) > 1) "them" else "it",
        ":\n  models = c(", paste0("\"", ok, "\"", collapse = ", "), ")"),
      "\nIf several models fail, check for predictors that do not vary or that perfectly predict the outcome.",
      call. = FALSE)
  }

  # Observed proportions and average predicted probabilities of each count
  k <- 0:maxcount
  obs <- sapply(k, function(j) mean(y == j))
  zprob <- function(m) colMeans(predict(m, type = "prob", at = 0:max(1, maxcount)))[seq_along(k)]  # pscl errors with one count
  probfun <- list(
    PRM  = function(m) sapply(k, function(j) mean(dpois(j, fitted(m)))),
    NBRM = function(m) sapply(k, function(j) mean(dnbinom(j, size = m$theta, mu = fitted(m)))),
    ZIP  = zprob, ZINB = zprob)
  pred <- matrix(unlist(lapply(sel, function(m) probfun[[m]](mods[[m]]))),
    nrow = length(k), dimnames = list(NULL, sel))
  dif <- obs - pred

  compare <- data.frame(
    MaxDiff = apply(dif, 2, function(d) d[which.max(abs(d))]),
    AtValue = apply(dif, 2, function(d) k[which.max(abs(d))]),
    MeanAbsDiff = colMeans(abs(dif)), row.names = sel)

  # Fit statistics and pairwise comparisons; dif = first model - second model
  bic <- sapply(mods, BIC)
  aic <- sapply(mods, AIC)
  ll <- sapply(mods, function(m) as.numeric(logLik(m)))
  fit <- data.frame(BIC = bic, AIC = aic, row.names = sel)
  tests <- lrtests <- NULL
  if (length(sel) > 1) {
    pairs <- combn(sel, 2)
    a <- pairs[1, ]
    b <- pairs[2, ]
    tests <- data.frame(Comparison = paste(a, "vs", b),
      BICdif = round(bic[a] - bic[b], 3), BICprefer = ifelse(bic[a] > bic[b], b, a),
      Evidence = cut(abs(bic[a] - bic[b]), c(0, 2, 6, 10, Inf), include.lowest = TRUE,
        labels = c("Weak", "Positive", "Strong", "Very strong")),
      AICdif = round(aic[a] - aic[b], 3), AICprefer = ifelse(aic[a] > aic[b], b, a))
    nest <- paste(a, b) %in% c("PRM NBRM", "ZIP ZINB")
    if (any(nest)) {
      lr <- pmax(2 * (ll[b[nest]] - ll[a[nest]]), 0)
      lrtests <- data.frame(Comparison = paste(a[nest], "vs", b[nest]),
        LRchi2 = round(lr, 3), p = round(pchisq(lr, df = 1, lower.tail = FALSE) / 2, 3))
    }
  }

  pdat <- data.frame(count = k, dif = as.vector(dif),
    model = factor(rep(sel, each = length(k)), levels = sel))
  pic <- ggplot2::ggplot(pdat, ggplot2::aes(count, dif, color = model)) +
    ggplot2::geom_hline(yintercept = 0) + ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::scale_x_continuous(breaks = k) +
    ggplot2::labs(x = "Count", y = "Observed - Predicted", color = NULL,
      caption = "Positive deviations show underpredictions") +
    ggplot2::theme_minimal()

  out <- structure(list(models = mods, compare = compare, fit = fit, tests = tests, lrtests = lrtests,
    probs = data.frame(Count = k, Actual = obs, pred), pic = pic, N = length(y)),
    class = "countfit")
  print(out)
  invisible(out)
}

print.countfit <- function(x, ...) {
  cat("\nComparison of Mean Observed and Predicted Count (N = ", x$N, ")\n\n", sep = "")
  print(round(x$compare, 3))
  cat("\nTests and Fit Statistics\n\n")
  print(format(round(x$fit, 3), nsmall = 3))
  if (!is.null(x$tests)) {
    cat("\n")
    print(x$tests, row.names = FALSE)
    cat("\nEvidence is based on the BIC difference\n")
  }
  if (!is.null(x$lrtests)) {
    cat("\nLikelihood ratio tests\n\n")
    print(format(x$lrtests, nsmall = 3), row.names = FALSE)
  }
  invisible(x)
}
