############################################################
# Count model helper functions
#
# countfit(): compares Poisson (PRM), negative binomial (NBRM),
# zero-inflated Poisson (ZIP), and zero-inflated negative binomial
# (ZINB) models, like Stata's countfit. The zero-inflated models
# use all predictors in both the count and inflation parts.
#
#   source("https://raw.githubusercontent.com/tdmize/Rfunctions/main/Count_helper_functions.R")
#   cf <- countfit(y ~ x1 + x2, data = mydata)
#   cf$pic   # observed - predicted graph (a ggplot object)
#
# Required packages: MASS, pscl, ggplot2
############################################################

countfit <- function(formula, data, maxcount = 9) {
  y <- model.response(model.frame(formula, data))
  mods <- list(
    PRM  = glm(formula, family = poisson, data = data),
    NBRM = MASS::glm.nb(formula, data = data),
    ZIP  = pscl::zeroinfl(formula, data = data, dist = "poisson"),
    ZINB = pscl::zeroinfl(formula, data = data, dist = "negbin"))

  # Observed proportions and average predicted probabilities of each count
  k <- 0:maxcount
  obs <- sapply(k, function(j) mean(y == j))
  pred <- cbind(
    PRM  = sapply(k, function(j) mean(dpois(j, fitted(mods$PRM)))),
    NBRM = sapply(k, function(j) mean(dnbinom(j, size = mods$NBRM$theta, mu = fitted(mods$NBRM)))),
    ZIP  = colMeans(predict(mods$ZIP, type = "prob", at = k)),
    ZINB = colMeans(predict(mods$ZINB, type = "prob", at = k)))
  dif <- obs - pred

  compare <- data.frame(
    MaxDiff = apply(dif, 2, function(d) d[which.max(abs(d))]),
    AtValue = apply(dif, 2, function(d) k[which.max(abs(d))]),
    MeanAbsDiff = colMeans(abs(dif)))

  # Fit statistics and pairwise comparisons; dif = first model - second model
  bic <- sapply(mods, BIC)
  aic <- sapply(mods, AIC)
  ll <- sapply(mods, function(m) as.numeric(logLik(m)))
  fit <- data.frame(BIC = bic, AIC = aic)
  a <- c("PRM", "PRM", "PRM", "NBRM", "NBRM", "ZIP")
  b <- c("NBRM", "ZIP", "ZINB", "ZIP", "ZINB", "ZINB")
  lr <- ifelse(paste(a, b) %in% c("PRM NBRM", "ZIP ZINB"), 2 * (ll[b] - ll[a]), NA)
  tests <- data.frame(Comparison = paste(a, "vs", b),
    BICdif = round(bic[a] - bic[b], 3), AICdif = round(aic[a] - aic[b], 3),
    Prefer = ifelse(bic[a] > bic[b], b, a),
    Evidence = cut(abs(bic[a] - bic[b]), c(0, 2, 6, 10, Inf),
      labels = c("Weak", "Positive", "Strong", "Very strong")),
    LRchi2 = ifelse(is.na(lr), "", sprintf("%.3f", lr)),
    p = ifelse(is.na(lr), "", sprintf("%.3f", pchisq(lr, df = 1, lower.tail = FALSE) / 2)))

  cat("\nComparison of Mean Observed and Predicted Count\n\n")
  print(round(compare, 3))
  cat("\nTests and Fit Statistics\n\n")
  print(round(fit, 3))
  cat("\n")
  print(tests, row.names = FALSE)
  cat("\nEvidence is based on the BIC difference; LR tests are for nested models\n")

  pdat <- data.frame(count = k, dif = as.vector(dif),
    model = factor(rep(colnames(dif), each = length(k)), levels = colnames(dif)))
  pic <- ggplot2::ggplot(pdat, ggplot2::aes(count, dif, color = model)) +
    ggplot2::geom_hline(yintercept = 0) + ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::scale_x_continuous(breaks = k) +
    ggplot2::labs(x = "Count", y = "Observed - Predicted", color = NULL,
      caption = "Positive deviations show underpredictions") +
    ggplot2::theme_minimal()

  invisible(list(models = mods, compare = compare, fit = fit, tests = tests,
    probs = data.frame(Count = k, Actual = obs, pred), pic = pic))
}
