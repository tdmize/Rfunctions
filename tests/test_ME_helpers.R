# Run from the repository root with Rscript tests/test_ME_helpers.R.
library(marginaleffects) # Marginal effects and hypotheses
library(nnet)            # Multinomial logit
library(MASS)            # Ordinal and negative binomial models
library(sandwich)        # Robust covariance matrices
source(Sys.getenv("ME_HELPER_FILE", "ME_helper_functions.R"))
set.seed(10012026)
options(width = 130, marginaleffects_parallel = FALSE)

results <- data.frame(test = character(), status = character(), detail = character())
check <- function(name, expr) {
  filter <- Sys.getenv("ME_TEST_FILTER", "")
  if (nzchar(filter) && !grepl(filter, name)) return(invisible(NULL))
  ans <- tryCatch({ force(expr); "PASS" }, error = function(e) conditionMessage(e))
  ok <- identical(ans, "PASS")
  results[nrow(results) + 1L, ] <<- list(name, if (ok) "PASS" else "FAIL", if (ok) "" else ans)
  cat(if (ok) "PASS" else "FAIL", name, if (ok) "" else paste(":", ans), "\n")
  flush.console()
}
near <- function(x, y, tol = 2e-6) {
  if (length(x) != length(y) || any(!is.finite(x)) || max(abs(x - y)) > tol)
    stop(sprintf("Got %s; expected %s", paste(signif(x, 8), collapse = ", "), paste(signif(y, 8), collapse = ", ")))
}
must_error <- function(expr, pattern = NULL) {
  ans <- tryCatch({ force(expr); NULL }, error = identity)
  if (is.null(ans)) stop("Expected rejection, but the helper returned a result.")
  if (!is.null(pattern) && !grepl(pattern, conditionMessage(ans), fixed = TRUE))
    stop("Unexpected error: ", conditionMessage(ans))
}

# Independent oracle: model matrices and link functions, with central-difference
# derivatives of the final statistic. Does not call either helper or comparisons().
parameters <- function(m) {
  if (inherits(m, "multinom")) return(as.vector(t(coef(m))))
  if (inherits(m, "polr")) return(c(coef(m), m$zeta))
  if (inherits(m, "svyolr")) return(c(m$coefficients, m$zeta))
  if (inherits(m, "clm")) return(c(m$beta, m$alpha))
  coef(m)
}
predict_manual <- function(m, d, beta) {
  X <- model.matrix(delete.response(terms(m)), d, contrasts.arg = m$contrasts, xlev = m$xlevels)
  if (inherits(m, "multinom")) {
    B <- matrix(beta, nrow = length(m$lev) - 1L, byrow = TRUE)
    eta <- cbind(0, X %*% t(B))
    eta <- eta - apply(eta, 1, max)
    p <- exp(eta); p <- p / rowSums(p); colnames(p) <- m$lev
    return(p)
  }
  if (inherits(m, c("polr", "svyolr", "clm"))) {
    coefficients <- if (inherits(m, "clm")) m$beta else m$coefficients
    X <- X[, names(coefficients), drop = FALSE]
    eta <- as.vector(X %*% beta[seq_along(coefficients)])
    zeta <- beta[-seq_along(coefficients)]
    link <- if (inherits(m, "clm")) m$link else m$method
    F <- if (link == "probit") pnorm else plogis
    C <- vapply(zeta, function(z) F(z - eta), numeric(nrow(d)))
    p <- cbind(C, 1) - cbind(0, C)
    colnames(p) <- if (inherits(m, "clm")) m$y.levels else m$lev
    return(p)
  }
  eta <- as.vector(X %*% beta)
  off <- model.offset(model.frame(delete.response(terms(m)), d, xlev = m$xlevels))
  if (is.null(off) && !is.null(m$call$offset)) off <- eval(m$call$offset, d, environment(formula(m)))
  if (!is.null(off)) eta <- eta + off
  if (inherits(m, "glm")) eta <- m$family$linkinv(eta)
  matrix(eta, ncol = 1L, dimnames = list(NULL, ""))
}
manual_summary <- function(m, d, var, total = FALSE, weighted = TRUE,
  shares_data = model.frame(m), shares_weights = NULL, avg_weights = NULL,
  by = NULL, beta = parameters(m)) {
  x <- shares_data[[var]]
  if (is.factor(x) || is.character(x) || is.logical(x)) {
    lev <- if (is.factor(x)) levels(droplevels(x)) else sort(unique(x))
    pairs <- combn(lev, 2, simplify = FALSE)
    p <- if (is.null(shares_weights)) prop.table(table(factor(x, levels = lev))) else
      tapply(shares_weights, factor(x, levels = lev), sum) / sum(shares_weights)
    pw <- if (weighted) vapply(pairs, function(a) sum(p[a]) / (length(lev) - 1), numeric(1)) else
      rep(1 / length(pairs), length(pairs))
  } else {
    pairs <- list(mean(d[[var]]) + c(-0.5, 0.5) * sd(d[[var]])); pw <- 1
  }
  rows <- if (is.null(by)) list(all = seq_len(nrow(d))) else split(seq_len(nrow(d)), d[[by]])
  out <- c()
  for (stratum in names(rows)) {
    ix <- rows[[stratum]]; dd <- d[ix, , drop = FALSE]
    effects <- lapply(pairs, function(pair) {
      lo <- hi <- dd
      if (is.factor(dd[[var]])) {
        lo[[var]] <- factor(pair[1], levels = levels(dd[[var]]), ordered = is.ordered(dd[[var]]))
        hi[[var]] <- factor(pair[2], levels = levels(dd[[var]]), ordered = is.ordered(dd[[var]]))
      } else { lo[[var]] <- pair[1]; hi[[var]] <- pair[2] }
      delta <- predict_manual(m, hi, beta) - predict_manual(m, lo, beta)
      if (is.null(avg_weights)) colMeans(delta) else apply(delta, 2, weighted.mean, w = avg_weights[ix])
    })
    E <- do.call(cbind, effects)
    val <- if (total) sum(abs(E) %*% pw) / 2 else as.vector(abs(E) %*% pw)
    names(val) <- if (total) stratum else paste(rownames(E), stratum, sep = "|")
    out <- c(out, val)
  }
  out
}
delta_manual <- function(fun, b, V) {
  J <- vapply(seq_along(b), function(j) {
    h <- 1e-5 * max(1, abs(b[j])); lo <- hi <- b; lo[j] <- lo[j] - h; hi[j] <- hi[j] + h
    (fun(hi) - fun(lo)) / (2 * h)
  }, numeric(length(fun(b))))
  J <- matrix(J, nrow = length(fun(b)))
  sqrt(pmax(0, diag(J %*% V %*% t(J))))
}
validate <- function(m, d, vars = "f", total = FALSE, weighted = TRUE,
  model_weights = TRUE, shares_data = model.frame(m), shares_weights = NULL,
  avg_weights = NULL, by = NULL, V = vcov(m)) {
  vv <- lapply(vars, function(v) if (is.numeric(d[[v]])) "sd" else "pairwise"); names(vv) <- vars
  hh <- if (total) totalme(m, weighted, model_weights) else meineq(m, weighted, model_weights)
  args <- list(model = m, variables = vv, newdata = d, hypothesis = hh, vcov = V, numderiv = "fdcenter")
  if (inherits(m, "svyolr")) args$vcov <- me_vcov(m)
  if (!is.null(by)) args$by <- unique(c(if (inherits(m, c("multinom", "polr", "clm", "svyolr"))) "group", "term", "contrast", by))
  if (!is.null(avg_weights)) args$wts <- avg_weights
  got <- do.call(avg_comparisons, args)
  exp <- se <- c()
  for (v in sort(vars)) {
    fun <- function(b) manual_summary(m, d, v, total, weighted, shares_data,
      if (model_weights) shares_weights else NULL, avg_weights, by, b)
    b <- parameters(m)
    covariance <- if (inherits(m, c("clm", "svyolr"))) V[names(b), names(b), drop = FALSE] else V
    exp <- c(exp, fun(b)); se <- c(se, delta_manual(fun, b, covariance))
  }
  # Package sorts outcome groups before subgroup columns. Sort both explicitly.
  got_key <- paste(got$term, if (is.null(got$group)) "" else as.character(got$group),
    if (is.null(by)) "all" else as.character(got[[by]]), sep = "|")
  expected_key <- c()
  for (v in sort(vars)) {
    groups <- if (total) "" else colnames(predict_manual(m, d, parameters(m)))
    strata <- if (is.null(by)) "all" else sort(unique(as.character(d[[by]])))
    expected_key <- c(expected_key, unlist(lapply(strata, function(s) paste(v, groups, s, sep = "|"))))
  }
  idx <- match(got_key, expected_key)
  if (anyNA(idx) || nrow(got) != length(expected_key)) stop("Result keys or row counts differ from the oracle.")
  near(got$estimate, exp[idx]); near(got$std.error, se[idx], tol = 8e-6)
  invisible(got)
}

n <- 600L
d <- data.frame(f = factor(sample(c("A", "B", "C", "D"), n, TRUE, c(.55, .25, .15, .05))),
  b = factor(sample(c("No", "Yes"), n, TRUE)), z = rnorm(n),
  s = factor(rep(c("First", "Second"), each = n / 2)), id = seq_len(n))
d$q <- ordered(sample(c("Low", "Medium", "High"), n, TRUE), levels = c("Low", "Medium", "High"))
d$w <- c(A = 1, B = 2, C = 4, D = 9)[d$f]
d$exposure <- exp(d$z / 4)
eta <- c(A = 0, B = .8, C = -.4, D = 1.5)[d$f] + .45 * d$z + .4 * (d$b == "Yes")
d$y <- eta + rnorm(n)
d$yb <- rbinom(n, 1, plogis(eta - 1))
d$yc <- rpois(n, exp(eta / 2))
d$ynb <- rnbinom(n, mu = exp(eta / 2), size = 2)
P <- cbind(0, .5 + .6 * eta, -.5 - .4 * eta, -.2 + .2 * eta)
P <- exp(P); P <- P / rowSums(P)
d$ym <- factor(vapply(seq_len(n), function(i) sample(c("Poor", "Fair", "Good", "Excellent"), 1, prob = P[i, ]), character(1)),
  levels = c("Poor", "Fair", "Good", "Excellent"))
d$yo <- ordered(cut(eta + rlogis(n), breaks = c(-Inf, -.5, .5, 1.5, Inf), labels = c("Poor", "Fair", "Good", "Excellent")))

lm1 <- lm(y ~ f + b + z + q, d)
check("me_vcov retains the native covariance for other model classes", {
  stopifnot(identical(me_vcov(lm1), vcov(lm1)))
})
check("lm weighted inequality and independent delta SE", validate(lm1, d))
check("lm unweighted inequality and independent delta SE", validate(lm1, d, weighted = FALSE))
check("multiple predictors including ordered and binary", validate(lm1, d, c("f", "b", "q")))
check("HC3 inequality standard error", validate(lm1, d, V = vcovHC(lm1, type = "HC3")))
check("clustered inequality standard error", validate(lm1, d, V = vcovCL(lm1, cluster = rep(1:60, each = 10))))
lmby <- lm(y ~ f * s + b + z, d)
check("by subgroup retains labels and effects", validate(lmby, d, by = "s"))
check("reverse pairwise and row order", {
  a <- avg_comparisons(lm1, variables = list(f = "pairwise"), vcov = FALSE)
  h <- meineq(lm1); near(h(a)$estimate, h(a[sample(nrow(a)), ])$estimate)
  b <- avg_comparisons(lm1, variables = list(f = "revpairwise"), hypothesis = meineq(lm1), vcov = FALSE)
  near(b$estimate, h(a)$estimate)
})
check("named exported weights and variable arguments", {
  a <- meineq_weights(lm1, f); v <- "f"
  near(a, meineq_weights(lm1, "f")); near(a, meineq_weights(lm1, v))
  p <- prop.table(table(d$f)); pair <- combn(names(p), 2)
  e <- setNames((p[pair[1, ]] + p[pair[2, ]]) / 3, paste(pair[2, ], "-", pair[1, ]))
  near(a, e[names(a)]); near(sum(a), 1)
})
wlm <- lm(y ~ f + z, d, weights = w)
check("weighted model prevalence shares", validate(wlm, d, shares_weights = d$w))
check("ignore model weights for prevalence", validate(wlm, d, model_weights = FALSE, shares_weights = d$w))
check("weighted averaging and weighted prevalence", validate(wlm, d, shares_weights = d$w, avg_weights = d$w))
dn <- d; dn$z[dn$f == "A" & dn$id %% 3 == 0] <- NA
mn <- lm(y ~ f + z, dn, subset = id > 80, na.action = na.exclude)
dn_used <- dn[rownames(model.frame(mn)), , drop = FALSE]
check("missingness subset and na.exclude use estimation sample", validate(mn, dn_used))
for (link in c("logit", "probit", "cloglog")) {
  gm <- glm(yb ~ f + b + z, d, family = binomial(link))
  check(paste("binary", link, "inequality and SE"), validate(gm, d, c("f", "b")))
}
gm <- glm(yc ~ f + z, d, family = poisson())
check("Poisson inequality and SE", validate(gm, d))
gmo <- glm(yc ~ f + z + offset(log(exposure)), d, family = poisson())
check("Poisson formula offset inequality and SE", validate(gmo, d))
nb <- glm.nb(ynb ~ f + z, d)
check("negative binomial inequality and SE", validate(nb, d))
mm <- multinom(ym ~ f + b + z + q + s, d, trace = FALSE, Hess = TRUE)
check("multinomial per-outcome inequalities and SE", validate(mm, d))
check("multinomial total ME inequalities and SE", validate(mm, d, c("f", "q"), total = TRUE))
check("multinomial unweighted total ME inequalities and SE", validate(mm, d, total = TRUE, weighted = FALSE))
wmm <- multinom(ym ~ f + b + z, d, weights = w, trace = FALSE, Hess = TRUE)
check("weighted multinomial total inequality and SE", validate(wmm, d, total = TRUE, shares_weights = d$w))
dmn <- dn; dmn$ym[dmn$id %% 9 == 0] <- NA
mmna <- multinom(ym ~ f + z, dmn, trace = FALSE, Hess = TRUE, model = TRUE, subset = id > 80)
mmused <- dmn[rownames(model.frame(mmna)), , drop = FALSE]
check("multinomial missingness subset prevalence", validate(mmna, mmused, total = TRUE))
check("multinomial continuous binary and nominal totals together", validate(mm, d, c("z", "b", "f"), total = TRUE))
check("multinomial totals by subgroup", validate(mm, d, c("f", "b", "z"), total = TRUE, by = "s"))
check("total identity: half sum of outcome inequalities", {
  i <- avg_comparisons(mm, variables = list(f = "pairwise"), hypothesis = meineq(mm), vcov = FALSE)
  t <- avg_comparisons(mm, variables = list(f = "pairwise"), hypothesis = totalme(mm), vcov = FALSE)
  near(t$estimate, sum(i$estimate) / 2)
})
for (link in c("logistic", "probit")) {
  om <- polr(yo ~ f + b + z, d, method = link, Hess = TRUE)
  check(paste("ordered", link, "inequality and SE"), validate(om, d))
  check(paste("ordered", link, "total and SE"), validate(om, d, c("f", "b", "z"), total = TRUE))
}
if (requireNamespace("ordinal", quietly = TRUE)) {
  for (link in c("logit", "probit")) {
    cm <- ordinal::clm(yo ~ f + b + z, data = d, link = link, model = TRUE)
    check(paste("clm", link, "inequality and mixed-predictor total SE"), {
      validate(cm, d)
      validate(cm, d, c("f", "b", "z"), total = TRUE, by = "s")
    })
  }
}
check("reject reference contrasts for multilevel nominal inequality", {
  must_error(avg_comparisons(lm1, variables = list(f = "reference"), hypothesis = meineq(lm1)))
})
check("reject scalar-outcome total", {
  must_error(avg_comparisons(lm1, variables = list(f = "pairwise"), hypothesis = totalme(lm1)))
})

# Safety checks: these should reject incomplete or non-probability input.
check("reject total after dropping outcome categories", {
  raw <- avg_comparisons(mm, variables = list(f = "pairwise"), vcov = FALSE)
  partial <- raw[raw$group == levels(raw$group)[1], ]
  must_error(totalme(mm)(partial))
})
check("selected categorical contrast total is calculated correctly", {
  raw <- avg_comparisons(mm, variables = list(f = c("A", "B")), vcov = FALSE)
  near(totalme(mm)(raw)$estimate, sum(abs(raw$estimate)) / 2)
})
check("reject duplicated pair replacing omitted pair", {
  raw <- avg_comparisons(lm1, variables = list(f = "pairwise"), vcov = FALSE)
  raw <- raw[c(1, 1, 3:6), ]
  must_error(meineq(lm1, weighted = FALSE)(raw))
})
check("reject reversed duplicate of a missing category pair", {
  raw <- avg_comparisons(lm1, variables = list(f = "pairwise"), vcov = FALSE)
  raw <- raw[c(1, 1, 3:6), ]
  pieces <- strsplit(raw$contrast[2], " - ", fixed = TRUE)[[1]]
  raw$contrast[2] <- paste(rev(pieces), collapse = " - ")
  raw$estimate[2] <- -raw$estimate[2]
  must_error(meineq(lm1, weighted = FALSE)(raw), "each category pair")
})
check("reject duplicated outcome replacing an omitted outcome", {
  raw <- avg_comparisons(mm, variables = list(b = "pairwise"), vcov = FALSE)
  raw <- raw[c(1, 1, 3:4), ]
  must_error(totalme(mm)(raw), "exactly one contrast for every outcome")
})
check("reject incomplete outcomes in only one subgroup", {
  raw <- avg_comparisons(mm, variables = list(f = "pairwise"),
    by = c("group", "term", "contrast", "s"), vcov = FALSE)
  must_error(totalme(mm)(raw[-1, ]), "exactly one contrast for every outcome")
})
check("reject non-probability differences and nonfinite input", {
  raw <- avg_comparisons(mm, variables = list(b = "pairwise"), vcov = FALSE)
  raw$estimate[1] <- raw$estimate[1] + .1
  must_error(totalme(mm)(raw), "sum to zero")
  raw$estimate[1] <- NA_real_
  must_error(totalme(mm)(raw), "finite marginal-effect contrasts")
})
check("reject zero-sum effects outside probability-difference bounds", {
  raw <- avg_comparisons(mm, variables = list(b = "pairwise"), vcov = FALSE)
  raw$estimate <- c(-2, 0, 0, 2)
  must_error(totalme(mm)(raw), "between -1 and 1")
  raw$estimate <- c(-.9, -.9, .9, .9)
  must_error(totalme(mm)(raw), "one unit of total probability")
})
check("reject incomplete category pairs in only one subgroup", {
  raw <- avg_comparisons(lmby, variables = list(f = "pairwise"),
    by = c("term", "contrast", "s"), vcov = FALSE)
  must_error(meineq(lmby, weighted = FALSE)(raw[-1, ]), "each category pair")
})
check("explicit subgroup names tolerate shuffled table columns", {
  if (packageVersion("marginaleffects") >= "1.0.0") {
    raw <- avg_comparisons(lmby, variables = list(f = "pairwise"),
      by = c("term", "contrast", "s"), vcov = FALSE)
    h <- meineq(lmby)
    shuffled <- raw[, c("estimate", "s", "contrast", "term")]
    near(h(shuffled, by = "s")$estimate, h(raw, by = "s")$estimate)
  } else stopifnot(identical(names(formals(meineq(lmby))), "x"))
})
check("logical binary predictor inequality", {
  dl <- d; dl$bl <- dl$b == "Yes"
  ml <- lm(y ~ f + bl + z, dl)
  out <- avg_comparisons(ml, variables = list(bl = "pairwise"), hypothesis = meineq(ml), vcov = FALSE)
  near(out$estimate, abs(coef(ml)["blTRUE"]))
})
check("character categorical predictor inequality", {
  dc <- d; dc$f <- as.character(dc$f); mc <- lm(y ~ f + z, dc)
  validate(mc, dc)
})
check("factor in formula from numeric predictor", {
  df <- d; df$fn <- as.integer(df$f); mf <- lm(y ~ factor(fn) + z, df)
  out <- avg_comparisons(mf, variables = list(fn = "pairwise"), hypothesis = meineq(mf), vcov = FALSE)
  stopifnot(nrow(out) == 1L, is.finite(out$estimate))
})

if (requireNamespace("survey", quietly = TRUE)) {
  des <- survey::svydesign(ids = ~1, weights = ~w, data = d)
  sm <- survey::svyglm(yb ~ f + z, des, family = quasibinomial())
  check("survey binomial inequality and SE", validate(sm, d, shares_weights = d$w))
  sdes <- subset(des, id > 80)
  smsub <- survey::svyglm(yb ~ f + z, sdes, family = quasibinomial())
  check("survey subset prevalence and SE", validate(smsub, d[d$id > 80, ], shares_weights = d$w[d$id > 80]))
  repdes <- survey::as.svrepdesign(des, type = "bootstrap", replicates = 30)
  srep <- survey::svyglm(yb ~ f + z, repdes, family = quasibinomial())
  check("replicate-weight survey prevalence and SE", validate(srep, d, shares_weights = d$w))
  desna <- survey::svydesign(ids = ~1, weights = ~w, data = dn)
  smna <- survey::svyglm(yb ~ f + z, desna, family = quasibinomial())
  use <- complete.cases(dn[c("yb", "f", "z")])
  check("survey missingness prevalence and SE", validate(smna, dn[use, ], shares_weights = dn$w[use]))
  for (link in c("logistic", "probit")) {
    so <- survey::svyolr(yo ~ f + b + z, des, method = link)
    check(paste("survey ordered", link, "aligned-covariance inequality and subgroup total SE"), {
      validate(so, d, shares_data = d, shares_weights = d$w)
      validate(so, d, c("f", "b", "z"), total = TRUE, by = "s",
        shares_data = d, shares_weights = d$w, avg_weights = d$w)
    })
  }
  # survey 4.2.1's replicate svyolr passes weights=w inside model.frame;
  # a data column called w shadows that internal replicate-weight argument.
  ordered_repdes <- repdes
  ordered_repdes$variables$w <- NULL
  ro <- survey::svyolr(yo ~ f + z, ordered_repdes)
  stopifnot(any(diag(vcov(ro)) > 0))
  check("replicate-weight survey ordered aligned-covariance total and SE", validate(ro, d, total = TRUE,
    shares_data = d, shares_weights = d$w, avg_weights = d$w))
  check("survey ordered ambiguous cutpoint updates are rejected and renamed labels work", {
    labeled <- d
    levels(labeled$yo) <- c("Intercept: A|B", "C", "A", "B|C")
    ld <- survey::svydesign(ids = ~1, weights = ~w, data = labeled)
    lo <- survey::svyolr(yo ~ f + z, ld)
    must_error(me_vcov(lo), "Cutpoint labels are ambiguous")
    levels(labeled$yo) <- paste0("Outcome", 1:4)
    ld <- survey::svydesign(ids = ~1, weights = ~w, data = labeled)
    lo <- survey::svyolr(yo ~ f + z, ld)
    validate(lo, labeled, total = TRUE, shares_data = labeled, shares_weights = labeled$w)
  })
}

if (requireNamespace("suest", quietly = TRUE)) {
  m2 <- lm(y ~ f + b + z, d)
  fit <- suest::suest(lm1, m2, model_names = c("Full", "Reduced"))
  check("suest same-sample inequality estimates and joint SE", {
    got <- avg_comparisons(fit, variables = list(f = "pairwise"), newdata = d, hypothesis = meineq(fit), numderiv = "fdcenter")
    expected <- vapply(fit$models, function(m) manual_summary(m, d, "f"), numeric(1))
    near(got$estimate, expected[as.character(got$group)])
    fun <- function(beta) vapply(seq_along(fit$models), function(i)
      manual_summary(fit$models[[i]], d, "f", beta = beta[fit$index[[i]]][seq_along(coef(fit$models[[i]]))]), numeric(1))
    se <- delta_manual(fun, coef(fit), vcov(fit)); names(se) <- fit$model_names
    near(got$std.error, se[as.character(got$group)], 8e-6)
    diff <- hypotheses(got, hypothesis = difference ~ revpairwise)
    J <- numDeriv::grad(function(b) diff(fun(b)), coef(fit))
    near(diff$std.error, sqrt(as.numeric(t(J) %*% vcov(fit) %*% J)), 8e-6)
  })
  pre <- glm(yb ~ f + z, d, family = binomial(), subset = id <= 300)
  post <- glm(yb ~ f + z, d, family = binomial(), subset = id > 300)
  sf <- suest::suest(pre, post, model_names = c("Pre", "Post"))
  check("suest disjoint samples own category shares and SE", {
    nd <- suest::suest_newdata(sf)
    got <- avg_comparisons(sf, variables = list(f = "pairwise"), newdata = nd, hypothesis = meineq(sf), numderiv = "fdcenter")
    fun <- function(beta) vapply(seq_along(sf$models), function(i) {
      dd <- d[if (i == 1) d$id <= 300 else d$id > 300, ]
      manual_summary(sf$models[[i]], dd, "f", beta = beta[sf$index[[i]]])
    }, numeric(1))
    near(got$estimate, fun(coef(sf))); near(got$std.error, delta_manual(fun, coef(sf), vcov(sf)), 8e-6)
  })
  mm2 <- multinom(ym ~ f + b + z, d, trace = FALSE, Hess = TRUE)
  smm <- suest::suest(mm, mm2, model_names = c("Full", "Reduced"))
  check("suest multinomial total inequalities and joint SE", {
    got <- avg_comparisons(smm, variables = list(f = "pairwise"), newdata = d, hypothesis = totalme(smm), numderiv = "fdcenter")
    fun <- function(beta) vapply(seq_along(smm$models), function(i)
      manual_summary(smm$models[[i]], d, "f", total = TRUE, beta = beta[smm$index[[i]]]), numeric(1))
    near(got$estimate, fun(coef(smm))); near(got$std.error, delta_manual(fun, coef(smm), vcov(smm)), 8e-6)
  })
  mma <- multinom(ym ~ f + z, d, subset = id <= 300, trace = FALSE, Hess = TRUE)
  mmb <- multinom(ym ~ f + z, d, subset = id > 300, trace = FALSE, Hess = TRUE)
  sdis <- suest::suest(mma, mmb, model_names = c("Pre", "Post"))
  check("suest disjoint multinomial totals and own shares", {
    got <- avg_comparisons(sdis, variables = list(f = "pairwise"), newdata = suest::suest_newdata(sdis),
      hypothesis = totalme(sdis), numderiv = "fdcenter")
    fun <- function(beta) vapply(seq_along(sdis$models), function(i) {
      dd <- d[if (i == 1) d$id <= 300 else d$id > 300, ]
      manual_summary(sdis$models[[i]], dd, "f", total = TRUE, beta = beta[sdis$index[[i]]])
    }, numeric(1))
    near(got$estimate, fun(coef(sdis))); near(got$std.error, delta_manual(fun, coef(sdis), vcov(sdis)), 8e-6)
  })
  check("reject scalar model total in mixed suest system", {
    mixed <- suest::suest(pre, mma, model_names = c("Binary", "Multinomial"))
    must_error(avg_comparisons(mixed, variables = list(f = "pairwise"), newdata = d,
      hypothesis = totalme(mixed), vcov = FALSE))
  })
  check("suest model names containing double colon", {
    sn <- suest::suest(lm1, m2, model_names = c("Study", "Study::Reduced"))
    out <- avg_comparisons(sn, variables = list(f = "pairwise"), newdata = d, hypothesis = meineq(sn), vcov = FALSE)
    stopifnot(nrow(out) == 2L)
    near(out$estimate, vapply(sn$models, function(m) manual_summary(m, d, "f"), numeric(1)))
  })
  check("suest categorical names containing double colon and total difference SE", {
    sn <- suest::suest(mma, mmb, model_names = c("Study", "Study::Reduced"))
    out <- avg_comparisons(sn, variables = list(f = "pairwise"), newdata = suest::suest_newdata(sn),
      hypothesis = totalme(sn), numderiv = "fdcenter")
    fun <- function(beta) vapply(seq_along(sn$models), function(i) {
      dd <- d[if (i == 1) d$id <= 300 else d$id > 300, ]
      manual_summary(sn$models[[i]], dd, "f", total = TRUE, beta = beta[sn$index[[i]]])
    }, numeric(1))
    near(out$estimate, fun(coef(sn)))
    diff <- hypotheses(out, hypothesis = difference ~ revpairwise)
    J <- numDeriv::grad(function(b) diff(fun(b)), coef(sn))
    near(diff$std.error, sqrt(as.numeric(t(J) %*% vcov(sn) %*% J)), 8e-6)
  })
  check("suest nested model names and outcome delimiters use exact category shares", {
    da <- d[d$id <= 300, ]; levels(da$ym) <- c("b::first", "b::second", "other", "last")
    ml <- multinom(ym ~ f + z, da, weights = w, trace = FALSE, Hess = TRUE)
    sn <- suest::suest(ml, mmb, model_names = c("A", "A::b"), weight_type = "pweight")
    for (h in list(meineq, totalme)) {
      got <- avg_comparisons(sn, variables = list(f = "pairwise"), newdata = suest::suest_newdata(sn),
        hypothesis = h(sn), vcov = FALSE)
      for (i in seq_along(sn$models)) {
        dd <- if (i == 1L) da else d[d$id > 300, ]
        expected <- manual_summary(sn$models[[i]], dd, "f", total = identical(h, totalme),
          shares_weights = if (i == 1L) dd$w else NULL)
        group <- if (identical(h, totalme)) sn$model_names[i] else
          paste0(sn$model_names[i], "::", sub("\\|all$", "", names(expected)))
        idx <- match(group, as.character(got$group))
        stopifnot(!anyNA(idx))
        near(got$estimate[idx], expected)
      }
    }
  })
  check("weighted overlapping multinomial and ordered subgroup summaries and joint SE", {
    da <- d[d$id <= 400, ]; db <- d[d$id > 200, ]
    ma <- multinom(ym ~ f + z + s, da, weights = w, trace = FALSE, Hess = TRUE)
    mb <- polr(yo ~ f + z + s, db, weights = w, Hess = TRUE)
    sn <- suest::suest(ma, mb, model_names = c("Nominal", "Ordered"), weight_type = "pweight")
    for (tot in c(FALSE, TRUE)) {
      fun <- function(beta) unlist(lapply(seq_along(sn$models), function(i) {
        dd <- if (i == 1L) da else db
        manual_summary(sn$models[[i]], dd, "f", total = tot,
          shares_weights = dd$w, avg_weights = dd$w, by = "s", beta = beta[sn$index[[i]]])
      }), use.names = FALSE)
      expected_keys <- unlist(lapply(seq_along(sn$models), function(i) {
        strata <- c("First", "Second")
        groups <- if (tot) sn$model_names[i] else paste0(sn$model_names[i], "::", sn$category_levels[[i]])
        unlist(lapply(strata, function(s) paste(groups, s, sep = "|")))
      }), use.names = FALSE)
      got <- avg_comparisons(sn, variables = list(f = "pairwise"), newdata = suest::suest_newdata(sn),
        by = c("group", "term", "contrast", "s"), wts = ".suest_weight",
        hypothesis = if (tot) totalme(sn) else meineq(sn), numderiv = "fdcenter")
      idx <- match(paste(got$group, got$s, sep = "|"), expected_keys)
      stopifnot(!anyNA(idx), nrow(got) == length(expected_keys))
      near(got$estimate, fun(coef(sn))[idx])
      near(got$std.error, delta_manual(fun, coef(sn), vcov(sn))[idx], 8e-6)
    }
  })
}

cat("\nAudit summary:", sum(results$status == "PASS"), "passed;", sum(results$status == "FAIL"), "failed.\n")
output <- Sys.getenv("ME_TEST_OUTPUT", ".")
dir.create(output, showWarnings = FALSE, recursive = TRUE)
write.csv(results, file.path(output, "test_results.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output, "sessionInfo.txt"))
print(results[results$status == "FAIL", ], row.names = FALSE)
if (any(results$status == "FAIL")) stop("ME helper validation failed.", call. = FALSE)
