############################################################
# Count model helper functions
#
# countfit(): compares Poisson (PRM), negative binomial (NBRM),
# zero-inflated Poisson (ZIP), and zero-inflated negative binomial
# (ZINB) models, with the same output as Stata's countfit. By default
# the zero-inflated models use all predictors in both the count and
# inflation parts; use inflate = ~ x1 + x2 to choose the inflation
# predictors, models = c("PRM", "NBRM", "ZIP") to compare only some
# of the models, and forcevuong = TRUE to add Vuong tests. As in Stata,
# noestimates, nodifferences, noprtable, and nofit = TRUE skip sections.
#
#   source("https://raw.githubusercontent.com/tdmize/Rfunctions/main/Count_helper_functions.R")
#   cf <- countfit(y ~ x1 + x2, data = mydata)
#   cf$pic   # observed - predicted graph (a ggplot object)
#
# Required packages: MASS, pscl, ggplot2
############################################################

countfit <- function(formula, data, inflate = NULL, maxcount = 9,
                     models = c("PRM", "NBRM", "ZIP", "ZINB"), forcevuong = FALSE,
                     noestimates = FALSE, nodifferences = FALSE, noprtable = FALSE, nofit = FALSE) {
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
  opts <- list(forcevuong = forcevuong, noestimates = noestimates, nodifferences = nodifferences,
    noprtable = noprtable, nofit = nofit)
  for (o in names(opts)) if (!isTRUE(opts[[o]]) && !isFALSE(opts[[o]]))
    stop(o, " must be TRUE or FALSE", call. = FALSE)
  alias <- c(prm = "PRM", poisson = "PRM", nbrm = "NBRM", nbreg = "NBRM", negbin = "NBRM",
    zip = "ZIP", zinb = "ZINB")
  if (!is.character(models) || length(models) == 0 || anyNA(alias[tolower(models)]))
    stop("models can include \"PRM\", \"NBRM\", \"ZIP\", and \"ZINB\"", call. = FALSE)
  sel <- intersect(c("PRM", "NBRM", "ZIP", "ZINB"), alias[tolower(models)])

  # One estimation sample for all models: complete cases on every model variable
  data <- as.data.frame(data)
  zformula <- sformula <- formula
  if (is.null(inflate)) {  # inflation part: the count predictors, without any offset
    tt <- terms(formula, data = data)
    zrhs <- if (length(attr(tt, "term.labels"))) reformulate(attr(tt, "term.labels"),
      intercept = attr(tt, "intercept") == 1)[[2]] else 1
  } else {
    zrhs <- inflate[[2]]
    sformula[[3]] <- call("+", formula[[3]], inflate[[2]])
  }
  zformula[[3]] <- call("|", formula[[3]], zrhs)
  vlab <- lapply(data, function(v) attr(v, "label", exact = TRUE))  # subsetting drops some labels
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
  for (f in c(list(formula), if (any(c("ZIP", "ZINB") %in% sel)) list(as.formula(call("~", zrhs))))) {
    mf <- model.frame(f, data, drop.unused.levels = TRUE)
    M <- model.matrix(attr(mf, "terms"), mf)
    q <- qr(M)
    if (q$rank < ncol(M))
      stop("These predictors do not vary or are perfectly collinear with other predictors: ",
        paste(colnames(M)[q$pivot[(q$rank + 1):ncol(M)]], collapse = ", "),
        ". Remove them and run countfit() again.", call. = FALSE)
  }

  # Fit each model; stop with a clear message if any fail or do not converge
  zfit <- function(dist, EM = FALSE, start = NULL, reltol = 1e-12) suppressWarnings(pscl::zeroinfl(zformula, data = data,
    dist = dist, na.action = na.omit, control = pscl::zeroinfl.control(reltol = reltol, EM = EM, start = start)))
  zem <- function(dist) {  # EM to a looser tolerance (it is slow near the boundary), then polish
    s <- zfit(dist, EM = TRUE, reltol = 1e-6)
    zfit(dist, start = list(count = coef(s, "count"), zero = coef(s, "zero"), theta = s$theta))
  }
  nbstart <- function(th) tryCatch(suppressWarnings(if (is.null(th))
    MASS::glm.nb(formula, data = data, na.action = na.omit) else
    MASS::glm.nb(formula, data = data, na.action = na.omit, init.theta = th)), error = function(e) NULL)
  nbfit <- function() {
    m <- nbstart(NULL)
    if (is.null(m) || !is.null(m$th.warn) || m$theta > 1000) {  # glm.nb's default start can fail; try others
      fits <- Filter(Negate(is.null), c(list(m), lapply(c(1, 0.1, 0.01), nbstart)))
      if (!length(fits)) stop("estimation failed from every starting value")
      m <- fits[[which.max(sapply(fits, function(x) as.numeric(logLik(x))))]]
    }
    if (m$theta > 1000) {  # no overdispersion: the estimate is at the Poisson boundary (alpha = 0)
      p <- glm(formula, family = poisson, data = data, na.action = na.omit)
      if (as.numeric(logLik(p)) >= as.numeric(logLik(m)))
        m <- structure(c(p, list(theta = Inf, SE.theta = NA, twologlik = 2 * as.numeric(logLik(p)))),
          class = c("negbin", class(p)))
    }
    m
  }
  fitters <- list(
    PRM  = function() glm(formula, family = poisson, data = data, na.action = na.omit),
    NBRM = nbfit,
    ZIP  = function() zfit("poisson"),
    ZINB = function() zfit("negbin"))
  mods <- list()
  why <- character()
  for (m in sel) {
    fit <- tryCatch(fitters[[m]](), error = function(e) conditionMessage(e))
    if (is.character(fit)) why[m] <- fit
    else if (isFALSE(fit$converged) || !is.finite(logLik(fit)) || !all(is.finite(na.omit(coef(fit)))))
      why[m] <- "the estimation stopped before converging"
    else mods[[m]] <- fit
  }
  # A zero-inflated fit no better than the simpler model inside it may have stalled; refit from the
  # simpler model's estimates (inflation part off) and from EM starting values, and keep the best fit
  nested <- list(ZIP = "PRM", ZINB = c("NBRM", "ZIP"))
  for (m in intersect(c("ZIP", "ZINB"), names(mods))) {
    base <- sapply(intersect(nested[[m]], names(mods)), function(b) as.numeric(logLik(mods[[b]])))
    if (length(base) && as.numeric(logLik(mods[[m]])) - max(base) >= .01) next
    dist <- if (m == "ZIP") "poisson" else "negbin"
    kz <- ncol(model.matrix(mods[[m]], "zero"))
    off0 <- replace(rep(0, kz), colnames(model.matrix(mods[[m]], "zero")) == "(Intercept)", -10)
    starts <- list()
    if (m == "ZIP" && "PRM" %in% names(mods)) starts <- c(starts, list(list(count = coef(mods$PRM), zero = off0)))
    if (m == "ZINB" && "NBRM" %in% names(mods)) starts <- c(starts, list(list(count = coef(mods$NBRM), zero = off0,
      theta = min(mods$NBRM$theta, 1e4))))
    if (m == "ZINB" && "ZIP" %in% names(mods)) starts <- c(starts, list(list(count = coef(mods$ZIP, "count"),
      zero = coef(mods$ZIP, "zero"), theta = 100)))
    for (st in c(starts, list("EM"))) {
      alt <- tryCatch(if (identical(st, "EM")) zem(dist) else zfit(dist, start = st), error = function(e) NULL)
      if (!is.null(alt) && isTRUE(alt$converged) && as.numeric(logLik(alt)) > as.numeric(logLik(mods[[m]])) + 1e-6)
        mods[[m]] <- alt
    }
  }
  # A ZINB near the Poisson boundary is nearly a ZIP; if it fits better, the ZIP stopped short, so refit it from there
  if (all(c("ZIP", "ZINB") %in% names(mods)) && mods$ZINB$theta > 1000 &&
      as.numeric(logLik(mods$ZINB)) > as.numeric(logLik(mods$ZIP)) + 1e-6) {
    alt <- tryCatch(zfit("poisson", start = list(count = coef(mods$ZINB, "count"), zero = coef(mods$ZINB, "zero"))),
      error = function(e) NULL)
    if (!is.null(alt) && isTRUE(alt$converged) && as.numeric(logLik(alt)) > as.numeric(logLik(mods$ZIP)) + 1e-6) {
      mods$ZIP <- alt
      alt <- tryCatch(zfit("negbin", start = list(count = coef(alt, "count"), zero = coef(alt, "zero"), theta = 1e4)),
        error = function(e) NULL)
      if (!is.null(alt) && isTRUE(alt$converged) && as.numeric(logLik(alt)) > as.numeric(logLik(mods$ZINB)) + 1e-6)
        mods$ZINB <- alt
    }
  }
  # A model that fits worse than the simpler model nested in it did not reach its maximum
  ll <- sapply(mods, function(m) as.numeric(logLik(m)))
  for (pr in list(c("PRM", "NBRM"), c("ZIP", "ZINB"), c("PRM", "ZIP"), c("NBRM", "ZINB"))) {
    if (all(pr %in% names(mods)) && ll[pr[2]] < ll[pr[1]] - 1) {
      why[pr[2]] <- "fit worse than the simpler nested model"
      mods[[pr[2]]] <- NULL
    }
  }
  why <- why[intersect(sel, names(why))]
  if (length(why)) {
    ok <- setdiff(sel, names(why))
    stop("The ", paste(names(why), collapse = " and "), " model", if (length(why) > 1) "s",
      " did not converge (", paste0(names(why), ": ", why, collapse = "; "), ").",
      if (length(ok)) paste0(" Run countfit() again without ", if (length(why) > 1) "them" else "it",
        ":\n  models = c(", paste0("\"", ok, "\"", collapse = ", "), ")"),
      "\nIf several models fail, check for predictors that do not vary or that perfectly predict the outcome.",
      call. = FALSE)
  }

  # Estimates: exponentiated coefficients and z values, with standard errors from the
  # observed information matrix (as Stata reports them)
  est <- lapply(sel, function(m) .cf_estimates(mods[[m]], m, y))
  names(est) <- sel
  cm <- if (any(c("PRM", "NBRM") %in% sel)) mods[[intersect(c("PRM", "NBRM"), sel)[1]]] else NULL
  zm <- if (any(c("ZIP", "ZINB") %in% sel)) mods[[intersect(c("ZIP", "ZINB"), sel)[1]]]
  rows <- list(count = if (is.null(cm)) .cf_rows(model.matrix(zm, "count"), zm$terms$count, data, vlab) else
    .cf_rows(model.matrix(cm), terms(cm), data, vlab))
  if (!is.null(zm)) rows$zero <- .cf_rows(model.matrix(zm, "zero"), zm$terms$zero, data, vlab, TRUE)

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

  # Largest difference and where it occurs, with Stata's rules for ties
  compare <- data.frame(
    MaxDiff = apply(dif, 2, function(d) if (max(abs(d)) == abs(max(d))) max(d) else min(d)),
    AtValue = apply(dif, 2, function(d) k[max(which(abs(d) > max(abs(d)) - 1e-5))]),
    MeanAbsDiff = colMeans(abs(dif)), row.names = sel)

  # Fit statistics and pairwise comparisons; dif = first model - second model
  bic <- sapply(mods, BIC)
  aic <- sapply(mods, AIC)
  fit <- data.frame(BIC = bic, AIC = aic, row.names = sel)
  tests <- lrtests <- vuong <- NULL
  if (length(sel) > 1) {
    pairs <- combn(sel, 2)
    a <- pairs[1, ]
    b <- pairs[2, ]
    bd <- bic[a] - bic[b]
    tests <- data.frame(Comparison = paste(a, "vs", b),
      BICdif = round(bd, 3), BICprefer = ifelse(abs(bd) <= 1e-10, "", ifelse(bd < 0, a, b)),
      BICover = ifelse(bd < 0 | abs(bd) <= 1e-10, b, a),
      Evidence = ifelse(abs(bd) <= 1e-10, "no", as.character(cut(abs(bd), c(0, 2, 6, 10, Inf),
        labels = c("Weak", "Positive", "Strong", "Very strong")))),
      AICdif = round(aic[a] - aic[b], 3), AICprefer = ifelse(aic[a] < aic[b], a, b),
      AICover = ifelse(aic[a] < aic[b], b, a))
    nest <- paste(a, b) %in% c("PRM NBRM", "ZIP ZINB")
    if (any(nest)) {
      lr <- pmax(2 * (ll[b[nest]] - ll[a[nest]]), 0)
      p <- pchisq(lr, df = 1, lower.tail = FALSE) / 2
      lrtests <- data.frame(Comparison = paste(a[nest], "vs", b[nest]),
        LRchi2 = round(lr, 3), p = round(p, 3),
        Prefer = ifelse(p < .05, b[nest], a[nest]), Over = ifelse(p < .05, a[nest], b[nest]))
    }
    vpair <- paste(a, b) %in% c("PRM ZIP", "NBRM ZINB")
    if (forcevuong && any(vpair)) {
      v <- mapply(function(m1, m2) {
        d <- est[[m2]]$llobs - est[[m1]]$llobs
        sqrt(length(d)) * mean(d) / sd(d)
      }, a[vpair], b[vpair])
      vuong <- data.frame(Comparison = paste(a[vpair], "vs", b[vpair]),
        Vuong = round(v, 3), p = round(pnorm(-abs(v)), 3),
        Prefer = ifelse(v > 0, b[vpair], a[vpair]), Over = ifelse(v > 0, a[vpair], b[vpair]))
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

  out <- structure(list(models = mods, estimates = lapply(est, function(e) e[names(e) != "llobs"]),
    rows = rows, outcome = deparse(formula[[2]]), compare = compare, fit = fit, ll = ll,
    tests = tests, lrtests = lrtests, vuong = vuong, opts = opts,
    probs = data.frame(Count = k, Actual = obs, pred), pic = pic, N = length(y)),
    class = "countfit")
  print(out)
  invisible(out)
}

# Exponentiated coefficients, z values (observed-information standard errors), and
# each observation's log likelihood for one fitted model
.cf_estimates <- function(m, type, y) {
  off <- function(o) if (is.null(o)) 0 else o
  hess <- function(gr, p, h) {
    H <- sapply(seq_along(p), function(j) {
      e <- replace(numeric(length(p)), j, h[j])
      (gr(p + e) - gr(p - e)) / (2 * h[j])
    })
    (H + t(H)) / 2
  }
  # Hessian by differencing the analytic gradient: first with steps scaled to each
  # predictor, then with steps scaled to each standard error
  sefrom <- function(gr, p, scale) {
    se <- rep(NA_real_, length(p))
    h <- 1e-5 / scale
    for (pass in 1:2) {
      V <- tryCatch(solve(-hess(gr, p, h)), error = function(e) NULL)
      if (is.null(V)) break
      se <- suppressWarnings(sqrt(diag(V)))
      if (!all(is.finite(se) & se > 0)) break
      h <- 1e-5 * se
    }
    se
  }
  rms <- function(M) pmax(sqrt(colMeans(M^2)), 1e-8)
  nbll <- function(y, mu, t) t * (digamma(y + t) - digamma(t) + log(t) - log(t + mu) + 1 - (y + t) / (t + mu))
  if (type %in% c("PRM", "NBRM")) {
    b <- coef(m)
    keep <- !is.na(b)
    X <- model.matrix(m)[, keep, drop = FALSE]
    o <- off(model.offset(m$model))
    mu <- fitted(m)
    if (type == "PRM") {
      se <- sqrt(diag(vcov(m)))[names(b)[keep]]
      llobs <- dpois(y, mu, log = TRUE)
    } else if (is.infinite(m$theta)) {
      se <- sqrt(diag(vcov(m)))[names(b)[keep]]
      lnalpha <- c(b = 0, z = NA)
      llobs <- dpois(y, mu, log = TRUE)
    } else {
      kk <- ncol(X) + 1
      gr <- function(p) {
        t <- exp(p[kk])
        mu <- exp(drop(X %*% p[-kk]) + o)
        c(colSums(X * ((y - mu) * t / (t + mu))), sum(nbll(y, mu, t)))
      }
      se <- sefrom(gr, c(b[keep], log(m$theta)), c(rms(X), 1))
      lnalpha <- c(b = 1 / m$theta, z = -log(m$theta) / se[kk])
      se <- se[-kk]
      llobs <- dnbinom(y, size = m$theta, mu = mu, log = TRUE)
    }
    z <- setNames(rep(NA_real_, length(b)), names(b))
    z[keep] <- b[keep] / se
    etac <- log(mu)
    etaz <- NULL
    out <- list(count = list(b = exp(b), z = z))
    if (type == "NBRM") out$lnalpha <- lnalpha
  } else {
    X <- model.matrix(m, "count")
    Z <- model.matrix(m, "zero")
    oc <- off(m$offset$count)
    oz <- off(m$offset$zero)
    kx <- ncol(X)
    kz <- ncol(Z)
    nb <- type == "ZINB"
    pieces <- function(p) {
      mu <- exp(drop(X %*% p[1:kx]) + oc)
      pz <- plogis(drop(Z %*% p[kx + 1:kz]) + oz)
      t <- if (nb) exp(p[kx + kz + 1]) else NA
      p0 <- if (nb) (t / (t + mu))^t else exp(-mu)
      list(mu = mu, pz = pz, t = t, p0 = p0, f0 = pz + (1 - pz) * p0)
    }
    gr <- function(p) {
      q <- pieces(p)
      z0 <- y == 0
      dp0 <- if (nb) -q$p0 * q$t * q$mu / (q$t + q$mu) else -q$mu * q$p0
      deta <- ifelse(z0, (1 - q$pz) * dp0 / q$f0, if (nb) q$t * (y - q$mu) / (q$t + q$mu) else y - q$mu)
      dzeta <- ifelse(z0, (1 - q$p0) * q$pz * (1 - q$pz) / q$f0, -q$pz)
      g <- c(colSums(X * deta), colSums(Z * dzeta))
      if (nb) g <- c(g, sum(ifelse(z0,
        (1 - q$pz) * q$t * q$p0 * (log(q$t / (q$t + q$mu)) + q$mu / (q$t + q$mu)) / q$f0,
        nbll(y, q$mu, q$t))))
      g
    }
    p <- c(coef(m, "count"), coef(m, "zero"), if (nb) log(m$theta))
    if (nb && m$theta > 1e4) {  # no overdispersion: alpha at 0, so hold it fixed
      se <- c(sefrom(function(b) gr(c(b, p[kx + kz + 1]))[1:(kx + kz)], p[1:(kx + kz)], c(rms(X), rms(Z))), NA)
    } else {
      se <- sefrom(gr, p, c(rms(X), rms(Z), if (nb) 1))
    }
    q <- pieces(p)
    etac <- log(q$mu)
    etaz <- qlogis(q$pz)
    llobs <- ifelse(y == 0, log(q$f0), log(1 - q$pz) +
      (if (nb) dnbinom(y, size = q$t, mu = q$mu, log = TRUE) else dpois(y, q$mu, log = TRUE)))
    out <- list(count = list(b = exp(p[1:kx]), z = p[1:kx] / se[1:kx]),
      zero = list(b = exp(p[kx + 1:kz]), z = p[kx + 1:kz] / se[kx + 1:kz]))
    names(out$count$b) <- names(out$count$z) <- colnames(X)
    names(out$zero$b) <- names(out$zero$z) <- colnames(Z)
    if (nb) out$lnalpha <- c(b = 1 / m$theta, z = -log(m$theta) / se[kx + kz + 1])
  }
  # Degenerate fits: predictions of essentially exactly 0 or 1 (perfect prediction), or huge
  # coefficients with tiny z values
  z <- c(out$count$z[!is.na(out$count$b)], out$zero$z)
  far <- function(e, cut) mean(e > cut) >= .01 || any(e > 25)  # a group of cases, not a few extreme ones
  extreme <- far(-etac, 14) || (!is.null(etaz) && far(abs(etaz), 15)) ||
    any(abs(log(c(out$count$b[!is.na(out$count$b)], out$zero$b))) > 10 & (is.na(z) | abs(z) < 1))
  if (extreme) {
    warning("countfit(): the ", type, " model has extreme estimates: some of its predictions are essentially ",
      "exactly 0 or 1. This usually means a predictor perfectly predicts the outcome or always-zero status, ",
      "or that the model's zero-inflation part is not needed; interpret that model with caution.", call. = FALSE)
  } else if (anyNA(z)) {
    warning("countfit(): some z values could not be computed for the ", type, " model and are shown as .", call. = FALSE)
  }
  out$llobs <- llobs
  out
}

# Row labels for the estimates table: variable labels where available, factor
# levels under their variable, the constant last
.cf_rows <- function(mm, tt, data, vlab, inflate = FALSE) {
  asg <- attr(mm, "assign")
  tl <- attr(tt, "term.labels")
  lab <- function(v) {
    l <- vlab[[gsub("^`|`$", "", v)]]
    if (is.character(l) && length(l) == 1 && nzchar(l)) l else v
  }
  out <- data.frame(label = character(), coef = character())
  add <- function(label, coef = NA_character_) out[nrow(out) + 1, ] <<- list(label, coef)
  for (t in setdiff(unique(asg), 0)) {
    cols <- colnames(mm)[asg == t]
    term <- tl[t]
    v <- data[[gsub("^`|`$", "", term)]]
    isfac <- is.factor(v) || is.character(v) || is.logical(v)
    if ((isfac || length(cols) > 1 || grepl(":", term)) && (nrow(out) > 0 || inflate)) add("")  # as Stata
    if (isfac) {
      add(lab(term))
      for (cc in cols) add(paste0(substring(cc, nchar(term) + 1), " "), cc)
    } else if (length(cols) == 1 && cols == term) {
      add(lab(term), cols)
    } else {
      for (cc in cols) add(cc, cc)
    }
  }
  if (0 %in% asg) add("Constant", colnames(mm)[asg == 0])
  out
}

print.countfit <- function(x, ...) {
  sel <- names(x$models)
  o <- x$opts
  fw <- function(v, d, w) {  # Stata's %w.df display, switching to e-format when too wide
    s <- ifelse(!is.finite(v), ".", sprintf(paste0("%.", d, "f"), v))
    long <- is.finite(v) & nchar(s) > w
    s[long] <- sprintf(if (w >= 9) "%.2e" else "%.1e", v[long])
    s
  }
  lab30 <- function(l) {
    sp <- grepl(" $", l)
    l <- sub(" $", "", l)
    l <- ifelse(nchar(l) > 30 - sp, paste0(substr(l, 1, 27 - sp), "..."), l)
    sprintf("%30s |", paste0(l, ifelse(sp, " ", "")))
  }
  cells <- function(s) paste0(sprintf("%10s  ", s), collapse = "")

  # Table of estimates (Stata: estimates table, eform b(%9.3f) t(%7.2f) label varwidth(30))
  if (!o$noestimates) {
    w <- 32 + 12 * length(sel)
    rule <- paste0(strrep("-", 31), "+", strrep("-", w - 32), "\n")
    cat("\n", strrep("-", w), "\n", sprintf("%30s |", "Variable"),
      paste0(sprintf("%7s     ", sel), collapse = ""), "\n", sep = "")
    eqtable <- function(name, rows, eq) {
      cat(rule, sprintf("%-31s|", name), "\n", sep = "")
      for (i in seq_len(nrow(rows))) {
        cf <- rows$coef[i]
        if (is.na(cf)) {
          cat(lab30(rows$label[i]), "\n", sep = "")
          next
        }
        get <- function(what, d, w) sapply(sel, function(m) {
          e <- x$estimates[[m]][[eq]]
          if (is.null(e) || !cf %in% names(e[[what]])) "" else fw(e[[what]][[cf]], d, w)
        })
        cat(lab30(rows$label[i]), cells(get("b", 3, 9)), "\n",
          sprintf("%30s |", ""), cells(get("z", 2, 7)), "\n", sep = "")
      }
    }
    eqtable(x$outcome, x$rows$count, "count")
    if (any(c("NBRM", "ZINB") %in% sel)) {
      a <- function(i, d, w) sapply(sel, function(m) {
        e <- x$estimates[[m]]$lnalpha
        if (is.null(e)) "" else fw(e[[i]], d, w)
      })
      cat(rule, sprintf("%-31s|", "lnalpha"), "\n", lab30("Constant"), cells(a("b", 3, 9)), "\n",
        sprintf("%30s |", ""), cells(a("z", 2, 7)), "\n", sep = "")
    }
    if (!is.null(x$rows$zero)) eqtable("inflate", x$rows$zero, "zero")
    stat <- function(name, s) cat(lab30(name), cells(s), "\n", sep = "")
    cat(rule, sprintf("%-31s|", "Statistics"), "\n", sep = "")
    if ("NBRM" %in% sel) stat("alpha", ifelse(sel == "NBRM", fw(x$estimates$NBRM$lnalpha[["b"]], 3, 9), ""))
    stat("N", rep(x$N, length(sel)))
    stat("ll", fw(x$ll[sel], 3, 9))
    stat("bic", fw(x$fit$BIC, 3, 9))
    stat("aic", fw(x$fit$AIC, 3, 9))
    cat(strrep("-", w), "\n", sprintf("%*s", w, "Legend: b/t"), "\n", sep = "")
  }

  # Comparison of mean observed and predicted count
  if (!o$nodifferences) {
    cat("Comparison of Mean Observed and Predicted Count\n\n",
      "            Maximum       At      Mean\n",
      "Model     Difference    Value    |Diff|\n", strrep("-", 45), "\n", sep = "")
    for (m in sel) cat(sprintf("%-9s%8.3f%7s%-6s%8.3f\n", m, x$compare[m, "MaxDiff"], "",
      paste0("  ", x$compare[m, "AtValue"]), x$compare[m, "MeanAbsDiff"]))
  }

  # Predicted and actual probabilities
  if (!o$noprtable) {
    for (m in sel) {
      act <- x$probs$Actual
      prd <- x$probs[[m]]
      ad <- abs(act - prd)
      pear <- ifelse(prd > 0, x$N * ad^2 / prd, NA)
      cat("\n", m, ": Predicted and actual probabilities\n\n",
        "Count   Actual    Predicted    |Diff|   Pearson\n", strrep("-", 48), "\n", sep = "")
      cat(sprintf("%-6s%8.3f    %8.3f   %8.3f  %8s\n", x$probs$Count, act, prd, ad, fw(pear, 3, 8)), sep = "")
      cat(strrep("-", 48), "\n", sprintf("%-6s%8.3f    %8.3f   %8.3f  %8s\n", "Sum", sum(act), sum(prd),
        sum(ad), fw(if (anyNA(pear)) NA else sum(pear), 3, 8)), sep = "")
    }
  }

  # Tests and fit statistics
  if (!o$nofit) {
    bic <- setNames(x$fit$BIC, sel)
    aic <- setNames(x$fit$AIC, sel)
    line <- strrep("-", 73)
    header <- function(a) cat(sprintf("%-15sBIC=%10.3f  AIC=%10.3f  Prefer  Over  Evidence\n", a, bic[a], aic[a]))
    versus <- function(a, b) {
      tr <- x$tests[x$tests$Comparison == paste(a, "vs", b), ]
      cat(line, "\n", sprintf("  vs %-10sBIC=%10.3f  dif=%10.3f  %-8s%-6s%s\n", b, bic[b], tr$BICdif,
        tr$BICprefer, tr$BICover, tr$Evidence), sprintf("%15sAIC=%10.3f  dif=%10.3f  %-8s%s\n", "",
        aic[b], tr$AICdif, tr$AICprefer, tr$AICover), sep = "")
      if (paste(a, b) %in% c("PRM NBRM", "ZIP ZINB")) {
        lr <- x$lrtests[x$lrtests$Comparison == paste(a, "vs", b), ]
        cat(sprintf("%15sLRX2=%9.3f  prob=%9.3f  %-8s%-6sp=%.3f\n", "", lr$LRchi2, lr$p, lr$Prefer, lr$Over, lr$p))
      }
      if (paste(a, b) %in% c("PRM ZIP", "NBRM ZINB")) {
        vu <- if (o$forcevuong) x$vuong[x$vuong$Comparison == paste(a, "vs", b), ]
        if (!is.null(vu)) cat(sprintf("%15sVuong=%8.3f  prob=%9.3f  %-8s%-6sp=%.3f\n", "", vu$Vuong, vu$p,
          vu$Prefer, vu$Over, vu$p)) else cat(sprintf("%15sVuong=%8s  prob=%9s  %-8s%-6sp=.\n", "", ".", ".", b, a))
      }
    }
    cat("\nTests and Fit Statistics\n\n")
    if ("PRM" %in% sel) {
      header("PRM")
      for (b in intersect(c("NBRM", "ZIP", "ZINB"), sel)) versus("PRM", b)
    }
    if ("NBRM" %in% sel) {
      cat(line, "\n", sep = "")
      header("NBRM")
      for (b in intersect(c("ZIP", "ZINB"), sel)) versus("NBRM", b)
    }
    if ("ZIP" %in% sel) {
      cat(line, "\n", sep = "")
      header("ZIP")
      if ("ZINB" %in% sel) {
        versus("ZIP", "ZINB")
        cat(line, "\n", sep = "")
      }
      cat("\nVuong test is not appropriate for testing zero-inflation. To force the\n",
        "computation of the test, use forcevuong = TRUE.\n", sep = "")
    }
  }
  invisible(x)
}
