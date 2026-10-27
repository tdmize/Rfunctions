# Tests for countfit() in Count_helper_functions.R
# Run from the repository root: source("tests/test_countfit.R")
# Set CDA_DATA_DIR to a local folder holding cda_gss.dta, cda_ah4.dta, and cda_ah5.dta to run offline

suppressMessages({library(haven); library(ggplot2)})
source("Count_helper_functions.R")
ddir <- Sys.getenv("CDA_DATA_DIR", "https://tdmize.github.io/data/data")

res <- data.frame(test = character(), pass = logical())
check <- function(name, ok) {
  ok <- isTRUE(all(ok))
  res[nrow(res) + 1, ] <<- list(name, ok)
  cat(if (ok) "PASS  " else "FAIL  ", name, "\n", sep = "")
}
cfq <- function(...) { out <- NULL; capture.output(out <- countfit(...)); out }
err <- function(name, expr, pattern) {
  m <- tryCatch({ capture.output(expr); "no error" }, error = function(e) conditionMessage(e))
  check(paste0(name, " [", m, "]"), grepl(pattern, m))
}

# Data
gss_raw <- read_dta(file.path(ddir, "cda_gss.dta"))
gss <- gss_raw[c("mntlhlth", "woman", "married", "age", "faminc")]
gss <- gss[complete.cases(gss), ]
gss[c("woman", "married")] <- lapply(gss[c("woman", "married")], as_factor)
ah4 <- read_dta(file.path(ddir, "cda_ah4.dta"))
ah4$depscale <- round(ah4$depressSS * 10)
ah4 <- ah4[c("drinktime", "depscale", "woman", "age", "college", "race", "income",
  "spouserole", "parrole", "role")]
ah4[c("woman", "college", "race", "spouserole", "parrole")] <-
  lapply(ah4[c("woman", "college", "race", "spouserole", "parrole")], as_factor)
ah5 <- read_dta(file.path(ddir, "cda_ah5.dta"))

f_gss <- mntlhlth ~ married + age + woman + faminc
f_dr <- drinktime ~ woman + age + college + race + income + spouserole + parrole
f_dep <- depscale ~ woman + poly(role, 2, raw = TRUE) + income + age + race


# **********************************************************
# A - Matches Stata's countfit (countfit 3.2.2 in Stata 19.5; logs of 26 Oct 2026)
# **********************************************************
stata_check <- function(label, cf, compare, bic, aic, bicdif, aicdif, bicprefer, aicprefer, lr, probs = NULL) {
  check(paste(label, "comparison table"), all(abs(as.matrix(round(cf$compare, 3)) - matrix(compare, 4)) < 1e-9))
  check(paste(label, "BIC and AIC"), all(abs(cf$fit$BIC - bic) < .0015, abs(cf$fit$AIC - aic) < .0015))
  check(paste(label, "BIC and AIC differences"), all(abs(cf$tests$BICdif - bicdif) < .002,
    abs(cf$tests$AICdif - aicdif) < .002))
  check(paste(label, "preferred models (BIC and AIC)"), identical(cf$tests$BICprefer, bicprefer) &&
    identical(cf$tests$AICprefer, aicprefer))
  check(paste(label, "LR tests"), all(abs(cf$lrtests$LRchi2 - lr) < .002) &&
    all(cf$lrtests$p == 0))
  if (!is.null(probs)) check(paste(label, "actual and predicted probabilities"),
    all(abs(round(as.matrix(cf$probs[-1]), 3) - probs) < 1e-9))
}
pref <- c("NBRM", "ZIP", "ZINB", "NBRM", "ZINB", "ZINB")

cf <- cfq(f_gss, data = gss)
check("A1 GSS: N = 10405", cf$N == 10405)
stata_check("A1 GSS:", cf,
  c(.506, -.059, .073, .029, 0, 1, 2, 5, .114, .018, .038, .015),
  c(110666.341, 43281.780, 65418.038, 42800.494), c(110630.090, 43238.280, 65345.537, 42720.743),
  c(67384.561, 45248.303, 67865.847, -22136.258, 481.286, 22617.544),
  c(67391.811, 45284.553, 67909.347, -22107.258, 517.537, 22624.794), pref, pref, c(67393.811, 22626.794),
  cbind(c(.544, .052, .079, .056, .030, .059, .006, .020, .004, .001),
    c(.038, .108, .167, .185, .167, .129, .088, .055, .031, .017),
    c(.529, .112, .064, .045, .034, .027, .022, .018, .016, .013),
    c(.544, .002, .006, .014, .025, .039, .050, .057, .058, .054),
    c(.544, .073, .054, .043, .036, .030, .026, .022, .019, .017)))
check("A1 GSS: evidence", all(cf$tests$Evidence == "Very strong"))

cf_dr <- cfq(f_dr, data = ah4)
fit_dr <- list(bic = c(23916.196, 19024.594, 19367.272, 18418.232), aic = c(23851.058, 18952.943, 19236.996, 18281.443),
  bicdif = c(4891.602, 4548.924, 5497.964, -342.677, 606.362, 949.039),
  aicdif = c(4898.115, 4614.062, 5569.615, -284.053, 671.500, 955.553))
stata_check("A2 Add Health drinks:", cf_dr, c(.259, -.098, .037, .032, 0, 1, 2, 2, .060, .023, .015, .011),
  fit_dr$bic, fit_dr$aic, fit_dr$bicdif, fit_dr$aicdif, pref, pref, c(4900.115, 957.553),
  cbind(c(.418, .104, .152, .113, .070, .045, .032, .009, .017, .003),
    c(.159, .257, .234, .162, .095, .050, .024, .011, .004, .002),
    c(.385, .203, .126, .083, .056, .039, .028, .020, .015, .011),
    c(.418, .083, .115, .115, .096, .069, .046, .027, .015, .008),
    c(.420, .131, .120, .095, .070, .050, .035, .024, .017, .012)))

cf_dr15 <- cfq(f_dr, data = ah4, maxcount = 15)
stata_check("A3 drinks, maxcount = 15:", cf_dr15, c(.259, -.098, .037, .032, 0, 1, 2, 2, .039, .016, .011, .009),
  fit_dr$bic, fit_dr$aic, fit_dr$bicdif, fit_dr$aicdif, pref, pref, c(4900.115, 957.553))

cf_dep <- cfq(f_dep, data = ah4)
stata_check("A4 depressive symptoms (squared term):", cf_dep,
  c(-.068, -.010, .074, .009, 6, 3, 1, 1, .042, .004, .037, .004),
  c(32243.930, 27106.561, 31061.395, 27161.049), c(32185.465, 27041.599, 30944.464, 27037.621),
  c(5137.369, 1182.535, 5082.882, -3954.834, -54.488, 3900.347),
  c(5143.865, 1241.001, 5147.843, -3902.865, 3.978, 3906.843),
  c("NBRM", "ZIP", "ZINB", "NBRM", "NBRM", "ZINB"), pref, c(5145.865, 3908.843),
  cbind(c(.059, .090, .097, .098, .105, .095, .076, .073, .056, .051),
    c(.005, .022, .056, .097, .132, .148, .144, .123, .095, .068),
    c(.052, .088, .105, .109, .103, .093, .082, .069, .058, .048),
    c(.060, .015, .042, .078, .112, .134, .136, .123, .099, .074),
    c(.060, .081, .100, .106, .103, .094, .083, .071, .059, .049)))

cf_inf <- cfq(f_dr, data = ah4, inflate = ~ woman + age)
stata_check("A5 drinks, inflate = ~ woman + age:", cf_inf, c(.259, -.098, .041, .033, 0, 1, 2, 2, .060, .023, .017, .011),
  c(23916.196, 19024.594, 19671.285, 18734.419), c(23851.058, 18952.943, 19586.606, 18643.226),
  c(4891.602, 4244.911, 5181.777, -646.691, 290.175, 936.866),
  c(4898.115, 4264.452, 5207.832, -633.663, 309.717, 943.380), pref, pref, c(4900.115, 945.380))


# **********************************************************
# B - Independent calculation from model matrices and coefficients
# **********************************************************
oracle <- function(cf, y, maxcount = 9) {
  k <- 0:maxcount
  n <- length(y)
  m <- cf$models
  off <- function(o) if (is.null(o)) 0 else o
  pred <- df <- NULL
  for (nm in names(m)) {
    x <- m[[nm]]
    if (nm %in% c("PRM", "NBRM")) {
      b <- coef(x)
      keep <- !is.na(b)
      mu <- exp(model.matrix(x)[, keep, drop = FALSE] %*% b[keep] + off(model.offset(x$model)))
      p <- sapply(k, function(j) mean(if (nm == "PRM") dpois(j, mu) else dnbinom(j, size = x$theta, mu = mu)))
      d <- sum(keep) + (nm == "NBRM")
    } else {
      mu <- exp(model.matrix(x, "count") %*% coef(x, "count") + off(x$offset$count))
      pz <- plogis(model.matrix(x, "zero") %*% coef(x, "zero") + off(x$offset$zero))
      p <- sapply(k, function(j) mean(pz * (j == 0) + (1 - pz) *
        (if (nm == "ZIP") dpois(j, mu) else dnbinom(j, size = x$theta, mu = mu))))
      d <- length(coef(x)) + (nm == "ZINB")
    }
    pred <- cbind(pred, p)
    df <- c(df, d)
  }
  colnames(pred) <- names(m)
  obs <- as.numeric(table(factor(y, levels = k))) / n
  dif <- obs - pred
  ll <- sapply(m, function(x) as.numeric(logLik(x)))
  lr <- c(if (all(c("PRM", "NBRM") %in% names(m))) 2 * (ll["NBRM"] - ll["PRM"]),
    if (all(c("ZIP", "ZINB") %in% names(m))) 2 * (ll["ZINB"] - ll["ZIP"]))
  list(probs = cbind(Actual = obs, pred), bic = -2 * ll + df * log(n), aic = -2 * ll + 2 * df,
    maxdiff = apply(dif, 2, function(d) d[which.max(abs(d))]), meanabs = colMeans(abs(dif)), lr = pmax(lr, 0))
}
run_oracle <- function(label, cf, y, maxcount = 9) {
  o <- oracle(cf, y, maxcount)
  check(paste(label, "probabilities"), max(abs(as.matrix(cf$probs[-1]) - o$probs)) < 1e-8)
  check(paste(label, "BIC and AIC"), max(abs(cf$fit$BIC - o$bic), abs(cf$fit$AIC - o$aic)) < 1e-6)
  check(paste(label, "comparison table"), max(abs(cf$compare$MaxDiff - o$maxdiff),
    abs(cf$compare$MeanAbsDiff - o$meanabs)) < 1e-10)
  if (!is.null(cf$lrtests)) check(paste(label, "LR tests"), max(abs(cf$lrtests$LRchi2 - o$lr)) < .001)
  check(paste(label, "same N in every model"), all(sapply(cf$models, function(x) length(fitted(x))) == cf$N))
}

run_oracle("B1 GSS mental health:", cf, gss$mntlhlth)
dr <- ah4[complete.cases(ah4[all.vars(f_dr)]), ]
run_oracle("B2 Add Health drinks (4-category factor):", cf_dr, dr$drinktime)
run_oracle("B3 drinks, maxcount = 15:", cf_dr15, dr$drinktime, 15)
dep <- ah4[complete.cases(ah4[all.vars(f_dep)]), ]
run_oracle("B4 depressive symptoms (squared term):", cf_dep, dep$depscale)
run_oracle("B5 drinks, inflate = ~ woman + age:", cf_inf, dr$drinktime)
check("B5 inflation part has only woman and age",
  identical(names(coef(cf_inf$models$ZIP, "zero")), c("(Intercept)", "womanWomen", "age")) ||
  length(coef(cf_inf$models$ZIP, "zero")) == 3)
f_cig <- cigdays ~ male + race + edyrs + employed
cf_cig <- cfq(f_cig, data = ah5)
dcig <- as.data.frame(ah5)[complete.cases(as.data.frame(ah5)[all.vars(f_cig)]), ]
p0 <- glm(f_cig, family = poisson, data = dcig)
nll <- function(par) -sum(dnbinom(dcig$cigdays, size = exp(par[1]), mu = exp(model.matrix(p0) %*% par[-1]), log = TRUE))
o <- suppressWarnings(optim(c(0, coef(p0)), nll, method = "BFGS", control = list(maxit = 1000, reltol = 1e-12)))
check("B6 negative binomial reaches its maximum when glm.nb's default start fails (AH5 cigarette days)",
  abs(as.numeric(logLik(cf_cig$models$NBRM)) + o$value) < .01 && abs(cf_cig$models$NBRM$theta - exp(o$par[1])) < 1e-3)


# **********************************************************
# C - Input handling
# **********************************************************
gss_lab <- gss_raw[c("mntlhlth", "woman", "married", "age", "faminc")]
cf_raw <- cfq(f_gss, data = gss_lab)
check("C1 drops missing data itself (raw data with NAs gives N = 10405)", cf_raw$N == 10405)
check("C2 labelled numeric predictors match factor version",
  isTRUE(all.equal(cf_raw$compare, cf$compare)) && isTRUE(all.equal(cf_raw$fit, cf$fit)))
cf_df <- cfq(f_gss, data = as.data.frame(gss))
check("C3 plain data.frame matches tibble", isTRUE(all.equal(cf_df$fit, cf$fit)))
cf_dot <- cfq(mntlhlth ~ ., data = gss)
check("C4 formula with . matches", isTRUE(all.equal(cf_dot$fit, cf$fit)))
gss_inf <- gss_raw[c("mntlhlth", "woman", "married", "age", "faminc", "edyrs")]
cf_samp <- cfq(f_gss, data = gss_inf, inflate = ~ woman + edyrs)
check("C5 inflation-only variable with missing data: one sample for all models",
  all(sapply(cf_samp$models, function(x) length(fitted(x))) == cf_samp$N) &&
  cf_samp$N == sum(complete.cases(gss_inf)))
cf0 <- cfq(f_gss, data = gss, maxcount = 0)
check("C6 maxcount = 0 runs", nrow(cf0$probs) == 1 && all(cf0$compare$AtValue == 0))
cf40 <- cfq(f_gss, data = gss, maxcount = 40)
check("C7 maxcount above largest count runs", nrow(cf40$probs) == 41 && all(cf40$probs$Actual[32:41] == 0))
cf_int <- cfq(mntlhlth ~ 1, data = gss)
check("C8 intercept-only model runs", cf_int$N == 10405)
small <- gss[round(seq(1, nrow(gss), length.out = 300)), ]
cf_small <- suppressWarnings(cfq(f_gss, data = small))
check("C9 small sample (N = 300) runs", cf_small$N == 300)
gss_lby <- gss; gss_lby$mntlhlth <- gss_raw$mntlhlth[complete.cases(gss_raw[names(gss)])]
check("C10 labelled count outcome runs", inherits(cfq(f_gss, data = gss_lby), "countfit"))
check("C11 sums of predicted probabilities over 0-1000 are 1",
  all(abs(colSums(as.matrix(cfq(f_gss, data = gss, maxcount = 1000)$probs[-(1:2)])) - 1) < 1e-6))
noasian <- ah4[ah4$race != "Asian", ]
cf_lev <- cfq(f_dr, data = noasian)
check("C12 empty factor level gives same results as dropping it",
  isTRUE(all.equal(cf_lev$fit, cfq(f_dr, data = droplevels(noasian))$fit)))
set.seed(2); sim <- data.frame(x = rnorm(2000), ch = sample(c("a", "b", "c"), 2000, TRUE))
sim$y <- rpois(2000, exp(.2 + .3 * sim$x))
cf_pois <- suppressWarnings(cfq(y ~ x + ch, data = sim))
check("C13 Poisson data: prefers PRM, LR test of PRM vs NBRM not significant, character predictor works",
  cf_pois$tests$BICprefer[1] == "PRM" && cf_pois$lrtests$p[1] > .05)

# Output
pr <- capture.output(print(cf))
check("C14 printing the result shows the two tables, not the list",
  any(grepl("Comparison of Mean Observed", pr)) && any(grepl("Tests and Fit", pr)) && !any(grepl("\\$models", pr)))
check("C15 result is returned invisibly",
  { v <- NULL; capture.output(v <- withVisible(countfit(f_gss, data = gss))); !v$visible })
p <- cf$pic + labs(title = "Edited")
check("C16 graph is an editable ggplot", inherits(p, "ggplot") && p$labels$title == "Edited" &&
  nrow(ggplot_build(p)$data[[2]]) == 40)


# **********************************************************
# D - Errors give clear messages
# **********************************************************
g <- gss
g$frac <- g$mntlhlth + .5; err("D1 non-integer outcome", countfit(frac ~ age, data = g), "must be a count")
g$neg <- g$mntlhlth - 1; err("D2 negative outcome", countfit(neg ~ age, data = g), "must be a count")
g$nozero <- g$mntlhlth + 1; err("D3 no zeros", countfit(nozero ~ age, data = g), "no zeros.*models = c")
g$allzero <- 0; err("D4 all zeros", countfit(allzero ~ age, data = g), "0 for every")
g$fac <- factor(g$mntlhlth); err("D5 factor outcome", countfit(fac ~ age, data = g), "numeric count")
err("D6 | in formula", countfit(mntlhlth ~ age | woman, data = gss), "inflate =")
err("D7 fitted model instead of formula", countfit(cf$models$PRM, data = gss), "needs a formula")
err("D8 one-sided formula", countfit(~ age, data = gss), "needs a formula")
err("D9 inflate as text", countfit(f_gss, data = gss, inflate = "woman"), "one-sided formula")
err("D10 two-sided inflate", countfit(f_gss, data = gss, inflate = mntlhlth ~ woman), "one-sided formula")
err("D11 negative maxcount", countfit(f_gss, data = gss, maxcount = -1), "maxcount")
err("D12 fractional maxcount", countfit(f_gss, data = gss, maxcount = 2.5), "maxcount")
err("D13 maxcount as text", countfit(f_gss, data = gss, maxcount = "9"), "maxcount")
err("D14 misspelled variable", countfit(mntlhlth ~ agee, data = gss), "agee")
g$allna <- NA_real_; err("D15 no complete cases", countfit(mntlhlth ~ allna, data = g), "No observations")
g$const <- 1; err("D16 predictor that does not vary", countfit(mntlhlth ~ woman + const, data = g), "ZIP and ZINB models did not converge.*models = c\\(\"PRM\", \"NBRM\"\\).*do not vary")


# **********************************************************
# M - Choosing models with models =
# **********************************************************
cf3 <- cfq(f_gss, data = gss, models = c("PRM", "NBRM", "ZIP"))
check("M1 three models: same results as the four-model run, three comparisons",
  isTRUE(all.equal(cf3$fit, cf$fit[1:3, ])) && isTRUE(all.equal(cf3$compare, cf$compare[1:3, ])) &&
  nrow(cf3$tests) == 3 && nrow(cf3$lrtests) == 1 && identical(names(cf3$models), c("PRM", "NBRM", "ZIP")))
cf3b <- cfq(f_gss, data = gss, models = c("zip", "poisson", "nbreg"))
check("M2 Stata-style or lowercase names, any order", isTRUE(all.equal(cf3b$fit, cf3$fit)))
cf1 <- cfq(f_gss, data = gss, models = "NBRM")
check("M3 one model runs with no comparisons", is.null(cf1$tests) && nrow(cf1$compare) == 1 &&
  !any(grepl("Evidence", capture.output(print(cf1)))))
check("M4 LR test shown only when both nested models are included",
  is.null(cfq(f_gss, data = gss, models = c("NBRM", "ZINB"))$lrtests))
cf_nz <- cfq(nozero ~ age, data = g, models = c("PRM", "NBRM"))
check("M5 outcome with no zeros works without the zero-inflated models", cf_nz$N == 10405)
err("M6 unknown model name", countfit(f_gss, data = gss, models = c("PRM", "ZINBB")), "models can include")
err("M7 models not text", countfit(f_gss, data = gss, models = 1:3), "models can include")
suppressMessages(trace("zeroinfl", where = asNamespace("pscl"), print = FALSE,
  tracer = quote(if (dist == "negbin") stop("simulated failure"))))
err("M8 one model fails: message names it and gives the models = line",
  countfit(f_gss, data = gss), "^The ZINB model did not converge.*models = c\\(\"PRM\", \"NBRM\", \"ZIP\"\\)")
cf_rerun <- cfq(f_gss, data = gss, models = c("PRM", "NBRM", "ZIP"))
suppressMessages(untrace("zeroinfl", where = asNamespace("pscl")))
check("M9 rerunning with the suggested models works", identical(names(cf_rerun$models), c("PRM", "NBRM", "ZIP")))
check("M10 constant predictor: the two models that work run with models =",
  identical(names(cfq(mntlhlth ~ woman + const, data = g, models = c("PRM", "NBRM"))$models), c("PRM", "NBRM")))
for (tc in list(
  list("M11 ZINB stops before converging", "zeroinfl.control(maxit = 1, reltol = 1e-12)", "the estimation stopped before converging"),
  list("M12 ZINB stops at a poor solution: caught because it fits worse than NBRM", "zeroinfl.control(reltol = 0.5)",
    "fit worse than the simpler nested model"))) {
  suppressMessages(trace("zeroinfl", where = asNamespace("pscl"), print = FALSE,
    tracer = str2lang(sprintf("if (dist == \"negbin\") control <- %s", tc[[2]]))))
  err(tc[[1]], suppressWarnings(countfit(f_gss, data = gss)),
    paste0("^The ZINB model did not converge \\(ZINB: ", tc[[3]], "\\).*models = c\\(\"PRM\", \"NBRM\", \"ZIP\"\\)"))
  suppressMessages(untrace("zeroinfl", where = asNamespace("pscl")))
}


# **********************************************************
# F - A range of outcomes and predictors, checked against the independent calculation
# **********************************************************
gssf <- gss_raw
gssf$race3f <- as_factor(gssf$race3); gssf$regionf <- as_factor(gssf$region)
gssf$womanf <- as_factor(gssf$woman); gssf$marriedf <- as_factor(gssf$married)
gssf$polord <- factor(gssf$polviews, ordered = TRUE); gssf$yearf <- factor(gssf$year)
for (v in c("male", "race", "married", "college", "smoker")) ah5[[paste0(v, "f")]] <- as_factor(ah5[[v]])
specs <- list(
  list("F1 GSS children (squared age, 3-category race):", numchild ~ womanf + age + I(age^2) + edyrs + race3f, gssf, list()),
  list("F2 GSS TV hours (factor x continuous interaction, 9-category region, few zeros):", tvhours ~ womanf * age + edyrs + regionf, gssf, list()),
  list("F3 GSS physical health days (orthogonal cubic, log term):", physhlth ~ poly(age, 3) + womanf + log(faminc) + race3f, gssf, list()),
  list("F4 GSS male partners (long tail, maxcount = 20):", nummen ~ womanf + age + edyrs + marriedf, gssf, list(maxcount = 20)),
  list("F5 GSS sex partners (continuous x continuous interaction):", sexpartners ~ age * edyrs + womanf, gssf, list()),
  list("F6 GSS church attendance (ordered factor):", attend ~ polord + age + womanf, gssf, list()),
  list("F7 GSS friends (N = 913):", numfrend ~ age + womanf + edyrs, gssf, list()),
  list("F8 GSS mental health (year as a 30-level factor):", mntlhlth ~ yearf + womanf + age, gssf, list()),
  list("F9 AH5 e-cigarette days (94% zeros):", ecigmon ~ malef + age + edyrs + racef, ah5, list()),
  list("F10 AH5 golf/basketball (85% zeros, custom inflation):", golfbball ~ malef + income + edyrs, ah5, list(inflate = ~ malef)),
  list("F11 AH5 children (little overdispersion):", numchilds ~ marriedf + age + edyrs + malef, ah5, list()),
  list("F12 AH5 sodas (10 predictors, scattered missing data):",
    sodas ~ malef + age + edyrs + income + marriedf + collegef + racef + polviews + bmi + smokerf, ah5, list()),
  list("F13 AH5 drinks per sitting (offset term):", numdrinksit ~ malef + age + income + offset(log(edyrs)), ah5, list()))
for (s in specs) {
  cfs <- suppressWarnings(do.call(cfq, c(list(s[[2]], data = s[[3]]), s[[4]])))
  vars <- union(all.vars(s[[2]]), if (!is.null(s[[4]]$inflate)) all.vars(s[[4]]$inflate))
  d <- as.data.frame(s[[3]])
  d <- d[complete.cases(d[intersect(vars, names(d))]), ]
  check(paste(s[[1]], "N equals complete cases"), cfs$N == nrow(d))
  run_oracle(s[[1]], cfs, as.numeric(d[[all.vars(s[[2]])[1]]]), if (is.null(s[[4]]$maxcount)) 9 else s[[4]]$maxcount)
}


# **********************************************************
# E - Works in a fresh session without MASS or pscl loaded
# **********************************************************
fresh <- tempfile(fileext = ".R")
writeLines(c(sprintf('source("%s")', normalizePath("Count_helper_functions.R", winslash = "/")),
  sprintf('d <- haven::read_dta("%s")', file.path(ddir, "cda_gss.dta")),
  'cf <- countfit(mntlhlth ~ married + age + woman + faminc, data = d)',
  'cat("\\nFRESHOK", cf$N, !any(c("package:MASS", "package:pscl") %in% search()), "\\n")'), fresh)
out <- system2(file.path(R.home("bin"), "Rscript"), fresh, stdout = TRUE, stderr = TRUE)
check("E1 fresh session, packages not attached", any(grepl("FRESHOK 10405 TRUE", out)))


# **********************************************************
cat("\n", sum(res$pass), " of ", nrow(res), " tests passed\n", sep = "")
if (!all(res$pass)) stop("Some countfit() tests failed: ", paste(res$test[!res$pass], collapse = "; "))
