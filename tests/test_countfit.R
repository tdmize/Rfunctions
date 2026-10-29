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

cf_dep <- suppressWarnings(cfq(f_dep, data = ah4))  # inflation coefficient for Asian respondents runs to -Inf
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
g$const <- 1; err("D16 predictor that does not vary", countfit(mntlhlth ~ woman + const, data = g), "do not vary or are perfectly collinear with other predictors: const\\.")
g$age2 <- g$age * 2; err("D17 perfectly collinear predictor", countfit(mntlhlth ~ age + age2, data = g), "perfectly collinear with other predictors: age2\\.")


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
  !any(grepl("  vs ", capture.output(print(cf1)))))
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
err("M10 collinear predictor is caught for any choice of models", countfit(mntlhlth ~ age + age2, data = g,
  models = c("PRM", "NBRM")), "perfectly collinear with other predictors: age2\\.")
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
# S - Output matches Stata's countfit (countfit.ado 3.2.2)
# **********************************************************
# S1: the printed output for the GSS model, line for line. Allowed differences: line 7, a z value
# Stata shows as -4.24 (R: -4.2349, confirmed by a separate calculation); lines 97 and 103, Pearson
# values that differ by 0.001; line 157, the note gives R's syntax
stata_gss <- strsplit(r"---(--------------------------------------------------------------------------------
                      Variable |    PRM        NBRM         ZIP        ZINB
-------------------------------+------------------------------------------------
mntlhlth                       |
             Currently Married |
                      Married  |     0.775       0.765       0.885       0.847
                               |    -23.00       -5.96      -11.00       -4.24
             age of respondent |     0.991       0.990       0.999       0.999
                               |    -23.98       -6.02       -2.17       -0.65
                               |
                       Woman R |
                        Women  |     1.291       1.281       0.974       0.962
                               |     24.81        5.77       -2.55       -1.02
        Family Income in 1000s |     0.996       0.997       0.997       0.996
                               |    -21.23       -5.42      -19.17       -6.59
                      Constant |     6.139       6.293      10.329       8.383
                               |    104.06       23.74      130.79       31.19
-------------------------------+------------------------------------------------
lnalpha                        |
                      Constant |                 4.449                   1.571
                               |                 78.37                   10.09
-------------------------------+------------------------------------------------
inflate                        |
                               |
             Currently Married |
                      Married  |                             1.278       1.283
                               |                              5.73        4.25
             age of respondent |                             1.016       1.021
                               |                             10.57        9.80
                               |
                       Woman R |
                        Women  |                             0.588       0.494
                               |                            -13.16      -12.30
        Family Income in 1000s |                             1.001       1.000
                               |                              1.04       -0.51
                      Constant |                             0.691       0.375
                               |                             -5.20       -8.89
-------------------------------+------------------------------------------------
Statistics                     |
                         alpha |                 4.449
                             N |     10405       10405       10405       10405
                            ll | -5.53e+04   -2.16e+04   -3.27e+04   -2.13e+04
                           bic |  1.11e+05   43281.780   65418.038   42800.494
                           aic |  1.11e+05   43238.280   65345.537   42720.743
--------------------------------------------------------------------------------
                                                                     Legend: b/t
Comparison of Mean Observed and Predicted Count

            Maximum       At      Mean
Model     Difference    Value    |Diff|
---------------------------------------------
PRM         0.506         0      0.114
NBRM       -0.059         1      0.018
ZIP         0.073         2      0.038
ZINB        0.029         5      0.015

PRM: Predicted and actual probabilities

Count   Actual    Predicted    |Diff|   Pearson
------------------------------------------------
0        0.544       0.038      0.506   7.0e+04
1        0.052       0.108      0.056   301.567
2        0.079       0.167      0.088   477.827
3        0.056       0.185      0.130   941.939
4        0.030       0.167      0.136  1158.285
5        0.059       0.129      0.069   386.307
6        0.006       0.088      0.082   796.725
7        0.020       0.055      0.035   233.471
8        0.004       0.031      0.027   242.103
9        0.001       0.017      0.016   159.483
------------------------------------------------
Sum      0.852       0.985      1.144   7.4e+04

NBRM: Predicted and actual probabilities

Count   Actual    Predicted    |Diff|   Pearson
------------------------------------------------
0        0.544       0.529      0.015     4.452
1        0.052       0.112      0.059   326.578
2        0.079       0.064      0.015    37.394
3        0.056       0.045      0.011    29.256
4        0.030       0.034      0.003     3.415
5        0.059       0.027      0.033   413.697
6        0.006       0.022      0.016   121.297
7        0.020       0.018      0.001     1.024
8        0.004       0.016      0.011    84.722
9        0.001       0.013      0.013   125.840
------------------------------------------------
Sum      0.852       0.879      0.178  1147.676

ZIP: Predicted and actual probabilities

Count   Actual    Predicted    |Diff|   Pearson
------------------------------------------------
0        0.544       0.544      0.000     0.000
1        0.052       0.002      0.051   1.5e+04
2        0.079       0.006      0.073  9485.500
3        0.056       0.014      0.042  1316.476
4        0.030       0.025      0.005    10.330
5        0.059       0.039      0.021   118.904
6        0.006       0.050      0.044   404.133
7        0.020       0.057      0.038   255.836
8        0.004       0.058      0.054   521.968
9        0.001       0.054      0.054   549.816
------------------------------------------------
Sum      0.852       0.849      0.381   2.8e+04

ZINB: Predicted and actual probabilities

Count   Actual    Predicted    |Diff|   Pearson
------------------------------------------------
0        0.544       0.544      0.000     0.000
1        0.052       0.073      0.021    61.714
2        0.079       0.054      0.025   119.158
3        0.056       0.043      0.013    37.549
4        0.030       0.036      0.005     8.146
5        0.059       0.030      0.029   298.459
6        0.006       0.026      0.020   157.634
7        0.020       0.022      0.002     2.755
8        0.004       0.019      0.015   119.812
9        0.001       0.017      0.016   160.449
------------------------------------------------
Sum      0.852       0.865      0.146   965.676

Tests and Fit Statistics

PRM            BIC=110666.341  AIC=110630.090  Prefer  Over  Evidence
-------------------------------------------------------------------------
  vs NBRM      BIC= 43281.780  dif= 67384.561  NBRM    PRM   Very strong
               AIC= 43238.280  dif= 67391.811  NBRM    PRM
               LRX2=67393.811  prob=    0.000  NBRM    PRM   p=0.000
-------------------------------------------------------------------------
  vs ZIP       BIC= 65418.038  dif= 45248.303  ZIP     PRM   Very strong
               AIC= 65345.537  dif= 45284.553  ZIP     PRM
               Vuong=       .  prob=        .  ZIP     PRM   p=.
-------------------------------------------------------------------------
  vs ZINB      BIC= 42800.494  dif= 67865.847  ZINB    PRM   Very strong
               AIC= 42720.743  dif= 67909.347  ZINB    PRM
-------------------------------------------------------------------------
NBRM           BIC= 43281.780  AIC= 43238.280  Prefer  Over  Evidence
-------------------------------------------------------------------------
  vs ZIP       BIC= 65418.038  dif=-22136.258  NBRM    ZIP   Very strong
               AIC= 65345.537  dif=-22107.258  NBRM    ZIP
-------------------------------------------------------------------------
  vs ZINB      BIC= 42800.494  dif=   481.286  ZINB    NBRM  Very strong
               AIC= 42720.743  dif=   517.537  ZINB    NBRM
               Vuong=       .  prob=        .  ZINB    NBRM  p=.
-------------------------------------------------------------------------
ZIP            BIC= 65418.038  AIC= 65345.537  Prefer  Over  Evidence
-------------------------------------------------------------------------
  vs ZINB      BIC= 42800.494  dif= 22617.544  ZINB    ZIP   Very strong
               AIC= 42720.743  dif= 22624.794  ZINB    ZIP
               LRX2=22626.794  prob=    0.000  ZINB    ZIP   p=0.000
-------------------------------------------------------------------------

Vuong test is not appropriate for testing zero-inflation. To force the
the computation of the test, use option -forcevuong-.)---", "\n")[[1]]
r_gss <- sub("[[:space:]]+$", "", capture.output(print(cf)))[-1]
check("S1 GSS model: printed output matches Stata line for line (157 lines)", length(r_gss) == 157 &&
  length(stata_gss) == 157 && identical(which(r_gss != stata_gss), c(7L, 97L, 103L, 157L)))

# S2: GSS 2002-2014 model with squared terms and forcevuong (Stata log of 28 Oct 2026)
g12 <- gss_raw[gss_raw$year >= 2002 & gss_raw$year <= 2014, ]
g12$notrelig <- ifelse(g12$relig %in% c(1, 2, 11), 0, ifelse(g12$relig == 4, 1, NA))
g12 <- g12[complete.cases(g12[c("mntlhlth", "physhlth", "woman", "race3", "year", "college", "parent",
  "faminc", "metro", "notrelig", "married")]), ]
g12[c("woman", "married", "race3", "college", "parent")] <- lapply(g12[c("woman", "married", "race3", "college", "parent")], as_factor)
g12$notrelig <- factor(g12$notrelig, 0:1, c("Christian", "Not Religious"))
g12$year <- factor(g12$year)
cf12 <- cfq(mntlhlth ~ woman + married + poly(age, 2, raw = TRUE) + poly(faminc, 2, raw = TRUE) +
    race3 + college + parent + notrelig + year, data = g12,
  inflate = ~ woman + married + poly(age, 2, raw = TRUE) + faminc + race3 + college + parent + notrelig + year,
  forcevuong = TRUE)
r2 <- function(e, eq, what) round(unname(e[[eq]][[what]]), c(b = 3, z = 2)[what])
same_num <- function(a, b) length(a) == length(b) && max(abs(a - b)) < 1e-9
check("S2 forcevuong model: N, log likelihoods, BIC, AIC", cf12$N == 5036 &&
  all(abs(cf12$fit$BIC - c(52903.251, 19473.344, 28702.860, 19137.134)) < .0015) &&
  all(abs(cf12$fit$AIC - c(52805.385, 19368.954, 28513.653, 18941.403)) < .0015) &&
  all(abs(cf12$ll[c("NBRM", "ZINB")] - c(-9668.477, -9440.701)) < .0015))
check("S2 forcevuong model: Vuong tests", all(abs(cf12$vuong$Vuong - c(32.061, 12.249)) < .0015) &&
  identical(cf12$vuong$Prefer, c("ZIP", "ZINB")) && all(cf12$vuong$p == 0))
check("S2 forcevuong model: exponentiated coefficients and z values, all four models",
  same_num(r2(cf12$estimates$PRM, "count", "b"), c(6.630, 1.310, 0.778, 0.999, 1.000, 0.988, 1.000, 0.681,
    0.900, 0.781, 1.205, 1.188, 0.758, 0.951, 0.863)) &&
  same_num(r2(cf12$estimates$NBRM, "count", "z"), c(5.55, 4.23, -3.33, -0.01, -0.68, -3.37, 2.09, -3.65,
    -1.05, -3.80, 1.47, 1.85, -3.44, -0.64, -1.51)) &&
  same_num(r2(cf12$estimates$ZIP, "count", "z"), c(41.70, -5.59, -2.98, -7.63, 7.45, -12.02, 8.88, -1.33,
    0.73, -16.39, 5.34, 1.67, -9.45, -5.73, -2.23)) &&
  same_num(r2(cf12$estimates$ZINB, "zero", "b"), c(2.407, 0.469, 1.405, 0.946, 1.001, 1.004, 1.974, 1.305,
    0.745, 0.844, 0.736, 1.125, 0.800, 1.214)) &&
  same_num(r2(cf12$estimates$ZINB, "zero", "z"), c(2.39, -10.06, 4.15, -3.09, 4.38, 3.08, 6.47, 1.99,
    -3.26, -1.87, -3.26, 1.25, -2.08, 1.92)) &&
  same_num(unname(round(c(cf12$estimates$NBRM$lnalpha, cf12$estimates$ZINB$lnalpha), c(3, 2, 3, 2))),
    c(5.404, 58.61, 1.331, 4.68)))

# S3: GSS 2010 on, number of children (Stata log of 28 Oct 2026). Stata's countfit prints the
# Prefer columns and p= of the ZIP vs ZINB line from the PRM vs NBRM test (countfit.ado lines 564-574);
# R reports that test's own result
gch <- gss_raw[gss_raw$year >= 2010, ]
gch[c("woman", "college")] <- lapply(gch[c("woman", "college")], as_factor)
cfch <- cfq(numchild ~ woman + college + faminc, data = gch, forcevuong = TRUE)
check("S3 number of children: N, BIC, AIC, comparison table", cfch$N == 14072 &&
  all(abs(cfch$fit$BIC - c(50523.434, 49541.455, 48412.931, 48420.957)) < .0015) &&
  all(abs(cfch$fit$AIC - c(50493.226, 49503.695, 48352.516, 48352.989)) < .0015) &&
  all(abs(as.matrix(round(cfch$compare, 3)) - matrix(c(-.133, -.117, .055, .056, 1, 1, 2, 2, .027, .027, .013, .013), 4)) < 1e-9))
check("S3 number of children: LR and Vuong tests", all(abs(cfch$lrtests$LRchi2 - c(991.532, 1.527)) < .0015) &&
  same_num(round(cfch$lrtests$p, 3), c(0, .108)) && identical(cfch$lrtests$Prefer, c("NBRM", "ZIP")) &&
  all(abs(cfch$vuong$Vuong - c(22.532, 18.536)) < .0015))
check("S3 number of children: printed ZIP vs ZINB LR line",
  any(capture.output(print(cfch)) == "               LRX2=    1.527  prob=    0.108  ZIP     ZINB  p=0.108"))
check("S3 number of children: z values (ZIP and ZINB)",
  same_num(r2(cfch$estimates$ZIP, "count", "z"), c(68.30, -2.82, -7.84, -8.77)) &&
  same_num(r2(cfch$estimates$ZINB, "zero", "z"), c(-11.07, -13.60, 13.52, -12.01)) &&
  same_num(round(unname(cfch$estimates$ZINB$lnalpha), c(3, 2)), c(0.009, -5.73)))

# S4: no overdispersion (Add Health male partners; Stata's countfit stops with "convergence not achieved")
ah5m <- ah5
ah5m$male <- as_factor(ah5m$male)
w <- character()
cfm <- withCallingHandlers(cfq(maleyr ~ male + age + edyrs, data = ah5m),
  warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
check("S4 no overdispersion: negative binomial at the Poisson boundary (alpha = 0, LR = 0, BIC higher by log N)",
  is.infinite(cfm$models$NBRM$theta) && abs(cfm$ll[["NBRM"]] - cfm$ll[["PRM"]]) < 1e-9 &&
  abs(cfm$fit["NBRM", "BIC"] - cfm$fit["PRM", "BIC"] - log(cfm$N)) < 1e-9 && cfm$lrtests$LRchi2[1] == 0 &&
  cfm$lrtests$Prefer[1] == "PRM")
check("S4 inflation coefficient running to infinity: warning for the zero-inflated models",
  sum(grepl("extreme estimates", w)) == 2)

# S5: options that skip sections, as in Stata
po <- function(...) capture.output(cfq(f_gss, data = gss, ...))
has <- function(out) c(est = any(grepl("Legend: b/t", out)), dif = any(grepl("Comparison of Mean", out)),
  prt = any(grepl("Predicted and actual", out)), fit = any(grepl("Tests and Fit", out)))
check("S5 noestimates, nodifferences, noprtable, and nofit each skip only their section",
  identical(unname(has(po(noestimates = TRUE))), c(FALSE, TRUE, TRUE, TRUE)) &&
  identical(unname(has(po(nodifferences = TRUE))), c(TRUE, FALSE, TRUE, TRUE)) &&
  identical(unname(has(po(noprtable = TRUE))), c(TRUE, TRUE, FALSE, TRUE)) &&
  identical(unname(has(po(nofit = TRUE))), c(TRUE, TRUE, TRUE, FALSE)))
err("S5 option that is not TRUE or FALSE", countfit(f_gss, data = gss, noprtable = "yes"), "noprtable must be TRUE or FALSE")


# **********************************************************
# R - Robustness fixes (stress testing, 28 Oct 2026)
# **********************************************************
# R1: an offset applies to the count part only, as with Stata's exposure()
cfo <- cfq(numdrinksit ~ malef + age + income + offset(log(edyrs)), data = ah5)
check("R1 offset in the count part only", !is.null(cfo$models$ZIP$offset$count) && is.null(cfo$models$ZIP$offset$zero) &&
  !is.null(cfo$models$ZINB$offset$count) && is.null(cfo$models$ZINB$offset$zero))
# R2: rescaling predictors changes no fit or z value (the zero-inflated negative binomial stalled at the
# plain negative binomial when age and education were divided by 10, before the EM refit was added)
gsp <- gss_raw[c("sexpartners", "age", "edyrs", "woman")]
gsp$woman <- as_factor(gsp$woman); gsp$age10 <- gsp$age / 10; gsp$edyrs10 <- gsp$edyrs / 10
cfa <- suppressWarnings(cfq(sexpartners ~ age * edyrs + woman, data = gsp, models = c("NBRM", "ZINB")))
cfb <- suppressWarnings(cfq(sexpartners ~ age10 * edyrs10 + woman, data = gsp, models = c("NBRM", "ZINB")))
zv <- function(cf) c(cf$estimates$ZINB$count$z, cf$estimates$ZINB$zero$z, cf$estimates$ZINB$lnalpha[["z"]])
check("R2 rescaled predictors: same zero-inflated negative binomial fit and z values",
  abs(cfa$ll[["ZINB"]] - cfb$ll[["ZINB"]]) < 1e-4 && max(abs(zv(cfa) - zv(cfb))) < .005)
g2 <- gss
g2$famdollars <- as.numeric(g2$faminc) * 1000
cfd <- cfq(mntlhlth ~ married + age + woman + famdollars, data = g2)
check("R2 income in dollars instead of thousands: same fits and z values",
  max(abs(cfd$ll - cf$ll)) < 1e-4 && max(abs(unlist(lapply(cfd$estimates, function(e) c(e$count$z, e$zero$z))) -
    unlist(lapply(cf$estimates, function(e) c(e$count$z, e$zero$z))))) < .002)
# R3: variable names with spaces still get their labels; overflowed coefficients print as . (as Stata)
g3 <- gss
names(g3)[names(g3) == "age"] <- "age years"
check("R3 variable name with a space keeps its label",
  any(grepl("age of respondent", capture.output(cfq(mntlhlth ~ married + `age years`, data = g3, nofit = TRUE)))))
g3$fam_b <- as.numeric(g3$faminc) / 1e6
o3 <- capture.output(suppressWarnings(cfq(mntlhlth ~ married + fam_b, data = g3, models = c("ZIP"), noprtable = TRUE)))
check("R3 an exponentiated coefficient too large to show prints as .", !any(grepl("Inf", o3)))
# R4: the warning for degenerate fits fires for perfect prediction (S4) but not for a smooth curve
# that predicts near-zero inflation only at the most extreme ages
g4 <- gss_raw[c("numchild", "woman", "age", "edyrs", "race3")]
g4[c("woman", "race3")] <- lapply(g4[c("woman", "race3")], as_factor)
w4 <- character()
invisible(withCallingHandlers(cfq(numchild ~ woman + poly(age, 2, raw = TRUE) + edyrs + race3, data = g4, models = "ZINB"),
  warning = function(x) { w4 <<- c(w4, conditionMessage(x)); invokeRestart("muffleWarning") }))
check("R4 no warning for a smooth age curve with a few extreme predictions", length(w4) == 0)
# R5: a zero-inflated negative binomial whose optimizer stalls below the zero-inflated Poisson
# (simulated data, n = 100, no real zero inflation); countfit refits from the simpler model's estimates
set.seed(1012)
d5 <- data.frame(age = runif(100, 18, 80), female = factor(rbinom(100, 1, .5), 0:1, c("Man", "Woman")),
  income = round(rlnorm(100, log(40), .8), 1), group = factor(sample(c("A", "B", "C"), 100, TRUE, prob = c(.6, .3, .1))))
d5$y <- rnbinom(100, size = 2, mu = exp(log(.5) - .01 * (d5$age - 45) + .3 * (d5$female == "Woman") -
  .004 * (d5$income - 40) + .2 * (d5$group == "B")))
cf5 <- suppressWarnings(cfq(y ~ poly(age, 2, raw = TRUE) + income, data = d5))
check("R5 stalled zero-inflated fit is refit: it is at least as good as the models nested in it",
  cf5$ll[["ZINB"]] >= max(cf5$ll[c("NBRM", "ZIP")]) - 1e-4 && cf5$ll[["ZIP"]] >= cf5$ll[["PRM"]] - 1e-4)
# R6: a zero-inflated Poisson that stops short of its maximum, found because the zero-inflated
# negative binomial reaches the Poisson boundary (alpha near 0) with a higher log likelihood
set.seed(1016)
d6 <- data.frame(age = runif(100, 18, 80), female = factor(rbinom(100, 1, .5), 0:1, c("Man", "Woman")),
  income = round(rlnorm(100, log(40), .8), 1), group = factor(sample(c("A", "B", "C"), 100, TRUE, prob = c(.6, .3, .1))))
d6$y <- rnbinom(100, size = 2, mu = exp(log(.5) - .01 * (d6$age - 45) + .3 * (d6$female == "Woman") -
  .004 * (d6$income - 40) + .2 * (d6$group == "B")))
d6$y[runif(100) < plogis(qlogis(.5) + .5 * (d6$female == "Man"))] <- 0
cf6 <- suppressWarnings(cfq(y ~ age + female + income + group, data = d6))
check("R6 zero-inflated negative binomial at alpha = 0 and zero-inflated Poisson have the same fit",
  cf6$models$ZINB$theta > 1000 && abs(cf6$ll[["ZINB"]] - cf6$ll[["ZIP"]]) < .01)


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
