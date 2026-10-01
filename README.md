# R helper functions

Source `ME_helper_functions.R` to use `meineq()`, `totalme()`, `meineq_weights()`, and `me_vcov()`.
Examples and definitions are on the [ME inequality](https://www.trentonmize.com/software/meinequality_r/)
and [Total ME](https://www.trentonmize.com/software/totalme_r/) pages.

```r
library(marginaleffects) # Marginal effects and hypotheses
source("ME_helper_functions.R")
avg_comparisons(model, variables = list(race4 = "pairwise"), hypothesis = meineq(model))
avg_comparisons(model, variables = list(age = "sd", race4 = "pairwise"), hypothesis = totalme(model))
```

`meineq()` requires every category pair once per predictor, outcome, and subgroup.
`totalme()` requires every fitted outcome category once per contrast and subgroup,
on the response-probability scale with `comparison = "difference"`. It also accepts
a single selected contrast, such as a continuous change or a chosen category pair.
Its completeness, zero-sum, and probability-bound checks establish necessary
invariants; arbitrary input tables cannot establish which prediction scale or
comparison formula produced their numbers.

Category shares come from each model's whole estimation sample, including when
`by` is used. Shares use fitting or survey weights by default; `model_weights = FALSE`
ignores those weights for the shares. Weighting the averaging in `avg_comparisons()`
is a separate choice, supplied through its `wts` argument. Standard errors treat
the category shares as fixed.

For `survey::svyolr()` models, supply `vcov = me_vcov(model)`. The native covariance
orders slopes before cutpoints, while the `insight` parameter table used by
`marginaleffects` can put cutpoints first and add `Intercept:` to their labels.
With the tested versions, relying on the default covariance can produce incorrect
standard errors; supplying an unaligned native matrix can also cause an error.
`me_vcov()` matches parameters by name and returns the explicitly aligned matrix.
It works with ordinary and replicate-weight survey designs; for other model
classes it returns `stats::vcov(model)` unchanged.

For replicate-weight ordered models under `survey` 4.2.1, also avoid a source
data column named `w`: the package's internal `weights = w` argument can resolve
to that column instead of the replicate weights, making every replicate fit
identical. Rename the source column before constructing the design, or remove
it from the design's variables after its weights have been stored. This is a
model-fitting issue; covariance alignment cannot repair an already affected fit.

Some unusual outcome labels containing `Intercept:` and `|` can also make a
native cutpoint label equal a different prefixed cutpoint label. The tested
`marginaleffects` parameter-update method can then perturb the wrong cutpoint.
`me_vcov()` rejects that ambiguity: rename the outcome categories and refit before
requesting standard errors. It cannot repair an upstream parameter-update method.

```r
avg_comparisons(model, variables = list(race4 = "pairwise"),
  hypothesis = totalme(model), vcov = me_vcov(model), wts = TRUE)
```

The summaries use absolute values. Normal-Wald confidence intervals and p-values
can be unreliable when component effects are zero or close to zero. Numerical
agreement of estimates and delta-method standard errors does not establish coverage
or calibration in those settings. A test of zero inequality can instead start
with a joint test that the underlying category effects are equal; this does not
provide a general interval method for differences of summary measures.

## Validation

Run from the repository root:

```r
source("tests/test_ME_helpers.R")
```

The independent numerical oracle uses model matrices, link functions, central
differences, and coefficient covariance matrices. It checks linear, binary,
count, multinomial, ordered, weighted, survey, subgroup, and combined-model
calculations, as well as incomplete-input rejection. Install `marginaleffects`,
`insight`, `nnet`, `MASS`, `sandwich`, and `numDeriv`; install `survey`, `ordinal`, and
the updated `suest` to run their integration cases too. The script writes `test_results.csv`
and `sessionInfo.txt` to `ME_TEST_OUTPUT` (the current directory by default), and
exits with an error if any scenario fails.

The GitHub workflow runs the complete helper suite on Windows and Linux with
`marginaleffects` 0.32.0 and 1.0.0. Install or merge the `suest` categorical sample
alignment fix first. A manual workflow run can choose its `suest_ref` explicitly;
automatic runs use `tdmize/suest` main. Workflow configuration is not evidence
that those platform runs have passed.
