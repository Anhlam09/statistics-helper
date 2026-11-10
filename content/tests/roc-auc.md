# ROC curve, AUC, sensitivity & specificity

> Evaluates how well a continuous test or model score separates patients with and without a disease, by plotting sensitivity against 1 − specificity across every possible cutoff and summarising discrimination as the area under the curve (AUC).

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Judge how well a continuous test or score separates diseased from healthy | Continuous or ordinal score; binary reference standard | Valid reference standard; independent patients; representative sample | High-sensitivity troponin to detect myocardial infarction |

## Learning objectives

- Explain sensitivity, specificity, PPV, and NPV, and how the ROC curve plots them over all cutoffs
- Interpret the AUC as a probability and read its DeLong CI
- Run ROC analysis in R and Python, and choose a cutoff while understanding the Youden index's limits
- Report accuracy with CIs and prevalence, and recognise that a cutoff chosen in the same data is optimistic

## What it is

A diagnostic test or prediction score is usually continuous (troponin, a risk-model probability), but decisions are binary. For any cutoff, patients scoring at or above it are called "positive", which gives a **sensitivity** (proportion of diseased patients called positive) and a **specificity** (proportion of non-diseased called negative). The **receiver operating characteristic (ROC) curve** plots sensitivity (y) against 1 − specificity (x) as the cutoff moves from high to low. A useless test follows the diagonal; a perfect test hugs the top-left corner.

The **AUC** is the probability that a randomly chosen diseased patient has a higher score than a randomly chosen non-diseased patient (ties count as one half). It equals the Mann–Whitney U statistic divided by the number of case–control pairs, so AUC = 0.5 means no discrimination and 1.0 perfect discrimination.

- **Test statistic:** AUC with standard error from the **DeLong** method; z = (AUC − 0.5) / SE tests H₀: AUC = 0.5
- **Effect size:** **AUC with 95% CI** (primary). At a chosen cutoff, report **sensitivity, specificity, PPV and NPV, each with a 95% CI**. A common (not the only) way to choose the cutoff is the **Youden index**, J = sensitivity + specificity − 1, maximised over cutoffs.
- **Definitions at a cutoff:** sensitivity = TP / (TP + FN); specificity = TN / (TN + FP); PPV = TP / (TP + FP); NPV = TN / (TN + FN). PPV and NPV depend on **prevalence** in the study sample; sensitivity and specificity do not.

## The intuition

Think of a smoke alarm. Set it very sensitive and it catches every fire (high sensitivity) but also rings for toast (low specificity). Set it dull and toast is ignored, but so, eventually, is a fire. **Every cutoff is a trade-off**, and the ROC curve shows all the trade-offs at once by sliding the cutoff from high to low.

Here the signal is the gap between the troponin of patients with and without MI; the noise is how much the two groups overlap. The AUC summarises the overlap as a single probability: pick one MI patient and one non-MI patient at random, and AUC is the chance the MI patient has the higher troponin (0.82 here). Perfect separation would give 1.0, complete overlap 0.5.

<div class="viz" data-viz="dots" data-unit="troponin, ng/L" data-groups='{"No MI":[9.7,4.9,12.8,5.9,5.8,6.1,22.3,3.7,67.8,20.3,4.5,22.3,8.8,19.3,2.9,5.2,4.7,34.0,8.4,19.3,9.7,3.9,7.4,21.2,9.4,12.4,11.5,12.9,7.4,2.2,5.5,29.5,54.9,8.5,4.4],"MI":[22.2,22.9,44.9,37.5,11.0,27.9,15.7,11.5,38.4,23.0,32.6,19.8,20.0,48.3,24.1,13.8,10.7,46.7,13.7,95.6,137.3,18.4,30.9,6.3,35.3]}' data-caption="Notice that MI patients sit higher (median 23.0 vs 8.8 ng/L) but the two groups overlap in the 10-30 ng/L zone. A cutoff near 10 ng/L misses only one MI patient (high sensitivity) but flags 14 of 35 patients without MI (modest specificity). Moving the cutoff right trades sensitivity for specificity."></div>

The AUC's DeLong standard error is about 0.054, so under H₀ (AUC = 0.5) the z statistic is (0.818 − 0.5) / 0.054 ≈ 5.8. This is far into the tail of the normal curve below. (This value is derived from the printed CI, not from a separate output line.)

<div class="viz" data-viz="z" data-stat="5.84"></div>

A far-out z only says the test is better than chance. How *good* it is for patients is a separate question that the AUC, sensitivity, specificity and predictive values answer.

## When to use it

Use ROC analysis when:

- You have a **continuous or ordinal** test result or model score and a **binary reference-standard diagnosis** (MI yes/no, cancer yes/no)
- You want to describe the test's **discrimination** across all cutoffs, compare tests, or choose a cutoff
- The reference standard was applied independently of the test result (to avoid verification bias)

**Clinical scenario.** In 60 emergency-department patients with chest pain, high-sensitivity troponin (ng/L) is measured at presentation. The final adjudicated diagnosis shows 25 had a myocardial infarction (MI). *How well does troponin discriminate MI, and what cutoff gives a good balance between sensitivity and specificity?*

**Don't use it when:**

- The test result is already binary → report sensitivity/specificity directly (with an [exact binomial CI](#/test/one-proportion)) or compare using [McNemar's test](#/test/mcnemar) for paired tests
- You need to quantify the **adjusted effect** of risk factors rather than discrimination → use [logistic regression](#/test/logistic-regression)
- The outcome is time to an event with censoring → use time-dependent ROC or the concordance index from [Cox regression](#/test/cox-ph)
- You want agreement between two raters or methods → use [Cohen's kappa](#/test/cohens-kappa) or [Bland–Altman analysis](#/test/bland-altman)

## Assumptions

1. **A valid binary reference standard** applied to all patients (or a random subset). A flawed gold standard bounds the apparent accuracy.
2. **Independent patients.** One result per patient; clustered or repeated data need cluster-aware methods.
3. **A representative sample.** Sensitivity and specificity may vary with disease spectrum; PPV and NPV depend on prevalence, so they transfer only to settings with similar prevalence.
4. **Enough patients in each group.** There is no magic number, but with only a few dozen per group the CIs are wide and the cutoff is unstable. The optimal cutoff chosen in the same data is **optimistic** and should be validated externally or by bootstrap/cross-validation.
5. **Direction of the test is known.** Decide beforehand whether higher values indicate disease (`direction` in pROC).
6. **For a prediction model:** the score should come from a model assessed on new patients, not only on the development data. AUC measures discrimination, not calibration.
7. **Model diagnostics.** When the score is a fitted logistic model, check its own assumptions (linearity in the logit, multicollinearity, influential points, enough events per variable; the traditional guide is about 10, but it is only a rough guide) as in [logistic regression](#/test/logistic-regression).

## Hypotheses

- **H₀:** AUC = 0.5 (the test cannot distinguish diseased from non-diseased patients)
- **H₁:** AUC ≠ 0.5 (two-sided). For comparing two correlated AUCs: H₀: AUC₁ = AUC₂ (DeLong test).

The cutoff, sensitivity and specificity are estimated descriptively; no hypothesis test is attached to them.

## Example data

Simulated data for 60 patients: `troponin` (ng/L) and `mi` (1 = MI, 0 = no MI; 25 cases, 35 controls). The first five rows:

| Patient | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| Troponin (ng/L) | 9.7 | 4.9 | 22.2 | 12.8 | 5.9 |
| MI (1 = yes) | 0 | 0 | 1 | 0 | 0 |

All 60 rows are typed inline in each code block below.

## How to run it

### R

```r
# Troponin (ng/L) at presentation and final diagnosis of MI (1 = yes), 60 patients
dat <- data.frame(
  troponin = c(9.7, 4.9, 22.2, 12.8, 5.9, 5.8, 6.1, 22.3, 3.7, 67.8, 22.9,
      20.3, 4.5, 44.9, 22.3, 8.8, 37.5, 11, 19.3, 27.9, 2.9, 5.2, 15.7, 4.7,
      34, 8.4, 19.3, 9.7, 3.9, 11.5, 38.4, 7.4, 21.2, 9.4, 23, 32.6, 19.8,
      12.4, 11.5, 12.9, 7.4, 20, 2.2, 48.3, 5.5, 24.1, 29.5, 13.8, 54.9, 10.7,
      46.7, 13.7, 95.6, 137.3, 18.4, 30.9, 6.3, 35.3, 8.5, 4.4),
  mi = c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0,
      1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0)
)

library(pROC)
# 1. ROC curve; direction = "<" means higher troponin -> more likely MI
r <- roc(dat$mi, dat$troponin, levels = c(0, 1), direction = "<")
auc(r)
ci.auc(r, method = "delong")               # DeLong 95% CI

# 2. Optimal cutoff by the Youden index (sensitivity + specificity - 1)
best <- coords(r, "best", best.method = "youden",
               ret = c("threshold", "sensitivity", "specificity", "ppv", "npv"))
best

# 3. Confusion table at that cutoff (positive test: troponin > threshold)
pos <- dat$troponin > best$threshold
table(Test = pos, MI = dat$mi)

# 4. Exact (Clopper-Pearson) 95% CIs for sensitivity and specificity
binom.test(sum(pos & dat$mi == 1), sum(dat$mi == 1))$conf.int
binom.test(sum(!pos & dat$mi == 0), sum(dat$mi == 0))$conf.int
# plot(r, print.thres = "best")
```

**What the code does**

- `roc(response, predictor, levels = c(0, 1), direction = "<")` builds the ROC curve. `levels` gives the control and case labels, and `direction = "<"` means controls have lower values. Left on its default (`"auto"`), pROC guesses the direction, which silently gives AUC < 0.5 when the guess is wrong.
- `auc()` and `ci.auc(method = "delong")` return the AUC and its DeLong 95% CI (the default for `ci.auc` with a single curve).
- `coords(r, "best", best.method = "youden", ret = ...)` finds the cutoff maximising sensitivity + specificity − 1 and returns the requested measures. pROC's threshold is a midpoint between two observed values, and a positive test is `value > threshold`.
- `table(Test = pos, MI = ...)` gives the 2 × 2 confusion table, so the numbers can be checked by hand.
- `binom.test(...)$conf.int` gives exact (Clopper–Pearson) CIs for sensitivity and specificity.

**Output**

```text
Area under the curve: 0.8177
95% CI: 0.7111-0.9244 (DeLong)
          threshold sensitivity specificity       ppv       npv
threshold      10.2        0.96         0.6 0.6315789 0.9545455
       MI
Test     0  1
  FALSE 21  1
  TRUE  14 24
[1] 0.7964831 0.9989878
attr(,"conf.level")
[1] 0.95
[1] 0.4211177 0.7612919
attr(,"conf.level")
[1] 0.95
```

### Python

```python
import numpy as np, pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve

# Troponin (ng/L) at presentation and final diagnosis of MI (1 = yes), 60 patients
dat = pd.DataFrame({
  'troponin': [9.7, 4.9, 22.2, 12.8, 5.9, 5.8, 6.1, 22.3, 3.7, 67.8, 22.9,
      20.3, 4.5, 44.9, 22.3, 8.8, 37.5, 11, 19.3, 27.9, 2.9, 5.2, 15.7, 4.7,
      34, 8.4, 19.3, 9.7, 3.9, 11.5, 38.4, 7.4, 21.2, 9.4, 23, 32.6, 19.8,
      12.4, 11.5, 12.9, 7.4, 20, 2.2, 48.3, 5.5, 24.1, 29.5, 13.8, 54.9, 10.7,
      46.7, 13.7, 95.6, 137.3, 18.4, 30.9, 6.3, 35.3, 8.5, 4.4],
  'mi': [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0,
      1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0]
})
y, x = dat["mi"].to_numpy(), dat["troponin"].to_numpy()

# 1. AUC (sklearn) and DeLong 95% CI (computed from structural components)
print("AUC = %.4f" % roc_auc_score(y, x))
pos, neg = x[y == 1], x[y == 0]
psi = lambda a, b: (a[:, None] > b[None, :]) + 0.5 * (a[:, None] == b[None, :])
v10 = psi(pos, neg).mean(axis=1)           # per-case placement values
v01 = psi(pos, neg).mean(axis=0)           # per-control placement values
auc = v10.mean()
se = np.sqrt(v10.var(ddof=1) / len(pos) + v01.var(ddof=1) / len(neg))
lo, hi = auc + np.array([-1, 1]) * stats.norm.ppf(0.975) * se
print("DeLong SE = %.4f, 95%% CI (%.4f, %.4f)" % (se, lo, hi))

# 2. Youden-optimal cutoff from roc_curve (rule: positive if troponin >= threshold)
fpr, tpr, thr = roc_curve(y, x)
j = np.argmax(tpr - fpr)                   # first maximum, like pROC
cut = thr[j]
lower = x[x < cut].max()                   # next lower observed value
print("Youden cutoff: x >= %.1f (pROC midpoint convention: > %.2f)" % (cut, (cut + lower) / 2))

# 3. Sensitivity, specificity, PPV, NPV at that cutoff
test = x >= cut
tp, fp = (test & (y == 1)).sum(), (test & (y == 0)).sum()
fn, tn = (~test & (y == 1)).sum(), (~test & (y == 0)).sum()
print(f"TP={tp} FP={fp} FN={fn} TN={tn}")
print("Sens=%.3f Spec=%.3f PPV=%.3f NPV=%.3f" % (tp/(tp+fn), tn/(tn+fp), tp/(tp+fp), tn/(tn+fn)))

# 4. Exact (Clopper-Pearson) 95% CIs for sensitivity and specificity
for name, k, n in [("Sens", tp, tp + fn), ("Spec", tn, tn + fp)]:
    c = stats.binomtest(int(k), int(n)).proportion_ci(confidence_level=0.95, method="exact")
    print(f"{name} 95% CI ({c.low:.3f}, {c.high:.3f})")
```

**What the code does**

- `roc_auc_score(y, x)` returns the AUC. `y` must be 0/1 with 1 the disease, and `x` the score (higher = more likely positive).
- scikit-learn has no CI for the AUC. The block computes the **DeLong** SE from the placement values (`v10`, `v01`) so it matches pROC; a bootstrap CI is an alternative.
- `roc_curve(y, x)` returns false- and true-positive rates and the thresholds (a positive test is `x >= threshold`). `argmax(tpr − fpr)` selects the Youden-optimal point. scikit-learn uses observed values as thresholds, while pROC uses the midpoint between observed values, so the number differs (10.7 versus 10.2) but the classification is identical.
- Sensitivity, specificity, PPV, and NPV are calculated from the 2 × 2 counts; `stats.binomtest(...).proportion_ci(method="exact")` gives Clopper–Pearson CIs (SciPy ≥ 1.7).

**Output**

```text
AUC = 0.8177
DeLong SE = 0.0544, 95% CI (0.7111, 0.9244)
Youden cutoff: x >= 10.7 (pROC midpoint convention: > 10.20)
TP=24 FP=14 FN=1 TN=21
Sens=0.960 Spec=0.600 PPV=0.632 NPV=0.955
Sens 95% CI (0.796, 0.999)
Spec 95% CI (0.421, 0.761)
```

### SAS

```sas
/* Troponin (ng/L) and final diagnosis of MI (1 = yes), 60 patients */
data trop;
  input troponin mi;
  datalines;
  9.7 0
  4.9 0
  22.2 1
  12.8 0
  5.9 0
  5.8 0
  6.1 0
  22.3 0
  3.7 0
  67.8 0
  22.9 1
  20.3 0
  4.5 0
  44.9 1
  22.3 0
  8.8 0
  37.5 1
  11 1
  19.3 0
  27.9 1
  2.9 0
  5.2 0
  15.7 1
  4.7 0
  34 0
  8.4 0
  19.3 0
  9.7 0
  3.9 0
  11.5 1
  38.4 1
  7.4 0
  21.2 0
  9.4 0
  23 1
  32.6 1
  19.8 1
  12.4 0
  11.5 0
  12.9 0
  7.4 0
  20 1
  2.2 0
  48.3 1
  5.5 0
  24.1 1
  29.5 0
  13.8 1
  54.9 0
  10.7 1
  46.7 1
  13.7 1
  95.6 1
  137.3 1
  18.4 1
  30.9 1
  6.3 1
  35.3 1
  8.5 0
  4.4 0
;
run;

/* 1. ROC curve, AUC and DeLong 95% CI.
      The ROC statement prints "ROC Association Statistics" (AUC, SE, 95% CI). */
ods graphics on;
proc logistic data=trop plots(only)=roc(id=prob);
  model mi(event='1') = troponin / outroc=rocdat;
  roc 'Troponin' troponin;
run;

/* 2. Youden index from the OUTROC data set: J = sensitivity + specificity - 1 */
data youden;
  set rocdat;
  j = _sensit_ - _1mspec_;
run;
proc sort data=youden; by descending j; run;
proc print data=youden(obs=1);
  var _prob_ _sensit_ _1mspec_ j;
run;

/* 3. Sens, spec, PPV, NPV at the chosen cutoff (troponin > 10.2 ng/L, from R) */
data cut;
  set trop;
  test = (troponin > 10.2);
run;
/* 2x2 table: Sens = TRUE row among MI=1, Spec = FALSE row among MI=0,
   PPV = MI=1 among test TRUE, NPV = MI=0 among test FALSE */
proc freq data=cut;
  tables test*mi / nopercent;
run;

/* Exact (Clopper-Pearson) 95% CIs: sensitivity, then specificity */
proc freq data=cut;
  where mi = 1;
  tables test / binomial(level='1' exact);
run;
proc freq data=cut;
  where mi = 0;
  tables test / binomial(level='0' exact);
run;
```

**What the code does**

- `model mi(event='1') = troponin` fits a logistic model with troponin as the only predictor and models P(MI = 1). `outroc=rocdat` saves the ROC points (`_PROB_`, `_SENSIT_`, `_1MSPEC_`).
- The `roc 'Troponin' troponin;` statement tells PROC LOGISTIC to compute the ROC curve and AUC for troponin directly. It prints the **ROC Association Statistics** table: AUC, its **DeLong** standard error, and the 95% Wald CI. The table may list the fitted model's own curve ('Model') as well as 'Troponin'; they are the same curve here, and the AUC should match R (0.818). Check your SAS version's documentation for the exact table layout. `plots(only)=roc` draws the curve.
- Because logistic regression with one predictor is a monotone transformation of troponin, the AUC is the same as for the raw marker. Cutoffs in the OUTROC data are on the predicted-probability scale; convert them back to troponin units (the cutoff above is the one found in R).
- The `youden` data step computes J = `_SENSIT_` − `_1MSPEC_` and sorts to display the optimum.
- `proc freq` tabulates the 2 × 2 table at the chosen cutoff. `binomial(level=... exact)` gives exact CIs for sensitivity (restricted to MI = 1) and for specificity (restricted to MI = 0).
- To compare two markers use several `roc` statements and `roccontrast` (DeLong test for correlated AUCs).

## Interpreting the output

- **AUC = 0.818 (95% CI 0.711 to 0.924).** There is an 82% chance that a random MI patient has a higher troponin than a random non-MI patient. This is good, not excellent, discrimination. The CI excludes 0.5, so the test discriminates better than chance, but the lower limit (0.71) shows the uncertainty with only 25 cases. R and Python give the same AUC and DeLong CI.
- **Youden cutoff = troponin > 10.2 ng/L** (Python reports x ≥ 10.7, which classifies patients identically). This cutoff maximises sensitivity + specificity − 1 (J = 0.56).
- **Sensitivity = 0.96 (24 of 25; exact 95% CI 0.80 to 1.00).** The test catches nearly all MI cases, with one false negative.
- **Specificity = 0.60 (21 of 35; 95% CI 0.42 to 0.76).** Of the 35 patients without MI, 14 are false positives. The data are compatible with specificity anywhere from about 40% to 75%.
- **PPV = 0.63 (24 of 38) and NPV = 0.95 (21 of 22).** A positive result means about a 63% chance of MI, a negative result about a 95% chance of no MI. These depend on the sample prevalence (25 of 60 = 42%), which is higher than in most chest-pain populations, so PPV would be lower and NPV higher in a lower-risk population.
- The cutoff was chosen **using the same data**, so these accuracy estimates are optimistic (the Youden cutoff is also unstable in samples this small, and it weighs false positives and false negatives equally, which may not match the clinical aim). Validate in a new sample or with a bootstrap. Clinically, the cutoff is often chosen to guarantee a high sensitivity (e.g., 95% for rule-out) rather than to maximise Youden's J.

**How to report it:**
> Troponin discriminated MI with an AUC of 0.82 (95% CI 0.71 to 0.92; DeLong). At the Youden-optimal cutoff of 10.2 ng/L, sensitivity was 96% (95% CI 80 to 100%), specificity 60% (95% CI 42 to 76%), PPV 63%, and NPV 95% (prevalence 42%, n = 60).

## Check your understanding

<details class="check"><summary>1. The troponin test has sensitivity 96% and specificity 60% in this sample, where prevalence of MI is 42%. In a general practice setting prevalence is 5%. Will the PPV still be 63%?</summary>

No. PPV and NPV depend on prevalence, whereas sensitivity and specificity do not (to a first approximation). At 5% prevalence, with sensitivity 0.96 and specificity 0.60, PPV = (0.96 × 0.05) / (0.96 × 0.05 + 0.40 × 0.95) = 0.048 / 0.428 = about 11%. Most positive results would then be false positives, while NPV would be above 99%.
</details>

<details class="check"><summary>2. A paper reports "the optimal cutoff of 10.2 ng/L had sensitivity 96% and specificity 60%." The cutoff was chosen from the same 60 patients. What is the concern?</summary>

Choosing the cutoff that gives the best performance in the same data capitalises on chance, so the reported sensitivity and specificity are optimistic. The Youden cutoff is also an unstable choice in small samples (a few patients can move it). Validate in a new sample or use a bootstrap or cross-validation optimism correction. And note that Youden's J weighs false positives and false negatives equally, which may not match the clinical aim: a rule-out test usually fixes a high sensitivity (for example 95%+) instead.
</details>

<details class="check"><summary>3. Two markers have AUC 0.82 and 0.79 on the same 60 patients. Can you declare the first better because 0.82 > 0.79?</summary>

Not from the point estimates alone. The AUCs come from the same patients and are correlated, so compare them with the DeLong test for correlated AUCs and look at the CI for the difference. With only 25 cases, an AUC CI of roughly ±0.1 means a 0.03 gap could easily be noise (the CI for the paired difference is narrower than for each AUC, so test the difference directly).
</details>

<details class="check"><summary>4. What would happen to the AUC if all troponin values were multiplied by 10 (a different unit)?</summary>

Nothing. The AUC depends only on the ranks of the scores (it is the Mann–Whitney U divided by n₁n₂). Any monotone transformation of the marker, such as multiplying by 10 or taking logs, leaves AUC and the ROC curve unchanged. The cutoff value changes, of course.
</details>

## Common mistakes

- **Choosing the cutoff and reporting its performance in the same data**, without validation or bootstrap correction.
- **Quoting PPV/NPV without prevalence**, or applying them to a population with different prevalence.
- **Wrong direction.** Letting software auto-detect the direction can give AUC < 0.5 or a flipped curve; set it explicitly.
- **Equating a high AUC with a well-calibrated or clinically useful model.** AUC ignores calibration and the consequences of errors; consider decision-curve analysis.
- **Comparing AUCs by eye or with independent tests.** AUCs from the same patients are correlated; use the DeLong test for paired curves.
- **Dichotomising too early.** Reporting only sensitivity and specificity at a single arbitrary cutoff discards the curve.
- **Reporting AUC with no CI**, especially in small studies, where the CI is wide.
- **Treating ordinal categories as continuous** without thinking about ties; ties are counted as one half in the AUC.
- **Overfitting the model** that produced the score (too few events per variable); use cross-validation or an external validation set.

## Related tests

- [Logistic regression](#/test/logistic-regression): builds the predicted probability that the ROC curve evaluates
- [Mann–Whitney U test](#/test/mann-whitney): the AUC is the Mann–Whitney U divided by n₁·n₂
- [McNemar's test](#/test/mcnemar): comparing two binary tests on the same patients
- [One-sample proportion (exact binomial) test](#/test/one-proportion): exact CIs for sensitivity and specificity
- [Cohen's kappa](#/test/cohens-kappa): agreement between two classifications
- [Chi-square test of independence](#/test/chi-square-independence): association of the 2 × 2 table

## Key takeaways

- ROC curve: sensitivity vs 1 − specificity across all cutoffs. AUC = probability a random case scores higher than a random control (0.5 = chance).
- Sensitivity and specificity are properties of the test. PPV and NPV also depend on **prevalence**.
- The Youden cutoff is only one way to choose, and it weights errors equally. A cutoff picked in the same data is optimistic: validate it.
- Report the AUC with a CI, plus sensitivity, specificity, PPV, and NPV with CIs at the chosen cutoff.
- AUC measures discrimination, not calibration or clinical usefulness. Compare correlated AUCs with DeLong's test.

## References

- Robin X, Turck N, Hainard A, et al. pROC: an open-source package for R and S+ to analyze and compare ROC curves. *BMC Bioinformatics.* 2011;12:77.
- DeLong ER, DeLong DM, Clarke-Pearson DL. Comparing the areas under two or more correlated receiver operating characteristic curves: a nonparametric approach. *Biometrics.* 1988;44:837–845.
- Bossuyt PM, Reitsma JB, Bruns DE, et al. STARD 2015: an updated list of essential items for reporting diagnostic accuracy studies. *BMJ.* 2015;351:h5527.
- R: `?pROC::roc`, `?pROC::coords`; scikit-learn: `sklearn.metrics.roc_auc_score`, `roc_curve`; SAS/STAT User's Guide: The LOGISTIC Procedure (ROC, ROCCONTRAST statements)
