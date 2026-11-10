# Paired t-test

> Compares two measurements taken on the **same subjects** — before vs after, left eye vs right eye, drug A vs drug B in a crossover — by testing whether the mean of the within-subject differences is zero.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Two measurements on the same subjects (before/after, crossover, matched pairs) | Continuous outcome | Pairs independent of each other; differences roughly normal | Systolic BP before vs after 8 weeks of a new drug |

## Learning objectives

- Recognize when data are **paired** and why that matters
- Explain the paired t-test as a one-sample test on the differences
- Run it in R and Python and read the mean difference, 95% CI, t, and p
- Report the result correctly and avoid the classic paired-data mistakes

## What it is

The paired t-test reduces each subject's two measurements to a single **difference** (e.g., `after − before`) and then runs a one-sample t-test on those differences against 0. Because each subject serves as their own control, between-person variability cancels out, which usually makes this test far more powerful than comparing two independent groups.

- **Test statistic:** t = d̄ / (s_d / √n), where d̄ is the mean difference, s_d its standard deviation, and n the number of pairs
- **Degrees of freedom:** n − 1
- **Effect size:** mean difference with 95% CI (primary), Cohen's d_z = d̄ / s_d (secondary)

## The intuition

Every test statistic is a **signal-to-noise ratio**. Here the signal is the average change, d̄. The noise is how much that average would bounce around from sample to sample by chance, its standard error s_d/√n. If the average change is large compared with that chance wobble, the change is unlikely to be a fluke.

Pairing is what makes the noise small. People differ a lot from each other: one patient's baseline SBP is 145 and another's is 162. That between-person spread is irrelevant when you look at each patient's *own* change. Compare the two pictures below. The raw before and after values overlap heavily, but the differences are all clearly below zero.

<div class="viz" data-viz="dots" data-unit="mmHg" data-groups='{"Baseline":[152,148,160,145,158,150,162,155,149,157],"Week 8":[144,145,151,140,150,147,152,149,146,148]}' data-caption="Raw values: the two groups overlap heavily because patients differ from each other."></div>

<div class="viz" data-viz="dots" data-unit="mmHg, week 8 − baseline" data-groups='{"Change":[-8,-3,-9,-5,-8,-3,-10,-6,-3,-9]}' data-zero="1" data-caption="Within-patient changes: every single patient went down."></div>

Our example gives t = −7.34 with 9 degrees of freedom. If the drug had no effect, t values would follow the curve below. Ours is far out in the tail, so the shaded p-value area is tiny.

<div class="viz" data-viz="t" data-stat="-7.34" data-df="9"></div>

## When to use it

Use the paired t-test when **all** of these are true:

- The outcome is **continuous** (blood pressure, HbA1c, pain on a 0–100 VAS, tumor size)
- There are **exactly two** measurements per subject, or two subjects that are deliberately matched (matched case–control pairs, twins)
- The **differences** are approximately normally distributed (the raw values do not need to be)

**Clinical scenario.** A pilot study enrols 10 patients with stage 1 hypertension. Systolic blood pressure (SBP, mmHg) is measured at baseline and again after 8 weeks on a new antihypertensive. *Did SBP change, and by how much?* Each patient contributes a before and an after value, so the observations are paired.

**Don't use it when:**

- The two groups contain **different people** → use the [independent t-test](#/test/independent-t)
- There are **three or more** time points → use [repeated-measures ANOVA](#/test/repeated-measures-anova) or a [mixed model](#/test/mixed-models)
- The differences are clearly skewed or have outliers with a small n → use the [Wilcoxon signed-rank test](#/test/wilcoxon-signed-rank)
- The outcome is binary (responder yes/no) → use [McNemar's test](#/test/mcnemar)

## Assumptions

1. **Pairs are independent of each other.** Patient 3's difference tells you nothing about patient 7's.
2. **Differences are approximately normal.** Check with a histogram or Q–Q plot of the differences, or a [Shapiro–Wilk test](#/test/shapiro-wilk). With n ≥ 30 pairs, the t-test is robust to moderate non-normality.
3. **The outcome is measured on an interval or ratio scale.**

## Hypotheses

- **H₀:** μ_d = 0 (the mean within-subject change is zero)
- **H₁:** μ_d ≠ 0 (two-sided; the default and almost always what regulators expect)

## Example data

| Patient | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| Baseline SBP | 152 | 148 | 160 | 145 | 158 | 150 | 162 | 155 | 149 | 157 |
| Week 8 SBP | 144 | 145 | 151 | 140 | 150 | 147 | 152 | 149 | 146 | 148 |
| Difference | −8 | −3 | −9 | −5 | −8 | −3 | −10 | −6 | −3 | −9 |

## How to run it

### R

```r
# Systolic BP (mmHg) for 10 patients, same order in both vectors
before <- c(152, 148, 160, 145, 158, 150, 162, 155, 149, 157)
after  <- c(144, 145, 151, 140, 150, 147, 152, 149, 146, 148)

# 1. Check the normality of the DIFFERENCES (not the raw values)
d <- after - before
shapiro.test(d)

# 2. Run the paired t-test
t.test(after, before, paired = TRUE)

# 3. Effect size: Cohen's d_z
mean(d) / sd(d)
```

**What the code does**

- `before` and `after` must be in the **same patient order** — R pairs them by position.
- `d <- after - before` computes each patient's change; negative means SBP went down.
- `shapiro.test(d)` checks normality of the differences, which is the actual assumption.
- `t.test(..., paired = TRUE)` is the key argument. Without it, R runs a Welch independent-samples test and gives a very different (wrong) answer.
- The order `after, before` means the reported difference is after − before.

**Output**

```text
	Shapiro-Wilk normality test

data:  d
W = 0.87233, p-value = 0.1064

	Paired t-test

data:  after and before
t = -7.3413, df = 9, p-value = 4.368e-05
alternative hypothesis: true mean difference is not equal to 0
95 percent confidence interval:
 -8.372103 -4.427897
sample estimates:
mean difference
           -6.4

> mean(d) / sd(d)
[1] -2.321524
```

### Python

```python
import numpy as np
from scipy import stats

# Systolic BP (mmHg) for 10 patients, same order in both arrays
before = np.array([152, 148, 160, 145, 158, 150, 162, 155, 149, 157])
after  = np.array([144, 145, 151, 140, 150, 147, 152, 149, 146, 148])

# 1. Normality of the differences
d = after - before
sw = stats.shapiro(d)
print(f"Shapiro-Wilk W = {sw.statistic:.4f}, p = {sw.pvalue:.4f}")

# 2. Paired t-test
res = stats.ttest_rel(after, before)
ci = res.confidence_interval(confidence_level=0.95)
print(f"t = {res.statistic:.3f}, df = {res.df}, p = {res.pvalue:.2e}")
print(f"Mean difference = {d.mean():.2f}, 95% CI ({ci.low:.2f}, {ci.high:.2f})")

# 3. Effect size: Cohen's d_z (ddof=1 gives the sample SD)
print(f"Cohen's d_z = {d.mean() / d.std(ddof=1):.2f}")
```

**What the code does**

- `stats.ttest_rel` is SciPy's **related-samples** (paired) t-test. `ttest_ind` is the independent version — a common mix-up.
- `res.confidence_interval()` returns the 95% CI for the mean difference (SciPy ≥ 1.10).
- `d.std(ddof=1)` gives the sample standard deviation. NumPy's default `ddof=0` divides by n, not n − 1, and would understate the SD.

**Output**

```text
Shapiro-Wilk W = 0.8723, p = 0.1064
t = -7.341, df = 9, p = 4.37e-05
Mean difference = -6.40, 95% CI (-8.37, -4.43)
Cohen's d_z = -2.32
```

### SAS

```sas
/* One row per patient, before and after side by side */
data bp;
  input patient before after;
  datalines;
1 152 144
2 148 145
3 160 151
4 145 140
5 158 150
6 150 147
7 162 152
8 155 149
9 149 146
10 157 148
;
run;

/* Paired t-test: the PAIRED statement computes after - before */
proc ttest data=bp alpha=0.05;
  paired after*before;
run;

/* Optional: normality of the differences */
data bp_diff;
  set bp;
  diff = after - before;
run;

proc univariate data=bp_diff normal;
  var diff;
  qqplot diff / normal(mu=est sigma=est);
run;
```

**What the code does**

- The data are in **wide format**: one row per patient, with both measurements on the same row.
- `paired after*before;` tells PROC TTEST to analyse the difference `after − before`. Reversing the order flips the sign of the difference.
- PROC TTEST prints the mean difference, its 95% CI, t, DF, and Pr > |t|, along with a histogram and Q–Q plot of the differences when ODS Graphics is on.
- `proc univariate ... normal` adds Shapiro–Wilk and other normality tests on `diff`.

## Interpreting the output

- **Mean difference = −6.4 mmHg.** On average, SBP fell by 6.4 mmHg after 8 weeks.
- **95% CI: −8.4 to −4.4 mmHg.** The data are consistent with a true mean reduction of roughly 4 to 8 mmHg. The CI excludes 0, which agrees with the p-value.
- **t(9) = −7.34, p < 0.001.** If the drug truly had no effect, a mean change this large would be very unlikely.
- **Shapiro–Wilk on the differences: p = 0.11.** There is no strong evidence against normality. With only 10 pairs, also look at a Q–Q plot.
- **Cohen's d_z = −2.3.** This is a very large standardized effect. For clinical reporting, the mean difference in mmHg is more meaningful.

**How to report it:**
> Systolic blood pressure decreased from baseline to week 8 (mean change −6.4 mmHg, 95% CI −8.4 to −4.4; paired t(9) = −7.34, p < 0.001).

⚠️ **Important for trials:** a before/after change in a *single arm* does not prove the drug caused it. Regression to the mean, placebo effects, and time trends can all produce a drop. To estimate a causal treatment effect, you need a control group, typically analysed with [ANCOVA](#/test/ancova).

## Check your understanding

<details class="check"><summary>1. Twenty patients each have their left and right eye treated with a different eye drop. Which test fits: independent or paired?</summary>

Paired. Each patient contributes both eyes, so the two measurements are linked. Analyse the within-patient differences.
</details>

<details class="check"><summary>2. The Shapiro–Wilk test on the raw "before" values gives p = 0.01. Does that rule out the paired t-test?</summary>

No. The paired t-test assumes the **differences** are roughly normal, not the raw values. Check the differences instead, ideally with a Q–Q plot.
</details>

<details class="check"><summary>3. A paired t-test gives p = 0.21. A colleague writes "the drug has no effect on blood pressure." What's wrong?</summary>

A non-significant p-value is not evidence of no effect. Look at the 95% CI. If it runs from, say, −9 to +2 mmHg, a clinically important reduction is still compatible with the data. To claim "no meaningful difference", you need an [equivalence design](#/test/tost-equivalence).
</details>

<details class="check"><summary>4. Suppose the average change stayed at −6.4 mmHg, but the changes were twice as spread out (SD 5.5 instead of 2.8 mmHg). What happens to t?</summary>

The signal is unchanged, but s_d doubles, so the standard error doubles and t halves: −6.4 / (5.51/√10) ≈ −3.7. The noise grew, so the signal-to-noise ratio shrank, and p gets larger.
</details>

## Common mistakes

- **Running an independent t-test on paired data.** This ignores the pairing, inflates the variance, and loses power. In R, check that `paired = TRUE` is set.
- **Testing normality of the raw values instead of the differences.** Only the differences need to be normal.
- **Mismatched order.** If the rows of `before` and `after` are not the same patients in the same order, the results are meaningless.
- **Dropping patients with one missing value without saying so.** The test uses only complete pairs. Report how many pairs were analysed and consider a [mixed model](#/test/mixed-models) if missingness is substantial.
- **Reporting only the p-value.** Always give the mean difference and its 95% CI.

## Related tests

- [Wilcoxon signed-rank test](#/test/wilcoxon-signed-rank): nonparametric alternative
- [One-sample t-test](#/test/one-sample-t): the paired t-test is a one-sample t-test on the differences
- [Repeated-measures ANOVA](#/test/repeated-measures-anova): three or more time points
- [McNemar's test](#/test/mcnemar): paired binary outcomes
- [Bland–Altman analysis](#/test/bland-altman): when the question is *agreement* between two methods, not a change

## Key takeaways

- Paired data mean **the same subjects measured twice** (or matched pairs). Analyse the differences.
- The paired t-test is a one-sample t-test on those differences: t = d̄ / (s_d/√n), df = n − 1.
- Only the **differences** need to be roughly normal. If they are skewed and n is small, use the Wilcoxon signed-rank test.
- Report the mean difference with its 95% CI, not just p.
- A before/after change in one group is not proof that the treatment caused it. You need a control group.

## References

- Altman DG. *Practical Statistics for Medical Research.* Chapman & Hall; 1991. Ch. 9.
- R documentation: `?t.test`
- SciPy documentation: `scipy.stats.ttest_rel`
- SAS/STAT User's Guide: The TTEST Procedure — PAIRED statement
