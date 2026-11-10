# One-sample t-test

> Tests whether the mean of a single group differs from a fixed reference value, such as a treatment target, a historical control, or a normal-range midpoint.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| One group compared with a fixed reference value set in advance | Continuous outcome | Independent observations; data roughly normal; μ₀ pre-specified | Mean HbA1c of 12 patients vs the 7.0% target |

## Learning objectives

- Recognize when a single group is being compared with a fixed reference value
- Explain t as the distance of the sample mean from μ₀ in standard errors
- Run the test in R and Python and read the mean, 95% CI, t, and p
- Interpret a non-significant result without claiming the mean equals the target

## What it is

The one-sample t-test compares the sample mean with a number you specify **in advance**. It measures how many standard errors the sample mean lies from that reference value. Because the population SD is unknown, it is estimated from the data, and the t distribution accounts for that extra uncertainty.

- **Test statistic:** t = (x̄ − μ₀) / (s / √n), where x̄ is the sample mean, μ₀ the reference value, s the sample SD, and n the sample size
- **Degrees of freedom:** n − 1
- **Effect size:** the mean (or mean difference from μ₀) with its 95% CI (primary); Cohen's d = (x̄ − μ₀) / s (secondary)

## The intuition

The test statistic is a **signal-to-noise ratio**. The signal is how far the sample mean sits from the reference value, x̄ − μ₀. The noise is the standard error s/√n: how much a mean of 12 patients would wobble by chance alone. A mean that lands 1 standard error from the target is unremarkable; one that lands 4 standard errors away is not.

An everyday analogy: a bathroom scale reads 71.0 kg, and your doctor's records say 70.0 kg. Is that real weight gain? It depends on how much the scale jitters from reading to reading. A 1 kg gap is huge for a precise scale and meaningless for a noisy one. The standard error plays the role of the jitter.

The plot below shows each patient's HbA1c minus the 7.0% target. The dots straddle zero: six patients are below target and six above, and the mean sits only 0.125 points to the right of the line, well within the scatter of the individual values.

<div class="viz" data-viz="dots" data-unit="HbA1c %, patient minus 7.0 target" data-groups='{"HbA1c minus 7.0":[0.4,0.4,0.3,-0.1,0.3,-0.2,-0.3,-0.2,0.8,-0.5,-0.1,0.7]}' data-zero="1" data-caption="Patients scatter on both sides of the target. The mean is only a little to the right of zero relative to the spread."></div>

Our example gives t = 1.05 with 11 degrees of freedom. If the true mean were exactly 7.0%, t values this large (in either direction) would happen about a third of the time.

<div class="viz" data-viz="t" data-stat="1.05" data-df="11"></div>

**What to notice:** the observed t is near the middle of the curve, and the shaded tail area (p = 0.32) is large. The data are not surprising under H₀.

## When to use it

Use the one-sample t-test when **all** of these are true:

- The outcome is **continuous** (HbA1c, blood pressure, hemoglobin, tumor size)
- There is **one group** of independent subjects
- You have a **meaningful reference value** μ₀ fixed before looking at the data
- The data are approximately normal, or n is reasonably large

**Clinical scenario.** After 6 months of a lifestyle programme, HbA1c (%) is measured in 12 patients with type 2 diabetes. The clinical target is a mean of 7.0%. *Is the mean HbA1c in this group different from 7.0%?*

**Don't use it when:**

- You are comparing **two groups** → use the [independent t-test](#/test/independent-t)
- The data are **before/after on the same people** → use the [paired t-test](#/test/paired-t)
- The data are clearly skewed or have outliers with a small n → use the [Wilcoxon signed-rank test](#/test/wilcoxon-signed-rank)
- The outcome is a binary event rate compared with a reference → use the [one-sample proportion test](#/test/one-proportion)

## Assumptions

1. **Observations are independent.** One value per patient; no repeated measurements pooled together.
2. **The data are approximately normal.** Check with a Q–Q plot or the [Shapiro–Wilk test](#/test/shapiro-wilk). With n ≥ 30 the test is robust to moderate non-normality.
3. **The outcome is on an interval or ratio scale.**
4. **μ₀ was chosen before analysis**, not after seeing the mean.

## Hypotheses

- **H₀:** μ = μ₀ (the population mean equals the reference value, here 7.0%)
- **H₁:** μ ≠ μ₀ (two-sided)

## Example data

HbA1c (%) in 12 patients at month 6:

| Patient | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HbA1c (%) | 7.4 | 7.4 | 7.3 | 6.9 | 7.3 | 6.8 | 6.7 | 6.8 | 7.8 | 6.5 | 6.9 | 7.7 |

## How to run it

### R

```r
# HbA1c (%) in 12 patients after 6 months of lifestyle intervention
a1c <- c(7.4, 7.4, 7.3, 6.9, 7.3, 6.8, 6.7, 6.8, 7.8, 6.5, 6.9, 7.7)

# Normality check
shapiro.test(a1c)

# One-sample t-test against the 7.0% treatment target
t.test(a1c, mu = 7.0)

# Effect size: Cohen's d
(mean(a1c) - 7.0) / sd(a1c)
```

**What the code does**

- `t.test(a1c, mu = 7.0)` tests the mean against 7.0. If you leave out `mu`, R silently tests against **0**, which is rarely what you want.
- The default is two-sided with a 95% CI. Use `alternative = "less"` or `"greater"` only for a pre-specified one-sided hypothesis.
- The reported CI is for the **mean** (not the difference from μ₀). Subtract 7.0 from its limits to get the CI for the difference.
- `sd(a1c)` uses n − 1 in the denominator.

**Output**

```text
	Shapiro-Wilk normality test

data:  a1c
W = 0.9388, p-value = 0.4827


	One Sample t-test

data:  a1c
t = 1.0467, df = 11, p-value = 0.3177
alternative hypothesis: true mean is not equal to 7
95 percent confidence interval:
 6.862156 7.387844
sample estimates:
mean of x 
    7.125 

[1] 0.3021613
```

### Python

```python
import numpy as np
from scipy import stats

# HbA1c (%) in 12 patients after 6 months of lifestyle intervention
a1c = np.array([7.4, 7.4, 7.3, 6.9, 7.3, 6.8, 6.7, 6.8, 7.8, 6.5, 6.9, 7.7])

# Normality check
sw = stats.shapiro(a1c)
print(f"Shapiro-Wilk W = {sw.statistic:.4f}, p = {sw.pvalue:.4f}")

# One-sample t-test against the 7.0% target
res = stats.ttest_1samp(a1c, popmean=7.0)
ci = res.confidence_interval(confidence_level=0.95)
print(f"t = {res.statistic:.3f}, df = {res.df}, p = {res.pvalue:.4f}")
print(f"Mean = {a1c.mean():.3f}, 95% CI ({ci.low:.3f}, {ci.high:.3f})")

# Effect size: Cohen's d (sample SD)
print(f"Cohen's d = {(a1c.mean() - 7.0) / a1c.std(ddof=1):.2f}")
```

**What the code does**

- `stats.ttest_1samp(a1c, popmean=7.0)` needs the reference value explicitly. It is two-sided by default.
- `res.confidence_interval()` (SciPy ≥ 1.10) returns the CI for the **mean**, matching R's output.
- `a1c.std(ddof=1)` is the sample SD. NumPy's default `ddof=0` divides by n and would overstate d.

**Output**

```text
Shapiro-Wilk W = 0.9388, p = 0.4827
t = 1.047, df = 11, p = 0.3177
Mean = 7.125, 95% CI (6.862, 7.388)
Cohen's d = 0.30
```

### SAS

```sas
data a1c;
  input patient a1c @@;
  datalines;
1 7.4 2 7.4 3 7.3 4 6.9 5 7.3 6 6.8
7 6.7 8 6.8 9 7.8 10 6.5 11 6.9 12 7.7
;
run;

/* H0=7 sets the reference value; default is 0 */
proc ttest data=a1c h0=7 alpha=0.05 sides=2;
  var a1c;
run;

/* Optional: Shapiro-Wilk normality test */
proc univariate data=a1c normal;
  var a1c;
run;
```

**What the code does**

- Each value is read as a patient number and an HbA1c value. `@@` lets several pairs sit on one `datalines` line.
- `h0=7` sets μ₀. As in R, leaving it out tests against 0.
- `sides=2` and `alpha=0.05` are the defaults; they are shown for clarity. `sides=L` or `sides=U` gives one-sided tests.
- PROC TTEST prints the mean, SD, 95% CI for the mean, t, DF, and Pr > |t|.
- `proc univariate ... normal` adds the Shapiro–Wilk test.

## Interpreting the output

- **Mean HbA1c = 7.125%.** This is 0.125 percentage points above the 7.0% target.
- **95% CI: 6.86% to 7.39%.** The data are compatible with a true mean anywhere from slightly below to slightly above target. The CI contains 7.0, which agrees with the p-value.
- **t(11) = 1.05, p = 0.32.** The data do not show a difference from 7.0%. This is not evidence of no difference: the mean could plausibly be anywhere in the CI, which is fairly wide with 12 patients.
- **Shapiro–Wilk p = 0.48.** No evidence against normality.
- **Cohen's d = 0.30.** A small standardized difference by Cohen's rough conventions, which are not universal benchmarks. The HbA1c units are more meaningful clinically.

**How to report it:**
> Mean HbA1c at month 6 was 7.13% (95% CI 6.86 to 7.39), which did not differ significantly from the 7.0% target (one-sample t(11) = 1.05, p = 0.32).

## Check your understanding

<details class="check"><summary>1. A trial reports that mean LDL in 40 patients "differs significantly from 0 mg/dL" (p < 0.001). What is the likely error?</summary>

The software was left at its default reference value of 0 (`mu` in R, `popmean` in SciPy, `h0` in SAS). Almost any positive outcome is "significantly different from zero". The reference value must be a meaningful, pre-specified target.
</details>

<details class="check"><summary>2. Here p = 0.32. Can we say the programme brought mean HbA1c to 7.0%?</summary>

No. p > 0.05 means the data are not inconsistent with 7.0%, not that the mean equals 7.0%. The 95% CI runs from 6.86 to 7.39, so true means a little above or below target are compatible with the data. To show closeness to a target, use an [equivalence (TOST) design](#/test/tost-equivalence).
</details>

<details class="check"><summary>3. Suppose the same mean (7.125) and SD (0.41) had come from 48 patients instead of 12. What happens to t?</summary>

The signal is unchanged but the standard error halves (√48 is twice √12), so t roughly doubles to about 2.1 and the result becomes significant at 0.05. The CI also narrows. This shows that "significance" depends on sample size as well as on the size of the difference.
</details>

<details class="check"><summary>4. The Shapiro–Wilk test gives p = 0.48 with n = 12. Does that prove the data are normal?</summary>

No. With 12 observations the test has little power to detect non-normality. Use it with a Q–Q plot and with subject-matter knowledge. For a one-sample test the raw values are what should be roughly normal.
</details>

## Common mistakes

- **Forgetting `mu` / `popmean` / `h0`.** The software then tests against 0, which gives a meaningless "significant" result for almost any positive outcome.
- **Picking μ₀ after seeing the data.** The reference value must come from the protocol or the literature.
- **Concluding "equal to target" from p > 0.05.** To show a mean is close to a target, use an equivalence approach ([TOST](#/test/tost-equivalence)).
- **Using it on paired data.** Two measurements per patient need a [paired t-test](#/test/paired-t).
- **Ignoring skewness in small samples.** Check the data before relying on t.

## Related tests

- [Paired t-test](#/test/paired-t): a one-sample t-test on the within-patient differences
- [Independent two-sample t-test](#/test/independent-t): two separate groups
- [Wilcoxon signed-rank test](#/test/wilcoxon-signed-rank): nonparametric alternative
- [Shapiro–Wilk test](#/test/shapiro-wilk): check normality
- [Equivalence (TOST) tests](#/test/tost-equivalence): when you want to show closeness to a target

## Key takeaways

- t = (x̄ − μ₀)/(s/√n): the distance from the reference value measured in standard errors, with df = n − 1.
- μ₀ must be fixed **before** looking at the data, and you must actually supply it (the default is 0).
- p > 0.05 means "not inconsistent with μ₀", never "equal to μ₀". Look at the CI.
- A larger n shrinks the standard error, so the same difference becomes easier to detect.
- Report the mean and its 95% CI, not just p.

## References

- Altman DG. *Practical Statistics for Medical Research.* Chapman & Hall; 1991. Ch. 9.
- Student. The probable error of a mean. *Biometrika.* 1908;6:1–25.
- R documentation: `?t.test`; SciPy documentation: `scipy.stats.ttest_1samp`
- SAS/STAT User's Guide: The TTEST Procedure — H0= option
