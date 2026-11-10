# One-sample z-test

> Compares a sample mean with a reference value when the population standard deviation σ is **truly known**, which is rare; when σ must be estimated from the data, use the [one-sample t-test](#/test/one-sample-t) instead.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| One sample vs a reference value when σ is genuinely known | Continuous outcome | σ known from outside the sample; independent; mean approximately normal | Mean height-for-age z-score vs the WHO reference of 0 |

## Learning objectives

- Explain what makes a z-test different from a t-test (a known σ) and why "n ≥ 30" is not the rule
- Compute z, the p-value, and a 95% CI by hand in R and Python
- Spot the statsmodels `ztest` pitfall, which silently uses the sample SD
- Report a result as a mean difference with its CI, not only a p-value

## What it is

The one-sample z-test measures how many **standard errors** the sample mean lies from a reference value μ₀, where the standard error is built from a σ that is known in advance and not estimated from your sample. Under H₀ the statistic follows the standard normal distribution exactly (if the data are normal) or approximately (large n, by the central limit theorem).

- **Test statistic:** z = (x̄ − μ₀) / (σ / √n), where x̄ is the sample mean, μ₀ the reference value, σ the known population SD, and n the sample size
- **Reference distribution:** standard normal N(0, 1); there are no degrees of freedom
- **Effect size:** the mean difference x̄ − μ₀ with its 95% CI, x̄ ± 1.96 σ/√n (primary); the standardized difference (x̄ − μ₀) / σ (secondary)

**An honest warning about "use z when n ≥ 30".** Many textbooks say to use a z-test for large samples and a t-test for small ones. That rule is a historical shortcut from the era of printed tables, when the normal table was short and the t table was long. It is not a statistical principle. If you **estimate** the SD from your data, the t-test is always the correct choice, whatever n is. For large n the t distribution becomes practically identical to the normal, so the two tests agree anyway and you lose nothing by using t. A z-test is justified only when σ comes from **outside** your sample.

## The intuition

A z-score is a **ruler measured in standard errors**. The signal is how far the sample mean sits from the reference value (x̄ − μ₀). The noise is the standard error σ/√n: how much a mean of 60 children would wobble from sample to sample purely by chance. Think of a bathroom scale whose accuracy is guaranteed by the manufacturer. If you weigh yourself and the reading is 3 guaranteed-error-widths off from your usual weight, you believe the change. With a known σ, the "error width" is a fact, not a guess, so the ruler is exactly a standard normal curve.

The dots below are the 60 HAZ values (split into two halves so each half is easy to read). The dashed no-difference line is the reference mean 0. Most children sit below it, and the mean of each half is below 0, even though individual children scatter widely on both sides.

<div class="viz" data-viz="dots" data-unit="HAZ" data-groups='{"Children 1-30": [0.63, 0.12, -0.46, -0.56, 0.81, 0.94, 0.18, -0.48, -1.57, -1.47, -2.02, 0.12, 0.48, -0.06, -3.62, -0.5, -1.28, 0.1, -1.96, -2.36, 0.79, -0.03, 0.69, -0.56, 0.05, -1.89, -0.47, 0.33, -2.41, -0.46], "Children 31-60": [-1.79, -0.68, -0.86, 0.85, -0.85, 0.04, -0.01, -0.84, -0.92, 0.08, -1.22, -0.35, 0.15, 1.62, -1.55, -0.22, -0.02, 0.32, 0.29, 1.59, -0.9, -1.12, -1, -1.38, -1.03, -1.78, 0.55, -0.59, 0.65, -0.38]}' data-zero="1" data-caption="Notice the individual scatter (SD about 1) is much bigger than the shift of the mean (about -0.47). The mean can be measured precisely only because it averages 60 children."></div>

The standard error is 1/√60 = 0.129, so the mean of −0.471 lies 3.65 standard errors below 0. If the programme made no difference, z values would follow the standard normal curve below. Ours falls far out in the left tail, so the shaded two-sided p-value area is tiny.

<div class="viz" data-viz="z" data-stat="-3.65"></div>

Notice that more children shrink the noise (σ/√n) and so make the same mean shift look more convincing. That is why a very large study can find a tiny, clinically meaningless shift "significant".

## When to use it

Use the one-sample z-test when **all** of these are true:

- The outcome is **continuous**
- There is **one group** of independent subjects and a reference value μ₀ fixed before looking at the data
- The population SD σ is **known** from a source other than your sample (not "estimated from a large previous study" with a lot of uncertainty)
- The data are approximately normal, or n is large enough that the sample mean is approximately normal

**Where σ is genuinely known in clinical work:**

- **Standardized scores.** Growth-chart z-scores (WHO height-for-age, weight-for-age, BMI-for-age) are defined so that the reference population has mean 0 and SD 1. IQ-like scales are built to have mean 100 and SD 15.
- **Well-calibrated assays.** A laboratory assay with a long, documented calibration history may have a stable measurement SD (for example from thousands of control runs), so the SD of replicate measurements of a reference material is effectively known.

**Clinical scenario.** A nutrition programme enrols 60 children. Their height-for-age z-score (HAZ) is computed against the WHO reference population, for which μ₀ = 0 and σ = 1 by construction. *Is the mean HAZ of children in the programme different from the reference population?* The reference SD of 1 is a property of the scale, not an estimate from these 60 children.

**Don't use it when:**

- σ is **estimated from your own sample** (the usual case) → use the [one-sample t-test](#/test/one-sample-t)
- You are comparing **two groups** → use the [independent t-test](#/test/independent-t) or the [paired t-test](#/test/paired-t)
- The outcome is a **binary** event rate → use the [one-sample proportion test](#/test/one-proportion), or the [two-proportion z-test](#/test/two-proportion-z) for two groups
- The data are strongly skewed with a small n → use the [Wilcoxon signed-rank test](#/test/wilcoxon-signed-rank)

**z-tests for proportions.** The z-test you meet most often in practice is for **proportions**, where it is legitimate because under H₀ the variance is fixed by the proportion itself: σ² = p₀(1 − p₀). See the [one-sample proportion test](#/test/one-proportion) and the [two-proportion z-test](#/test/two-proportion-z).

## Assumptions

1. **σ is known** and applies to the population you are sampling from. If the programme children could plausibly have a different spread than the reference population, the "known σ" is not valid for them, and the z-test will be too liberal or too conservative. Compare your sample SD with σ as a sanity check (it should be close, but do not substitute it).
2. **Observations are independent.** One value per child; siblings or clustered clinics need a different analysis.
3. **The mean is approximately normal.** Check with a Q–Q plot or [Shapiro–Wilk test](#/test/shapiro-wilk). With n ≥ 30 the sample mean is usually close to normal even if the data are not.
4. **μ₀ and σ were fixed before analysis.**

## Hypotheses

- **H₀:** μ = μ₀ (the programme's mean HAZ equals the reference mean, here 0)
- **H₁:** μ ≠ μ₀ (two-sided)

## Example data

Height-for-age z-scores (HAZ) in 60 children enrolled in the programme (simulated; reference μ₀ = 0, σ = 1):

| Children | HAZ values |
|---|---|
| 1–10 | 0.63, 0.12, −0.46, −0.56, 0.81, 0.94, 0.18, −0.48, −1.57, −1.47 |
| 11–20 | −2.02, 0.12, 0.48, −0.06, −3.62, −0.50, −1.28, 0.10, −1.96, −2.36 |
| 21–30 | 0.79, −0.03, 0.69, −0.56, 0.05, −1.89, −0.47, 0.33, −2.41, −0.46 |
| 31–40 | −1.79, −0.68, −0.86, 0.85, −0.85, 0.04, −0.01, −0.84, −0.92, 0.08 |
| 41–50 | −1.22, −0.35, 0.15, 1.62, −1.55, −0.22, −0.02, 0.32, 0.29, 1.59 |
| 51–60 | −0.90, −1.12, −1.00, −1.38, −1.03, −1.78, 0.55, −0.59, 0.65, −0.38 |

## How to run it

### R

```r
# HAZ in 60 children; reference population mean 0, SD 1
haz <- c(0.63, 0.12, -0.46, -0.56, 0.81, 0.94, 0.18, -0.48, -1.57, -1.47,
         -2.02, 0.12, 0.48, -0.06, -3.62, -0.5, -1.28, 0.1, -1.96, -2.36,
         0.79, -0.03, 0.69, -0.56, 0.05, -1.89, -0.47, 0.33, -2.41, -0.46,
         -1.79, -0.68, -0.86, 0.85, -0.85, 0.04, -0.01, -0.84, -0.92, 0.08,
         -1.22, -0.35, 0.15, 1.62, -1.55, -0.22, -0.02, 0.32, 0.29, 1.59,
         -0.9, -1.12, -1, -1.38, -1.03, -1.78, 0.55, -0.59, 0.65, -0.38)
mu0 <- 0; sigma <- 1            # KNOWN population SD
n <- length(haz); xbar <- mean(haz)

# z-test by hand
se <- sigma / sqrt(n)
z  <- (xbar - mu0) / se
p  <- 2 * pnorm(-abs(z))
ci <- xbar + c(-1, 1) * qnorm(0.975) * se
cat(sprintf("mean = %.3f, SE = %.4f, z = %.3f, p = %.5f\n", xbar, se, z, p))
cat(sprintf("95%% CI (known sigma): %.3f to %.3f\n", ci[1], ci[2]))

# For comparison: one-sample t-test (SD estimated from data)
cat(sprintf("sample SD = %.3f\n", sd(haz)))
t.test(haz, mu = mu0)
```

**What the code does**

- Base R has **no** `z.test()`, so the test is computed by hand. This is only three lines and makes the role of σ explicit. (The `BSDA` package offers `BSDA::z.test(x, mu = 0, sigma.x = 1)` if you prefer a ready-made function.)
- `sigma <- 1` is the **known** population SD, supplied by you, never `sd(haz)`.
- `pnorm(-abs(z))` is the one-sided tail area; doubling it gives the two-sided p-value.
- `qnorm(0.975)` is 1.96, the multiplier for a 95% CI. The CI is based on the known σ.
- `t.test(haz, mu = 0)` is run on the same data for comparison. It estimates the SD from the sample (1.020), uses a t distribution with n − 1 = 59 df, and therefore gives a slightly wider CI and larger p-value.

**Output**

```text
mean = -0.471, SE = 0.1291, z = -3.650, p = 0.00026
95% CI (known sigma): -0.724 to -0.218
sample SD = 1.020

	One Sample t-test

data:  haz
t = -3.5783, df = 59, p-value = 0.0006991
alternative hypothesis: true mean is not equal to 0
95 percent confidence interval:
 -0.7346465 -0.2076869
sample estimates:
 mean of x 
-0.4711667 
```

### Python

```python
import numpy as np
from scipy import stats
from statsmodels.stats.weightstats import ztest

# HAZ in 60 children; reference population mean 0, SD 1
haz = np.array([0.63, 0.12, -0.46, -0.56, 0.81, 0.94, 0.18, -0.48, -1.57, -1.47,
                -2.02, 0.12, 0.48, -0.06, -3.62, -0.5, -1.28, 0.1, -1.96, -2.36,
                0.79, -0.03, 0.69, -0.56, 0.05, -1.89, -0.47, 0.33, -2.41, -0.46,
                -1.79, -0.68, -0.86, 0.85, -0.85, 0.04, -0.01, -0.84, -0.92, 0.08,
                -1.22, -0.35, 0.15, 1.62, -1.55, -0.22, -0.02, 0.32, 0.29, 1.59,
                -0.9, -1.12, -1, -1.38, -1.03, -1.78, 0.55, -0.59, 0.65, -0.38])
mu0, sigma = 0.0, 1.0            # sigma is KNOWN
n, xbar = len(haz), haz.mean()

# z-test by hand with the known sigma
se = sigma / np.sqrt(n)
z = (xbar - mu0) / se
p = 2 * stats.norm.sf(abs(z))
lo, hi = xbar + np.array([-1, 1]) * stats.norm.ppf(0.975) * se
print(f"mean = {xbar:.3f}, SE = {se:.4f}, z = {z:.3f}, p = {p:.5f}")
print(f"95% CI (known sigma): {lo:.3f} to {hi:.3f}")

# statsmodels ztest: uses the SAMPLE SD (ddof=1), not your sigma
zs, ps = ztest(haz, value=mu0)
print(f"statsmodels ztest: z = {zs:.3f}, p = {ps:.5f}  (sample SD = {haz.std(ddof=1):.3f})")

# One-sample t-test for comparison
r = stats.ttest_1samp(haz, mu0)
ci = r.confidence_interval()
print(f"t-test: t = {r.statistic:.3f}, df = {r.df}, p = {r.pvalue:.5f}, "
      f"95% CI {ci.low:.3f} to {ci.high:.3f}")
```

**What the code does**

- The manual calculation with `scipy.stats.norm` mirrors the R code: `norm.sf(abs(z))` is the upper-tail probability (more accurate than `1 - cdf` for small p), and `norm.ppf(0.975)` is 1.96.
- `sigma = 1.0` is the **known** σ and is the only place it enters.
- ⚠️ **Pitfall:** `statsmodels.stats.weightstats.ztest(x, value=mu0)` does **not** accept a known σ. It computes the standard error from the **sample SD** (`ddof=1`) and then refers the statistic to the normal distribution. It is therefore a "t-statistic with a normal p-value": z = −3.578 here, identical to the t statistic, instead of the correct −3.650 from σ = 1. The p-value (0.00035) is too small compared with the t-test's 0.00070 because the normal distribution ignores the uncertainty in the estimated SD. Use `ztest` only if you accept this approximation (large n), and never describe its result as a known-σ z-test.
- `stats.ttest_1samp` gives the proper t-test on the same data. `res.confidence_interval()` needs SciPy ≥ 1.10.

**Output**

```text
mean = -0.471, SE = 0.1291, z = -3.650, p = 0.00026
95% CI (known sigma): -0.724 to -0.218
statsmodels ztest: z = -3.578, p = 0.00035  (sample SD = 1.020)
t-test: t = -3.578, df = 59, p = 0.00070, 95% CI -0.735 to -0.208
```

## Interpreting the output

- **Mean HAZ = −0.471.** On average, children in the programme are about half a reference standard deviation shorter for their age than the reference population.
- **Standard error = 1/√60 = 0.129.** It uses the known σ = 1.
- **z = −3.65, p = 0.00026.** If the programme's true mean HAZ were 0, a sample mean this far from 0 would occur in about 3 of 10,000 samples.
- **95% CI (known σ): −0.72 to −0.22.** The data are consistent with a true mean HAZ between about 0.2 and 0.7 SD below the reference. The CI excludes 0, in agreement with the p-value.
- **Comparison with the t-test.** The sample SD is 1.020, close to the assumed 1, so the results are similar: t(59) = −3.58, p = 0.0007, 95% CI −0.73 to −0.21. The t result is a little more conservative because it uses the sample SD (1.020, slightly above 1) and a t distribution that allows for estimating it. With n = 60 the conclusion does not change, but the t-test remains the safe default whenever σ is not truly known.
- **statsmodels `ztest` gives z = −3.578**, the same as t, because it silently used the sample SD.

**How to report it:**
> The mean height-for-age z-score of the 60 children was −0.47 (95% CI −0.72 to −0.22), significantly below the WHO reference mean of 0 (one-sample z-test with known σ = 1, z = −3.65, p < 0.001).

⚠️ A mean HAZ below 0 describes the group, not each child. The proportion of children with HAZ < −2 (stunting) is often the more clinically meaningful quantity; analyse it with a [proportion test](#/test/one-proportion).

## Check your understanding

<details class="check"><summary>1. A clinic measures the HbA1c of 50 patients and wants to compare the mean with 7.0%. It uses the sample SD of 0.9 and calls the result a "z-test because n ≥ 30". Is that right?</summary>

No. The SD was estimated from the sample, so the one-sample [t-test](#/test/one-sample-t) is the correct choice at any n. With n = 50 the answers will be nearly identical, but "n ≥ 30 so use z" is a printed-table shortcut, not a rule. A z-test needs a σ that comes from outside the data.
</details>

<details class="check"><summary>2. The z-test gives p = 0.0003 for a mean HAZ of −0.47. A colleague says "the programme is harming children's growth." What is missing?</summary>

Three things. (1) A mean difference from a reference says nothing about *cause*: the children may have started below 0 before the programme. (2) The effect size matters: −0.47 SD with 95% CI −0.72 to −0.22 is moderate, and clinical importance should be judged on that, not on p. (3) The mean describes the group, not every child; the proportion with HAZ < −2 may be more informative.
</details>

<details class="check"><summary>3. What would happen to z if the same mean (−0.471) had been observed in only 15 children instead of 60?</summary>

The signal is unchanged but the standard error doubles (1/√15 = 0.258 instead of 0.129), so z halves to about −1.82 and the two-sided p-value rises to roughly 0.07. Same effect, less information, weaker evidence. The wider CI (−0.98 to 0.03) would show that the data are also compatible with no shift.
</details>

<details class="check"><summary>4. Python's statsmodels `ztest` gives z = −3.578 for these data but the hand calculation gives −3.650. Which is the real z-test with known σ = 1, and why do they differ?</summary>

The hand calculation (−3.650) is the z-test with σ = 1. `statsmodels.stats.weightstats.ztest` ignores your σ and uses the sample SD (1.020), so its statistic equals the t statistic (−3.578), only referred to the normal curve. Always check what a software "z-test" uses for the SD.
</details>

## Common mistakes

- **Using the z-test because n ≥ 30.** The rule is obsolete. If σ is estimated, use the t-test; it gives the right answer at any n.
- **Plugging in the sample SD and calling it a z-test.** This is what `statsmodels.ztest` does. Unless σ is known, it is a t-test with the wrong reference distribution.
- **Treating a published SD as "known".** An SD estimated from a previous study is itself an estimate, with its own error. Only a σ that is defined by the scale (z-scores, IQ) or established by extensive calibration counts.
- **Assuming the reference σ applies to your population.** A malnourished population may well be more variable than the reference. Check the sample SD, and if it differs a lot, use the t-test.
- **Reporting only the p-value.** Give the mean difference and its 95% CI.
- **Testing many reference values after looking at the data.** Fix μ₀ in advance.

## Related tests

- [One-sample t-test](#/test/one-sample-t): the default when σ is estimated
- [One-sample proportion test](#/test/one-proportion): z-test for a single binary rate
- [Two-proportion z-test](#/test/two-proportion-z): compares two binary rates
- [Wilcoxon signed-rank test](#/test/wilcoxon-signed-rank): nonparametric alternative
- [Shapiro–Wilk test](#/test/shapiro-wilk): normality check

## Key takeaways

- A z-test measures (x̄ − μ₀) in units of the standard error σ/√n, with σ **known from outside the sample**.
- If σ is estimated from your data, use the t-test, whatever n is. "n ≥ 30 so use z" is a myth.
- The statsmodels `ztest` uses the sample SD, so it is a t statistic with a normal p-value.
- Under H₀, z follows N(0, 1). |z| > 1.96 gives p < 0.05 (two-sided).
- Report the mean difference with its 95% CI. A small p does not say how large or how important the effect is.

## References

- Altman DG. *Practical Statistics for Medical Research.* Chapman & Hall; 1991. Ch. 9.
- WHO Multicentre Growth Reference Study Group. *WHO Child Growth Standards: Length/height-for-age, weight-for-age, weight-for-length, weight-for-height and body mass index-for-age.* Geneva: WHO; 2006.
- statsmodels documentation: `statsmodels.stats.weightstats.ztest`
- SciPy documentation: `scipy.stats.norm`, `scipy.stats.ttest_1samp`; R documentation: `?pnorm`, `?qnorm`
