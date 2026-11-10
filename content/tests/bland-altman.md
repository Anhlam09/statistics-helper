# Bland–Altman analysis

> Assesses agreement between two methods that measure the same quantity on the same patients by estimating the average difference (bias) and the range within which 95% of differences are expected to fall (the limits of agreement).

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Judge whether two methods can be used interchangeably | Continuous, both methods on the same patients | Independent pairs; differences roughly normal with constant spread | Point-of-care vs laboratory HbA1c in 30 patients |

## Learning objectives

- Explain why a high correlation does not show that two methods agree
- Compute the bias and 95% limits of agreement from the paired differences
- Read a Bland–Altman plot for bias, trends, and fanning
- Judge the limits against a clinically acceptable difference set in advance

## What it is

Suppose a new point-of-care device and the reference laboratory both measure HbA1c on the same patients. A high correlation does not tell you whether they *agree*; it only says that one rises as the other rises. Bland and Altman's approach looks directly at the **differences**:

1. For each patient, compute the difference d = (new method − reference) and the mean of the two readings.
2. Estimate the **bias**: the mean difference d̄ (a systematic over- or under-reading).
3. Estimate the **95% limits of agreement (LoA)**: d̄ ± 1.96 × s_d, where s_d is the SD of the differences. About 95% of the differences between the two methods are expected to fall between these limits.
4. Plot d against the mean of the two methods (the Bland–Altman plot) with the bias line and the two limit lines, and look for patterns.

The statistics say how big the disagreement is in the **units of measurement**. Whether the limits are acceptable is a **clinical decision** made *before* the study (for example, differences within ±0.5% HbA1c would not change management).

- **Test statistic:** none is needed for agreement itself. The bias has a one-sample (paired) t-test, t = d̄ / (s_d / √n), but that tests whether the bias is zero, not whether agreement is adequate.
- **Effect size:** bias with 95% CI; 95% limits of agreement, each with its own 95% CI (the limits are estimated, especially with small n, and their CIs are wide)

**Why correlation is NOT agreement.** Pearson's r measures the strength of a linear relationship, not closeness to the line of equality. A device that always reads 2 units too high, or reads exactly twice the true value, correlates perfectly (r = 1) with the reference but agrees with it poorly. Correlation also depends on the spread of the patients: a wide range of true values gives a high r almost regardless of measurement error. In the example below r = 0.975, which looks excellent, yet the device reads on average 0.28 percentage points high and individual differences range from −0.6 to +1.0.

## The intuition

Agreement is about the **differences**, so Bland–Altman strips the problem down to one number per patient: the POC reading minus the laboratory reading. The **signal** is the typical offset (bias, +0.28) and the **noise** is how much individual patients scatter around it (SD 0.41). Two methods can be used interchangeably only if both the offset and the scatter are small compared with what matters clinically.

Analogy: two bathroom scales. If one always reads 1 kg heavier, they correlate perfectly, but you cannot swap them for tracking a 0.5 kg change. Correlation tells you that heavier people read higher on both scales; it says nothing about the numbers matching.

The dot plot shows the 30 differences. Notice the line at 0 (perfect agreement): most dots lie to its right (POC reads higher), and the mean difference is +0.28. The spread from −0.6 to +1.0 is the real clinical story: the limits of agreement −0.52 to +1.09 are where about 95% of individual differences are expected to fall, assuming the differences are roughly normal.

<div class="viz" data-viz="dots" data-unit="HbA1c %, POC minus laboratory" data-groups='{"Difference":[0.2,0.0,-0.6,0.3,0.1,0.5,0.5,0.7,0.4,0.5,0.1,0.7,0.7,0.0,-0.1,0.3,0.6,1.0,0.6,0.7,0.5,-0.2,-0.5,-0.4,0.2,0.7,-0.1,0.2,0.9,0.0]}' data-zero="1" data-caption="The mean difference (bias) is about +0.28. Nine of 30 patients differ by more than +0.5."></div>

The paired t-test on these differences (t = 3.78, df = 29) tests H₀: the true mean difference (bias) is 0, and the curve below is the t distribution under that H₀. It says the bias is not zero. That answers a different question from "is agreement good enough?".

<div class="viz" data-viz="t" data-stat="3.78" data-df="29"></div>

**What to notice:** the t statistic is far in the tail, so the bias is statistically detectable. But the bias is not the whole story: even if you subtracted 0.28 from every POC value, individual differences would still scatter by about ±0.8. Statistical significance of the bias does not tell you whether the scatter is acceptable.

## When to use it

Use Bland–Altman analysis when **all** of these are true:

- The measurement is **continuous**, and the **same patients** are measured by **both methods** (paired)
- The question is whether the two methods can be **used interchangeably** (or whether a new, cheaper, or less invasive method can replace the old)
- The differences are approximately normally distributed, with spread that does not change with the size of the measurement (otherwise transform or use a regression-based version)

**Clinical scenario.** A clinic wants to replace central-laboratory HbA1c with a point-of-care (POC) device that gives results during the visit. HbA1c (%) is measured on 30 patients by both methods. *By how much does the POC device differ from the laboratory, and within what range will most individual results differ?*

**Don't use it when:**

- You want the reliability of repeated or multi-rater measurements, with variance components → use the [ICC](#/test/icc)
- The ratings are categorical → use [Cohen's kappa](#/test/cohens-kappa)
- You only want to test whether the mean measurement differs between conditions → use the [paired t-test](#/test/paired-t)
- You want to relate the measurement to another variable → use [Pearson](#/test/pearson) or [linear regression](#/test/linear-regression)
- One method is a gold standard being used to classify disease → use [ROC analysis](#/test/roc-auc)
- Each patient has several readings by each method → use extended Bland–Altman methods for repeated measures (not covered here)

## Assumptions

1. **The pairs are independent** (one pair per patient, or properly handled repeated measures).
2. **The differences are approximately normal.** Check with a histogram, Q–Q plot, or [Shapiro–Wilk test](#/test/shapiro-wilk). Otherwise use percentile-based limits or log-transform the readings.
3. **The bias and spread are constant across the measurement range.** On the plot, the points should form a horizontal band without a trend or fan shape. If the spread grows with the mean, plot ratios or log differences. If the differences trend with the mean (proportional bias), fit a regression of d on the mean.
4. **The range of the data covers the clinically relevant range.** Limits of agreement apply only within the range studied.
5. **Adequate sample size.** The CIs of the LoA are wide with fewer than about 50 pairs, so plan the sample size around the precision of the limits.

## Hypotheses

Bland–Altman analysis is descriptive and estimation-focused; there is no single H₀ for agreement. What is estimated:

- **Bias** μ_d = the mean of (method A − method B), with 95% CI. The associated test is H₀: μ_d = 0 against H₁: μ_d ≠ 0.
- **Limits of agreement** μ_d ± 1.96σ_d, with a CI for each limit.
- **Decision rule** (pre-specified): the methods are interchangeable if both limits, including their outer CI limits, lie within the clinically acceptable difference.

## Example data

HbA1c (%) in 30 patients by central laboratory (reference) and by the point-of-care (POC) device.

| Patient | Lab | POC | Patient | Lab | POC | Patient | Lab | POC |
|---|---|---|---|---|---|---|---|---|
| 1 | 6.2 | 6.4 | 11 | 6.6 | 6.7 | 21 | 10.3 | 10.8 |
| 2 | 9.1 | 9.1 | 12 | 7.9 | 8.6 | 22 | 9.3 | 9.1 |
| 3 | 10.5 | 9.9 | 13 | 6.9 | 7.6 | 23 | 6.3 | 5.8 |
| 4 | 6.7 | 7 | 14 | 8.4 | 8.4 | 24 | 6.4 | 6 |
| 5 | 5.6 | 5.7 | 15 | 6.6 | 6.5 | 25 | 5.8 | 6 |
| 6 | 9.2 | 9.7 | 16 | 6.2 | 6.5 | 26 | 7.9 | 8.6 |
| 7 | 8.2 | 8.7 | 17 | 7.3 | 7.9 | 27 | 7.6 | 7.5 |
| 8 | 9.8 | 10.5 | 18 | 10.3 | 11.3 | 28 | 10.8 | 11 |
| 9 | 10.7 | 11.1 | 19 | 8.3 | 8.9 | 29 | 5.9 | 6.8 |
| 10 | 5.7 | 6.2 | 20 | 10.1 | 10.8 | 30 | 10.7 | 10.7 |

## How to run it

### R

```r
# HbA1c (%) in 30 patients: central laboratory (reference) vs point-of-care (POC) device
dat <- read.table(text = "
patient lab poc
1 6.2 6.4
2 9.1 9.1
3 10.5 9.9
4 6.7 7
5 5.6 5.7
6 9.2 9.7
7 8.2 8.7
8 9.8 10.5
9 10.7 11.1
10 5.7 6.2
11 6.6 6.7
12 7.9 8.6
13 6.9 7.6
14 8.4 8.4
15 6.6 6.5
16 6.2 6.5
17 7.3 7.9
18 10.3 11.3
19 8.3 8.9
20 10.1 10.8
21 10.3 10.8
22 9.3 9.1
23 6.3 5.8
24 6.4 6
25 5.8 6
26 7.9 8.6
27 7.6 7.5
28 10.8 11
29 5.9 6.8
30 10.7 10.7
", header = TRUE)

# 1. Differences and means
d <- dat$poc - dat$lab                   # POC minus laboratory
m <- (dat$poc + dat$lab) / 2
n <- length(d); bias <- mean(d); s <- sd(d)

# 2. Bias and 95% limits of agreement
loa <- bias + c(-1, 1) * qnorm(0.975) * s

# 3. 95% CIs for the bias and for each limit (Bland & Altman 1986, t-based)
tq      <- qt(0.975, n - 1)
ci_bias <- bias + c(-1, 1) * tq * s / sqrt(n)
se_loa  <- s * sqrt(1 / n + qnorm(0.975)^2 / (2 * (n - 1)))
ci_lo   <- loa[1] + c(-1, 1) * tq * se_loa
ci_hi   <- loa[2] + c(-1, 1) * tq * se_loa
cat(sprintf("n = %d, bias = %.3f (95%% CI %.3f to %.3f), SD of differences = %.3f\n",
            n, bias, ci_bias[1], ci_bias[2], s))
cat(sprintf("Limits of agreement: %.3f (CI %.3f to %.3f) to %.3f (CI %.3f to %.3f)\n",
            loa[1], ci_lo[1], ci_lo[2], loa[2], ci_hi[1], ci_hi[2]))

# 4. Why correlation is not agreement, and a check on normality of the differences
cat(sprintf("Pearson r = %.3f\n", cor(dat$lab, dat$poc)))
print(shapiro.test(d))

# 5. Bland-Altman plot
plot(m, d, pch = 19, xlab = "Mean of POC and laboratory HbA1c (%)",
     ylab = "Difference (POC - laboratory, % HbA1c)", ylim = c(-1.2, 1.2))
abline(h = bias, lwd = 2)
abline(h = loa, lty = 2, col = "red")
abline(h = 0, col = "grey60")
```

**What the code does**

- `d <- poc - lab` fixes the direction: a positive bias means the new method reads higher. Always say which way the difference was taken.
- `loa` is d̄ ± 1.96 s_d, using the normal quantile (`qnorm(0.975)`). Some authors use 2 instead of 1.96; the difference is small.
- The CIs follow Bland & Altman (1986): the bias has SE = s_d/√n, and each limit has SE = s_d × √(1/n + 1.96²/(2(n−1))). Both use a t quantile with n − 1 df.
- `cor()` and `shapiro.test(d)` demonstrate that correlation is high even though the agreement is imperfect, and check the normality assumption on the **differences**.
- The plot uses base R only: the points, a solid line at the bias, dashed red lines at the limits, and a grey reference line at 0 (perfect agreement). The plot is drawn on screen (or into `Rplots.pdf` when run with `Rscript`).

**Output**

```text
n = 30, bias = 0.283 (95% CI 0.130 to 0.437), SD of differences = 0.410
Limits of agreement: -0.521 (CI -0.786 to -0.256) to 1.087 (CI 0.823 to 1.352)
Pearson r = 0.975

	Shapiro-Wilk normality test

data:  d
W = 0.96614, p-value = 0.4397
```

### Python

```python
import io
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

# HbA1c (%) in 30 patients: central laboratory (reference) vs point-of-care (POC) device
dat = pd.read_csv(io.StringIO("""
patient lab poc
1 6.2 6.4
2 9.1 9.1
3 10.5 9.9
4 6.7 7
5 5.6 5.7
6 9.2 9.7
7 8.2 8.7
8 9.8 10.5
9 10.7 11.1
10 5.7 6.2
11 6.6 6.7
12 7.9 8.6
13 6.9 7.6
14 8.4 8.4
15 6.6 6.5
16 6.2 6.5
17 7.3 7.9
18 10.3 11.3
19 8.3 8.9
20 10.1 10.8
21 10.3 10.8
22 9.3 9.1
23 6.3 5.8
24 6.4 6
25 5.8 6
26 7.9 8.6
27 7.6 7.5
28 10.8 11
29 5.9 6.8
30 10.7 10.7
"""), sep=r"\s+")

# 1. Differences and means
d = (dat["poc"] - dat["lab"]).to_numpy()          # POC minus laboratory
m = ((dat["poc"] + dat["lab"]) / 2).to_numpy()
n = len(d); bias = d.mean(); s = d.std(ddof=1)    # ddof=1: sample SD

# 2. Bias and 95% limits of agreement
z = stats.norm.ppf(0.975)
loa = bias + np.array([-1, 1]) * z * s

# 3. 95% CIs for the bias and for each limit (Bland & Altman 1986, t-based)
tq = stats.t.ppf(0.975, n - 1)
ci_bias = bias + np.array([-1, 1]) * tq * s / np.sqrt(n)
se_loa = s * np.sqrt(1 / n + z**2 / (2 * (n - 1)))
ci_lo = loa[0] + np.array([-1, 1]) * tq * se_loa
ci_hi = loa[1] + np.array([-1, 1]) * tq * se_loa
print(f"n = {n}, bias = {bias:.3f} (95% CI {ci_bias[0]:.3f} to {ci_bias[1]:.3f}), SD of differences = {s:.3f}")
print(f"Limits of agreement: {loa[0]:.3f} (CI {ci_lo[0]:.3f} to {ci_lo[1]:.3f}) "
      f"to {loa[1]:.3f} (CI {ci_hi[0]:.3f} to {ci_hi[1]:.3f})")

# 4. Why correlation is not agreement, and a check on normality of the differences
print(f"Pearson r = {np.corrcoef(dat['lab'], dat['poc'])[0, 1]:.3f}")
sw = stats.shapiro(d)
print(f"Shapiro-Wilk W = {sw.statistic:.4f}, p = {sw.pvalue:.4f}")

# 5. Bland-Altman plot
fig, ax = plt.subplots(figsize=(6, 4))
ax.scatter(m, d, color="black")
ax.axhline(bias, color="black", lw=2, label=f"Bias {bias:.2f}")
ax.axhline(loa[0], color="red", ls="--", label="95% limits of agreement")
ax.axhline(loa[1], color="red", ls="--")
ax.axhline(0, color="grey", lw=1)
ax.set_xlabel("Mean of POC and laboratory HbA1c (%)")
ax.set_ylabel("Difference (POC - laboratory, % HbA1c)")
ax.legend(); fig.tight_layout(); fig.savefig("bland_altman.png", dpi=150)
```

**What the code does**

- `d.std(ddof=1)` is the sample SD (n − 1 in the denominator). NumPy's default (`ddof=0`) would give a slightly narrower interval.
- `stats.norm.ppf(0.975)` is 1.96 and `stats.t.ppf(0.975, n - 1)` is the t quantile used for the CIs.
- The CI formulas are identical to the R code, so the numbers match exactly.
- `matplotlib` draws the Bland–Altman plot and `fig.savefig("bland_altman.png")` writes it to a file; use `plt.show()` to display it interactively.

**Output**

```text
n = 30, bias = 0.283 (95% CI 0.130 to 0.437), SD of differences = 0.410
Limits of agreement: -0.521 (CI -0.786 to -0.256) to 1.087 (CI 0.823 to 1.352)
Pearson r = 0.975
Shapiro-Wilk W = 0.9661, p = 0.4397
```

### SAS

```sas
data hba1c;
  input patient lab poc;
  datalines;
1 6.2 6.4
2 9.1 9.1
3 10.5 9.9
4 6.7 7
5 5.6 5.7
6 9.2 9.7
7 8.2 8.7
8 9.8 10.5
9 10.7 11.1
10 5.7 6.2
11 6.6 6.7
12 7.9 8.6
13 6.9 7.6
14 8.4 8.4
15 6.6 6.5
16 6.2 6.5
17 7.3 7.9
18 10.3 11.3
19 8.3 8.9
20 10.1 10.8
21 10.3 10.8
22 9.3 9.1
23 6.3 5.8
24 6.4 6
25 5.8 6
26 7.9 8.6
27 7.6 7.5
28 10.8 11
29 5.9 6.8
30 10.7 10.7
;
run;

/* 1. Differences and means (POC minus laboratory) */
data ba;
  set hba1c;
  diff = poc - lab;
  avg  = (poc + lab) / 2;
run;

/* 2. Bias, SD, 95% limits of agreement and their CIs */
proc means data=ba noprint;
  var diff;
  output out=st n=n mean=bias std=sd;
run;

data st;
  set st;
  z  = quantile('normal', 0.975);
  tq = quantile('t', 0.975, n - 1);
  loa_lo = bias - z * sd;
  loa_hi = bias + z * sd;
  se_bias = sd / sqrt(n);
  se_loa  = sd * sqrt(1 / n + z**2 / (2 * (n - 1)));
  bias_lcl = bias - tq * se_bias;   bias_ucl = bias + tq * se_bias;
  lo_lcl   = loa_lo - tq * se_loa;  lo_ucl   = loa_lo + tq * se_loa;
  hi_lcl   = loa_hi - tq * se_loa;  hi_ucl   = loa_hi + tq * se_loa;
  call symputx('bias',   bias);
  call symputx('loa_lo', loa_lo);
  call symputx('loa_hi', loa_hi);
run;

proc print data=st noobs;
  var n bias bias_lcl bias_ucl sd loa_lo lo_lcl lo_ucl loa_hi hi_lcl hi_ucl;
  format bias bias_lcl bias_ucl sd loa_lo lo_lcl lo_ucl loa_hi hi_lcl hi_ucl 8.3;
run;

/* 3. Correlation (for contrast) and normality of the differences */
proc corr data=ba; var lab poc; run;
proc univariate data=ba normal; var diff; run;

/* 4. Bland-Altman plot */
proc sgplot data=ba;
  scatter x=avg y=diff;
  refline &bias   / axis=y lineattrs=(thickness=2);
  refline &loa_lo &loa_hi / axis=y lineattrs=(pattern=dash color=red);
  refline 0 / axis=y lineattrs=(color=gray);
  xaxis label="Mean of POC and laboratory HbA1c (%)";
  yaxis label="Difference (POC - laboratory, % HbA1c)";
run;
```

**What the code does**

- The first `DATA` step creates the difference and the pairwise mean. Because the new method is first in `poc - lab`, a positive bias means the POC device reads higher.
- `PROC MEANS ... output out=st` stores n, the mean difference (the bias) and the SD of the differences.
- The next `DATA` step computes the limits of agreement, `quantile('normal', 0.975)` = 1.96, and the same Bland–Altman CI formulas as in R and Python, using `quantile('t', 0.975, n-1)` for the t value. `call symputx` saves the bias and limits as macro variables for the plot.
- `PROC CORR` gives Pearson's r for contrast with the agreement statistics. `PROC UNIVARIATE ... normal` provides the Shapiro–Wilk test of the differences.
- `PROC SGPLOT` draws the scatter of the difference against the mean, with `REFLINE` statements for the bias (solid), the two limits (dashed red), and zero.
- All numbers should agree with R and Python, because the same formulas are used.

## Interpreting the output

- **Bias = +0.283 percentage points of HbA1c** (95% CI 0.130 to 0.437). The POC device reads about 0.3 percentage points higher than the laboratory on average. The CI excludes 0, so the bias is statistically detectable.
- **SD of differences = 0.410.**
- **95% limits of agreement: −0.521 to +1.087.** For a new patient, the POC reading is expected to lie between 0.52 below and 1.09 above the laboratory reading in 95% of cases.
- **CIs of the limits:** lower limit −0.786 to −0.256; upper limit 0.823 to 1.352. These show the uncertainty from having only 30 patients.
- **Clinical judgement.** If the pre-specified acceptable difference were ±0.5 percentage points, the upper limit (1.09, even its lower CI limit of 0.82) exceeds it, so the devices would not be interchangeable without a bias correction. A constant bias like this one could be corrected, which would shift the limits to about ±0.8 around zero; whether that is acceptable depends on the clinical use.
- **Pearson r = 0.975**, an excellent correlation that hides this disagreement. This is the key reason not to use correlation as a measure of agreement.
- **Shapiro–Wilk on the differences: p = 0.44**, so there is no evidence against normality of the differences. Check the plot as well: the points should form a horizontal band with no trend or fanning.
- R and Python agree exactly to three decimals.

**How to report it:**
> Agreement between point-of-care and laboratory HbA1c was assessed in 30 patients using Bland–Altman analysis. The point-of-care device read higher by a mean of 0.28 percentage points (95% CI 0.13 to 0.44). The 95% limits of agreement were −0.52 (95% CI −0.79 to −0.26) to +1.09 (95% CI 0.82 to 1.35) percentage points, exceeding the pre-specified acceptable difference of ±0.5.

## Check your understanding

<details class="check"><summary>1. A new device and the reference method have Pearson r = 0.99. Does this show they agree?</summary>

No. r measures how well one is a linear function of the other, not whether the values match. A device reading a constant 2 units too high, or twice the true value, has r = 1. Look at the differences: the bias and the limits of agreement.
</details>

<details class="check"><summary>2. The paired t-test on the differences gives p = 0.40 in a small study. A colleague concludes the methods agree. What is wrong?</summary>

p > 0.05 does not show agreement. A small sample can fail to detect a bias, and the test says nothing about the scatter of individual differences. Report the bias with its CI and the limits of agreement, with their CIs, and compare them with a pre-specified acceptable difference.
</details>

<details class="check"><summary>3. If the acceptable difference is ±0.5 HbA1c and the limits are −0.52 to +1.09, can the devices be swapped?</summary>

Not as they are. The upper limit (and even the lower end of its CI, 0.82) exceeds +0.5. Correcting the constant bias of 0.28 would shift the limits to about ±0.8, which still exceeds ±0.5. So the POC device fails the criterion in this study.
</details>

<details class="check"><summary>4. What if the spread of differences were wider at high HbA1c values (a fan shape on the plot)?</summary>

A single bias and fixed limits would be misleading: too wide at low values and too narrow at high ones. Use log-transformed data (ratios), or a regression-based Bland–Altman analysis where the limits depend on the mean.
</details>

## Common mistakes

- **Using correlation (or a paired t-test alone) to claim agreement.** High r says nothing about bias or the size of individual differences, and a non-significant paired t-test with a small n can reflect low power.
- **Not defining acceptable limits in advance.** Bland–Altman gives numbers, not a verdict. State the clinically acceptable difference before the study.
- **Omitting the CIs of the limits of agreement.** With small n the limits are imprecise; report them.
- **Confusing LoA with a CI for the mean.** The LoA describe individual differences (about 95% of patients), not the uncertainty of the mean.
- **Ignoring proportional bias or heteroscedasticity.** If differences grow with the mean, a single bias and fixed limits mislead; use log-transformed data, ratios, or a regression-based LoA.
- **Plotting the difference against one method instead of the mean.** Against the reference value, an artificial correlation appears. Use the mean (or the reference only if it is error-free).
- **Not stating the direction of subtraction.** "Bias = 0.28" is ambiguous without "POC minus laboratory".
- **Treating repeated measures on the same patient as independent pairs.** This narrows the limits too much.

## Related tests

- [Intraclass correlation coefficient (ICC)](#/test/icc): reliability, with variance components
- [Paired t-test](#/test/paired-t): tests whether the mean difference (bias) is zero
- [Pearson correlation](#/test/pearson): association, not agreement
- [Cohen's kappa](#/test/cohens-kappa): agreement for categorical ratings
- [Equivalence (TOST) & non-inferiority tests](#/test/tost-equivalence): formal test that the mean difference lies within a margin

## Key takeaways

- Bland–Altman looks at the **differences** between methods: bias = mean difference, limits of agreement = bias ± 1.96 SD.
- High correlation does not mean agreement.
- The limits of agreement describe individual differences, and need their own CIs (wide when n is small).
- Decide the clinically acceptable difference **before** the study; the method gives numbers, not a verdict.
- Always state the direction of subtraction and check for trend or fanning on the plot.

## References

- Bland JM, Altman DG. Statistical methods for assessing agreement between two methods of clinical measurement. *Lancet.* 1986;i:307–310.
- Bland JM, Altman DG. Measuring agreement in method comparison studies. *Statistical Methods in Medical Research.* 1999;8:135–160.
- Kottner J, et al. Guidelines for Reporting Reliability and Agreement Studies (GRRAS). *Journal of Clinical Epidemiology.* 2011;64:96–106.
- SAS/STAT and ODS Graphics documentation: the SGPLOT procedure (REFLINE statement)
