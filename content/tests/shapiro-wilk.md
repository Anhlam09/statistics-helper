# Shapiro–Wilk test

> Tests whether a sample could plausibly have come from a normal distribution, by checking how closely the ordered data line up with the straight line expected under normality.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Screen one sample (or model residuals) for departure from normality | Continuous, ideally few ties | Independent observations; few tied values | Fasting triglycerides in 12 patients |

## Learning objectives

- Explain what W measures and why a **small** p-value signals non-normality
- Run Shapiro–Wilk in R and Python and read W and p alongside a Q–Q plot
- Recognize why p-values from this test depend heavily on sample size
- Avoid using the test as an automatic gatekeeper between parametric and nonparametric tests

## What it is

The Shapiro–Wilk test compares your sorted observations with the values you would expect from a normal distribution with the same mean and spread. If the points track the expected values closely, W is near 1. If the data are skewed, heavy-tailed, or contain outliers, W drops and the p-value shrinks.

- **Test statistic:** W, between 0 and 1. Values near 1 indicate approximate normality; small values indicate departure from it.
- **Sample size:** works for 3 ≤ n ≤ 5000 in R and SciPy. It is most informative for small to moderate n.
- **Effect size:** there is no conventional one. Report W, the p-value, and a **Q–Q plot**, and describe the departure (skewness, outliers) rather than relying on the p-value alone.

## The intuition

Think of W as a **fit score** between your sorted data and the straight line a normal sample would follow. There is no "signal vs noise" here: the question is how far your data bend away from a bell shape. W = 1 means a perfect match; the more the data are skewed or have outliers, the lower W falls.

An everyday analogy: you lay a ruler along a row of fence posts. If the posts follow a straight line, the ruler fits. If a few posts at one end stick out, the ruler no longer fits. W is a number summarizing how well the ruler fits.

The dot plot below shows our 12 triglyceride values. Most patients sit between 88 and 172 mg/dL, but a few stretch far to the right. That long right tail is what pulls W down to 0.83.

<div class="viz" data-viz="dots" data-unit="mg/dL" data-groups='{"Raw TG":[88,95,102,110,118,125,134,150,172,210,265,340]}' data-caption="Raw triglycerides: a cluster at the low end and a long tail to the right (skewed)."></div>

After a log transformation the same patients look much more evenly spread, with no isolated far-away points. W rises to 0.93.

<div class="viz" data-viz="dots" data-unit="natural log of TG" data-groups='{"Log TG":[4.48,4.55,4.62,4.7,4.77,4.83,4.9,5.01,5.15,5.35,5.58,5.83]}' data-caption="Same patients on the log scale: the right tail has largely disappeared."></div>

**Sample size matters.** With a handful of patients, only dramatic departures are detected (low power). With thousands, even trivial, harmless bumps give p < 0.001. So the p-value tells you how *detectable* a departure is, not how *important* it is. Always look at the shape.

## When to use it

Use the Shapiro–Wilk test when:

- You want a quick, objective **supporting check for normality** alongside a Q–Q plot, for a method that assumes it (t-tests, ANOVA, linear regression)
- The variable is **continuous**
- The sample is small to moderate, where visual checks alone are hard to judge

**Clinical scenario.** A lipid clinic records fasting triglycerides (mg/dL) in 12 patients and plans to compare them with a reference value using a t-test. Triglycerides are well known to be right-skewed. *Is it reasonable to treat the raw values as normal, or should they be log-transformed first?*

**Don't use it when:**

- You want to **prove** normality. A non-significant result only means the test did not detect a departure; with small n it has little power. Combine it with a Q–Q plot.
- The sample is **very large** (hundreds or more). The test flags tiny, clinically irrelevant departures. Rely on the Q–Q plot and the robustness of the planned test.
- You are testing the **raw outcome** in a model with several groups or covariates. The assumption applies to the residuals; test those instead (see [one-way ANOVA](#/test/one-way-anova) and [linear regression](#/test/linear-regression)).
- Your question is about equal variances → use [Levene's test](#/test/levene)

## Assumptions

1. **Observations are independent.** Repeated measures on the same patient violate this. Check the study design.
2. **The variable is continuous**, with few ties. Many tied or rounded values distort W. Check the number of distinct values.
3. **The data are a random sample** from the population of interest, with no gross data-entry errors. Check the minimum and maximum before testing.

## Hypotheses

- **H₀:** the data come from a normal distribution
- **H₁:** the data do not come from a normal distribution

Note the direction: a **small** p-value is evidence of non-normality, which is the opposite of how most tests are read.

## Example data

Fasting triglycerides (mg/dL) in 12 patients, sorted:

| Patient | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| TG (mg/dL) | 88 | 95 | 102 | 110 | 118 | 125 | 134 | 150 | 172 | 210 | 265 | 340 |

## How to run it

### R

```r
# Fasting triglycerides (mg/dL), 12 patients
tg <- c(88, 95, 102, 110, 118, 125, 134, 150, 172, 210, 265, 340)

# Shapiro-Wilk on the raw values
shapiro.test(tg)

# Shapiro-Wilk after a log transformation
shapiro.test(log(tg))

# Visual check: Q-Q plot
qqnorm(tg); qqline(tg)
```

**What the code does**

- `shapiro.test(tg)` returns W and the p-value. It needs between 3 and 5000 non-missing values.
- `log(tg)` is the natural log. Right-skewed lab values often become roughly symmetric on the log scale.
- `qqnorm()` and `qqline()` draw a Q–Q plot. Points that bend away from the line, especially at the upper end, show right skew.

**Output**

```text
	Shapiro-Wilk normality test

data:  tg
W = 0.83327, p-value = 0.02296


	Shapiro-Wilk normality test

data:  log(tg)
W = 0.92584, p-value = 0.3381
```

### Python

```python
import numpy as np
from scipy import stats

# Fasting triglycerides (mg/dL), 12 patients
tg = np.array([88, 95, 102, 110, 118, 125, 134, 150, 172, 210, 265, 340])

# Shapiro-Wilk on the raw values
raw = stats.shapiro(tg)
print(f"Raw:  W = {raw.statistic:.4f}, p = {raw.pvalue:.4f}")

# Shapiro-Wilk after a log transformation
lg = stats.shapiro(np.log(tg))
print(f"Log:  W = {lg.statistic:.4f}, p = {lg.pvalue:.4f}")

# Skewness as a descriptive companion
print(f"Skewness (raw) = {stats.skew(tg, bias=False):.2f}")
```

**What the code does**

- `stats.shapiro` returns an object with `.statistic` (W) and `.pvalue`. It gives the same W and p as R's `shapiro.test`.
- `np.log` is the natural log, matching R's `log()`.
- `stats.skew(..., bias=False)` gives the sample skewness with the small-sample correction. A value well above 1 means a pronounced right tail. The default `bias=True` gives a slightly smaller number.

**Output**

```text
Raw:  W = 0.8333, p = 0.0230
Log:  W = 0.9258, p = 0.3381
Skewness (raw) = 1.48
```

### SAS

```sas
data tg;
  input tg @@;
  logtg = log(tg);
  datalines;
88 95 102 110 118 125 134 150 172 210 265 340
;
run;

/* NORMAL requests Shapiro-Wilk and other normality tests */
proc univariate data=tg normal;
  var tg logtg;
  qqplot tg logtg / normal(mu=est sigma=est);
run;
```

**What the code does**

- `input tg @@;` reads several values per line, which keeps the data compact.
- `logtg = log(tg);` creates the natural-log variable in the DATA step.
- `normal` in PROC UNIVARIATE prints the "Tests for Normality" table: Shapiro–Wilk, Kolmogorov–Smirnov, Cramér–von Mises, and Anderson–Darling. Use the Shapiro–Wilk line (W and Pr < W); it matches R and Python.
- `qqplot ... / normal(mu=est sigma=est)` draws a normal Q–Q plot with a reference line from the sample mean and SD.

## Interpreting the output

- **Raw TG: W = 0.833, p = 0.023.** There is evidence against normality at the 5% level. The sample skewness of 1.48 confirms a long right tail (a few patients with very high TG).
- **Log TG: W = 0.926, p = 0.338.** After the log transformation there is no evidence against normality. This does not prove the log data are normal; with n = 12 the test has limited power.
- **R and Python agree** to all printed digits.
- **What to do next.** Analyse log-TG with a parametric test and back-transform the results (a difference of logs becomes a ratio of geometric means), or use a rank-based test such as the [Wilcoxon signed-rank test](#/test/wilcoxon-signed-rank) or [Mann–Whitney U test](#/test/mann-whitney).

**How to report it:**
> Fasting triglycerides were right-skewed (Shapiro–Wilk W = 0.83, p = 0.023) and were therefore analysed on the natural-log scale, where no departure from normality was detected (W = 0.93, p = 0.34).

## Check your understanding

<details class="check"><summary>1. You have 400 patients and Shapiro–Wilk gives p = 0.002, but the histogram looks almost bell-shaped. Should you abandon the t-test?</summary>

Not automatically. With n = 400 the test detects tiny departures that barely affect a t-test, which is robust at that size thanks to the central limit theorem. Trust the Q–Q plot and the size of the departure (strong skew, heavy outliers) more than the p-value.
</details>

<details class="check"><summary>2. A colleague says "Shapiro–Wilk gave p = 0.4, so my data are normal and I can use the t-test." What is wrong with that?</summary>

A large p-value means no departure was *detected*, not that normality is proven. With small n the test has little power. Also, choosing the test after looking at the data distorts error rates. Pre-specify the analysis plan, and use the Q–Q plot and subject-matter knowledge (lab values are often right-skewed) to justify the choice.
</details>

<details class="check"><summary>3. You compare a biomarker between three treatment arms. Where should you apply Shapiro–Wilk?</summary>

To the **residuals** of the ANOVA model (each value minus its group mean), or within each group if groups are large enough. Pooling the raw values across arms mixes group differences with within-group shape and can look non-normal even when the model assumption holds.
</details>

<details class="check"><summary>4. What if the largest triglyceride value (340 mg/dL) were 190 mg/dL instead? What would you expect to happen to W?</summary>

W would move closer to 1 and the p-value would increase. The extreme outlier is what drives the long right tail; remove it and the sorted data track the normal reference line much better.
</details>

## Common mistakes

- **Treating p > 0.05 as proof of normality.** It only means no departure was detected. With n = 12 only large departures are detectable.
- **Testing the wrong thing.** For ANOVA and regression, test the **residuals**, not the pooled raw outcome. For a paired t-test, test the **differences**.
- **Using it as a hard gatekeeper.** Choosing between a t-test and a nonparametric test based on a Shapiro–Wilk result distorts the overall error rate of the two-stage procedure, and the test is least reliable exactly when n is small. Pre-specify the plan in the statistical analysis plan (ICH E9), and use the test only as supporting evidence next to a Q–Q plot and subject-matter knowledge.
- **Ignoring the plot.** Always look at a Q–Q plot or histogram. It tells you *how* the data depart from normal.
- **Over-reacting with large n.** In a big sample, trivial skewness gives p < 0.001 even though the t-test is robust.

## Related tests

- [Levene's test](#/test/levene): checks the other key ANOVA assumption, equal variances
- [One-sample t-test](#/test/one-sample-t): assumes normality of the data
- [Independent two-sample t-test](#/test/independent-t): assumes normality in each group
- [One-way ANOVA](#/test/one-way-anova): assumes normal residuals
- [Wilcoxon signed-rank test](#/test/wilcoxon-signed-rank) and [Mann–Whitney U test](#/test/mann-whitney): rank-based alternatives

## Key takeaways

- Shapiro–Wilk gives W (near 1 = looks normal). A **small** p-value is evidence of non-normality, the reverse of how most tests are read.
- p > 0.05 does not prove normality; with small n the test has little power.
- With large n it flags trivial departures. Judge the **size** of the departure with a Q–Q plot.
- Test the **residuals** (ANOVA, regression) or the **differences** (paired designs), not the pooled raw outcome.
- Do not use it as an automatic switch between a t-test and a nonparametric test. Pre-specify the plan.

## References

- Shapiro SS, Wilk MB. An analysis of variance test for normality (complete samples). *Biometrika.* 1965;52:591–611.
- Altman DG. *Practical Statistics for Medical Research.* Chapman & Hall; 1991. Ch. 4.
- R documentation: `?shapiro.test`
- SciPy documentation: `scipy.stats.shapiro`; SAS/STAT User's Guide: The UNIVARIATE Procedure — NORMAL option
