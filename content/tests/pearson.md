# Pearson correlation

> Measures the strength and direction of the **linear** relationship between two continuous variables, from −1 (perfect negative) through 0 (none) to +1 (perfect positive).

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Strength and direction of a linear relationship between two variables | Two continuous variables | Independent pairs; roughly linear; no influential outliers | BMI and systolic blood pressure in 30 adults |

## Learning objectives

- Interpret r as a measure of **linear** association, with its sign and size
- Use a scatterplot to check linearity and outliers before trusting r
- Run Pearson correlation in R and Python and read r, 95% CI, t, and p
- Report r with its CI and avoid reading it as a slope, as agreement, or as causation

## What it is

Pearson's correlation coefficient r summarises how closely a scatterplot of two continuous measurements follows a straight line. It is the average product of the two variables' standardised (z-score) values, so it does not depend on the units: BMI in kg/m² and blood pressure in mmHg can be compared directly. Correlation describes association only; it says nothing about which variable influences the other.

- **Test statistic:** t = r √(n − 2) / √(1 − r²), with n − 2 degrees of freedom
- **Confidence interval:** from Fisher's z-transformation, z = arctanh(r), with standard error 1 / √(n − 3), back-transformed to the r scale
- **Effect size:** r itself with its 95% CI (primary); r² (the proportion of variance shared) as a secondary summary. As rough guides only, |r| ≈ 0.1, 0.3, and 0.5 are often called small, moderate, and large, but clinical context matters more than any cut-off.

## The intuition

Pearson's r answers: **when one variable is above its average, is the other usually above its average too?** Each person is turned into two z-scores (how many SDs from the mean for BMI, and for SBP). r is the average product of these. Pairs that are both high or both low add a positive amount; pairs that disagree subtract.

Signal vs noise: the signal is the tendency of the two variables to move together in a straight line. The noise is the scatter around that line. A tight band of points gives r near ±1; a shapeless cloud gives r near 0. The t statistic compares r with how much r would bounce around by chance with only n pairs.

We cannot draw a scatterplot here, but we can see the same trend by splitting the 30 adults into thirds by BMI and plotting SBP in each third. Systolic pressure climbs from lowest-BMI to higher-BMI groups (means about 128, 139, and 140 mmHg), yet the groups overlap a lot. That overlap is why r is moderate (0.56) rather than close to 1.

<div class="viz" data-viz="dots" data-unit="mmHg, SBP" data-groups='{"Lowest third of BMI":[133,141,120,133,116,121,125,134,134,127],"Middle third":[145,137,132,143,132,136,155,133,134,145],"Highest third":[137,140,143,139,146,131,138,139,140,148]}' data-caption="SBP by BMI third (BMI 19.8 to 33.2 kg/m2). Notice the upward drift in the means and the heavy overlap between groups."></div>

Our test gives t = 3.57 with 28 df. Under H₀ (ρ = 0), t follows the curve below, and ours is far out in the tail.

<div class="viz" data-viz="t" data-stat="3.57" data-df="28"></div>

**What matters most.** r only captures *straight-line* association, and a single extreme point can swing it. So the scatterplot checks (is it roughly linear? is any point driving the result?) matter more than a formal normality test. Normality of each variable is a minor concern, mostly for the CI and p-value in small samples.

## When to use it

Use Pearson correlation when **all** of these are true:

- Both variables are **continuous** (BMI, blood pressure, HbA1c, cholesterol)
- Each subject contributes **one pair** of values, and subjects are independent
- The relationship looks **roughly linear** on a scatterplot
- There are no extreme outliers. Approximate normality helps the CI and p-value in small samples, but linearity and outliers matter far more than normality

**Clinical scenario.** A cardiology clinic wants to know whether body mass index (BMI, kg/m²) tracks with resting systolic blood pressure (SBP, mmHg) in 30 adults without antihypertensive treatment. *How strong is the linear association, and how precisely is it estimated?*

**Don't use it when:**

- The relationship is monotonic but curved, the data are ordinal, or there are outliers → use [Spearman rank correlation](#/test/spearman)
- You want to **predict** SBP from BMI or adjust for other variables → use [linear regression](#/test/linear-regression)
- The two variables are two methods measuring the same quantity and you care about agreement → use [Bland–Altman analysis](#/test/bland-altman) or the [ICC](#/test/icc)
- The observations are repeated measures on the same patients → use a [mixed model](#/test/mixed-models)

## Assumptions

1. **Independent observations.** Each pair comes from a different person. Repeated measures or family clusters violate this.
2. **Linear relationship.** Inspect a scatterplot. A U-shaped or saturating pattern can give r near 0 despite a strong association.
3. **No influential outliers.** One extreme point can create or destroy a correlation. Re-compute r without it as a sensitivity check.
4. **Approximate bivariate normality** (only for the p-value and CI, mainly when n < 30). r itself is a valid descriptive summary of linear association without it. Check with Q–Q plots; a [Shapiro–Wilk test](#/test/shapiro-wilk) can support this but is not a gate. If strongly skewed, use Spearman or a bootstrap CI as a sensitivity analysis.
5. **Homoscedasticity.** The spread of one variable should be similar across the range of the other (the scatter should look like a band, not a funnel).
6. **No restriction of range.** Studying only a narrow slice of one variable (e.g., only obese patients) shrinks r.

## Hypotheses

- **H₀:** ρ = 0 (no linear association in the population)
- **H₁:** ρ ≠ 0 (two-sided)

## Example data

Simulated data for 30 adults; each row is one person. The first five rows:

| Patient | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| BMI (kg/m²) | 31.9 | 29.9 | 27.6 | 27.1 | 32.6 |
| SBP (mmHg) | 131 | 137 | 132 | 145 | 140 |

The full 30 pairs are typed inline in each code block below.

## How to run it

### R

```r
# BMI (kg/m²) and systolic BP (mmHg) in 30 adults
dat <- data.frame(
  bmi = c(31.9, 29.9, 27.6, 27.1, 32.6, 33.2, 30.1, 27.5, 23.1, 23.5, 21.3,
      29.9, 31.3, 29.2, 22.9, 27.4, 24.3, 29.8, 21.6, 20, 32.6, 29.3, 32.2,
      27.1, 29.6, 21.8, 27.5, 30.7, 19.8, 27.6),
  sbp = c(131, 137, 132, 145, 140, 148, 143, 132, 125, 134, 120, 140, 146,
      155, 121, 137, 134, 145, 133, 141, 139, 133, 138, 127, 134, 116, 143,
      139, 133, 136)
)

# 1. Look at the data and check linearity / outliers
plot(dat$bmi, dat$sbp, xlab = "BMI", ylab = "SBP")

# 2. Normality of each variable (bivariate normality matters for the CI/p-value)
shapiro.test(dat$bmi); shapiro.test(dat$sbp)

# 3. Pearson correlation with 95% CI
cor.test(dat$bmi, dat$sbp, method = "pearson")

# 4. Effect size: r² = share of SBP variance shared with BMI
cor(dat$bmi, dat$sbp)^2
```

**What the code does**

- `data.frame(...)` holds the 30 pairs; BMI and SBP in the same position belong to the same person.
- `plot()` gives the scatterplot. Always look before you compute r.
- `shapiro.test()` is run on each variable. It tests normality separately; it is a screening aid rather than a proof of bivariate normality.
- `cor.test(..., method = "pearson")` returns r, t, df, the two-sided p-value, and a **95% CI based on Fisher's z**. `method = "pearson"` is the default, but writing it out avoids confusion with `"spearman"`.
- `cor()^2` gives r², the proportion of variance shared.

**Output**

```text

	Shapiro-Wilk normality test

data:  dat$bmi
W = 0.92562, p-value = 0.03762


	Shapiro-Wilk normality test

data:  dat$sbp
W = 0.97701, p-value = 0.7417


	Pearson's product-moment correlation

data:  dat$bmi and dat$sbp
t = 3.5664, df = 28, p-value = 0.001326
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
 0.2487092 0.7651106
sample estimates:
      cor 
0.5588997 

[1] 0.3123689
```

### Python

```python
import numpy as np, pandas as pd
from scipy import stats

# BMI (kg/m²) and systolic BP (mmHg) in 30 adults
dat = pd.DataFrame({
  'bmi': [31.9, 29.9, 27.6, 27.1, 32.6, 33.2, 30.1, 27.5, 23.1, 23.5, 21.3,
      29.9, 31.3, 29.2, 22.9, 27.4, 24.3, 29.8, 21.6, 20, 32.6, 29.3, 32.2,
      27.1, 29.6, 21.8, 27.5, 30.7, 19.8, 27.6],
  'sbp': [131, 137, 132, 145, 140, 148, 143, 132, 125, 134, 120, 140, 146,
      155, 121, 137, 134, 145, 133, 141, 139, 133, 138, 127, 134, 116, 143,
      139, 133, 136]
})

# 1. Normality of each variable
print("Shapiro BMI: p = %.3f; SBP: p = %.3f" %
      (stats.shapiro(dat.bmi).pvalue, stats.shapiro(dat.sbp).pvalue))

# 2. Pearson r, p-value and 95% CI (CI needs SciPy >= 1.9)
res = stats.pearsonr(dat.bmi, dat.sbp)
ci = res.confidence_interval(confidence_level=0.95)
print(f"r = {res.statistic:.3f}, p = {res.pvalue:.2e}, df = {len(dat) - 2}")
print(f"95% CI ({ci.low:.3f}, {ci.high:.3f})")

# 3. Effect size: r squared
print(f"r^2 = {res.statistic**2:.3f}")
```

**What the code does**

- `stats.pearsonr` returns r and the two-sided p-value. In SciPy ≥ 1.9 the result object also has `.confidence_interval()` (Fisher z, same as R).
- `stats.shapiro(...)` checks normality of each variable.
- Degrees of freedom are not returned by `pearsonr`; they are n − 2.
- r² is just `r ** 2`.

**Output**

```text
Shapiro BMI: p = 0.038; SBP: p = 0.742
r = 0.559, p = 1.33e-03, df = 28
95% CI (0.249, 0.765)
r^2 = 0.312
```

### SAS

```sas
/* BMI (kg/m2) and systolic BP (mmHg) in 30 adults */
data bpbmi;
  input bmi sbp;
  datalines;
  31.9 131
  29.9 137
  27.6 132
  27.1 145
  32.6 140
  33.2 148
  30.1 143
  27.5 132
  23.1 125
  23.5 134
  21.3 120
  29.9 140
  31.3 146
  29.2 155
  22.9 121
  27.4 137
  24.3 134
  29.8 145
  21.6 133
  20 141
  32.6 139
  29.3 133
  32.2 138
  27.1 127
  29.6 134
  21.8 116
  27.5 143
  30.7 139
  19.8 133
  27.6 136
;
run;

/* Pearson r with 95% CI (Fisher z) and a scatterplot */
proc corr data=bpbmi pearson fisher(biasadj=no alpha=0.05) plots=scatter(alpha=.05);
  var bmi sbp;
run;

/* Normality of each variable */
proc univariate data=bpbmi normal;
  var bmi sbp;
run;
```

**What the code does**

- `proc corr ... pearson` is the default method; it is spelled out for clarity.
- `fisher(biasadj=no alpha=0.05)` adds the **Fisher z-based 95% CI** and a p-value for H₀: ρ = 0. By default (`biasadj=yes`) SAS centres the interval on a bias-adjusted estimate, which differs slightly from r and from R and SciPy; `biasadj=no` matches them. Without the `fisher` option, PROC CORR prints only r and the p-value.
- `plots=scatter(alpha=.05)` draws the scatterplot with a prediction ellipse (requires ODS Graphics).
- `proc univariate ... normal` prints Shapiro–Wilk and other normality tests for each variable.
- Missing values are dropped pairwise by default; add `nomiss` to use only complete rows.

## Interpreting the output

- **r = 0.559.** A positive, moderate-to-strong linear association: patients with higher BMI tend to have higher SBP.
- **95% CI: 0.25 to 0.77.** The data are compatible with anything from a weak-to-moderate to a strong correlation. The interval is wide because n is only 30; plan sample size for the precision you need, not just for p < 0.05.
- **t(28) = 3.57, p = 0.0013.** There is strong evidence that ρ ≠ 0. The p-value says the correlation is unlikely to be exactly zero; it does not say it is large.
- **r² = 0.31.** About 31% of the variation in SBP in this sample is linearly associated with BMI; the remaining 69% is due to other factors and noise.
- **Shapiro–Wilk:** SBP shows no evidence against normality (p = 0.74), while BMI shows some evidence of non-normality (p = 0.038). With n = 30 this is a mild departure that matters less than linearity and outliers. Check the scatterplot and Q–Q plot, and consider a [Spearman](#/test/spearman) analysis as a sensitivity check.

**How to report it:**
> BMI and systolic blood pressure were positively correlated (Pearson r = 0.56, 95% CI 0.25 to 0.77, n = 30, p = 0.001).

⚠️ Correlation is not causation, and r does not tell you the **slope**. To state how many mmHg SBP rises per 1 kg/m² of BMI, fit a [linear regression](#/test/linear-regression).

## Check your understanding

<details class="check"><summary>1. A study finds r = 0.02 between age and a biomarker, but the scatterplot shows a clear U-shape. What does r = 0.02 tell you?</summary>

Only that there is no *straight-line* trend. The U-shape is a strong, non-linear association. Always plot the data first; consider a non-linear model, not a correlation.
</details>

<details class="check"><summary>2. Shapiro–Wilk on BMI gives p = 0.038. Does that invalidate Pearson's r?</summary>

Not by itself. With n = 30 it signals mild non-normality. Look at the Q–Q plot and scatterplot, check for outliers and curvature (more important for r), and run [Spearman](#/test/spearman) or a bootstrap CI as a sensitivity analysis. Do not switch methods solely because a normality p-value crossed 0.05.
</details>

<details class="check"><summary>3. Two lab methods measure creatinine, and r = 0.96. A colleague concludes the methods agree. Is that justified?</summary>

No. High correlation means the values move together, but one method could read 20% higher everywhere and still give r = 0.96. Use [Bland–Altman analysis](#/test/bland-altman) or the [ICC](#/test/icc) to judge agreement.
</details>

<details class="check"><summary>4. What if we studied only obese patients (BMI 30 to 33)? What would likely happen to r?</summary>

r would shrink, a **restriction of range** effect. With a narrow BMI span, the BMI-related variation in SBP is small compared with other sources of variation, so the scatter looks less linear even though the underlying relationship has not changed.
</details>

## Common mistakes

- **Not plotting the data.** Very different patterns (curves, clusters, one outlier) can share the same r.
- **Reading r as a percentage or as the slope.** r = 0.56 does not mean "56% agreement" or "SBP rises 0.56 mmHg per unit".
- **Reporting only p.** A tiny r can be "significant" in a huge sample and a large r "non-significant" in a small one. Give r and its CI.
- **Pooling groups.** Combining men and women, or treated and untreated patients, can create or hide a correlation (Simpson's paradox).
- **Using correlation to assess agreement.** Two methods can correlate at r = 0.95 yet differ systematically; use [Bland–Altman](#/test/bland-altman).
- **Running many correlations and reporting the significant ones.** Pre-specify the pairs, or apply a [multiple-comparison correction](#/test/multiple-comparisons).

## Related tests

- [Spearman rank correlation](#/test/spearman): rank-based, robust to outliers and curved monotonic trends
- [Linear regression](#/test/linear-regression): slope, prediction, and adjustment for covariates
- [Bland–Altman analysis](#/test/bland-altman): agreement between two measurement methods
- [Intraclass correlation coefficient (ICC)](#/test/icc): reliability and agreement
- [Shapiro–Wilk test](#/test/shapiro-wilk): normality check

## Key takeaways

- Pearson r (−1 to +1) measures **linear** association; it is unit-free.
- Always look at the scatterplot first. Curvature, clusters, and outliers matter more than a normality test.
- Report r with its 95% CI. r² is the share of variance linearly shared, and r is not the slope.
- Correlation does not mean causation or agreement (use regression, Bland–Altman, or ICC for those).
- Restricted range shrinks r; pooled subgroups can create or hide it.

## References

- Altman DG. *Practical Statistics for Medical Research.* Chapman & Hall; 1991. Ch. 11.
- Bland JM, Altman DG. Statistics notes: Correlation, regression, and repeated data. *BMJ.* 1994;308:896.
- R documentation: `?cor.test`
- SciPy documentation: `scipy.stats.pearsonr`; SAS/STAT User's Guide: The CORR Procedure (FISHER option)
