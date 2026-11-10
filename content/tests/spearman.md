# Spearman rank correlation

> Measures the strength and direction of a **monotonic** relationship (as one variable rises, the other consistently rises or falls) by correlating the *ranks* of two variables, so it tolerates skewed data, outliers, and ordinal scales.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Measure monotonic association between two variables, robust to skew and outliers | Continuous or ordinal, paired observations | Independent subjects; monotonic relationship; few ties | CRP vs 6-minute walk distance in COPD |

## Learning objectives

- Explain how ranking turns a correlation into a measure of monotonic association
- Run Spearman correlation in R and Python and report ρ with its 95% CI
- Distinguish monotonic from linear relationships, and know when ρ ≈ 0 is misleading
- Avoid reading ρ as a slope, a size of effect, or agreement

## What it is

Spearman's ρ (rho) replaces each value by its rank within its own variable (smallest = 1) and then computes the Pearson correlation of the ranks. Because only the ordering matters, any increasing transformation of the data (log, square root) leaves ρ unchanged, and a single extreme value can move a point by at most a few ranks. Ties receive the average of the ranks they span.

- **Test statistic:** without ties, the classic form is ρ = 1 − 6 Σd² / (n(n² − 1)), where d is the rank difference for each subject. The p-value comes from the exact permutation distribution (small n, no ties) or from t = ρ √(n − 2) / √(1 − ρ²) with n − 2 df (large n, ties).
- **Confidence interval:** approximate, from Fisher's z-transformation with a standard error of √((1 + ρ²/2) / (n − 3)) (Bonett–Wright); a bootstrap CI is a good alternative
- **Effect size:** ρ with its 95% CI. ρ² is the share of *rank* variance explained, but it is less often reported than for Pearson r.

## The intuition

Suppose you line up all 24 patients twice: once from lowest to highest CRP, and once from shortest to longest walk. If inflammation and walking ability go together, the two line-ups should be mirror images, with the patient at the top of one list near the bottom of the other. ρ measures how closely the two line-ups agree, from −1 (perfect mirror) through 0 (unrelated) to +1 (identical). Because only each patient's *position in line* matters, one patient with a CRP of 200 counts only as "highest", not as a number that can drag a whole fit around.

- **Signal:** how consistently high ranks on one variable pair with low ranks on the other.
- **Noise:** how much the rank pairings would scatter if the two variables were unrelated.

First, see the problem that ranks solve. CRP is strongly right-skewed: most patients cluster below 10 mg/L and a few trail out to about 20.

<div class="viz" data-viz="dots" data-unit="mg/L, CRP" data-groups='{"CRP":[18.8,9.7,7.3,5.2,19.7,4,0.7,8.5,7,0.7,12.4,1.5,5.8,4.8,5.3,3.4,16.1,8.3,1.8,4.5,1.5,10.9,9.8,0.9]}' data-caption="CRP values (n = 24; mean 7.0, median 5.6). The skew is why a rank method is preferred."></div>

For example, patient 5 has the highest CRP (rank 24) and the lowest walk distance (rank 1), the extreme of a negative association. Patient 6 has CRP rank 8 and walk rank 19. Across all patients, ρ = −0.67.

To test it, convert ρ to a t statistic: t = ρ √(n − 2) / √(1 − ρ²) = −4.28 on 22 degrees of freedom. If there were no monotonic association, t would follow the curve below. Ours is far into the tail.

<div class="viz" data-viz="t" data-stat="-4.28" data-df="22"></div>

What to notice: the shaded area is tiny (p = 0.0003). But p only says the association is unlikely to be zero. The width of the CI (−0.86 to −0.33) tells you how well the strength is pinned down.

## When to use it

Use Spearman correlation when:

- Both variables are at least **ordinal** (continuous, counts, or ordered categories such as a 0–10 pain score)
- The relationship is **monotonic**, though not necessarily a straight line
- The data are **skewed**, contain outliers, or the sample is small and clearly non-normal
- Each subject contributes one pair, and subjects are independent

**Clinical scenario.** In 24 patients with stable COPD, serum C-reactive protein (CRP, mg/L) is measured together with the 6-minute walk distance (6MWD, metres). CRP is strongly right-skewed. *Is higher systemic inflammation associated with a shorter walking distance?*

**Don't use it when:**

- Both variables are roughly normal and the relationship is linear, and you want the most efficient estimate → use [Pearson correlation](#/test/pearson)
- The relationship is not monotonic (e.g., U-shaped, such as J-shaped mortality versus BMI) → rank correlation will be near 0; model the shape with [linear regression](#/test/linear-regression) with a non-linear term
- You need a slope, prediction, or covariate adjustment → use [linear regression](#/test/linear-regression)
- The question is agreement between two measurement methods → use [Bland–Altman analysis](#/test/bland-altman) or the [ICC](#/test/icc)
- You have **many ties** in two short ordinal scales (e.g., a 3-point scale) → consider Kendall's τ-b (`cor.test(method = "kendall")`)

## Assumptions

1. **Independent observations.** Each pair belongs to a different subject.
2. **Monotonic relationship.** Check with a scatterplot (or a plot of the ranks). The relationship may be curved but should not change direction.
3. **At least ordinal measurement.** The ordering of values must be meaningful.
4. **Few ties.** Heavy ties reduce the usable information and break exact p-values; software then applies a correction or an asymptotic approximation.
5. **No distributional assumption.** Unlike Pearson, normality is *not* required. However, ρ is not a measure of agreement and, like any correlation, is affected by restriction of range.

## Hypotheses

- **H₀:** ρₛ = 0 (no monotonic association in the population)
- **H₁:** ρₛ ≠ 0 (two-sided)

## Example data

Simulated data for 24 COPD patients. The first six rows:

| Patient | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| CRP (mg/L) | 18.8 | 9.7 | 7.3 | 5.2 | 19.7 | 4.0 |
| 6MWD (m) | 353 | 417 | 384 | 397 | 337 | 446 |

The full 24 pairs are typed inline in each code block below. Several CRP and walk values are tied (for example, CRP 0.7 and 1.5 each appear twice).

## How to run it

### R

```r
# CRP (mg/L, right-skewed) and 6-minute walk distance (m) in 24 COPD patients
dat <- data.frame(
  crp = c(18.8, 9.7, 7.3, 5.2, 19.7, 4, 0.7, 8.5, 7, 0.7, 12.4, 1.5, 5.8, 4.8,
      5.3, 3.4, 16.1, 8.3, 1.8, 4.5, 1.5, 10.9, 9.8, 0.9),
  walk = c(353, 417, 384, 397, 337, 446, 453, 365, 368, 520, 386, 389, 360,
      368, 404, 465, 432, 400, 490, 417, 417, 386, 370, 525)
)

# 1. Show why Pearson is a poor fit: CRP is skewed
shapiro.test(dat$crp)

# 2. Spearman rho (exact = FALSE because there are tied values)
fit <- cor.test(dat$crp, dat$walk, method = "spearman", exact = FALSE)
fit

# 3. Approximate 95% CI (Fisher z, Bonett-Wright standard error)
rho <- unname(fit$estimate); n <- nrow(dat)
se  <- sqrt((1 + rho^2 / 2) / (n - 3))
tanh(atanh(rho) + c(-1, 1) * qnorm(0.975) * se)

# 4. For comparison: Pearson on the raw (skewed) values
cor(dat$crp, dat$walk)
```

**What the code does**

- `cor.test(..., method = "spearman")` computes ρ and its p-value. The default is `"pearson"`, so forgetting the `method` argument silently gives a different analysis.
- `exact = FALSE` requests the asymptotic p-value. With tied values R cannot compute an exact p-value and would otherwise print a warning.
- `cor.test` does **not** give a CI for Spearman, so the CI is built by hand: `atanh(rho)` is Fisher's z, the standard error uses the Bonett–Wright formula, and `tanh()` converts back to the ρ scale.
- The last line shows the Pearson coefficient on the raw values for contrast.

**Output**

```text

	Shapiro-Wilk normality test

data:  dat$crp
W = 0.90432, p-value = 0.0266


	Spearman's rank correlation rho

data:  dat$crp and dat$walk
S = 3849.7, p-value = 0.0003068
alternative hypothesis: true rho is not equal to 0
sample estimates:
       rho 
-0.6737807 

[1] -0.8594917 -0.3309321
[1] -0.6160618
```

### Python

```python
import numpy as np, pandas as pd
from scipy import stats

# CRP (mg/L, right-skewed) and 6-minute walk distance (m) in 24 COPD patients
dat = pd.DataFrame({
  'crp': [18.8, 9.7, 7.3, 5.2, 19.7, 4, 0.7, 8.5, 7, 0.7, 12.4, 1.5, 5.8, 4.8,
      5.3, 3.4, 16.1, 8.3, 1.8, 4.5, 1.5, 10.9, 9.8, 0.9],
  'walk': [353, 417, 384, 397, 337, 446, 453, 365, 368, 520, 386, 389, 360,
      368, 404, 465, 432, 400, 490, 417, 417, 386, 370, 525]
})

# 1. CRP is skewed
print(f"Shapiro CRP: W = {stats.shapiro(dat.crp).statistic:.4f}, "
      f"p = {stats.shapiro(dat.crp).pvalue:.5f}")

# 2. Spearman rho and p-value (asymptotic t-approximation, handles ties)
res = stats.spearmanr(dat.crp, dat.walk)
n = len(dat)
print(f"rho = {res.statistic:.4f}, p = {res.pvalue:.5f}, n = {n}")

# 3. Approximate 95% CI: Fisher z with Bonett-Wright standard error
rho = res.statistic
se = np.sqrt((1 + rho**2 / 2) / (n - 3))
lo, hi = np.tanh(np.arctanh(rho) + np.array([-1, 1]) * 1.959964 * se)
print(f"95% CI ({lo:.4f}, {hi:.4f})")

# 4. For comparison: Pearson on the raw (skewed) values
print(f"Pearson r = {stats.pearsonr(dat.crp, dat.walk).statistic:.4f}")
```

**What the code does**

- `stats.spearmanr(x, y)` returns ρ and a two-sided p-value. With ties it uses the t-distribution approximation, which matches R's `exact = FALSE`.
- SciPy has no built-in CI for ρ, so the same Fisher-z / Bonett–Wright interval as in R is computed by hand (1.959964 is the 97.5th percentile of the standard normal).
- `stats.shapiro()` documents the skewness of CRP.
- `stats.pearsonr` is shown only for comparison. Argument order does not matter for either coefficient.

**Output**

```text
Shapiro CRP: W = 0.9043, p = 0.02660
rho = -0.6738, p = 0.00031, n = 24
95% CI (-0.8595, -0.3309)
Pearson r = -0.6161
```

### SAS

```sas
/* CRP (mg/L) and 6-minute walk distance (m) in 24 COPD patients */
data copd;
  input crp walk;
  datalines;
  18.8 353
  9.7 417
  7.3 384
  5.2 397
  19.7 337
  4 446
  0.7 453
  8.5 365
  7 368
  0.7 520
  12.4 386
  1.5 389
  5.8 360
  4.8 368
  5.3 404
  3.4 465
  16.1 432
  8.3 400
  1.8 490
  4.5 417
  1.5 417
  10.9 386
  9.8 370
  0.9 525
;
run;

/* Spearman rho; FISHER adds a CI and p-value based on Fisher's z.
   Pearson is requested too, for comparison. */
proc corr data=copd spearman pearson fisher(alpha=0.05)
          plots=scatter;
  var crp walk;
run;

/* Skewness check of CRP */
proc univariate data=copd normal;
  var crp;
  histogram crp / normal;
run;
```

**What the code does**

- `proc corr ... spearman pearson` prints both coefficients. Spearman ranks the data internally; the Spearman p-value is based on a t approximation.
- `fisher(alpha=0.05)` adds a Fisher-z confidence interval and p-value, labelled separately for the Pearson and Spearman coefficients in recent SAS/STAT releases (check your version's documentation for rank-based coefficients). SAS uses its own variance approximation for ρ and, by default, a bias-adjusted centre (add `biasadj=no` to turn that off), so its interval can differ slightly from the R and Python intervals above, while ρ and the p-value agree.
- `plots=scatter` draws the scatterplot (requires ODS Graphics).
- `proc univariate ... normal` and the histogram document the skewness of CRP.

## Interpreting the output

- **ρ = −0.67.** A moderate-to-strong negative monotonic association: patients with higher CRP tended to walk shorter distances.
- **95% CI: −0.86 to −0.33.** The data are compatible with a true rank correlation from moderate to very strong. The interval is wide because n is only 24.
- **p = 0.0003.** Strong evidence against no monotonic association. It does not tell you how large the association is; use ρ and the CI for that.
- **Shapiro–Wilk for CRP: W = 0.904, p = 0.027.** CRP departs from normality, which is why a rank method is preferred here.
- **Pearson r = −0.62 on the raw data.** It is somewhat weaker than ρ, as is typical when one variable is skewed. The two methods agree on the direction and approximate strength, which is reassuring.
- ρ says nothing about the *size* of the effect in metres. To estimate, say, the walk-distance change per doubling of CRP, fit a [linear regression](#/test/linear-regression) on log(CRP).

**How to report it:**
> CRP was inversely correlated with 6-minute walk distance (Spearman ρ = −0.67, 95% CI −0.86 to −0.33, n = 24, p < 0.001).

## Check your understanding

<details class="check"><summary>1. A study finds ρ = 0.02 between BMI and mortality across a wide BMI range. Does BMI have no relationship with mortality?</summary>

Not necessarily. Spearman only detects **monotonic** relationships. Mortality is often J-shaped or U-shaped against BMI, with higher risk at both extremes. Those opposing trends can cancel out and give ρ ≈ 0. Always look at a scatterplot.
</details>

<details class="check"><summary>2. Both variables are roughly normal and the scatterplot is a clean straight line. Should you still prefer Spearman because it is "safer"?</summary>

Not really. Pearson is more efficient for linear, normal data. Spearman is still valid, but it throws away information about the actual distances between values. Choose based on the shape of the relationship and the variable types, and say so in advance.
</details>

<details class="check"><summary>3. Applying a log transform to CRP before computing ρ gives −0.67 again. Why?</summary>

Spearman uses only ranks, and the log is an increasing transformation, which does not change any ranks. Pearson r *would* change, which is why log-transforming skewed variables matters for Pearson but not for Spearman.
</details>

<details class="check"><summary>4. A colleague reports ρ = −0.67 and says "each doubling of CRP shortens the walk by 67 metres". What is wrong, and what would you use instead?</summary>

ρ is a unit-free measure of how consistently ranks go in opposite directions. It is not a slope and carries no information about metres. To estimate the change in walk distance per doubling of CRP, fit a [linear regression](#/test/linear-regression) of walk distance on log₂(CRP) and report the slope with its CI.
</details>

## Common mistakes

- **Using Spearman for a non-monotonic relationship.** A U-shaped pattern gives ρ ≈ 0 even when the association is strong.
- **Thinking Spearman is "Pearson but safer" for everything.** When data are normal and linear, Pearson is more efficient; with a non-linear but monotonic trend, Spearman is the better descriptor.
- **Ignoring ties.** Exact p-values are not valid with ties; make sure your software uses the corrected or asymptotic version (`exact = FALSE` in R).
- **Reporting ρ without a CI.** Give the interval, or bootstrap it when n is small.
- **Reading ρ as a slope or as agreement.** It describes ordering, not the size of change, and not whether two methods give the same values.
- **Testing many pairs and reporting the best one.** Pre-specify the comparisons or apply a [multiple-comparison correction](#/test/multiple-comparisons).

## Related tests

- [Pearson correlation](#/test/pearson): the linear, parametric counterpart
- [Linear regression](#/test/linear-regression): slope and covariate adjustment (possibly after a log transform)
- [Kruskal–Wallis test](#/test/kruskal-wallis) and [Mann–Whitney U test](#/test/mann-whitney): other rank-based methods
- [Bland–Altman analysis](#/test/bland-altman): agreement between two methods
- [Shapiro–Wilk test](#/test/shapiro-wilk): normality check

## Key takeaways

- Spearman's ρ is Pearson's correlation applied to the ranks. It measures **monotonic**, not necessarily linear, association.
- It is robust to outliers and skew, and unchanged by any increasing transformation.
- ρ ≈ 0 does not mean no relationship: a U-shape can give ρ ≈ 0.
- Report ρ with a 95% CI and n. ρ is not a slope or an agreement measure.
- With many ties, use the tie-corrected p-value or consider Kendall's τ-b.

## References

- Altman DG. *Practical Statistics for Medical Research.* Chapman & Hall; 1991. Ch. 11.
- Bonett DG, Wright TA. Sample size requirements for estimating Pearson, Kendall and Spearman correlations. *Psychometrika.* 2000;65:23–28.
- R documentation: `?cor.test`; SciPy documentation: `scipy.stats.spearmanr`
- SAS/STAT User's Guide: The CORR Procedure (SPEARMAN and FISHER options)
