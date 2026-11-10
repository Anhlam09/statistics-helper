# Repeated-measures ANOVA

> Tests whether the mean of a continuous outcome changes across **three or more time points (or conditions) measured on the same subjects**, by separating stable between-patient differences from the within-patient effect of time.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Same subjects measured at 3+ time points or conditions | Continuous outcome, complete data | Independent subjects; roughly normal; sphericity | SBP at baseline and weeks 4, 8, 12 on a new drug |

## Learning objectives

- Explain how RM-ANOVA removes stable between-patient differences from the error term
- Recognize when sphericity matters and read Mauchly's test together with ε
- Run RM-ANOVA in R and Python and read F, df, and the Greenhouse–Geisser correction
- Follow a significant omnibus F with adjusted pairwise comparisons and report it properly

## What it is

Repeated-measures (RM) ANOVA is the extension of the [paired t-test](#/test/paired-t) to three or more measurements per subject. Every patient is measured at each visit, so the patient's overall level (some people simply run higher blood pressure than others) can be removed from the error term. What remains is the **time × patient** residual, against which the time effect is tested. This makes the test much more sensitive than a one-way ANOVA on the same data.

- **Test statistic:** F = MS_time / MS_error, where MS_error is the time × subject interaction mean square
- **Degrees of freedom:** (k − 1) and (k − 1)(n − 1) for k time points and n subjects; here 3 and 27
- **Sphericity:** the variances of all pairwise differences between time points must be equal. When this fails, the F test is too liberal, so the degrees of freedom are multiplied by an estimate ε (Greenhouse–Geisser or Huynh–Feldt)
- **Effect size:** the visit-to-visit mean differences with 95% CI (primary); partial η² = SS_time / (SS_time + SS_error) or generalized η² (secondary)

## The intuition

RM-ANOVA is the paired t-test's big sibling, and it uses the same trick: **let every patient be their own control**. The signal is how much the visit means differ from each other (MS_time). The noise is how much a patient's change fluctuates *after* removing the patient's overall level (MS_error). F is the ratio of the two.

Think of a school class tested four times in a year. Some pupils always score high, some always low. If you ignored that, the "pupil differences" would swamp any real progress. RM-ANOVA subtracts each pupil's personal level first, so only progress versus random wobble is left.

The dot plot shows our four visits. Each dot is a patient. The groups overlap at neighbouring visits because patients differ (baseline runs from 139 to 160), but the means (148.8, 141.6, 137.3, 137.4) step down and then plateau.

<div class="viz" data-viz="dots" data-unit="mmHg" data-groups='{"Baseline":[148,150,139,143,160,145,156,155,147,145],"Week 4":[138,144,134,137,151,138,151,145,143,135],"Week 8":[130,137,135,125,146,135,143,150,139,133],"Week 12":[135,142,126,127,147,137,144,141,141,134]}' data-caption="Notice that the spread within each visit (patients differ) is large compared with the drop in the means. Repeated-measures ANOVA works because it compares each patient with themselves, not with other patients. Also notice the plateau between week 8 and week 12."></div>

F = 35.93 with 3 and 27 degrees of freedom. If SBP did not change at all, F would follow the curve below, and values much above 3 would already be rare. Ours is far to the right, so the upper-tail p-value is minuscule (this is the uncorrected curve; the Greenhouse–Geisser correction shrinks the df, but the conclusion is the same).

<div class="viz" data-viz="f" data-stat="35.93" data-df1="3" data-df2="27"></div>

## When to use it

Use RM-ANOVA when **all** of these are true:

- The outcome is **continuous** (blood pressure, HbA1c, VAS pain score)
- The **same subjects** are measured at **three or more** occasions or conditions
- Every subject has **all** measurements (complete data), and the design is simple (one within-subject factor, perhaps one between-subject factor)

**Clinical scenario.** Ten patients with stage 1 hypertension start a new antihypertensive. Systolic blood pressure (SBP, mmHg) is measured at baseline and at weeks 4, 8, and 12. *Does mean SBP change over the course of treatment, and when does the effect level off?*

**Don't use it when:**

- There are only **two** time points → use the [paired t-test](#/test/paired-t)
- Some patients **miss visits** or visits are at irregular times → use a [mixed model for repeated measures](#/test/mixed-models) (RM-ANOVA drops every patient with any missing visit) or [GEE](#/test/gee)
- The outcome is skewed, ordinal, or n is very small → use the [Friedman test](#/test/friedman)
- The groups are **independent** people → use [one-way ANOVA](#/test/one-way-anova)
- You need to compare randomized arms while adjusting for baseline → use [ANCOVA](#/test/ancova) or a mixed model

## Assumptions

1. **Subjects are independent of each other.** Repeated measures within a subject are correlated by design; different subjects are not. Check from the study design.
2. **The outcome is approximately normal at each time point** (more precisely, the residuals are). Check with Q–Q plots or [Shapiro–Wilk](#/test/shapiro-wilk) on the residuals or on each visit; RM-ANOVA is fairly robust when n is moderate.
3. **Sphericity.** The variances of the differences between every pair of time points are equal (compound symmetry is a sufficient special case). Check with **Mauchly's test**, but remember that it has little power when n is small and rejects trivially when n is large. In practice many analysts report the Greenhouse–Geisser-corrected result whenever ε is clearly below 1, regardless of the Mauchly p-value. With two time points, sphericity is automatic.
4. **No missing visits.** Complete cases only; otherwise see the mixed-model alternative above.

## Hypotheses

- **H₀:** μ_baseline = μ_wk4 = μ_wk8 = μ_wk12 (the mean SBP is the same at every time point)
- **H₁:** at least one time point has a different mean

A significant omnibus F says only that the means differ somewhere. Locate the differences with pre-specified pairwise comparisons, adjusted for multiplicity (see [multiple-comparison corrections](#/test/multiple-comparisons)).

## Example data

SBP (mmHg) for 10 patients:

| Patient | Baseline | Week 4 | Week 8 | Week 12 |
|---|---|---|---|---|
| 1 | 148 | 138 | 130 | 135 |
| 2 | 150 | 144 | 137 | 142 |
| 3 | 139 | 134 | 135 | 126 |
| 4 | 143 | 137 | 125 | 127 |
| 5 | 160 | 151 | 146 | 147 |
| 6 | 145 | 138 | 135 | 137 |
| 7 | 156 | 151 | 143 | 144 |
| 8 | 155 | 145 | 150 | 141 |
| 9 | 147 | 143 | 139 | 141 |
| 10 | 145 | 135 | 133 | 134 |
| **Mean** | **148.8** | **141.6** | **137.3** | **137.4** |

## How to run it

### R

```r
library(car)
# Systolic BP (mmHg): 10 patients measured at baseline, week 4, 8 and 12
sbp <- data.frame(
  id = 1:10,
  wk0  = c(148,150,139,143,160,145,156,155,147,145),
  wk4  = c(138,144,134,137,151,138,151,145,143,135),
  wk8  = c(130,137,135,125,146,135,143,150,139,133),
  wk12 = c(135,142,126,127,147,137,144,141,141,134))

# Multivariate linear model with one column per time point
fit   <- lm(cbind(wk0, wk4, wk8, wk12) ~ 1, data = sbp)
idata <- data.frame(time = factor(c("wk0","wk4","wk8","wk12"),
                                  levels = c("wk0","wk4","wk8","wk12")))
res <- Anova(fit, idata = idata, idesign = ~time)

# Univariate F test + Mauchly + GG/HF corrections
summary(res, multivariate = FALSE)

# Mean SBP at each visit
colMeans(sbp[, -1])
```

**What the code does**

- The data are in **wide format** (one column per visit). `cbind(...) ~ 1` fits a multivariate model with only an intercept; the within-subject structure is declared separately.
- `idata` describes the within-subject factor (`time`) and `idesign = ~time` tells `car::Anova` to test it. This route is used because it reports Mauchly's test and the corrections automatically.
- `summary(res, multivariate = FALSE)` prints the **univariate** (classical) RM-ANOVA table assuming sphericity, then Mauchly's test, then the Greenhouse–Geisser and Huynh–Feldt ε and corrected p-values. The note about "only an intercept" is expected and harmless.
- The `(Intercept)` row tests whether the grand mean is zero; ignore it.
- The corrected p-values come from multiplying **both** numerator and denominator df by ε; the F statistic itself does not change.

**Output**

```text
Univariate Type III Repeated-Measures ANOVA Assuming Sphericity

            Sum Sq num Df Error SS den Df  F value    Pr(>F)    
(Intercept) 798345      1  1449.22      9 4957.895 1.188e-13 ***
time           875      3   219.27     27   35.933 1.436e-09 ***

Mauchly Tests for Sphericity

     Test statistic p-value
time        0.46038 0.31061

Greenhouse-Geisser and Huynh-Feldt Corrections
 for Departure from Sphericity

      GG eps Pr(>F[GG])    
time 0.66988  4.901e-07 ***

        HF eps   Pr(>F[HF])
time 0.8629276 1.607569e-08

  wk0   wk4   wk8  wk12 
148.8 141.6 137.3 137.4 
```

### Python

```python
import pandas as pd, pingouin as pg

# Systolic BP (mmHg): 10 patients measured at baseline, week 4, 8 and 12
wide = pd.DataFrame({
    "id":   range(1, 11),
    "wk0":  [148,150,139,143,160,145,156,155,147,145],
    "wk4":  [138,144,134,137,151,138,151,145,143,135],
    "wk8":  [130,137,135,125,146,135,143,150,139,133],
    "wk12": [135,142,126,127,147,137,144,141,141,134]})

# pingouin needs long format: one row per patient x visit
long = wide.melt(id_vars="id", var_name="time", value_name="sbp")

# Mauchly's test of sphericity
sph = pg.sphericity(long, dv="sbp", within="time", subject="id", method="mauchly")
print(f"Mauchly W = {sph.W:.4f}, chi2 = {sph.chi2:.3f}, df = {sph.dof}, p = {sph.pval:.4f}")

# RM-ANOVA with Greenhouse-Geisser correction (also reports uncorrected p)
aov = pg.rm_anova(long, dv="sbp", within="time", subject="id",
                  correction=True, detailed=True)
print(aov[["Source","SS","DF","MS","F","p_unc","p_GG_corr","eps"]]
      .to_string(float_format=lambda x: f"{x:.4g}"))
print(wide.drop(columns="id").mean())
```

**What the code does**

- `statsmodels.AnovaRM` fits the plain RM-ANOVA but offers **no** sphericity test or correction, so this page uses `pingouin`, which does.
- `melt` converts wide to **long** format (one row per patient per visit), the layout pingouin expects. Python does not know the visit order, so check that the means are sensible rather than relying on row order.
- `pg.sphericity(..., method="mauchly")` runs Mauchly's test and returns W, χ², df, and the p-value (Mauchly is the default method; it is written out here so the match with R and SAS is explicit).
- `pg.rm_anova(..., correction=True)` always computes the Greenhouse–Geisser ε and corrected p-value (`p_GG_corr`) next to the uncorrected `p_unc`. The default, `correction="auto"`, would compute them only when the sphericity test rejects, which is not what you want with a low-power test.
- `detailed=True` shows the SS and MS of the time and error rows.

**Output**

```text
Mauchly W = 0.4604, chi2 = 5.990, df = 5, p = 0.3106
  Source    SS  DF    MS     F     p_unc  p_GG_corr    eps
0   time 875.5   3 291.8 35.93 1.436e-09  4.901e-07 0.6699
1  Error 219.3  27 8.121   NaN       NaN        NaN    NaN
wk0     148.8
wk4     141.6
wk8     137.3
wk12    137.4
dtype: float64
```

### SAS

```sas
/* One row per patient, one column per visit (wide format) */
data sbp;
  input id wk0 wk4 wk8 wk12;
  datalines;
1 148 138 130 135
2 150 144 137 142
3 139 134 135 126
4 143 137 125 127
5 160 151 146 147
6 145 138 135 137
7 156 151 143 144
8 155 145 150 141
9 147 143 139 141
10 145 135 133 134
;
run;

/* Repeated-measures ANOVA with sphericity test and epsilon corrections */
proc glm data=sbp;
  model wk0 wk4 wk8 wk12 = / nouni;
  repeated time 4 (0 4 8 12) / printe summary;
run;
quit;
```

**What the code does**

- The `MODEL` statement lists all four repeated measures on the left and **no** between-subject effects on the right; `nouni` suppresses the four separate univariate analyses.
- `repeated time 4 (0 4 8 12)` names the within-subject factor `time`, with 4 levels at the spacing 0, 4, 8, 12 weeks (this spacing is used for the polynomial contrasts).
- `printe` prints the error covariance matrix and **Mauchly's test of sphericity** (SAS reports it on the orthogonalized components).
- The output section "Univariate Tests of Hypotheses for Within Subject Effects" gives the unadjusted F test and, in the columns **Adj Pr > F: G-G** and **H-F-L**, the Greenhouse–Geisser and Huynh–Feldt corrected p-values, with both epsilons listed underneath.
- `summary` adds per-contrast ANOVA tables (linear, quadratic, cubic trend over time).
- For unbalanced or incomplete data, use `PROC MIXED` with a `repeated` statement and an appropriate covariance structure instead.

## Interpreting the output

- **Means:** 148.8, 141.6, 137.3, 137.4 mmHg at baseline, weeks 4, 8, and 12. SBP falls by about 11 mmHg in total, mostly by week 8, and then plateaus.
- **Sphericity:** Mauchly W = 0.460, p = 0.31. There is no formal evidence against sphericity, but with n = 10 this test has low power, and ε = 0.67 (well below 1) suggests a real departure. The cautious choice is to report the Greenhouse–Geisser result.
- **Omnibus test (uncorrected):** F(3, 27) = 35.93, p = 1.4 × 10⁻⁹. The sums of squares are SS_time = 875.5 and SS_error = 219.3, so MS_time / MS_error = 291.8 / 8.12.
- **Greenhouse–Geisser corrected:** ε = 0.670, which gives df of 3 × 0.670 ≈ 2.0 and 27 × 0.670 ≈ 18.1, and p = 4.9 × 10⁻⁷. The conclusion is unchanged, because the effect is large.
- **Huynh–Feldt:** ε = 0.863, p = 1.6 × 10⁻⁸ (less conservative than GG).
- **Effect size:** partial η² = 875.5 / (875.5 + 219.3) = 0.80. As always, the mean changes in mmHg matter more clinically; follow up with adjusted pairwise comparisons.

**How to report it:**
> Mean systolic blood pressure fell from 148.8 mmHg at baseline to 141.6, 137.3, and 137.4 mmHg at weeks 4, 8, and 12 (repeated-measures ANOVA with Greenhouse–Geisser correction, ε = 0.67: F(2.0, 18.1) = 35.93, p < 0.001).

⚠️ **Important for trials:** RM-ANOVA applied to a single arm shows that the outcome changed, not that the treatment caused it. In randomized trials the preferred analysis of repeated measures is a [mixed model (MMRM)](#/test/mixed-models), which uses all available data under a missing-at-random assumption.

## Check your understanding

<details class="check"><summary>1. A trial measures HbA1c at baseline, 3, 6 and 12 months. Three of the 40 patients miss the month-6 visit. Is RM-ANOVA still the right tool?</summary>

It will run, but it silently drops all three patients, which wastes data and can bias results if the missingness is related to outcome. A [mixed model (MMRM)](#/test/mixed-models) or [GEE](#/test/gee) uses every available measurement and is the preferred approach for this design.
</details>

<details class="check"><summary>2. Mauchly's test gives p = 0.31, but ε = 0.67. A colleague says "p > 0.05, so sphericity holds and we can ignore the correction." Comment.</summary>

A non-significant Mauchly test is not proof that sphericity holds; with n = 10 it has very little power. An ε of 0.67 (the minimum possible here is 1/3, the maximum is 1) suggests a real departure. The safe habit is to report the Greenhouse–Geisser result whenever ε is clearly below 1: F(2.0, 18.1) = 35.93, p = 4.9 × 10⁻⁷.
</details>

<details class="check"><summary>3. The omnibus F is significant (p < 0.001). Can you conclude that SBP at week 12 is lower than at week 8?</summary>

No. The omnibus test says only that the four means are not all equal. Here the means at weeks 8 and 12 are practically identical (137.3 vs 137.4), so most of the effect is from baseline versus later. Use pre-specified pairwise comparisons with a Holm or Bonferroni adjustment to locate the differences.
</details>

<details class="check"><summary>4. What would happen to F if the patient-to-patient differences in how much they respond to treatment were much larger (the same visit means, but noisier time × patient residuals)?</summary>

MS_error would grow while MS_time stays the same, so F would shrink and p would rise. The means pattern is the signal; the inconsistency of individual trajectories is the noise. The baseline differences between patients, however, do not add noise, because RM-ANOVA removes them.
</details>

## Common mistakes

- **Ignoring sphericity.** Reading only the uncorrected F and p can give false positives when ε is low. Report the Greenhouse–Geisser (or Huynh–Feldt when ε > 0.75) result.
- **Relying on Mauchly's p alone.** It has low power in small samples. Look at ε as well.
- **Running a one-way ANOVA on repeated data.** This treats the visits as independent groups, inflates the error term, and violates independence.
- **Deleting patients with missing visits silently.** RM-ANOVA is a complete-case analysis and can bias results if dropout is related to the outcome. Use a [mixed model](#/test/mixed-models) and report how many patients were analysed.
- **Unadjusted pairwise tests after a significant F.** Use a Bonferroni or Holm adjustment for the follow-up comparisons.
- **Wrong data layout.** R's `car::Anova` and SAS PROC GLM want wide data; pingouin and `aov()` with an `Error()` term want long data.

## Related tests

- [Paired t-test](#/test/paired-t): the two-time-point special case
- [Friedman test](#/test/friedman): nonparametric alternative
- [Linear mixed-effects models (MMRM)](#/test/mixed-models): handles missing visits and flexible covariance structures
- [Generalized estimating equations (GEE)](#/test/gee): population-averaged repeated-measures models, including non-normal outcomes
- [Multiple-comparison corrections](#/test/multiple-comparisons): for the follow-up pairwise comparisons
- [One-way ANOVA](#/test/one-way-anova): the independent-groups analogue

## Key takeaways

- RM-ANOVA tests whether means differ across 3+ measurements on the same subjects: F = MS_time / MS_error, df = (k − 1) and (k − 1)(n − 1).
- Between-patient differences are removed from the error term, which makes it more powerful than one-way ANOVA.
- Sphericity (equal variances of all pairwise differences) is the key extra assumption. Report Greenhouse–Geisser when ε is clearly below 1, even if Mauchly's p is large.
- A significant F says only "the means differ somewhere". Follow up with adjusted pairwise comparisons.
- Patients with any missing visit are dropped, so with missing data use a mixed model (MMRM).

## References

- Maxwell SE, Delaney HD, Kelley K. *Designing Experiments and Analyzing Data: A Model Comparison Perspective.* 3rd ed. Routledge; 2017. Chapters on within-subjects designs.
- ICH E9 *Statistical Principles for Clinical Trials* (1998), on handling repeated measurements and missing data.
- R documentation: `?car::Anova` and `?summary.Anova.mlm`
- pingouin documentation: `pingouin.rm_anova` and `pingouin.sphericity`
- SAS/STAT User's Guide: The GLM Procedure, REPEATED statement
