# MANOVA

> Tests whether two or more groups differ on a **set of correlated continuous outcomes taken together** (for example SBP, DBP and heart rate), using one omnibus test instead of one ANOVA per outcome.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Compare 2+ groups on several correlated continuous outcomes at once | Continuous outcomes (2–3 typical), categorical group | Independent subjects; multivariate normal residuals; equal covariance matrices | Drug A, Drug B, placebo compared on SBP, DBP and heart rate |

## Learning objectives

- Explain what MANOVA tests that a set of separate ANOVAs does not
- Run a one-way MANOVA in R and Python and read Pillai's trace, Wilks' lambda, F and p
- Follow a significant MANOVA with univariate ANOVAs, and know what that follow-up does and does not show
- Explain why MANOVA is rare in clinical trials and name the better-suited alternatives

## What it is

A one-way ANOVA asks whether group means differ on **one** outcome. MANOVA asks whether the groups differ on the **vector of means** of several outcomes at once, i.e. whether the point (SBP, DBP, heart rate) sits in a different place for each group. It uses the correlation among the outcomes, so a pattern that is small on every single outcome but consistent across them can still be detected.

Four classic test statistics exist. All are built from the between-group (H) and within-group (E) sums-of-squares-and-cross-products matrices, and all give an (approximate) F test:

- **Pillai's trace:** V = Σ λᵢ/(1 + λᵢ), where λᵢ are the eigenvalues of H·E⁻¹. It ranges from 0 to min(p, g − 1) and is the most **robust** to unequal covariance matrices and non-normality. It is R's default.
- **Wilks' lambda:** Λ = |E| / |H + E|, the share of generalized variance **not** explained by group. Smaller means stronger separation. It is the traditional choice and is closely tied to likelihood-ratio ideas.
- **Hotelling–Lawley trace** and **Roy's largest root:** more powerful when group differences lie along a single dimension, but less robust (Roy's is the most fragile).
- **Effect size:** partial η² from Pillai (V / s, where s = min(p, g − 1)), plus the **raw group mean differences in clinical units** for each outcome with 95% CIs.

With two groups, MANOVA reduces to Hotelling's T² test. With one outcome, it reduces to [one-way ANOVA](#/test/one-way-anova).

## The intuition

Think of each patient as a point on a map with SBP on one axis and DBP on another (plus heart rate as a third). Separate ANOVAs look at the map one axis at a time. MANOVA asks whether the **clouds of points** for the three groups sit in different places. Because SBP and DBP move together (r = 0.78 here), a patient whose SBP is a bit high but DBP is unusually low is odd in a way neither single outcome reveals. MANOVA accounts for that shared movement.

The dot plot below shows one slice of the picture, SBP by arm. Both active arms sit about 14 mmHg below placebo, with some overlap. DBP shows the same pattern at a smaller scale; heart rate shows none. MANOVA pools all three.

<div class="viz" data-viz="dots" data-unit="mmHg (SBP)" data-groups='{"Placebo":[136,170,146,162,128,149,140,145,141,157,162,152,149,138,161],"Drug A":[145,132,143,142,136,121,137,130,144,135,134,137,137,140,115],"Drug B":[147,148,134,143,143,122,139,124,136,130,136,124,124,141,133]}' data-caption="SBP by arm (n = 15 each). The active arms are shifted down by about 14 mmHg; MANOVA tests the shift on SBP, DBP and heart rate together."></div>

Pillai's trace converts to an F statistic. Ours is F(6, 82) = 2.97. Under H₀ (identical mean vectors in all groups), F values follow the curve below, and ours falls in the shaded upper tail.

<div class="viz" data-viz="f" data-stat="2.97" data-df1="6" data-df2="82"></div>

## When to use it

Use MANOVA when **all** of these are true:

- There are **2–3 (occasionally a few more) continuous outcomes** that are conceptually related and **correlated**
- You want **one omnibus claim** ("these groups differ on the cardiovascular profile") rather than separate claims per outcome
- Subjects are independent, and you have **more patients per group than outcomes**, comfortably (at least about 20 per group is a common guide)

**Clinical scenario.** A phase 2 study randomizes 45 adults with stage 1 hypertension to Drug A, Drug B or placebo (15 per arm). At week 8 the team records SBP (mmHg), DBP (mmHg) and resting heart rate (bpm). *Do the three arms differ on the cardiovascular profile, and if so, on which measurements?*

**Don't use it when:**

- There is **one** outcome → use [one-way ANOVA](#/test/one-way-anova) (or [Welch's ANOVA](#/test/welch-anova) for unequal variances)
- The outcome is measured **repeatedly over time** → use a [mixed model](#/test/mixed-models), [GEE](#/test/gee) or [repeated-measures ANOVA](#/test/repeated-measures-anova)
- You need to **adjust for baseline** (almost always true in trials) → use [ANCOVA](#/test/ancova) on the pre-specified primary endpoint
- One outcome is clearly **primary** and the others are secondary → test the primary and control multiplicity with a [hierarchical procedure](#/test/multiple-comparisons)
- Groups are tiny relative to the number of outcomes, or outcomes are heavily skewed → consider permutation methods, or [Kruskal–Wallis](#/test/kruskal-wallis) per outcome

## Assumptions

1. **Independent observations.** One row per patient; no clustering or repeated measures.
2. **Multivariate normality of the residuals** (within groups). Check each outcome's residuals with a Q–Q plot or [Shapiro–Wilk](#/test/shapiro-wilk), and look at pairwise scatterplots for outliers. With about 20+ per group and balanced arms, Pillai's trace is fairly tolerant.
3. **Homogeneity of covariance matrices.** The spread **and** the correlations among outcomes should be similar across groups. Box's M test checks this but is very sensitive; prefer Pillai's trace and equal group sizes instead of relying on it.
4. **Linearity among outcomes** and **no severe multicollinearity.** If two outcomes correlate above about 0.9, they are essentially duplicates; drop or combine one.
5. **Adequate sample size.** Each group needs more subjects than there are outcomes; power drops quickly when outcomes are weakly related.

## Hypotheses

- **H₀:** μ(Placebo) = μ(Drug A) = μ(Drug B), where each μ is the vector (SBP, DBP, heart rate) of means. The groups have identical mean vectors.
- **H₁:** at least one group's mean vector differs from the others, on at least one outcome.

## Example data

Simulated data, 15 patients per arm, three correlated outcomes (SBP and DBP correlate at 0.78; heart rate is unrelated to them). Group means (SD):

| Arm | SBP (mmHg) | DBP (mmHg) | Heart rate (bpm) |
|---|---|---|---|
| Placebo | 149.1 (11.6) | 88.7 (8.3) | 73.8 (5.9) |
| Drug A | 135.2 (8.3) | 80.9 (7.0) | 75.4 (5.9) |
| Drug B | 134.9 (8.7) | 81.1 (7.6) | 73.5 (4.7) |

The full 45 rows are typed inline in each code block below.

## How to run it

### R

```r
# SBP, DBP (mmHg) and heart rate (bpm); 15 patients per arm, order = Drug A, Drug B, Placebo
v <- function(s) scan(text = s, quiet = TRUE)
dat <- data.frame(
  grp = factor(rep(c("Drug A", "Drug B", "Placebo"), each = 15),
               levels = c("Placebo", "Drug A", "Drug B")),
  sbp = v("145 132 143 142 136 121 137 130 144 135 134 137 137 140 115 147
      148 134 143 143 122 139 124 136 130 136 124 124 141 133 136 170
      146 162 128 149 140 145 141 157 162 152 149 138 161"),
  dbp = v("89 77 88 85 80 65 79 84 83 85 84 76 88 82 68 97 81 78 90 86 79 79
      79 93 78 78 74 77 81 67 71 98 86 98 86 91 89 90 82 92 101 80 87 81
      99"),
  hr  = v("83 75 76 69 76 65 71 65 82 74 84 81 75 77 78 78 76 72 76 66 72 64
      80 72 74 77 72 81 70 73 73 69 83 71 67 67 68 73 76 72 78 71 81 86
      72"))

# 1. Correlation among the outcomes (the reason to consider MANOVA)
round(cor(dat[, c("sbp", "dbp", "hr")]), 2)

# 2. One-way MANOVA: Pillai's trace (default) and Wilks' lambda
fit <- manova(cbind(sbp, dbp, hr) ~ grp, data = dat)
summary(fit, test = "Pillai")
summary(fit, test = "Wilks")

# 3. Follow-up univariate ANOVAs (only after a significant omnibus test)
summary.aov(fit)
```

**What the code does**

- `scan(text = ..., quiet = TRUE)` reads the compact data; the arms are ordered Drug A, Drug B, Placebo in each vector. Placebo is set as the reference level.
- `cor()` shows how related the outcomes are. If they are all near 0, MANOVA gains little over separate ANOVAs.
- `manova(cbind(sbp, dbp, hr) ~ grp)` fits the multivariate linear model. `cbind()` on the left stacks the outcomes.
- `summary(fit, test = "Pillai")` is the default; `test = "Wilks"` gives Wilks' lambda. Hotelling–Lawley and Roy are also available.
- `summary.aov(fit)` prints one ordinary ANOVA per outcome: the usual **follow-up** after a significant omnibus test.

**Output**

```text
      sbp   dbp    hr
sbp  1.00  0.78 -0.03
dbp  0.78  1.00 -0.01
hr  -0.03 -0.01  1.00
          Df  Pillai approx F num Df den Df  Pr(>F)  
grp        2 0.35666   2.9661      6     82 0.01135 *
Residuals 42                                         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
          Df   Wilks approx F num Df den Df   Pr(>F)   
grp        2 0.65045   3.1989      6     80 0.007248 **
Residuals 42                                           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
 Response sbp :
            Df Sum Sq Mean Sq F value    Pr(>F)    
grp          2 1960.5  980.27  10.513 0.0001988 ***
Residuals   42 3916.3   93.24                      
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

 Response dbp :
            Df  Sum Sq Mean Sq F value  Pr(>F)  
grp          2  598.58  299.29  5.1256 0.01019 *
Residuals   42 2452.40   58.39                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

 Response hr :
            Df  Sum Sq Mean Sq F value Pr(>F)
grp          2   30.58  15.289  0.4987 0.6109
Residuals   42 1287.73  30.660               
```

### Python

```python
import numpy as np, pandas as pd
import statsmodels.api as sm
from statsmodels.multivariate.manova import MANOVA
from statsmodels.formula.api import ols

# SBP, DBP (mmHg) and heart rate (bpm); 15 patients per arm, order = Drug A, Drug B, Placebo
v = lambda s: np.array(s.split(), dtype=float)
dat = pd.DataFrame({
    "grp": np.repeat(["Drug A", "Drug B", "Placebo"], 15),
    "sbp": v("145 132 143 142 136 121 137 130 144 135 134 137 137 140 "
            "115 147 148 134 143 143 122 139 124 136 130 136 124 124 "
            "141 133 136 170 146 162 128 149 140 145 141 157 162 152 "
            "149 138 161"),
    "dbp": v("89 77 88 85 80 65 79 84 83 85 84 76 88 82 68 97 81 78 90 "
            "86 79 79 79 93 78 78 74 77 81 67 71 98 86 98 86 91 89 90 "
            "82 92 101 80 87 81 99"),
    "hr":  v("83 75 76 69 76 65 71 65 82 74 84 81 75 77 78 78 76 72 76 "
            "66 72 64 80 72 74 77 72 81 70 73 73 69 83 71 67 67 68 73 "
            "76 72 78 71 81 86 72")})

# 1. Correlation among the outcomes
print(dat[["sbp", "dbp", "hr"]].corr().round(2))

# 2. One-way MANOVA: mv_test() prints Wilks, Pillai, Hotelling-Lawley, Roy
fit = MANOVA.from_formula("sbp + dbp + hr ~ C(grp, Treatment('Placebo'))", data=dat)
print(fit.mv_test())

# 3. Follow-up univariate ANOVAs
for y in ["sbp", "dbp", "hr"]:
    a = sm.stats.anova_lm(ols(f"{y} ~ C(grp)", data=dat).fit(), typ=1)
    print(y, f"F = {a.loc['C(grp)', 'F']:.2f}, p = {a.loc['C(grp)', 'PR(>F)']:.4f}")
```

**What the code does**

- `MANOVA.from_formula("sbp + dbp + hr ~ ...")` puts the outcomes on the left, joined by `+`. `Treatment('Placebo')` sets the reference arm.
- `mv_test()` prints all four statistics at once for the intercept and the group term. Only the **group** block matters; the intercept test (F ≈ 1967) just says the mean vector differs from zero, which is meaningless here.
- The loop runs the follow-up univariate ANOVAs with `ols` and `anova_lm`.

**Output**

```text
      sbp   dbp    hr
sbp  1.00  0.78 -0.03
dbp  0.78  1.00 -0.01
hr  -0.03 -0.01  1.00
                    Multivariate linear model
==================================================================
                                                                  
------------------------------------------------------------------
         Intercept         Value   Num DF  Den DF  F Value  Pr > F
------------------------------------------------------------------
            Wilks' lambda   0.0067 3.0000 40.0000 1967.3987 0.0000
           Pillai's trace   0.9933 3.0000 40.0000 1967.3987 0.0000
   Hotelling-Lawley trace 147.5549 3.0000 40.0000 1967.3987 0.0000
      Roy's greatest root 147.5549 3.0000 40.0000 1967.3987 0.0000
------------------------------------------------------------------
                                                                  
------------------------------------------------------------------
 C(grp, Treatment('Placebo')) Value  Num DF  Den DF F Value Pr > F
------------------------------------------------------------------
                Wilks' lambda 0.6504 6.0000 80.0000  3.1989 0.0072
               Pillai's trace 0.3567 6.0000 82.0000  2.9662 0.0114
       Hotelling-Lawley trace 0.5265 6.0000 51.5932  3.4687 0.0059
          Roy's greatest root 0.5048 3.0000 41.0000  6.8990 0.0007
==================================================================

sbp F = 10.51, p = 0.0002
dbp F = 5.13, p = 0.0102
hr F = 0.50, p = 0.6109
```

## Interpreting the output

- **Correlations:** SBP and DBP correlate at 0.78, heart rate is unrelated (≈ 0). So SBP and DBP carry overlapping information.
- **Pillai's trace = 0.357, F(6, 82) = 2.97, p = 0.011.** The mean vectors differ across arms. The maximum possible Pillai is 2 here (s = min(3 outcomes, 3 − 1 groups)), so 0.357/2 ≈ 0.18 is a moderate multivariate effect (partial η²).
- **Wilks' lambda = 0.650, F(6, 80) = 3.20, p = 0.007.** About 65% of the generalized variance is *not* explained by arm. R and Python agree to the printed digits. The two statistics give slightly different F and df because they approximate differently; with balanced groups and a clear effect, the conclusions agree. If they ever disagree, prefer Pillai.
- **Follow-up ANOVAs:** SBP F(2, 42) = 10.51, p = 0.0002; DBP F(2, 42) = 5.13, p = 0.010; heart rate F(2, 42) = 0.50, p = 0.61. The multivariate signal comes from the two pressure measures. Heart rate contributes nothing.
- **Multiplicity:** the three univariate tests are protected only loosely by the omnibus test. With a Bonferroni correction over three outcomes, SBP stays significant (p = 0.0006) and DBP (p = 0.031) just survives.
- **What MANOVA does not tell you:** which arms differ from which. For that, run post-hoc comparisons on the outcomes of interest ([Tukey](#/test/tukey-hsd) or [Dunnett vs placebo](#/test/dunnett)), reporting mean differences with 95% CIs.

**How to report it:**
> A one-way MANOVA showed that the arms differed on the combined outcomes SBP, DBP and heart rate (Pillai's trace = 0.36, F(6, 82) = 2.97, p = 0.011). Follow-up ANOVAs with Bonferroni correction indicated differences in SBP (F(2, 42) = 10.51, adjusted p < 0.001) and DBP (F(2, 42) = 5.13, adjusted p = 0.031) but not heart rate (F(2, 42) = 0.50, p > 0.99). Mean SBP was 135.2 mmHg (Drug A), 134.9 (Drug B) and 149.1 (placebo).

## Why MANOVA is rare in clinical trials

Despite its elegance, MANOVA is seldom the primary analysis of a clinical trial:

- **ICH E9 wants one pre-specified primary endpoint** (or a clearly structured set). "The vector of outcomes differs" is not a treatment effect a clinician or regulator can act on or put on a label.
- **The result is not interpretable in clinical units.** A significant Pillai's trace does not say how many mmHg the drug lowers SBP.
- **The test can be significant for the wrong reason.** A combination of tiny differences, or a difference in correlation structure, can drive it, and a direction like "SBP down, heart rate up" is hidden.
- **Follow-up ANOVAs recreate the multiplicity problem** that MANOVA seemed to solve.
- **Trials have baselines, covariates and repeated visits.** MANOVA does not handle these elegantly; ANCOVA and mixed models do.

**What to use instead:**

- A **single pre-specified primary endpoint** analysed with [ANCOVA](#/test/ancova), with the others as secondary under a hierarchical (gatekeeping) or [multiplicity-adjusted](#/test/multiple-comparisons) plan
- A [**mixed model**](#/test/mixed-models) (or [GEE](#/test/gee)) for outcomes measured at several visits, which can also model SBP and DBP jointly as a multivariate response with an outcome-type factor
- A **composite or global test** (e.g., an O'Brien OLS/GLS global statistic) when you must claim benefit on "any of several correlated outcomes"
- MANOVA still fits in **exploratory** or **laboratory** work with several measures on the same units and no obvious primary one.

## Check your understanding

<details class="check"><summary>1. A trial measures SBP, DBP and heart rate on 90 patients in 3 arms. A colleague proposes running three separate ANOVAs at α = 0.05 each. What is the problem, and how does MANOVA help?</summary>

Three tests at 0.05 each give an overall false-positive chance of about 14% if there is truly no effect (1 − 0.95³). MANOVA gives **one** omnibus test at α = 0.05 for "any difference on the vector." It also uses the correlation between SBP and DBP. But a follow-up univariate ANOVA on each outcome again needs a multiplicity strategy.
</details>

<details class="check"><summary>2. MANOVA gives p = 0.011, so a colleague writes: "Drug A lowers blood pressure more than Drug B." Is that supported?</summary>

No. A significant MANOVA only says that **at least one** group's mean vector differs. It does not say which outcome drove the difference, or which pairs of groups differ. In our data, both active arms were about 14 mmHg below placebo on SBP, with nearly identical means (135.2 vs 134.9). To compare Drug A and Drug B, use planned or post-hoc contrasts on the outcome of interest, with CIs.
</details>

<details class="check"><summary>3. Suppose heart rate were highly variable and unrelated to the other outcomes, while SBP and DBP still showed the same differences. What would you expect to happen to the MANOVA p-value and to the follow-up ANOVAs?</summary>

The follow-up ANOVAs for SBP and DBP would be unchanged (each outcome is tested alone). The omnibus MANOVA would typically get **weaker**: adding an outcome that carries only noise adds degrees of freedom and noise without adding signal, so the overall F can shrink. This is why you should include only outcomes that you expect to be affected.
</details>

<details class="check"><summary>4. The follow-up ANOVA for DBP gives p = 0.0102. Apply a Bonferroni correction for the three outcomes. Is it still significant at 0.05?</summary>

Adjusted p = 0.0102 × 3 = 0.031, which is still below 0.05. SBP (0.0002 × 3 = 0.0006) also remains significant, and heart rate (0.61 × 3, capped at 1) does not.
</details>

## Common mistakes

- **Running MANOVA when there is a clear primary endpoint.** Analyse the primary endpoint by itself (ANCOVA) and treat the rest as secondary.
- **Interpreting the omnibus result as an effect on every outcome.** Significance can come from just one outcome.
- **Skipping the follow-up logic.** Univariate ANOVAs after a significant MANOVA still need multiplicity control and effect sizes in clinical units.
- **Including too many outcomes for the sample size.** Groups need many more subjects than outcomes; otherwise power and the stability of the covariance matrix collapse.
- **Relying on Box's M to "prove" equal covariances.** It is hypersensitive to non-normality. Use balanced groups and Pillai's trace.
- **Adding uncorrelated or irrelevant outcomes** to "improve" the test. It usually reduces power.

## Related tests

- [One-way ANOVA](#/test/one-way-anova): the single-outcome version
- [ANCOVA](#/test/ancova): the usual trial analysis for one outcome with baseline adjustment
- [Mixed models](#/test/mixed-models): repeated outcomes over time, multiple correlated responses
- [Multiple comparisons](#/test/multiple-comparisons): controlling the false-positive rate across outcomes
- [Tukey HSD](#/test/tukey-hsd) and [Dunnett](#/test/dunnett): post-hoc group comparisons on one outcome

## Key takeaways

- MANOVA tests whether groups differ on the **vector of means** of several correlated outcomes, with one omnibus test.
- **Pillai's trace** is the most robust statistic and the safest default; Wilks' lambda is the traditional alternative. They usually agree.
- A significant MANOVA does not say **which outcome** or **which groups** differ; follow up carefully, with multiplicity control and CIs.
- In clinical trials, prefer **one pre-specified primary endpoint** (ANCOVA) or a **mixed model**; MANOVA's result is hard to translate into a clinical effect.

## References

- Johnson RA, Wichern DW. *Applied Multivariate Statistical Analysis.* 6th ed. Pearson; 2007. Ch. 6.
- ICH E9: *Statistical Principles for Clinical Trials.* 1998 (primary endpoint and multiplicity).
- R documentation: `?manova`, `?summary.manova`
- statsmodels documentation: `statsmodels.multivariate.manova.MANOVA`
