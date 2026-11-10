# Cohen's kappa

> Measures how well two raters agree on a categorical rating after removing the agreement expected by chance alone, with a weighted version for ordinal grades.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Agreement between two raters on a categorical rating, beyond chance | Nominal or ordinal categories | Same categories; independent ratings; same two raters for all patients | Two radiologists grading chest X-ray severity |

## Learning objectives

- Explain why raw percent agreement overstates reliability
- Compute kappa from observed and chance-expected agreement
- Choose between unweighted and weighted kappa for ordinal ratings
- Report kappa with its 95% CI and interpret benchmarks cautiously

## What it is

Two raters (two radiologists, two pathologists, or a clinician and an algorithm) classify the same patients into the same categories. The simplest summary is the **observed agreement** pₒ, the proportion of patients given the same category. But some agreement happens by luck, especially when one category is very common. **Cohen's kappa** corrects for this:

κ = (pₒ − pₑ) / (1 − pₑ)

where pₑ is the agreement expected if the two raters classified independently, given each rater's own category frequencies. κ = 1 is perfect agreement, κ = 0 is agreement no better than chance, and negative values mean agreement worse than chance.

For **ordinal** categories (mild < moderate < severe), a mild-vs-severe disagreement is worse than a mild-vs-moderate one. **Weighted kappa** gives partial credit for near misses. *Linear* weights penalise a disagreement in proportion to the distance between categories; *quadratic* (Fleiss–Cohen) weights penalise by the squared distance and give a value close to the ICC.

- **Test statistic:** κ̂ with an asymptotic standard error; z = κ̂ / SE tests H₀: κ = 0 (rarely informative, because almost any real rating process beats chance)
- **Effect size:** κ itself with its 95% CI. Rough bands (Landis & Koch): < 0.20 slight, 0.21–0.40 fair, 0.41–0.60 moderate, 0.61–0.80 substantial, > 0.80 almost perfect. These are widely quoted conventions, not rules: kappa depends on how common each category is and how many categories there are, so the same level of rater skill can produce different kappas in different settings. Judge κ and its CI against what the clinical use requires.

## The intuition

Two people who each say "yes" 90% of the time will agree most of the time even if they pay no attention at all. **Kappa asks how much of the agreement beyond luck was actually achieved.** The idea is a ratio: (what you got − what chance gives) / (what perfect would give − what chance gives).

Think of two weather forecasters in a desert who both predict "no rain" almost every day. They agree 95% of the time, but that tells you little about skill, because "no rain" is almost always right by default. Kappa subtracts out that built-in agreement.

In our data, the radiologists agreed on 59 of 80 X-rays (73.8%). But from each radiologist's own category totals, we would expect about 27 agreements by chance alone, 33.8% of the patients. The table shows where those chance agreements come from (row total × column total ÷ 80):

| Grade | Observed agreements | Expected by chance |
|---|---|---|
| Mild | 22 | 29 × 27 / 80 = 9.8 |
| Moderate | 20 | 30 × 30 / 80 = 11.3 |
| Severe | 17 | 21 × 23 / 80 = 6.0 |
| **Total** | **59 (0.738)** | **27.1 (0.338)** |

Kappa = (0.738 − 0.338) / (1 − 0.338) = 0.60. The extra agreement the radiologists achieved is 60% of the maximum extra agreement possible.

<div class="viz" data-viz="dots" data-unit="grade steps (B minus A), disagreements only" data-groups='{"Disagreements":[-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,1,1,2]}' data-zero="1" data-caption="The 21 disagreements among 80 X-rays (the 59 agreements at 0 are not drawn). Almost all are one grade apart; only one is two grades apart."></div>

**What to notice:** chance agreement is large (a third of the patients), so percent agreement alone would flatter the raters. Also, most disagreements are one grade apart (mild/moderate, moderate/severe). Weighted kappa gives partial credit for those near misses, which is why it comes out higher (0.68 linear, 0.76 quadratic).

**A word on benchmarks.** Labels such as "moderate" or "substantial" are rough conventions. Kappa depends on how common each category is and on the number of categories, so the same skill can give different kappas in different populations. Judge κ and its confidence interval against what your clinical use requires.

## When to use it

Use Cohen's kappa when **all** of these are true:

- The rating is **categorical** (nominal or ordinal) with the **same categories** for both raters
- **Exactly two** raters (or two methods) rate **every** patient
- The same two raters are used for all patients, and the patients were rated independently

**Clinical scenario.** Two radiologists independently grade the severity of 80 chest X-rays as mild, moderate or severe. Before the grading is used as a trial entry criterion, the study team wants to know how well the two readers agree. *Is their agreement better than chance, and how strong is it?*

**Don't use it when:**

- There are **three or more raters**, or different raters for different patients → use Fleiss' kappa or Krippendorff's alpha (not covered here)
- The rating is **continuous** (tumour size, blood pressure) → use the [ICC](#/test/icc) or [Bland–Altman analysis](#/test/bland-altman)
- You want to know whether the raters have **different overall rates** of a positive rating (marginal homogeneity) → use [McNemar's test](#/test/mcnemar) for a binary rating
- The goal is accuracy against a gold standard rather than agreement between raters → use [ROC analysis](#/test/roc-auc) or sensitivity/specificity
- The categories are very unbalanced (for example 95% negative) and you rely on kappa alone, since it can be low despite high raw agreement; report pₒ and the category frequencies as well

## Assumptions

1. **Both raters use the same set of categories** for every patient. Check that the table is square and that the categories are in the same order.
2. **Ratings are independent between raters.** Each rater must be blinded to the other's rating.
3. **Patients are independent**, and the same two raters rate all of them.
4. **Weighted kappa only: the categories are truly ordered**, and the chosen weights (linear or quadratic) are pre-specified. Different weights give different values.
5. **The sample is large enough for the asymptotic CI.** With small n or a kappa near 1, use a bootstrap CI. Check n and the table for empty cells.

## Hypotheses

- **H₀:** κ = 0 (agreement is no better than chance)
- **H₁:** κ ≠ 0

The more useful question is usually one of estimation: is κ's lower confidence limit above the level of agreement the application requires (for example 0.6)?

## Example data

80 chest X-rays, severity graded by two radiologists (rows = Radiologist A, columns = Radiologist B). The diagonal (22 + 20 + 17 = 59) is the number of agreements.

| A \ B | Mild | Moderate | Severe | Total |
|---|---|---|---|---|
| **Mild** | 22 | 6 | 1 | 29 |
| **Moderate** | 5 | 20 | 5 | 30 |
| **Severe** | 0 | 4 | 17 | 21 |
| **Total** | 27 | 30 | 23 | 80 |

## How to run it

### R

```r
library(DescTools)

# Cross-tabulation of 80 chest X-rays: rows = Radiologist A, columns = Radiologist B
# Severity: 1 = mild, 2 = moderate, 3 = severe
tab <- matrix(c(22,  6,  1,
                 5, 20,  5,
                 0,  4, 17), nrow = 3, byrow = TRUE,
              dimnames = list(A = c("mild","moderate","severe"),
                              B = c("mild","moderate","severe")))
n <- sum(tab); po <- sum(diag(tab)) / n
cat("n =", n, " observed agreement =", round(po, 3), "\n")

# Unweighted kappa: any disagreement counts equally
CohenKappa(tab, weights = "Unweighted", conf.level = 0.95)

# Weighted kappa: linear ("Equal-Spacing") and quadratic ("Fleiss-Cohen") weights
CohenKappa(tab, weights = "Equal-Spacing", conf.level = 0.95)
CohenKappa(tab, weights = "Fleiss-Cohen",  conf.level = 0.95)

# Same result from patient-level data with the irr package
long <- as.data.frame(as.table(tab))
pat  <- long[rep(seq_len(nrow(long)), long$Freq), c("A", "B")]
irr::kappa2(pat, weight = "unweighted")$value
```

**What the code does**

- The data are a 3 × 3 matrix of counts. `CohenKappa()` from **DescTools** takes the table directly and returns κ with a confidence interval.
- `weights = "Unweighted"` treats every disagreement alike. `"Equal-Spacing"` is **linear** weighting and `"Fleiss-Cohen"` is **quadratic** weighting. Both weighted versions assume the categories are ordered as they appear in the table.
- The `irr::kappa2()` lines show the same unweighted kappa from patient-level data (one row per patient); `irr` needs that format and has no CI. Use `weight = "equal"` or `"squared"` for the weighted versions.
- `conf.level = 0.95` requests the 95% CI; without it `CohenKappa()` returns only the point estimate.

**Output**

```text
n = 80  observed agreement = 0.738 
    kappa    lwr.ci    upr.ci 
0.6032121 0.4564622 0.7499620 
    kappa    lwr.ci    upr.ci 
0.6781273 0.5514384 0.8048162 
    kappa    lwr.ci    upr.ci 
0.7580645 0.6449538 0.8711752 
[1] 0.6032121
```

### Python

```python
import numpy as np
from sklearn.metrics import cohen_kappa_score
from statsmodels.stats.inter_rater import cohens_kappa

# 80 chest X-rays: rows = Radiologist A, columns = Radiologist B (1 = mild, 2 = moderate, 3 = severe)
tab = np.array([[22,  6,  1],
                [ 5, 20,  5],
                [ 0,  4, 17]])
print(f"n = {tab.sum()}, observed agreement = {np.trace(tab)/tab.sum():.3f}")

# Expand the table to one row per patient for sklearn
a, b = np.nonzero(tab)
ra = np.repeat(a, tab[a, b]) + 1
rb = np.repeat(b, tab[a, b]) + 1

# Point estimates: unweighted, linear-weighted, quadratic-weighted
for w in (None, "linear", "quadratic"):
    print(f"sklearn kappa ({w or 'unweighted'}): {cohen_kappa_score(ra, rb, weights=w):.3f}")

# statsmodels adds the standard error and a 95% CI
for name, wt in (("unweighted", None), ("linear", "linear"), ("quadratic", "quadratic")):
    r = cohens_kappa(tab, wt=wt)
    print(f"{name:10s} kappa = {r.kappa:.3f}, SE = {r.std_kappa:.3f}, "
          f"95% CI ({r.kappa_low:.3f}, {r.kappa_upp:.3f})")
```

**What the code does**

- `np.repeat` expands the count table back to one row per patient, because `cohen_kappa_score` takes two vectors of ratings (it has no table input).
- `weights=None` is unweighted; `"linear"` and `"quadratic"` are the weighted versions. The category order is taken from the sorted labels, so the codes 1, 2, 3 must be in clinical order.
- `sklearn` returns only the point estimate. `statsmodels.stats.inter_rater.cohens_kappa` works directly on the table and also returns the SE and CI (`kappa_low`, `kappa_upp`); use `wt="linear"` or `"quadratic"` there.

**Output**

```text
n = 80, observed agreement = 0.738
sklearn kappa (unweighted): 0.603
sklearn kappa (linear): 0.678
sklearn kappa (quadratic): 0.758
unweighted kappa = 0.603, SE = 0.075, 95% CI (0.456, 0.750)
linear     kappa = 0.678, SE = 0.065, 95% CI (0.551, 0.805)
quadratic  kappa = 0.758, SE = 0.058, 95% CI (0.645, 0.871)
```

### SAS

```sas
/* One row per cell of the cross-table; count is the cell frequency */
data xray;
  input a b count;
  datalines;
1 1 22
1 2 6
1 3 1
2 1 5
2 2 20
2 3 5
3 1 0
3 2 4
3 3 17
;
run;

/* Simple kappa and weighted kappa (Cicchetti-Allison linear weights by default) */
proc freq data=xray;
  tables a*b / agree;
  weight count / zeros;
  test kappa wtkap;
run;

/* Fleiss-Cohen (quadratic) weights */
proc freq data=xray;
  tables a*b / agree(wt=fc);
  weight count / zeros;
run;
```

**What the code does**

- The data have one row per cell with a frequency `count`. `weight count / zeros` tells PROC FREQ to use `count` as the frequency; the `zeros` option keeps the empty cell (3,1) in the table, so the table stays 3 × 3.
- `tables a*b / agree` requests the McNemar-type symmetry test (Bowker's test for tables larger than 2 × 2), **simple kappa** and **weighted kappa** with 95% confidence limits.
- By default, weighted kappa uses **Cicchetti–Allison** (linear) weights. `agree(wt=fc)` switches to Fleiss–Cohen (quadratic) weights.
- `test kappa wtkap` adds the asymptotic tests of H₀: κ = 0 for the simple and weighted kappa.
- The kappa values should match R and Python. SAS's confidence limits can differ in the third decimal, because the asymptotic SE formula differs slightly between programs.

## Interpreting the output

- **Observed agreement = 59/80 = 0.738.** The radiologists give the same grade for 74% of patients.
- **Expected by chance: 0.338.** Given the row and column totals, about 34% agreement would occur by chance, so kappa rescales the extra agreement beyond that.
- **Unweighted κ = 0.60** (95% CI 0.46 to 0.75; SE 0.075). By the common Landis–Koch labels this is "moderate" to "substantial", but those labels are only rough conventions. The CI is wide, so the true agreement could be anywhere from moderate to substantial; the lower limit (0.46) is what you would compare with the level your clinical use requires.
- **Linear weighted κ = 0.68** (95% CI 0.55 to 0.80) and **quadratic weighted κ = 0.76** (95% CI 0.64 to 0.87). Weighted values are higher because most disagreements are by one grade (mild vs moderate or moderate vs severe); only one patient was rated mild by A and severe by B, and none the other way round.
- R and Python give identical estimates, SEs and CIs to the digits shown.
- **Which one to report?** For an ordinal scale, pre-specify the weighting; quadratic weights are the most common in clinical studies and are close to the ICC.

**How to report it:**
> Two radiologists independently graded severity (mild, moderate, severe) on 80 chest X-rays. Observed agreement was 73.8%. Cohen's kappa was 0.60 (95% CI 0.46 to 0.75); the linearly weighted kappa was 0.68 (95% CI 0.55 to 0.80) and the quadratically weighted kappa was 0.76 (95% CI 0.64 to 0.87).

## Check your understanding

<details class="check"><summary>1. Two pathologists agree on 95% of slides, but 95% of slides are benign and kappa is only 0.47. Is that a contradiction?</summary>

No. When one category is very common, chance agreement is already high (here about 90%), so little of the 95% is "extra". This is the kappa paradox. Report observed agreement, each rater's category frequencies, and kappa together.
</details>

<details class="check"><summary>2. A paper says "kappa = 0.60, so agreement is substantial" and stops. What is missing?</summary>

The confidence interval (0.46 to 0.75 here, which spans "moderate" to "substantial"), and the context. The Landis–Koch labels are conventions, not thresholds: κ = 0.6 may suit a screening aid but not a primary endpoint. Judge whether the **lower confidence limit** meets the level your application needs.
</details>

<details class="check"><summary>3. Severity is graded mild, moderate, severe. Why is the quadratic weighted kappa (0.76) higher than the unweighted (0.60)?</summary>

Unweighted kappa counts every disagreement as complete failure. Weighted kappa gives partial credit for near misses. Most disagreements here are one grade apart, and only one patient was graded mild by one rater and severe by the other, so weighting raises the value. Pre-specify the weights; they are not a free choice after seeing the data.
</details>

<details class="check"><summary>4. What if the p-value for H₀: κ = 0 is below 0.001? Does that show the raters agree well?</summary>

No. It only shows agreement is better than chance, which almost any real rating process achieves. The useful question is how large κ is and how precisely, so report κ with its CI.
</details>

## Common mistakes

- **Reporting only the p-value for H₀: κ = 0.** Nearly any pair of raters beats chance. Report κ and its CI.
- **Using unweighted kappa for ordinal ratings.** It treats a mild-vs-severe disagreement the same as a mild-vs-moderate one and understates agreement.
- **Choosing the weights after seeing the data.** Weighted kappa depends on the weights; pre-specify them.
- **Ignoring prevalence and bias.** With very common or very rare categories, kappa can be low even with high raw agreement (the "kappa paradox"). Also report pₒ and the marginal distributions of each rater.
- **Treating the Landis–Koch labels as thresholds.** κ = 0.6 may be good enough for screening and not for a primary endpoint.
- **Using kappa for continuous measurements** that were cut into categories. Analyse the continuous data with the [ICC](#/test/icc) or [Bland–Altman](#/test/bland-altman) instead, since categorising loses information.
- **Mismatched category order or levels** between the two raters, which silently changes the table and the weights.

## Related tests

- [Intraclass correlation coefficient (ICC)](#/test/icc): reliability of continuous ratings (quadratic weighted kappa is approximately an ICC)
- [Bland–Altman analysis](#/test/bland-altman): agreement between two continuous methods
- [McNemar's test](#/test/mcnemar): tests whether two raters differ in their rate of a binary rating
- [Chi-square test of independence](#/test/chi-square-independence): association, which is not the same as agreement
- [ROC curve, AUC, sensitivity & specificity](#/test/roc-auc): accuracy versus a reference standard

## Key takeaways

- κ = (pₒ − pₑ)/(1 − pₑ): agreement **beyond chance**, from 0 (chance) to 1 (perfect); it is negative if agreement is worse than chance.
- Raw percent agreement overstates reliability, especially when one category is very common.
- Use weighted kappa (pre-specified weights) for ordinal ratings; unweighted kappa treats all disagreements alike.
- Report κ with its 95% CI, plus pₒ and category frequencies.
- Benchmarks like "moderate" or "substantial" are rough conventions; judge κ against what the clinical use requires.

## References

- Cohen J. A coefficient of agreement for nominal scales. *Educational and Psychological Measurement.* 1960;20:37–46.
- Landis JR, Koch GG. The measurement of observer agreement for categorical data. *Biometrics.* 1977;33:159–174.
- Kottner J, et al. Guidelines for Reporting Reliability and Agreement Studies (GRRAS). *Journal of Clinical Epidemiology.* 2011;64:96–106.
- R documentation: `?DescTools::CohenKappa`, `?irr::kappa2`; scikit-learn: `cohen_kappa_score`; SAS/STAT User's Guide: The FREQ Procedure (AGREE option)
