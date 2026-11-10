# Friedman test

> Compares **three or more related measurements** on the same subjects (repeated time points, or conditions in a crossover) using within-subject ranks; it is the nonparametric alternative to one-way repeated-measures ANOVA.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Same subjects measured under 3 or more conditions, nonparametric | Ordinal or non-normal continuous | Independent subjects; complete blocks; no carry-over | Pain VAS under three analgesics in a crossover |

## Learning objectives

- Recognize repeated-measures designs with three or more conditions where ranks suit the data
- Explain Friedman as a comparison of within-subject mean ranks
- Run it in R and Python and read Q, p, Kendall's W, and mean ranks
- Follow up a significant omnibus test with adjusted pairwise comparisons

## What it is

The Friedman test treats each subject as a **block**. Within every subject, the k measurements are ranked from 1 to k, and then the ranks are averaged for each condition across subjects. If the conditions have no effect, each condition should receive about the same mean rank, (k + 1)/2, because the order within a subject is then random. The Friedman statistic Q measures how far the mean ranks are from that expectation, and is compared with a χ² distribution with k − 1 degrees of freedom (a reasonable approximation with about 10 or more blocks; use exact tables for smaller designs).

Ranking within each subject removes between-subject differences, so the test uses only the *ordering* of conditions inside each person.

- **Test statistic:** Q (also written χ²_F), tie-corrected, with k − 1 df
- **Effect size:** **Kendall's W** = Q / (n(k − 1)), the coefficient of concordance from 0 (no agreement in the within-subject order) to 1 (every subject orders the conditions identically); plus mean ranks and medians per condition

## The intuition

The Friedman test works with **ranks inside each person**. For each patient we ask only: which regimen gave the lowest pain, which the middle, which the highest? Between-patient differences (one person always scores high, another low) vanish because ranking happens within the person.

Signal vs noise: if the regimens did nothing, the order within each patient would be a random shuffle, and each regimen would earn a mean rank near 2 (the middle of 1, 2, 3). The **signal** is how far the mean ranks drift from 2. The **noise** is how far they would drift by chance from random shuffling.

An everyday picture: ten judges each rank three wines. If the wines are equal, no wine should consistently win. If most judges put wine C first and A last, the mean ranks pull apart.

In our data, the mean ranks are A 2.8, B 2.0, C 1.2. In the plot, most patients have their highest VAS in regimen A and their lowest in C, even though the raw values overlap across patients (patient 4 on C scored 69, higher than patient 3 on A scored 55).

<div class="viz" data-viz="dots" data-unit="mm, pain VAS" data-groups='{"Regimen A":[62,70,55,81,67,73,58,77,64,69],"Regimen B":[50,62,50,72,51,76,52,63,60,71],"Regimen C":[44,54,32,69,52,49,47,72,58,52]}' data-caption="Raw values overlap because patients differ in baseline pain. The Friedman test looks only at the order within each patient."></div>

Our statistic is Q = 12.8 with 2 df. Under H₀ it follows the χ² curve below. Ours sits well into the upper tail.

<div class="viz" data-viz="chisq" data-stat="12.8" data-df="2"></div>

## When to use it

Use the Friedman test when **all** of these are true:

- The **same subjects** are measured under **three or more** conditions or time points (complete blocks, no missing values)
- The outcome is **ordinal or continuous but non-normal**
- Subjects are independent of each other
- The conditions are exchangeable under H₀ (for time points, there is no assumption that the effect is linear)

**Clinical scenario.** In a three-period crossover study, 10 patients with chronic pain receive three analgesic regimens (A, B and C) in randomised order, with a washout between periods. After each period they rate their pain on a 0–100 mm visual analogue scale (VAS). Pain scores are not normally distributed. *Do the three regimens differ in pain, and which pairs differ?*

**Don't use it when:**

- The data are approximately normal and sphericity is acceptable → use [repeated-measures ANOVA](#/test/repeated-measures-anova)
- There are only **two** conditions → use the [Wilcoxon signed-rank test](#/test/wilcoxon-signed-rank)
- The groups are **independent** → use the [Kruskal–Wallis test](#/test/kruskal-wallis)
- Some subjects have **missing** measurements or the visits are unbalanced → use a [mixed model (MMRM)](#/test/mixed-models), which uses all available data
- The outcome is binary (responder yes/no) at each time → use Cochran's Q test, or [GEE](#/test/gee)

## Assumptions

1. **Subjects (blocks) are independent of each other.** Check the design.
2. **Each subject contributes one value under every condition**, with no missing cells. Check with a data listing; incomplete subjects are dropped by R's `friedman.test`.
3. **The outcome is at least ordinal**, so that values can be ranked within a subject.
4. **No carry-over or period effects** (for a crossover) that systematically favour one condition. Check the design, randomised order and washout.
5. **Enough blocks** for the χ² approximation (about 10 or more). With fewer, use exact p-values.

## Hypotheses

- **H₀:** the k conditions have the same distribution within subjects (the mean ranks are equal)
- **H₁:** at least one condition tends to give higher or lower values than another

## Example data

Pain VAS (mm), 10 patients, each on three regimens.

| Patient | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| Regimen A | 62 | 70 | 55 | 81 | 67 | 73 | 58 | 77 | 64 | 69 |
| Regimen B | 50 | 62 | 50 | 72 | 51 | 76 | 52 | 63 | 60 | 71 |
| Regimen C | 44 | 54 | 32 | 69 | 52 | 49 | 47 | 72 | 58 | 52 |

## How to run it

### R

```r
# Pain VAS (mm) in 10 patients, each tried three analgesic regimens in a crossover
pain <- matrix(c(62, 50, 44,
                 70, 62, 54,
                 55, 50, 32,
                 81, 72, 69,
                 67, 51, 52,
                 73, 76, 49,
                 58, 52, 47,
                 77, 63, 72,
                 64, 60, 58,
                 69, 71, 52), ncol = 3, byrow = TRUE,
               dimnames = list(paste0("P", 1:10), c("A", "B", "C")))

# Friedman test: rows = blocks (patients), columns = treatments
ft <- friedman.test(pain)
ft

# Kendall's W = Q / (n * (k - 1))
unname(ft$statistic) / (nrow(pain) * (ncol(pain) - 1))

# Mean ranks within patient
colMeans(t(apply(pain, 1, rank)))

# Post hoc: paired Wilcoxon tests, exact, Holm-adjusted
long <- data.frame(y = as.vector(pain), trt = rep(colnames(pain), each = 10))
pairwise.wilcox.test(long$y, long$trt, paired = TRUE, p.adjust.method = "holm")
```

**What the code does**

- `friedman.test(pain)` takes a **matrix with one row per subject (block) and one column per condition**. A common error is to supply the transpose. You can also use the formula form `friedman.test(y ~ trt | patient, data = long)` on long data.
- R ranks within each row, averages ties, and applies the **tie correction** to Q, the same as SciPy and SAS.
- `Q / (n(k − 1))` converts the statistic to **Kendall's W**.
- `apply(pain, 1, rank)` ranks within each patient; `colMeans` of the transposed result gives the mean rank per regimen (1 = lowest pain).
- The omnibus test does not say which regimens differ. `pairwise.wilcox.test(..., paired = TRUE, p.adjust.method = "holm")` runs paired Wilcoxon signed-rank tests for each pair, with **Holm** adjustment. This works because the patient order is the same within each regimen. It reports "exact" tests here (no ties among differences). For a rank-based all-pairs procedure that uses the Friedman ranks directly, use the **Nemenyi** test (`PMCMRplus::frdAllPairsNemenyiTest`) or **Conover** test.

**Output**

```text

	Friedman rank sum test

data:  pain
Friedman chi-squared = 12.8, df = 2, p-value = 0.001662

[1] 0.64
  A   B   C 
2.8 2.0 1.2 

	Pairwise comparisons using Wilcoxon signed rank exact test 

data:  long$y and long$trt 

  A      B     
B 0.0195 -     
C 0.0059 0.0488

P value adjustment method: holm 
```

### Python

```python
import numpy as np
from itertools import combinations
from scipy import stats
from statsmodels.stats.multitest import multipletests

# Pain VAS (mm) in 10 patients, each tried three analgesic regimens in a crossover
A = np.array([62, 70, 55, 81, 67, 73, 58, 77, 64, 69])
B = np.array([50, 62, 50, 72, 51, 76, 52, 63, 60, 71])
C = np.array([44, 54, 32, 69, 52, 49, 47, 72, 58, 52])

# Friedman test: one array per treatment, same patient order
res = stats.friedmanchisquare(A, B, C)
n, k = len(A), 3
print(f"Friedman chi2 = {res.statistic:.2f}, df = {k - 1}, p = {res.pvalue:.5f}")
print(f"Kendall's W = {res.statistic / (n * (k - 1)):.2f}")

# Mean ranks within patient
ranks = stats.rankdata(np.column_stack([A, B, C]), axis=1)
print("Mean ranks (A, B, C):", ranks.mean(axis=0))

# Post hoc: paired exact Wilcoxon tests, Holm-adjusted
data = {"A": A, "B": B, "C": C}
pairs = list(combinations(data, 2))
raw = [stats.wilcoxon(data[a], data[b], method="exact").pvalue for a, b in pairs]
adj = multipletests(raw, method="holm")[1]
for (a, b), r, p in zip(pairs, raw, adj):
    print(f"{a} vs {b}: raw p = {r:.4f}, Holm p = {p:.4f}")
```

**What the code does**

- `stats.friedmanchisquare(A, B, C)` takes **one array per condition**, each in the same subject order. It requires at least 3 conditions. It handles ties with a correction and returns the χ² statistic and an asymptotic p-value.
- Kendall's W is computed by hand as `Q / (n × (k − 1))`.
- `stats.rankdata(..., axis=1)` ranks within each row (patient), giving the mean ranks.
- The post-hoc step runs exact paired `stats.wilcoxon` tests for every pair and adjusts them with `multipletests(method="holm")` from statsmodels. This is the same approach as R's `pairwise.wilcox.test`. `scikit_posthocs.posthoc_nemenyi_friedman` provides the Nemenyi alternative.

**Output**

```text
Friedman chi2 = 12.80, df = 2, p = 0.00166
Kendall's W = 0.64
Mean ranks (A, B, C): [2.8 2.  1.2]
A vs B: raw p = 0.0098, Holm p = 0.0195
A vs C: raw p = 0.0020, Holm p = 0.0059
B vs C: raw p = 0.0488, Holm p = 0.0488
```

### SAS

```sas
/* Wide format: one row per patient, pain VAS (mm) under regimens A, B, C */
data wide;
  input patient a b c;
  datalines;
1 62 50 44
2 70 62 54
3 55 50 32
4 81 72 69
5 67 51 52
6 73 76 49
7 58 52 47
8 77 63 72
9 64 60 58
10 69 71 52
;
run;

/* Long format: one row per patient x regimen */
data long;
  set wide;
  length trt $1;
  trt = 'A'; pain = a; output;
  trt = 'B'; pain = b; output;
  trt = 'C'; pain = c; output;
  keep patient trt pain;
run;

/* Rank the three values within each patient (ties get mid-ranks) */
proc rank data=long out=ranked;
  by patient;
  var pain;
  ranks rk;
run;

/* Friedman test = CMH "row mean scores" statistic, patients as strata */
proc freq data=ranked;
  tables patient*trt*rk / cmh2 scores=table noprint;
run;
```

**What the code does**

- PROC FREQ is the SAS route to the Friedman test. The data must be in **long format**, with one row per patient and regimen, which the second DATA step creates from the wide data.
- `proc rank ... by patient` ranks the three pain values **within each patient** (`ranks rk`; ties get mid-ranks).
- `tables patient*trt*rk / cmh2 scores=table` builds a three-way table with **patient as the stratum**, and `cmh2` requests the CMH statistics for nonzero correlation and **row mean scores differ**. With patient-wise ranks as scores, the **"Row Mean Scores Differ"** statistic (df = k − 1 = 2) is the **Friedman chi-square**; it should reproduce Q = 12.8 and p ≈ 0.0017. The other CMH line (nonzero correlation, df = 1) is not the Friedman test; ignore it.
- `noprint` suppresses the individual stratum tables (one per patient); the CMH summary tables are still printed.
- For post-hoc comparisons, compute paired differences and run `PROC UNIVARIATE` (signed-rank test) for each pair, then adjust the three p-values with `PROC MULTTEST ... holm`.

## Interpreting the output

- **Friedman χ² = 12.80, df = 2, p = 0.0017.** There is evidence that pain differs across the three regimens. R and SciPy agree.
- **Mean ranks: A 2.8, B 2.0, C 1.2** (1 = lowest pain, 3 = highest). Regimen A tends to give the most pain and regimen C the least; B is in the middle.
- **Kendall's W = 0.64.** Concordance among patients is substantial: most patients order the regimens similarly, but not all (patients 6 and 10 had B higher than A, and patients 5 and 8 had C higher than B).
- **Pairwise Wilcoxon, Holm-adjusted:** A vs C p = 0.0059 and A vs B p = 0.0195 are significant; B vs C p = 0.0488 is borderline (just under 0.05 after adjustment). Pain ordering A > B > C is supported, but the B vs C difference is the least certain.
- Follow up with the **size** of the differences (for example the Hodges–Lehmann median paired difference in mm with a CI for each pair), not just p-values.

**How to report it:**
> Pain differed across the three regimens (Friedman χ²(2) = 12.80, p = 0.002; Kendall's W = 0.64; mean ranks 2.8, 2.0 and 1.2 for A, B and C). Pairwise Wilcoxon signed-rank tests with Holm adjustment showed lower pain with C than A (p = 0.006) and with B than A (p = 0.020); C vs B was borderline (p = 0.049).

## Check your understanding

<details class="check"><summary>1. Twelve patients each try four diets in random order and a 0 to 10 satisfaction score is recorded. Which test fits?</summary>

The Friedman test: same subjects, four conditions, an ordinal score. Independent groups would call for Kruskal–Wallis; normal continuous data with sphericity would allow repeated-measures ANOVA.
</details>

<details class="check"><summary>2. The Friedman test gives p = 0.002. Does that tell you which regimens differ?</summary>

No. It only says at least one regimen differs. Follow up with pairwise comparisons (paired Wilcoxon with Holm adjustment, or Nemenyi), and report effect sizes with CIs for each pair.
</details>

<details class="check"><summary>3. Patient 5 is missing the regimen C measurement. What does R's friedman.test do, and what should you consider?</summary>

R drops the whole patient. With substantial missingness this wastes data and can bias results, so consider a [mixed model](#/test/mixed-models) that uses all available measurements.
</details>

<details class="check"><summary>4. What if every patient ranked the regimens identically (A highest, then B, then C)? What would Kendall's W be?</summary>

W = 1, perfect agreement among patients. Our observed W = 0.64 means substantial but imperfect agreement. If patients ordered the regimens at random, W would be near 0.
</details>

## Common mistakes

- **Wrong matrix orientation in R.** Rows must be subjects and columns conditions. Transposing silently gives a different (wrong) test.
- **Missing values.** R drops entire subjects with any missing value; SciPy does not handle missing values and will give meaningless results. Use a [mixed model](#/test/mixed-models) when data are incomplete.
- **Running unadjusted pairwise tests.** With k conditions, there are k(k − 1)/2 pairs; adjust the p-values.
- **Treating order effects as treatment effects.** In a crossover, check period and carry-over effects.
- **Using Friedman with independent groups.** Use [Kruskal–Wallis](#/test/kruskal-wallis) instead.
- **Reporting only the p-value.** Give mean ranks, Kendall's W, and pairwise effect sizes with CIs.

## Related tests

- [Repeated-measures ANOVA](#/test/repeated-measures-anova): parametric counterpart
- [Wilcoxon signed-rank test](#/test/wilcoxon-signed-rank): the two-condition case
- [Kruskal–Wallis test](#/test/kruskal-wallis): independent groups
- [Linear mixed-effects models (MMRM)](#/test/mixed-models): handles missing visits
- [Multiple-comparison corrections](#/test/multiple-comparisons): Holm, Bonferroni, Benjamini–Hochberg

## Key takeaways

- Friedman ranks the k conditions **within each subject** and compares mean ranks across conditions.
- Under H₀ each condition's mean rank is about (k + 1)/2; Q is compared with χ² (k − 1 df).
- Kendall's W = Q / (n(k − 1)) measures agreement among subjects, from 0 to 1.
- A significant result needs adjusted pairwise follow-ups (Holm, Nemenyi).
- Subjects with any missing value are dropped, so consider a mixed model.

## References

- Friedman M. The use of ranks to avoid the assumption of normality implicit in the analysis of variance. *J Am Stat Assoc.* 1937;32:675–701.
- Hollander M, Wolfe DA, Chicken E. *Nonparametric Statistical Methods.* 3rd ed. Wiley; 2014. Ch. 7.
- Documentation: R `?friedman.test`; SciPy `scipy.stats.friedmanchisquare`; SAS/STAT User's Guide, The FREQ Procedure — CMH statistics
