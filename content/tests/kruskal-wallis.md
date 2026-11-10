# Kruskal–Wallis test

> Compares **three or more independent groups** by ranks, testing whether at least one group tends to give higher or lower values than the others; it is the rank-based alternative to one-way ANOVA.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Compare 3+ independent groups when data are skewed or ordinal | Ordinal or skewed continuous | Independent observations; similar shapes if you want to talk about medians | IL-6 across placebo, low dose and high dose |

## Learning objectives

- Explain how Kruskal–Wallis works on ranks and what its H statistic measures
- State correctly what it tests (stochastic dominance), and when it is a comparison of medians
- Run it in R and Python, follow up with Dunn's post-hoc test, and read the adjusted p-values
- Report H, ε², and mean ranks without calling it a "median test"

## What it is

Kruskal–Wallis extends the [Mann–Whitney test](#/test/mann-whitney) to k groups. All observations are pooled and ranked, and the mean rank of each group is compared with the overall mean rank (N + 1)/2. If the groups come from the same distribution, the group mean ranks should be close to each other. The statistic H grows as the group mean ranks spread apart, and, for reasonably sized groups, H follows a χ² distribution with k − 1 degrees of freedom under H₀.

Like Mann–Whitney, it is a test of **stochastic dominance**, and it is a test of medians only if the group distributions have the same shape. It is an **omnibus** test: a small p-value says that *some* group differs, not which one. Follow it with a post-hoc procedure such as **Dunn's test**.

- **Test statistic:** H, tie-corrected, compared with χ²(k − 1)
- **Effect size:** epsilon-squared ε² = H / (N − 1) (proportion of rank variability explained; η²_H = (H − k + 1)/(N − k) is a common alternative), plus group medians and mean ranks. For pairwise contrasts, report Hodges–Lehmann shifts or probabilities of superiority.

## The intuition

Imagine lining up all 24 patients from lowest to highest IL-6 and handing out rank numbers 1 to 24. If the drug did nothing, the three arms would be shuffled randomly through the line, and each arm's **average rank** would be near the overall average of 12.5. The signal is how far the arms' average ranks sit from 12.5; the noise is how far they would wander by luck when 24 numbers are dealt into three groups of 8. H is that signal-to-noise ratio, and under H₀ it follows a χ² curve with k − 1 = 2 degrees of freedom.

Ranks matter because they tame outliers. One placebo patient has 22.5 pg/mL, which would drag a mean and inflate an ANOVA's variance. On the rank scale that patient is just "24th", one step above the next, no matter how extreme the value.

<div class="viz" data-viz="dots" data-unit="IL-6, pg/mL" data-groups='{"Placebo":[12.4,8.9,15.2,9.8,22.5,11.1,14.0,9.5],"Low dose":[8.1,6.4,9.9,5.8,7.2,10.4,6.9,8.8],"High dose":[4.9,6.1,3.8,5.5,7.0,4.2,6.1,5.0]}' data-caption="Notice that the placebo values sit highest and the high-dose values lowest, with almost no overlap between those two arms. The 22.5 outlier would dominate a mean but counts only as one rank. The arms also have different spreads, so we describe the result as a tendency toward higher values, not as a difference in medians."></div>

Our H = 17.12 on 2 degrees of freedom. The curve below shows what H does under H₀ (all arms the same). Ours is far into the upper tail, and only the upper tail counts for this test.

<div class="viz" data-viz="chisq" data-stat="17.12" data-df="2"></div>

## When to use it

Use Kruskal–Wallis when **all** of these are true:

- There are **three or more independent groups** (dose arms, disease stages, study sites)
- The outcome is **ordinal, or continuous but skewed** (cytokine levels, symptom scores, length of stay)
- Observations are independent, with one per subject
- Groups are not so small that the χ² approximation is poor (there is no sharp cutoff; with only a handful of observations per group, prefer an exact or permutation version)

**Clinical scenario.** In a phase 2 trial, 24 patients with rheumatoid arthritis are randomised to placebo, a low dose or a high dose of an anti-inflammatory drug (8 per arm). Serum interleukin-6 (IL-6, pg/mL) is measured at week 12. IL-6 is strongly right-skewed (one placebo patient has 22.5 pg/mL). *Does IL-6 differ across arms, and which arms differ?*

**Don't use it when:**

- The outcome is approximately normal with similar variances → use [one-way ANOVA](#/test/one-way-anova)
- There are only **two** groups → use the [Mann–Whitney U test](#/test/mann-whitney)
- The same patients are measured under every condition → use the [Friedman test](#/test/friedman)
- You need to adjust for covariates or have more than one factor → use [ANCOVA](#/test/ancova) or a regression model
- The groups are naturally ordered (dose-response) and you want a trend test → use a Jonckheere–Terpstra test, which has more power against ordered alternatives

## Assumptions

1. **Observations are independent**, within and between groups. Check the design (no repeated measures, no clustering).
2. **The outcome is at least ordinal.**
3. **For the "difference in medians" interpretation:** the group distributions have **similar shapes and spreads**. Check with side-by-side box plots. Otherwise interpret the result as "some group tends to give larger values".
4. **Group sizes are not tiny.** The χ² approximation can be inaccurate when groups are very small (roughly five or fewer observations is a common warning sign, not a hard threshold); use exact or permutation methods then. Heavy ties lower power slightly and are handled by the tie correction.

## Hypotheses

- **H₀:** all k groups come from the same distribution (no group tends to give larger values than another)
- **H₁:** at least one group tends to give larger or smaller values than at least one other group

## Example data

Serum IL-6 (pg/mL) at week 12, 8 patients per arm.

| Arm | IL-6 (pg/mL) |
|---|---|
| Placebo | 12.4, 8.9, 15.2, 9.8, 22.5, 11.1, 14.0, 9.5 |
| Low dose | 8.1, 6.4, 9.9, 5.8, 7.2, 10.4, 6.9, 8.8 |
| High dose | 4.9, 6.1, 3.8, 5.5, 7.0, 4.2, 6.1, 5.0 |

## How to run it

### R

```r
# Serum IL-6 (pg/mL) after 12 weeks, three arms, n = 8 each
il6 <- c(12.4, 8.9, 15.2, 9.8, 22.5, 11.1, 14.0, 9.5,    # placebo
          8.1, 6.4, 9.9, 5.8, 7.2, 10.4, 6.9, 8.8,        # low dose
          4.9, 6.1, 3.8, 5.5, 7.0, 4.2, 6.1, 5.0)         # high dose
arm <- factor(rep(c("Placebo", "Low", "High"), each = 8),
              levels = c("Placebo", "Low", "High"))

# Omnibus Kruskal-Wallis test (tie-corrected H, chi-square approximation)
kt <- kruskal.test(il6 ~ arm)
kt

# Mean rank per arm
tapply(rank(il6), arm, mean)

# Epsilon-squared effect size = H / (n - 1)
unname(kt$statistic) / (length(il6) - 1)

# Dunn's post-hoc test with Holm adjustment
library(DescTools)
DunnTest(il6 ~ arm, method = "holm")
```

**What the code does**

- `kruskal.test(il6 ~ arm)` uses a formula with the outcome on the left and the grouping factor on the right. `arm` must be a factor (or character); if the arms were coded 1, 2, 3 and stored as numbers, you would get the wrong grouping.
- R's statistic is **tie-corrected** by default; there is no `correct` argument to set. SciPy and SAS use the same correction, so H agrees across languages.
- The p-value comes from the χ² approximation with k − 1 df; there is no exact option in base R. For small groups, use a permutation version (for example the `coin` package).
- `tapply(rank(il6), arm, mean)` shows the **mean rank** per arm, which is what the test actually compares.
- `H / (N − 1)` is ε².
- `DunnTest(..., method = "holm")` from **DescTools** runs **Dunn's post-hoc test**: z-tests on the differences between mean ranks, using the pooled ranking and a tie-corrected variance, with p-values adjusted by Holm's method (`"bonferroni"` and `"BH"` are alternatives). The `FSA::dunnTest` and `dunn.test::dunn.test` functions do the same.

**Output**

```text

	Kruskal-Wallis rank sum test

data:  il6 by arm
Kruskal-Wallis chi-squared = 17.122, df = 2, p-value = 0.0001914

Placebo     Low    High 
 19.750  12.625   5.125 
[1] 0.7444541

 Dunn's test of multiple comparisons using rank sums : holm  

             mean.rank.diff    pval    
Low-Placebo          -7.125 0.06771 .  
High-Placebo        -14.625 0.00011 ***
High-Low             -7.500 0.06771 .  
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Python

```python
import numpy as np
from scipy import stats
from itertools import combinations

# Serum IL-6 (pg/mL) after 12 weeks, three arms, n = 8 each
groups = {
    "Placebo": [12.4, 8.9, 15.2, 9.8, 22.5, 11.1, 14.0, 9.5],
    "Low":     [8.1, 6.4, 9.9, 5.8, 7.2, 10.4, 6.9, 8.8],
    "High":    [4.9, 6.1, 3.8, 5.5, 7.0, 4.2, 6.1, 5.0],
}

# Omnibus Kruskal-Wallis (tie-corrected H, chi-square approximation)
H, p = stats.kruskal(*groups.values())
N = sum(len(g) for g in groups.values())
print(f"H = {H:.4f}, df = {len(groups) - 1}, p = {p:.5f}")
print(f"epsilon^2 = {H / (N - 1):.3f}")

# Dunn's post-hoc test (tie-corrected), Holm-adjusted
allv = np.concatenate(list(groups.values()))
ranks = stats.rankdata(allv)
_, counts = np.unique(allv, return_counts=True)
tie = (counts**3 - counts).sum() / (12 * (N - 1))
sigma2 = N * (N + 1) / 12 - tie
meanrank, start = {}, 0
for k, g in groups.items():
    meanrank[k] = ranks[start:start + len(g)].mean(); start += len(g)
print({k: round(float(v), 3) for k, v in meanrank.items()})
rows = []
for a, b in combinations(groups, 2):
    z = (meanrank[a] - meanrank[b]) / np.sqrt(sigma2 * (1/len(groups[a]) + 1/len(groups[b])))
    rows.append((a, b, z, 2 * stats.norm.sf(abs(z))))
raw = np.array([r[3] for r in rows])
order = np.argsort(raw)
holm = np.empty(len(raw)); run = 0
for rank_i, idx in enumerate(order):
    run = max(run, (len(raw) - rank_i) * raw[idx]); holm[idx] = min(1, run)
for (a, b, z, pr), ph in zip(rows, holm):
    print(f"{a} vs {b}: z = {z:.3f}, raw p = {pr:.5f}, Holm p = {ph:.4f}")
```

**What the code does**

- `stats.kruskal(*groups.values())` takes one array per group and returns H and the p-value. It applies the **tie correction** automatically, so it matches R and SAS.
- ε² is computed by hand as `H / (N − 1)`.
- SciPy has no Dunn's test. The code computes it directly from the pooled ranks: z = (R̄ᵢ − R̄ⱼ) / √(σ² (1/nᵢ + 1/nⱼ)), where σ² = N(N + 1)/12 − Σ(t³ − t)/(12(N − 1)) is the tie-corrected variance of a single rank. Two-sided p-values come from the standard normal, and the **Holm** adjustment is applied by hand. In practice, you can use `scikit_posthocs.posthoc_dunn(..., p_adjust="holm")` for the same result.
- Adjust p-values for the number of pairwise comparisons (here 3). Without adjustment, Placebo vs Low (p = 0.044) and Low vs High (p = 0.034) would look significant at 0.05; after Holm correction they are not.

**Output**

```text
H = 17.1224, df = 2, p = 0.00019
epsilon^2 = 0.744
{'Placebo': 19.75, 'Low': 12.625, 'High': 5.125}
Placebo vs Low: z = 2.016, raw p = 0.04383, Holm p = 0.0677
Placebo vs High: z = 4.137, raw p = 0.00004, Holm p = 0.0001
Low vs High: z = 2.122, raw p = 0.03386, Holm p = 0.0677
```

### SAS

```sas
/* IL-6 (pg/mL) at week 12: 8 patients per arm */
data trial;
  input arm $ il6 @@;
  datalines;
Placebo 12.4 Placebo 8.9  Placebo 15.2 Placebo 9.8
Placebo 22.5 Placebo 11.1 Placebo 14.0 Placebo 9.5
Low 8.1  Low 6.4  Low 9.9  Low 5.8
Low 7.2  Low 10.4 Low 6.9  Low 8.8
High 4.9 High 6.1 High 3.8 High 5.5
High 7.0 High 4.2 High 6.1 High 5.0
;
run;

/* WILCOXON on 3+ groups gives the Kruskal-Wallis test (tie-corrected).
   DSCF adds pairwise Dwass-Steel-Critchlow-Fligner comparisons. */
proc npar1way data=trial wilcoxon dscf;
  class arm;
  var il6;
run;
```

**What the code does**

- The data are in **long format**: one row per patient (`arm`, `il6`). The `@@` in the INPUT statement lets several observations sit on one data line.
- With three or more classes, `proc npar1way ... wilcoxon` reports the **Kruskal–Wallis test**: a table of rank-score sums and mean scores per class, then "Kruskal-Wallis Test" with the chi-square and its df and p-value. SAS applies the tie correction.
- `dscf` adds pairwise two-sample comparisons by the **Dwass–Steel–Critchlow–Fligner** method, a rank-based all-pairs procedure that controls the family-wise error rate. It is closely related to, but not identical with, Dunn's test (Dunn's uses the pooled ranks from all groups; DSCF re-ranks each pair). SAS/STAT in recent releases also offers a `dunn` option on the PROC NPAR1WAY statement; check the documentation for your version. Otherwise, run pairwise Wilcoxon tests and adjust with `PROC MULTTEST`.
- Add `exact wilcoxon;` for an exact p-value in small samples (computation can be slow).

## Interpreting the output

- **H = 17.12, df = 2, p = 0.0002.** There is strong evidence that IL-6 is not distributed identically across the three arms. R and SciPy agree to all printed digits.
- **Mean ranks: Placebo 19.75, Low 12.63, High 5.13** (among 24 patients, so the overall mean rank is 12.5). Higher ranks mean higher IL-6; the ordering placebo > low > high suggests a dose-response.
- **ε² = 0.74.** About three quarters of the rank variability is associated with arm, a very large effect. As a rough guide, 0.01 is small, 0.06 medium and 0.14 large.
- **Dunn's test, Holm-adjusted:**
  - High vs Placebo: mean-rank difference 14.6, z = 4.14, adjusted p = 0.0001, so a clear difference.
  - Low vs Placebo: z = 2.02, unadjusted p = 0.044, **adjusted p = 0.068**, not significant after correction.
  - High vs Low: z = 2.12, unadjusted p = 0.034, **adjusted p = 0.068**, not significant after correction.
  - The three-way pattern shows the high dose is separated from placebo, while the low-dose arm lies between the two and cannot be clearly distinguished from either with n = 8 per arm.

**How to report it:**
> IL-6 differed across arms (Kruskal–Wallis H(2) = 17.12, p < 0.001, ε² = 0.74; mean ranks 19.8, 12.6 and 5.1 for placebo, low and high dose). In Dunn's post-hoc tests with Holm adjustment, only high dose vs placebo was significant (adjusted p < 0.001; low vs placebo p = 0.068; high vs low p = 0.068).

## Check your understanding

<details class="check"><summary>1. A study compares pain scores (0-10) in four independent groups of patients. Which test fits: Kruskal–Wallis or Friedman?</summary>

Kruskal–Wallis. The groups contain different patients. Friedman is for the same patients measured under three or more conditions.
</details>

<details class="check"><summary>2. A report says "Kruskal–Wallis shows the group medians differ (p = 0.0002)." Is that correct?</summary>

Not strictly. Kruskal–Wallis tests whether some group tends to yield larger values than another (stochastic dominance). It is a test of medians only when the groups have the same distribution shape and spread. Here the placebo arm is much more spread out than the high-dose arm, so the safer wording is "IL-6 tended to be higher on placebo", with medians reported as descriptive statistics.
</details>

<details class="check"><summary>3. Dunn's test gives an unadjusted p = 0.034 for high vs low dose, but the Holm-adjusted p is 0.068. Which should be reported, and what does it mean?</summary>

Report the adjusted p-value, because three pairwise comparisons were made. After accounting for multiplicity, there is not enough evidence to separate the low and high doses with 8 patients per arm. This does not mean they are the same; the trial is too small to tell.
</details>

<details class="check"><summary>4. What would happen to H if the single 22.5 pg/mL placebo value were replaced by 150 pg/mL?</summary>

Nothing. 22.5 is already the largest value (rank 24), and 150 is still the largest, so the ranks and H = 17.12 are unchanged. That robustness to outliers is the main reason to use a rank-based test. A one-way ANOVA would be strongly affected, because the outlier would inflate the placebo variance and mean.
</details>

## Common mistakes

- **Stopping at the omnibus p-value.** Report which groups differ and by how much, using a pairwise procedure that controls for multiplicity.
- **Running unadjusted pairwise Mann–Whitney tests.** With k groups there are k(k − 1)/2 comparisons, which inflates the type I error. Adjust (Holm, Bonferroni) or use Dunn/DSCF.
- **Calling it "a comparison of medians".** It compares rank distributions. Medians are valid descriptors only if the group shapes are similar.
- **Using it on paired or repeated data.** Use the [Friedman test](#/test/friedman) or a mixed model.
- **Coding groups numerically in R** without converting to a factor. Always check `levels(arm)`.
- **Using it with very small groups** without an exact or permutation p-value. Five per arm is a warning sign, not a magic cutoff.

## Related tests

- [One-way ANOVA](#/test/one-way-anova): parametric alternative
- [Mann–Whitney U test](#/test/mann-whitney): the two-group case
- [Friedman test](#/test/friedman): three or more related samples
- [Multiple-comparison corrections](#/test/multiple-comparisons): Holm, Bonferroni, Benjamini–Hochberg
- [Tukey's HSD](#/test/tukey-hsd): the parametric post-hoc counterpart

## Key takeaways

- Kruskal–Wallis ranks all observations together and compares the arms' mean ranks with the overall mean rank (N + 1)/2.
- H follows χ²(k − 1) under H₀. A small p says some group differs, not which one: follow with Dunn's test (adjusted).
- It tests stochastic dominance. It is a comparison of medians only if the shapes and spreads are similar.
- It is not the same as a "median test", and the ranks make it robust to outliers.
- Report H, df, p, ε² and mean ranks, plus adjusted pairwise results.

## References

- Kruskal WH, Wallis WA. Use of ranks in one-criterion variance analysis. *J Am Stat Assoc.* 1952;47:583–621.
- Dunn OJ. Multiple comparisons using rank sums. *Technometrics.* 1964;6:241–252.
- Documentation: R `?kruskal.test` and DescTools `DunnTest`; SciPy `scipy.stats.kruskal`; SAS/STAT User's Guide, The NPAR1WAY Procedure
