# Mann–Whitney U (Wilcoxon rank-sum) test

> Compares two **independent** groups by ranking all observations together and asking whether values from one group tend to be larger than values from the other; it is the standard alternative to the independent t-test for skewed or ordinal outcomes.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Two independent groups with skewed, ordinal, or outlier-prone outcomes | Ordinal or skewed continuous | Independent observations; same-shape distributions only if you interpret it as a shift | Length of stay after hip replacement: ERAS vs usual care |

## Learning objectives

- Explain the Mann–Whitney test as a comparison of ranks, i.e. "does one group tend to be higher?"
- State why it is not automatically a test of medians
- Run it in R and Python and read U or W, p, the probability of superiority, and the Hodges–Lehmann shift
- Report the result with an effect size and a confidence interval

## What it is

Pool the observations from both groups, rank them from smallest to largest (tied values share the average rank), and add up the ranks that belong to each group. If the groups come from the same distribution, the rank sums are close to what chance would give. The Mann–Whitney U statistic and the Wilcoxon rank-sum statistic W are the same test written differently: U = W − n₁(n₁ + 1)/2.

**It is not simply a test of medians.** What it really tests is whether a random patient from one group tends to have a larger value than a random patient from the other: H₀ is P(X > Y) + ½P(X = Y) = ½. It reduces to "the medians differ" only when the two distributions have the same shape and spread and differ by a pure shift. With different spreads or skewness, two groups can have identical medians yet a very small p-value, or different medians and a large one. Say "stochastically larger" or "tends to be higher", and show the shape of the data.

- **Test statistic:** U (or W); the p-value comes from the exact null distribution (small samples, no ties) or a normal approximation with tie correction
- **Effect size:** the **probability of superiority** (also called the common-language effect size, or AUC) = U/(n₁n₂); the **Hodges–Lehmann shift** (median of all n₁n₂ pairwise differences) with a CI, which is a valid location-shift estimate if the shapes are similar; rank-biserial correlation r = 2·U/(n₁n₂) − 1

## The intuition

Ranks turn the data into a **signal-to-noise** problem that does not care about outliers. Pool all 20 patients, order them by length of stay, and give each a rank. If the two arms are no different, high and low ranks should be sprinkled evenly between them. The signal is how lopsided the rank sums are (the standard-care group collects most of the high ranks); the noise is how lopsided they could be by chance alone.

Analogy: a race between two running clubs. Rather than comparing average finish times (one runner who got lost for three hours would wreck the mean), you line everyone up by finishing place and ask whether one club's runners tend to finish ahead. A runner who finishes last counts as "last", whether they were 10 minutes or 3 hours behind.

In the plot, notice the 21-day and 14-day stays in the standard-care group. They would drag a mean upward but count only as the two highest ranks. Most ERAS values cluster between 4 and 9 days while standard care sits higher. About 81% of random patient pairs (one per arm) have the standard-care patient staying longer. That 81% is the probability of superiority.

<div class="viz" data-viz="dots" data-unit="length of stay, days" data-groups='{"Standard care":[7,9,6,12,8,10,14,7,9,21],"ERAS":[5,6,4,8,7,5,9,6,11,6]}' data-caption="Standard care stays tend to be longer, and the right tail is long. Ranks tame that tail."></div>

With 10 per group and ties, the p-value uses a normal approximation. The observed standardized statistic is about z = 2.32.

<div class="viz" data-viz="z" data-stat="2.32"></div>

**What to notice:** z = 2.32 falls in the shaded tails, giving p = 0.020. Be careful what this says: it says one arm's values tend to be larger. Here standard care is also more spread out (the 21-day stay widens it), so the arms are not a clean shift of one another and the 3-day Hodges–Lehmann shift is only a rough summary; the probability of superiority is the safer description.

## When to use it

Use the Mann–Whitney test when **all** of these are true:

- There are **two independent groups** (treatment vs control, men vs women)
- The outcome is **ordinal or continuous but skewed** (length of stay, hospital costs, symptom scores, viral load)
- Observations are independent, with one observation per subject
- The sample is small or outliers make the t-test unreliable

**Clinical scenario.** A hospital compares length of stay (days) after hip replacement between patients on usual care (n = 10) and patients on an enhanced-recovery pathway (ERAS, n = 10). Stay lengths are right-skewed: most patients leave in 5–10 days, but a few have long, complicated admissions. *Do ERAS patients tend to be discharged sooner?*

**Don't use it when:**

- The outcome is roughly normal in each group and the samples are not tiny → the [independent t-test](#/test/independent-t) gives a mean difference with a CI
- The same patients are measured twice → use the [Wilcoxon signed-rank test](#/test/wilcoxon-signed-rank)
- There are **three or more** groups → use the [Kruskal–Wallis test](#/test/kruskal-wallis)
- You must adjust for baseline covariates → use [ANCOVA](#/test/ancova) or a regression model
- You need to claim a difference in **medians** specifically, and the group shapes differ → compare medians directly (for example with quantile regression or a bootstrap CI)

## Assumptions

1. **Observations are independent**, both within and between groups. Check the design (no repeated patients, no clustering by site without accounting for it).
2. **The outcome is at least ordinal**, so ranking makes sense.
3. **For the "shift in location" / "difference in medians" interpretation only:** the two distributions have the **same shape and spread**. Check with side-by-side box plots or histograms. Without this assumption, interpret the result as stochastic dominance, P(X > Y) ≠ ½.
4. **Ties are limited.** Many ties (typical with integer scores) reduce power slightly and force the normal approximation in R and SciPy; the tie-corrected variance handles this.

## Hypotheses

- **H₀:** P(X > Y) + ½P(X = Y) = ½: values from either group are equally likely to be the larger one (under the usual shift model: no difference in location)
- **H₁:** P(X > Y) + ½P(X = Y) ≠ ½ (two-sided; the default for regulators)

## Example data

Length of stay (days), 10 patients per arm.

| Arm | Length of stay (days) | Median |
|---|---|---|
| Standard care | 7, 9, 6, 12, 8, 10, 14, 7, 9, 21 | 9 |
| ERAS | 5, 6, 4, 8, 7, 5, 9, 6, 11, 6 | 6 |

Several values appear more than once (ties), so exact p-values are not available and the normal approximation is used in all three languages.

## How to run it

### R

```r
# Length of hospital stay (days) after hip replacement
standard <- c(7, 9, 6, 12, 8, 10, 14, 7, 9, 21)   # usual care, n = 10
eras     <- c(5, 6, 4, 8, 7, 5, 9, 6, 11, 6)      # enhanced recovery, n = 10

# Wilcoxon rank-sum (Mann-Whitney) test; ties force the normal approximation
wt <- wilcox.test(standard, eras, exact = FALSE, correct = TRUE, conf.int = TRUE)
wt

# Probability of superiority P(standard > ERAS) + 0.5 P(tie)
unname(wt$statistic) / (length(standard) * length(eras))

# Group summaries: medians differ, but the test is about the whole distribution
sapply(list(standard = standard, eras = eras), median)
```

**What the code does**

- `wilcox.test(x, y)` with two vectors and **no** `paired = TRUE` runs the rank-sum test. R calls it "Wilcoxon rank sum test", and its statistic `W` is the same as Mann–Whitney's U for the **first** sample.
- `exact = FALSE` selects the normal approximation. R's default is exact when both samples have fewer than 50 observations **and no ties**; with ties it switches (with a warning) anyway, so stating `exact = FALSE` makes the choice visible.
- `correct = TRUE` is R's default continuity correction (±0.5 in the numerator). SciPy's `use_continuity=True` matches it; SAS PROC NPAR1WAY also applies it by default.
- `conf.int = TRUE` returns the **Hodges–Lehmann** "difference in location" (median of the pairwise differences) and a CI. With the normal approximation the CI is found numerically, so you see 2.999979 rather than exactly 3 and limits of about 1.0 and 6.0.
- `W / (n1 * n2)` converts the statistic into the probability of superiority.
- The group medians are shown only to illustrate the point above: the test is *about* the whole distribution, not the medians.

**Output**

```text

	Wilcoxon rank sum test with continuity correction

data:  standard and eras
W = 81, p-value = 0.02017
alternative hypothesis: true location shift is not equal to 0
95 percent confidence interval:
 0.9999611 5.9999575
sample estimates:
difference in location 
              2.999979 

[1] 0.81
standard     eras 
       9        6 
```

### Python

```python
import numpy as np
from scipy import stats

# Length of hospital stay (days) after hip replacement
standard = np.array([7, 9, 6, 12, 8, 10, 14, 7, 9, 21])   # usual care, n = 10
eras     = np.array([5, 6, 4, 8, 7, 5, 9, 6, 11, 6])      # enhanced recovery, n = 10

# Mann-Whitney U: normal approximation with continuity correction (tie-corrected)
res = stats.mannwhitneyu(standard, eras, alternative="two-sided",
                         method="asymptotic", use_continuity=True)
print(f"U (standard) = {res.statistic:.1f}, p = {res.pvalue:.4f}")

# Probability of superiority = U / (n1 * n2)
print(f"P(standard > ERAS) + 0.5*P(tie) = {res.statistic / (len(standard) * len(eras)):.3f}")

# Hodges-Lehmann shift: median of all n1*n2 pairwise differences
diffs = (standard[:, None] - eras[None, :]).ravel()
print(f"Hodges-Lehmann shift = {np.median(diffs):.1f} days")
print("Medians:", np.median(standard), np.median(eras))
```

**What the code does**

- `stats.mannwhitneyu(x, y)` returns **U for the first sample** (here `standard`). `U₂ = n₁n₂ − U₁`; if you swap the arguments you get U₂ instead. For a two-sided test SciPy ≥ 1.7 returns U₁, not min(U₁, U₂).
- `method="asymptotic"` is the normal approximation, tie-corrected. SciPy's default `method="auto"` would use the exact method only when one of the samples has fewer than 8 observations and there are no ties, so here it would already fall back to this approximation. Stating it keeps the code explicit.
- `use_continuity=True` is the default and matches R's `correct = TRUE`. Setting `False` here (and `correct = FALSE` in R) would give identical, slightly smaller p-values in both languages.
- SciPy does not return a Hodges–Lehmann CI; the point estimate is the median of all 100 pairwise differences. If a CI is needed, use R's `conf.int = TRUE`, SAS `HL`, or a bootstrap.

**Output**

```text
U (standard) = 81.0, p = 0.0202
P(standard > ERAS) + 0.5*P(tie) = 0.810
Hodges-Lehmann shift = 3.0 days
Medians: 9.0 6.0
```

### SAS

```sas
/* One row per patient: arm and length of stay (days) */
data los;
  input arm $ days @@;
  datalines;
Standard 7  Standard 9  Standard 6  Standard 12 Standard 8
Standard 10 Standard 14 Standard 7  Standard 9  Standard 21
ERAS 5  ERAS 6  ERAS 4  ERAS 8  ERAS 7
ERAS 5  ERAS 9  ERAS 6  ERAS 11 ERAS 6
;
run;

/* Wilcoxon rank-sum test with Hodges-Lehmann shift estimate and 95% CI */
proc npar1way data=los wilcoxon hl alpha=0.05;
  class arm;
  var days;
run;
```

**What the code does**

- The data are in **long format**: one row per patient with a group variable (`arm`) and the outcome (`days`). `@@` in the INPUT statement lets several observations share one data line.
- `proc npar1way ... wilcoxon` performs the Wilcoxon rank-sum test. The output has the sum of scores (rank sums) for each class and a **normal approximation Z with a continuity correction of 0.5** (the default), plus a t approximation. The two-sided p-value should be about 0.020, as in R and SciPy. For small samples, add `exact wilcoxon;` to get the exact p-value.
- SAS lists the classes **alphabetically**, so `ERAS` is the first class. The reported rank-sum statistic is therefore for ERAS (74; Standard has 136, and 74 + 136 = 210 = 20·21/2), and the Z statistic is negative (about −2.3).
- `hl` adds the **Hodges–Lehmann estimator of shift and its asymptotic confidence limits**, computed as first class minus second class, i.e. **ERAS − Standard ≈ −3.0 days**, the opposite sign of R's (standard − ERAS = +3.0).
- Differences in default (continuity correction, exact vs approximate) between packages are the usual reason for small mismatches. If your numbers differ in the third decimal, compare these options first.

## Interpreting the output

- **W = U = 81 for standard care** (out of a maximum of 10 × 10 = 100). Standard-care patients had the longer stay in most of the 100 patient-to-patient comparisons.
- **Probability of superiority = 0.81.** If you pick one patient at random from each arm, there is an 81% chance (ties counted as half) that the usual-care patient stayed longer. A value of 0.5 would mean no difference.
- **p = 0.020 (R 0.02017; SciPy 0.0202).** Under the null hypothesis of no difference between the distributions, a rank separation this large would be uncommon. This p-value uses the normal approximation with continuity correction because of the ties.
- **Hodges–Lehmann shift = 3.0 days, 95% CI ≈ 1.0 to 6.0 days (R).** The typical patient on usual care stays about 3 days longer, with plausible values from about 1 to 6 days. The CI excludes 0 and agrees with the p-value. The sample medians (9 vs 6 days) happen to give the same 3-day difference, but that is not guaranteed in general.

**How to report it:**
> Length of stay was shorter with the ERAS pathway than with usual care (median 6 vs 9 days; Mann–Whitney U = 81, p = 0.020; Hodges–Lehmann difference 3.0 days, 95% CI 1.0 to 6.0; probability of superiority 0.81).

## Check your understanding

<details class="check"><summary>1. Two groups have the same median (10 days) but one is much more spread out and skewed. Could Mann–Whitney still give a small p-value?</summary>

Yes. The test compares whole distributions (is P(X > Y) = ½?). Equal medians with different shapes can still give P(X > Y) ≠ ½. A "difference in medians" reading is only valid when the shapes and spreads are the same.
</details>

<details class="check"><summary>2. Twelve patients have BP measured before and after a drug. Is Mann–Whitney the right test?</summary>

No. The two measurements are on the same patients, so the data are paired. Use the [Wilcoxon signed-rank test](#/test/wilcoxon-signed-rank). Treating paired data as independent wastes the pairing and loses power.
</details>

<details class="check"><summary>3. The probability of superiority is 0.81. Explain it to a clinician in one sentence.</summary>

If you pick one patient from each arm at random, the usual-care patient will have the longer stay about 81% of the time (ties count as half). 0.5 would mean no difference.
</details>

<details class="check"><summary>4. What if the 21-day stay were instead 40 days? What happens to U, p, and to a t-test?</summary>

U and p do not change: 21 and 40 are both the highest rank. A t-test would see a larger mean and a much larger SD, so its result could change a lot. This robustness to extremes is the main attraction of rank tests.
</details>

## Common mistakes

- **Describing the result as "the medians differ".** The test compares whole distributions. Use medians only as descriptive statistics and the Hodges–Lehmann shift as a location estimate under the equal-shape assumption.
- **Mixing up paired and independent data.** Independent groups need the rank-sum test (no `paired = TRUE` in R); paired data need the signed-rank test.
- **Forgetting which group the U statistic belongs to.** R's W and SciPy's U refer to the first sample, SAS reports the first class (alphabetical). Swapping arguments flips U to n₁n₂ − U and the sign of the shift.
- **Ignoring ties and software defaults.** Continuity correction and exact-vs-asymptotic choices differ between R, SciPy and SAS, so set them explicitly and report them.
- **Using it for clustered data** (e.g. several eyes or several visits per patient). The independence assumption fails; use a mixed model or a cluster-adjusted rank method.
- **Reporting only the p-value.** Add the probability of superiority and/or the Hodges–Lehmann shift with a CI.

## Related tests

- [Independent two-sample t-test](#/test/independent-t): parametric alternative that compares means
- [Wilcoxon signed-rank test](#/test/wilcoxon-signed-rank): paired-data analogue
- [Kruskal–Wallis test](#/test/kruskal-wallis): three or more independent groups
- [ROC curve, AUC](#/test/roc-auc): the AUC is the same quantity as the probability of superiority
- [Shapiro–Wilk test](#/test/shapiro-wilk): check normality before choosing the t-test

## Key takeaways

- Mann–Whitney ranks all observations together and asks whether one group tends to have higher values.
- It tests P(X > Y) = ½, **not** "medians are equal" unless the two distributions have the same shape.
- Report the probability of superiority and/or the Hodges–Lehmann shift with a CI, not just p.
- Ties and software defaults (continuity correction, exact vs approximate) change p slightly; state your choices.
- Paired data need the Wilcoxon signed-rank test instead.

## References

- Hollander M, Wolfe DA, Chicken E. *Nonparametric Statistical Methods.* 3rd ed. Wiley; 2014. Ch. 4.
- Fay MP, Proschan MA. Wilcoxon–Mann–Whitney or t-test? On assumptions for hypothesis tests and multiple interpretations of decision rules. *Statistics Surveys.* 2010;4:1–39.
- Documentation: R `?wilcox.test`; SciPy `scipy.stats.mannwhitneyu`; SAS/STAT User's Guide, The NPAR1WAY Procedure
