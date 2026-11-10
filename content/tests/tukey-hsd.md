# Tukey's HSD

> Compares every pair of group means after a one-way ANOVA while keeping the overall (family-wise) type I error at 5%, and gives a simultaneous confidence interval for each difference.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| All pairwise comparisons among 3+ groups with family-wise error control | Continuous outcome, grouped | Independent observations; normal residuals; roughly equal variances | Week-8 SBP across four diets, 6 patients each |

## Learning objectives

- Explain why many separate t-tests inflate the chance of a false positive
- Describe how Tukey's HSD widens intervals to control the family-wise error rate
- Run Tukey's HSD in R and Python and read the difference, simultaneous CI, and adjusted p
- Report pairwise results without treating non-significant pairs as "no difference"

## What it is

When you compare k groups pairwise, there are k(k − 1)/2 comparisons, and running separate t-tests inflates the chance of at least one false positive. Tukey's Honestly Significant Difference (HSD) controls this by using the **studentized range distribution**, which accounts for the fact that you are comparing the largest and smallest of several means. It uses the pooled within-group variance (MS_within) from the ANOVA.

- **Test statistic:** q = |x̄ᵢ − x̄ⱼ| / √(MS_within / 2 · (1/nᵢ + 1/nⱼ)), compared with the studentized range distribution with k groups and N − k df (the Tukey–Kramer version handles unequal n)
- **Output:** for each pair, the mean difference, a simultaneous 95% CI, and an **adjusted p-value**
- **Effect size:** the mean difference with its simultaneous CI; the CIs are wider than unadjusted CIs, which is the price of multiplicity control

## The intuition

Every pairwise test is a **signal-to-noise** comparison: the signal is a difference between two group means, the noise is its standard error. The new problem is that with four groups you look at six pairs. If each pair has a 5% chance of a false alarm, the chance that at least one pair falsely lights up is about 26%, not 5%.

Analogy: a screening panel of six blood tests on a healthy person. Each test flags 5% of healthy people, so more than a quarter of healthy people get at least one flag. To keep the panel as a whole honest, each test's cut-off must be stricter. Tukey does the same for pairs of means: it asks how big the *largest* gap among four means could be by chance alone, and sets the bar there. That is why Tukey intervals are wider (about ±14.7 mmHg here) than ordinary ones.

The plot shows week-8 SBP for each diet. The control group sits highest and DASH + low-sodium lowest, but the groups overlap a lot and each has only 6 patients.

<div class="viz" data-viz="dots" data-unit="systolic BP, mmHg" data-groups='{"Control":[135,133,158,135,130,144],"Low-sodium":[125,120,142,124,140,128],"DASH":[132,125,139,119,139,123],"DASH + low-sodium":[122,128,131,125,109,131]}' data-caption="The best-vs-worst gap (about 15 mmHg) is modest compared with the spread within each diet."></div>

The omnibus ANOVA gives F(3, 20) = 2.78. Here is where that falls on the null F curve.

<div class="viz" data-viz="f" data-stat="2.78" data-df1="3" data-df2="20"></div>

**What to notice:** the observed F is only just into the tail (p = 0.068). The biggest pairwise gap, DASH + low-sodium vs control, barely clears the Tukey bar (adjusted p = 0.047), and its CI nearly touches zero.

## When to use it

Use Tukey's HSD when:

- You have **three or more independent groups** and a [one-way ANOVA](#/test/one-way-anova) as the framework
- You want **all pairwise comparisons**, with no single control group
- Group variances are roughly equal

**Clinical scenario.** A trial randomises 24 patients with hypertension to four diets (6 per group): control, low-sodium, DASH, and DASH plus low-sodium. The outcome is **systolic BP (mmHg) at week 8**. *Which diets differ from which, and by how much?*

**Don't use it when:**

- You only compare each arm with **one control** → use [Dunnett's test](#/test/dunnett), which has more power for that question
- Variances are clearly unequal → use a Games–Howell test, or Welch t-tests with a [multiplicity correction](#/test/multiple-comparisons)
- The data are skewed and groups are small → use the [Kruskal–Wallis test](#/test/kruskal-wallis) with a rank-based post-hoc test
- You have a small, pre-specified set of contrasts → use [Bonferroni, Holm, or BH corrections](#/test/multiple-comparisons)

## Assumptions

1. **Independent observations**, as in [one-way ANOVA](#/test/one-way-anova).
2. **Residuals are approximately normal.** Check a Q–Q plot or [Shapiro–Wilk](#/test/shapiro-wilk) on the residuals.
3. **Equal variances across groups.** Check group SDs or [Levene's test](#/test/levene). Tukey is somewhat robust with balanced groups.
4. **The pairwise comparisons are of interest.** If only some are, the correction is wasteful.

## Hypotheses

For each pair of groups i and j:

- **H₀:** μᵢ = μⱼ
- **H₁:** μᵢ ≠ μⱼ

All pairs are tested together at a **family-wise** α = 0.05, so the chance of any false positive across the whole set is at most 5%.

## Example data

Systolic BP (mmHg) at week 8, 6 patients per diet:

| Diet | Values | Mean |
|---|---|---|
| Control | 135, 133, 158, 135, 130, 144 | 139.17 |
| Low-sodium | 125, 120, 142, 124, 140, 128 | 129.83 |
| DASH | 132, 125, 139, 119, 139, 123 | 129.50 |
| DASH + low-sodium | 122, 128, 131, 125, 109, 131 | 124.33 |

## How to run it

### R

```r
# Systolic BP (mmHg) under four diets, n = 6 per group
bp <- data.frame(
  sbp = c(135,133,158,135,130,144,    # Control
          125,120,142,124,140,128,    # Low-sodium
          132,125,139,119,139,123,    # DASH
          122,128,131,125,109,131),   # DASH + low-sodium
  diet = factor(rep(c("Control","LowSodium","DASH","DASH_LowSodium"), each = 6),
                levels = c("Control","LowSodium","DASH","DASH_LowSodium"))
)

# Step 1: omnibus ANOVA
fit <- aov(sbp ~ diet, data = bp)
summary(fit)

# Step 2: Tukey HSD on all 6 pairwise comparisons
TukeyHSD(fit, conf.level = 0.95)
```

**What the code does**

- `TukeyHSD()` needs an object from **`aov()`**, not `lm()`. For `lm` objects, use `emmeans` or `multcomp`.
- `diet` must be a factor. The order of levels sets the direction of each difference (later minus earlier, e.g., `LowSodium-Control`).
- `conf.level = 0.95` gives **simultaneous** CIs: they are wider than the usual per-comparison CIs.
- The `p adj` column is already adjusted for all 6 comparisons; do not adjust it again.

**Output**

```text
            Df Sum Sq Mean Sq F value Pr(>F)  
diet         3  686.5  228.82   2.779 0.0677 .
Residuals   20 1646.5   82.33                 

  Tukey multiple comparisons of means
    95% family-wise confidence level

$diet
                                diff       lwr        upr     p adj
LowSodium-Control         -9.3333333 -23.99550  5.3288361 0.3107723
DASH-Control              -9.6666667 -24.32884  4.9955028 0.2824595
DASH_LowSodium-Control   -14.8333333 -29.49550 -0.1711639 0.0467621
DASH-LowSodium            -0.3333333 -14.99550 14.3288361 0.9999042
DASH_LowSodium-LowSodium  -5.5000000 -20.16217  9.1621694 0.7227769
DASH_LowSodium-DASH       -5.1666667 -19.82884  9.4955028 0.7587532
```

### Python

```python
import numpy as np
from scipy import stats

# Systolic BP (mmHg) under four diets, n = 6 per group
control = [135,133,158,135,130,144]
lowsod  = [125,120,142,124,140,128]
dash    = [132,125,139,119,139,123]
dashlow = [122,128,131,125,109,131]

# Step 1: omnibus ANOVA
f, p = stats.f_oneway(control, lowsod, dash, dashlow)
print(f"ANOVA: F(3, 20) = {f:.3f}, p = {p:.4f}")

# Step 2: Tukey HSD (SciPy >= 1.8); pairwise differences are group i - group j
res = stats.tukey_hsd(control, lowsod, dash, dashlow)
ci = res.confidence_interval(0.95)
names = ["Control", "LowSodium", "DASH", "DASH_LowSodium"]
for i in range(4):
    for j in range(i + 1, 4):
        # statistic[i][j] = mean_i - mean_j, so flip signs to report j - i
        print(f"{names[j]} - {names[i]}: diff = {-res.statistic[i][j]:.2f}, "
              f"95% CI ({-ci.high[i][j]:.2f}, {-ci.low[i][j]:.2f}), adj p = {res.pvalue[i][j]:.4f}")
```

**What the code does**

- `stats.tukey_hsd` takes each group as an argument and handles unequal n (Tukey–Kramer) automatically.
- **Sign convention differs from R.** SciPy reports `mean_i − mean_j` in its matrices (earlier minus later), while R reports later minus earlier. The loop flips the signs so the output matches R.
- `res.confidence_interval()` returns matrices `.low` and `.high`; `res.pvalue` holds the adjusted p-values.
- `statsmodels.stats.multicomp.pairwise_tukeyhsd` is an alternative that works on long-format data.

**Output**

```text
ANOVA: F(3, 20) = 2.779, p = 0.0677
LowSodium - Control: diff = -9.33, 95% CI (-24.00, 5.33), adj p = 0.3108
DASH - Control: diff = -9.67, 95% CI (-24.33, 5.00), adj p = 0.2825
DASH_LowSodium - Control: diff = -14.83, 95% CI (-29.50, -0.17), adj p = 0.0468
DASH - LowSodium: diff = -0.33, 95% CI (-15.00, 14.33), adj p = 0.9999
DASH_LowSodium - LowSodium: diff = -5.50, 95% CI (-20.16, 9.16), adj p = 0.7228
DASH_LowSodium - DASH: diff = -5.17, 95% CI (-19.83, 9.50), adj p = 0.7588
```

### SAS

```sas
data diets;
  input diet $ sbp @@;
  datalines;
Control 135 Control 133 Control 158 Control 135 Control 130 Control 144
LowSod 125 LowSod 120 LowSod 142 LowSod 124 LowSod 140 LowSod 128
DASH 132 DASH 125 DASH 139 DASH 119 DASH 139 DASH 123
DASHLow 122 DASHLow 128 DASHLow 131 DASHLow 125 DASHLow 109 DASHLow 131
;
run;

/* ANOVA plus Tukey-Kramer comparisons with confidence limits */
proc glm data=diets;
  class diet;
  model sbp = diet;
  means diet / tukey cldiff alpha=0.05;
run;
quit;
```

**What the code does**

- `class diet;` and `model sbp = diet;` fit the one-way ANOVA, whose F test is printed first.
- `means diet / tukey cldiff` requests Tukey's HSD (Tukey–Kramer when n is unequal) and prints **confidence limits** for each pairwise difference. Without `cldiff`, SAS prints only letter groupings.
- `alpha=0.05` sets the family-wise level (the default).
- For **adjusted p-values**, use `lsmeans diet / pdiff=all adjust=tukey cl;` instead. Note that `LSMEANS` labels differences by order of the levels.
- SAS orders `diet` levels alphabetically (Control, DASH, DASHLow, LowSod), so the sign and order of the differences will not match R's order unless you reorder the levels (e.g., with a numeric coding and a format).

## Interpreting the output

- **Omnibus ANOVA: F(3, 20) = 2.78, p = 0.068.** The overall test is not significant at 5%. If your analysis plan requires a significant omnibus test first (a "protected" approach), you would stop here. Tukey's procedure does not need that gate, but note that the two results can disagree, as they do here.
- **DASH + low-sodium vs control: −14.8 mmHg, 95% CI −29.5 to −0.2, adjusted p = 0.047.** This is the only pair whose simultaneous CI excludes 0. The CI is wide, and its upper limit is almost 0, so the result is fragile.
- **Other pairs:** low-sodium vs control (−9.3 mmHg, adj p = 0.31) and DASH vs control (−9.7 mmHg, adj p = 0.28) have CIs spanning 0. The differences among the three intervention diets are small (−0.3 to −5.5 mmHg) with adjusted p ≥ 0.72.
- **Why the adjusted p-values are large.** With 6 comparisons, each CI is about ±14.7 mmHg here, so only large differences stand out in a study of 6 per group.
- **R and Python agree** for every difference, CI, and adjusted p.

**How to report it:**
> After Tukey adjustment for six pairwise comparisons, systolic BP was 14.8 mmHg lower with DASH plus low-sodium than with the control diet (95% CI −29.5 to −0.2; adjusted p = 0.047); no other pairwise difference was statistically significant.

## Check your understanding

<details class="check"><summary>1. A trial compares four arms and you run all six pairwise t-tests at 0.05 without adjustment. Roughly what is the chance of at least one false positive if all true means are equal?</summary>

About 1 − 0.95⁶ ≈ 26% if the tests were independent (they are not exactly, but it is close). Tukey's HSD brings the family-wise rate back to 5%.
</details>

<details class="check"><summary>2. Three of four arms are active drugs and one is placebo, and you only care about each drug vs placebo. Is Tukey the right tool?</summary>

Not ideal. You only need 3 comparisons, and [Dunnett's test](#/test/dunnett) is designed for many-vs-one, spending less of the error budget and giving more power.
</details>

<details class="check"><summary>3. Low-sodium vs control has adjusted p = 0.31 and a CI of −24.0 to +5.3 mmHg. Does low-sodium diet not work?</summary>

Not shown. The CI includes reductions of up to 24 mmHg, which would be clinically important. A non-significant pair in a study of 6 per group mainly says the study is small.
</details>

<details class="check"><summary>4. What would happen to the width of the Tukey CIs if we compared only two diets rather than four?</summary>

They would narrow. With two groups there is a single comparison and Tukey reduces to an ordinary (pooled) t-test CI. The more groups you compare, the bigger the multiplicity penalty.
</details>

## Common mistakes

- **Reading unadjusted and adjusted p-values together.** Report the Tukey-adjusted p-values and simultaneous CIs only. Do not apply Bonferroni on top.
- **Treating a non-significant pair as "no difference".** The CIs are wide; absence of significance is not equivalence.
- **Using Tukey with strongly unequal variances.** Switch to Games–Howell or Welch tests with a correction.
- **Sign confusion across software.** R reports later − earlier; SciPy's matrices report row − column. Always label the direction of each difference.
- **Using Tukey when only comparing with a control.** [Dunnett's test](#/test/dunnett) is more powerful for that structure.
- **Data dredging.** Pairwise testing should be pre-specified in the analysis plan (ICH E9), not chosen after seeing the means.

## Related tests

- [One-way ANOVA](#/test/one-way-anova): the omnibus test that usually accompanies Tukey (not required for it)
- [Dunnett's test](#/test/dunnett): each group vs a single control
- [Multiple-comparison corrections](#/test/multiple-comparisons): Bonferroni, Holm, and BH
- [Kruskal–Wallis test](#/test/kruskal-wallis): nonparametric alternative
- [Levene's test](#/test/levene): check equal variances

## Key takeaways

- Tukey's HSD compares **all pairs** while holding the family-wise error rate at 5%.
- It uses the studentized range distribution, so its CIs are wider than unadjusted ones.
- Report the adjusted p-value and simultaneous CI; never adjust again.
- It assumes normal residuals and roughly equal variances; otherwise use Games–Howell.
- For many-vs-control questions use Dunnett; a non-significant pair is not evidence of no difference.

## References

- Tukey JW. Comparing individual means in the analysis of variance. *Biometrics.* 1949;5:99–114.
- Kramer CY. Extension of multiple range tests to group means with unequal numbers of replications. *Biometrics.* 1956;12:307–310.
- R documentation: `?TukeyHSD`; SciPy documentation: `scipy.stats.tukey_hsd`
- SAS/STAT User's Guide: The GLM Procedure — MEANS statement (TUKEY, CLDIFF)
