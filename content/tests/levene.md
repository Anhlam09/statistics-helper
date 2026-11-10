# Levene's test

> Tests whether two or more groups have the same variance (spread), by running an ANOVA on how far each observation lies from its group's centre.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Check whether two or more independent groups have equal spread | Continuous outcome, 2+ groups | Independent observations; continuous, not heavily tied | LDL-C change at three clinics before pooling |

## Learning objectives

- Explain what Levene's test compares and why equal variances matter
- Distinguish the mean-centred (Levene) and median-centred (Brown–Forsythe) versions
- Run it in R and Python and read W, df, and p alongside the group SDs
- Recognize why a non-significant result is not proof of equal variances, and why pre-testing is not a good gatekeeper

## What it is

Many tests, including the classic Student t-test and one-way ANOVA, assume that the groups are equally variable. Levene's test checks that assumption. It replaces each value with its **absolute deviation from the group centre** and asks, with a one-way ANOVA, whether those deviations differ between groups. If one group is far more spread out, its deviations are larger and the test is significant.

- **Test statistic:** W (an F statistic) with k − 1 and N − k degrees of freedom, for k groups and N observations
- **Variants:** centring on the **mean** (original Levene) or on the **median** (Brown–Forsythe), which is more robust to skewed data. Median centring is the default in R's `car::leveneTest` and in SciPy.
- **Effect size:** the ratio of group SDs (largest ÷ smallest) is the most useful descriptive measure; report it alongside the p-value.

## The intuition

Think of three clinics measuring the same thing. If one clinic's patients scatter widely around their clinic average while another's cluster tightly, the clinics differ in **spread**, not just in average. Levene's trick is to turn a question about spread into a question about averages. Replace each value with its **distance from its own group's centre** (always positive). A widely scattered group now has large distances, a tight group has small ones. Then run an ordinary one-way ANOVA on those distances: if the average distance differs between groups, the spreads differ.

Here the signal is the difference in average distance between clinics, and the noise is how much distances vary within each clinic. Look at the raw data first. Clinic C reaches from about −22 to +3, while clinic A stays between −18 and −7.

<div class="viz" data-viz="dots" data-unit="mg/dL, LDL-C change" data-groups='{"A":[-14,-12,-18,-17,-7,-16,-7,-10],"B":[-10,-19,-17,-13,-24,-12,-20,-10],"C":[-15,3,-20,-21,-22,-11,-18,-5]}' data-caption="Clinic C is visibly wider than A (SD 8.8 vs 4.3), but with 8 patients per clinic the eye can easily be fooled. Note the single +3 value in C."></div>

The Brown–Forsythe statistic is W = 1.27 on 2 and 21 degrees of freedom. If all three variances were truly equal, W values would follow the F curve below. Ours sits near the middle of the bulk of the curve, so the shaded upper-tail area (the p-value) is large, about 0.30.

<div class="viz" data-viz="f" data-stat="1.27" data-df1="2" data-df2="21"></div>

Notice the tension: the SD ratio is about 2, yet p is not small. That is a statement about **low power with small groups**, not evidence that the spreads are equal.

## When to use it

Use Levene's test when:

- The outcome is **continuous** and you have **two or more independent groups**
- You want to check the **equal-variance assumption** of Student's t-test or one-way ANOVA
- The data may not be perfectly normal (Levene is far less sensitive to non-normality than Bartlett's test)

**Clinical scenario.** A multi-centre study measures the change in LDL-C (mg/dL) after 12 weeks of treatment at three clinics, with 8 patients per clinic. Before pooling the clinics in a one-way ANOVA, the analyst asks: *is the variability in response similar at all three sites?*

**Don't use it when:**

- You want to **decide between the Student and Welch t-test** or between ANOVA and Welch ANOVA. Pre-testing and then choosing the analysis changes the error rate of the final test, and Levene has little power exactly when n is small and the choice matters. Use Welch by default and skip the pre-test (see [independent t-test](#/test/independent-t)). Levene is better used descriptively, or to *report* heterogeneity.
- The groups are **paired or repeated** measures on the same patients → compare variances with a method that accounts for correlation (e.g., a [mixed model](#/test/mixed-models) with a variance structure)
- The question is about **normality** → use the [Shapiro–Wilk test](#/test/shapiro-wilk)

## Assumptions

1. **Observations are independent**, within and between groups. Check the study design.
2. **The outcome is continuous.** Very coarse or heavily tied data weaken the test.
3. **Groups are reasonably sized.** With very small n per group the test has little power, so a non-significant result is weak reassurance. Look at the group SDs too.

## Hypotheses

- **H₀:** σ₁² = σ₂² = … = σ_k² (all group variances are equal)
- **H₁:** at least one group variance differs

## Example data

Change in LDL-C (mg/dL), 8 patients per clinic (negative = reduction):

| Clinic | Values | SD |
|---|---|---|
| A | −14, −12, −18, −17, −7, −16, −7, −10 | 4.34 |
| B | −10, −19, −17, −13, −24, −12, −20, −10 | 5.15 |
| C | −15, 3, −20, −21, −22, −11, −18, −5 | 8.81 |

## How to run it

### R

```r
# Change in LDL-C (mg/dL) at three clinics, n = 8 each
ldl <- data.frame(
  change = c(-14,-12,-18,-17,-7,-16,-7,-10,
             -10,-19,-17,-13,-24,-12,-20,-10,
             -15,3,-20,-21,-22,-11,-18,-5),
  clinic = factor(rep(c("A","B","C"), each = 8))
)

# SDs by clinic
tapply(ldl$change, ldl$clinic, sd)

# Brown-Forsythe version (center = median) -- the car default
library(car)
leveneTest(change ~ clinic, data = ldl)

# Classic Levene (center = mean)
leveneTest(change ~ clinic, data = ldl, center = mean)
```

**What the code does**

- `clinic` must be a **factor**. If it is numeric, `leveneTest` warns and converts it.
- `tapply(..., sd)` gives each group's SD, the descriptive counterpart of the test.
- `leveneTest(change ~ clinic)` uses `center = median` by default, which is the **Brown–Forsythe** test. Add `center = mean` to get the original Levene test.
- The output is an ANOVA-style table: `Df` is k − 1 and N − k, then the F value and p-value.

**Output**

```text
       A        B        C 
4.340425 5.153016 8.814558 

Levene's Test for Homogeneity of Variance (center = median)
      Df F value Pr(>F)
group  2  1.2693 0.3017
      21               

Levene's Test for Homogeneity of Variance (center = mean)
      Df F value Pr(>F)
group  2   2.407 0.1145
      21               
```

### Python

```python
import numpy as np
from scipy import stats

# Change in LDL-C (mg/dL) at three clinics, n = 8 each
A = [-14,-12,-18,-17,-7,-16,-7,-10]
B = [-10,-19,-17,-13,-24,-12,-20,-10]
C = [-15,3,-20,-21,-22,-11,-18,-5]

# SDs by clinic
print("SDs:", [f"{np.std(g, ddof=1):.2f}" for g in (A, B, C)])

# scipy's default is center='median' (Brown-Forsythe)
bf = stats.levene(A, B, C)
print(f"Brown-Forsythe: W = {bf.statistic:.4f}, p = {bf.pvalue:.4f}")

# Classic Levene (center = mean)
lv = stats.levene(A, B, C, center='mean')
print(f"Levene (mean):  W = {lv.statistic:.4f}, p = {lv.pvalue:.4f}")
```

**What the code does**

- `stats.levene(A, B, C)` takes each group as a separate argument. The default is `center='median'`, the same as R's `car`.
- `center='mean'` switches to the original Levene test. SciPy also offers `'trimmed'` for heavy-tailed data.
- `np.std(g, ddof=1)` is the sample SD. NumPy's default `ddof=0` would understate it.

**Output**

```text
SDs: ['4.34', '5.15', '8.81']
Brown-Forsythe: W = 1.2693, p = 0.3017
Levene (mean):  W = 2.4070, p = 0.1145
```

### SAS

```sas
data ldl;
  input clinic $ change @@;
  datalines;
A -14 A -12 A -18 A -17 A -7 A -16 A -7 A -10
B -10 B -19 B -17 B -13 B -24 B -12 B -20 B -10
C -15 C 3 C -20 C -21 C -22 C -11 C -18 C -5
;
run;

/* HOVTEST=BF is the Brown-Forsythe (median-centred) test */
proc glm data=ldl;
  class clinic;
  model change = clinic;
  means clinic / hovtest=bf;
run;
quit;
```

**What the code does**

- `input clinic $ change @@;` reads a clinic label and a value as a pair, repeatedly along each line.
- `hovtest=bf` in the `MEANS` statement requests the Brown–Forsythe test, which matches the R and Python defaults.
- `hovtest=levene` is also available, but SAS's default Levene uses **squared** deviations (`type=square`). To reproduce the mean-centred absolute-deviation version used by R and SciPy, write `hovtest=levene(type=abs)`.
- PROC GLM also prints the usual ANOVA table; `HOVTEST` adds a separate "Brown and Forsythe's Test for Homogeneity" table.

## Interpreting the output

- **Group SDs: 4.3, 5.2, and 8.8 mg/dL.** Clinic C is about twice as variable as clinic A (SD ratio ≈ 2.0). This is a noticeable difference in the raw data.
- **Brown–Forsythe: W(2, 21) = 1.27, p = 0.30.** No statistically significant evidence that variances differ.
- **Classic Levene (mean-centred): W = 2.41, p = 0.11.** Also not significant, but smaller, because the mean-centred version reacts more to the one extreme value (+3) in clinic C.
- **R and Python agree** exactly for both versions.
- **Caution.** With 8 patients per clinic the test has low power. A non-significant result does not show the variances are equal; the SD ratio of 2 is a reason to consider a Welch ANOVA or a variance-robust method as a sensitivity analysis.

**How to report it:**
> Variability in LDL-C change did not differ significantly between clinics (Brown–Forsythe W(2, 21) = 1.27, p = 0.30; SDs 4.3, 5.2, and 8.8 mg/dL).

## Check your understanding

<details class="check"><summary>1. A colleague runs Levene's test (p = 0.40) and then uses Student's t-test "because variances are equal". What is wrong with this workflow?</summary>

Two problems. First, a non-significant Levene test is not evidence of equal variances, especially with small groups. Second, choosing the analysis based on a pre-test distorts the type I error of the final test. Welch's test loses almost nothing when variances are equal, so use it by default.
</details>

<details class="check"><summary>2. In our example the Brown–Forsythe p is 0.30 but the mean-centred Levene p is 0.11. Why do they differ?</summary>

The mean-centred version measures distances from the group mean, and one extreme value (+3 in clinic C) drags that mean and creates large distances. The median-centred version is more robust to such outliers and skewness. Neither is significant here, but the gap is a hint that one point is influential. Always state which version you ran.
</details>

<details class="check"><summary>3. Clinic SDs are 4.3, 5.2, and 8.8 mg/dL and p = 0.30. A reviewer says "variances are homogeneous". Do you agree?</summary>

Not necessarily. The SD ratio of about 2 is large and n = 8 per clinic gives the test little power. The data are compatible with equal variances but also with a real 2-fold difference. A Welch ANOVA, or a model with separate variances, is a safer primary or sensitivity analysis.
</details>

<details class="check"><summary>4. What if clinic C had 80 patients with the same SD instead of 8?</summary>

The test would have far more power, so the same SD ratio would likely give a small p-value. The statistic depends on both the size of the spread difference and the sample size, so the same descriptive SDs can lead to very different conclusions at different n.
</details>

## Common mistakes

- **Mixing up the variants.** R and SciPy default to the median; SAS's `hovtest=levene` defaults to squared deviations around the mean. Always state which version you used.
- **Using the p-value as a gatekeeper.** Pre-testing and then switching methods affects the overall error rate. Prefer a method that is robust to unequal variances by default.
- **Reading non-significance as equality.** With small groups, even a 2-fold SD ratio can be non-significant. Report the SDs.
- **Using Bartlett's test on skewed data.** Bartlett is very sensitive to non-normality; Levene/Brown–Forsythe is safer.
- **Forgetting that groups must be independent.** It is not valid for paired or repeated measurements.

## Related tests

- [Shapiro–Wilk test](#/test/shapiro-wilk): the companion check for normality
- [Independent two-sample t-test](#/test/independent-t): Student (equal variances) vs Welch (unequal)
- [One-way ANOVA](#/test/one-way-anova): the equal-variance assumption applies here too
- [Kruskal–Wallis test](#/test/kruskal-wallis): rank-based alternative when assumptions fail
- [Linear mixed-effects models](#/test/mixed-models): model unequal variances directly

## Key takeaways

- Levene's test runs a one-way ANOVA on each value's absolute distance from its group centre.
- Median centring (Brown–Forsythe) is the robust default in R and SciPy; mean centring reacts more to outliers.
- p > 0.05 does not show variances are equal. Look at the SD ratio and the group sizes.
- Do not use it as a gatekeeper between Student and Welch. Use Welch by default.

## References

- Levene H. Robust tests for equality of variances. In: Olkin I, ed. *Contributions to Probability and Statistics.* Stanford University Press; 1960:278–292.
- Brown MB, Forsythe AB. Robust tests for the equality of variances. *J Am Stat Assoc.* 1974;69:364–367.
- R documentation: `car::leveneTest`; SciPy documentation: `scipy.stats.levene`
- SAS/STAT User's Guide: The GLM Procedure — MEANS statement, HOVTEST option
