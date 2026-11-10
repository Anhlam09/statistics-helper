# Independent two-sample t-test (Student & Welch)

> Compares the means of two separate groups of patients; the Student version assumes equal variances, the Welch version does not and is the recommended default.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Compare the means of two separate groups of different people | Continuous outcome | Independent observations; roughly normal data (or large groups) | Drug vs placebo change in SBP in a parallel-group trial |

## Learning objectives

- Recognize two-group designs with **different** people in each group
- Explain why Welch's t-test is the recommended default over Student's
- Run it in R and Python and read the mean difference, 95% CI, t, df, and p
- Report the result and avoid the classic pitfalls (software defaults, variance pre-tests)

## What it is

The independent t-test asks whether the difference between two group means is larger than chance would produce. It divides the observed difference by its standard error. The two versions differ only in how that standard error is built.

- **Student's t-test** pools the two variances into one estimate. It is exact when the groups really have equal variances and normal data.
- **Welch's t-test** keeps the two variances separate and adjusts the degrees of freedom (Satterthwaite). It stays valid when variances or group sizes differ, and loses almost nothing when they are equal. **Use Welch by default**; this avoids having to pre-test variances.
- **Test statistic (Welch):** t = (x̄₁ − x̄₂) / √(s₁²/n₁ + s₂²/n₂), with fractional df
- **Test statistic (Student):** t = (x̄₁ − x̄₂) / (s_p √(1/n₁ + 1/n₂)), df = n₁ + n₂ − 2
- **Effect size:** the mean difference with its 95% CI (primary); Cohen's d or Hedges' g (secondary)

## The intuition

The t statistic is a **signal-to-noise ratio**. The signal is the gap between the two group means. The noise is the standard error of that gap: how much it would wobble by chance if you repeated the trial. A big gap relative to the wobble is hard to blame on luck.

Think of two clinics whose waiting times you compare. If one averages 12 minutes less than the other, and waits in each clinic vary by only a few minutes, the gap is convincing. If waits in each clinic range from 5 to 90 minutes, the same 12-minute gap could easily be chance. Spread matters as much as the gap.

The plot below shows exactly this. The drug group's mean (−15.75) is far below placebo (−3.90), but the drug group is also much more spread out, from +5 to −40 mmHg. The placebo group is tightly bunched. Welch's test builds the noise term from each group's own spread, so it does not pretend the two groups are equally variable.

<div class="viz" data-viz="dots" data-unit="mmHg, change in SBP" data-groups='{"Drug":[-16,5,-21,-18,-34,-14,-3,-24,-3,-18,-40,-3],"Placebo":[-1,-4,-1,-6,-9,-3,-4,-2,-6,-3]}' data-zero="1" data-caption="Drug: lower on average but widely spread. Placebo: tightly clustered near a small reduction. Note that some drug patients overlap with placebo."></div>

Our Welch result is t = −3.02 with 11.94 degrees of freedom. Under H₀ (no difference in means), t values follow the curve below. Ours lands in the tail, so the shaded two-sided p-value area is small (about 0.011).

<div class="viz" data-viz="t" data-stat="-3.02" data-df="11.94"></div>

Notice the df: 11.94, not 20. Welch lowers the degrees of freedom when variances differ, which makes the tails a little heavier and the test more honest.

## When to use it

Use the independent t-test when **all** of these are true:

- The outcome is **continuous** (change in blood pressure, LDL-C, length of stay in days)
- There are **exactly two groups** of **different** people (e.g., drug vs placebo in a parallel-group trial)
- Observations are independent
- Data are approximately normal in each group, or each group is reasonably large

**Clinical scenario.** In a parallel-group trial, 12 patients receive a new antihypertensive and 10 receive placebo. The outcome is the **change in systolic BP (mmHg) from baseline to week 12**. The drug group is much more variable than the placebo group. *How much more does the drug lower SBP, and how precise is that estimate?*

**Don't use it when:**

- The same patients are measured twice → use the [paired t-test](#/test/paired-t)
- There are **three or more** groups → use [one-way ANOVA](#/test/one-way-anova)
- The data are skewed or have outliers and groups are small → use the [Mann–Whitney U test](#/test/mann-whitney)
- You need to adjust for baseline → use [ANCOVA](#/test/ancova)
- The outcome is binary → use the [two-proportion z-test](#/test/two-proportion-z) or [Fisher's exact test](#/test/fisher-exact)

## Assumptions

1. **The two groups are independent**, and observations within each group are independent. Check the design: no patient appears in both groups.
2. **Each group is approximately normal.** Check with dot plots or Q–Q plots; [Shapiro–Wilk](#/test/shapiro-wilk) can support this but is underpowered in small groups and over-sensitive in large ones, so do not use it to choose between tests. Robust for n ≳ 30 per group unless the data are heavily skewed.
3. **Equal variances (Student only).** Compare SDs, or use [Levene's test](#/test/levene). Welch does not need this assumption.
4. **Continuous outcome** on an interval or ratio scale.

## Hypotheses

- **H₀:** μ₁ = μ₂ (the difference in population means is zero)
- **H₁:** μ₁ ≠ μ₂ (two-sided)

## Example data

Change in SBP (mmHg) from baseline to week 12 (negative = reduction):

| Group | n | Values | Mean | SD |
|---|---|---|---|---|
| Drug | 12 | −16, 5, −21, −18, −34, −14, −3, −24, −3, −18, −40, −3 | −15.75 | 13.31 |
| Placebo | 10 | −1, −4, −1, −6, −9, −3, −4, −2, −6, −3 | −3.90 | 2.51 |

## How to run it

### R

```r
# Change in SBP (mmHg) from baseline to week 12
drug    <- c(-16, 5, -21, -18, -34, -14, -3, -24, -3, -18, -40, -3)   # n = 12
placebo <- c(-1, -4, -1, -6, -9, -3, -4, -2, -6, -3)                   # n = 10

# Compare variances first (descriptive)
sd(drug); sd(placebo)

# Welch (R's default): var.equal = FALSE
t.test(drug, placebo)

# Student (pooled variance)
t.test(drug, placebo, var.equal = TRUE)
```

**What the code does**

- `t.test(drug, placebo)` runs **Welch** by default, because `var.equal = FALSE` is the default. This is the recommended analysis.
- `var.equal = TRUE` switches to Student's pooled test.
- The reported difference is first group − second group (drug − placebo). The CI is for that difference.
- The two groups are different lengths, which is fine for an independent test (it would be an error for a paired test).

**Output**

```text
[1] 13.3084
[1] 2.514403

	Welch Two Sample t-test

data:  drug and placebo
t = -3.0205, df = 11.936, p-value = 0.01071
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -20.403075  -3.296925
sample estimates:
mean of x mean of y 
   -15.75     -3.90 

	Two Sample t-test

data:  drug and placebo
t = -2.764, df = 20, p-value = 0.01197
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -20.793054  -2.906946
sample estimates:
mean of x mean of y 
   -15.75     -3.90 
```

### Python

```python
import numpy as np
from scipy import stats

# Change in SBP (mmHg) from baseline to week 12
drug    = np.array([-16, 5, -21, -18, -34, -14, -3, -24, -3, -18, -40, -3])  # n = 12
placebo = np.array([-1, -4, -1, -6, -9, -3, -4, -2, -6, -3])                  # n = 10

print(f"SD drug = {drug.std(ddof=1):.2f}, SD placebo = {placebo.std(ddof=1):.2f}")
diff = drug.mean() - placebo.mean()

# Welch: equal_var=False (SciPy's default is Student, i.e. True!)
w = stats.ttest_ind(drug, placebo, equal_var=False)
wci = w.confidence_interval(0.95)
print(f"Welch:   diff = {diff:.2f}, t = {w.statistic:.3f}, df = {w.df:.2f}, p = {w.pvalue:.4f}, "
      f"95% CI ({wci.low:.2f}, {wci.high:.2f})")

# Student: pooled variance
s = stats.ttest_ind(drug, placebo, equal_var=True)
sci = s.confidence_interval(0.95)
print(f"Student: diff = {diff:.2f}, t = {s.statistic:.3f}, df = {s.df:.0f}, p = {s.pvalue:.4f}, "
      f"95% CI ({sci.low:.2f}, {sci.high:.2f})")
```

**What the code does**

- **Defaults differ.** SciPy's `ttest_ind` defaults to `equal_var=True` (Student), whereas R's `t.test` defaults to Welch. Always set `equal_var=False` explicitly to get Welch.
- `res.df` and `res.confidence_interval()` need SciPy ≥ 1.10. The CI is for the mean difference (first − second).
- The mean difference is computed by hand because the result object does not store it.

**Output**

```text
SD drug = 13.31, SD placebo = 2.51
Welch:   diff = -11.85, t = -3.020, df = 11.94, p = 0.0107, 95% CI (-20.40, -3.30)
Student: diff = -11.85, t = -2.764, df = 20, p = 0.0120, 95% CI (-20.79, -2.91)
```

### SAS

```sas
data bp;
  input group $ change @@;
  datalines;
Drug -16 Drug 5 Drug -21 Drug -18 Drug -34 Drug -14
Drug -3 Drug -24 Drug -3 Drug -18 Drug -40 Drug -3
Plac -1 Plac -4 Plac -1 Plac -6 Plac -9 Plac -3
Plac -4 Plac -2 Plac -6 Plac -3
;
run;

/* Reports Pooled (Student) AND Satterthwaite (Welch) in one run */
proc ttest data=bp alpha=0.05 sides=2;
  class group;
  var change;
run;
```

**What the code does**

- The data are in **long format**: one row per patient with a group label and the outcome.
- `class group;` names the two-level grouping variable. SAS orders levels alphabetically, so the difference is **Drug − Plac**, as in R.
- PROC TTEST prints both tests in the "Equality of Variances" and "T-Tests" tables: the **Pooled** row is Student and the **Satterthwaite** row is Welch. Read the Satterthwaite row.
- The "Folded F" equality-of-variances test is shown automatically. Use it as supporting information, not as a gate for choosing the test.
- The "Difference (1-2)" CI table gives the CI for both methods.

## Interpreting the output

- **Mean difference = −11.85 mmHg** (drug −15.75 vs placebo −3.90). The drug lowered SBP by about 12 mmHg more than placebo.
- **Welch 95% CI: −20.40 to −3.30 mmHg.** The data are consistent with a true extra reduction between roughly 3 and 20 mmHg. The CI excludes 0.
- **Welch t(11.94) = −3.02, p = 0.011.** The fractional df (11.94, not 20) reflects the large SD difference: the drug SD (13.3) is more than five times the placebo SD (2.5).
- **Student: t(20) = −2.76, p = 0.012, CI −20.79 to −2.91.** Here Student's and Welch's conclusions agree, but the Student CI is slightly wider and its p-value larger. With unequal group sizes and unequal variances, Student's error rate can be wrong in either direction; here the larger group has the larger variance, which makes Student *conservative*.
- **R and Python agree** for both tests.

**How to report it:**
> Systolic BP fell by 15.8 mmHg (SD 13.3) with the drug and 3.9 mmHg (SD 2.5) with placebo. The mean difference was −11.9 mmHg (95% CI −20.4 to −3.3; Welch t(11.9) = −3.02, p = 0.011).

## Check your understanding

<details class="check"><summary>1. A trial measures HbA1c in 40 patients on drug A and 40 different patients on drug B. Which t-test, and which version?</summary>

The independent two-sample t-test, using **Welch** by default. The groups contain different people, so the paired test does not apply.
</details>

<details class="check"><summary>2. A colleague runs Levene's test, gets p = 0.30, and switches to Student's t-test "because the variances are equal." Is this good practice?</summary>

Not really. Choosing the test based on a variance pre-test distorts error rates, and p = 0.30 does not prove variances are equal, especially with small groups. Welch loses almost nothing when variances are equal, so use it from the start and pre-specify that in the analysis plan.
</details>

<details class="check"><summary>3. The Shapiro–Wilk test is non-significant in each of two groups of n = 8. Does that confirm the t-test is valid?</summary>

No. With n = 8 the test has little power, so it would miss moderate skewness. Look at dot plots or Q–Q plots, think about whether the outcome is naturally skewed (lab values often are), and consider a transformation or the [Mann–Whitney U test](#/test/mann-whitney) as pre-specified sensitivity analysis. Deciding between tests from a normality p-value is itself discouraged.
</details>

<details class="check"><summary>4. What if the placebo group had the same large spread as the drug group (SD about 13 instead of 2.5), with the same means? What happens to t?</summary>

The signal (−11.85 mmHg) is unchanged but the standard error grows: roughly √(13.3²/12 + 13²/10) ≈ 5.6 versus 3.9 now. So |t| falls to about 2.1, the p-value rises, and the CI widens. More noise means a weaker signal-to-noise ratio.
</details>

## Common mistakes

- **Using Student by default.** If variances differ and group sizes are unequal, Student's p-values and CIs can be misleading. Welch costs almost nothing when variances are equal.
- **Pre-testing variances to choose.** Choosing between the tests based on a Levene or F test distorts error rates. Pre-specify Welch.
- **Software defaults.** R defaults to Welch; SciPy and many other tools default to Student.
- **Using it on paired data.** Two measurements per patient need the [paired t-test](#/test/paired-t).
- **Comparing change scores without considering baseline.** In randomised trials, [ANCOVA](#/test/ancova) is usually more efficient.
- **Reporting only the p-value.** Always give both group means and the difference with its CI.

## Related tests

- [Paired t-test](#/test/paired-t): same patients measured twice
- [Mann–Whitney U test](#/test/mann-whitney): nonparametric alternative
- [One-way ANOVA](#/test/one-way-anova): three or more groups
- [Levene's test](#/test/levene): compare variances
- [ANCOVA](#/test/ancova): adjust for baseline
- [Equivalence (TOST) tests](#/test/tost-equivalence): when the goal is to show similarity

## Key takeaways

- The independent t-test compares two means from **different** people: t = mean difference / standard error.
- **Use Welch by default.** It is safe whether or not variances are equal, and needs no variance pre-test.
- R defaults to Welch; SciPy defaults to Student. Set `equal_var=False`.
- Report both group means, the mean difference, and its 95% CI, not just p.
- If the same people appear twice, use the paired t-test instead.

## References

- Welch BL. The generalization of "Student's" problem when several different population variances are involved. *Biometrika.* 1947;34:28–35.
- Delacre M, Lakens D, Leys C. Why psychologists should by default use Welch's t-test instead of Student's t-test. *Int Rev Soc Psychol.* 2017;30:92–101.
- R documentation: `?t.test`; SciPy documentation: `scipy.stats.ttest_ind`
- SAS/STAT User's Guide: The TTEST Procedure
