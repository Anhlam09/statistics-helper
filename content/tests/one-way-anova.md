# One-way ANOVA

> Tests whether the means of three or more independent groups are all equal, by comparing the variation between group means with the variation within groups.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Compare means of three or more independent groups | Continuous outcome, one categorical factor | Independent observations; residuals roughly normal; similar spreads | Pain reduction on placebo, low dose, or high dose |

## Learning objectives

- Explain the F statistic as between-group variance divided by within-group variance
- Run one-way ANOVA in R and Python and read the ANOVA table
- Report an effect size (η², ideally also ω²) together with group means
- Explain why a significant F needs a multiplicity-adjusted follow-up, and when to prefer Welch's ANOVA

## What it is

One-way ANOVA splits the total variability in the outcome into two parts: the variability **between** group means and the variability **within** groups (the noise). If the groups truly have the same mean, the between-group variability should be about the same size as the within-group variability. The ratio of the two is the F statistic.

- **Test statistic:** F = MS_between / MS_within, where MS = SS / df
- **Degrees of freedom:** k − 1 (between) and N − k (within), for k groups and N observations
- **Effect size:** group means with 95% CIs (primary); η² = SS_between / SS_total (the proportion of variance explained). η² is **biased upward**, especially with small samples, because it counts chance differences between group means as explained variance. **ω²** = (SS_between − (k − 1) · MS_within) / (SS_total + MS_within) is a less biased alternative; report it (or partial ω²) when n is small.
- **Important:** a significant F says *some* group differs, not *which* one. Follow up with [Tukey's HSD](#/test/tukey-hsd) or [Dunnett's test](#/test/dunnett).

## The intuition

ANOVA asks one question: *are the group averages further apart than I would expect from ordinary patient-to-patient scatter?* Picture three clinics' waiting-room crowds. If the average height differs by 1 cm between crowds but individuals within each crowd vary by 15 cm, the clinics look alike. If the crowds differ by 15 cm and individuals vary by only 2 cm, they clearly differ. ANOVA formalises that comparison.

- **Signal (between-group variance):** how far the group means sit from the overall mean.
- **Noise (within-group variance):** how far individual patients sit from *their own* group mean.
- **F = signal / noise.** If H₀ is true, F is about 1. A large F means the group means are spread far more than the within-group scatter can explain.

In the plot below, look at the horizontal marks for the three group means (3.8, 4.4, 6.1) and compare their spread with the vertical scatter of the dots around each mean. The means are separated, but the dots overlap a lot, so the signal is moderate relative to the noise.

<div class="viz" data-viz="dots" data-unit="VAS points, pain reduction" data-groups='{"Placebo":[5.4,3.3,1.9,6.7,3.0,0.3,4.0,5.7],"Low dose":[3.6,3.2,3.9,5.9,5.5,4.6,3.7,5.1],"High dose":[7.7,6.2,5.0,6.7,4.4,6.6,5.9,6.5]}' data-caption="Group means differ (3.8, 4.4, 6.1), yet within-group scatter is wide. F compares these two sources of variation."></div>

Here MS_between = 11.65 and MS_within = 2.17, so F = 5.36 with 2 and 21 degrees of freedom. Under H₀ the F curve is skewed right with its bulk near 1. Our F falls in the shaded upper tail, which is why p = 0.013.

<div class="viz" data-viz="f" data-stat="5.36" data-df1="2" data-df2="21"></div>

Notice that the F test is **one-sided by construction**: only large F counts as evidence, because both a larger signal and a smaller noise push F upward.

## When to use it

Use one-way ANOVA when **all** of these are true:

- The outcome is **continuous**
- There are **three or more independent groups** defined by one factor (dose level, diet, clinic)
- Observations are independent
- Residuals are approximately normal and group variances are similar

**Clinical scenario.** A dose-ranging study randomises 24 patients with chronic pain to placebo, low-dose, or high-dose analgesic (8 per arm). The outcome is the **reduction in pain VAS (0–10 points) after 4 weeks**. *Is there any difference in mean pain reduction across the three arms?*

**Don't use it when:**

- There are only **two** groups → use the [independent t-test](#/test/independent-t) (Welch)
- The same patients are measured at several times → use [repeated-measures ANOVA](#/test/repeated-measures-anova) or a [mixed model](#/test/mixed-models)
- There are **two factors** (e.g., drug and sex) → use [two-way ANOVA](#/test/two-way-anova)
- You need to adjust for baseline → use [ANCOVA](#/test/ancova)
- Data are strongly skewed or ordinal → use the [Kruskal–Wallis test](#/test/kruskal-wallis)

## Assumptions

1. **Independent observations.** Each patient contributes one value and belongs to one group.
2. **Residuals are approximately normal.** Check a Q–Q plot of the residuals or run [Shapiro–Wilk](#/test/shapiro-wilk) on them. The test is fairly robust with balanced groups and n ≥ 10–15 per group.
3. **Equal variances across groups (homoscedasticity).** Look at the group SDs (a ratio above about 2 is a warning sign). [Levene's test](#/test/levene) can be reported, but it should not decide the method: it has little power in small samples, and choosing the test after a pre-test distorts error rates. If in doubt, or if group sizes are unequal, use Welch's ANOVA (`oneway.test` in R) as the primary analysis.
4. **Continuous outcome.**

## Hypotheses

- **H₀:** μ₁ = μ₂ = … = μ_k (all group means are equal)
- **H₁:** at least one group mean differs from the others

## Example data

Reduction in pain VAS (points), 8 patients per arm:

| Arm | Values | Mean |
|---|---|---|
| Placebo | 5.4, 3.3, 1.9, 6.7, 3.0, 0.3, 4.0, 5.7 | 3.79 |
| Low dose | 3.6, 3.2, 3.9, 5.9, 5.5, 4.6, 3.7, 5.1 | 4.44 |
| High dose | 7.7, 6.2, 5.0, 6.7, 4.4, 6.6, 5.9, 6.5 | 6.12 |

## How to run it

### R

```r
# Pain reduction (VAS points, 0-10) after 4 weeks at three doses of an analgesic
pain <- data.frame(
  reduction = c(5.4,3.3,1.9,6.7,3.0,0.3,4.0,5.7,    # placebo
                3.6,3.2,3.9,5.9,5.5,4.6,3.7,5.1,    # low dose
                7.7,6.2,5.0,6.7,4.4,6.6,5.9,6.5),   # high dose
  group = factor(rep(c("Placebo","Low","High"), each = 8),
                 levels = c("Placebo","Low","High"))
)

# Fit the model and get the ANOVA table
fit <- aov(reduction ~ group, data = pain)
summary(fit)

# Effect size: eta-squared = SS_between / SS_total
ss <- summary(fit)[[1]][["Sum Sq"]]
ss[1] / sum(ss)

# Check assumptions on the residuals
shapiro.test(residuals(fit))
library(car)
leveneTest(reduction ~ group, data = pain)   # median-centred (Brown-Forsythe)
```

**What the code does**

- `group` must be a **factor**. If it were numeric (e.g., 0/1/2), `aov` would treat it as a continuous slope with 1 df, not a 3-group comparison.
- `levels = c("Placebo","Low","High")` sets the order of groups; the first level is the reference in later contrasts.
- `summary(fit)` prints the ANOVA table: Df, Sum Sq, Mean Sq, F value, Pr(>F).
- `ss[1] / sum(ss)` is η²: the between-group SS divided by the total SS.
- The assumption checks use the **residuals** (`residuals(fit)`) for normality and Levene's test for equal variances.

**Output**

```text
            Df Sum Sq Mean Sq F value Pr(>F)  
group        2  23.29  11.645   5.358 0.0132 *
Residuals   21  45.64   2.173                 
---
[1] 0.3378748

	Shapiro-Wilk normality test

data:  residuals(fit)
W = 0.98707, p-value = 0.9842

Levene's Test for Homogeneity of Variance (center = median)
      Df F value  Pr(>F)  
group  2  2.9869 0.07218 .
      21                  
```

### Python

```python
import numpy as np
from scipy import stats

# Pain reduction (VAS points, 0-10) after 4 weeks at three doses
placebo = [5.4,3.3,1.9,6.7,3.0,0.3,4.0,5.7]
low     = [3.6,3.2,3.9,5.9,5.5,4.6,3.7,5.1]
high    = [7.7,6.2,5.0,6.7,4.4,6.6,5.9,6.5]

# One-way ANOVA (classic, equal-variance F-test)
f, p = stats.f_oneway(placebo, low, high)
print(f"F(2, 21) = {f:.3f}, p = {p:.4f}")

# Eta-squared = SS_between / SS_total
allv = np.concatenate([placebo, low, high])
gm = allv.mean()
ss_b = sum(len(g) * (np.mean(g) - gm) ** 2 for g in (placebo, low, high))
ss_t = ((allv - gm) ** 2).sum()
print(f"SS_between = {ss_b:.2f}, SS_total = {ss_t:.2f}, eta^2 = {ss_b/ss_t:.3f}")

# Assumption checks
resid = np.concatenate([np.array(g) - np.mean(g) for g in (placebo, low, high)])
print(f"Shapiro-Wilk on residuals: p = {stats.shapiro(resid).pvalue:.4f}")
print(f"Levene (median): p = {stats.levene(placebo, low, high).pvalue:.4f}")
print("Means:", [f"{np.mean(g):.2f}" for g in (placebo, low, high)])
```

**What the code does**

- `stats.f_oneway` takes each group as a separate argument and runs the classic equal-variance F-test. It returns only F and p, not the full table, so η² is computed by hand.
- Residuals in a one-way design are each value minus its group mean.
- `stats.levene` uses the median-centred (Brown–Forsythe) form by default, matching R's `car::leveneTest`.
- `statsmodels` (`ols` + `anova_lm`) gives a full ANOVA table if you prefer.

**Output**

```text
F(2, 21) = 5.358, p = 0.0132
SS_between = 23.29, SS_total = 68.93, eta^2 = 0.338
Shapiro-Wilk on residuals: p = 0.9842
Levene (median): p = 0.0722
Means: ['3.79', '4.44', '6.12']
```

### SAS

```sas
data pain;
  input group $ reduction @@;
  datalines;
Placebo 5.4 Placebo 3.3 Placebo 1.9 Placebo 6.7
Placebo 3.0 Placebo 0.3 Placebo 4.0 Placebo 5.7
Low 3.6 Low 3.2 Low 3.9 Low 5.9 Low 5.5 Low 4.6 Low 3.7 Low 5.1
High 7.7 High 6.2 High 5.0 High 6.7 High 4.4 High 6.6 High 5.9 High 6.5
;
run;

/* One-way ANOVA with Brown-Forsythe test of equal variances */
proc glm data=pain plots=diagnostics;
  class group;
  model reduction = group;
  means group / hovtest=bf;
run;
quit;
```

**What the code does**

- `class group;` declares the categorical factor (without it, GLM would fit a slope).
- `model reduction = group;` fits the one-way ANOVA. The output has the model F test, the Type I/III SS table, and **R-Square, which equals η²** for a one-way design.
- `means group / hovtest=bf;` prints group means and SDs and the Brown–Forsythe test of equal variances.
- `plots=diagnostics` produces residual plots (Q–Q plot, residual-vs-predicted) for assumption checks.
- For a Welch ANOVA when variances differ, add `welch` to the `MEANS` statement: `means group / welch;`.

## Interpreting the output

- **F(2, 21) = 5.36, p = 0.013.** There is evidence that mean pain reduction is not the same in all three arms.
- **Group means: placebo 3.79, low dose 4.44, high dose 6.12 VAS points.** The pattern is dose-ordered, with high dose highest. The ANOVA does not say which pairs differ; follow up with [Tukey's HSD](#/test/tukey-hsd) for all pairs or [Dunnett's test](#/test/dunnett) for each dose vs placebo.
- **η² = 0.338.** About 34% of the variance in pain reduction is associated with treatment arm; a large effect by conventional benchmarks. η² overstates the population effect in small samples. The less biased ω², computed from the table above as (23.29 − 2 × 2.173) / (68.93 + 2.173), is about 0.27.
- **Residual normality: Shapiro–Wilk p = 0.98.** No evidence against normality.
- **Levene (median) p = 0.072.** Not significant at 5%, but that does not show the variances are equal; it is borderline with only 8 per arm. A sensitivity analysis with Welch's ANOVA (`oneway.test(reduction ~ group, data = pain)` in R) is reasonable.
- **R and Python agree** on F, p, and η².

**How to report it:**
> Mean pain reduction differed across arms (one-way ANOVA F(2, 21) = 5.36, p = 0.013, η² = 0.34): placebo 3.8, low dose 4.4, and high dose 6.1 VAS points.

## Check your understanding

<details class="check"><summary>1. A study compares HbA1c in three diets, with 40 patients per diet. A colleague proposes three t-tests (A vs B, A vs C, B vs C) at α = 0.05 each. What is the problem?</summary>

Three tests give roughly a 14% chance of at least one false positive when all diets are truly equal. Run one ANOVA, then a multiplicity-adjusted follow-up such as [Tukey's HSD](#/test/tukey-hsd), or pre-specify the comparisons you care about.
</details>

<details class="check"><summary>2. The ANOVA gives p = 0.013. Can you say the high dose is better than placebo?</summary>

Not from the F test alone. It says at least one mean differs. To claim high vs placebo you need a pairwise comparison with multiplicity control ([Tukey](#/test/tukey-hsd) or [Dunnett](#/test/dunnett)), reported with a CI.
</details>

<details class="check"><summary>3. Another trial gets F(2, 57) = 1.1, p = 0.34. The abstract says "the doses do not differ." Fair?</summary>

No. A non-significant F means the data did not show a difference; it does not show the means are equal. Report the group means with CIs and check whether clinically important differences are still compatible with the data.
</details>

<details class="check"><summary>4. What if the three group means stayed at 3.8, 4.4, and 6.1 but the within-group scatter were half as large?</summary>

MS_within would drop to about a quarter (SD halves, so variance falls to one quarter), while MS_between stays the same. F would be about four times bigger (roughly 21), and p would be far smaller. Same signal, less noise.
</details>

## Common mistakes

- **Running several t-tests instead.** Three pairwise t-tests inflate the type I error. Use ANOVA, then a multiplicity-adjusted post-hoc test.
- **Stopping at the F test.** Report group means with CIs and appropriate pairwise contrasts; a lone p-value does not tell the reader where the differences are.
- **Coding the group as a number.** In R and SAS, make sure the variable is a factor / `class` variable.
- **Testing normality of the raw outcome rather than residuals.** The assumption concerns residuals.
- **Ignoring unequal variances with unequal n.** In that situation the classic F-test can be badly miscalibrated; use Welch's ANOVA.
- **Skipping baseline adjustment in trials.** Where baseline values exist, [ANCOVA](#/test/ancova) is usually preferred.

## Related tests

- [Tukey's HSD](#/test/tukey-hsd): all pairwise comparisons after a significant ANOVA
- [Dunnett's test](#/test/dunnett): each arm vs a single control
- [Kruskal–Wallis test](#/test/kruskal-wallis): nonparametric alternative
- [Two-way ANOVA](#/test/two-way-anova): two factors
- [ANCOVA](#/test/ancova): adjust for a covariate such as baseline
- [Levene's test](#/test/levene) and [Shapiro–Wilk test](#/test/shapiro-wilk): assumption checks
- [Multiple-comparison corrections](#/test/multiple-comparisons)

## Key takeaways

- One-way ANOVA compares three or more means with F = MS_between / MS_within, where "between" is signal and "within" is noise.
- A significant F means *some* mean differs; follow up with adjusted pairwise or pre-specified contrasts.
- η² is biased upward in small samples; ω² is less biased. Report means and CIs too.
- Use Welch's ANOVA when variances or group sizes are unequal; do not let a Levene pre-test choose for you.
- p > 0.05 is not proof that the groups are equal.

## References

- Altman DG. *Practical Statistics for Medical Research.* Chapman & Hall; 1991. Ch. 9.
- ICH E9: Statistical Principles for Clinical Trials (1998) — multiplicity and pre-specification.
- R documentation: `?aov`, `?oneway.test`; SciPy documentation: `scipy.stats.f_oneway`
- SAS/STAT User's Guide: The GLM Procedure — MEANS statement
