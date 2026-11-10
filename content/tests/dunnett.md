# Dunnett's test

> Compares each of several treatment arms with a single control arm, while keeping the family-wise error rate at 5% across all of those comparisons.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Compare each of several treatment arms with one control | Continuous outcome, independent groups | Normal residuals; equal variances; comparisons pre-specified | Three LDL-C doses each vs placebo |

## Learning objectives

- Explain why several unadjusted comparisons against one control inflate the false-positive rate
- Recognize when Dunnett applies (treatment vs control) and when Tukey or Holm is needed instead
- Run Dunnett's test in R and Python and read adjusted p-values and simultaneous CIs
- Report each arm versus control with a family-wise adjusted CI

## What it is

When a trial has one control (placebo or standard of care) and k active arms, the natural questions are the k comparisons *"each arm vs control"*. Running k separate t-tests inflates the chance of at least one false positive. Dunnett's test fixes this by using the **joint distribution** of the k test statistics: the comparisons all share the same control group and therefore are positively correlated, so the exact critical value (a multivariate t quantile) is smaller than the one Bonferroni would use. The result is more power than Bonferroni or Tukey for this specific set of comparisons.

- **Test statistic:** for each arm i, t_i = (ȳ_i − ȳ_control) / (s_p · √(1/n_i + 1/n_control)), where s_p² is the pooled within-group variance (the ANOVA mean square error)
- **Reference distribution:** the multivariate t with the pooled error degrees of freedom; the p-values and CIs are **adjusted** so that the family-wise error rate is controlled at α
- **Effect size:** the mean difference from control for each arm with a **simultaneous** 95% CI (primary); Cohen's d against the pooled SD (secondary)

## The intuition

Each comparison with placebo is a signal-to-noise ratio: the signal is a dose's mean difference from placebo, the noise is the standard error of that difference (2.09 mg/dL here). The twist is **multiplicity**. Test three doses at α = 0.05 each and the chance that at least one is a false positive is about 1 − 0.95³ = 14%. Think of three smoke detectors that each have a 5% false-alarm rate: with all three in one house, false alarms become routine.

Dunnett's fix is to raise the bar for "signal is big enough". Because all three comparisons share the same placebo group, the three t statistics are correlated, so the exact bar can be lower than Bonferroni's. Here the 95% family-wise critical value is 2.48. An ordinary single t-test would use 2.05. That is the price of looking at three doses at once.

<div class="viz" data-viz="dots" data-unit="mg/dL change in LDL-C" data-groups='{"Placebo":[-7.8,-4.5,-1.7,-8.8,-2.0,-2.8,-2.6,2.6],"10 mg":[-12.1,0.3,-9.7,-11.7,-9.6,-4.7,-5.2,-7.5],"20 mg":[-13.8,-12.2,-2.9,-8.0,-11.9,-13.7,-10.0,-17.3],"40 mg":[-13.4,-14.7,-5.2,-5.9,-11.4,-16.7,-6.5,-6.7]}' data-zero="1" data-caption="Notice how the 20 mg and 40 mg dots sit clearly lower than placebo, while the 10 mg dots overlap a lot with placebo. The differences from placebo (-4.1, -7.8, -6.6 mg/dL) are compared with the spread within arms (pooled SD 4.2)."></div>

What to notice: the 10 mg arm is about 4 mg/dL below placebo, which is less than 2 standard errors, so even before adjustment it is weak evidence (adjusted p = 0.149). The 20 mg arm is 3.7 standard errors away, which is beyond 2.48, so it still counts as significant after adjustment (adjusted p = 0.0026).

## When to use it

Use Dunnett's test when **all** of these are true:

- The outcome is **continuous**
- There is **one control group** and **two or more** treatment groups, and the planned comparisons are *each treatment vs control* only
- The groups contain **different subjects**

**Clinical scenario.** A phase 2 dose-ranging trial randomizes 32 patients to placebo or to 10, 20, or 40 mg of a lipid-lowering drug (8 per arm). The outcome is the change in LDL-C (mg/dL) after 12 weeks. *Which doses reduce LDL-C more than placebo?* The three comparisons are pre-specified, so Dunnett's adjustment applies.

**Don't use it when:**

- You want **all** pairwise comparisons (including dose vs dose) → use [Tukey's HSD](#/test/tukey-hsd)
- You have only **two** groups → use the [independent t-test](#/test/independent-t)
- The comparisons are not all vs a single control, or involve different endpoints → use [Bonferroni, Holm, or Benjamini–Hochberg corrections](#/test/multiple-comparisons)
- The variances are clearly unequal or the data are strongly skewed → use a rank-based approach such as [Kruskal–Wallis](#/test/kruskal-wallis) followed by corrected [Mann–Whitney tests](#/test/mann-whitney), or a robust variant
- The same patients are measured in several conditions → use [repeated-measures ANOVA](#/test/repeated-measures-anova) or a [mixed model](#/test/mixed-models)

## Assumptions

1. **Independent observations.** Different patients in each arm.
2. **Residuals are approximately normal** within each group. Check with a Q–Q plot of residuals or [Shapiro–Wilk](#/test/shapiro-wilk); with n = 8 per arm, plots are more informative than tests.
3. **Equal variances across arms.** Dunnett's test uses the pooled variance. Check with [Levene's test](#/test/levene) or by comparing group SDs.
4. **Pre-specified comparisons.** The family is "each treatment vs the control". Do not add other comparisons after the fact; they are not covered by this adjustment.
5. **A true shared control group** with similar sample size, ideally; unequal sizes are allowed but the control arm is often made larger (about √k times the treatment arms) to optimize power.

## Hypotheses

For each treatment arm i (i = 1, …, k):

- **H₀ᵢ:** μᵢ = μ_control
- **H₁ᵢ:** μᵢ ≠ μ_control (two-sided; one-sided versions exist when a direction is pre-specified)

The family-wise error rate (the chance of rejecting at least one true null) is controlled at α = 0.05.

## Example data

Change in LDL-C (mg/dL) after 12 weeks; negative values are reductions. n = 8 per arm:

| Arm | Individual values | Mean |
|---|---|---|
| Placebo | −7.8, −4.5, −1.7, −8.8, −2.0, −2.8, −2.6, 2.6 | −3.45 |
| 10 mg | −12.1, 0.3, −9.7, −11.7, −9.6, −4.7, −5.2, −7.5 | −7.53 |
| 20 mg | −13.8, −12.2, −2.9, −8.0, −11.9, −13.7, −10.0, −17.3 | −11.23 |
| 40 mg | −13.4, −14.7, −5.2, −5.9, −11.4, −16.7, −6.5, −6.7 | −10.06 |

## How to run it

### R

```r
library(multcomp)
# Change in LDL-C (mg/dL) after 12 weeks; placebo plus three doses, n = 8 each
dose <- factor(rep(c("Placebo", "D10", "D20", "D40"), each = 8),
               levels = c("Placebo", "D10", "D20", "D40"))   # control FIRST
y <- c(-7.8,-4.5,-1.7,-8.8,-2.0,-2.8,-2.6,2.6,     # Placebo
       -12.1,0.3,-9.7,-11.7,-9.6,-4.7,-5.2,-7.5,   # D10
       -13.8,-12.2,-2.9,-8.0,-11.9,-13.7,-10.0,-17.3, # D20
       -13.4,-14.7,-5.2,-5.9,-11.4,-16.7,-6.5,-6.7)   # D40
dat <- data.frame(dose, y)

# Omnibus one-way ANOVA first (optional; Dunnett does not require it)
fit <- aov(y ~ dose, data = dat)
summary(fit)

# Dunnett: each dose vs the first factor level (Placebo)
dn <- glht(fit, linfct = mcp(dose = "Dunnett"))
summary(dn)
confint(dn)
```

**What the code does**

- The factor `levels` list puts **Placebo first**. `mcp(dose = "Dunnett")` always treats the **first level as the control**, so level order decides which arm is the reference.
- `aov(y ~ dose)` fits the one-way ANOVA; its residual mean square (17.48, df = 28) is the pooled variance Dunnett uses.
- `glht(fit, linfct = mcp(dose = "Dunnett"))` from the `multcomp` package builds the three treatment-minus-control contrasts and computes **single-step** adjusted p-values from the multivariate t distribution.
- `summary(dn)` gives the estimates, SEs, t-values, and **adjusted** p-values; `confint(dn)` gives **simultaneous** 95% CIs (they are wider than ordinary CIs, using a critical value of 2.48 instead of 2.05).
- The adjusted p-values come from numerical integration, so the last digits can differ slightly between runs or between software packages.

**Output**

```text
            Df Sum Sq Mean Sq F value  Pr(>F)   
dose         3  284.5   94.84   5.427 0.00453 **
Residuals   28  489.3   17.48                   

	 Simultaneous Tests for General Linear Hypotheses

Multiple Comparisons of Means: Dunnett Contrasts

Linear Hypotheses:
                   Estimate Std. Error t value Pr(>|t|)   
D10 - Placebo == 0   -4.075      2.090  -1.950  0.14866   
D20 - Placebo == 0   -7.775      2.090  -3.720  0.00255 **
D40 - Placebo == 0   -6.612      2.090  -3.164  0.01013 * 
(Adjusted p values reported -- single-step method)

	 Simultaneous Confidence Intervals

Quantile = 2.4832
95% family-wise confidence level

Linear Hypotheses:
                   Estimate lwr      upr     
D10 - Placebo == 0  -4.0750  -9.2652   1.1152
D20 - Placebo == 0  -7.7750 -12.9652  -2.5848
D40 - Placebo == 0  -6.6125 -11.8027  -1.4223
```

### Python

```python
import numpy as np
from scipy.stats import dunnett

# Change in LDL-C (mg/dL) after 12 weeks; placebo plus three doses, n = 8 each
placebo = np.array([-7.8,-4.5,-1.7,-8.8,-2.0,-2.8,-2.6,2.6])
d10 = np.array([-12.1,0.3,-9.7,-11.7,-9.6,-4.7,-5.2,-7.5])
d20 = np.array([-13.8,-12.2,-2.9,-8.0,-11.9,-13.7,-10.0,-17.3])
d40 = np.array([-13.4,-14.7,-5.2,-5.9,-11.4,-16.7,-6.5,-6.7])

# Dunnett: treatment samples first, control passed by keyword
res = dunnett(d10, d20, d40, control=placebo, random_state=1)   # two-sided; fixed seed
ci = res.confidence_interval(confidence_level=0.95)
for name, diff, p, lo, hi in zip(["D10", "D20", "D40"],
        [x.mean() - placebo.mean() for x in (d10, d20, d40)],
        res.pvalue, ci.low, ci.high):
    print(f"{name} - Placebo: diff = {diff:6.2f}, 95% CI ({lo:6.2f}, {hi:6.2f}), adj p = {p:.4f}")
```

**What the code does**

- `scipy.stats.dunnett` (SciPy 1.11 or later) takes the **treatment** samples as positional arguments and the **control** as the `control=` keyword; the order is the reverse of R's "control is the first level".
- The default `alternative="two-sided"` matches R; use `"less"` or `"greater"` for a pre-specified one-sided test.
- `res.pvalue` holds the adjusted p-values (one per treatment) and `res.statistic` the t statistics. `confidence_interval()` returns **simultaneous** CIs for each treatment-minus-control difference.
- SciPy computes the multivariate t probabilities by numerical integration with a random component, so p-values can differ from R in the third or fourth significant digit (here 0.0025 vs 0.00255). The code sets `random_state` for reproducibility.
- The code takes only the treatment-versus-control comparisons; SciPy does not run the omnibus ANOVA for you (use `scipy.stats.f_oneway` if you want it).

**Output**

```text
D10 - Placebo: diff =  -4.07, 95% CI ( -9.27,   1.12), adj p = 0.1487
D20 - Placebo: diff =  -7.78, 95% CI (-12.97,  -2.58), adj p = 0.0025
D40 - Placebo: diff =  -6.61, 95% CI (-11.80,  -1.42), adj p = 0.0102
```

### SAS

```sas
data ldl;
  input dose $ y @@;
  datalines;
Placebo -7.8  Placebo -4.5  Placebo -1.7  Placebo -8.8
Placebo -2.0  Placebo -2.8  Placebo -2.6  Placebo 2.6
D10 -12.1  D10 0.3  D10 -9.7  D10 -11.7
D10 -9.6  D10 -4.7  D10 -5.2  D10 -7.5
D20 -13.8  D20 -12.2  D20 -2.9  D20 -8.0
D20 -11.9  D20 -13.7  D20 -10.0  D20 -17.3
D40 -13.4  D40 -14.7  D40 -5.2  D40 -5.9
D40 -11.4  D40 -16.7  D40 -6.5  D40 -6.7
;
run;

/* One-way ANOVA with Dunnett comparisons against the Placebo arm */
proc glm data=ldl;
  class dose;
  model y = dose;
  lsmeans dose / pdiff=control('Placebo') adjust=dunnett cl alpha=0.05;
run;
quit;
```

**What the code does**

- `input dose $ y @@;` reads several (dose, y) pairs per line (the double trailing `@@` keeps the same data line active), which keeps the datalines compact.
- `model y = dose;` fits the one-way ANOVA; the omnibus F test appears first in the output.
- `lsmeans dose / pdiff=control('Placebo') adjust=dunnett` compares each level with the named control level and applies Dunnett's adjustment. The control level label is **case-sensitive** and in quotes.
- `cl` adds **simultaneous** confidence limits for the differences, adjusted in the same way (without `adjust=`, they would be unadjusted).
- `alpha=0.05` is the default; it is shown for clarity.
- SAS labels each difference by the pair of levels compared; check the sign in the output so that you report treatment − control, as R and Python do.
- You can run `proc glm ... means dose / dunnett('Placebo') cldiff;` to get the same comparisons from the `MEANS` statement.

## Interpreting the output

- **Omnibus ANOVA:** F(3, 28) = 5.43, p = 0.005, so the arms do not all have the same mean. Dunnett's test does not require a significant omnibus test.
- **10 mg vs placebo:** the mean difference is −4.08 mg/dL (simultaneous 95% CI −9.27 to 1.12), adjusted p = 0.149. This is not significant, and the CI includes 0.
- **20 mg vs placebo:** −7.78 mg/dL (95% CI −12.97 to −2.58), adjusted p = 0.00255 (0.0025 in SciPy; the difference reflects numerical integration).
- **40 mg vs placebo:** −6.61 mg/dL (95% CI −11.80 to −1.42), adjusted p = 0.010.
- **SE of each difference:** 2.09 mg/dL, from the pooled SD √17.48 = 4.18 and n = 8 per arm. Simultaneous CIs use t = 2.48 instead of 2.05, which is the price of testing three doses.
- **Dose-response:** the 40 mg arm does not reduce LDL-C more than the 20 mg arm here. Dunnett's test does not compare doses with each other; use a pre-specified dose-response model if that is the question.

**How to report it:**
> Compared with placebo, LDL-C change was −4.1 mg/dL with 10 mg (95% CI −9.3 to 1.1; adjusted p = 0.15), −7.8 mg/dL with 20 mg (95% CI −13.0 to −2.6; adjusted p = 0.003), and −6.6 mg/dL with 40 mg (95% CI −11.8 to −1.4; adjusted p = 0.010). Dunnett's test was used to control the family-wise error rate at 5% across the three dose comparisons.

## Check your understanding

<details class="check"><summary>1. A trial compares four doses with placebo on LDL-C and also wants to know whether 40 mg beats 20 mg. Is Dunnett's test enough?</summary>

No. Dunnett's guarantee covers only the "each dose vs control" family. Adding dose-vs-dose comparisons changes the family, so use [Tukey's HSD](#/test/tukey-hsd) for all pairs, or pre-specify a hierarchy or a [Holm/Bonferroni](#/test/multiple-comparisons) adjustment that covers every planned comparison.
</details>

<details class="check"><summary>2. For 10 mg vs placebo the adjusted p is 0.149 and the CI is −9.3 to +1.1. A report says "10 mg has no effect on LDL-C." What would you say?</summary>

A non-significant result is not evidence of no effect. The simultaneous CI is compatible with a reduction as big as 9.3 mg/dL and with a small increase. With only 8 patients per arm the trial is simply too small to settle the question for 10 mg. Say "no significant difference was detected; the estimate was −4.1 mg/dL (CI −9.3 to 1.1)".
</details>

<details class="check"><summary>3. In a different three-dose trial, a colleague runs three ordinary t-tests against placebo, gets p = 0.034, 0.002 and 0.18, then claims two doses are significant. What is the problem?</summary>

Without adjustment the family-wise error rate is about 14%, not 5%, and with a cherry-picked "significant" arm the false-positive risk is higher still. Unadjusted pairwise t-tests also use only two groups' variances, not the pooled variance. Use Dunnett's test (adjusted p-values and simultaneous CIs) instead.
</details>

<details class="check"><summary>4. Suppose a fifth arm (60 mg) is added. Would the critical value 2.48 stay the same?</summary>

No, it would rise. With four comparisons instead of three, the family-wise critical value gets larger (and the error df change), so each comparison needs more signal to be declared significant. More comparisons mean a higher price for protecting the family-wise error rate.
</details>

## Common mistakes

- **Wrong control level.** In R the **first factor level** is the control; if the factor is alphabetical, "D10" would become the control. Always set `levels = ` explicitly.
- **Running unadjusted t-tests against the control** and calling the result "significant at 0.05". Three comparisons carry roughly a 14% chance of at least one false positive.
- **Reporting unadjusted CIs with adjusted p-values.** Use simultaneous CIs so that they agree with the decision.
- **Using Dunnett for all pairwise comparisons.** It covers only treatment-vs-control; add dose-vs-dose comparisons and the guarantee is lost. Use [Tukey's HSD](#/test/tukey-hsd).
- **Waiting for a significant omnibus F.** Dunnett's single-step procedure controls the error rate without this gatekeeping step; requiring one is allowed but unnecessary.
- **Not pre-specifying the family.** Under ICH E9 the comparisons and the multiplicity strategy belong in the statistical analysis plan.

## Related tests

- [Tukey's HSD](#/test/tukey-hsd): all pairwise comparisons with family-wise control
- [One-way ANOVA](#/test/one-way-anova): the omnibus test that precedes or accompanies Dunnett
- [Multiple-comparison corrections (Bonferroni, Holm, BH)](#/test/multiple-comparisons): general-purpose adjustments for any set of p-values
- [Independent t-test](#/test/independent-t): the unadjusted two-group comparison
- [Kruskal–Wallis test](#/test/kruskal-wallis): nonparametric omnibus alternative
- [Levene's test](#/test/levene): checking the equal-variance assumption

## Key takeaways

- Dunnett's test compares each treatment arm with one control while holding the family-wise error rate at 5%.
- Three unadjusted comparisons give about a 14% chance of at least one false positive.
- It uses the pooled variance and the correlation induced by the shared control, so its critical value (about 2.48 for 3 arms here) is lower than Bonferroni's.
- It covers treatment-vs-control only. For all pairs use Tukey; for other families use Holm or Bonferroni.
- Report mean differences with simultaneous (adjusted) CIs. A non-significant arm is not shown to be ineffective.

## References

- Dunnett CW. A multiple comparison procedure for comparing several treatments with a control. *J Am Stat Assoc.* 1955;50:1096–1121.
- Hothorn T, Bretz F, Westfall P. Simultaneous inference in general parametric models. *Biometrical Journal.* 2008;50:346–363 (the `multcomp` package).
- ICH E9 *Statistical Principles for Clinical Trials* (1998), section on multiplicity.
- SciPy documentation: `scipy.stats.dunnett`; SAS/STAT User's Guide: The GLM Procedure — LSMEANS statement, ADJUST=DUNNETT
