# ANCOVA

> Compares group means on a continuous outcome after adjusting for a covariate, most often the baseline value of that outcome, and is the standard analysis of a randomized trial with a baseline measurement.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Compare groups on an outcome while adjusting for a baseline covariate | Continuous outcome plus a continuous covariate | Linear covariate effect; parallel slopes; roughly normal residuals | Week-12 SBP, active vs control, adjusted for baseline SBP |

## Learning objectives

- Explain how adjusting for baseline removes noise and corrects chance imbalance
- Run ANCOVA in R and Python and read the adjusted mean difference with its 95% CI
- Explain why change-from-baseline adjusted for baseline equals the follow-up model
- Check the parallel-slopes assumption without over-reading a non-significant interaction test

## What it is

Analysis of covariance (ANCOVA) is a linear regression of the follow-up outcome on **treatment group plus the baseline value**. Baseline values are strongly correlated with follow-up values, so adjusting for them removes much of the between-patient noise. The treatment effect is then estimated more precisely (a narrower confidence interval or a smaller required sample size) and chance baseline imbalances between arms are corrected for.

- **Model:** Y_followup = β₀ + β₁ · baseline + β₂ · treatment + ε
- **Test statistic:** t (or F = t²) for the treatment coefficient β₂, with N − 3 degrees of freedom in a two-arm trial
- **Effect size:** the **adjusted mean difference** with 95% CI (primary). It is the difference between arms in follow-up values at the same baseline value.
- **Change from baseline:** analysing the change score (follow-up − baseline) as the outcome **and** adjusting for baseline gives exactly the same treatment effect as the model above, because baseline appears on both sides and only the baseline coefficient shifts by 1. Analysing the change score **without** adjusting for baseline is not equivalent and is less efficient unless the baseline–follow-up correlation is near 1.
- **Regulatory context:** EMA guidance on adjustment for baseline covariates (2015) and FDA guidance on covariate adjustment (2023) both support a **pre-specified** adjusted analysis, with baseline value of the outcome and stratification factors as covariates, and discourage selecting covariates based on observed imbalance.

## The intuition

Patients do not start in the same place. A patient who starts at 182 mmHg will probably still be high at week 12, whatever the drug does. That baseline-driven scatter is **noise** for the question "does the drug work?". ANCOVA uses the baseline to predict where each patient would be at week 12 anyway, and then compares the arms on what is *left over*. Think of a race with a handicap: you compare finishing times after accounting for each runner's starting position.

- **Signal:** the difference between arms at the *same* baseline.
- **Noise:** the scatter of week-12 values around the fitted baseline lines, which is smaller than the raw scatter because baseline explains much of it.
- **Bonus:** if one arm happened to start a bit higher (here control: 160.5 vs 159.0 mmHg), ANCOVA corrects for that chance imbalance.

The first plot shows the unadjusted change scores. Most active patients fall below zero, but there is considerable overlap with the control arm.

<div class="viz" data-viz="dots" data-unit="mmHg, week 12 − baseline" data-groups='{"Control":[-18,-8,-8,3,0,-5,9,-3,3,4,0,-16,-4,-7,6],"Active":[-15,-11,-5,-9,-20,-17,-2,-12,-13,-9,-8,-14,-16,-6,-5]}' data-zero="1" data-caption="Change from baseline by arm (means -2.9 and -10.8). The arms differ, but individual patients overlap a lot."></div>

ANCOVA's adjusted treatment effect is −8.12 mmHg with SE 2.34, giving t = −3.47 on 27 degrees of freedom. If the drug had no effect, t would follow the curve below. Ours lies far out in the tail on the negative side.

<div class="viz" data-viz="t" data-stat="-3.47" data-df="27"></div>

Notice that the adjusted t (−3.47) is slightly larger in size than the unadjusted change-score t (−3.30). The gain here is modest because baseline-to-follow-up slope (0.83) is already close to 1, which is the case where change scores do almost as well. With a weaker baseline correlation the advantage of ANCOVA over simple change scores is bigger.

## When to use it

Use ANCOVA when **all** of these are true:

- The outcome is **continuous** and measured at baseline and at follow-up
- Subjects are **independent** and, ideally, randomized to **two or more** groups
- The relationship between baseline and follow-up is approximately **linear**, and similar in each arm

**Clinical scenario.** A randomized trial assigns 30 adults with hypertension to an active drug or control (15 per arm). Systolic blood pressure (SBP, mmHg) is measured at baseline and at week 12. The primary estimand is the *difference in week-12 SBP between arms at a given baseline SBP*, equivalent to the adjusted difference in change from baseline.

**Don't use it when:**

- There is no baseline measurement, or the covariate is itself affected by treatment → use [independent t-test](#/test/independent-t) or [one-way ANOVA](#/test/one-way-anova)
- The outcome is measured at **several** follow-up visits with missing data → use a [mixed model for repeated measures (MMRM)](#/test/mixed-models)
- The outcome is binary or time-to-event → use [logistic regression](#/test/logistic-regression) or [Cox regression](#/test/cox-ph) with baseline covariates
- The slopes clearly differ between arms (a treatment × baseline interaction) → report the interaction and effects at several baseline values, or see [linear regression](#/test/linear-regression)
- The post-baseline outcome is highly skewed with small n → transform, or use a rank-based approach

## Assumptions

1. **Independent observations.** Different patients in each group; for clustering, use a [mixed model](#/test/mixed-models).
2. **Linear relationship between the covariate and the outcome.** Check with a scatterplot of follow-up vs baseline, by arm.
3. **Homogeneity of regression slopes** (parallel lines). Test the baseline × treatment interaction; a large interaction means the treatment effect varies with baseline and a single adjusted difference is misleading. A non-significant interaction test (p > 0.05) has low power in typical trial sizes and is not proof of parallel slopes, so also inspect the plot by arm.
4. **Residuals are approximately normal with constant variance.** Check a Q–Q plot and a residual-vs-fitted plot; see [Shapiro–Wilk](#/test/shapiro-wilk) and [Levene](#/test/levene).
5. **The covariate is measured before randomization and is not influenced by treatment.** Never adjust for post-randomization variables.
6. **The covariates are pre-specified** in the protocol or statistical analysis plan.

## Hypotheses

- **H₀:** μ_active,adj = μ_control,adj (the adjusted mean follow-up SBP is the same in both arms; the adjusted difference is zero)
- **H₁:** μ_active,adj ≠ μ_control,adj (two-sided)

## Example data

SBP (mmHg) at baseline and week 12; 15 patients per arm:

| Arm | Baseline (n = 15) | Week 12 |
|---|---|---|
| Control | 168, 165, 177, 147, 182, 164, 144, 151, 161, 160, 137, 168, 155, 162, 166 | 150, 157, 169, 150, 182, 159, 153, 148, 164, 164, 137, 152, 151, 155, 172 |
| Active | 175, 167, 171, 152, 156, 164, 160, 150, 147, 158, 167, 163, 149, 153, 153 | 160, 156, 166, 143, 136, 147, 158, 138, 134, 149, 159, 149, 133, 147, 148 |

Mean baseline: 160.5 (control) and 159.0 (active), a small chance imbalance. Mean change from baseline: −2.9 (control) and −10.8 (active).

## How to run it

### R

```r
library(emmeans)
# RCT: SBP (mmHg) at baseline and week 12; n = 15 per arm
trial <- data.frame(
  arm = factor(rep(c("Control", "Active"), each = 15), levels = c("Control", "Active")),
  base = c(168,165,177,147,182,164,144,151,161,160,137,168,155,162,166,
           175,167,171,152,156,164,160,150,147,158,167,163,149,153,153),
  wk12 = c(150,157,169,150,182,159,153,148,164,164,137,152,151,155,172,
           160,156,166,143,136,147,158,138,134,149,159,149,133,147,148))
trial$chg <- trial$wk12 - trial$base          # change from baseline

# Unadjusted comparison of change scores (for contrast)
t.test(chg ~ arm, data = trial, var.equal = TRUE)

# ANCOVA: week-12 outcome adjusted for baseline (pre-specified covariate)
fit <- lm(wk12 ~ base + arm, data = trial)
summary(fit)$coefficients
confint(fit)["armActive", ]

# Adjusted (LS) means at the average baseline, and the treatment contrast
emm <- emmeans(fit, ~ arm)
as.data.frame(emm)
pairs(emm, reverse = TRUE) |> confint()

# Same treatment effect if change from baseline is the outcome
coef(lm(chg ~ base + arm, data = trial))["armActive"]

# Check parallel slopes: baseline x arm interaction
anova(lm(wk12 ~ base * arm, data = trial))["base:arm", ]
```

**What the code does**

- `factor(..., levels = c("Control", "Active"))` makes Control the **reference** level, so the `armActive` coefficient is Active − Control.
- `t.test(chg ~ arm, var.equal = TRUE)` is the unadjusted comparison of change scores, shown only for contrast. R's `t.test` defaults to Welch; `var.equal = TRUE` gives the pooled test to match the Python and linear-model results. The group order in its output is Control then Active.
- `lm(wk12 ~ base + arm)` is the ANCOVA. The `armActive` row is the **adjusted mean difference** (Active − Control) at any common baseline; the `base` row is the pooled within-arm slope.
- `emmeans` returns the **adjusted (least-squares) means**, evaluated at the overall mean baseline (159.7 mmHg), and `pairs(..., reverse = TRUE)` reports Active − Control with a CI.
- The change-score model `chg ~ base + arm` returns the **identical** treatment coefficient, showing that "change from baseline adjusted for baseline" and "follow-up adjusted for baseline" are the same analysis.
- The last line tests the baseline × arm interaction (the parallel-slopes assumption).

**Output**

```text
	Two Sample t-test

data:  chg by arm
t = 3.2973, df = 28, p-value = 0.002659
95 percent confidence interval:
  2.979568 12.753766
sample estimates:
mean in group Control  mean in group Active 
            -2.933333            -10.800000 

              Estimate Std. Error   t value     Pr(>|t|)
(Intercept) 24.8847040 18.5985725  1.337990 1.920603e-01
base         0.8266429  0.1154458  7.160443 1.062306e-07
armActive   -8.1209238  2.3402072 -3.470173 1.763542e-03
     2.5 %     97.5 % 
-12.922632  -3.319215 
 arm     emmean       SE df lower.CL upper.CL
 Control 156.9271 1.652609 27 153.5362 160.3180
 Active  148.8062 1.652609 27 145.4153 152.1971

 contrast         estimate   SE df lower.CL upper.CL
 Active - Control    -8.12 2.34 27    -12.9    -3.32

armActive 
-8.120924 

         Df Sum Sq Mean Sq F value Pr(>F)
base:arm  1 61.482  61.482  1.5345 0.2265
```

### Python

```python
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm
from scipy import stats

# RCT: SBP (mmHg) at baseline and week 12; n = 15 per arm
trial = pd.DataFrame({
    "arm": ["Control"]*15 + ["Active"]*15,
    "base": [168,165,177,147,182,164,144,151,161,160,137,168,155,162,166,
             175,167,171,152,156,164,160,150,147,158,167,163,149,153,153],
    "wk12": [150,157,169,150,182,159,153,148,164,164,137,152,151,155,172,
             160,156,166,143,136,147,158,138,134,149,159,149,133,147,148]})
trial["chg"] = trial.wk12 - trial.base          # change from baseline

# Unadjusted comparison of change scores (pooled-variance t-test)
a = trial.loc[trial.arm == "Active", "chg"]; c = trial.loc[trial.arm == "Control", "chg"]
r = stats.ttest_ind(a, c)
print(f"Unadjusted change diff = {a.mean()-c.mean():.2f}, t = {r.statistic:.3f}, p = {r.pvalue:.4f}")

# ANCOVA: week-12 outcome adjusted for baseline (Control is the reference level)
fit = smf.ols("wk12 ~ base + C(arm, Treatment('Control'))", data=trial).fit()
print(fit.summary().tables[1])

# Adjusted (LS) means at the overall mean baseline
newd = pd.DataFrame({"arm": ["Control", "Active"], "base": trial.base.mean()})
print(fit.get_prediction(newd).summary_frame()[["mean", "mean_se"]].set_index(newd.arm))

# Same treatment effect with change from baseline as the outcome
fit2 = smf.ols("chg ~ base + C(arm, Treatment('Control'))", data=trial).fit()
print(f"Change-score model, arm effect = {fit2.params.iloc[1]:.3f}")

# Check parallel slopes: baseline x arm interaction
fit3 = smf.ols("wk12 ~ base * C(arm)", data=trial).fit()
row = anova_lm(fit3).loc["base:C(arm)"]
print(f"Slope x arm: F = {row.F:.3f}, p = {row['PR(>F)']:.4f}")
```

**What the code does**

- `stats.ttest_ind(a, c)` is the pooled-variance t-test (SciPy's default, `equal_var=True`), on Active − Control change scores.
- `smf.ols("wk12 ~ base + C(arm, Treatment('Control'))")` is the ANCOVA; `Treatment('Control')` sets the reference level, otherwise Python chooses it alphabetically (Active), which would flip the sign.
- `summary().tables[1]` gives the coefficient table with 95% CIs. The arm coefficient is the adjusted mean difference.
- `get_prediction(newd)` computes **adjusted means** by predicting each arm at the average baseline, equivalent to R's `emmeans` and SAS `LSMEANS`.
- The change-score model reproduces the same treatment effect (the arm is the first non-intercept coefficient there, hence `params.iloc[1]`).
- The last block tests the parallel-slopes assumption with a baseline × arm interaction.

**Output**

```text
Unadjusted change diff = -7.87, t = -3.297, p = 0.0027
==========================================================================================================
                                             coef    std err          t      P>|t|      [0.025      0.975]
----------------------------------------------------------------------------------------------------------
Intercept                                 24.8847     18.599      1.338      0.192     -13.276      63.046
C(arm, Treatment('Control'))[T.Active]    -8.1209      2.340     -3.470      0.002     -12.923      -3.319
base                                       0.8266      0.115      7.160      0.000       0.590       1.064
==========================================================================================================
               mean   mean_se
arm                          
Control  156.927129  1.652609
Active   148.806205  1.652609
Change-score model, arm effect = -8.121
Slope x arm: F = 1.535, p = 0.2265
```

### SAS

```sas
data trial;
  input arm $ base wk12;
  chg = wk12 - base;                      /* change from baseline */
  datalines;
Control 168 150
Control 165 157
Control 177 169
Control 147 150
Control 182 182
Control 164 159
Control 144 153
Control 151 148
Control 161 164
Control 160 164
Control 137 137
Control 168 152
Control 155 151
Control 162 155
Control 166 172
Active 175 160
Active 167 156
Active 171 166
Active 152 143
Active 156 136
Active 164 147
Active 160 158
Active 150 138
Active 147 134
Active 158 149
Active 167 159
Active 163 149
Active 149 133
Active 153 147
Active 153 148
;
run;

/* ANCOVA: Control is first in the data, so order=data makes it the reference */
proc glm data=trial order=data;
  class arm;
  model wk12 = base arm / solution clparm;
  lsmeans arm / pdiff cl stderr;         /* adjusted means at mean baseline */
run;

/* Check parallel slopes: baseline x arm interaction */
proc glm data=trial order=data;
  class arm;
  model wk12 = base arm base*arm / ss3;
run;
quit;
```

**What the code does**

- `chg = wk12 - base;` creates the change score inside the DATA step; it is not needed for the primary model but shows the equivalence (`model chg = base arm;` gives the same arm effect).
- `order=data` orders the levels of `arm` as they appear in the data (Control, then Active). Without it SAS sorts alphabetically (Active, then Control). Note that SAS then treats the **last** level as the reference in `solution`, so the `arm` estimate printed there is **Control − Active** (the sign is flipped relative to R and Python).
- `model wk12 = base arm` fits the ANCOVA; `solution` prints the coefficients and `clparm` adds their 95% CIs.
- `lsmeans arm / pdiff cl stderr` prints **adjusted means** at the mean baseline, the difference between them, and its CI. With `order=data`, the `LSMEANS` difference is also labelled Control − Active (about +8.12), the opposite sign to R and Python; swap the order or negate it to report Active − Control.
- The Type III table in the first run tests `base` and `arm`; the arm F equals the square of the t in R and Python.
- The second `PROC GLM` adds `base*arm`; a small F for that term supports parallel slopes.

## Interpreting the output

- **Unadjusted change scores:** −2.9 mmHg (control) vs −10.8 mmHg (active), a difference of −7.87 mmHg (95% CI −12.75 to −2.98, t = −3.30, p = 0.003; Active − Control, as in Python. R's `t.test` prints Control − Active, so its signs are reversed).
- **ANCOVA treatment effect:** the adjusted mean difference is **−8.12 mmHg** (95% CI −12.92 to −3.32; SE 2.34; t(27) = −3.47, p = 0.002). At the same baseline SBP, patients on the active drug have a week-12 SBP about 8 mmHg lower.
- **Adjusted means:** 156.9 mmHg (control) and 148.8 mmHg (active) at the average baseline of 159.7 mmHg, each with SE 1.65.
- **Baseline slope:** 0.83 mmHg of week-12 SBP per mmHg of baseline (p < 0.001). Patients with higher baseline stay higher, which is why adjustment helps. The SE for the treatment effect is 2.34, against 2.39 for the unadjusted change-score comparison (−7.87 / −3.30), a modest gain here because the slope is already close to 1.
- **Same effect, change scores:** the change-score model adjusted for baseline gives −8.121, identical to the follow-up model.
- **Parallel slopes:** baseline × arm F(1, 26) = 1.53, p = 0.23, so there is no evidence that the baseline relationship differs between arms (with only 30 patients this test is weak).
- **Why adjusted ≠ unadjusted:** the control arm happened to start slightly higher (160.5 vs 159.0 mmHg), which ANCOVA corrects for.

**How to report it:**
> At week 12, adjusted for baseline systolic blood pressure, the mean SBP was 8.1 mmHg lower with the active drug than with control (148.8 vs 156.9 mmHg; adjusted difference −8.12 mmHg, 95% CI −12.92 to −3.32; ANCOVA, p = 0.002).

⚠️ **Important for trials:** the covariates must be pre-specified in the protocol or statistical analysis plan (ICH E9). Choosing covariates after looking at the data, or adjusting for post-randomization variables, invalidates the type I error rate.

## Check your understanding

<details class="check"><summary>1. A trial has baseline and week-12 LDL-C. A colleague compares change scores with a t-test and says "that already accounts for baseline." Is that equivalent to ANCOVA?</summary>

No. Change scores subtract baseline with a fixed coefficient of 1. ANCOVA estimates the coefficient from the data. Unless the baseline-to-follow-up slope is about 1, ANCOVA is more precise, and if arms differ at baseline by chance the change-score analysis can be biased. Putting baseline *in the model* as well as using change as the outcome gives the same answer as ANCOVA.
</details>

<details class="check"><summary>2. The baseline × arm interaction has p = 0.23. Does that prove the slopes are parallel?</summary>

No. With 30 patients the interaction test has little power, so p = 0.23 is only a lack of evidence against parallel slopes. Plot follow-up against baseline by arm and judge whether the lines look reasonably parallel, and prefer pre-specifying the model rather than choosing it from this test.
</details>

<details class="check"><summary>3. Active arm patients with high baseline SBP responded best. After seeing the data, you decide to add "early response at week 4" as a covariate. Is that fine?</summary>

No. Week-4 response happens after randomization and is partly caused by treatment, so adjusting for it removes part of the treatment effect. Covariates must be measured before randomization and pre-specified.
</details>

<details class="check"><summary>4. What if the baseline-to-week-12 slope were 0.2 instead of 0.83? Would adjusting still help?</summary>

Much less. Adjustment removes variance in proportion to the squared baseline–outcome correlation, so a weak relationship removes little noise. ANCOVA would still be unbiased and never worse than the unadjusted comparison, but the precision gain would be small. Change scores would then perform *worse* than the unadjusted follow-up comparison.
</details>

## Common mistakes

- **Analysing change scores without adjusting for baseline.** This is less efficient and, in the presence of baseline imbalance, can be biased. Adjust for baseline in the model, whatever the outcome definition.
- **Testing baseline differences between randomized arms and adjusting only if "significant".** Any baseline difference in a randomized trial is due to chance; adjust for the pre-specified prognostic covariates regardless.
- **Using percent change from baseline as the outcome.** It is skewed, has unstable variance near zero baselines, and does not remove the baseline correlation.
- **Ignoring non-parallel slopes.** If the baseline × treatment interaction is large, a single adjusted difference does not describe the effect.
- **Sign and reference-level confusion.** R and Python use the first level as reference (here set explicitly); SAS `solution` uses the last. State the direction (Active − Control) in the report.
- **Adjusting for a covariate measured after randomization** (such as an early response), which can remove part of the treatment effect.

## Related tests

- [Independent t-test](#/test/independent-t): the unadjusted two-group comparison
- [One-way ANOVA](#/test/one-way-anova): the unadjusted comparison of 2+ groups
- [Linear regression](#/test/linear-regression): ANCOVA is a regression with a categorical and a continuous predictor
- [Linear mixed-effects models (MMRM)](#/test/mixed-models): baseline-adjusted analysis with multiple follow-up visits and missing data
- [Paired t-test](#/test/paired-t): single-arm before/after change, with no control group
- [Two-way ANOVA](#/test/two-way-anova): two categorical factors rather than a continuous covariate

## Key takeaways

- ANCOVA = regression of follow-up on treatment plus baseline. The treatment coefficient is the adjusted mean difference, reported with a 95% CI.
- Adjusting for baseline cuts noise and corrects chance imbalance; it is the standard analysis for randomized trials with a baseline measurement.
- Change score adjusted for baseline gives exactly the same treatment effect as the follow-up model.
- Covariates must be measured before randomization and pre-specified, never chosen by looking at the data.
- A non-significant slope-by-treatment test does not prove parallel slopes. Look at the plot.

## References

- European Medicines Agency. *Guideline on adjustment for baseline covariates in clinical trials* (EMA/CHMP/295050/2013). 2015.
- US FDA. *Adjusting for Covariates in Randomized Clinical Trials for Drugs and Biological Products: Guidance for Industry.* 2023.
- Vickers AJ, Altman DG. Analysing controlled trials with baseline and follow up measurements. *BMJ.* 2001;323:1123–4.
- R documentation: `?lm`, and the `emmeans` package vignette; SAS/STAT User's Guide: The GLM Procedure — LSMEANS statement
