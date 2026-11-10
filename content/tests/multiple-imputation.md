# Missing data and multiple imputation

> Handles missing values by filling them in several times with plausible draws, analysing each completed dataset, and combining the results with Rubin's rules so the uncertainty about the missing values is carried through.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Outcomes or covariates are missing and the reason is plausibly related to observed data | Any (continuous, binary, count) with missing values | Data missing at random (MAR); imputation model includes the predictors of missingness | Week-12 blood pressure missing in patients who dropped out after a poor week-6 result |

## Learning objectives

- Distinguish **MCAR, MAR, and MNAR** and explain why complete-case analysis can be biased
- Explain what multiple imputation does (and does not do) and how Rubin's rules combine the results
- Run MI in R (`mice`) and Python (`statsmodels`) and compare with a complete-case analysis
- Link missing-data handling to the **estimand** and plan sensitivity analyses

## What it is

Almost every real trial has missing data. **Multiple imputation (MI)** is a three-step recipe:

1. **Impute.** Create m completed datasets. Each missing value is replaced by a random draw from its predictive distribution given the observed data, so the m datasets differ where the data were missing.
2. **Analyse.** Fit the *same* analysis model (here a linear regression) to each completed dataset.
3. **Pool.** Combine the m estimates with **Rubin's rules**.

- **Pooled estimate:** Q̄ = the average of the m estimates
- **Within-imputation variance:** Ū = the average of the m squared standard errors
- **Between-imputation variance:** B = the variance of the m estimates
- **Total variance:** T = Ū + (1 + 1/m)·B. The pooled standard error is √T.
- **Fraction of missing information (FMI)** summarizes how much of the uncertainty is due to missingness.

MI does **not create information**. The imputed values are not data. The extra variance term B is exactly the price of not knowing the missing values, and it is what makes standard errors honest.

## The intuition

Three kinds of missingness, in plain words:

- **MCAR (missing completely at random):** whether a value is missing has nothing to do with anything: a lab sample dropped on the floor. Complete-case analysis is unbiased but loses power.
- **MAR (missing at random):** missingness depends on things you **observed**, not on the missing value itself once those are accounted for. Example: patients with a high week-6 blood pressure are more likely to drop out of a trial, and week-6 is recorded. MI and likelihood-based methods are valid here **if the imputation model includes those observed variables**.
- **MNAR (missing not at random):** missingness depends on the unobserved value itself, even after accounting for what you observed: patients stop coming because they feel worse, and nobody measured how they feel. No data-only method can fix this; you need assumptions and sensitivity analyses.

You can never prove MAR from the data. It is an assumption you justify clinically and then stress-test.

Why complete-case analysis can mislead: if the patients who drop out differ systematically from those who stay, the completers are no longer a random sample of the randomized groups. In our example, more placebo patients with a poor week-6 result dropped out, so the remaining placebo patients look better than the placebo group really was, and the treatment effect is **shrunk toward zero**.

The dot plot shows the observed week-12 values; fewer than two thirds of the patients are in it.

<div class="viz" data-viz="dots" data-unit="mmHg, week-12 SBP" data-groups='{"Placebo (24 of 50 observed)":[136,140,133,155,157,128,152,159,164,153,142,150,158,151,134,141,155,164,167,150,131,158,145,139],"Treatment (38 of 50 observed)":[138,149,157,156,141,152,141,128,123,162,153,151,173,135,133,127,145,116,151,133,142,148,144,144,151,141,147,136,136,146,143,132,152,151,153,157,150,141]}' data-caption="Only completers are shown. Only 24 of 50 placebo patients are observed, versus 38 of 50 treated patients. Placebo patients who dropped out had high blood pressure, so the observed placebo group looks better than the whole placebo group really was."></div>

## When to use it

- You have **missing values** in the outcome or covariates, and a MAR assumption is clinically plausible
- You can name the **variables that predict** both missingness and the outcome (interim measurements, baseline severity, centre) and include them
- The proportion missing is moderate (up to roughly 30–40%); with more, the results depend heavily on the assumptions

**Clinical scenario.** A 100-patient trial randomizes patients to a new antihypertensive or placebo and measures systolic blood pressure (SBP) at baseline, week 6, and week 12. Thirty-eight patients have no week-12 value, mostly placebo patients who dropped out after a high week-6 reading. *What is the treatment effect on week-12 SBP, adjusted for baseline?*

**The estimand link.** ICH E9(R1) asks you to define the **estimand** first: the treatment effect you want, including how you handle intercurrent events such as treatment discontinuation or rescue medication. Missing data are a *separate* problem: "how do I handle the values I did not observe, given the estimand?" For instance, with a "treatment policy" estimand you want outcomes even after discontinuation, so MI should impute what would have been observed had patients been followed up; a "hypothetical" estimand imputes values under the situation where the intercurrent event did not occur. MI is one tool; the chosen estimand tells you what you are imputing.

**Don't use it when:**

- Data are likely **MNAR** and you have no sensitivity analysis plan: MI under MAR may be as biased as complete cases
- Only a few outcomes are missing completely at random: a complete-case analysis or a [mixed model](#/test/mixed-models) (MMRM) that uses all observed repeated measures may be simpler
- You are tempted to fill gaps once with the mean or the last observation: these "single imputation" shortcuts understate uncertainty (see Common mistakes)
- The analysis is a pure comparison with no covariates and few missing: see [independent t-test](#/test/independent-t) or [ANCOVA](#/test/ancova) directly

## Assumptions

1. **Missing at random (MAR).** Given the variables in the imputation model, missingness does not depend on the missing value. Check by clinical reasoning and by comparing completers with dropouts on observed variables. This cannot be tested formally.
2. **The imputation model is adequate and "congenial" with the analysis model.** It must include every variable in the analysis model (including the outcome and treatment), and any variable that predicts missingness or the outcome. Interactions or non-linearities in the analysis model must be in the imputation model.
3. **Enough imputations.** Use m of about 20–100; a common rule is m at least equal to the percentage of incomplete cases (here 38%, so m = 20 is a teaching-size choice; a real analysis would use 40 or more).
4. **Imputed values are plausible.** Compare imputed and observed distributions, and check that imputations lie in a realistic range.

## Hypotheses

MI is a missing-data **method**, not a test. The hypotheses belong to the analysis model; here the same regression is used throughout.

- **Estimated:** the treatment effect β (difference in mean week-12 SBP, treatment − placebo, adjusted for baseline SBP)
- **H₀:** β = 0 vs **H₁:** β ≠ 0, tested with the pooled estimate and the Rubin standard error, using the Barnard–Rubin degrees of freedom

## Example data

A simulated trial of 100 patients (50 per arm). **True treatment effect on week-12 SBP is −8 mmHg.** Baseline and week-6 SBP are known for everyone. Week-12 SBP is missing for 38 patients, with dropout more likely after a high week-6 value, which is **MAR**. The 100 patients' values are typed into the code so that R and Python use the same data.

## How to run it

### R

```r
library(mice)

# Simulated trial, n = 100: SBP (mmHg) at baseline, week 6 (all observed),
# and week 12 (outcome; NA = dropped out). True treatment effect = -8 mmHg.
trt <- rep(0:1, each = 50)
base <- scan(text = "
168 165 177 147 182 164 144 151 161 160 137 168 155 162 166 175 167
171 152 156 164 160 150 147 158 167 163 149 153 153 142 156 160 169
176 161 178 161 174 175 161 145 160 158 181 162 165 177 164 147 154
178 158 156 166 170 160 151 169 172 155 169 139 143 178 168 155 179
156 144 162 149 165 142 165 161 162 168 153 163 171 164 157 158 161
151 156 149 152 173 157 166 163 161 171 154 175 159 160 139")
wk6 <- scan(text = "
159 143 173 141 178 163 136 153 163 161 124 170 140 156 160 172 155
166 145 153 158 152 146 145 154 170 154 140 146 157 140 161 154 173
174 167 167 157 160 164 162 137 150 160 173 148 165 167 158 147 148
165 141 147 161 165 157 137 163 156 147 160 129 125 166 163 157 175
136 133 160 131 153 120 152 149 152 161 147 144 162 152 148 147 149
142 148 141 143 167 152 161 147 150 162 153 163 157 153 140")
wk12 <- scan(text = "
NA 136 NA 140 NA NA 133 155 NA 157 128 NA NA 152 159 164 153 NA 142
150 158 151 NA NA NA NA NA 134 NA NA 141 NA 155 NA NA NA NA NA 164 167
150 131 NA 158 NA 145 NA NA NA 139 138 NA NA 149 157 156 NA 141 152 NA
141 NA 128 123 162 153 151 173 135 133 NA 127 145 116 151 133 142 148
144 144 NA NA 151 141 147 NA NA 136 136 NA 146 143 132 152 151 153 NA
157 150 141", na.strings = "NA")
d <- data.frame(trt, base, wk6, wk12)
cat("Missing outcomes:", sum(is.na(d$wk12)), "of 100\n")

# 1. Complete-case analysis: lm() silently drops rows with a missing outcome
cc <- lm(wk12 ~ trt + base, data = d)
print(round(summary(cc)$coef["trt", ], 3))

# 2. Multiple imputation: m = 20 datasets, normal linear model, then pool
imp <- mice(d, m = 20, method = "norm", seed = 1, printFlag = FALSE)
fit <- with(imp, lm(wk12 ~ trt + base))
print(summary(pool(fit), conf.int = TRUE)[2, c("term","estimate","std.error","conf.low","conf.high")])
print(round(pool(fit)$pooled[2, c("b", "ubar", "t", "fmi")], 4))
```

**What the code does**

- `scan(text = ...)` reads the numbers typed into the script; `NA` marks the missing outcomes.
- `lm()` drops rows with a missing outcome automatically: this is the **complete-case** analysis.
- `mice(..., m = 20, method = "norm")` creates 20 imputed datasets using Bayesian linear regression for the one incomplete variable. Every column in `d` (including `wk6`) is used as a predictor, which is what makes the MAR assumption work here. `seed` makes it reproducible.
- `with(imp, lm(...))` fits the analysis model in each dataset; `pool()` applies Rubin's rules.
- In the second table, `b` = between-imputation variance, `ubar` = within-imputation variance, `t` = total variance, and `fmi` = fraction of missing information.

**Output**

```text
Missing outcomes: 38 of 100
  Estimate Std. Error    t value   Pr(>|t|) 
    -5.056      1.831     -2.762      0.008 
  term  estimate std.error  conf.low conf.high
2  trt -6.510841  1.607863 -9.732486 -3.289196
       b   ubar      t    fmi
2 0.6539 1.8987 2.5852 0.2907
```

### Python

```python
import numpy as np, pandas as pd, statsmodels.api as sm
from statsmodels.imputation.mice import MICEData, MICE

# Same simulated trial as the R tab (typed inline): SBP in mmHg
trt = np.repeat([0, 1], 50)
base = np.array("""
168 165 177 147 182 164 144 151 161 160 137 168 155 162 166 175 167
171 152 156 164 160 150 147 158 167 163 149 153 153 142 156 160 169
176 161 178 161 174 175 161 145 160 158 181 162 165 177 164 147 154
178 158 156 166 170 160 151 169 172 155 169 139 143 178 168 155 179
156 144 162 149 165 142 165 161 162 168 153 163 171 164 157 158 161
151 156 149 152 173 157 166 163 161 171 154 175 159 160 139
""".split(), dtype=float)
wk6 = np.array("""
159 143 173 141 178 163 136 153 163 161 124 170 140 156 160 172 155
166 145 153 158 152 146 145 154 170 154 140 146 157 140 161 154 173
174 167 167 157 160 164 162 137 150 160 173 148 165 167 158 147 148
165 141 147 161 165 157 137 163 156 147 160 129 125 166 163 157 175
136 133 160 131 153 120 152 149 152 161 147 144 162 152 148 147 149
142 148 141 143 167 152 161 147 150 162 153 163 157 153 140
""".split(), dtype=float)
wk12 = np.array("""
NA 136 NA 140 NA NA 133 155 NA 157 128 NA NA 152 159 164 153 NA 142
150 158 151 NA NA NA NA NA 134 NA NA 141 NA 155 NA NA NA NA NA 164 167
150 131 NA 158 NA 145 NA NA NA 139 138 NA NA 149 157 156 NA 141 152 NA
141 NA 128 123 162 153 151 173 135 133 NA 127 145 116 151 133 142 148
144 144 NA NA 151 141 147 NA NA 136 136 NA 146 143 132 152 151 153 NA
157 150 141
""".replace("NA", "nan").split(), dtype=float)
d = pd.DataFrame(dict(trt=trt, base=base, wk6=wk6, wk12=wk12))
print("Missing outcomes:", int(d.wk12.isna().sum()), "of 100")

# 1. Complete-case analysis: drop rows with a missing outcome
cc = sm.OLS.from_formula("wk12 ~ trt + base", d, missing="drop").fit()
print(cc.summary2().tables[1].loc["trt", ["Coef.", "Std.Err.", "P>|t|"]].round(3))

# 2. MICE: 20 imputed datasets (predictive mean matching; k_pmm=5 donors)
imp = MICEData(d, k_pmm=5, rng=1)
mi = MICE("wk12 ~ trt + base", sm.OLS, imp).fit(n_burnin=10, n_imputations=20)
print(mi.summary())
```

**What the code does**

- The same 100 patients are typed in; `"NA"` is converted to `nan` so pandas treats it as missing.
- `OLS.from_formula(..., missing="drop")` is the complete-case analysis.
- `MICEData(d, k_pmm=5, rng=1)` sets up chained-equation imputation using **predictive mean matching** (each missing value is replaced by an observed value from a similar patient). It is a different algorithm from R's `method = "norm"` above.
- `MICE(...).fit(n_burnin=10, n_imputations=20)` runs burn-in cycles, then fits the model to 20 imputed datasets and pools them with Rubin's rules. `rng=1` seeds the imputation (statsmodels 0.14+), so reruns reproduce exactly.

**Output**

```text
Missing outcomes: 38 of 100
Coef.      -5.056
Std.Err.    1.831
P>|t|       0.008
Name: trt, dtype: float64
                         Results: MICE
================================================================
Method:                    MICE      Sample size:          100  
Model:                     OLS       Scale                 50.52
Dependent variable:        wk12      Num. imputations      20   
----------------------------------------------------------------
           Coef.  Std.Err.    t    P>|t|   [0.025  0.975]  FMI  
----------------------------------------------------------------
Intercept 23.0582  13.0765  1.7633 0.0778 -2.5713 48.6877 0.2558
trt       -6.2610   1.6076 -3.8948 0.0001 -9.4118 -3.1103 0.2155
base       0.8026   0.0819  9.8026 0.0000  0.6421  0.9631 0.2760
================================================================
```

**The R and Python MI results do not match, and that is expected.** MI is a random procedure and the two programs use different imputation algorithms (Bayesian normal regression in `mice(method = "norm")`, predictive mean matching in statsmodels), different random number streams, and different pooling details. Compare them qualitatively: both give an estimate **further from zero than complete-case**, both have SEs similar to or smaller than complete-case, and both have confidence intervals that include the true effect. Do not expect digits to agree; do expect conclusions to agree. Here they happen to land close together (−6.5 vs −6.3 mmHg), but a different seed or dataset would move each by a few tenths of a mmHg.

## Interpreting the output

- **Complete-case:** treatment effect −5.06 mmHg (SE 1.83, p = 0.008). Only 62 patients contribute, and because dropout depended on week-6 SBP and differed by arm, this estimate is pulled toward zero.
- **MI in R (m = 20):** −6.51 mmHg (SE 1.61, 95% CI −9.73 to −3.29). Closer to the −8 mmHg truth, with a narrower SE than complete-case because the 38 incomplete patients still contribute their baseline and week-6 information.
- **Rubin's rules, checked by hand.** Between-imputation variance b = 0.654, within ū = 1.899. Total t = 1.899 + (1 + 1/20) × 0.654 = 2.585, so SE = √2.585 = 1.61. ✔
- **FMI = 0.29.** About 29% of the uncertainty about the treatment effect is due to the missing outcomes.
- **MI in Python:** −6.26 mmHg (SE 1.61, 95% CI −9.41 to −3.11). Same story as R: further from zero than complete-case and near the full-data estimate, but not identical digits, as explained above.
- **Sensitivity analyses.** Because MAR is untestable, always add a **delta-adjustment** (tipping-point) analysis: impute under MAR, then shift imputed values in the dropout group by δ mmHg (for example, assume dropouts in the treatment arm are 5 mmHg worse than MAR predicts) and see how large δ must be before the conclusion changes. Reference-based imputation ("jump to reference") is another common regulatory choice. Report these as pre-specified sensitivity analyses.

**How to report it:**
> Week-12 SBP was missing for 38 of 100 patients. Assuming data were missing at random given baseline and week-6 SBP, we imputed missing outcomes using chained equations (m = 20) and pooled the estimates with Rubin's rules. The adjusted treatment effect was −6.5 mmHg (95% CI −9.7 to −3.3). A complete-case analysis gave −5.1 mmHg (95% CI −8.7 to −1.4). Sensitivity analyses under MNAR (delta adjustment) are reported in the supplement.

## Check your understanding

<details class="check"><summary>1. "Multiple imputation invents data, so it makes the results look more certain than they are." True or false?</summary>

False. MI does not add information; the imputed values are not treated as real observations. Because the m completed datasets differ, the between-imputation variance B is added to the standard error (T = Ū + (1 + 1/m)B), so MI **widens** uncertainty relative to filling in each gap once. The mistake of overconfidence belongs to *single* imputation (mean or last value carried forward), not to MI done properly.
</details>

<details class="check"><summary>2. Patients with severe symptoms are more likely to miss the follow-up visit, and symptom severity was recorded at baseline. Is this MCAR, MAR, or MNAR? What if severity was not recorded?</summary>

With severity recorded: **MAR**, because missingness depends on an observed variable, and MI that includes severity is valid. If severity was never recorded and still drives dropout, the data are **MNAR** with respect to your dataset: no analysis of the observed data can fully correct it, and you need sensitivity analyses.
</details>

<details class="check"><summary>3. A colleague says complete-case analysis is "conservative" because it uses fewer patients. Using this lesson's numbers, what is wrong?</summary>

Fewer patients makes the CI wider, but it does not guarantee it is centred in the right place. Here the complete-case effect (−5.06 mmHg) is biased toward zero because dropout was related to week-6 SBP and differed between arms. Losing data can lower power **and** introduce bias. In a non-inferiority trial, such bias toward "no difference" can even favour the new drug.
</details>

<details class="check"><summary>4. We use m = 20 imputations. If the between-imputation variance had been 0, what would the total variance and FMI be?</summary>

T = Ū + (1 + 1/20) × 0 = Ū (1.899), so the SE would be the same as in a complete dataset and FMI would be 0. That would mean the missing values are perfectly predictable from the observed data, which almost never happens. A large B is the signal that the missing values matter.
</details>

## Common mistakes

- **Mean imputation or last observation carried forward** as a single fill-in. This shrinks variance, distorts correlations, and gives overconfident p-values.
- **Leaving the outcome or treatment out of the imputation model.** This biases the treatment effect toward zero. The imputation model must contain everything in the analysis model.
- **Imputing, then analysing the average dataset.** Analyse each dataset separately and pool; do not average the data first.
- **Too few imputations** (m = 3–5) for modest to heavy missingness, giving noisy standard errors.
- **Claiming MAR has been "tested".** It is an assumption. Justify it and run MNAR sensitivity analyses.
- **Ignoring the estimand.** Handling dropout differently from what the estimand implies answers a different clinical question.

## Related tests

- [Linear mixed-effects models (MMRM)](#/test/mixed-models): likelihood-based analysis that is valid under MAR without imputation
- [ANCOVA](#/test/ancova): the usual baseline-adjusted analysis model that MI feeds
- [Generalized estimating equations](#/test/gee): requires MCAR unless weighted
- [Propensity scores](#/test/propensity-scores): covariate adjustment for confounding, a different problem often combined with MI
- [Linear regression](#/test/linear-regression): the analysis model used here

## Key takeaways

- **MCAR** is harmless but rare, **MAR** can be handled by MI or mixed models, and **MNAR** needs sensitivity analyses.
- Complete-case analysis can be biased, not just less powerful, when dropout depends on observed outcomes.
- MI = impute m times, analyse each, pool with Rubin's rules: T = Ū + (1 + 1/m)B.
- MI does not create information; it carries uncertainty about the missing values into the standard error.
- Define the estimand (ICH E9(R1)) first, then choose the missing-data method and pre-specify sensitivity analyses.

## References

- Rubin DB. *Multiple Imputation for Nonresponse in Surveys.* Wiley; 1987.
- van Buuren S, Groothuis-Oudshoorn K. mice: Multivariate imputation by chained equations in R. *J Stat Softw.* 2011;45(3).
- ICH E9(R1). *Addendum on estimands and sensitivity analysis in clinical trials.* 2019.
- White IR, Royston P, Wood AM. Multiple imputation using chained equations. *Statistics in Medicine.* 2011;30:377–399.
