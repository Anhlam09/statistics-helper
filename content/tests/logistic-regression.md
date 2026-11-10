# Logistic regression

> Models the probability of a binary outcome (complication yes/no, responder yes/no) as a function of predictors, and reports each predictor's effect as an **odds ratio** with a confidence interval, adjusted for the others.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Binary outcome with one or more predictors, adjusted effects wanted | Binary outcome; continuous or categorical predictors | Independent patients; linear in the logit; enough events; no separation | Post-operative complication vs age, diabetes and smoking |

## Learning objectives

- Explain log-odds, odds ratios, and why exp(β) is the effect measure
- Fit a logistic model in R and Python and read ORs with 95% CIs, the Wald z, and the LR test
- Check the main assumptions (events per variable, linearity in the logit, collinearity, influence)
- Report ORs correctly: as odds ratios, with absolute risks alongside

## What it is

Logistic regression models the **log-odds** of the outcome as a linear function of the predictors: log(p / (1 − p)) = β₀ + β₁X₁ + β₂X₂ + … Coefficients are estimated by maximum likelihood. Exponentiating a coefficient gives an **odds ratio (OR) = exp(β)**: the multiplicative change in the odds of the outcome for a one-unit increase in the predictor (or for the exposed versus unexposed group), holding the other predictors constant.

- **Test statistic:** Wald z = β̂ / SE(β̂) for each coefficient; the likelihood-ratio (LR) χ² (null deviance − residual deviance, df = number of predictors) for the whole model
- **Effect size:** **OR = exp(β) with 95% CI** (primary), preferably for a clinically meaningful unit (per 10 years, per 1 SD). For discrimination, report the C-statistic/[AUC](#/test/roc-auc); for calibration, a calibration plot or Hosmer–Lemeshow test.
- Odds ratios exaggerate risk ratios when the outcome is common (above about 10%); describe them as odds, not risks.

## The intuition

Linear regression models a continuous outcome directly. A probability cannot be modelled that way (it would eventually go below 0 or above 1), so logistic regression models the **log-odds** instead. Each unit of a predictor adds a fixed amount to the log-odds, which is the same as multiplying the odds by a fixed factor: the odds ratio. Think of odds as the bookmaker's language: odds of 1 (50:50) become 2 (2:1 on) when multiplied by an OR of 2.

The signal for each predictor is its coefficient β̂ (the log-OR); the noise is its standard error. Their ratio, the Wald z, is read against a standard normal curve. For age, β̂ = 0.0772 and SE = 0.0285, so z = 2.71.

The dot plot shows the raw reason age matters here: patients who had a complication tend to be older (mean 71.5 vs 63.2 years), although the ranges overlap a lot. Logistic regression turns "older tends to go with complications" into a number, OR = 1.08 per year, after accounting for diabetes and smoking.

<div class="viz" data-viz="dots" data-unit="age, years" data-groups='{"No complication":[69,53,55,78,65,79,63,84,51,59,47,71,82,65,55,65,50,66,50,57,75,55,59,58,45,76,83,55,70,56],"Complication":[84,71,60,85,78,72,49,71,73,79,76,77,77,68,49,51,64,71,66,82,63,82,62,75,79,72,74,81,76,78]}' data-caption="Notice the complication group is shifted toward older ages (mean 71.5 vs 63.2), but there is heavy overlap: age is a risk factor, not a deterministic rule, which is why the OR per year is modest (1.08) and the CI is not narrow."></div>

If age were unrelated to complications (β = 0), the Wald z would follow the curve below. Our z = 2.71 is in the tail, giving p = 0.007.

<div class="viz" data-viz="z" data-stat="2.71"></div>

An OR of 1.08 per year sounds small, but ORs compound: 1.08¹⁰ = 2.16 per decade, which is why the lesson also reports the OR per 10 years.

## When to use it

Use logistic regression when:

- The outcome is **binary** (30-day complication, remission, death by 1 year with complete follow-up)
- You want **adjusted** odds ratios, or predicted probabilities for individual patients
- Observations are independent
- You have enough events: a traditional rough guide is **about 10 events per variable (EPV)**, counting the *smaller* outcome category (it is a guide, not a law; more complex or less certain models need more)

**Clinical scenario.** A surgical team reviews 60 patients over 45 years old after elective surgery. Thirty developed a post-operative complication. They want to know whether age, diabetes, and current smoking independently predict complications. *What are the adjusted odds ratios?*

**Don't use it when:**

- The outcome is continuous → use [linear regression](#/test/linear-regression)
- The outcome is a count or event rate with varying follow-up → use [Poisson or negative binomial regression](#/test/poisson-regression)
- Follow-up times differ and the time to the event matters → use [Cox regression](#/test/cox-ph) or the [log-rank test](#/test/log-rank)
- Observations are clustered or repeated within patients → use [GEE](#/test/gee) or a [mixed model](#/test/mixed-models)
- There is a single binary predictor and no adjustment → a [chi-square test](#/test/chi-square-independence) or [Fisher's exact test](#/test/fisher-exact) is enough
- There are very few events or complete separation (a predictor perfectly predicts the outcome) → consider Firth's penalised likelihood

## Assumptions

1. **Binary outcome and independent observations.** One row per patient; no repeated measures.
2. **Linearity in the logit.** Each continuous predictor must be linearly related to the *log-odds*. Check: add a Box–Tidwell term (x · log x), or compare with a spline fit; plot binned log-odds.
3. **Sufficient events per variable (EPV about 10 or more).** With too few events, coefficients are biased away from 0 and CIs are unreliable. EPV of 10 is a traditional rough guide rather than a hard threshold, and newer sample-size approaches also consider the expected shrinkage. Here EPV = 30 / 3 = 10, at the lower limit; count each dummy variable as a predictor.
4. **No multicollinearity.** Check VIFs (above about 5 is worrying) or the correlation matrix of predictors.
5. **No influential outliers.** Check Cook's distance, Pearson/deviance residuals, and leverage; refit without flagged patients as a sensitivity analysis.
6. **No complete or quasi-complete separation.** Very large coefficients and standard errors signal a problem.
7. **Adequate calibration and a sensible model.** Use a calibration plot or Hosmer–Lemeshow test (`lackfit` in SAS); pre-specify predictors rather than selecting them by p-values.

## Hypotheses

For each predictor j:

- **H₀:** βⱼ = 0, i.e. OR = 1 (no association with the outcome, after adjustment)
- **H₁:** βⱼ ≠ 0, i.e. OR ≠ 1 (two-sided)

Overall model: H₀: all slope coefficients are 0 (LR test).

## Example data

Simulated data for 60 patients (age 45 to 85 years; about 40% with diabetes and 35% smokers). The first five rows:

| Patient | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| Age (years) | 69 | 53 | 84 | 71 | 60 |
| Diabetes (1 = yes) | 0 | 1 | 0 | 1 | 1 |
| Smoker (1 = yes) | 0 | 0 | 0 | 0 | 0 |
| Complication (1 = yes) | 0 | 0 | 1 | 1 | 1 |

The full 60 rows are typed inline in each code block below.

## How to run it

### R

```r
# Post-operative complication (1 = yes) in 60 patients: age, diabetes, smoking
dat <- data.frame(
  age = c(69, 53, 84, 71, 60, 85, 78, 55, 72, 78, 49, 71, 65, 73, 79, 79, 63,
      84, 51, 76, 77, 77, 59, 47, 68, 71, 49, 51, 82, 64, 65, 55, 65, 50, 66,
      71, 50, 66, 57, 75, 82, 63, 82, 55, 59, 62, 58, 75, 79, 45, 72, 74, 76,
      83, 55, 81, 70, 76, 78, 56),
  diabetes = c(0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1,
      1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0,
      0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0),
  smoker = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1,
      0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0),
  event = c(0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1,
      0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1,
      0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0)
)

# 1. Events per variable: 30 events / 3 predictors = 10
sum(dat$event); sum(dat$event) / 3

# 2. Fit the logistic model
fit <- glm(event ~ age + diabetes + smoker, data = dat, family = binomial)
summary(fit)$coefficients

# 3. Odds ratios with Wald 95% CIs (exp of the log-odds scale)
round(exp(cbind(OR = coef(fit), confint.default(fit))), 3)
exp(10 * coef(fit)["age"])                 # OR per 10 years of age

# 4. Overall LR test and diagnostics
lr <- fit$null.deviance - fit$deviance       # likelihood-ratio chi-square
c(chisq = lr, df = 3, p = pchisq(lr, df = 3, lower.tail = FALSE))
car::vif(fit)                              # collinearity
max(cooks.distance(fit))                   # influential observations
# Box-Tidwell: is the age effect linear on the logit scale?
summary(glm(event ~ age + I(age * log(age)) + diabetes + smoker,
            data = dat, family = binomial))$coefficients[3, ]
```

**What the code does**

- `glm(..., family = binomial)` fits a logistic model with the logit link (the default for `binomial`). The outcome must be 0/1 (or a factor whose *second* level is the event).
- `confint.default(fit)` gives **Wald** intervals on the log-odds scale; `exp()` converts coefficients and CIs to odds ratios. `confint(fit)` would give profile-likelihood intervals, which are slightly better in small samples and are the default in recent R versions.
- `exp(10 * coef(fit)["age"])` re-expresses the age effect per 10 years, which is more interpretable than per year.
- The LR χ² compares the model with the intercept-only model via the deviance drop.
- `car::vif()` gives VIFs; `cooks.distance()` flags influential observations.
- The final `glm` adds the Box–Tidwell term `age·log(age)`; its p-value (row 3 of the coefficient table) tests non-linearity on the logit scale.

**Output**

```text
[1] 30
[1] 10
               Estimate Std. Error    z value    Pr(>|z|)
(Intercept) -6.21610110 2.08957431 -2.9748170 0.002931634
age          0.07722116 0.02849085  2.7103844 0.006720528
diabetes     1.32916839 0.60397634  2.2006961 0.027757544
smoker       0.63675946 0.64316028  0.9900479 0.322150728
               OR 2.5 % 97.5 %
(Intercept) 0.002 0.000  0.120
age         1.080 1.022  1.142
diabetes    3.778 1.156 12.341
smoker      1.890 0.536  6.668
     age 
2.164548 
       chisq           df            p 
14.455685390  3.000000000  0.002346155 
     age diabetes   smoker 
1.067055 1.023490 1.071640 
[1] 0.1558613
   Estimate  Std. Error     z value    Pr(>|z|) 
-0.02337341  0.35672477 -0.06552226  0.94775818
```

### Python

```python
import numpy as np, pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import variance_inflation_factor as vif

# Post-operative complication (1 = yes) in 60 patients: age, diabetes, smoking
dat = pd.DataFrame({
  'age': [69, 53, 84, 71, 60, 85, 78, 55, 72, 78, 49, 71, 65, 73, 79, 79, 63,
      84, 51, 76, 77, 77, 59, 47, 68, 71, 49, 51, 82, 64, 65, 55, 65, 50, 66,
      71, 50, 66, 57, 75, 82, 63, 82, 55, 59, 62, 58, 75, 79, 45, 72, 74, 76,
      83, 55, 81, 70, 76, 78, 56],
  'diabetes': [0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1,
      1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0,
      0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0],
  'smoker': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1,
      0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0],
  'event': [0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1,
      0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1,
      0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0]
})

# 1. Events per variable
print("events =", dat.event.sum(), " EPV =", dat.event.sum() / 3)

# 2. Fit the logistic model
fit = smf.logit("event ~ age + diabetes + smoker", data=dat).fit(disp=False)
print(fit.summary2().tables[1].round(4))

# 3. Odds ratios with Wald 95% CIs
ci = fit.conf_int()
print(pd.DataFrame({"OR": np.exp(fit.params), "2.5%": np.exp(ci[0]),
                    "97.5%": np.exp(ci[1])}).round(3))
print("OR per 10 years of age: %.3f" % np.exp(10 * fit.params["age"]))

# 4. Overall LR test and diagnostics
print("LR chi2 = %.2f, df = %d, p = %.5f" % (fit.llr, fit.df_model, fit.llr_pvalue))
X = fit.model.exog
print("VIF:", np.round([vif(X, i) for i in range(1, X.shape[1])], 2).tolist())
print("max Cook's D = %.3f" % fit.get_influence().cooks_distance[0].max())
bt = smf.logit("event ~ age + I(age * np.log(age)) + diabetes + smoker",
               data=dat).fit(disp=False)
print("Box-Tidwell term p = %.3f" % bt.pvalues.iloc[2])
```

**What the code does**

- `smf.logit(...)` fits logistic regression by maximum likelihood; `disp=False` hides the optimiser output. (Alternatively, `smf.glm(..., family=sm.families.Binomial())`.) The outcome must be 0/1.
- `fit.conf_int()` returns Wald CIs on the log-odds scale; exponentiating gives the OR CIs, which agree with R's `confint.default`.
- `fit.llr` and `fit.llr_pvalue` are the LR χ² and its p-value against the intercept-only model.
- VIFs are computed on the design matrix, skipping the intercept column.
- `get_influence().cooks_distance` gives Cook's distance.
- `I(age * np.log(age))` in the formula is the Box–Tidwell term.

**Output**

```text
events = 30  EPV = 10.0
            Coef.  Std.Err.       z   P>|z|   [0.025  0.975]
Intercept -6.2161    2.0897 -2.9747  0.0029 -10.3118 -2.1204
age        0.0772    0.0285  2.7103  0.0067   0.0214  0.1331
diabetes   1.3292    0.6040  2.2006  0.0278   0.1454  2.5130
smoker     0.6368    0.6432  0.9900  0.3222  -0.6238  1.8974
              OR   2.5%   97.5%
Intercept  0.002  0.000   0.120
age        1.080  1.022   1.142
diabetes   3.778  1.156  12.342
smoker     1.890  0.536   6.668
OR per 10 years of age: 2.165
LR chi2 = 14.46, df = 3, p = 0.00235
VIF: [1.04, 1.02, 1.04]
max Cook's D = 0.156
Box-Tidwell term p = 0.948
```

### SAS

```sas
/* Post-operative complication (event=1) in 60 patients */
data surg;
  input age diabetes smoker event;
  datalines;
  69 0 0 0
  53 1 0 0
  84 0 0 1
  71 1 0 1
  60 1 0 1
  85 1 0 1
  78 1 0 1
  55 0 0 0
  72 1 0 1
  78 1 0 0
  49 1 0 1
  71 1 0 1
  65 0 0 0
  73 1 1 1
  79 1 1 0
  79 1 0 1
  63 0 0 0
  84 0 0 0
  51 0 1 0
  76 1 0 1
  77 1 0 1
  77 1 0 1
  59 1 0 0
  47 1 0 0
  68 1 0 1
  71 1 0 0
  49 0 0 1
  51 1 1 1
  82 0 0 0
  64 1 1 1
  65 1 0 0
  55 0 1 0
  65 0 0 0
  50 0 1 0
  66 1 1 0
  71 0 1 1
  50 1 0 0
  66 0 1 1
  57 1 1 0
  75 0 0 0
  82 1 0 1
  63 1 1 1
  82 1 1 1
  55 0 1 0
  59 0 0 0
  62 0 1 1
  58 1 0 0
  75 1 0 1
  79 1 0 1
  45 0 1 0
  72 0 0 1
  74 1 1 1
  76 0 0 0
  83 0 1 0
  55 1 1 0
  81 1 1 1
  70 1 0 0
  76 0 1 1
  78 0 0 1
  56 0 0 0
;
run;

/* Logistic regression: model P(event=1); odds ratios with Wald CIs */
proc logistic data=surg;
  model event(event='1') = age diabetes smoker / clodds=wald lackfit rsquare
                                                 influence iplots;
  units age = 10;                  /* extra OR per 10 years */
run;

/* Collinearity: VIF from PROC REG on the same predictors */
proc reg data=surg;
  model event = age diabetes smoker / vif;
run;
quit;
```

**What the code does**

- `event(event='1')` makes SAS model the probability that `event = 1`. **Without it (or `descending` on the PROC statement), SAS models the probability of the lowest sorted value (0), and every OR is inverted**, a very common error. Check the "Response Profile" table in the output.
- `clodds=wald` prints odds ratios with Wald 95% CIs; `units age = 10;` adds the OR per 10 years.
- `lackfit` gives the Hosmer–Lemeshow test and `rsquare` the generalised R²; the C-statistic is printed in the "Association of Predicted Probabilities" table by default.
- `influence iplots` outputs case-wise diagnostics (leverage, DFBETAs, Cook-type C and CBAR) and plots.
- `proc reg ... / vif` uses a linear model purely to obtain collinearity diagnostics; VIFs depend only on the predictors, so they remain valid for the logistic model.
- For categorical predictors with more than two levels, add `class sex (ref='F') / param=ref;` so the odds ratios compare each level with the reference level.

## Interpreting the output

- **Events per variable:** 30 events and 3 predictors give EPV = 10, adequate but only just. No more predictors should be added.
- **Age: OR = 1.08 per year (95% CI 1.02 to 1.14, p = 0.007).** Each additional year multiplies the odds of a complication by 1.08, holding diabetes and smoking constant. Per 10 years, OR = 2.16 (95% CI 1.24 to 3.78), i.e. the odds roughly double per decade.
- **Diabetes: OR = 3.78 (95% CI 1.16 to 12.34, p = 0.028).** Diabetic patients have almost four times the odds of a complication, but the CI is very wide (small sample): the data are compatible with a 16% to 12-fold increase.
- **Smoker: OR = 1.89 (95% CI 0.54 to 6.67, p = 0.32).** The point estimate suggests increased odds, but the interval includes 1. This is inconclusive, not evidence of no effect.
- **Intercept: OR = 0.002.** The odds of a complication for a 0-year-old non-diabetic non-smoker; not interpretable. Centre age to make it meaningful.
- **LR χ²(3) = 14.46, p = 0.002.** The model as a whole is better than the intercept-only model.
- **VIF = 1.04, 1.02, 1.04.** No collinearity. **Max Cook's D = 0.16.** No observation dominates the fit. **Box–Tidwell p = 0.95.** No evidence that the age effect is non-linear on the logit scale.
- In **SAS**, check the Hosmer–Lemeshow (`lackfit`) p-value for calibration (a large p-value, e.g. > 0.05, shows no evidence of poor fit); this test has low power for n = 60.

**How to report it:**
> In multivariable logistic regression (n = 60, 30 events), the odds of a post-operative complication increased with age (adjusted OR 1.08 per year, 95% CI 1.02 to 1.14; 2.16 per 10 years, 95% CI 1.24 to 3.78) and with diabetes (adjusted OR 3.78, 95% CI 1.16 to 12.34). Smoking was not clearly associated (OR 1.89, 95% CI 0.54 to 6.67).

## Check your understanding

<details class="check"><summary>1. A researcher wants to predict 30-day readmission (yes/no) using 12 candidate predictors in a cohort of 90 patients, of whom 20 were readmitted. What is the concern?</summary>

Only 20 events for 12 predictors is under 2 events per variable. The model will overfit: ORs will be inflated and predictions will fail in new patients. The traditional guide of about 10 events per variable is a rough guide, not a law (newer work bases the required size on the expected shrinkage and the anticipated model performance), but 20 events for 12 predictors is far outside any sensible range. Cut the predictors, use penalisation, or collect more data.
</details>

<details class="check"><summary>2. In the output the smoking OR is 1.89 (95% CI 0.54 to 6.67, p = 0.32). A colleague concludes that smoking does not affect complications. Respond.</summary>

p = 0.32 is not evidence of no effect. The CI runs from a 46% reduction in odds to a 6.7-fold increase, so the data are compatible with a harmful effect that this small study simply could not detect. Say "smoking was not clearly associated", and keep the CI in the report.
</details>

<details class="check"><summary>3. Diabetes has an adjusted OR of 3.78. Can you say that diabetic patients are 3.78 times as likely to have a complication?</summary>

No. That is the odds ratio, not the risk ratio. The outcome is common here (50% overall), so the OR exaggerates the RR. To describe risk, compute predicted probabilities for the two groups at the same age and smoking status (for example, about 23% vs 53% for a non-smoker aged 65), or report the risk difference.
</details>

<details class="check"><summary>4. What would happen to the age OR if age were measured in decades instead of years?</summary>

The model is the same but the coefficient is multiplied by 10, so the OR becomes exp(10 × 0.0772) = 2.16 (per decade), with CI 1.24 to 3.78. The p-value and the Wald z do not change. Choose a unit that is clinically meaningful.
</details>

## Common mistakes

- **Reading odds ratios as risk ratios.** With a common outcome (here 50%), the OR is much larger than the RR. Report absolute risks (predicted probabilities) as well.
- **Too few events for the number of predictors.** Overfitted models give inflated ORs and poor prediction in new patients. Count events, not patients.
- **Wrong event coding.** In SAS, forgetting `event='1'` (or `descending`); in R, an outcome factor whose levels are in the wrong order.
- **Interpreting a non-significant OR as "no effect".** Look at the CI width.
- **Dichotomising continuous predictors** (age ≥ 65) to simplify the model. It discards information.
- **Variable selection by p-value** (stepwise). It biases ORs and invalidates CIs; pre-specify the model.
- **Ignoring separation.** Huge ORs with huge standard errors indicate sparse cells.
- **Assessing fit only by p-values.** Also evaluate discrimination ([AUC](#/test/roc-auc)) and calibration, preferably in a validation sample.

## Related tests

- [ROC curve, AUC, sensitivity & specificity](#/test/roc-auc): discrimination of the fitted model
- [Chi-square test of independence](#/test/chi-square-independence) and [Fisher's exact test](#/test/fisher-exact): unadjusted association of binary variables
- [Cochran–Mantel–Haenszel test](#/test/cmh): stratified 2 × 2 tables
- [Cox proportional hazards regression](#/test/cox-ph): time-to-event outcome
- [GEE](#/test/gee): repeated or clustered binary outcomes
- [Linear regression](#/test/linear-regression): continuous outcome

## Key takeaways

- Logistic regression models the log-odds. exp(β) is the odds ratio for a one-unit increase, adjusted for the other predictors.
- An OR is not a risk ratio. For common outcomes (more than about 10%) it overstates the RR, so also report predicted probabilities.
- About 10 events per variable is only a rough guide. Count events (the smaller outcome group) and avoid overfitting.
- A wide CI that includes 1 means "inconclusive", not "no effect".
- Check linearity in the logit, collinearity, influence, discrimination (AUC) and calibration.

## References

- Hosmer DW, Lemeshow S, Sturdivant RX. *Applied Logistic Regression.* 3rd ed. Wiley; 2013.
- Peduzzi P, Concato J, Kemper E, Holford TR, Feinstein AR. A simulation study of the number of events per variable in logistic regression analysis. *J Clin Epidemiol.* 1996;49:1373–1379.
- R documentation: `?glm`, `?confint.default`; statsmodels: `statsmodels.formula.api.logit`
- SAS/STAT User's Guide: The LOGISTIC Procedure (EVENT= and DESCENDING options, CLODDS=, LACKFIT)
