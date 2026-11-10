# Linear regression

> Models a continuous outcome as a straight-line function of one or more predictors, giving each predictor's effect in the outcome's own units (mmHg, mg/dL) with a confidence interval, adjusted for the other predictors.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Continuous outcome, one or more predictors; adjusted effects or prediction | Continuous outcome; continuous or binary predictors | Linear relationship; independent errors; constant variance; roughly normal residuals | SBP vs age, BMI, and sex in 50 adults |

## Learning objectives

- Explain a regression coefficient as an adjusted change in the outcome per unit of a predictor
- Fit a multiple linear regression in R and Python and read coefficients, CIs, R², and residual SE
- Check the assumptions using residual plots and influence diagnostics, not the raw outcome
- Report adjusted effects per clinically meaningful unit without implying causation

## What it is

Multiple linear regression fits the equation Y = β₀ + β₁X₁ + β₂X₂ + … + ε by ordinary least squares (OLS), which chooses the coefficients that minimise the sum of squared residuals. Each coefficient β is the expected change in Y for a one-unit increase in that predictor **holding the other predictors constant**. The residual ε is the part of Y the model does not explain, assumed to be random noise with constant variance.

- **Test statistic:** for each coefficient, t = β̂ / SE(β̂) with n − p − 1 degrees of freedom (p = number of predictors); for the whole model, F = (explained variance / p) / (residual variance / (n − p − 1))
- **Effect size:** the **coefficients with 95% CIs** on the clinical scale (primary). Report R² (adjusted R² for model comparison) and the residual standard error (typical prediction error in the outcome's units) as model-level summaries.

## The intuition

Regression is a **signal-to-noise** story at the level of each coefficient. The signal is the slope: how many mmHg SBP rises per year of age. The noise is the standard error of that slope, which shrinks when the predictor varies widely, when the scatter around the line is small, and when n is large. The t statistic is slope ÷ standard error. For age it is 0.742 / 0.092 = 8.04, so the slope is eight standard errors away from zero.

Analogy: a growth chart. A child's height rises roughly steadily with age, but individual children scatter above and below the average line. Regression draws the average line and tells you how steep it is; the residual standard error (7.5 mmHg here) is the typical distance of a patient from the line.

Notice in the plot how SBP climbs across age bands (mean about 127, 135 and 149 mmHg), but patients within a band still differ by 20 mmHg or more. Regression uses all that information at once and adjusts for BMI and sex, which the plot ignores.

<div class="viz" data-viz="dots" data-unit="systolic BP, mmHg" data-groups='{"Under 45":[114,129,125,123,128,131,115,128,125,145,131,123,133,135,123],"45 to 59":[135,116,139,127,135,130,153,130,138,132,136,139,143,126,142,153,131,123,136,135,134,141],"60 and over":[156,151,148,157,149,142,158,145,150,136,155,163,127]}' data-caption="Unadjusted SBP by age band: a clear upward trend, but wide scatter within each band."></div>

The slope for age has t = 8.04 with 46 residual degrees of freedom. If age had no linear relationship with SBP, a t this large would essentially never occur.

<div class="viz" data-viz="t" data-stat="8.04" data-df="46"></div>

**What to notice:** the age statistic sits far out in the tail, so the shaded area is invisible. Compare with the male coefficient (t = 1.92, p = 0.062), which falls just short of the usual 0.05 cut-off: the point estimate is 4.2 mmHg, but its standard error (2.2 mmHg) is about half of it.

## When to use it

Use linear regression when:

- The outcome is **continuous** (blood pressure, LDL cholesterol, change in HbA1c, FEV₁)
- You want an **adjusted** effect (e.g., of BMI after accounting for age and sex), a **prediction**, or a **dose–response slope**
- Observations are independent, and the relationship with each predictor is approximately linear
- You have enough subjects: a common rule of thumb is at least 10–20 observations per predictor

**Clinical scenario.** A hypertension clinic records systolic blood pressure (SBP, mmHg), age (years), BMI (kg/m²) and sex in 50 adults. *By how many mmHg does SBP differ per additional year of age and per 1 kg/m² of BMI, independent of sex?*

**Don't use it when:**

- The outcome is **binary** → use [logistic regression](#/test/logistic-regression)
- The outcome is a **count** or a rate → use [Poisson or negative binomial regression](#/test/poisson-regression)
- The outcome is time to an event → use [Cox proportional hazards regression](#/test/cox-ph)
- The same patients are measured repeatedly → use a [mixed model](#/test/mixed-models) or [GEE](#/test/gee)
- You only need the strength of association between two variables → use [Pearson](#/test/pearson) or [Spearman](#/test/spearman) correlation
- You are comparing treatment groups with a baseline covariate in a trial → [ANCOVA](#/test/ancova) is the same model in specialised form

## Assumptions

1. **Linearity.** The mean of Y is a straight-line function of each predictor. Check: residuals-versus-fitted and residuals-versus-predictor plots should show no curve; add splines or polynomial terms if they do.
2. **Independent errors.** Observations (and residuals) are independent; no clustering or repeated measures. Check by design; for time-ordered data, a Durbin–Watson test or residual autocorrelation plot.
3. **Constant variance (homoscedasticity).** The residual spread is the same across fitted values. Check: scale–location plot, Breusch–Pagan test. If violated, use robust (sandwich) standard errors or transform the outcome.
4. **Approximately normal residuals.** Needed for the small-sample t-tests and CIs, not for unbiased coefficients. Check: Q–Q plot of residuals, [Shapiro–Wilk test](#/test/shapiro-wilk). With n ≥ 30–50 moderate departures matter little.
5. **No problematic multicollinearity.** Highly correlated predictors inflate standard errors. Check: variance inflation factors (VIF); values above about 5 (some use 10) deserve attention.
6. **No unduly influential points.** A few observations should not drive the fit. Check: Cook's distance (compare with 4/n), studentized residuals beyond ±3, and leverage; refit without flagged points as a sensitivity analysis, and do not delete them silently.
7. **Correct model specification.** All important confounders are included; the model is not overfitted. Variable selection should be pre-specified, not driven by p-values.

## Hypotheses

For each predictor j (the overall F-test has the analogous global form):

- **H₀:** βⱼ = 0 (no linear association with the outcome, after adjusting for the other predictors)
- **H₁:** βⱼ ≠ 0 (two-sided)

The F-test tests H₀: all slope coefficients are 0.

## Example data

Simulated data for 50 adults (age 35–75 years, BMI about 27 kg/m², about half male). The first five rows:

| Patient | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| Age (years) | 46 | 67 | 40 | 41 | 39 |
| BMI (kg/m²) | 27.1 | 33.8 | 23.6 | 22.2 | 29.1 |
| Male (1 = yes) | 0 | 1 | 0 | 1 | 1 |
| SBP (mmHg) | 135 | 156 | 114 | 129 | 125 |

The full 50 rows are typed inline in each code block below.

## How to run it

### R

```r
# Systolic BP (mmHg) vs age (years), BMI (kg/m²) and sex (male = 1), n = 50
dat <- data.frame(
  age = c(46, 67, 40, 41, 39, 57, 70, 55, 45, 70, 55, 73, 46, 54, 36, 42, 74,
      46, 63, 72, 51, 41, 68, 72, 51, 37, 50, 51, 39, 41, 52, 68, 45, 67, 55,
      41, 41, 55, 40, 52, 46, 47, 63, 35, 42, 37, 57, 64, 45, 58),
  bmi = c(27.1, 33.8, 23.6, 22.2, 29.1, 28.3, 25.8, 29.2, 22, 23.8, 34, 30.5,
      29.9, 33.9, 27.6, 30.5, 31.2, 22.5, 26.8, 36.6, 32, 23.9, 15.2, 28.3,
      25.8, 24.1, 28.2, 34, 27, 27, 33.2, 25.5, 24.5, 27.5, 20.1, 31.9, 35.8,
      25.5, 30.7, 25, 27.8, 28.6, 31.6, 32, 27.7, 29.9, 25.5, 29.8, 28, 29.7),
  male = c(0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1,
      1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0,
      0, 0, 1, 1),
  sbp = c(135, 156, 114, 129, 125, 116, 151, 139, 127, 148, 135, 157, 130,
      153, 123, 128, 149, 130, 142, 158, 138, 131, 145, 150, 132, 115, 136,
      139, 128, 125, 143, 136, 126, 155, 142, 145, 131, 153, 123, 131, 123,
      136, 163, 133, 135, 123, 135, 127, 134, 141)
)

# 1. Fit the model; coefficients are in mmHg per unit of each predictor
fit <- lm(sbp ~ age + bmi + male, data = dat)
summary(fit)
confint(fit)

# 2. Diagnostics
library(car); library(lmtest)
vif(fit)                                # multicollinearity
shapiro.test(resid(fit))                # normality of residuals
bptest(fit)                             # constant variance (Breusch-Pagan)
max(cooks.distance(fit))                # influential points (look > 4/n = 0.08)
which(abs(rstudent(fit)) > 3)           # outliers

# 3. Sensitivity analysis: refit without the one outlier (patient 6)
confint(update(fit, data = dat[-6, ]))["age", ]
summary(update(fit, data = dat[-6, ]))$r.squared
# plot(fit, which = 1:4)                # residuals vs fitted, Q-Q, scale-location, Cook
```

**What the code does**

- `lm(sbp ~ age + bmi + male)` fits OLS. `male` is coded 0/1, so its coefficient is the adjusted male-minus-female difference. For a multi-level factor, use `factor()`; R then creates dummy variables against the first level.
- `confint(fit)` gives 95% CIs based on the t-distribution; the `summary()` table shows only standard errors.
- `vif()` (package **car**) reports variance inflation factors.
- `shapiro.test(resid(fit))` and `bptest()` (package **lmtest**) test normality of residuals and constant variance. Tests are only a supplement to the diagnostic plots (the commented `plot(fit, which = 1:4)` line).
- `cooks.distance()` and `rstudent()` flag influential points and outliers. Here one patient (row 6) is flagged.
- The last block refits without that patient to see whether conclusions change.

**Output**

```text

Call:
lm(formula = sbp ~ age + bmi + male, data = dat)

Residuals:
     Min       1Q   Median       3Q      Max 
-25.9628  -3.6533   0.3643   3.2583  14.5503 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept) 75.01302    8.82247   8.503 5.43e-11 ***
age          0.74162    0.09227   8.038 2.60e-10 ***
bmi          0.72496    0.25891   2.800  0.00745 ** 
male         4.16134    2.17139   1.916  0.06153 .  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 7.484 on 46 degrees of freedom
Multiple R-squared:  0.6272,	Adjusted R-squared:  0.6029 
F-statistic:  25.8 on 3 and 46 DF,  p-value: 6.123e-10

                 2.5 %     97.5 %
(Intercept) 57.2543164 92.7717186
age          0.5558912  0.9273399
bmi          0.2037937  1.2461231
male        -0.2094453  8.5321273
     age      bmi     male 
1.002578 1.027350 1.025349 

	Shapiro-Wilk normality test

data:  resid(fit)
W = 0.93106, p-value = 0.006061


	studentized Breusch-Pagan test

data:  fit
BP = 1.8022, df = 3, p-value = 0.6145

[1] 0.1255552
6 
6 
    2.5 %    97.5 % 
0.6004195 0.9215810 
[1] 0.7109612
```

### Python

```python
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.outliers_influence import variance_inflation_factor as vif

# Systolic BP (mmHg) vs age (years), BMI (kg/m²) and sex (male = 1), n = 50
dat = pd.DataFrame({
  'age': [46, 67, 40, 41, 39, 57, 70, 55, 45, 70, 55, 73, 46, 54, 36, 42, 74,
      46, 63, 72, 51, 41, 68, 72, 51, 37, 50, 51, 39, 41, 52, 68, 45, 67, 55,
      41, 41, 55, 40, 52, 46, 47, 63, 35, 42, 37, 57, 64, 45, 58],
  'bmi': [27.1, 33.8, 23.6, 22.2, 29.1, 28.3, 25.8, 29.2, 22, 23.8, 34, 30.5,
      29.9, 33.9, 27.6, 30.5, 31.2, 22.5, 26.8, 36.6, 32, 23.9, 15.2, 28.3,
      25.8, 24.1, 28.2, 34, 27, 27, 33.2, 25.5, 24.5, 27.5, 20.1, 31.9, 35.8,
      25.5, 30.7, 25, 27.8, 28.6, 31.6, 32, 27.7, 29.9, 25.5, 29.8, 28, 29.7],
  'male': [0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1,
      1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0,
      0, 0, 1, 1],
  'sbp': [135, 156, 114, 129, 125, 116, 151, 139, 127, 148, 135, 157, 130,
      153, 123, 128, 149, 130, 142, 158, 138, 131, 145, 150, 132, 115, 136,
      139, 128, 125, 143, 136, 126, 155, 142, 145, 131, 153, 123, 131, 123,
      136, 163, 133, 135, 123, 135, 127, 134, 141]
})

# 1. Fit the model; coefficients are in mmHg per unit of each predictor
fit = smf.ols("sbp ~ age + bmi + male", data=dat).fit()
print(fit.summary().tables[1])
print(f"R2 = {fit.rsquared:.3f}, adj R2 = {fit.rsquared_adj:.3f}, "
      f"sigma = {np.sqrt(fit.scale):.2f}, F = {fit.fvalue:.2f}, p = {fit.f_pvalue:.2e}")

# 2. Diagnostics
X = fit.model.exog                      # design matrix incl. intercept
vifs = [vif(X, i) for i in range(1, X.shape[1])]   # skip the intercept
print("VIF:", dict(zip(fit.model.exog_names[1:], np.round(vifs, 2).tolist())))
print("Shapiro resid p = %.3f" % stats.shapiro(fit.resid).pvalue)
from statsmodels.stats.diagnostic import het_breuschpagan
print("Breusch-Pagan p = %.3f" % het_breuschpagan(fit.resid, X)[1])
infl = fit.get_influence()
print("max Cook's D = %.3f" % infl.cooks_distance[0].max())
print("n |studentized resid| > 3:", int((np.abs(infl.resid_studentized_external) > 3).sum()))

# 3. Sensitivity analysis: refit without the one outlier (patient 6, row index 5)
fit2 = smf.ols("sbp ~ age + bmi + male", data=dat.drop(index=5)).fit()
print("Without patient 6: age CI", fit2.conf_int().loc["age"].round(3).tolist(), "R2 = %.3f" % fit2.rsquared)
```

**What the code does**

- `smf.ols("sbp ~ age + bmi + male", data=dat)` uses an R-style formula and adds the intercept automatically. `.fit()` runs OLS.
- `fit.summary().tables[1]` prints coefficients with SEs, t, p, and 95% CIs; the CI columns match R's `confint()`. Printed p-values of `0.000` mean p < 0.0005.
- `variance_inflation_factor(X, i)` needs the design matrix including the intercept column, so the loop starts at column 1 to skip the intercept's own VIF.
- `het_breuschpagan` returns (LM statistic, LM p-value, F statistic, F p-value); index `[1]` is the LM p-value, which is the same studentized test as `bptest()` in R.
- `get_influence()` gives Cook's distance and externally studentized residuals (R's `rstudent`).
- The final refit drops row index 5 (patient 6; Python counts from 0).

**Output**

```text
==============================================================================
                 coef    std err          t      P>|t|      [0.025      0.975]
------------------------------------------------------------------------------
Intercept     75.0130      8.822      8.503      0.000      57.254      92.772
age            0.7416      0.092      8.038      0.000       0.556       0.927
bmi            0.7250      0.259      2.800      0.007       0.204       1.246
male           4.1613      2.171      1.916      0.062      -0.209       8.532
==============================================================================
R2 = 0.627, adj R2 = 0.603, sigma = 7.48, F = 25.80, p = 6.12e-10
VIF: {'age': 1.0, 'bmi': 1.03, 'male': 1.03}
Shapiro resid p = 0.006
Breusch-Pagan p = 0.614
max Cook's D = 0.126
n |studentized resid| > 3: 1
Without patient 6: age CI [0.6, 0.922] R2 = 0.711
```

### SAS

```sas
/* SBP (mmHg) vs age (years), BMI (kg/m2) and sex (male=1), n = 50 */
data bp;
  input age bmi male sbp;
  datalines;
  46 27.1 0 135
  67 33.8 1 156
  40 23.6 0 114
  41 22.2 1 129
  39 29.1 1 125
  57 28.3 1 116
  70 25.8 1 151
  55 29.2 0 139
  45 22 0 127
  70 23.8 1 148
  55 34 0 135
  73 30.5 0 157
  46 29.9 1 130
  54 33.9 0 153
  36 27.6 0 123
  42 30.5 0 128
  74 31.2 0 149
  46 22.5 1 130
  63 26.8 1 142
  72 36.6 0 158
  51 32 1 138
  41 23.9 1 131
  68 15.2 1 145
  72 28.3 0 150
  51 25.8 1 132
  37 24.1 0 115
  50 28.2 1 136
  51 34 0 139
  39 27 1 128
  41 27 1 125
  52 33.2 1 143
  68 25.5 1 136
  45 24.5 0 126
  67 27.5 1 155
  55 20.1 1 142
  41 31.9 1 145
  41 35.8 1 131
  55 25.5 1 153
  40 30.7 0 123
  52 25 0 131
  46 27.8 1 123
  47 28.6 1 136
  63 31.6 1 163
  35 32 1 133
  42 27.7 0 135
  37 29.9 0 123
  57 25.5 0 135
  64 29.8 0 127
  45 28 1 134
  58 29.7 1 141
;
run;

/* Fit, 95% CI for coefficients, VIF, and influence statistics */
ods graphics on;
proc reg data=bp plots(only label)=(diagnostics cooksd);
  model sbp = age bmi male / clb vif;
  output out=diag r=resid student=stdres cookd=cooks;
run;
quit;

/* Residual normality (Shapiro-Wilk) */
proc univariate data=diag normal;
  var resid;
run;

/* Any large Cook's D? (rule of thumb: > 4/n = 0.08) */
proc print data=diag;
  where cooks > 0.08;
  var sbp age bmi male resid stdres cooks;
run;

/* Heteroscedasticity test: add SPEC to the MODEL statement of PROC REG */
proc reg data=bp;
  model sbp = age bmi male / spec;
run;
quit;
```

**What the code does**

- `proc reg` fits OLS. `clb` requests 95% confidence limits for the coefficients and `vif` the variance inflation factors.
- `plots(only label)=(diagnostics cooksd)` produces the residual, Q–Q, leverage, and Cook's distance plots (requires ODS Graphics; `label` names the observation numbers).
- `output out=diag r= student= cookd=` saves residuals, studentized residuals, and Cook's D for further checks.
- `proc univariate ... normal` gives Shapiro–Wilk for the residuals.
- `proc print` lists observations with Cook's D above 4/n = 0.08.
- The `spec` option on a second `proc reg` runs White's test of constant variance (it differs from the Breusch–Pagan test in R and Python, so p-values will not match exactly).
- For a categorical predictor with more than two levels, use `proc glm` or `proc genmod` with a `class` statement. In PROC REG, `male` must already be a 0/1 numeric variable.

## Interpreting the output

- **Age: 0.74 mmHg per year (95% CI 0.56 to 0.93, p < 0.001).** Holding BMI and sex constant, each extra year of age is associated with a 0.74 mmHg higher SBP, or about 7.4 mmHg per decade.
- **BMI: 0.72 mmHg per kg/m² (95% CI 0.20 to 1.25, p = 0.007).** Each unit of BMI is associated with a 0.72 mmHg higher SBP, adjusted for age and sex.
- **Male: 4.2 mmHg (95% CI −0.2 to 8.5, p = 0.062).** Men averaged 4.2 mmHg higher SBP than women of the same age and BMI, but the data are also compatible with no difference or a modest difference in either direction. Absence of significance is not evidence of no effect.
- **Intercept = 75.0 mmHg.** The predicted SBP for a 0-year-old woman with BMI 0, which has no clinical meaning. Centre age and BMI (e.g., at 55 years and 27 kg/m²) if you want an interpretable intercept.
- **R² = 0.627, adjusted R² = 0.603; residual SE = 7.48 mmHg; F(3, 46) = 25.8, p < 0.001.** The model explains about 63% of the variation in SBP; typical prediction error is about ±7.5 mmHg (a 95% prediction interval is roughly ±15 mmHg).
- **VIF = 1.00, 1.03, 1.03.** No collinearity concern.
- **Residual diagnostics.** Breusch–Pagan p = 0.61 (no evidence of non-constant variance). The Shapiro–Wilk test rejects normality of residuals (p = 0.006); the reason is one outlier, patient 6 (a 57-year-old man with BMI 28.3 and SBP 116 mmHg, about 26 mmHg below the model's prediction of 142, with a studentized residual beyond ±3). Max Cook's D = 0.126, above the 4/n = 0.08 guide, so this patient is somewhat influential.
- **Sensitivity analysis without patient 6.** The age slope is unchanged in substance (95% CI 0.60 to 0.92), and R² rises to 0.71. The conclusions are robust, but the outlier should be investigated (data-entry error? true physiology?) and reported, not simply removed.

**How to report it:**
> In a multiple linear regression adjusted for BMI and sex (n = 50), each additional year of age was associated with a 0.74 mmHg higher systolic blood pressure (95% CI 0.56 to 0.93; p < 0.001). Each 1 kg/m² higher BMI was associated with 0.72 mmHg higher SBP (95% CI 0.20 to 1.25; p = 0.007). The model explained 63% of the variance (adjusted R² = 0.60). One outlier was identified; excluding it did not change the conclusions.

## Check your understanding

<details class="check"><summary>1. The Shapiro–Wilk test on the raw SBP values gives p = 0.30. A reviewer asks whether the outcome is normal. What should you check instead?</summary>

The **residuals** of the fitted model, with a Q–Q plot. Regression assumes normal errors, not a normal outcome. SBP can be skewed overall simply because age, BMI, and sex differ between patients.
</details>

<details class="check"><summary>2. The coefficient for male is 4.2 mmHg with p = 0.062 and a 95% CI of −0.2 to 8.5. A colleague says sex has no effect on SBP. Is that justified?</summary>

No. The CI includes differences of up to about 8.5 mmHg, which would matter clinically. p > 0.05 means the data cannot rule out zero, not that the true effect is zero. In 50 patients, the study is simply too small for a clear answer.
</details>

<details class="check"><summary>3. What would you expect to happen to the age coefficient's standard error if the same study enrolled only patients aged 50 to 60?</summary>

The SE would grow. A narrow age range means less variation in the predictor, so the slope is estimated less precisely. Wide predictor ranges give more signal for the same noise.
</details>

<details class="check"><summary>4. Using the model, how much higher is predicted SBP for a patient 10 years older, with the same BMI and sex?</summary>

About 10 × 0.742 = 7.4 mmHg (95% CI roughly 5.6 to 9.3 mmHg). Reporting per decade is often more meaningful clinically than per year.
</details>

## Common mistakes

- **Interpreting coefficients causally.** In observational data, adjusted associations can still be confounded; a regression coefficient is not a treatment effect unless the design supports it.
- **Overfitting.** Too many predictors for the sample size gives unstable coefficients and optimistic R². Specify the model beforehand; avoid stepwise selection.
- **Checking only the outcome's distribution.** Normality applies to the **residuals**, not to Y or to the predictors.
- **Dichotomising a continuous predictor** (e.g., BMI ≥ 30) to "simplify". This throws away information and power.
- **Extrapolating beyond the observed range** of the predictors (e.g., predicting for ages 20 or 90).
- **Ignoring influential points or deleting them silently.** Report a sensitivity analysis.
- **Reporting standardized coefficients only.** Clinical readers need the effect per meaningful unit (per 10 years, per 5 kg/m²).
- **Including a predictor that is a consequence of the exposure** (a mediator) when you want the total effect.
- **Treating repeated measures as independent.** Use a [mixed model](#/test/mixed-models).

## Related tests

- [Logistic regression](#/test/logistic-regression): binary outcome
- [Poisson & negative binomial regression](#/test/poisson-regression): count outcome
- [ANCOVA](#/test/ancova): treatment comparison adjusted for a continuous baseline covariate
- [Pearson correlation](#/test/pearson): unadjusted linear association
- [Linear mixed-effects models (MMRM)](#/test/mixed-models): repeated or clustered continuous outcomes
- [Cox proportional hazards regression](#/test/cox-ph): time-to-event outcome

## Key takeaways

- Each coefficient is the expected change in Y per one-unit change in X, **holding the other predictors constant**.
- t = estimate / SE; always report the estimate with its 95% CI on the clinical scale.
- Check **residuals** (not the raw outcome) for normality, constant variance, and influence.
- R² is the share of variance explained; the residual SE gives typical prediction error in outcome units.
- Adjusted associations from observational data are not necessarily causal.

## References

- Harrell FE. *Regression Modeling Strategies.* 2nd ed. Springer; 2015.
- Kutner MH, Nachtsheim CJ, Neter J, Li W. *Applied Linear Statistical Models.* 5th ed. McGraw-Hill; 2005.
- ICH E9 Statistical Principles for Clinical Trials (1998), section 5.7 on covariate adjustment.
- R documentation: `?lm`, `?confint`, `car::vif`; statsmodels: `statsmodels.formula.api.ols`; SAS/STAT User's Guide: The REG Procedure
