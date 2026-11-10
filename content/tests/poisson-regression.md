# Poisson & negative binomial regression

> Models event **counts** (exacerbations, seizures, infections) as a function of predictors and reports **rate ratios**, using person-time as an offset; the negative binomial version handles the extra variability (overdispersion) that real clinical counts nearly always show.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Model event counts and compare event rates between groups, with unequal follow-up | Non-negative integer counts plus person-time | Independent patients; log-linear mean; offset for follow-up; check overdispersion | COPD exacerbations on new drug vs control |

## Learning objectives

- Explain why counts with unequal follow-up need a rate model with a person-time offset
- Interpret exp(β) as a rate ratio, with its 95% CI
- Detect overdispersion and explain why ignoring it produces over-confident results
- Choose between Poisson, negative binomial, quasi-Poisson, and robust standard errors

## What it is

Poisson regression is a generalised linear model for counts with a log link: log(μ) = log(t) + β₀ + β₁X₁ + …, where μ is the expected number of events and t the person-time at risk (follow-up). The term log(t) is the **offset**, a predictor whose coefficient is fixed at 1, so the model is really for the event *rate* μ / t. Exponentiating a coefficient gives a **rate ratio (RR) = exp(β)**: the multiplicative change in the event rate per unit of the predictor, or for treated versus control patients.

The Poisson distribution assumes the variance equals the mean. Clinical counts usually vary more than that (some patients are simply more exacerbation-prone), which is **overdispersion**. It leaves the rate ratios roughly right but makes the standard errors too small, so CIs are too narrow and p-values too optimistic. The **negative binomial (NB2)** model adds a dispersion parameter: Var(Y) = μ + αμ² (α = 1/θ), and corrects the standard errors.

- **Test statistic:** Wald z = β̂ / SE(β̂); likelihood-ratio χ² for model comparison
- **Effect size:** **rate ratio exp(β) with 95% CI** (primary); also report the observed event rate per person-year in each group. For trials, pre-specify the count model, the offset and how overdispersion will be handled (ICH E9).
- **Overdispersion check:** Pearson χ² / residual df; values well above 1 indicate overdispersion (about 1 is expected under the Poisson model)

## The intuition

Counting exacerbations is like counting buses at a stop. If you watch for two hours you expect more buses than if you watch for 30 minutes, so what matters is buses **per hour**, not the raw count. That is the role of the offset: each patient's count is compared with their own follow-up time, and the model estimates a *rate*. The rate ratio then says how many times higher or lower one group's rate is.

- **Signal:** the log rate ratio, here log(0.48) ≈ −0.74.
- **Noise:** how much counts vary beyond the average rate. A pure Poisson process assumes the variance equals the mean. Real patients differ in how exacerbation-prone they are, so counts bunch into "none" and "lots", and the variance exceeds the mean. That extra noise is **overdispersion**. If the model ignores it, the noise is underestimated and the signal looks more convincing than it is.

Look at each patient's own crude rate (exacerbations per person-year) below. Most dots sit low, but a few patients in both arms reach 3 or more, and a number sit at exactly zero, so the spread is wide. The drug arm has more zeros and fewer high rates.

<div class="viz" data-viz="dots" data-unit="exacerbations per person-year" data-groups='{"Control":[2.0,3.3,5.0,0.6,2.2,1.2,0.0,0.8,1.4,0.5,0.8,0.0,0.5,0.9,0.7,1.6,1.7,0.0,3.6,1.6,0.9,3.7,1.1,0.0,0.0],"Drug":[0.0,1.8,3.0,0.0,0.9,0.0,0.7,0.5,0.0,0.0,0.7,0.0,1.0,0.0,0.7,3.2,0.0,0.0,1.2,0.5,0.0,0.0,0.0,0.0,3.0]}' data-caption="Per-patient exacerbation rates (rounded to 0.1). Pooled arm rates (total events / total person-years) are 1.27 and 0.62; the means of the individual rates shown are 1.36 and 0.69. Patients vary a lot."></div>

Overdispersion is visible in the raw counts too: in the control arm the mean count is 1.68 but the variance is 2.56 (ratio 1.5), which would be 1.0 under a pure Poisson model. This matches the Pearson χ²/df of 1.59.

For the negative binomial model, the drug effect has log(RR) = −0.74 and a standard error of about 0.32 (back-calculated from the CI in the output), so z ≈ −2.31. Under H₀ (RR = 1), z follows the standard normal curve below. The observed z is in the tail, but it is not as extreme as a naive Poisson fit would suggest, which gives z ≈ −2.65.

<div class="viz" data-viz="z" data-stat="-2.31"></div>

What to notice: allowing for overdispersion widens the CI (0.26 to 0.89 versus 0.28 to 0.83) and moves z toward the middle, without changing the estimated rate ratio. That is the cost of honest standard errors.

## When to use it

Use Poisson or negative binomial regression when:

- The outcome is a **non-negative integer count** (number of COPD exacerbations, hospitalisations, relapses, falls)
- Patients are followed for **different lengths of time**, so the count must be related to person-time (offset)
- Events are reasonably independent within a patient (no heavy clustering in time)
- You want **rate ratios** adjusted for covariates

**Clinical scenario.** A randomised trial allocates 50 patients with COPD to a new inhaled drug (n = 25) or control (n = 25). The number of exacerbations is recorded over follow-up of 0.5 to 2 years (person-years), and patients are classified as having severe baseline disease or not. *By how much does the drug reduce the exacerbation rate?*

**Don't use it when:**

- The outcome is **binary** (any exacerbation yes/no, with fixed follow-up) → use [logistic regression](#/test/logistic-regression)
- The outcome is the **time to the first event** → use [Cox regression](#/test/cox-ph) or the [log-rank test](#/test/log-rank)
- The outcome is **continuous** → use [linear regression](#/test/linear-regression)
- The same patient is measured repeatedly (e.g., counts per month) → use [GEE](#/test/gee) or a [mixed model](#/test/mixed-models) with a count distribution
- There are **excess zeros** (many patients can never have the event) → consider zero-inflated or hurdle models

## Assumptions

1. **Independent observations.** One count per patient (or cluster-robust methods for repeated counts).
2. **Correct mean model.** The log of the expected rate is linear in the predictors; check with residual plots or by adding splines for continuous predictors.
3. **Offset specified.** The log of person-time enters with a fixed coefficient of 1. Omitting it models counts rather than rates, so unequal follow-up biases the results.
4. **Variance = mean (Poisson only).** This is the assumption most often violated, and ignoring it is the classic error. Check the Pearson χ² / df and compare Poisson with negative binomial by AIC or a likelihood-ratio test. If dispersion is above 1, use negative binomial, quasi-Poisson, or Poisson with robust (sandwich) standard errors. In a trial, decide this in the analysis plan, not after seeing which gives the smaller p-value. A dispersion check can itself be unreliable in small samples.
5. **Events occur at a constant rate over time** (within groups), and follow-up is non-informative. Dropout because of exacerbations violates this.
6. **No multicollinearity and no influential points.** Check VIFs and Cook's distance, and refit without flagged patients as a sensitivity analysis.

## Hypotheses

For each predictor j:

- **H₀:** βⱼ = 0, i.e. RR = 1 (the event rate does not differ, after adjustment)
- **H₁:** βⱼ ≠ 0, i.e. RR ≠ 1 (two-sided)

## Example data

Simulated data for 50 patients (25 per arm; 16 with severe disease). Columns: `treat` (1 = new drug), `severe` (1 = yes), `py` (person-years of follow-up), `count` (exacerbations). The first five rows:

| Patient | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| Treat | 0 | 0 | 0 | 0 | 0 |
| Severe | 0 | 0 | 1 | 0 | 0 |
| Person-years | 1.98 | 0.61 | 0.8 | 1.62 | 1.35 |
| Exacerbations | 4 | 2 | 4 | 1 | 3 |

Overall, the control arm had 42 exacerbations in 33.1 person-years (1.27 per person-year) and the drug arm 19 in 30.6 person-years (0.62 per person-year). All 50 rows are typed inline in each code block below.

## How to run it

### R

```r
# COPD exacerbations over follow-up (count), by treatment, disease severity, person-years
dat <- data.frame(
  treat = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1),
  severe = c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0,
      1, 0, 0, 1),
  py = c(1.98, 0.61, 0.8, 1.62, 1.35, 1.66, 1.11, 1.3, 1.41, 1.99, 1.27, 1.21,
      1.91, 1.16, 1.48, 0.63, 0.58, 0.53, 1.69, 1.22, 1.17, 1.35, 1.83, 1.96,
      1.26, 1.93, 1.66, 0.66, 1.72, 1.06, 0.69, 1.46, 1.88, 0.64, 1.56, 1.39,
      1.01, 1.03, 1.15, 1.42, 0.63, 0.83, 1.14, 0.81, 1.97, 0.77, 1.68, 0.83,
      1.36, 1.33),
  count = c(4, 2, 4, 1, 3, 2, 0, 1, 2, 1, 1, 0, 1, 1, 1, 1, 1, 0, 6, 2, 1, 5,
      2, 0, 0, 0, 3, 2, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 2, 0, 0, 1, 1, 0,
      0, 0, 0, 4)
)
dat$logpy <- log(dat$py)                     # offset: log of person-years

# 1. Poisson rate-ratio model with a log(person-time) offset
pois <- glm(count ~ treat + severe + offset(logpy), family = poisson, data = dat)
round(exp(cbind(RR = coef(pois), confint.default(pois))), 3)

# 2. Check for overdispersion: Pearson chi-square / residual df (should be ~1)
sum(residuals(pois, type = "pearson")^2) / df.residual(pois)

# 3. Negative binomial model (MASS::glm.nb estimates theta; variance = mu + mu^2/theta)
library(MASS)
nb <- glm.nb(count ~ treat + severe + offset(logpy), data = dat)
round(exp(cbind(RR = coef(nb), confint.default(nb))), 3)
c(theta = nb$theta, alpha = 1 / nb$theta)
AIC(pois, nb)

# 4. Multicollinearity and influence
car::vif(nb)
max(cooks.distance(nb))

# 5. Rates per person-year by arm (descriptive)
tapply(dat$count, dat$treat, sum) / tapply(dat$py, dat$treat, sum)
```

**What the code does**

- `offset(logpy)` in the formula supplies log(person-years) with a fixed coefficient of 1. Equivalent: `glm(..., offset = logpy)`. Forgetting it is a common error.
- `family = poisson` uses the log link by default. `exp(cbind(coef, confint.default))` turns coefficients and Wald CIs into rate ratios.
- The **dispersion statistic** is the sum of squared Pearson residuals divided by the residual df. A value near 1 is consistent with Poisson.
- `MASS::glm.nb` fits the negative binomial. It estimates θ (`nb$theta`); the variance is μ + μ²/θ, so the α reported by other software is 1/θ. Note that `glm.nb` treats θ as known when computing standard errors, so its CIs are slightly narrower than those from a joint-likelihood fit (e.g., Python or SAS).
- `AIC(pois, nb)` compares the models; a lower AIC is better (a difference of 2 or more is usually meaningful).
- `car::vif` and `cooks.distance` are the standard multicollinearity and influence checks. The last line gives the crude event rate per person-year in each arm.

**Output**

```text
               RR 2.5 % 97.5 %
(Intercept) 1.142 0.798  1.636
treat       0.481 0.280  0.827
severe      1.363 0.815  2.280
[1] 1.589229
               RR 2.5 % 97.5 %
(Intercept) 1.166 0.758  1.794
treat       0.478 0.256  0.894
severe      1.368 0.734  2.548
    theta     alpha 
3.3511426 0.2984057 
     df      AIC
pois  3 153.6450
nb    4 152.8259
   treat   severe 
1.005575 1.005575 
[1] 0.2039313
        0         1 
1.2696493 0.6207122
```

### Python

```python
import numpy as np, pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import variance_inflation_factor as vif

# COPD exacerbations over follow-up (count), by treatment, disease severity, person-years
dat = pd.DataFrame({
  'treat': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1],
  'severe': [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0,
      1, 0, 0, 1],
  'py': [1.98, 0.61, 0.8, 1.62, 1.35, 1.66, 1.11, 1.3, 1.41, 1.99, 1.27, 1.21,
      1.91, 1.16, 1.48, 0.63, 0.58, 0.53, 1.69, 1.22, 1.17, 1.35, 1.83, 1.96,
      1.26, 1.93, 1.66, 0.66, 1.72, 1.06, 0.69, 1.46, 1.88, 0.64, 1.56, 1.39,
      1.01, 1.03, 1.15, 1.42, 0.63, 0.83, 1.14, 0.81, 1.97, 0.77, 1.68, 0.83,
      1.36, 1.33],
  'count': [4, 2, 4, 1, 3, 2, 0, 1, 2, 1, 1, 0, 1, 1, 1, 1, 1, 0, 6, 2, 1, 5,
      2, 0, 0, 0, 3, 2, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 2, 0, 0, 1, 1, 0,
      0, 0, 0, 4]
})

# 1. Poisson rate-ratio model; exposure= adds log(person-years) as the offset
pois = smf.glm("count ~ treat + severe", data=dat, family=sm.families.Poisson(),
               exposure=dat["py"]).fit()
ci = pois.conf_int()
print(pd.DataFrame({"RR": np.exp(pois.params), "2.5%": np.exp(ci[0]),
                    "97.5%": np.exp(ci[1])}).round(3))

# 2. Overdispersion: Pearson chi-square / residual df (should be ~1)
print("Pearson chi2/df = %.3f" % (pois.pearson_chi2 / pois.df_resid))

# 3. Negative binomial (NB2) fitted by maximum likelihood; alpha = 1/theta
nb = smf.negativebinomial("count ~ treat + severe", data=dat,
                          exposure=dat["py"]).fit(disp=False)
ci = nb.conf_int()
print(pd.DataFrame({"RR": np.exp(nb.params), "2.5%": np.exp(ci[0]),
                    "97.5%": np.exp(ci[1])}).round(3).iloc[:3])
alpha = nb.params["alpha"]
print("alpha = %.3f, theta = 1/alpha = %.3f" % (alpha, 1 / alpha))
print("AIC Poisson = %.2f, AIC NB = %.2f" % (pois.aic, nb.aic))

# 4. Multicollinearity (VIF on the predictors) and influence
X = pois.model.exog
print("VIF:", np.round([vif(X, i) for i in range(1, X.shape[1])], 2).tolist())

# Cook's distance for the NB model, computed by hand (same formula as R)
mu = nb.predict(dat, exposure=dat["py"]).to_numpy()
w = mu / (1 + alpha * mu)                       # NB working weights
Xw = X * np.sqrt(w)[:, None]
h = np.diag(Xw @ np.linalg.inv(Xw.T @ Xw) @ Xw.T)   # hat values
pr = (dat["count"] - mu) / np.sqrt(mu + alpha * mu**2)   # Pearson residuals
print("max Cook's D = %.3f" % (pr**2 * h / (3 * (1 - h)**2)).max())

# 5. Rates per person-year by arm (descriptive)
g = dat.groupby("treat")[["count", "py"]].sum()
print((g["count"] / g["py"]).round(3).to_dict())
```

**What the code does**

- `exposure=dat["py"]` tells statsmodels to add log(exposure) as the offset; no need to take the log yourself. (You can also pass `offset=np.log(py)`.) For `GLM` and `NegativeBinomial` alike, this is the equivalent of R's `offset()`.
- `sm.families.Poisson()` is the family; the default link is log. `pois.pearson_chi2 / pois.df_resid` is the dispersion statistic.
- `smf.negativebinomial` fits the NB2 model by maximum likelihood and estimates `alpha` jointly with the coefficients (θ = 1/α). It prints `alpha` as an extra parameter, which is why the OR table is trimmed with `.iloc[:3]`. Its CIs are slightly wider than R's because the uncertainty in α enters the covariance matrix.
- `sm.families.NegativeBinomial` (the GLM family) requires you to fix α yourself; that is not the same as estimating it.
- The influence calculation is done by hand because the built-in influence class does not reproduce R's weights for a GLM with an offset; it matches R's `cooks.distance` here.

**Output**

```text
              RR   2.5%  97.5%
Intercept  1.142  0.798  1.636
treat      0.481  0.280  0.827
severe     1.363  0.815  2.280
Pearson chi2/df = 1.589
              RR   2.5%  97.5%
Intercept  1.166  0.755  1.801
treat      0.478  0.254  0.899
severe     1.368  0.731  2.559
alpha = 0.298, theta = 1/alpha = 3.351
AIC Poisson = 153.64, AIC NB = 152.83
VIF: [1.01, 1.01]
max Cook's D = 0.204
{0: 1.27, 1: 0.621}
```

### SAS

```sas
/* COPD exacerbations: count, treat (1 = new drug), severe (1 = yes), person-years */
data exac;
  input treat severe py count;
  datalines;
  0 0 1.98 4
  0 0 0.61 2
  0 1 0.8 4
  0 0 1.62 1
  0 0 1.35 3
  0 0 1.66 2
  0 0 1.11 0
  0 0 1.3 1
  0 0 1.41 2
  0 0 1.99 1
  0 1 1.27 1
  0 0 1.21 0
  0 1 1.91 1
  0 0 1.16 1
  0 1 1.48 1
  0 0 0.63 1
  0 0 0.58 1
  0 0 0.53 0
  0 1 1.69 6
  0 0 1.22 2
  0 1 1.17 1
  0 0 1.35 5
  0 1 1.83 2
  0 0 1.96 0
  0 0 1.26 0
  1 0 1.93 0
  1 0 1.66 3
  1 0 0.66 2
  1 0 1.72 0
  1 0 1.06 1
  1 0 0.69 0
  1 0 1.46 1
  1 1 1.88 1
  1 1 0.64 0
  1 0 1.56 0
  1 0 1.39 1
  1 1 1.01 0
  1 1 1.03 1
  1 0 1.15 0
  1 1 1.42 1
  1 0 0.63 2
  1 0 0.83 0
  1 1 1.14 0
  1 1 0.81 1
  1 0 1.97 1
  1 0 0.77 0
  1 1 1.68 0
  1 0 0.83 0
  1 0 1.36 0
  1 1 1.33 4
;
run;

data exac2;
  set exac;
  logpy = log(py);                        /* offset: log person-time */
run;

/* 1. Poisson rate-ratio model. ESTIMATE ... / EXP gives the rate ratio with CI.
      The "Criteria For Assessing Goodness Of Fit" table shows
      Pearson Chi-Square / DF: the overdispersion check. */
proc genmod data=exac2;
  model count = treat severe / dist=poisson link=log offset=logpy type3;
  estimate 'Rate ratio, new drug vs control' treat 1 / exp;
run;

/* 2. Negative binomial: DIST=NEGBIN also estimates the dispersion
      parameter (labelled "Dispersion"; it equals alpha = 1/theta). */
proc genmod data=exac2;
  model count = treat severe / dist=negbin link=log offset=logpy type3;
  estimate 'Rate ratio, new drug vs control' treat 1 / exp;
run;
```

**What the code does**

- `offset=logpy` is the SAS equivalent of R's `offset()`; the variable must be created beforehand as log(person-time).
- `dist=poisson link=log` fits Poisson regression. The goodness-of-fit table reports Pearson Chi-Square and its DF; divide to get the dispersion statistic (about 1.59 here, same as R and Python). Adding `scale=pearson` to the MODEL statement would fit a quasi-Poisson model that inflates the standard errors by √dispersion.
- `estimate '...' treat 1 / exp;` prints the treatment log-rate-ratio and, with `exp`, the rate ratio and its 95% CI.
- `type3` gives likelihood-ratio tests for each predictor.
- `dist=negbin` fits NB2 and reports the dispersion parameter (α = 1/θ) with its own CI.
- `proc genmod` has no VIF option; compute collinearity with `proc reg ... / vif` if needed. Case-level residuals and leverage are available with `/ obstats residuals` on the MODEL statement.

## Interpreting the output

- **Poisson, drug effect: RR = 0.48 (95% CI 0.28 to 0.83).** After adjusting for severity, patients on the drug had about half the exacerbation rate of controls (a 52% reduction).
- **Overdispersion: Pearson χ² / df = 1.59.** Counts vary more than a Poisson model allows (1.0 is expected), so the Poisson CI is probably too narrow. A modest value like this is typical; values of 2 or more would be more serious.
- **Negative binomial, drug effect: RR = 0.48 (95% CI 0.26 to 0.89 in R; 0.25 to 0.90 in Python).** The point estimate is essentially unchanged, but the interval is wider, which is the honest answer. The dispersion is α = 0.30 (θ = 3.35). The NB model fits slightly better (AIC 152.8 versus 153.6).
- **Severe disease: RR = 1.37 (95% CI 0.73 to 2.55 under NB).** The rate was about 37% higher in severe disease, but the interval is compatible with no difference.
- **Baseline rate (intercept): 1.17 exacerbations per person-year** for control patients with non-severe disease (NB model). This is the rate when all other predictors are 0.
- **VIF = 1.01, max Cook's D = 0.20.** No collinearity, and no observation is a concern.
- **Crude rates: 1.27 per person-year (control) versus 0.62 (drug).** The crude ratio is 0.49, close to the adjusted RR, as expected in a randomised trial.

**How to report it:**
> Over a median follow-up of 1.29 person-years per patient, the exacerbation rate was 0.62 per person-year on the new drug versus 1.27 on control. In a negative binomial model with log person-time offset adjusted for baseline severity, the rate ratio was 0.48 (95% CI 0.26 to 0.89), a 52% reduction. Overdispersion was present (Pearson χ²/df = 1.59), so negative binomial standard errors were used.

## Check your understanding

<details class="check"><summary>1. Patient A has 2 exacerbations in 0.5 years of follow-up and patient B has 3 in 2 years. Who has the higher rate?</summary>

Patient A: 2 / 0.5 = 4 per person-year, versus 3 / 2 = 1.5 for B. Comparing raw counts would point the wrong way, which is why the offset (log person-time) is essential.
</details>

<details class="check"><summary>2. A Poisson model gives RR = 0.48 with CI 0.28 to 0.83, and Pearson χ²/df = 1.59. A colleague says the Poisson model is fine because the RR is "so significant". What is the problem?</summary>

The dispersion is above 1, so Poisson standard errors are too small and the CI too narrow. The estimate of 0.48 is roughly right, but the confidence is overstated. A negative binomial model, or robust standard errors, gives a more honest CI (here 0.26 to 0.89).
</details>

<details class="check"><summary>3. What if you had forgotten the offset in this trial, and drug-arm patients happened to be followed for half as long as control patients?</summary>

You would be modelling counts, not rates. The drug arm would look like it had fewer events even if the true rate were identical, because shorter follow-up means fewer events. The rate ratio would be biased towards benefit.
</details>

<details class="check"><summary>4. The rate ratio is 0.48. Does that mean 48% of patients on the drug had an exacerbation?</summary>

No. It is a ratio of event *rates* per unit of person-time: the drug group had about 0.48 times the control rate of exacerbations, a 52% reduction. It says nothing about the proportion of patients with at least one event, which would be a risk or odds ratio.
</details>

## Common mistakes

- **Forgetting the offset** when follow-up differs between patients. The model then compares counts, not rates.
- **Ignoring overdispersion.** Poisson standard errors are too small when variance exceeds the mean; the result can be a false-positive finding. Always check, and report the model used.
- **Treating the rate ratio as a risk ratio or odds ratio.** An RR from a count model is a ratio of event *rates*, not of the proportion of patients with an event.
- **Using `scale=pearson` or `quasi-Poisson` and then reporting AIC or likelihood-based comparisons.** Quasi-likelihood models have no true likelihood.
- **Counting events per patient but modelling per visit (or vice versa)** without adjusting for clustering.
- **Ignoring excess zeros.** If many patients have no events and cannot have any, a zero-inflated model may fit better.
- **Taking the log of zero person-time.** Each patient needs positive follow-up.
- **Informative dropout.** Patients who leave after an event bias the rate estimate.

## Related tests

- [Logistic regression](#/test/logistic-regression): binary outcome
- [Linear regression](#/test/linear-regression): continuous outcome
- [Cox proportional hazards regression](#/test/cox-ph): time to first event
- [Generalized estimating equations (GEE)](#/test/gee): repeated or clustered counts
- [Log-rank test](#/test/log-rank): comparing time-to-event curves
- [Linear mixed-effects models (MMRM)](#/test/mixed-models): longitudinal continuous outcomes

## Key takeaways

- Poisson regression models counts with a log link. exp(β) is a **rate ratio**, and log person-time enters as an offset.
- Poisson assumes variance = mean. Real clinical counts are usually overdispersed, which makes Poisson CIs too narrow.
- Check Pearson χ²/df; if it is above 1, use negative binomial, quasi-Poisson, or robust standard errors.
- Report rates per person-year in each arm and the rate ratio with its 95% CI.
- Forgetting the offset, or ignoring overdispersion, are the two classic errors.

## References

- Cameron AC, Trivedi PK. *Regression Analysis of Count Data.* 2nd ed. Cambridge University Press; 2013.
- Keene ON, Jones MRK, Lane PW, Anderson J. Analysis of exacerbation rates in asthma and chronic obstructive pulmonary disease: example from the TRISTAN study. *Pharm Stat.* 2007;6:89–97.
- ICH E9 Statistical Principles for Clinical Trials (1998).
- R documentation: `?glm`, `MASS::glm.nb`; statsmodels: `GLM`, `NegativeBinomial`; SAS/STAT User's Guide: The GENMOD Procedure (DIST=NEGBIN, OFFSET=)
