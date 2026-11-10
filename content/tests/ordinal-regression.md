# Ordinal logistic (proportional odds) regression

> Models an ordered outcome with several categories, such as the 90-day modified Rankin Scale, using a **single common odds ratio** that describes the shift of the whole distribution, and adjusts for covariates such as age.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Ordered outcome with 3+ categories; one summary effect with covariate adjustment | Ordinal (ordered categories) | Proportional odds; independent observations; linear in logit | 90-day mRS in a stroke trial, adjusted for age |

## Learning objectives

- Recognize ordered outcomes where a "shift" analysis beats dichotomising
- Explain the **common odds ratio** and the proportional odds assumption
- Fit the model in R and Python and read the OR, CI, and thresholds
- Check the sign convention and the proportional odds assumption before reporting

## What it is

Ordinal logistic regression, also called the **proportional odds model** or **cumulative logit model**, is used when the outcome has more than two categories that have a natural order but unequal spacing. Instead of collapsing the scale to one cut-point, the model looks at *every* cut-point at once. For an outcome Y with categories 0, 1, …, 6 it models the cumulative probabilities P(Y ≤ j) for j = 0 … 5 and assumes that the covariate effect on the log-odds is the **same** at every cut-point:

logit P(Y ≤ j) = ζⱼ − (β₁·treatment + β₂·age)

- **Intercepts (thresholds)** ζⱼ: one per cut-point, increasing in j. They describe the control-group distribution at the reference covariate values.
- **β:** the log of the **common odds ratio**. With the sign convention above, **OR = exp(β) > 1 means higher odds of being in a higher (worse) mRS category**; OR < 1 means a favourable shift toward better outcomes.
- **Test statistic:** Wald z = β̂ / SE(β̂) for each coefficient (a likelihood-ratio test is also available).
- **Effect size to report:** the common OR with 95% CI. In stroke trials this is often called the **shift analysis** (an "ordinal shift" of the mRS distribution).

**Sign conventions across software (read this before comparing outputs).**

| Software | Model it fits | Treatment coefficient here | Common OR reported on this page |
|---|---|---|---|
| R `MASS::polr` | logit P(Y ≤ j) = ζⱼ − η | −0.851 | exp(−0.851) = 0.43 |
| Python `OrderedModel` | logit P(Y ≤ j) = ζⱼ − η | −0.851 | 0.43 |
| SAS `PROC LOGISTIC` (default) | logit P(Y ≤ j) = αⱼ **+** η | +0.851 | 2.34 (this is the OR for a *lower* score, i.e. 1/0.43) |
| SAS with `DESCENDING` | models P(Y ≥ j) with +η | −0.851 | 0.43 |

All three languages give the same model. Only the labelling of the sign differs. This page always reports the OR for **a higher (worse) mRS category**, so use `DESCENDING` in SAS to get the same sign as R and Python.

## The intuition

Imagine the mRS scale as a staircase from 0 (no symptoms) to 6 (death). Dichotomising asks one question at one step: "is the patient on the good side of step 2?" The proportional odds model asks the same question at **every** step (0 vs 1–6, 0–1 vs 2–6, and so on) and assumes the drug has the *same odds ratio* at each one. The common OR is one number that summarizes how far the whole distribution shifts.

Signal vs noise: the signal is the estimated shift on the log-odds scale (β = −0.851). The noise is its standard error (0.365). The Wald z is their ratio: −0.851 / 0.365 = −2.33.

An everyday analogy: a drug that moves every patient one notch down the staircase helps everyone, including those who go from 5 to 4, which a "good outcome" cut at 2 would ignore entirely. The shift analysis uses that information, which is why it typically has more power.

Here are the observed counts. Notice how the drug arm piles up at mRS 0 (19 vs 6 patients), while the two arms look similar further up the scale.

| mRS | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|
| Control (n = 50) | 6 | 8 | 9 | 11 | 7 | 8 | 1 |
| Drug (n = 50) | 19 | 7 | 7 | 6 | 6 | 3 | 2 |
| Cumulative P(mRS ≤ j), control | 0.12 | 0.28 | 0.46 | 0.68 | 0.82 | 0.98 | 1.00 |
| Cumulative P(mRS ≤ j), drug | 0.38 | 0.52 | 0.66 | 0.78 | 0.90 | 0.96 | 1.00 |

The drug's cumulative probability is higher at nearly every cut-point: more patients at or below each level. That consistent gap is the "shift" the common OR describes.

Under H₀ (OR = 1), the Wald z for treatment follows a standard normal curve. Ours is −2.33, in the lower tail.

<div class="viz" data-viz="z" data-stat="-2.33"></div>

## When to use it

Use ordinal logistic regression when **all** of these are true:

- The outcome is **ordinal** with three or more ordered categories (mRS 0–6, pain none/mild/moderate/severe, NYHA class, Likert items)
- You want a **single summary effect** for a treatment or exposure, with adjustment for covariates
- The effect is **roughly constant across cut-points** (proportional odds)
- Observations are independent (one row per patient)

**Clinical scenario.** A randomised trial of 100 patients with acute ischaemic stroke compares a new neuroprotective drug (n = 50) with standard care (n = 50). The primary outcome is the **modified Rankin Scale at 90 days**: 0 = no symptoms, 1–5 = increasing disability, 6 = dead. *Does the drug shift the whole mRS distribution toward better outcomes, after adjusting for age?* The analysis uses mRS as a seven-level ordinal outcome (a "shift analysis"), with treatment and age as covariates.

**Dichotomising versus the shift analysis.** The traditional primary outcome is "good outcome", mRS 0–2 versus 3–6. Dichotomising throws away information: a patient moving from mRS 5 to 3 is "no change", but a patient moving from 3 to 2 is a "success". The proportional odds model uses every category, which typically gives **more statistical power** (or the same power with a smaller sample) and answers a broader question. The price is the proportional odds assumption and a less immediate interpretation, so many trials report both.

**Don't use it when:**

- The outcome is **binary** → use [logistic regression](#/test/logistic-regression)
- The outcome is **continuous** → use [linear regression](#/test/linear-regression)
- The outcome is a **count** → use [Poisson regression](#/test/poisson-regression)
- You only compare two or more groups with no covariates and are happy with a rank test → the [Mann–Whitney test](#/test/mann-whitney) or [Kruskal–Wallis test](#/test/kruskal-wallis) is enough
- Categories have **no natural order** (type of stroke, blood group) → use multinomial logistic regression
- The proportional odds assumption clearly fails for the main effect → use a **partial proportional odds** model, or report the cut-point-specific ORs

## Assumptions

1. **Proportional odds (parallel lines).** The OR for each covariate is the same at every cut-point of the outcome. Check with the **score test** (SAS prints this automatically as the "Score Test for the Proportional Odds Assumption"), with `ordinal::nominal_test()` in R (a likelihood-ratio test per covariate), with a Brant test, or informally by fitting separate binary logistic models for P(Y > j) and comparing the ORs. A significant test means the effect varies across cut-points. These tests are over-sensitive in large samples and under-powered in small ones, so look at the size of the differences too.
2. **Independent observations.** One row per patient. For repeated measures or clusters use [GEE](#/test/gee) or a [mixed model](#/test/mixed-models) with an ordinal link.
3. **Ordered categories with a sensible scale.** The coding must follow the clinical order. mRS 6 (death) is treated as just "worse than 5", which is itself a modelling choice.
4. **Linearity in the logit for continuous covariates.** Here, the log-odds change linearly with age. Check with splines or categories.
5. **No perfect separation and adequate counts per category.** Very sparse categories (for example only one patient with mRS 6 in an arm) make the thresholds unstable. Consider merging adjacent categories if the scale is sparse.

## Hypotheses

For the treatment coefficient β (log common odds ratio):

- **H₀:** β = 0, that is OR = 1 (the distribution of mRS is the same in both arms, adjusting for age)
- **H₁:** β ≠ 0, that is OR ≠ 1 (two-sided)

For the proportional odds check the null hypothesis is that the **effect is the same at every cut-point**, so here a *large* p-value is the reassuring one.

## Example data

The data are generated once in R (`set.seed(2026)`, a logistic latent-variable model) and typed in compactly in all three languages: 100 patients, treatment = 0 (control, first 50) or 1 (drug, last 50), age in years (range 45–90), and mRS at 90 days. Age is rescaled to **age10 = (age − 68) / 10**, so the age coefficient is per 10 years.

| mRS | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|
| Control (n = 50) | 6 | 8 | 9 | 11 | 7 | 8 | 1 |
| Drug (n = 50) | 19 | 7 | 7 | 6 | 6 | 3 | 2 |

mRS 0–2 (the usual "good outcome" cut) occurs in 23/50 control and 33/50 drug patients.

## How to run it

### R

```r
library(MASS)
trt <- rep(c(0, 1), each = 50)                 # 0 = control, 1 = drug
age <- c(73,58,69,67,62,45,61,59,69,64,64,61,66,66,45,80,74,70,61,74,65,67,55,81,68,85,84,69,74,84,63,69,66,71,70,78,73,60,73,61,58,75,57,71,60,81,74,64,75,72,57,66,60,72,62,66,57,59,69,59,67,70,66,65,52,66,52,78,59,69,73,68,80,60,72,59,87,74,65,75,58,82,58,54,68,76,68,61,79,90,50,81,64,77,52,66,59,64,58,71)
mrs <- c(3,2,3,1,4,0,3,1,4,4,5,3,1,4,0,5,3,2,3,3,5,2,0,6,1,5,5,3,2,5,2,4,4,5,1,5,1,2,2,3,1,3,0,0,4,0,2,1,3,2,1,6,0,0,2,4,4,0,1,2,0,0,4,0,1,0,2,5,1,4,3,2,5,1,3,0,3,0,0,0,2,5,0,1,6,3,2,0,0,3,4,4,0,0,0,2,0,3,1,0)
d <- data.frame(trt, age10 = (age - 68) / 10,
                mrs = factor(mrs, levels = 0:6, ordered = TRUE))

# 1. Proportional odds model; Hess = TRUE is needed for SEs and summary()
fit <- polr(mrs ~ trt + age10, data = d, Hess = TRUE, method = "logistic")
summary(fit)

# 2. Common odds ratios: profile-likelihood CI (default of confint) and Wald CI
exp(cbind(OR = coef(fit), confint(fit)))
exp(cbind(OR = coef(fit), confint.default(fit)))

# 3. Proportional odds check (likelihood-ratio test per covariate)
library(ordinal)
nominal_test(clm(mrs ~ trt + age10, data = d))

# 4. For comparison: the dichotomised outcome mRS >= 3 (poor outcome)
poor <- as.numeric(as.character(d$mrs)) >= 3
g <- glm(poor ~ trt + age10, data = d, family = binomial)
exp(cbind(OR = coef(g), confint.default(g)))
```

**What the code does**

- `mrs` must be an **ordered factor** with the levels in clinical order. If it is a plain number or an unordered factor, `polr` errors or fits the wrong thing.
- `polr(..., Hess = TRUE)` fits the proportional odds model. Without `Hess = TRUE` the Hessian is not stored and `summary()` cannot report standard errors.
- **Sign convention:** `polr` models logit P(Y ≤ j) = ζⱼ − η. A **negative** coefficient therefore means a shift toward **lower** (better) mRS, and exp(coef) is the OR for being in a **higher** category.
- `summary(fit)` reports no p-values. Use z = Value / Std. Error (`t value` in the output, but it is a z statistic).
- `confint(fit)` uses **profile likelihood** (it prints "Waiting for profiling to be done..."), which is more accurate than the Wald intervals from `confint.default(fit)`. Both are shown here.
- `nominal_test()` from the `ordinal` package re-fits the model letting each covariate have cut-point-specific effects and gives a likelihood-ratio test for each. If a covariate fails, the partial proportional odds model is `clm(mrs ~ trt, nominal = ~ age10, data = d)`.
- The dichotomised logistic model uses the usual "mRS ≥ 3" cut. Compare its CI with the common OR above.

**Output**

```text
Coefficients:
        Value Std. Error t value
trt   -0.8512     0.3651  -2.331
age10  0.8089     0.2134   3.791

Intercepts:
        Value   Std. Error t value
0|1    -1.7872  0.3211    -5.5667
1|2    -0.9893  0.2837    -3.4873
2|3    -0.2109  0.2659    -0.7931
3|4     0.6459  0.2773     2.3295
4|5     1.5652  0.3303     4.7382
5|6     3.3656  0.6112     5.5063

Residual Deviance: 346.3353 
AIC: 362.3353 

Profile-likelihood CI:
             OR     2.5 %   97.5 %
trt   0.4269003 0.2068704 0.868395
age10 2.2453747 1.4849466 3.441517

Wald CI:
             OR     2.5 %    97.5 %
trt   0.4269003 0.2087015 0.8732276
age10 2.2453747 1.4779951 3.4111800

Tests of nominal effects
       Df  logLik    AIC     LRT Pr(>Chi)   
<none>    -173.17 362.34                    
trt     5 -170.69 367.39  4.9484 0.422210   
age10   5 -164.19 354.38 17.9543 0.003004 **

Dichotomised (mRS >= 3):
                   OR     2.5 %   97.5 %
trt         0.4477591 0.1866995 1.073855
age10       2.6896798 1.5500276 4.667257
```

### Python

```python
import numpy as np, pandas as pd
from statsmodels.miscmodels.ordinal_model import OrderedModel

trt = np.repeat([0, 1], 50)                    # 0 = control, 1 = drug
age = np.array([73,58,69,67,62,45,61,59,69,64,64,61,66,66,45,80,74,70,61,74,65,67,55,81,68,85,84,69,74,84,63,69,66,71,70,78,73,60,73,61,58,75,57,71,60,81,74,64,75,72,57,66,60,72,62,66,57,59,69,59,67,70,66,65,52,66,52,78,59,69,73,68,80,60,72,59,87,74,65,75,58,82,58,54,68,76,68,61,79,90,50,81,64,77,52,66,59,64,58,71])
mrs = np.array([3,2,3,1,4,0,3,1,4,4,5,3,1,4,0,5,3,2,3,3,5,2,0,6,1,5,5,3,2,5,2,4,4,5,1,5,1,2,2,3,1,3,0,0,4,0,2,1,3,2,1,6,0,0,2,4,4,0,1,2,0,0,4,0,1,0,2,5,1,4,3,2,5,1,3,0,3,0,0,0,2,5,0,1,6,3,2,0,0,3,4,4,0,0,0,2,0,3,1,0])

X = pd.DataFrame({"trt": trt, "age10": (age - 68) / 10})
y = pd.Series(pd.Categorical(mrs, categories=range(7), ordered=True))

# Proportional odds (logit) model; no constant column in X
res = OrderedModel(y, X, distr="logit").fit(method="bfgs", disp=False, maxiter=500)
print(res.summary())

# Common odds ratios with Wald 95% CI (exponentiate coefficient and limits)
b = res.params[["trt", "age10"]]
ci = res.conf_int().loc[["trt", "age10"]]
print(pd.DataFrame({"OR": np.exp(b), "lo": np.exp(ci[0]), "hi": np.exp(ci[1])}).round(3))
```

**What the code does**

- `y` is an **ordered categorical**. Passing the integer codes 0–6 also works, but an ordered `Categorical` makes the category order explicit and keeps empty categories visible.
- `OrderedModel(..., distr="logit")` is the proportional odds model. Do **not** add a constant to `X`; the thresholds play the role of intercepts.
- **Sign convention:** the same as R's `polr`, logit P(Y ≤ j) = ζⱼ − η. The coefficients (−0.851 and 0.809) and the common ORs are identical to R. Negative means a shift to lower (better) mRS.
- **Threshold parametrisation differs.** The first threshold printed (`0/1`) is the first cut-point itself. The later thresholds are the **log of the increments** between successive cut-points, which keeps them ordered. For example, 1/2 = log(−0.9893 − (−1.7872)) = −0.226. To get R's ζⱼ, cumulate: ζ₂ = ζ₁ + exp(`1/2`).
- Only Wald intervals are available. `maxiter` and `method="bfgs"` are set because the default optimiser sometimes stops early.
- statsmodels has no built-in proportional odds test. Fit separate binary logits for P(Y > j) and compare the ORs, or use R or SAS for the formal test.

**Output**

```text
Log-Likelihood: -173.17   AIC: 362.3
                 coef    std err          z      P>|z|      [0.025      0.975]
------------------------------------------------------------------------------
trt           -0.8512      0.365     -2.331      0.020      -1.567      -0.136
age10          0.8089      0.213      3.791      0.000       0.391       1.227
0/1           -1.7873      0.321     -5.567      0.000      -2.417      -1.158
1/2           -0.2258      0.242     -0.933      0.351      -0.700       0.249
2/3           -0.2505      0.234     -1.071      0.284      -0.709       0.208
3/4           -0.1545      0.226     -0.684      0.494      -0.597       0.288
4/5           -0.0841      0.263     -0.320      0.749      -0.599       0.431
5/6            0.5880      0.306      1.920      0.055      -0.012       1.188

          OR     lo     hi
trt    0.427  0.209  0.873
age10  2.245  1.478  3.411
```

### SAS

```sas
/* trt = 0 control, 1 drug; age in years; mrs = 90-day modified Rankin Scale */
data stroke;
  input trt age mrs @@;
  age10 = (age - 68) / 10;
  datalines;
  0 73 3  0 58 2  0 69 3  0 67 1  0 62 4  0 45 0  0 61 3  0 59 1  0 69 4  0 64 4
  0 64 5  0 61 3  0 66 1  0 66 4  0 45 0  0 80 5  0 74 3  0 70 2  0 61 3  0 74 3
  0 65 5  0 67 2  0 55 0  0 81 6  0 68 1  0 85 5  0 84 5  0 69 3  0 74 2  0 84 5
  0 63 2  0 69 4  0 66 4  0 71 5  0 70 1  0 78 5  0 73 1  0 60 2  0 73 2  0 61 3
  0 58 1  0 75 3  0 57 0  0 71 0  0 60 4  0 81 0  0 74 2  0 64 1  0 75 3  0 72 2
  1 57 1  1 66 6  1 60 0  1 72 0  1 62 2  1 66 4  1 57 4  1 59 0  1 69 1  1 59 2
  1 67 0  1 70 0  1 66 4  1 65 0  1 52 1  1 66 0  1 52 2  1 78 5  1 59 1  1 69 4
  1 73 3  1 68 2  1 80 5  1 60 1  1 72 3  1 59 0  1 87 3  1 74 0  1 65 0  1 75 0
  1 58 2  1 82 5  1 58 0  1 54 1  1 68 6  1 76 3  1 68 2  1 61 0  1 79 0  1 90 3
  1 50 4  1 81 4  1 64 0  1 77 0  1 52 0  1 66 2  1 59 0  1 64 3  1 58 1  1 71 0
;
run;

/* Proportional odds model. DESCENDING makes the sign match R and Python:
   the OR is then for a HIGHER (worse) mRS category. */
proc logistic data=stroke descending;
  model mrs = trt age10 / link=clogit clodds=pl;   /* clogit = cumulative logit (default) */
run;

/* Dichotomised outcome for comparison: poor outcome mRS >= 3 */
data stroke2;
  set stroke;
  poor = (mrs >= 3);
run;
proc logistic data=stroke2 descending;
  model poor = trt age10;
run;
```

**What the code does**

- `input trt age mrs @@;` reads several patients per line (the double trailing `@` holds the line), then `age10` is computed in the DATA step.
- With an ordinal outcome (more than two levels), `PROC LOGISTIC` fits the **cumulative logit** model by default, so `link=clogit` is optional and written here only for clarity. `clodds=pl` adds profile-likelihood CIs for the odds ratios (`clodds=wald` for Wald).
- **Sign convention:** by default SAS orders the response ascending and models logit P(Y ≤ j) = αⱼ **+** β′x, so the printed β for treatment would be +0.851 (OR 2.34 for a *lower* score). **`DESCENDING`** reverses the ordering, so β = −0.851 and the OR is 0.43 for a higher score, matching R and Python. The intercepts also change sign and meaning with `DESCENDING`, so do not compare them directly with R.
- SAS automatically prints **"Score Test for the Proportional Odds Assumption"** (a single global chi-square with df = (number of cut-points − 1) × number of covariates, here 10). It is a test of both covariates jointly. A small p-value means at least one covariate violates the assumption. SAS prints no per-covariate test, so use the R `nominal_test` or fit separate `proc logistic` models for P(Y > j) for diagnosis. The partial proportional odds model is available through the `UNEQUALSLOPES` option on the MODEL statement (for example `unequalslopes=age10`).
- The "Odds Ratio Estimates" table gives the common OR for `trt` and `age10`.

## Interpreting the output

- **Treatment coefficient β = −0.851, common OR = 0.43 (95% CI 0.21 to 0.87 profile, 0.21 to 0.87 Wald).** The odds of being in a *higher* (worse) mRS category rather than a lower one are about 57% lower with the drug than with control, at the same age. Equivalently, the odds of a *better* outcome at any cut-point are about 2.3 times higher (1/0.43). The CI excludes 1, in agreement with z = −2.33 and p = 0.020.
- **Age (per 10 years): OR = 2.25 (95% CI 1.48 to 3.44).** Each additional 10 years of age multiplies the odds of a worse mRS category by about 2.25, with treatment held fixed.
- **Thresholds ζⱼ (−1.79 to 3.37).** These are the log-odds of mRS ≤ j for a control patient aged 68. They are rarely interpreted directly, but they must increase with j. For example the control-group probability of mRS ≤ 2 at age 68 is 1/(1 + exp(0.211)) = 0.45 (the observed control proportion is 23/50 = 0.46).
- **R and Python agree** on coefficients, SEs and log-likelihood (−173.17). Only the threshold parametrisation differs (see above).
- **Proportional odds check.** For treatment the LR test is non-significant (p = 0.42), so a single common treatment OR is reasonable. For age it is significant (LRT 17.95 on 5 df, p = 0.003), so the age effect is not constant over the scale. Treatment is the quantity of interest, and the partial proportional odds model that frees only age gives a treatment OR of 0.44 (95% CI 0.21 to 0.89), almost unchanged. In a real trial, this check would be pre-specified and the sensitivity analysis reported.
- **Separate binary logits** for P(mRS > j) give treatment ORs of 0.37, 0.45, 0.63, 0.48 and 2.17 for j = 1 to 5. They are all below 1 except the last, which rests on only 3 patients at mRS 6. This pattern is consistent with a common OR of about 0.4 to 0.5.
- **Dichotomised analysis (mRS ≥ 3):** OR 0.45 (95% CI 0.19 to 1.07), which is no longer significant. The point estimate is similar, but the CI is about 20% wider on the log scale (log-width 1.75 versus 1.43). This is the typical gain in precision from the shift analysis.

**How to report it:**
> Treatment shifted the 90-day modified Rankin Scale distribution toward better outcomes (common odds ratio for a higher mRS category 0.43, 95% CI 0.21 to 0.87, adjusted for age; proportional odds model). The proportional odds assumption was acceptable for treatment (likelihood-ratio p = 0.42).

## Check your understanding

<details class="check"><summary>1. In R, polr gives a treatment coefficient of −0.851. Does the drug help or harm?</summary>

It helps. In `polr`, a **negative** coefficient means a shift toward lower (better) mRS. exp(−0.851) = 0.43 is the odds ratio for being in a *higher* (worse) category. Always state which direction your OR refers to.
</details>

<details class="check"><summary>2. A trial reports an ordinal shift OR = 0.43 and also a dichotomised mRS ≥ 3 analysis with OR = 0.45 (CI 0.19 to 1.07). Why might the first be significant and the second not?</summary>

The point estimates are similar, but the shift analysis uses all seven categories, so its estimate is more precise (narrower CI). Dichotomising discards how far patients move within each side of the cut, which costs power.
</details>

<details class="check"><summary>3. The proportional odds test for age is significant (p = 0.003), but for treatment p = 0.42. Is the treatment OR invalid?</summary>

No. The assumption is checked per covariate. The treatment effect looks constant across cut-points, so a single OR is reasonable. A partial proportional odds model freeing only age gave almost the same treatment OR (0.44). These tests are also over-sensitive in large samples and under-powered in small ones, so look at the size of the differences too.
</details>

<details class="check"><summary>4. What if the drug helped greatly at low mRS cut-points but not at all at high ones (a non-constant effect)? What goes wrong with a single common OR?</summary>

The common OR would average over very different effects and describe none of them well. Check the proportional odds assumption, and if it fails, use a partial proportional odds model or report cut-point-specific ORs.
</details>

## Common mistakes

- **Misreading the sign.** `polr` and statsmodels use ζⱼ − η, SAS defaults to αⱼ + η. Check which direction your OR refers to and say so in the manuscript ("OR for a worse outcome" or "for a better outcome").
- **Forgetting `Hess = TRUE`** in `polr`, or not ordering the factor levels correctly.
- **Interpreting the common OR as the OR for each cut-point** without checking proportional odds. If the assumption fails, report a partial proportional odds model or the cut-point-specific ORs.
- **Relying on a single global test.** The SAS score test is global; a failure may be driven by one covariate (here age) and not the treatment.
- **Treating mRS as continuous** (linear regression or t-test on the 0–6 score) when the steps are not equally spaced.
- **Adding a constant** to `X` in statsmodels, or comparing the intercepts between software without noting the parametrisation.
- **Merging sparse categories after seeing the data** without saying so. Pre-specify the collapsing rule.
- **Reporting only the p-value.** Always give the common OR with its 95% CI and the direction.

## Related tests

- [Logistic regression](#/test/logistic-regression): the binary special case (two categories)
- [Mann–Whitney test](#/test/mann-whitney): rank-based comparison of two groups without covariates
- [Kruskal–Wallis test](#/test/kruskal-wallis): rank-based comparison of three or more groups
- [Poisson regression](#/test/poisson-regression): count outcomes
- [GEE](#/test/gee): marginal ordinal models for repeated or clustered data

## Key takeaways

- Ordinal logistic regression fits one **common OR** across all cut-points of an ordered outcome.
- OR sign depends on software: R and Python give the OR for a *higher* category; SAS default gives it for a *lower* one.
- Shift analysis uses the whole scale, so it usually has more power than dichotomising.
- Check proportional odds per covariate, and look at effect sizes, not only p-values.
- If it fails, use a partial proportional odds model or report cut-point-specific ORs.

## References

- Agresti A. *Categorical Data Analysis.* 3rd ed. Wiley; 2013. Ch. 8 (cumulative logit models).
- Saver JL. Novel end point analytic techniques and interpreting shifts across the entire range of outcome scales in acute stroke trials. *Stroke.* 2007;38:3055–62.
- R documentation: `MASS::polr`, `ordinal::clm`, `ordinal::nominal_test`; statsmodels: `OrderedModel`
- SAS/STAT User's Guide: The LOGISTIC Procedure (cumulative logit model, score test for the proportional odds assumption, `UNEQUALSLOPES`)
