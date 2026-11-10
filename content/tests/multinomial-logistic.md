# Multinomial logistic regression

> Models an **unordered categorical outcome with three or more levels** (home, rehab, nursing facility) as a function of predictors, reporting each predictor's effect as a **relative risk ratio** against a chosen reference category.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Outcome has 3+ unordered categories and you want adjusted effects per category | Nominal outcome; continuous or categorical predictors | Independent observations; IIA; linear logit; enough events per category | Stroke discharge destination (home / rehab / nursing) by age and NIHSS |

## Learning objectives

- Recognize an unordered multi-category outcome and choose multinomial rather than ordinal or binary regression
- Run `nnet::multinom` in R and `MNLogit` in Python, and compute Wald z and p in R
- Interpret relative risk ratios (exp β) against the reference category
- Explain how the model differs from ordinal (proportional-odds) regression and why that matters

## What it is

A multinomial logistic model is a **set of binary logistic regressions run simultaneously**, each comparing one outcome category with a common **reference category**. For an outcome with J categories there are J − 1 equations:

log( P(category j) / P(reference) ) = β₀ⱼ + β₁ⱼ·X₁ + β₂ⱼ·X₂ + …

Each category has its **own** intercept and its **own** slope for every predictor. Coefficients are estimated together by maximum likelihood so that the predicted probabilities across all categories add to 1.

- **Test statistics:** Wald z = β̂ / SE(β̂) for each coefficient; the likelihood-ratio χ² (null deviance − residual deviance, df = (J − 1) × number of predictors) for the whole model
- **Effect size:** **exp(β) = relative risk ratio (RRR)**: the multiplicative change, per unit of the predictor, in the **odds of category j versus the reference**, with a 95% CI. Despite the name, it is not a risk ratio in the epidemiologic sense; it is an odds-type ratio between two outcome categories.

## The intuition

Imagine a triage nurse sorting stroke patients into three exits: home, rehab, or a nursing facility. Older age and a more severe deficit each push patients toward the busier exits, but **they do not necessarily push them equally hard** toward rehab and toward nursing care. A multinomial model gives each exit its own slope, so the story can be "NIHSS barely matters for rehab versus home, but it strongly separates nursing from home."

The dot plot shows admission NIHSS by destination. Home patients are clustered at low scores, rehab sits slightly higher, and the nursing-facility group has the highest scores and the widest spread. The model turns that visible tilt into ratios with confidence intervals.

<div class="viz" data-viz="dots" data-unit="NIHSS at admission" data-groups='{"Home":[5,12,2,7,8,3,3,6,8,7,23,15,7,5,5,5,7,3,13,11,6,9,8,8,12,4,11,5,4,6],"Rehab":[5,6,6,11,17,3,9,6,9,6,7,2,6,4,11,8,5,9,10,10,8,4,6,10,9,8,11,11,3,11],"Nursing":[20,7,5,10,8,22,11,2,11,9,7,13,11,9,12,21,25,7,7,5,13]}' data-caption="Admission NIHSS by discharge destination (first 30 patients per group; all 21 in the nursing group). Mean NIHSS: 7.2 home, 7.6 rehab, 11.2 nursing (all patients)."></div>

The overall likelihood-ratio test of the model against an intercept-only model gives χ² = 24.33 on 4 df. It sits far in the upper tail of the χ² curve:

<div class="viz" data-viz="chisq" data-stat="24.33" data-df="4"></div>

## When to use it

Use multinomial logistic regression when:

- The outcome has **three or more categories with no natural order** (discharge destination, type of infection, tumor subtype, drug class prescribed)
- Or the outcome is ordered but the **proportional-odds assumption clearly fails** (see below)
- You want **adjusted** effects for each category against a clinically sensible reference
- Observations are independent, with **enough patients in the smallest category** (about 10 events per predictor per non-reference category is a common guide)

**Clinical scenario.** A stroke unit reviews 150 consecutive patients and records age, admission NIHSS score, and discharge destination: home (n = 94), inpatient rehabilitation (n = 35) or a nursing facility (n = 21). *How do age and stroke severity change the chance of each destination compared with going home?*

**Don't use it when:**

- The outcome is **binary** → use [logistic regression](#/test/logistic-regression)
- The categories are **ordered** (mild / moderate / severe) and the effect of a predictor is similar across the cut-points → use [ordinal regression](#/test/ordinal-regression), which is more powerful and more parsimonious
- Patients contribute **repeated or clustered** outcomes → use [GEE](#/test/gee) or a [mixed model](#/test/mixed-models)
- The outcome is a **count** → use [Poisson regression](#/test/poisson-regression)
- You only have one categorical predictor and no adjustment → a [chi-square test](#/test/chi-square-independence) is enough

## Assumptions

1. **Independent observations.** One row per patient.
2. **Outcome categories are mutually exclusive and exhaustive.** Each patient falls into exactly one category.
3. **Linearity in the logit.** Each continuous predictor must be linear with the log-odds of each category versus the reference. Check with splines or by adding a squared term.
4. **Independence of irrelevant alternatives (IIA).** The odds of rehab versus home should not depend on whether nursing care is available as an option. If one category is a near-substitute for another (e.g., two similar facilities), the model can mislead. Hausman–McFadden or Small–Hsiao tests exist, but they are unreliable; judge IIA on clinical grounds.
5. **Enough events in each category** and no (quasi-)complete separation. Sparse categories give huge SEs; consider merging categories or [Firth's penalized regression](#/test/firth-logistic).
6. **No severe multicollinearity** among predictors.

## Hypotheses

For each predictor j in each non-reference category k:

- **H₀:** βⱼₖ = 0, i.e. RRR = 1 (the predictor does not change the odds of category k versus the reference)
- **H₁:** βⱼₖ ≠ 0, i.e. RRR ≠ 1 (two-sided)

Overall model: H₀: all slope coefficients are 0 (likelihood-ratio test, df = 2 categories × 2 predictors = 4).

## Example data

Simulated data for 150 stroke patients. Counts by destination, with mean age and NIHSS:

| Destination | n | Mean age (years) | Mean NIHSS |
|---|---|---|---|
| Home (reference) | 94 | 73.4 | 7.2 |
| Rehab | 35 | 77.1 | 7.6 |
| Nursing facility | 21 | 77.8 | 11.2 |

All 150 rows (age, NIHSS and the destination code H / R / N) are typed inline in each code block below.

## How to run it

### R

```r
library(nnet)
# 150 stroke patients: age (years), admission NIHSS (0-42), discharge destination
v <- function(s) scan(text = s, quiet = TRUE)
dat <- data.frame(
  age   = v("73 64 79 83 77 77 75 82 62 84 66 77 81 76 60 74 77 83 79 79 84
        68 68 74 77 87 66 66 81 70 72 79 68 67 72 70 68 82 70 75 75 57
        77 88 74 67 86 76 83 76 80 87 80 77 68 79 70 94 77 67 89 70 60
        65 74 70 78 72 57 71 90 70 69 73 89 66 81 68 71 71 81 79 84 84
        93 75 51 69 67 71 84 69 67 78 80 72 78 78 77 69 81 54 80 78 92
        76 75 85 81 74 74 75 64 85 69 71 80 70 73 88 68 72 88 71 69 76
        80 68 73 73 65 74 89 82 67 77 72 73 82 79 82 78 67 70 63 80 86
        59 75 79"),
  nihss = v("5 20 12 5 2 6 7 8 3 3 6 8 7 7 23 15 7 5 6 5 11 5 7 17 5 3 13 10
        11 6 9 8 8 8 12 4 3 9 11 6 5 4 6 1 9 2 3 8 6 22 2 11 3 2 2 7 7 2
        11 3 9 6 4 8 11 1 8 5 3 5 7 5 9 3 13 10 11 9 10 2 8 3 16 4 6 7
        13 10 7 12 12 7 7 21 4 25 9 11 8 13 11 14 11 4 11 7 13 2 4 12 6
        3 15 11 10 11 7 4 4 4 8 5 6 5 7 8 10 3 14 4 9 7 7 4 5 6 10 4 5 8
        4 13 4 13 8 20 4 8 12 5"),
  dest  = strsplit(gsub("\\s", "", "HNHRHRHHHHHHHNHHHHRHRHHRNHHNHHHHNHHHRRHRHHHHRHHHRNHNHHNRHRNHNRRH
        RHRHHRNHRHNRNNRHRHHRRHHRHHNHHNHNRHRHRHHHRHHHHHHRHRRHHHHRRRHHHRHH
        HHHNNHHHHHNHHNHHHHHHHH"), "")[[1]])
dat$dest <- factor(dat$dest, levels = c("H", "R", "N"),
                   labels = c("Home", "Rehab", "Nursing"))   # Home = reference
dat$age10 <- dat$age / 10                                     # effects per 10 years
table(dat$dest)

# 1. Fit the multinomial model (trace = FALSE hides the iteration log)
fit  <- multinom(dest ~ age10 + nihss, data = dat, trace = FALSE)
fit0 <- multinom(dest ~ 1, data = dat, trace = FALSE)

# 2. Wald z and p (multinom does not print them), then RRRs with 95% CIs
b <- coef(fit); se <- summary(fit)$standard.errors
z <- b / se; p <- 2 * pnorm(-abs(z))
tab <- data.frame(vs_Home = rep(rownames(b), ncol(b)), term = rep(colnames(b), each = nrow(b)),
                  est = c(b), se = c(se), z = c(z), p = c(p),
                  RRR = exp(c(b)), lo = exp(c(b - 1.96 * se)), hi = exp(c(b + 1.96 * se)))
print(tab[tab$term != "(Intercept)", ], digits = 3, row.names = FALSE)

# 3. Overall likelihood-ratio test against the intercept-only model
anova(fit0, fit)
```

**What the code does**

- `factor(..., levels = c("Home", "Rehab", "Nursing"))` makes **Home the reference**, because the first level is the reference. Change the order to change the comparison.
- `age10 <- age / 10` expresses the age effect **per 10 years**, which is easier to interpret than per year.
- `multinom()` from the `nnet` package fits the model with a neural-network optimizer. `trace = FALSE` hides its iteration log. It prints coefficients and SEs but **no z or p values**, so we compute them: z = β / SE and p = 2·pnorm(−|z|).
- `exp(est ± 1.96·se)` converts coefficients and Wald limits to RRRs with 95% CIs.
- `anova(fit0, fit)` is the likelihood-ratio test of the whole model against the intercept-only model.

**Output**

```text

   Home   Rehab Nursing 
     94      35      21 
 vs_Home  term    est     se     z        p  RRR   lo   hi
   Rehab age10 0.6715 0.2749 2.442 0.014593 1.96 1.14 3.35
 Nursing age10 1.0164 0.3644 2.789 0.005282 2.76 1.35 5.64
   Rehab nihss 0.0476 0.0506 0.941 0.346637 1.05 0.95 1.16
 Nursing nihss 0.2043 0.0555 3.679 0.000235 1.23 1.10 1.37
          Model Resid. df Resid. Dev   Test    Df LR stat.      Pr(Chi)
1             1       298   272.3069           NA       NA           NA
2 age10 + nihss       294   247.9765 1 vs 2     4 24.33041 6.857204e-05
```

### Python

```python
import numpy as np, pandas as pd
import statsmodels.api as sm

# 150 stroke patients: age (years), admission NIHSS, discharge destination
v = lambda s: np.array(s.split(), dtype=float)
age   = v("73 64 79 83 77 77 75 82 62 84 66 77 81 76 60 74 77 83 79 79 "
          "84 68 68 74 77 87 66 66 81 70 72 79 68 67 72 70 68 82 70 75 "
          "75 57 77 88 74 67 86 76 83 76 80 87 80 77 68 79 70 94 77 67 "
          "89 70 60 65 74 70 78 72 57 71 90 70 69 73 89 66 81 68 71 71 "
          "81 79 84 84 93 75 51 69 67 71 84 69 67 78 80 72 78 78 77 69 "
          "81 54 80 78 92 76 75 85 81 74 74 75 64 85 69 71 80 70 73 88 "
          "68 72 88 71 69 76 80 68 73 73 65 74 89 82 67 77 72 73 82 79 "
          "82 78 67 70 63 80 86 59 75 79")
nihss = v("5 20 12 5 2 6 7 8 3 3 6 8 7 7 23 15 7 5 6 5 11 5 7 17 5 3 13 "
          "10 11 6 9 8 8 8 12 4 3 9 11 6 5 4 6 1 9 2 3 8 6 22 2 11 3 2 "
          "2 7 7 2 11 3 9 6 4 8 11 1 8 5 3 5 7 5 9 3 13 10 11 9 10 2 8 "
          "3 16 4 6 7 13 10 7 12 12 7 7 21 4 25 9 11 8 13 11 14 11 4 11 "
          "7 13 2 4 12 6 3 15 11 10 11 7 4 4 4 8 5 6 5 7 8 10 3 14 4 9 "
          "7 7 4 5 6 10 4 5 8 4 13 4 13 8 20 4 8 12 5")
dest  = np.array(list("HNHRHRHHHHHHHNHHHHRHRHHRNHHNHHHHNHHHRRHRHHHHRHHHRNHNHHNRHRNH "
          "NRRHRHRHHRNHRHNRNNRHRHHRRHHRHHNHHNHNRHRHRHHHRHHHHHHRHRRHHHHR "
          "RRHHHRHHHHHNNHHHHHNHHNHHHHHHHH".replace(" ", "")))
y = pd.Series(dest).map({"H": 0, "R": 1, "N": 2}).values   # 0 = Home (reference)
X = sm.add_constant(pd.DataFrame({"age10": age / 10, "nihss": nihss}))
print(pd.Series(y).value_counts().sort_index().to_string())

# 1. Fit the multinomial logit; MNLogit uses the lowest code (0) as reference
fit = sm.MNLogit(y, X).fit(disp=0)
print(fit.summary())

# 2. Relative risk ratios with 95% CIs (each vs Home)
rrr = pd.concat({"RRR": np.exp(fit.params), "lo": np.exp(fit.params - 1.96 * fit.bse),
                 "hi": np.exp(fit.params + 1.96 * fit.bse)}, axis=1)
rrr.columns = pd.MultiIndex.from_product([["RRR", "lo", "hi"], ["Rehab", "Nursing"]])
print(rrr.drop(index='const').round(3))

# 3. Overall likelihood-ratio test against the intercept-only model
print(f"LR chi2 = {fit.llr:.2f}, df = {int(fit.df_model)}, p = {fit.llr_pvalue:.2e}")
```

**What the code does**

- The destination is coded 0 = Home, 1 = Rehab, 2 = Nursing. **`MNLogit` uses the lowest code as the reference.**
- `sm.add_constant` adds the intercept; `MNLogit` does not add one automatically.
- The summary prints one block per non-reference category (`y=1` is Rehab, `y=2` is Nursing) with coefficients, SEs, z and p already computed.
- RRRs and CIs are `exp()` of the coefficients and Wald limits. `fit.llr` and `fit.llr_pvalue` give the overall likelihood-ratio test.

**Output**

```text
0    94
1    35
2    21
                          MNLogit Regression Results                          
==============================================================================
Dep. Variable:                      y   No. Observations:                  150
Model:                        MNLogit   Df Residuals:                      144
Method:                           MLE   Df Model:                            4
Date:                Sat, 17 Oct 2026   Pseudo R-squ.:                 0.08935
Time:                        16:31:33   Log-Likelihood:                -123.99
converged:                       True   LL-Null:                       -136.15
Covariance Type:            nonrobust   LLR p-value:                 6.857e-05
==============================================================================
       y=1       coef    std err          z      P>|z|      [0.025      0.975]
------------------------------------------------------------------------------
const         -6.3958      2.200     -2.908      0.004     -10.707      -2.085
age10          0.6717      0.275      2.443      0.015       0.133       1.211
nihss          0.0476      0.051      0.941      0.347      -0.052       0.147
------------------------------------------------------------------------------
       y=2       coef    std err          z      P>|z|      [0.025      0.975]
------------------------------------------------------------------------------
const        -11.0189      3.019     -3.650      0.000     -16.936      -5.102
age10          1.0164      0.364      2.789      0.005       0.302       1.731
nihss          0.2043      0.056      3.679      0.000       0.095       0.313
==============================================================================
         RRR             lo             hi        
       Rehab Nursing  Rehab Nursing  Rehab Nursing
age10  1.958   2.763  1.142   1.353  3.355   5.644
nihss  1.049   1.227  0.950   1.100  1.158   1.368
LR chi2 = 24.33, df = 4, p = 6.86e-05
```

## Interpreting the output

Each row compares a destination with **going home**, holding the other predictor fixed.

- **Age, rehab vs home: RRR = 1.96 per 10 years (95% CI 1.14 to 3.35), p = 0.015.** Each additional decade roughly doubles the odds of going to rehab rather than home.
- **Age, nursing vs home: RRR = 2.76 per 10 years (1.35 to 5.64), p = 0.005.** Age has an even stronger association with a nursing-facility discharge.
- **NIHSS, rehab vs home: RRR = 1.05 per point (0.95 to 1.16), p = 0.35.** No clear association: the CI includes 1.
- **NIHSS, nursing vs home: RRR = 1.23 per point (1.10 to 1.37), p = 0.0002.** Each additional NIHSS point increases the odds of a nursing facility by 23%. Per 5 points: exp(5 × 0.2043) = 2.78.
- **Likelihood-ratio test: χ² = 24.33, df = 4, p < 0.001.** Age and NIHSS together improve the model over the intercept only.
- **R and Python agree** on all RRRs, CIs and p-values to the displayed precision. The tiny difference in coefficients (e.g., age, rehab: 0.6715 vs 0.6717) comes from the optimizers' stopping rules.
- **To compare rehab with nursing**, subtract the coefficients: for age, exp(1.0164 − 0.6715) = 1.41, so the odds of nursing (versus rehab) are 41% higher per decade. For a formal CI, re-fit with another reference category.
- **Odds, not probabilities.** RRRs exaggerate risk ratios for common outcomes. To communicate the absolute picture, present **predicted probabilities** for typical patients (R: `predict(fit, newdata, type = "probs")`).

**How to report it:**
> In a multinomial logistic model for discharge destination (reference: home), each additional 10 years of age was associated with higher odds of rehabilitation (RRR 1.96, 95% CI 1.14–3.35) and of a nursing facility (RRR 2.76, 1.35–5.64). Each additional NIHSS point was associated with higher odds of a nursing facility (RRR 1.23, 1.10–1.37) but not of rehabilitation (RRR 1.05, 0.95–1.16). The model improved on the null model (likelihood-ratio χ²(4) = 24.33, p < 0.001).

## How this differs from ordinal regression

| | Multinomial | Ordinal (proportional odds) |
|---|---|---|
| Outcome | Unordered categories | Ordered categories |
| Slopes per predictor | **One per non-reference category** (J − 1) | **One in total** |
| Intercepts | J − 1 | J − 1 cut-points |
| Reports | RRR for each category vs a reference | A single cumulative OR across all cut-points |
| Key assumption | IIA | Proportional odds |
| Parameters | More; needs more data | Fewer; more power if the ordering is real |

Our destinations could be argued to be ordered by dependency (home < rehab < nursing), but the slopes are **not** the same: NIHSS is strongly tied to nursing and barely to rehab. That is exactly the kind of pattern a proportional-odds model would flatten into one number, so the multinomial model is the honest choice here. See [ordinal regression](#/test/ordinal-regression).

## Check your understanding

<details class="check"><summary>1. Which model fits each outcome? (a) Pain rated none / mild / moderate / severe; (b) type of infection: bacterial / viral / fungal; (c) readmitted within 30 days: yes / no.</summary>

(a) is ordered → [ordinal regression](#/test/ordinal-regression). (b) is unordered with three categories → multinomial logistic regression. (c) is binary → [logistic regression](#/test/logistic-regression).
</details>

<details class="check"><summary>2. The RRR for age (nursing vs home) is 2.76 per 10 years. A colleague says "a 10-year-older patient is 2.76 times as likely to go to a nursing facility." What is wrong?</summary>

The RRR compares **odds of nursing versus home**, not the probability of nursing among all patients. It is 2.76 times the odds of nursing relative to home, and because it ignores the rehab category and the baseline odds, it cannot be read as "times as likely" unless the outcome is rare. Report predicted probabilities for concrete patients to give the absolute picture.
</details>

<details class="check"><summary>3. What changes if we switch the reference category from Home to Nursing?</summary>

The model fit (log-likelihood, LR χ² and predicted probabilities) is **identical**, but the coefficients change meaning. You would get Home vs Nursing and Rehab vs Nursing, and the RRRs would be inverted or re-based: the age RRR for Home vs Nursing would be 1 / 2.76 = 0.36, and Rehab vs Nursing would be 1 / 1.41 = 0.71 per decade.
</details>

<details class="check"><summary>4. A patient's NIHSS is 5 points higher than another's, same age. Using the output, what is the multiplier on the odds of nursing vs home?</summary>

exp(5 × 0.2043) = exp(1.02) ≈ 2.78. Odds of a nursing-facility discharge versus home are about 2.8 times higher. The multiplier for rehab versus home would be exp(5 × 0.0476) = 1.27, which is compatible with no effect given the wide CI.
</details>

## Common mistakes

- **Using ordinal regression when slopes differ by category,** or multinomial when the outcome is ordered and a single slope would do (losing power).
- **Forgetting the reference category.** Every RRR is relative to it; state it in the text and in tables.
- **Reading RRR as a risk ratio** or a probability, especially for common outcomes.
- **Running separate binary logistic regressions** on subsets (home vs rehab, home vs nursing). It is roughly similar but uses different subsets of patients and does not constrain probabilities to sum to 1.
- **Too many predictors for the smallest category.** With 21 nursing patients, two predictors is about the limit.
- **Trusting the IIA test** rather than thinking about whether categories are close substitutes.
- **Ignoring that `multinom` prints no p-values;** compute them from z, or use another package.

## Related tests

- [Logistic regression](#/test/logistic-regression): the binary special case (J = 2)
- [Ordinal regression](#/test/ordinal-regression): for ordered categories, with one slope per predictor
- [Chi-square test of independence](#/test/chi-square-independence): unadjusted association of two categorical variables
- [Firth logistic regression](#/test/firth-logistic): when a category is so sparse that estimates blow up
- [Poisson regression](#/test/poisson-regression): count outcomes

## Key takeaways

- Multinomial regression = J − 1 simultaneous logistic regressions against one **reference category**, each with its own coefficients.
- exp(β) is a **relative risk ratio**: the multiplicative change in the odds of category j versus the reference per unit of the predictor.
- Predictors can matter for one category and not another (NIHSS: nursing yes, rehab no).
- Unlike ordinal regression, there is **no single OR across categories**; use ordinal only if the order is real and slopes are similar.
- In R, `multinom` gives no p-values; compute Wald z, or check with a likelihood-ratio test.

## References

- Hosmer DW, Lemeshow S, Sturdivant RX. *Applied Logistic Regression.* 3rd ed. Wiley; 2013. Ch. 8.
- Venables WN, Ripley BD. *Modern Applied Statistics with S.* 4th ed. Springer; 2002 (the `nnet` package).
- statsmodels documentation: `statsmodels.discrete.discrete_model.MNLogit`
- R documentation: `?nnet::multinom`
