# Conditional logistic regression

> Analyses **matched case–control (or other matched-set) data** by comparing each case only with its own matched controls, giving odds ratios that respect the matching and are not distorted by the nuisance variables used to match.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Matched case–control sets (1:1, 1:2, 1:k) or other strata with one event per set | Binary case status, matched sets; any predictors | Correct matching; independent sets; linear logit; exposure measured similarly in cases and controls | GI bleeding cases, each matched to 2 controls, with NSAID use and diabetes as exposures |

## Learning objectives

- Recognize data from a matched design and why the matched structure must be analysed
- Run `survival::clogit` in R and `ConditionalLogit` in Python and read the matched odds ratios
- Explain why ordinary logistic regression is biased for matched data
- Interpret a matched OR correctly and say which sets contribute information

## What it is

In a matched case–control study, each case is paired with one or more controls who share key characteristics (age band, sex, hospital, calendar period). Each matched **set** has its own baseline risk that depends on those matching factors. Conditional logistic regression removes the unknown baseline of each set from the likelihood by **conditioning**: for a set with one case and k controls, it asks "given that exactly one of these k + 1 people is the case, what is the probability it was the person who actually was, given everyone's exposures?"

P(case is person i | one case in the set) = exp(β·Xᵢ) / Σⱼ exp(β·Xⱼ)

The set-specific intercepts cancel, so they are never estimated. Only the slopes β are, by maximum likelihood on this conditional likelihood.

- **Test statistic:** Wald z = β̂ / SE(β̂) per coefficient; likelihood-ratio, Wald or score χ² for the whole model
- **Effect size:** **matched OR = exp(β) with 95% CI**, the odds of exposure in cases versus controls **within the same matched set** (equivalently, the odds ratio for being the case)
- **Equivalence:** for 1:1 pairs and one binary exposure, the conditional OR equals the ratio of discordant pairs (the [McNemar](#/test/mcnemar) odds ratio).

## The intuition

Matching is a way of comparing like with like. A case aged 70 is compared with controls aged 70 of the same sex, not with the whole control pool. The analysis must keep that comparison intact. **Within each set, the question is whether the case was more exposed than their own controls.** Sets where everyone has the same exposure (all exposed or all unexposed) say nothing about the exposure effect, so they contribute no information.

In our example, 50 sets each contain 1 case and 2 controls. Of the 50 sets, 39 contain a mix of NSAID exposure (informative) and 11 do not. Crudely, 30 of 50 cases (60%) were exposed versus 25 of 100 controls (25%), but the right analysis compares within sets.

| NSAID use in the set | Case exposed | Case unexposed |
|---|---|---|
| Neither control exposed | 19 sets | 9 sets (uninformative) |
| One control exposed | 9 sets | 10 sets |
| Both controls exposed | 2 sets (uninformative) | 1 set |

In 28 of the 39 informative sets the case was exposed; in the other 11 the case was unexposed. Under H₀ (NSAIDs unrelated to bleeding), the Wald z statistic would follow the standard normal curve; ours is z = 3.52, well into the tail.

<div class="viz" data-viz="z" data-stat="3.52"></div>

## When to use it

Use conditional logistic regression when:

- The design is a **matched case–control study** (1:1, 1:2, up to about 1:5), or a **matched cohort**, **case–crossover** or **stratified** study with exactly one event per stratum
- You want **adjusted** matched ORs for exposures that were *not* used for matching (matching variables themselves cannot be analysed as predictors)
- Matched sets are independent of each other

**Clinical scenario.** A pharmacovigilance team studies upper GI bleeding. Fifty hospitalised patients with a confirmed bleed (cases) are each matched to two controls of the same sex, age within 3 years and the same hospital ward. They record **NSAID use in the preceding 30 days** and **diabetes**. *After accounting for the matching, how much higher are the odds of NSAID use among cases?*

**Don't use it when:**

- The study is **unmatched** → use ordinary [logistic regression](#/test/logistic-regression), adjusting for the confounders directly
- You matched **frequency** (by broad strata) and have many events per stratum → logistic regression with the stratifying variables as covariates is usually fine
- The outcome is a **time to event** with censoring, not a matched binary event → use [Cox regression](#/test/cox-ph)
- Observations are **clustered but not matched** (patients within clinics) → use [GEE](#/test/gee) or a [mixed model](#/test/mixed-models)
- Matched **binary pairs with a single exposure and no adjustment** → [McNemar's test](#/test/mcnemar) is enough
- You want to **control confounding without matching** → consider [propensity scores](#/test/propensity-scores) or [CMH stratification](#/test/cmh)

## Assumptions

1. **The matching is real and recorded.** Each row must carry the correct set ID; one case per set (or a known number).
2. **Sets are independent of each other** and members of a set are comparable apart from exposure.
3. **Linearity in the logit** for continuous exposures. Check with splines or categories.
4. **Matched factors are not entered as predictors.** Their effect is conditioned away. You can examine effect modification (an interaction between exposure and a matching factor) but not main effects.
5. **Exposure is measured in the same way for cases and controls** (no recall bias, same time window). Matching does not fix information bias.
6. **Enough informative sets.** The effective sample size is the number of sets with **discordant** exposure, not the number of patients. Few discordant sets give wide CIs.
7. **No overmatching.** Do not match on a factor that is a consequence of the exposure or on the exposure's causal pathway.

## Hypotheses

For each exposure j:

- **H₀:** βⱼ = 0, i.e. matched OR = 1 (within matched sets, cases and controls are equally likely to be exposed)
- **H₁:** βⱼ ≠ 0, i.e. matched OR ≠ 1 (two-sided)

## Example data

Simulated data: 50 matched sets, each with **1 case and 2 controls** (150 rows), with two binary exposures. Crude exposure by case status:

| | Controls (n = 100) | Cases (n = 50) |
|---|---|---|
| NSAID use in last 30 days | 25 (25%) | 30 (60%) |
| Diabetes | 24 (24%) | 22 (44%) |

The rows are in set order: rows 1–3 are set 1 (case first), rows 4–6 are set 2, and so on. All 150 rows of the two exposures are typed inline below.

## How to run it

### R

```r
library(survival)
# 50 matched sets, each 1 case + 2 controls (rows 1-3 = set 1, 4-6 = set 2, ...)
bits <- function(s) as.integer(strsplit(gsub("\\s", "", s), "")[[1]])
dat <- data.frame(
  set    = rep(1:50, each = 3),
  case   = rep(c(1, 0, 0), times = 50),
  nsaid  = bits("000100101100001101100100101000101001000110000100100100000010000100
      010101100000100010010010110100100010100000100100100111000010100011
      110101100100010111"),
  diab   = bits("110000100100000100001101000101100100010000010100011010000000100010
      111100100001000100001010000001101000100101000100000010000110000000
      001101010100000101"))

# 1. Exposure by case status (crude 2 x 2 table, ignoring matching)
table(case = dat$case, nsaid = dat$nsaid)

# 2. Conditional logistic regression: strata() removes each set's own intercept
fit <- clogit(case ~ nsaid + diab + strata(set), data = dat)
summary(fit)

# 3. For contrast: ordinary logistic regression that ignores the matching
naive <- glm(case ~ nsaid + diab, data = dat, family = binomial)
round(exp(cbind(OR = coef(naive), confint.default(naive)))[-1, ], 3)
```

**What the code does**

- `bits()` turns a string of 0/1 digits into a vector. The case indicator repeats `1, 0, 0` for each set (case first, then two controls).
- `table()` shows the crude 2 × 2 table, which **ignores** the matching and is not the valid analysis.
- `clogit(case ~ nsaid + diab + strata(set))` is the conditional logistic model. `strata(set)` is the essential part: it identifies the matched sets. Do not put `set` in as an ordinary predictor.
- `clogit` is a wrapper around `coxph` with the exact partial likelihood, which is why the printed call says `coxph(...method = "exact")` and the summary says "number of events" and prints a concordance. Ignore the Cox labels, and read the coefficient table.
- The last block fits an **ordinary** logistic model and prints its ORs for comparison only.

**Output**

```text
    nsaid
case  0  1
   0 75 25
   1 20 30
Call:
coxph(formula = Surv(rep(1, 150L), case) ~ nsaid + diab + strata(set), 
    data = dat, method = "exact")

  n= 150, number of events= 50 

        coef exp(coef) se(coef)     z Pr(>|z|)    
nsaid 1.3171    3.7326   0.3742 3.520 0.000431 ***
diab  0.7834    2.1889   0.3975 1.971 0.048750 *  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

      exp(coef) exp(-coef) lower .95 upper .95
nsaid     3.733     0.2679     1.793     7.772
diab      2.189     0.4568     1.004     4.771

Concordance= 0.715  (se = 0.06 )
Likelihood ratio test= 19.44  on 2 df,   p=6e-05
Wald test            = 15.8  on 2 df,   p=4e-04
Score (logrank) test = 19.12  on 2 df,   p=7e-05

         OR 2.5 % 97.5 %
nsaid 4.394 2.100  9.195
diab  2.383 1.104  5.145
```

### Python

```python
import numpy as np, pandas as pd
import statsmodels.api as sm
from statsmodels.discrete.conditional_models import ConditionalLogit

# 50 matched sets, each 1 case + 2 controls (rows 1-3 = set 1, 4-6 = set 2, ...)
bits = lambda s: np.array(list(s.replace(" ", "")), dtype=int)
dat = pd.DataFrame({
    "set":   np.repeat(np.arange(1, 51), 3),
    "case":  np.tile([1, 0, 0], 50),
    "nsaid": bits("000100101100001101100100101000101001000110000100100100000010 "
          "000100010101100000100010010010110100100010100000100100100111 "
          "000010100011110101100100010111"),
    "diab":  bits("110000100100000100001101000101100100010000010100011010000000 "
          "100010111100100001000100001010000001101000100101000100000010 "
          "000110000000001101010100000101")})

# 1. Crude 2 x 2 table, ignoring the matching
print(pd.crosstab(dat["case"], dat["nsaid"]))

# 2. Conditional logistic regression: groups= identifies the matched sets
fit = ConditionalLogit(dat["case"], dat[["nsaid", "diab"]], groups=dat["set"]).fit(disp=0)
print(fit.summary())
print(pd.DataFrame({"OR": np.exp(fit.params), "lo": np.exp(fit.conf_int()[0]),
                    "hi": np.exp(fit.conf_int()[1])}).round(3))

# 3. For contrast: ordinary logistic regression that ignores the matching
naive = sm.Logit(dat["case"], sm.add_constant(dat[["nsaid", "diab"]])).fit(disp=0)
print(pd.DataFrame({"OR": np.exp(naive.params), "lo": np.exp(naive.conf_int()[0]),
                    "hi": np.exp(naive.conf_int()[1])}).drop("const").round(3))
```

**What the code does**

- `ConditionalLogit(endog, exog, groups=...)` takes the case indicator, the exposures and the set IDs. There is **no intercept** (it cancels out), so do not add a constant.
- The summary reports coefficients, SEs, z and p. We exponentiate them for matched ORs and 95% CIs.
- The last block uses `sm.Logit` ignoring the sets, for contrast.

**Output**

```text
nsaid   0   1
case         
0      75  25
1      20  30
                  Conditional Logit Model Regression Results                  
==============================================================================
Dep. Variable:                   case   No. Observations:                  150
Model:               ConditionalLogit   No. groups:                         50
Log-Likelihood:               -45.211   Min group size:                      3
Method:                          BFGS   Max group size:                      3
Date:                Sat, 17 Oct 2026   Mean group size:                   3.0
Time:                        16:31:35                                         
==============================================================================
                 coef    std err          z      P>|z|      [0.025      0.975]
------------------------------------------------------------------------------
nsaid          1.3171      0.374      3.520      0.000       0.584       2.050
diab           0.7835      0.398      1.971      0.049       0.004       1.563
==============================================================================
          OR     lo     hi
nsaid  3.733  1.793  7.772
diab   2.189  1.004  4.771
          OR     lo     hi
nsaid  4.394  2.100  9.195
diab   2.383  1.104  5.145
```

## Interpreting the output

- **NSAID use: matched OR = 3.73 (95% CI 1.79 to 7.77), z = 3.52, p < 0.001.** Within matched sets, the odds of NSAID use are about 3.7 times higher in cases than in their matched controls.
- **Diabetes: matched OR = 2.19 (95% CI 1.00 to 4.77), z = 1.97, p = 0.049.** Suggestive, but the CI lower limit is right at 1; with only 50 sets this is imprecise. Do not over-interpret a p-value just under 0.05.
- **Each OR is adjusted** for the other exposure, and **conditioned on** the matching factors (sex, age band, ward). No coefficient exists for the matching factors.
- **R and Python agree:** coefficients 1.3171 and 0.7834 (R) vs 1.3171 and 0.7835 (Python); SEs 0.374 and 0.398. The R summary also prints Cox-style "likelihood ratio test = 19.44 on 2 df", which is the valid overall test for the model.
- **Ordinary logistic regression gave 4.39 (2.10 to 9.20) for NSAIDs and 2.38 (1.10 to 5.15) for diabetes.** These do not answer the matched question. They treat the 150 patients as unrelated and treat the matching factors as if they had no effect on who becomes a case. The direction and size of the bias depend on how strongly the matching factors relate to exposure and outcome: classically the unmatched OR is **biased toward 1** when the matching factors are associated with the exposure, and the standard errors are wrong in either direction. Here the gap (4.39 vs 3.73) happens to go the other way because this is one small simulated dataset, so do not read the pattern; the point is that the unconditional analysis targets the wrong quantity.
- **Why a separate set intercept does not rescue ordinary logistic regression.** Adding a dummy variable per set is tempting but gives a **biased** estimate when sets are small: with 1:1 pairs the estimated OR is roughly the *square* of the true OR, because the number of nuisance parameters grows with the sample size (the incidental-parameter problem). Conditioning avoids estimating them.

**How to report it:**
> In a matched case–control analysis (50 cases, 100 controls; 1:2 matching on sex, age and ward), NSAID use in the previous 30 days was associated with GI bleeding (conditional OR 3.73, 95% CI 1.79–7.77, p < 0.001) after adjustment for diabetes. The matched OR for diabetes was 2.19 (95% CI 1.00–4.77).

## What matched ORs mean

- A matched OR is **conditional on the set**: it compares cases with *their own* controls, so it is an odds ratio for exposure (or equivalently for being the case) within groups of similar people. It is not a population-averaged OR.
- **Interpretation:** "Holding the matching factors fixed, the odds of being the case are 3.7 times as high for NSAID users." In a case–control study the OR approximates a relative risk or rate ratio only if the disease is rare (or controls are sampled by density sampling).
- **Matching factors cannot be estimated.** You cannot say whether older people bleed more from this analysis; that information was removed on purpose.
- **Absolute risk cannot be estimated** from case–control data: the proportion of cases is fixed by the design.

## Check your understanding

<details class="check"><summary>1. A study matches each of 80 cases to one control by age and sex. Which analysis fits for a single binary exposure: ordinary logistic regression, conditional logistic regression, or McNemar's test?</summary>

Either **McNemar's test** (unadjusted, with the OR from discordant pairs) or **conditional logistic regression** if you also want to adjust for other exposures. Ordinary logistic regression ignores the pairing and is not appropriate.
</details>

<details class="check"><summary>2. A colleague adds `age` and `sex` (the matching variables) to the conditional model to "control for them" and gets no estimates for them. Is something wrong?</summary>

No. Within a matched set, age and sex are constant (or nearly so), so they cancel from the conditional likelihood. They are controlled for by the design and cannot be estimated. If you want to know whether the exposure effect differs by sex, add an exposure × sex **interaction**.
</details>

<details class="check"><summary>3. In 1:1 pairs with one binary exposure, 40 pairs have only the case exposed and 10 pairs have only the control exposed. All other pairs are concordant. What is the matched OR, and what do the concordant pairs contribute?</summary>

The matched OR is 40 / 10 = 4.0. Concordant pairs (both exposed or both unexposed) contribute **nothing** to the estimate; only the 50 discordant pairs carry information. That is why the effective sample size can be much smaller than the number of patients.
</details>

<details class="check"><summary>4. In our 1:2 data, 11 of 50 sets had all three members with the same NSAID exposure. What happens to those sets in the conditional analysis?</summary>

They are **uninformative for NSAID** and contribute nothing to its coefficient, because nothing in the set distinguishes the case from the controls. They stay in the data (they might be informative for diabetes) but the NSAID estimate rests on the 39 mixed sets, which is why the CI is wide (1.79 to 7.77).
</details>

## Common mistakes

- **Using ordinary logistic regression on matched data.** It ignores the design; the estimate targets the wrong quantity and the SEs are not valid.
- **Fitting a dummy variable per set** in an ordinary model. In small sets this gives a biased OR (about squared for pairs).
- **Entering matching variables as predictors** and expecting estimates.
- **Mishandling missing data:** if a control is lost, keep the remaining set members (a 1:1 set with one member missing is dropped, since it has no comparison; a 1:2 set becomes 1:1). Do not drop the whole set unnecessarily.
- **Ignoring uninformative sets** when quoting sample size. Report the number of informative sets.
- **Overmatching** on variables affected by the exposure, which biases the OR toward 1.
- **Misreading Cox labels** in the `clogit` output (it says "events" and "concordance"; the odds ratio is `exp(coef)`).

## Related tests

- [Logistic regression](#/test/logistic-regression): the unconditional model for unmatched data
- [McNemar's test](#/test/mcnemar): the unadjusted test for 1:1 matched binary pairs
- [Cochran–Mantel–Haenszel](#/test/cmh): stratified analysis of 2 × 2 tables
- [Propensity scores](#/test/propensity-scores): an alternative way to control confounding
- [Effect measures](#/test/effect-measures): odds ratio versus risk ratio

## Key takeaways

- Matched sets must be analysed as sets: **conditional logistic regression** compares each case with its own controls.
- Set-specific intercepts are conditioned away, so **matching variables get no coefficient**; only exposures do.
- exp(β) is a **matched OR**, conditional on the set; sets with no exposure discordance contribute nothing.
- Ordinary logistic regression on matched data (or with a dummy per set) gives the wrong estimand and can be biased.
- In R use `survival::clogit(... + strata(set))`; in Python `statsmodels` `ConditionalLogit(..., groups=set)`.

## References

- Breslow NE, Day NE. *Statistical Methods in Cancer Research, Vol. 1: The Analysis of Case-Control Studies.* IARC; 1980.
- Rothman KJ, Greenland S, Lash TL. *Modern Epidemiology.* 3rd ed. Lippincott; 2008. Ch. 11 (matching).
- R documentation: `?survival::clogit`
- statsmodels documentation: `statsmodels.discrete.conditional_models.ConditionalLogit`
