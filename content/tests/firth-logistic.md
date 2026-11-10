# Firth (penalized) logistic regression

> Fits a logistic model with a small **penalty on the likelihood** (Jeffreys prior) that removes small-sample bias and gives finite, sensible estimates even when a predictor perfectly or almost perfectly separates the outcome.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Binary outcome with few events or (quasi-)complete separation | Binary outcome; any predictors | Independent observations; correct model form; few events is the reason to use it | Adverse drug reaction in 40 patients, where every risk-allele carrier reacted |

## Learning objectives

- Recognize (quasi-)complete separation and sparse-data bias in ordinary logistic regression
- Explain how Firth's penalty works and what it changes
- Run `logistf` in R and reproduce its coefficients with a short Newton–Raphson in Python
- Report Firth results with profile-likelihood CIs and penalized likelihood-ratio p-values

## What it is

Ordinary logistic regression maximizes the likelihood. When a predictor lines up perfectly with the outcome (**complete separation**) or nearly so, with one outcome level in a cell having zero patients (**quasi-complete separation**), the likelihood keeps rising as the coefficient heads to infinity. The software stops at an arbitrary huge number with an enormous standard error, and sometimes warns that "fitted probabilities numerically 0 or 1 occurred".

Firth (1993) adds a penalty term, ½·log|I(β)|, where I(β) is the Fisher information matrix. This is the Jeffreys invariant prior. The penalized log-likelihood is

ℓ*(β) = ℓ(β) + ½ · log |X′WX|

The penalty goes to −∞ as any coefficient runs off to ±∞, so a **finite maximum always exists**. As a bonus, it removes the first-order small-sample bias of the maximum-likelihood estimate (which, in ordinary logistic regression, is biased away from zero when events are few).

- **Test statistic:** the **penalized likelihood-ratio** χ² for each coefficient (preferred); Wald z = β̂ / SE(β̂) also exists, but is less reliable in small samples
- **Confidence intervals:** **profile penalized-likelihood** CIs (preferred, asymmetric); Wald CIs are symmetric on the log-odds scale and can be poor near separation
- **Effect size:** **OR = exp(β)** with a 95% CI, as in ordinary logistic regression. The OR is shrunk toward 1 compared with the (infinite) maximum-likelihood value.

## The intuition

Suppose that 12 of 12 patients carrying a risk allele had a severe reaction, while only 3 of 28 non-carriers did. The ordinary-logistic question is "how large must the carrier coefficient be to make the observed data most likely?" The answer is "as large as possible: infinity", because a carrier probability of exactly 1 explains 12 of 12 perfectly. A zero cell has made the data look more certain than 12 patients can justify.

Firth's correction acts like adding half a pseudo-observation to each cell: a mild prior that says "no probability is truly exactly 0 or 1." The estimate becomes large but finite (here log OR = 5.16, OR ≈ 174), and the CI honestly shows that the odds ratio could be anywhere from about 16 up to many thousands.

Age does not separate anything here, as the dot plot shows: ages overlap widely among the carriers, the non-carriers who reacted, and those who did not. The separation comes entirely from the carrier status.

<div class="viz" data-viz="dots" data-unit="age (years)" data-groups='{"Carrier, reaction":[68,56,45,49,58,78,62,42,53,46,47,65],"Non-carrier, reaction":[54,81,54],"Non-carrier, no reaction":[69,60,59,75,53,61,92,54,76,42,61,53,43,52,56,52,69,51,60,61,54,75,56,80,56]}' data-caption="Age by carrier status and outcome. All 12 carriers reacted (quasi-complete separation); age overlaps in every group."></div>

The penalized likelihood-ratio test for the carrier effect gives χ² = 28.65 on 1 df. It sits far into the upper tail of the χ² curve:

<div class="viz" data-viz="chisq" data-stat="28.65" data-df="1"></div>

## When to use it

Use Firth regression when:

- There is **complete or quasi-complete separation**, or ordinary logistic regression gives enormous coefficients and SEs
- **Events are few** (for example, fewer than about 10 events per predictor), or the outcome or an exposure is **rare**
- You want **less biased** odds ratios from a small sample
- You need a **pre-specified, defensible** analysis for a rare adverse event, pharmacogenetic association or small pilot

**Clinical scenario.** A pharmacogenetic study follows 40 patients starting a new drug. The team records age, whether the patient carries a risk allele (12 carriers), and whether a severe cutaneous adverse reaction occurred (15 patients). *Does the allele predict the reaction, and by how much?* All 12 carriers reacted, so ordinary logistic regression cannot give a finite odds ratio.

**Don't use it when:**

- There are plenty of events and no separation → ordinary [logistic regression](#/test/logistic-regression) is fine (Firth gives nearly the same answer)
- You have a 2 × 2 table with one predictor and no adjustment → [Fisher's exact test](#/test/fisher-exact) may be enough
- The outcome has more than two unordered categories → [multinomial regression](#/test/multinomial-logistic) (Firth versions exist)
- Observations are clustered or matched → [conditional logistic regression](#/test/conditional-logistic), [GEE](#/test/gee) or a [mixed model](#/test/mixed-models)
- You want to overcome genuinely **too few events** by adding many predictors. Firth fixes separation and bias, but it cannot create information

## Assumptions

1. **Binary outcome and independent observations.** One row per patient.
2. **The logistic model form is right:** linear on the log-odds scale and with no important omitted predictors.
3. **Separation is a feature of sparse data, not a coding error.** Check that the perfectly predicting variable is not a leak (for instance, the outcome itself recorded as a predictor).
4. **Keep the model small.** Firth's penalty reduces bias but does not cure overfitting; limit the number of predictors, and pre-specify them.
5. **Use profile-likelihood CIs and LR tests** for inference. Wald intervals can be unreliable in precisely the situations where you need Firth.

## Hypotheses

For each predictor j:

- **H₀:** βⱼ = 0, i.e. OR = 1 (no association with the outcome, adjusted for the other predictors)
- **H₁:** βⱼ ≠ 0, i.e. OR ≠ 1 (two-sided), tested with the penalized likelihood-ratio test

## Example data

Simulated data for 40 patients: age, risk-allele carrier status (12 carriers) and a severe adverse reaction (15 reactions).

| | No reaction | Reaction | Reaction rate |
|---|---|---|---|
| Non-carrier (n = 28) | 25 | 3 | 10.7% |
| Carrier (n = 12) | 0 | 12 | 100% |

The zero in the carrier row is the problem. The 40 ages and outcomes are typed inline in each code block below (carriers are the first 12 patients).

## How to run it

### R

```r
library(logistf)
# 40 patients on a new drug: age, risk-allele carrier (1 = yes), severe adverse reaction (1 = yes)
bits <- function(s) as.integer(strsplit(gsub("\\s", "", s), "")[[1]])
dat <- data.frame(
  age     = scan(text = "68 56 45 49 58 78 62 42 53 46 47 65 54 69 60 59 75 53 61 92 54
        76 42 61 53 43 52 81 56 52 69 51 60 61 54 54 75 56 80 56", quiet = TRUE),
  carrier = rep(c(1, 0), c(12, 28)),
  adr     = bits("1111111111111000000000000001000000100000"))
table(carrier = dat$carrier, adr = dat$adr)    # every carrier had a reaction

# 1. Ordinary logistic regression: separation -> huge estimate and SE
fit <- glm(adr ~ carrier + age, data = dat, family = binomial)
summary(fit)$coefficients

# 2. Firth penalized logistic regression with profile-likelihood CIs
fth <- logistf(adr ~ carrier + age, data = dat)
summary(fth)
round(exp(cbind(OR = fth$coefficients, fth$ci.lower, fth$ci.upper)), 3)
```

**What the code does**

- The data are the same 40 patients as in Python. `bits()` converts a digit string to 0/1.
- `glm(..., family = binomial)` is the **ordinary** logistic regression. Look at the carrier row: the estimate is about 22 on the log-odds scale with an SE in the thousands. That is the signature of separation. (R may or may not print a "fitted probabilities numerically 0 or 1" warning; here it did not, so do not rely on warnings.)
- `logistf()` fits the Firth model. Its `summary()` shows the penalized coefficients, **profile-likelihood** 95% limits (`method 2`), penalized LR χ² and p-values.
- `exp(cbind(coef, ci.lower, ci.upper))` gives ORs with profile-likelihood CIs.

**Output**

```text
       adr
carrier  0  1
      0 25  3
      1  0 12
               Estimate   Std. Error      z value  Pr(>|z|)
(Intercept) -3.06770067 3.185122e+00 -0.963134579 0.3354800
carrier     21.78853021 3.098354e+03  0.007032293 0.9943891
age          0.01531214 4.982895e-02  0.307293985 0.7586196
logistf(formula = adr ~ carrier + age, data = dat)

Model fitted by Penalized ML
Coefficients:
                   coef   se(coef) lower 0.95 upper 0.95      Chisq
(Intercept) -3.00216970 2.63368106 -8.9120264  2.8769322  1.0800169
carrier      5.15633436 1.53679790  2.7698671 10.1373913 28.6471674
age          0.01796959 0.04066419 -0.0815504  0.1032844  0.1599247
                       p method
(Intercept) 2.986938e-01      2
carrier     8.684093e-08      2
age         6.892258e-01      2

Method: 1-Wald, 2-Profile penalized log-likelihood, 3-None

Likelihood ratio test=29.43935 on 2 df, p=4.0488e-07, n=40
Wald test = 17.11821 on 2 df, p = 0.0001917905                 OR                 
(Intercept)   0.050  0.000    17.760
carrier     173.527 15.957 25270.456
age           1.018  0.922     1.109
```

### Python

```python
import numpy as np, pandas as pd
from scipy import stats

# 40 patients on a new drug: age, risk-allele carrier (1 = yes), severe adverse reaction (1 = yes)
age = np.array("68 56 45 49 58 78 62 42 53 46 47 65 54 69 60 59 75 53 61 92 "
          "54 76 42 61 53 43 52 81 56 52 69 51 60 61 54 54 75 56 80 56".split(), dtype=float)
carrier = np.repeat([1.0, 0.0], [12, 28])
adr = np.array(list("1111111111111000000000000001000000100000".replace(" ", "")), dtype=float)
X = np.column_stack([np.ones(40), carrier, age])

def firth_logit(X, y, fix=(), maxit=100, tol=1e-8):
    """Firth (1993) logistic regression by Newton-Raphson. Coefficients in `fix` stay at 0."""
    free = [j for j in range(X.shape[1]) if j not in fix]
    beta = np.zeros(X.shape[1])
    for _ in range(maxit):
        p = 1 / (1 + np.exp(-X @ beta))
        W = p * (1 - p)
        info = X.T @ (X * W[:, None])                       # Fisher information
        h = W * np.einsum("ij,jk,ik->i", X, np.linalg.inv(info), X)   # hat-matrix diagonal
        U = X.T @ (y - p + h * (0.5 - p))                   # penalized score
        step = np.clip(np.linalg.solve(info[np.ix_(free, free)], U[free]), -5, 5)
        beta[free] += step
        if np.abs(step).max() < tol:
            break
    p = 1 / (1 + np.exp(-X @ beta)); W = p * (1 - p)
    info = X.T @ (X * W[:, None])
    pll = np.sum(y * np.log(p) + (1 - y) * np.log(1 - p)) + 0.5 * np.linalg.slogdet(info)[1]
    return beta, np.sqrt(np.diag(np.linalg.inv(info))), pll   # pll = penalized log-likelihood

# 1. Firth estimates, Wald z and Wald 95% CIs
b, se, pll = firth_logit(X, adr)
z = b / se; p = 2 * stats.norm.sf(abs(z))
out = pd.DataFrame({"coef": b, "se": se, "z": z, "p_wald": p,
                    "OR": np.exp(b), "lo": np.exp(b - 1.96 * se), "hi": np.exp(b + 1.96 * se)},
                   index=["Intercept", "carrier", "age"])
print(out.round(4))

# 2. Penalized likelihood-ratio test for carrier (carrier fixed at 0), as logistf does
_, _, pll0 = firth_logit(X, adr, fix=(1,))
lr = 2 * (pll - pll0)
print(f"carrier: penalized LR chi2 = {lr:.3f}, p = {stats.chi2.sf(lr, 1):.2e}")
```

**What the code does**

- No mainstream Python package offers Firth logistic regression with profile CIs, so `firth_logit()` implements the algorithm directly in about 25 lines.
- Each Newton–Raphson iteration computes the hat-matrix diagonal `h` and uses the **modified score** `X′(y − p + h(0.5 − p))`. The extra `h(0.5 − p)` term is the gradient of the Jeffreys penalty.
- `fix=(1,)` holds the carrier coefficient at 0 while the others are re-estimated, which gives the **penalized likelihood-ratio test** (twice the difference in penalized log-likelihoods, with the penalty always evaluated from the full model, as `logistf` does).
- The reported CIs here are **Wald** intervals (β ± 1.96·SE). They are shown for completeness; for profile-likelihood intervals use R `logistf` (or write a root-finding routine over the penalized log-likelihood).

**Output**

```text
             coef      se       z  p_wald        OR      lo         hi
Intercept -3.0022  2.7934 -1.0748  0.2825    0.0497  0.0002    11.8559
carrier    5.1563  1.6005  3.2217  0.0013  173.5272  7.5331  3997.2543
age        0.0180  0.0435  0.4130  0.6796    1.0181  0.9349     1.1088
carrier: penalized LR chi2 = 28.647, p = 8.68e-08
```

## Interpreting the output

- **Ordinary glm:** carrier coefficient = 21.8, SE = 3,098, p = 0.99. The estimate is arbitrary (it would grow without limit if the optimizer ran longer), and the Wald test is meaningless. Do not report it.
- **Firth carrier coefficient = 5.156, OR ≈ 174.** Both R and Python give 5.156 (to four decimals). Carriers have very high odds of a reaction. The estimate is finite because of the penalty, and it is deliberately smaller than the infinite ML value.
- **Profile-likelihood 95% CI for the carrier OR: 16.0 to 25,270 (R).** The lower limit is clearly above 1, so the association is real, but the data allow anything from a large to an astronomically large OR. Report the CI, not just the point estimate.
- **Penalized LR test: χ² = 28.65, p = 8.7 × 10⁻⁸** (R's `logistf`, and Python's `penalized LR chi2 = 28.647, p = 8.68e-08`). This test is trustworthy.
- **Python's Wald CI for the carrier OR is 7.5 to 3,997, and its Wald p is 0.0013.** The Wald interval is placed and shaped differently from the profile interval (16.0 to 25,270), and the Wald p-value understates the evidence compared with the LR test (p = 8.7 × 10⁻⁸). This is why profile-likelihood CIs and LR tests are preferred; use R's `logistf` for reporting.
- **Age: OR = 1.018 per year (profile CI 0.92 to 1.11), p = 0.69.** No evidence that age matters after accounting for carrier status. The age coefficient is 0.0180 in both languages; standard errors differ slightly (0.0407 R, 0.0435 Python) because the two outputs report different variance estimates.
- **Intercept = −3.002** in both languages: the log-odds of a reaction for a non-carrier at age zero, a number without clinical meaning. Center age if you need a usable intercept.

**How to report it:**
> Because all 12 allele carriers developed the reaction (quasi-complete separation), Firth penalized logistic regression was used. Carrying the risk allele was strongly associated with the reaction (OR 173.5, 95% profile-likelihood CI 16.0 to 25,270; penalized likelihood-ratio p < 0.001), adjusted for age (OR 1.02 per year, 95% CI 0.92 to 1.11).

## Check your understanding

<details class="check"><summary>1. An ordinary logistic regression reports an OR of 3 × 10⁹ for a genotype with a 95% CI from 0 to infinity, and no warning. What has happened, and what should you do?</summary>

This is (quasi-)complete separation: the genotype perfectly or nearly perfectly predicts the outcome, so the ML estimate does not exist. Check the cross-tabulation for a zero cell, then refit with Firth's penalized logistic regression (R `logistf`) and report profile-likelihood CIs. Do not rely on the glm warning alone, since glm may not warn.
</details>

<details class="check"><summary>2. Firth's estimate for the carrier OR is 174, with a CI from 16 to 25,270. A colleague says "the true OR is 174." Is that right?</summary>

No. The point estimate is only the penalized maximum; the CI shows how little the data pin the size down. A 12/12 result is compatible with the OR being 16 or being in the thousands. State the CI, and also report the absolute rates (12 of 12 carriers vs 3 of 28 non-carriers).
</details>

<details class="check"><summary>3. Suppose that two non-carriers who had reactions had instead not reacted, so only one non-carrier reacted. What would happen to the Firth estimate for the carrier, and to the ordinary glm?</summary>

The separation would get even more complete (12 of 12 carriers vs 1 of 28 non-carriers), so the glm would still diverge. The Firth OR would be expected to increase (a bigger contrast), but would remain finite, and the CI would stay very wide on the upper side. Separation is always handled by the penalty; it never makes Firth fail.
</details>

<details class="check"><summary>4. Python gives a Wald p = 0.0013 for the carrier effect, while the penalized LR p is 8.7 × 10⁻⁸. Which should you report, and why do they differ so much?</summary>

Report the **penalized likelihood-ratio** p. Near separation the log-likelihood is very asymmetric, so the normal (Wald) approximation, which assumes a symmetric bell shape around the estimate, is poor. The LR test uses the actual shape of the penalized likelihood. This is also why the Python code is useful for checking coefficients but R `logistf` is preferred for final inference.
</details>

## Common mistakes

- **Reporting the ordinary glm estimate with an OR of millions.** It is an optimizer artifact, not a result.
- **Waiting for a warning.** glm may fail silently; check SEs, cross-tabs and extreme fitted probabilities.
- **Using Wald CIs** with Firth near separation. Use profile penalized-likelihood CIs and LR tests.
- **Dropping the separating variable** because it "breaks the model". It is often the finding.
- **Adding many predictors** to a model with 15 events. Firth reduces bias, but it does not give you extra events.
- **Treating the OR as precise.** The CI is usually enormous; present the raw counts too.
- **Using Firth to "fix" a coding error** in which the outcome leaks into a predictor.

## Related tests

- [Logistic regression](#/test/logistic-regression): the ordinary maximum-likelihood model
- [Fisher's exact test](#/test/fisher-exact): the exact unadjusted test for a 2 × 2 table with zero cells
- [Multinomial logistic regression](#/test/multinomial-logistic): sparse categories also cause separation
- [Conditional logistic regression](#/test/conditional-logistic): matched data (a Firth version exists in `coxphf`)
- [Effect measures](#/test/effect-measures): odds ratio versus risk ratio

## Key takeaways

- **Separation** (a zero cell for a predictor) makes the ML odds ratio infinite; glm returns huge estimates and SEs, with or without a warning.
- **Firth's penalty** (½ log|I|) always gives finite estimates and also reduces small-sample bias.
- Use **profile penalized-likelihood CIs and LR tests**; Wald CIs and p-values can be misleading.
- The Firth OR is large but finite; the CI is the honest summary of how little the data determine it.
- R `logistf` is the reference; the Python Newton–Raphson here reproduces its coefficients and LR test, but only gives Wald CIs.

## References

- Firth D. Bias reduction of maximum likelihood estimates. *Biometrika.* 1993;80(1):27–38.
- Heinze G, Schemper M. A solution to the problem of separation in logistic regression. *Stat Med.* 2002;21(16):2409–2419.
- Puhr R, Heinze G, Nold M, Lusa L, Geroldinger A. Firth's logistic regression with rare events: accurate effect estimates and predictions? *Stat Med.* 2017;36(14):2302–2317.
- R documentation: `?logistf::logistf`
