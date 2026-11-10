# Propensity scores

> Summarizes many measured confounders into one number, the probability of receiving treatment given those covariates, and uses it (by matching, weighting, or stratification) to compare treated and untreated patients who look alike.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Observational comparison of treatments with several measured confounders | Binary treatment; any outcome | No unmeasured confounding; positivity (overlap); correct PS model | New therapy given to sicker patients, outcome function score at 6 months |

## Learning objectives

- Explain **confounding by indication** and why a naive comparison can have the wrong sign
- Estimate a propensity score and use it for 1:1 matching and inverse-probability weighting (IPTW)
- Check covariate balance with **standardized mean differences (SMDs)** before and after
- State what propensity scores cannot do: they do not fix unmeasured confounding, and need overlap (positivity)

## What it is

In a randomized trial, treatment is assigned by a coin, so the groups are alike on average. In observational data, doctors choose treatment, usually for reasons that also affect the outcome. This is **confounding by indication**: sicker patients get the new therapy and also do worse for reasons unrelated to it.

The **propensity score** is e(x) = P(treated | covariates x), usually estimated by logistic regression. Patients with the same propensity score have, on average, the same distribution of the measured covariates whether treated or not, so comparing them mimics a randomized comparison **with respect to the measured covariates**.

Common ways to use it:

- **Matching:** pair each treated patient with an untreated one with a similar score (here 1:1 nearest neighbour on the logit scale, with a **caliper** of 0.2 SD of the logit PS so poor matches are dropped). Estimates the effect in the treated (ATT).
- **Inverse-probability-of-treatment weighting (IPTW):** weight each patient by 1/e(x) if treated and 1/(1 − e(x)) if untreated, creating a pseudo-population in which treatment is unrelated to the covariates. Estimates the average effect in everyone (ATE).
- **Stratification or covariate adjustment** by the score are alternatives.
- **Balance diagnostic:** SMD = (mean₁ − mean₀) / pooled SD for each covariate. A common rule of thumb is |SMD| < 0.1 after adjustment. p-values are not used to judge balance.

## The intuition

Think of a clinic in which the sickest patients are steered toward the new treatment. If you compare average outcomes, the treated group looks worse, not because the treatment hurts, but because those patients were sicker to begin with. The propensity score answers: "how likely was this patient to be steered to treatment?" Two patients with the same answer, one treated and one not, differ mainly by chance, so comparing them is much fairer.

The picture below shows the estimated PS at 30 evenly spaced quantiles (minimum to maximum) of each group. The treated sit higher (they were more likely to be treated), but the two groups **overlap** over most of the range. That overlap is what makes the comparison possible. Where one group has no counterpart, no method can compare them.

<div class="viz" data-viz="dots" data-unit="estimated propensity score" data-groups='{"Untreated":[0.14,0.18,0.2,0.22,0.25,0.26,0.28,0.31,0.35,0.36,0.38,0.41,0.44,0.46,0.47,0.48,0.5,0.51,0.52,0.53,0.54,0.55,0.57,0.59,0.61,0.62,0.64,0.7,0.76,0.83],"Treated":[0.22,0.32,0.35,0.38,0.41,0.42,0.43,0.44,0.46,0.47,0.5,0.51,0.53,0.54,0.55,0.56,0.57,0.59,0.6,0.63,0.66,0.68,0.7,0.72,0.76,0.78,0.8,0.81,0.82,0.86]}' data-caption="Propensity scores (30 quantile points per group, 200 patients in total). Treated patients have higher scores on average, but the ranges overlap, so positivity holds. Untreated patients below about 0.25 have few treated counterparts."></div>

In our simulated data the true treatment effect is **+5 points** on a 6-month function score. A naive comparison gives **−4.0**: it suggests the treatment is harmful. Matching and weighting move the estimate to the positive side, close to the truth.

## When to use it

- The exposure is **binary** (treated vs not) and the data are **observational**
- You have many measured confounders, so outcome regression alone is awkward (especially with few outcome events)
- You can argue that the **important confounders were measured**

**Clinical scenario.** A registry of 200 patients with a chronic condition records age, diabetes, and a baseline severity score (0–100), whether the patient received a new rehabilitation therapy, and a function score at 6 months. Sicker patients were more likely to receive the therapy. *What is the effect of the therapy on the function score?*

**Don't use it when:**

- Treatment was randomized: a simple comparison or [ANCOVA](#/test/ancova) is the right analysis
- Key confounders were **not measured**: no propensity score can fix this
- Some patients have essentially **no chance** of treatment (or of no treatment); positivity is violated. Restrict the population instead.
- You only need to adjust for one or two covariates: [linear regression](#/test/linear-regression) or [logistic regression](#/test/logistic-regression) may be simpler
- You have a time-to-event outcome: combine weights with [Cox regression](#/test/cox-ph)

## Assumptions

1. **No unmeasured confounding (exchangeability).** Given the measured covariates, treatment assignment is as good as random. This is untestable: justify it with subject-matter knowledge and consider sensitivity analyses (such as the E-value).
2. **Positivity (overlap).** Every patient has a probability strictly between 0 and 1 of each treatment. Check PS distributions by group; extreme weights (say > 10) signal violations.
3. **Correct propensity-score model.** The PS model must include confounders (variables that affect both treatment and outcome). Check **balance**, not model fit.
4. **Consistency and no interference.** The treatment is well defined, and one patient's treatment does not affect another's outcome.

## Hypotheses

Propensity-score analysis is an estimation approach. The estimand is a causal effect, defined before analysis:

- **ATT:** mean outcome difference if the treated had been treated vs untreated (what matching estimates)
- **ATE:** the same difference in the whole population (what the IPTW weights above estimate)
- **H₀:** effect = 0 vs **H₁:** effect ≠ 0, tested with a robust standard error that accounts for the weighting or the matching

## Example data

A simulated registry of 200 patients (103 treated, 97 untreated). Outcome = function score at 6 months (higher is better). The **true therapy effect is +5** for everyone. Severity, age, and diabetes raise the chance of treatment and lower the outcome. Raw standardized mean differences between treated and untreated are large for severity (0.68), moderate for age (0.32) and diabetes (0.24).

## How to run it

### R

```r
library(MatchIt); library(sandwich); library(lmtest)

# 200 simulated patients; outcome y = 6-month function score; TRUE effect = +5
age <- scan(text = "
59 65 50 51 77 56 78 71 65 55 57 62 50 62 54 65 63 74 59 58 58 65 61
69 66 65 63 57 63 55 54 56 72 49 56 70 63 80 59 62 49 65 74 56 74 62
43 74 72 67 73 63 57 70 67 70 63 69 80 66 57 88 64 45 70 82 57 54 59
73 70 64 57 67 76 67 58 70 54 69 64 68 59 56 88 59 52 70 56 50 77 55
74 60 70 61 78 54 79 59 61 67 71 64 44 71 70 55 78 76 73 81 69 42 64
60 61 68 65 66 90 59 77 72 63 79 78 71 78 56 64 53 68 54 60 64 62 42
62 68 49 70 50 62 64 68 69 85 62 69 49 72 69 72 55 73 82 81 57 62 65
59 71 66 52 64 51 62 65 65 70 68 74 68 54 76 55 76 74 77 64 64 72 58
88 62 61 63 60 68 65 60 74 66 69 93 82 53 67 75", quiet = TRUE)
sev <- scan(text = "
30 44 39 77 84 59 62 42 77 54 24 51 35 41 64 45 51 67 44 52 58 41 39
47 53 38 47 67 57 54 33 43 16 48 54 62 21 32 81 16 22 59 59 39 67 41
22 43 51 48 55 46 61 51 49 64 60 54 48 60 39 75 72 55 30 51 48 55 41
26 60 75 34 23 78 45 63 29 46 48 41 48 72 61 65 61 74 69 52 67 77 57
48 84 58 56 57 27 53 48 69 35 70 8 64 60 66 33 29 84 34 71 80 66 51 43
41 37 55 59 34 59 57 63 59 36 49 38 53 52 52 61 65 25 59 56 27 36 48
57 35 69 55 50 74 71 60 76 35 50 25 79 34 84 67 52 59 46 24 23 64 48
50 54 40 46 23 31 47 43 51 24 60 59 55 44 59 84 58 61 52 74 34 53 83
71 80 46 47 38 61 51 51 28 24 63 74 65 86 85", quiet = TRUE)
y <- scan(text = "
87 68 80 51 33 74 63 83 55 84 86 80 83 82 59 80 78 57 67 56 74 85 76
70 58 77 71 71 75 66 89 74 101 95 60 54 74 99 72 84 84 65 72 77 70 83
99 77 66 78 75 66 79 68 67 49 69 78 77 73 51 62 65 84 84 58 70 69 79
95 64 52 88 88 64 75 62 78 82 65 87 63 44 67 48 64 53 52 72 59 53 76
70 61 65 69 61 99 74 72 65 76 59 106 60 64 64 98 74 58 89 52 56 76 74
80 76 92 64 69 75 69 55 72 76 74 57 75 62 60 74 63 57 93 76 69 81 80
73 70 87 62 66 74 43 59 70 38 81 77 91 52 83 57 76 56 56 80 96 90 55
83 65 75 75 80 96 87 79 90 73 69 72 79 79 74 72 38 72 56 81 64 79 73
48 46 38 71 69 70 64 71 82 66 87 55 52 69 48 42", quiet = TRUE)
chr <- function(s) as.integer(strsplit(gsub("\\s", "", s), "")[[1]])
diab  <- chr("
1110001010000000010011100100110100100000000000000000001100011000000000
0101100010101110101100010001010000111000011101000000000000001100000010
000011100000010010000000010001100110010100010110011101010100")
treat <- chr("
1100110111010101000001100101000111000110001000001101100111010111000001
1000001010000110110011010000011110011100011111100100011110111100110001
010011111101111011001101010010100000010110111101111010011111")
d <- data.frame(age, diab, sev, treat, y)

# 1. Naive comparison (confounded by indication)
print(round(coef(lm(y ~ treat, d))["treat"], 2))

# 2. Propensity score by logistic regression
ps <- fitted(glm(treat ~ age + diab + sev, binomial, d))
w_ate <- ifelse(d$treat == 1, 1/ps, 1/(1 - ps))        # IPTW (ATE) weights

# 3. 1:1 nearest-neighbour matching on the logit PS, caliper 0.2 SD
m <- matchit(treat ~ age + diab + sev, d, method = "nearest",
             distance = "glm", link = "linear.logit",
             caliper = 0.2, ratio = 1)
md <- match.data(m)
cat("Matched pairs:", sum(md$treat), "of", sum(d$treat), "treated\n")

# Standardized mean difference; denominator = pooled SD in the ORIGINAL data
smd <- function(v, dat, w = rep(1, nrow(dat))) {
  x <- dat[[v]]; t <- dat$treat
  sd0 <- sqrt((var(d[[v]][d$treat == 1]) + var(d[[v]][d$treat == 0])) / 2)
  (weighted.mean(x[t == 1], w[t == 1]) - weighted.mean(x[t == 0], w[t == 0])) / sd0
}
vars <- c("age", "diab", "sev")
print(round(cbind(raw     = sapply(vars, smd, dat = d),
                  matched = sapply(vars, smd, dat = md),
                  iptw    = sapply(vars, smd, dat = d, w = w_ate)), 3))

# 4. Effect estimates (matched: robust SE; IPTW: HC1 robust SE)
fm <- lm(y ~ treat, md, weights = weights)
print(round(coeftest(fm, vcovCL(fm, cluster = md$subclass))["treat", 1:2], 2))
fw <- lm(y ~ treat, d, weights = w_ate)
print(round(coeftest(fw, vcovHC(fw, "HC1"))["treat", 1:2], 2))
```

**What the code does**

- Data are typed in as compact vectors, shared with the Python tab, so both languages use the same patients.
- `glm(..., binomial)` gives the propensity score; `w_ate` are the IPTW weights.
- `matchit(..., method = "nearest", distance = "glm", caliper = 0.2)` does greedy 1:1 matching without replacement, dropping treated patients with no control within 0.2 SD of the distance. `match.data()` returns the matched sample with weights.
- `link = "linear.logit"` makes `MatchIt` match on the **logit** of the PS (the linear predictor). Without it, `distance = "glm"` matches on the **probability** scale, and the caliper becomes 0.2 SD of the probability instead. That gives different pairs.
- `smd()` computes standardized mean differences (raw, after matching, after weighting), always using the **original** pooled SD as the denominator.
- The matched estimate uses a regression with cluster-robust SE (pairs are clusters); the IPTW estimate uses a heteroskedasticity-robust (HC1) SE. Ordinary SEs are wrong for weighted data.

**Output**

```text
treat 
-4.01 
Matched pairs: 77 of 103 treated
       raw matched  iptw
age  0.324  -0.004 0.066
diab 0.240   0.056 0.000
sev  0.678   0.068 0.050
  Estimate Std. Error 
      3.32       1.38 
  Estimate Std. Error 
      2.84       2.04 
```

### Python

```python
import numpy as np, pandas as pd, statsmodels.api as sm
from sklearn.neighbors import NearestNeighbors

# Same 200 simulated patients as the R tab; TRUE effect = +5
nums = lambda s: np.array(s.split(), dtype=float)
age = nums("""
59 65 50 51 77 56 78 71 65 55 57 62 50 62 54 65 63 74 59 58 58 65 61
69 66 65 63 57 63 55 54 56 72 49 56 70 63 80 59 62 49 65 74 56 74 62
43 74 72 67 73 63 57 70 67 70 63 69 80 66 57 88 64 45 70 82 57 54 59
73 70 64 57 67 76 67 58 70 54 69 64 68 59 56 88 59 52 70 56 50 77 55
74 60 70 61 78 54 79 59 61 67 71 64 44 71 70 55 78 76 73 81 69 42 64
60 61 68 65 66 90 59 77 72 63 79 78 71 78 56 64 53 68 54 60 64 62 42
62 68 49 70 50 62 64 68 69 85 62 69 49 72 69 72 55 73 82 81 57 62 65
59 71 66 52 64 51 62 65 65 70 68 74 68 54 76 55 76 74 77 64 64 72 58
88 62 61 63 60 68 65 60 74 66 69 93 82 53 67 75
""")
sev = nums("""
30 44 39 77 84 59 62 42 77 54 24 51 35 41 64 45 51 67 44 52 58 41 39
47 53 38 47 67 57 54 33 43 16 48 54 62 21 32 81 16 22 59 59 39 67 41
22 43 51 48 55 46 61 51 49 64 60 54 48 60 39 75 72 55 30 51 48 55 41
26 60 75 34 23 78 45 63 29 46 48 41 48 72 61 65 61 74 69 52 67 77 57
48 84 58 56 57 27 53 48 69 35 70 8 64 60 66 33 29 84 34 71 80 66 51 43
41 37 55 59 34 59 57 63 59 36 49 38 53 52 52 61 65 25 59 56 27 36 48
57 35 69 55 50 74 71 60 76 35 50 25 79 34 84 67 52 59 46 24 23 64 48
50 54 40 46 23 31 47 43 51 24 60 59 55 44 59 84 58 61 52 74 34 53 83
71 80 46 47 38 61 51 51 28 24 63 74 65 86 85
""")
y = nums("""
87 68 80 51 33 74 63 83 55 84 86 80 83 82 59 80 78 57 67 56 74 85 76
70 58 77 71 71 75 66 89 74 101 95 60 54 74 99 72 84 84 65 72 77 70 83
99 77 66 78 75 66 79 68 67 49 69 78 77 73 51 62 65 84 84 58 70 69 79
95 64 52 88 88 64 75 62 78 82 65 87 63 44 67 48 64 53 52 72 59 53 76
70 61 65 69 61 99 74 72 65 76 59 106 60 64 64 98 74 58 89 52 56 76 74
80 76 92 64 69 75 69 55 72 76 74 57 75 62 60 74 63 57 93 76 69 81 80
73 70 87 62 66 74 43 59 70 38 81 77 91 52 83 57 76 56 56 80 96 90 55
83 65 75 75 80 96 87 79 90 73 69 72 79 79 74 72 38 72 56 81 64 79 73
48 46 38 71 69 70 64 71 82 66 87 55 52 69 48 42
""")
chr_ = lambda s: np.array([int(c) for c in s if c.isdigit()])
diab = chr_("""
1110001010000000010011100100110100100000000000000000001100011000000000
0101100010101110101100010001010000111000011101000000000000001100000010
000011100000010010000000010001100110010100010110011101010100
""")
treat = chr_("""
1100110111010101000001100101000111000110001000001101100111010111000001
1000001010000110110011010000011110011100011111100100011110111100110001
010011111101111011001101010010100000010110111101111010011111
""")
d = pd.DataFrame(dict(age=age, diab=diab, sev=sev, treat=treat, y=y))

# 1. Naive comparison (confounded by indication)
print("Naive:", round(sm.OLS.from_formula("y ~ treat", d).fit().params["treat"], 2))

# 2. Propensity score by logistic regression (no penalty), logit scale for matching
fit = sm.Logit.from_formula("treat ~ age + diab + sev", d).fit(disp=0)
ps = fit.predict(d).to_numpy()
lps = np.log(ps / (1 - ps))
w_ate = np.where(d.treat == 1, 1 / ps, 1 / (1 - ps))      # IPTW (ATE) weights

# 3. Greedy 1:1 nearest-neighbour matching without replacement, caliper 0.2 SD
caliper = 0.2 * lps.std(ddof=1)
treated = np.where(d.treat == 1)[0]
treated = treated[np.argsort(-lps[treated])]               # largest PS first
pool = list(np.where(d.treat == 0)[0]); pairs = []
for i in treated:
    nn = NearestNeighbors(n_neighbors=1).fit(lps[pool].reshape(-1, 1))
    dist, j = nn.kneighbors([[lps[i]]])
    if dist[0, 0] <= caliper:
        pairs.append((i, pool.pop(j[0, 0])))
ti, ci = map(np.array, zip(*pairs))
print("Matched pairs:", len(pairs), "of", len(treated), "treated")

# Standardized mean difference; denominator = pooled SD in the ORIGINAL data
def smd(v, idx=None, w=None):
    x, t = d[v].to_numpy(), d.treat.to_numpy()
    sd0 = np.sqrt((x[t == 1].var(ddof=1) + x[t == 0].var(ddof=1)) / 2)
    if idx is not None: return (x[idx[0]].mean() - x[idx[1]].mean()) / sd0
    w = np.ones(len(x)) if w is None else w
    return (np.average(x[t == 1], weights=w[t == 1]) - np.average(x[t == 0], weights=w[t == 0])) / sd0
bal = pd.DataFrame({v: [smd(v), smd(v, (ti, ci)), smd(v, w=w_ate)]
                    for v in ["age", "diab", "sev"]}, index=["raw", "matched", "iptw"]).T
print(bal.round(3))

# 4. Effect estimates (matched: paired differences; IPTW: WLS with HC1 robust SE)
dif = y[ti] - y[ci]
print(f"Matched: {dif.mean():.2f} (SE {dif.std(ddof=1)/np.sqrt(len(dif)):.2f})")
fw = sm.WLS.from_formula("y ~ treat", d, weights=w_ate).fit(cov_type="HC1")
print(f"IPTW: {fw.params['treat']:.2f} (SE {fw.bse['treat']:.2f})")
```

**What the code does**

- `statsmodels` `Logit` estimates the PS (sklearn's `LogisticRegression` penalizes by default, so we avoid it here).
- The loop is **hand-written greedy matching**: take treated patients from highest PS down, find the nearest unused control on the logit scale with `NearestNeighbors`, accept it only if within the caliper, then remove it from the pool.
- `smd()` handles both the matched and the weighted cases.
- The matched estimate is the mean of the paired outcome differences with its paired-SE; IPTW uses `WLS` with `cov_type="HC1"`.

**Output**

```text
Naive: -4.01
Matched pairs: 77 of 103 treated
        raw  matched   iptw
age   0.324   -0.004  0.066
diab  0.240    0.056 -0.000
sev   0.678    0.068  0.050
Matched: 3.32 (SE 1.37)
IPTW: 2.84 (SE 2.04)
```

**R and Python agree.** Both match greedily on the logit PS (largest score first) with a caliper of 0.2 SD of the logit, so they pick the same 77 pairs and give the same matched SMDs and estimate (3.32). Only the matched SE differs slightly (1.38 vs 1.37): R uses a cluster-robust SE with pairs as clusters, while Python uses the SE of the paired differences. The scale matters. If you leave out `link = "linear.logit"`, `MatchIt` matches on the probability scale instead and picks different pairs (here the age SMD after matching becomes −0.084 and the estimate 3.74). Neither scale is wrong, but they are different analyses: pre-specify the scale and the caliper, and say which you used.

## Interpreting the output

- **Naive difference: −4.01.** The treated group scored 4 points **lower** on average. This is confounding by indication: treated patients were sicker. The truth is +5.
- **Balance before:** SMD for severity 0.68, age 0.32, diabetes 0.24: the groups were not comparable.
- **Balance after matching:** all |SMD| are about 0.1 or below (R: age −0.08, diabetes 0.08, severity 0.10; Python: ≤ 0.07). Balance is adequate, though severity in R sits on the 0.1 threshold, so a more flexible PS model could be tried.
- **Balance after IPTW:** all |SMD| ≤ 0.07. Weighting balances all measured covariates well here.
- **Matched estimate (ATT): 3.32 (SE 1.38 in R, 1.37 in Python).** 95% CI about 0.6 to 6.0. The effect is positive and the CI includes the truth (+5). Only 77 of 103 treated patients were matched, so this estimates the effect in the *matchable* treated population.
- **IPTW estimate (ATE): 2.84 (SE 2.04).** 95% CI about −1.2 to 6.8. It is also in the right direction, but wide: IPTW weights range up to about 6, which inflates variance. The truth is inside the interval.
- **Why do matching and IPTW differ?** They target different populations: matching estimates the **ATT** (here, among the 77 matched treated patients), IPTW the **ATE** (everyone). If the effect varied by patient, the two estimands would differ. In this simulation the true effect is +5 for everyone, so ATT = ATE = +5, and the gap between 3.3 and 2.8 is sampling noise plus method differences (IPTW uses all 200 patients but with unstable weights), not a real difference between the estimands. To estimate the ATT with weighting, give treated patients weight 1 and untreated patients e(x)/(1 − e(x)).
- **Comparison with the truth (+5).** Naive −4.0 (wrong sign), matched 3.3, IPTW 2.8. Both PS methods fix the sign and get much closer. Neither hits +5 exactly: with 200 patients the sampling error is large (SEs 1.4–2.0), and the estimates differ from the truth by less than about 2 SEs.
- **What was not checked.** These results are only valid if no **unmeasured** confounder remains (here, by construction, none exists). In real data you never know.

**How to report it:**
> Of 200 patients, 103 received the therapy. Propensity scores were estimated by logistic regression on age, diabetes, and severity. After 1:1 nearest-neighbour matching on the logit propensity score (caliper 0.2 SD; 77 pairs), all standardized mean differences were below 0.1 (maximum 0.068). The therapy was associated with a 3.3-point higher function score (95% CI 0.6 to 6.0). IPTW gave a similar estimate (2.8 points; 95% CI −1.2 to 6.8). Results assume no unmeasured confounding.

## Check your understanding

<details class="check"><summary>1. "Propensity score adjustment controls for all confounding, measured or not." Is this true?</summary>

No. The propensity score is estimated **only from measured covariates**, so it can balance only those. Unmeasured confounders (for example, frailty that was never recorded) remain unbalanced and bias the result. Randomization is the only method that balances unmeasured factors on average. Use sensitivity analyses, such as the E-value, to see how strong an unmeasured confounder would have to be to explain away your result.
</details>

<details class="check"><summary>2. After matching, the logistic model for the PS has a significant Hosmer–Lemeshow test and a low AUC. Is the PS model bad?</summary>

Not necessarily. The purpose of the PS model is **balance**, not prediction. A high AUC can even signal a problem, because it means treated and untreated patients are easy to tell apart (poor overlap). Judge the model by SMDs (|SMD| < 0.1) and PS overlap, not by fit tests or AUC.
</details>

<details class="check"><summary>3. Some patients with very high severity never receive the untreated option in your data (PS ≈ 0.99). What problem is this and what can you do?</summary>

This is a **positivity violation**. There are no untreated patients comparable to them, so any method must extrapolate. IPTW will give them huge weights; matching will leave them unmatched. Options: restrict the analysis to the region of overlap (and redefine the target population accordingly), trim or truncate extreme weights, or report that no comparison is possible for them.
</details>

<details class="check"><summary>4. In our example the matched analysis used 77 of 103 treated patients. Does that make the result biased?</summary>

Not necessarily biased, but it changes the **target population**: the estimate applies to the treated patients who could be matched within the caliper, not to all 103. The loss of patients also widens the CI (SE 1.38). Always report how many were matched, and compare the characteristics of matched and unmatched treated patients.
</details>

## Common mistakes

- **Treating the PS as a cure for unmeasured confounding.** It does not. State the assumption openly.
- **Including variables affected by the treatment** (mediators or colliders) in the PS model. Use pre-treatment covariates only.
- **Judging balance with p-values.** p-values depend on the sample size, whereas SMDs do not. Use SMDs, and look at the whole distribution for continuous covariates.
- **Ignoring positivity.** Check overlap plots and extreme weights, and trim or restrict when needed.
- **Using ordinary standard errors** after weighting or matching. Use robust, cluster-robust, or bootstrap SEs.
- **Matching with a wide or no caliper**, which keeps bad matches, then claiming "balance".
- **Not checking which scale the software matches on.** In R, `MatchIt` with `distance = "glm"` matches on the probability unless you add `link = "linear.logit"`. The common recommendation (caliper 0.2 SD of the logit) assumes the logit scale.
- **Not defining the estimand (ATT vs ATE)** before choosing matching or weighting.

## Related tests

- [Logistic regression](#/test/logistic-regression): used to estimate the propensity score
- [ANCOVA](#/test/ancova) and [linear regression](#/test/linear-regression): outcome regression adjustment, an alternative or complement
- [Cochran–Mantel–Haenszel test](#/test/cmh): stratified analysis, the ancestor of PS stratification
- [Missing data and multiple imputation](#/test/multiple-imputation): often needed before estimating the PS
- [Cox proportional hazards regression](#/test/cox-ph): for time-to-event outcomes with weights or matched sets

## Key takeaways

- Confounding by indication can reverse the sign of a naive comparison.
- The propensity score is P(treatment | measured covariates). Matching or IPTW uses it to build comparable groups.
- Judge success by covariate balance (|SMD| < 0.1), not p-values or model fit.
- PS methods need **positivity** and **no unmeasured confounding**; neither can be tested from the data.
- Use robust SEs, and say whether the estimand is the ATT (matching) or the ATE (IPTW).

## References

- Rosenbaum PR, Rubin DB. The central role of the propensity score in observational studies for causal effects. *Biometrika.* 1983;70:41–55.
- Austin PC. An introduction to propensity score methods for reducing the effects of confounding in observational studies. *Multivariate Behav Res.* 2011;46:399–424.
- Ho DE, Imai K, King G, Stuart EA. MatchIt: nonparametric preprocessing for parametric causal inference. *J Stat Softw.* 2011;42(8).
- Hernán MA, Robins JM. *Causal Inference: What If.* Chapman & Hall/CRC; 2020.
