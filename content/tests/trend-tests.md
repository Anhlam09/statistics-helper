# Tests for trend (Cochran–Armitage, Jonckheere–Terpstra)

> Test whether an outcome rises or falls steadily across ordered groups, such as increasing drug doses: Cochran–Armitage for binary outcomes and Jonckheere–Terpstra for continuous or ordinal outcomes.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Groups have a natural order and you expect a monotonic trend | Binary (Cochran–Armitage); continuous/ordinal (Jonckheere–Terpstra) | Independent subjects; ordered groups defined before the analysis | Nausea rate or LDL-C change across 0, 10, 20, 40 mg |

## Learning objectives

- Recognize when an ordered-group design calls for a trend test rather than an omnibus test
- Explain why a trend test is more powerful against a monotonic alternative
- Run Cochran–Armitage and Jonckheere–Terpstra in R and in Python (by formula)
- Interpret and report a dose–response trend with an effect size

## What it is

An omnibus test (chi-square, Kruskal–Wallis) asks "do the groups differ in *any* way?" A trend test asks the narrower, more informative question: "does the outcome go steadily up (or down) as the group order increases?"

- **Cochran–Armitage:** a weighted linear-regression-style test on proportions. With scores sⱼ (usually the dose), T = Σ sⱼ(rⱼ − nⱼp̄), where rⱼ is the number of events and p̄ the pooled proportion. Divide by its SD to get z; z² follows χ² with **1 df**.
- **Jonckheere–Terpstra:** counts, over all pairs of groups i < j, how many pairs of observations are ordered in the expected direction (the later group's value is larger). JT is the total count; it is standardised to a z score using its mean and variance under H₀.
- **Effect size:** for Cochran–Armitage, the change in log-odds per unit dose (logistic slope) or the absolute rate in each group; for Jonckheere–Terpstra, group medians or the probability that a higher-dose patient has a larger value.

## The intuition

Think of a dose–response curve. If a drug works, more drug should mean more effect. An omnibus test spreads its "attention" over every possible way four groups could differ (group 2 high, group 3 low, and so on). A trend test puts all its attention on one pattern, a steady slope. When that pattern is really there, the signal is concentrated and the test is more powerful. In the omnibus χ² you "spend" 3 degrees of freedom; in the trend test only 1, so the same signal faces less noise.

Here are the LDL-C reductions at four ascending doses. The means climb steadily, with some overlap between neighbours.

<div class="viz" data-viz="dots" data-unit="mg/dL LDL-C reduction" data-groups='{"0 mg":[2,-3,5,1,-4,6,0,3],"10 mg":[8,4,12,1,9,15,6,3],"20 mg":[14,10,6,19,12,8,22,11],"40 mg":[20,26,15,31,18,24,12,29]}' data-caption="Each dose group sits higher than the previous one. A trend test uses this ordering."></div>

For the nausea data, the Cochran–Armitage statistic is χ² = 15.35 on 1 df. The omnibus χ² on the same table is 15.43, but on 3 df, so its p-value (0.0015) is much weaker than the trend p-value (0.00009).

<div class="viz" data-viz="chisq" data-stat="15.35" data-df="1"></div>

The observed value is far into the upper tail of the χ²(1) curve, so almost no area is left beyond it.

## When to use it

- The groups have a **natural order** (dose, severity stage, time since exposure) fixed **before** looking at data
- You expect (or want to check for) a **monotonic** trend
- Subjects are independent, one group each

**Clinical scenario.** A phase 2 trial randomises 200 patients to placebo (0 mg) or 10, 20, or 40 mg of a new drug (50 per arm). Nausea occurs in 5, 9, 14, and 21 patients respectively. *Does the nausea rate rise with dose?* In a second, smaller study, 8 patients per dose have LDL-C reduction measured. *Does cholesterol lowering increase with dose?*

**Don't use it when:**

- The groups are **unordered** (drug A, B, C) → use the [chi-square test](#/test/chi-square-independence) or [Kruskal–Wallis](#/test/kruskal-wallis)
- The expected pattern is **non-monotonic** (e.g., a U-shape or an inverted U) → use an omnibus test
- The same patients are measured at each dose (crossover) → use [Cochran's Q](#/test/cochrans-q) or the [Friedman test](#/test/friedman)
- You need to model the dose–response curve or adjust for covariates → use [logistic](#/test/logistic-regression) or [linear regression](#/test/linear-regression) with dose as a numeric predictor

## Assumptions

1. **Independent observations.** Each subject is in one dose group. Check the design.
2. **Meaningful order.** The ordering must be specified in advance. Choosing it after seeing the data invalidates the p-value.
3. **Monotonic alternative.** The test is powerful for a steady rise or fall, and weak for a non-monotonic pattern. Plot the group rates.
4. **Scores (Cochran–Armitage).** Results depend on the numeric scores. Equally spaced scores (1, 2, 3, 4) are the default; real doses (0, 10, 20, 40) are better when dose spacing is meaningful. State which one you used.
5. **Reasonably large counts** (Cochran–Armitage uses a normal approximation). Check expected counts or use an exact version.

## Hypotheses

- **H₀:** the outcome distribution is the same in all groups (no trend)
- **H₁:** the outcome increases (or decreases) monotonically with group order. A two-sided test detects either direction; a one-sided test is justified only if the direction is pre-specified.

## Example data

| Dose | 0 mg | 10 mg | 20 mg | 40 mg |
|---|---|---|---|---|
| Patients | 50 | 50 | 50 | 50 |
| Nausea | 5 (10%) | 9 (18%) | 14 (28%) | 21 (42%) |

LDL-C reduction (mg/dL), 8 patients per dose: see the dot plot above (medians 1.5, 7, 11.5, 22).

## How to run it

### R

```r
# Nausea by dose: 50 patients per arm
dose <- c(0, 10, 20, 40)          # mg
n    <- c(50, 50, 50, 50)
ae   <- c(5, 9, 14, 21)           # patients with nausea

# Cochran-Armitage test for trend (scores = actual doses)
prop.trend.test(ae, n, score = dose)

# Omnibus chi-square ignores the ordering
chisq.test(rbind(ae, n - ae))
```

```r
library(clinfun)
# LDL-C reduction (mg/dL), 8 patients per dose, doses in ascending order
y <- c( 2, -3,  5,  1, -4,  6,  0,  3,     # 0 mg (placebo)
        8,  4, 12,  1,  9, 15,  6,  3,     # 10 mg
       14, 10,  6, 19, 12,  8, 22, 11,     # 20 mg
       20, 26, 15, 31, 18, 24, 12, 29)     # 40 mg
dose <- rep(c(0, 10, 20, 40), each = 8)

# Jonckheere-Terpstra test (normal approximation because of ties)
jonckheere.test(y, dose, alternative = "two.sided")

# Kruskal-Wallis ignores the ordering
kruskal.test(y ~ factor(dose))
```

**What the code does**

- `prop.trend.test(events, totals, score = dose)` is base R's Cochran–Armitage test; `score` sets the numeric spacing of the groups.
- `chisq.test` on the 2×4 table is the omnibus comparison, for contrast.
- `clinfun::jonckheere.test(x, g)` takes the outcome and the **numeric group order**. `DescTools::JonckheereTerpstraTest` is an alternative.
- Because the data contain tied values, the p-value is based on the normal approximation, which is why R prints a warning. That warning is informational here.

**Output**

```text
	Chi-squared Test for Trend in Proportions

data:  ae out of n ,
 using scores: 0 10 20 40
X-squared = 15.354, df = 1, p-value = 8.912e-05

	Pearson's Chi-squared test

data:  rbind(ae, n - ae)
X-squared = 15.435, df = 3, p-value = 0.001481

	Jonckheere-Terpstra test

JT = 349, p-value = 1.317e-07
alternative hypothesis: two.sided

	Kruskal-Wallis rank sum test

data:  y by factor(dose)
Kruskal-Wallis chi-squared = 22.751, df = 3, p-value = 4.551e-05
```

### Python

```python
import numpy as np
from scipy import stats

# ---- Cochran-Armitage by formula ----
dose = np.array([0, 10, 20, 40.])
n    = np.array([50., 50., 50., 50.])
ae   = np.array([5., 9., 14., 21.])
N, pbar = n.sum(), ae.sum() / n.sum()
T = np.sum(dose * (ae - n * pbar))                       # score-weighted excess events
var = pbar * (1 - pbar) * (np.sum(n * dose**2) - np.sum(n * dose)**2 / N)
z = T / np.sqrt(var)
print(f"Cochran-Armitage: X2 = {z**2:.3f}, df = 1, p = {2*stats.norm.sf(abs(z)):.3e}")
chi2, p, df, _ = stats.chi2_contingency(np.array([ae, n - ae]))
print(f"Chi-square (omnibus): X2 = {chi2:.3f}, df = {df}, p = {p:.4f}")

# ---- Jonckheere-Terpstra by formula ----
groups = [[2,-3,5,1,-4,6,0,3], [8,4,12,1,9,15,6,3],
          [14,10,6,19,12,8,22,11], [20,26,15,31,18,24,12,29]]
JT = 0.0
for i in range(len(groups)):                  # count pairs ordered like the doses
    for j in range(i + 1, len(groups)):
        for x in groups[i]:
            for y in groups[j]:
                JT += (y > x) + 0.5 * (y == x)
ns = np.array([len(g) for g in groups]); Nn = ns.sum()
E = (Nn**2 - np.sum(ns**2)) / 4
V = (Nn*(Nn-1)*(2*Nn+5) - np.sum(ns*(ns-1)*(2*ns+5))) / 72   # no tie correction, as in clinfun
zj = (JT - E) / np.sqrt(V)
print(f"Jonckheere-Terpstra: JT = {JT:.0f}, E = {E:.0f}, z = {zj:.3f}, p = {2*stats.norm.sf(abs(zj)):.3e}")
kw = stats.kruskal(*groups)
print(f"Kruskal-Wallis: H = {kw.statistic:.3f}, df = 3, p = {kw.pvalue:.3e}")
```

**What the code does**

- SciPy has no built-in trend tests, so both are written from the formulas. The Cochran–Armitage block is the same arithmetic that `prop.trend.test` performs.
- `T` measures how strongly events are concentrated in the high-score groups; `var` is its variance under H₀.
- The loop adds 1 for every pair (lower-dose value, higher-dose value) that is in the expected order and 0.5 for ties. `E` and `V` are the mean and variance of JT under H₀.
- The variance here has no tie correction, matching `clinfun`. A tie-corrected variance would give a slightly different p-value (1.27 × 10⁻⁷), the same conclusion.

**Output**

```text
Cochran-Armitage: X2 = 15.354, df = 1, p = 8.912e-05
Chi-square (omnibus): X2 = 15.435, df = 3, p = 0.0015
Jonckheere-Terpstra: JT = 349, E = 192, z = 5.277, p = 1.317e-07
Kruskal-Wallis: H = 22.751, df = 3, p = 4.551e-05
```

## Interpreting the output

- **Cochran–Armitage: χ²(1) = 15.35, p = 0.00009.** Nausea increases steadily with dose (10%, 18%, 28%, 42%).
- **Omnibus χ²(3) = 15.43, p = 0.0015.** Same signal, but spread over 3 df. The test is significant but the p-value is roughly 17 times larger than the trend p-value.
- **Jonckheere–Terpstra: JT = 349, z = 5.28, p = 1.3 × 10⁻⁷.** Under H₀ we would expect JT = 192; we see 349, so LDL-C reduction rises with dose.
- **Kruskal–Wallis: H(3) = 22.75, p = 0.00005.** Also significant, but about 350 times larger p than the trend test.
- **Neither test gives the size of the effect.** Add group rates or medians (nausea 10% → 42%; median reduction 1.5 → 22 mg/dL), and consider a regression slope (e.g., OR per 10 mg).

**How to report it:**
> Nausea incidence rose with dose, from 10% (5/50) at 0 mg to 42% (21/50) at 40 mg (Cochran–Armitage trend test with dose scores, χ²(1) = 15.35, p < 0.001). LDL-C reduction also increased with dose (median 1.5, 7, 11.5, and 22 mg/dL; Jonckheere–Terpstra z = 5.28, p < 0.001).

⚠️ A significant trend says the relationship is monotonic, not that it is linear or that every adjacent pair differs. It says nothing about which dose to choose.

## Check your understanding

<details class="check"><summary>1. Patients with stage I, II, III, and IV cancer are tested for a biomarker (positive or negative). Which test fits the question "does positivity rise with stage"?</summary>

Cochran–Armitage. The outcome is binary and the groups (stage) are ordered. A plain chi-square would ignore the ordering and be less powerful.
</details>

<details class="check"><summary>2. The Cochran–Armitage test gives p = 0.30 across 4 doses. Does that show there is no dose–response?</summary>

No. It shows no *evidence of a steady linear-in-score trend*. The sample may be too small, or the pattern may be non-monotonic (e.g., a plateau or reversal at the highest dose), which a trend test is poor at detecting. Look at the group rates and the CI of the slope.
</details>

<details class="check"><summary>3. Suppose the nausea rates were 10%, 40%, 12%, 11% for the four doses. Which test is better, trend or chi-square, and why?</summary>

The omnibus chi-square. The pattern is not monotonic, so the trend test would see almost no steady slope and have low power, while the chi-square detects any difference in any group.
</details>

<details class="check"><summary>4. Using our numbers: the trend test has 1 df and the omnibus χ² has 3 df, with almost the same statistic (15.35 vs 15.43). Why is the trend p-value smaller?</summary>

For the same statistic value, a χ² with fewer degrees of freedom has a thinner upper tail, so p is smaller. The trend test puts essentially all the 15.4 units of "signal" into one df because the data really are monotonic. The omnibus test dilutes it across 3.
</details>

## Common mistakes

- **Choosing the order or scores after looking at the data.** Pre-specify them.
- **Treating a significant trend as proof of linearity** or as a way to pick a dose.
- **Using a trend test on unordered groups** (e.g., different drugs).
- **Applying a one-sided test after seeing the direction.** Use two-sided unless the direction was pre-specified.
- **Forgetting ties.** Jonckheere–Terpstra with many ties needs a tie-corrected variance or a permutation p-value.
- **Reporting only p.** Give group rates or medians and, if possible, a slope with its CI.

## Related tests

- [Chi-square test of independence](#/test/chi-square-independence): omnibus test for any difference in proportions
- [Kruskal–Wallis test](#/test/kruskal-wallis): omnibus rank test for unordered groups
- [Logistic regression](#/test/logistic-regression): model the dose slope on the log-odds
- [Spearman rank correlation](#/test/spearman) and [Kendall's tau](#/test/kendall-tau): ordinal associations
- [Ordinal logistic regression](#/test/ordinal-regression): ordered outcome categories

## Key takeaways

- Use a trend test when groups have a **pre-specified order** and you expect a **monotonic** pattern.
- Cochran–Armitage (binary outcome) has 1 df; Jonckheere–Terpstra (ordinal/continuous) uses a z statistic from the count of correctly ordered pairs.
- Against a monotonic alternative, trend tests beat chi-square and Kruskal–Wallis; against a non-monotonic one, they lose.
- A significant trend is not evidence of linearity; report group rates or medians and a slope.
- State the scores and sidedness you used.

## References

- Armitage P. Tests for linear trends in proportions and frequencies. *Biometrics.* 1955;11:375–86.
- Jonckheere AR. A distribution-free k-sample test against ordered alternatives. *Biometrika.* 1954;41:133–45.
- R documentation: `?prop.trend.test`; `clinfun::jonckheere.test`
- Agresti A. *Categorical Data Analysis.* 3rd ed. Wiley; 2013. Ch. 5.
