# Risk difference, relative risk, odds ratio, and NNT

> Four ways to summarise the difference in a binary outcome between two groups: the absolute gap (RD), the ratio of risks (RR), the ratio of odds (OR), and the number of patients you must treat to prevent one event (NNT).

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Quantify how much a treatment or exposure changes a binary outcome | Binary outcome, two groups (2×2 table) | Independent subjects; fixed follow-up; large enough counts for Wald CIs | Stroke within 1 year, drug vs placebo |

## Learning objectives

- Calculate RD, RR, OR, and NNT from a 2×2 table, with 95% CIs
- Explain why OR ≈ RR only when the outcome is rare
- Choose between relative and absolute measures when talking to patients
- Run the calculations in R and Python and report them correctly

## What it is

A test of two proportions answers "is there a difference?" Effect measures answer the more useful question, "how big is it?" All four measures come from the same 2×2 table, but they tell different stories.

- **Risk difference (RD)** = p₁ − p₀. An *absolute* measure, in percentage points.
- **Relative risk (RR)** = p₁ / p₀. A *ratio* of risks; 1 means no effect.
- **Odds ratio (OR)** = [p₁/(1−p₁)] / [p₀/(1−p₀)]. A ratio of odds; 1 means no effect. It is what logistic regression estimates.
- **Number needed to treat (NNT)** = 1 / |RD|. How many patients must be treated for one extra patient to benefit (or, for harm, the NNH).

RR and OR are given with CIs on the log scale; RD with a Wald CI on its natural scale.

## The intuition

Imagine two news headlines about the same trial: "Drug cuts stroke risk by 40%" and "Drug prevents 8 strokes per 100 patients treated." Both are true. The first is the RR (0.60); the second is the RD (−8 percentage points). The relative number sounds dramatic whatever the baseline risk; the absolute number tells you what a patient actually gains.

Odds are the awkward cousin. A 20% risk is odds of 0.25 (1 to 4). When events are rare, odds and risks are nearly the same number, so OR ≈ RR. When events are common, odds inflate faster than risks, and the OR runs further from 1 than the RR.

In our trial, the RD of −8 points has a standard error of 3.6 points, so the Wald z is −2.20. The RD is about 2.2 standard errors below zero.

<div class="viz" data-viz="z" data-stat="-2.20"></div>

The shaded tails show how often chance alone would give a z this far from 0. They are small, but not tiny, matching a CI for the RD that only just excludes 0 (−15.1 to −0.9 points). Notice what that does to the NNT interval below.

## When to use it

- The outcome is **binary** (event / no event) and there are **two groups**
- You want the **size** of the difference, not just a p-value
- Follow-up time is equal in both groups (otherwise use [survival methods](#/test/kaplan-meier))

**Clinical scenario.** A randomised trial assigns 400 patients with atrial fibrillation to a new anticoagulant or placebo (200 each). Within one year, 24 patients on the drug and 40 on placebo have a stroke. *How much does the drug reduce stroke risk, and how many patients must be treated to prevent one stroke?*

**Don't use it when:**

- You need to adjust for covariates → use [logistic regression](#/test/logistic-regression) (adjusted OR) or Poisson regression with robust errors (adjusted RR)
- Follow-up differs between patients or events are time-dependent → use the [log-rank test](#/test/log-rank) and [Cox regression](#/test/cox-ph) (hazard ratio)
- The same patients are measured twice → use [McNemar's test](#/test/mcnemar)
- Counts are very small (any expected cell < 5) → use [Fisher's exact test](#/test/fisher-exact) and exact CIs

## Assumptions

1. **Independent observations.** One row per patient; no clustering by site or family unless modelled. Check the design.
2. **Same follow-up in both groups.** The measures describe a *fixed* period (here, 1 year).
3. **Adequate counts for Wald CIs.** Roughly ≥ 5–10 events and non-events per group. Otherwise use exact or score-based CIs.
4. **Random assignment for causal reading.** In an observational study, RD, RR, and OR describe association, not effect.

## Hypotheses

These are estimation methods, so the main target is the value of each measure. The matching null values are:

- **H₀:** RD = 0, RR = 1, OR = 1 (all equivalent: the two risks are equal)
- **H₁:** RD ≠ 0, RR ≠ 1, OR ≠ 1

A 95% CI that excludes the null value (0 for RD; 1 for RR and OR) corresponds to p < 0.05.

## Example data

| | Stroke | No stroke | Total | Risk |
|---|---|---|---|---|
| Drug | 24 | 176 | 200 | 12% |
| Placebo | 40 | 160 | 200 | 20% |

## How to run it

### R

```r
# RCT: 24/200 strokes on drug, 40/200 on placebo
a <- 24; n1 <- 200      # drug: events, total
b <- 40; n0 <- 200      # placebo: events, total
p1 <- a / n1; p0 <- b / n0
zc <- qnorm(0.975)

# Risk difference (Wald CI) and its z statistic
rd    <- p1 - p0
se_rd <- sqrt(p1*(1-p1)/n1 + p0*(1-p0)/n0)
rd_ci <- rd + c(-1, 1) * zc * se_rd
cat(sprintf("RD = %.3f (%.3f, %.3f), z = %.2f\n", rd, rd_ci[1], rd_ci[2], rd/se_rd))

# Relative risk: CI on the log scale
rr    <- p1 / p0
rr_ci <- exp(log(rr) + c(-1, 1) * zc * sqrt(1/a - 1/n1 + 1/b - 1/n0))
cat(sprintf("RR = %.3f (%.3f, %.3f)\n", rr, rr_ci[1], rr_ci[2]))

# Odds ratio: CI on the log scale
or    <- (a/(n1-a)) / (b/(n0-b))
or_ci <- exp(log(or) + c(-1, 1) * zc * sqrt(1/a + 1/(n1-a) + 1/b + 1/(n0-b)))
cat(sprintf("OR = %.3f (%.3f, %.3f)\n", or, or_ci[1], or_ci[2]))

# NNT = 1/|RD|; CI limits come from inverting the RD limits (only valid if CI excludes 0)
cat(sprintf("NNT = %.1f (%.1f to %.1f)\n", 1/abs(rd), 1/abs(rd_ci[1]), 1/abs(rd_ci[2])))

# A common outcome: response 140/200 vs 80/200
cat(sprintf("Common outcome: RR = %.2f, OR = %.2f\n",
            (140/200)/(80/200), (140/60)/(80/120)))

# Cross-check with DescTools (rows = groups, columns = event / no event)
library(DescTools)
tab <- matrix(c(24, 176, 40, 160), nrow = 2, byrow = TRUE)
OddsRatio(tab, conf.level = 0.95, method = "wald")
```

**What the code does**

- Computes each measure from the four cell counts, using the standard log-scale standard errors for RR and OR.
- `1/abs(rd_ci[1])` and `1/abs(rd_ci[2])` turn the RD limits into NNT limits. The RD limit nearest 0 (−0.009) gives the *largest* NNT.
- The last block shows what happens with a common outcome: the OR (3.50) is far from the RR (1.75).
- `DescTools::OddsRatio` (and `RelRisk`) reproduce the OR as a cross-check. `RelRisk` uses a slightly different SE formula, so its CI can differ in the third decimal.

**Output**

```text
RD = -0.080 (-0.151, -0.009), z = -2.20
RR = 0.600 (0.376, 0.957)
OR = 0.545 (0.315, 0.945)
NNT = 12.5 (6.6 to 116.6)
Common outcome: RR = 1.75, OR = 3.50
odds ratio     lwr.ci     upr.ci 
 0.5454545  0.3148623  0.9449232 
```

### Python

```python
import numpy as np
from scipy.stats import norm
from statsmodels.stats.contingency_tables import Table2x2

a, n1, b, n0 = 24, 200, 40, 200          # events and totals: drug, placebo
p1, p0 = a / n1, b / n0
zc = norm.ppf(0.975)

# Risk difference (Wald CI) and z
rd = p1 - p0
se = np.sqrt(p1*(1-p1)/n1 + p0*(1-p0)/n0)
lo, hi = rd - zc*se, rd + zc*se
print(f"RD = {rd:.3f} ({lo:.3f}, {hi:.3f}), z = {rd/se:.2f}")

# RR and OR with CIs from statsmodels (rows = groups, columns = event / no event)
t = Table2x2(np.array([[24, 176], [40, 160]]))
rr_lo, rr_hi = t.riskratio_confint()
or_lo, or_hi = t.oddsratio_confint()
print(f"RR = {t.riskratio:.3f} ({rr_lo:.3f}, {rr_hi:.3f})")
print(f"OR = {t.oddsratio:.3f} ({or_lo:.3f}, {or_hi:.3f})")

# NNT from the RD limits (valid only because the RD CI excludes 0)
print(f"NNT = {1/abs(rd):.1f} ({1/abs(lo):.1f} to {1/abs(hi):.1f})")

# A common outcome: response 140/200 vs 80/200
c = Table2x2(np.array([[140, 60], [80, 120]]))
print(f"Common outcome: RR = {c.riskratio:.2f}, OR = {c.oddsratio:.2f}")
```

**What the code does**

- `Table2x2` expects rows = exposure groups and columns = event / no event, with the *exposed* group in the first row and the *event* in the first column.
- `riskratio_confint()` and `oddsratio_confint()` give Wald CIs on the log scale, identical to the R formulas.
- The RD has no built-in in `Table2x2`, so it is computed by hand.

**Output**

```text
RD = -0.080 (-0.151, -0.009), z = -2.20
RR = 0.600 (0.376, 0.957)
OR = 0.545 (0.315, 0.945)
NNT = 12.5 (6.6 to 116.6)
Common outcome: RR = 1.75, OR = 3.50
```

## Interpreting the output

- **RD = −8.0 points (95% CI −15.1 to −0.9).** The drug lowers one-year stroke risk from 20% to 12%, an absolute reduction of 8 points. The CI excludes 0.
- **RR = 0.60 (0.38 to 0.96).** The drug's risk is 60% of the placebo risk: a 40% *relative* reduction.
- **OR = 0.55 (0.31 to 0.95).** Further from 1 than the RR, because a 12–20% outcome is not that rare. Reading it as "45% fewer strokes" would overstate the benefit.
- **NNT = 12.5, so 13.** Treat about 13 patients for a year to prevent one stroke. The CI (7 to 117) is very wide and asymmetric: the NNT is the inverse of the RD, and the RD's upper limit is close to 0.
- **Common outcome check.** When 70% respond on drug and 40% on control, RR = 1.75 but OR = 3.50. They diverge, so never read an OR as a "times as likely" ratio unless the outcome is rare (roughly < 10%).

**How to report it:**
> Within one year, stroke occurred in 12.0% (24/200) of patients on the drug and 20.0% (40/200) on placebo (risk difference −8.0 percentage points, 95% CI −15.1 to −0.9; RR 0.60, 95% CI 0.38 to 0.96; OR 0.55, 95% CI 0.31 to 0.95). The NNT was 13 (95% CI 7 to 117).

**Communicating with patients.** Say "Out of 100 people like you, about 20 would have a stroke this year without the drug and about 12 with it." Absolute numbers prevent the common error of hearing "40% lower risk" as "my risk falls from 50% to 10%." The same RR of 0.60 is a big deal at a baseline risk of 20% (RD −8) and trivial at a baseline of 0.1% (RD −0.04 points, NNT 2,500). Always give the baseline risk.

## Check your understanding

<details class="check"><summary>1. A cohort study finds that an exposure doubles the odds of a rare cancer (OR = 2.0, outcome prevalence 0.5%). Is it fair to say the exposure doubles the risk?</summary>

Roughly yes. With an outcome this rare, the OR is almost identical to the RR, so "about twice as likely" is a fair summary. For a common outcome (say 40% prevalence) the same OR would correspond to a much smaller RR, and that wording would be wrong.
</details>

<details class="check"><summary>2. The RR is 0.60 with CI 0.38 to 0.96. A colleague says the drug is "highly effective." What else do you need to know?</summary>

The baseline (control) risk, and the absolute difference. A 40% relative reduction from 20% saves 8 per 100; from 0.1% it saves 0.04 per 100. The wide CI (anything from a 62% to a 4% reduction) also matters.
</details>

<details class="check"><summary>3. The RD CI had run from −15 to +2 points. Could you still report an NNT with a CI?</summary>

Not as a single interval. The RD CI includes 0, which means an NNT of infinity, so the CI for the NNT runs from 7 up through ∞ and then wraps around to a number needed to harm of about 50. Report it as "NNTB 7 to ∞ to NNTH 50", or simply report the RD with its CI.
</details>

<details class="check"><summary>4. Using our numbers, what is the NNT if the placebo risk were 10% and the RR stayed 0.60?</summary>

The drug risk would be 6%, so RD = −4 points and NNT = 1/0.04 = 25. The NNT doubled even though the RR did not change, because the NNT depends on the baseline risk.
</details>

## Common mistakes

- **Reading an OR as an RR for common outcomes.** It exaggerates the effect when the outcome is common.
- **Reporting only a relative measure.** Always add the absolute risks or the RD.
- **Rounding the NNT down.** Round *up* (12.5 becomes 13).
- **Giving a CI for the NNT when the RD CI crosses 0.** The interval then is not a simple range.
- **Flipping the reference group.** Say which group is the reference (placebo here); swapping it inverts the RR and OR and changes the sign of the RD.
- **Applying NNT across different follow-up periods.** An NNT of 13 is for 1 year, not forever.

## Related tests

- [Two-proportion z-test](#/test/two-proportion-z): the significance test for the RD
- [Fisher's exact test](#/test/fisher-exact): small counts
- [Chi-square test of independence](#/test/chi-square-independence): is there any association?
- [Logistic regression](#/test/logistic-regression): adjusted odds ratios
- [Cochran–Mantel–Haenszel test](#/test/cmh): a common OR across strata
- [Confidence intervals](#/test/confidence-intervals): how to read the intervals above

## Key takeaways

- RD = p₁ − p₀ (absolute), RR = p₁/p₀, OR = odds₁/odds₀, NNT = 1/|RD| (round up).
- OR ≈ RR only when the outcome is rare (under about 10%); for common outcomes the OR overstates the effect.
- Relative measures hide the baseline risk. Communicate with absolute risks and the NNT.
- Compute CIs for RR and OR on the log scale; give an NNT CI only when the RD CI excludes 0.
- Always state the reference group and the follow-up period.

## References

- Altman DG. Confidence intervals for the number needed to treat. *BMJ.* 1998;317:1309–12.
- Agresti A. *Categorical Data Analysis.* 3rd ed. Wiley; 2013. Ch. 2.
- DescTools documentation: `RelRisk`, `OddsRatio`, `BinomDiffCI`
- statsmodels documentation: `Table2x2`
