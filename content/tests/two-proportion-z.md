# Two-proportion z-test

> Compares the proportion of patients with an event (such as a response) between two independent groups, and reports the difference in proportions with a confidence interval.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Compare a binary outcome between two independent groups | Binary (yes/no) in each group | Independent patients; expected counts of at least 5 in every cell | Asthma response: 60/100 on drug vs 42/100 on placebo |

## Learning objectives

- Recognize a two-independent-group, binary-outcome question
- Explain the z statistic as a difference in rates divided by its chance wobble, and its link to the 2 × 2 χ²
- Run the test in R and Python and report the risk difference with its CI, plus RR and OR
- Choose between pooled and unpooled variance, and know when to switch to an exact test

## What it is

You have a success count and a sample size in each of two independent groups. The test asks whether the two true proportions are equal. Under H₀ both groups share one common proportion, which is estimated by pooling all patients. The z statistic divides the observed difference by its standard error under that pooled assumption. Squaring z gives exactly the Pearson χ² statistic of the 2 × 2 table, so the two tests are the same test.

- **Test statistic:** z = (p̂₁ − p̂₂) / √[ p̂(1 − p̂)(1/n₁ + 1/n₂) ], with p̂ the pooled proportion
- **Effect size (primary):** the **risk difference** RD = p̂₁ − p̂₂ with a 95% CI (unpooled Wald standard error)
- **Secondary effect sizes:** the **risk ratio** RR = p̂₁ / p̂₂ and the **odds ratio** OR = [p̂₁/(1 − p̂₁)] / [p̂₂/(1 − p̂₂)], each with a CI on the log scale

## The intuition

The z statistic is a **signal-to-noise ratio** again. The signal is the difference in response rates (60% − 42% = 18 percentage points). The noise is the standard error of that difference: with 100 patients per arm, a response rate is only known to within roughly ±5 points, so a gap of 18 points is about 2.5 standard errors wide.

Analogy: two coins. Flip each 100 times and compare the number of heads. A gap of 2 or 3 heads is easily luck; a gap of 18 is not. The test asks how often two *identical* coins (here, one common response rate of 51%, the pooled rate) would show a gap this large.

<div class="viz" data-viz="z" data-stat="2.55"></div>

**What to notice:** z = 2.55 (z² = 6.48, the same number as the χ² statistic) lies in the shaded tails, so under H₀ a gap this big in either direction happens about 1% of the time (p = 0.011). The test is two-sided: a gap of 18 points in favour of placebo would be just as surprising.

The table below shows what the data look like; the interval from 4.4 to 31.6 points is the range of true differences compatible with them.

| Arm | Responders | Non-responders | Response rate |
|---|---|---|---|
| Drug | 60 | 40 | 60% |
| Placebo | 42 | 58 | 42% |

Why pooling? Under H₀ both arms share one true rate, so the best estimate of the noise uses all 200 patients. The CI for the difference makes no such assumption, so it uses each arm's own rate. That is why a p-value and a CI can occasionally disagree near the boundary.

## When to use it

Use the two-proportion z-test when **all** of these are true:

- The outcome is **binary** (responder, event, infection)
- There are **two independent groups** (treatment vs placebo, exposed vs unexposed)
- The groups are large enough for the normal approximation: roughly **expected counts of at least 5** in every cell of the 2 × 2 table

**Clinical scenario.** A randomised trial assigns 200 patients with moderate asthma to a new inhaled drug or placebo (100 per arm). The outcome is "clinical response at week 12" (an improvement of at least 0.5 points on the Asthma Control Questionnaire). Sixty drug patients and 42 placebo patients respond. *Is the response rate different, and by how many percentage points?*

**Don't use it when:**

- Expected counts are small (below 5) → use [Fisher's exact test](#/test/fisher-exact)
- The same patients appear in both groups (before/after) → use [McNemar's test](#/test/mcnemar)
- Data come from several sites or strata and you want a stratified comparison → use the [Cochran–Mantel–Haenszel test](#/test/cmh)
- You need to adjust for baseline covariates → use [logistic regression](#/test/logistic-regression)
- The aim is to show non-inferiority or equivalence rather than a difference → see [equivalence and non-inferiority tests](#/test/tost-equivalence)
- There are more than two groups or outcome categories → use the [chi-square test of independence](#/test/chi-square-independence)

## Assumptions

1. **The two groups are independent, and patients within a group are independent.** Check from the design: no patient in both arms, no clustering (for example, families or sites) ignored.
2. **Each patient has a binary outcome measured the same way.** Check that the response definition is identical in both arms and fixed in advance.
3. **The sample is large enough for the normal approximation.** Compute expected counts n·p̂ and n·(1 − p̂) in each arm; all should be at least 5 (some authors use 10). Otherwise use Fisher's exact test.
4. **Complete follow-up.** Patients lost to follow-up are either excluded or counted as non-responders by a pre-specified rule. Report how many.

## Hypotheses

- **H₀:** p₁ = p₂ (equal response probabilities; risk difference = 0, RR = 1, OR = 1)
- **H₁:** p₁ ≠ p₂ (two-sided)

## Example data

| Arm | Responders | Non-responders | n | Response rate |
|---|---|---|---|---|
| Drug | 60 | 40 | 100 | 60% |
| Placebo | 42 | 58 | 100 | 42% |

## How to run it

### R

```r
# Responders and group sizes: drug first, placebo second
resp <- c(60, 42)
n    <- c(100, 100)

# 1. z test (shown as chi-square), NO continuity correction
prop.test(resp, n, correct = FALSE)

# 2. Same test with R's default Yates correction (p-value only)
prop.test(resp, n, correct = TRUE)$p.value

# 3. Risk ratio and odds ratio with Wald CIs on the log scale
rr <- (60/100) / (42/100)
or <- (60 * 58) / (40 * 42)
se_lrr <- sqrt(1/60 - 1/100 + 1/42 - 1/100)
se_lor <- sqrt(1/60 + 1/40 + 1/42 + 1/58)
cat("RR =", round(rr, 3), " 95% CI", round(exp(log(rr) + c(-1, 1) * qnorm(.975) * se_lrr), 3), "\n")
cat("OR =", round(or, 3), " 95% CI", round(exp(log(or) + c(-1, 1) * qnorm(.975) * se_lor), 3), "\n")
```

**What the code does**

- `prop.test(x, n)` takes the vector of successes and the vector of group sizes. The first group is the first element, so the difference is **drug − placebo**.
- **`correct = FALSE` matters.** For a 2 × 2 comparison R applies the **Yates continuity correction by default**, which makes p-values larger and the CI wider. Without it, R's X² equals z², and it matches the standard z test and statsmodels.
- The confidence interval printed by `prop.test` is the risk difference interval (unpooled Wald; with `correct = TRUE` it is widened by the correction).
- R reports X-squared rather than z; take the square root for z.
- RR and OR are computed by hand because `prop.test` does not give them. The standard errors are the usual large-sample log-scale formulas.

**Output**

```text
	2-sample test for equality of proportions without continuity correction

data:  resp out of n
X-squared = 6.4826, df = 1, p-value = 0.01089
alternative hypothesis: two.sided
95 percent confidence interval:
 0.04370145 0.31629855
sample estimates:
prop 1 prop 2 
  0.60   0.42 

[1] 0.01618823
RR = 1.429  95% CI 1.079 1.891 
OR = 2.071  95% CI 1.179 3.64 
```

### Python

```python
import numpy as np
from statsmodels.stats.proportion import proportions_ztest, confint_proportions_2indep

resp = np.array([60, 42])    # drug, placebo
n    = np.array([100, 100])

# 1. z test: pooled variance under H0, no continuity correction
z, p = proportions_ztest(resp, n, value=0)
print(f"z = {z:.4f}, z^2 = {z**2:.4f}, p = {p:.5f}")

# 2. Risk difference with Wald CI (unpooled SE)
rd = resp[0]/n[0] - resp[1]/n[1]
lo, hi = confint_proportions_2indep(resp[0], n[0], resp[1], n[1],
                                    method="wald", compare="diff")
print(f"Risk difference = {rd:.3f}, 95% CI ({lo:.3f}, {hi:.3f})")

# 3. Risk ratio and odds ratio with log-scale Wald CIs
rr_lo, rr_hi = confint_proportions_2indep(resp[0], n[0], resp[1], n[1],
                                          method="log", compare="ratio")
or_lo, or_hi = confint_proportions_2indep(resp[0], n[0], resp[1], n[1],
                                          method="logit", compare="odds-ratio")
print(f"RR = {(60/100)/(42/100):.3f}, 95% CI ({rr_lo:.3f}, {rr_hi:.3f})")
print(f"OR = {(60*58)/(40*42):.3f}, 95% CI ({or_lo:.3f}, {or_hi:.3f})")
```

**What the code does**

- `proportions_ztest(count, nobs, value=0)` tests the difference against 0 using the **pooled** standard error. It has **no continuity correction**, so it matches R's `prop.test(..., correct = FALSE)` and not R's default.
- `confint_proportions_2indep(..., method="wald", compare="diff")` returns the CI for the risk difference using the unpooled SE, the same interval R prints.
- The `compare` argument chooses the scale (`"diff"`, `"ratio"`, `"odds-ratio"`), and the allowed `method` depends on it: `"log"` for ratios and `"logit"` for odds ratios are the Wald-type options.
- SciPy's `chi2_contingency` can also run this test but applies Yates for 2 × 2 tables unless `correction=False` is passed.

**Output**

```text
z = 2.5461, z^2 = 6.4826, p = 0.01089
Risk difference = 0.180, 95% CI (0.044, 0.316)
RR = 1.429, 95% CI (1.079, 1.891)
OR = 2.071, 95% CI (1.179, 3.640)
```

### SAS

```sas
data trial;
  input arm $ response $ count;
  datalines;
Drug    Yes 60
Drug    No  40
Placebo Yes 42
Placebo No  58
;
run;

proc freq data=trial order=data;
  weight count;
  tables arm*response / chisq riskdiff(cl=wald) relrisk;
run;
```

**What the code does**

- `order=data` keeps categories in the order they appear: Drug before Placebo, and **Yes before No**. This matters because SAS computes risks for **column 1** (here "Yes") and the difference as row 1 − row 2.
- `chisq` gives the Pearson χ² (equal to z² with no continuity correction). For 2 × 2 tables SAS prints it separately from the continuity-adjusted χ², so there is no hidden Yates default.
- `riskdiff(cl=wald)` prints the risk difference with the Wald CI. Other CL options (Newcombe, Miettinen–Nurminen) are available.
- `relrisk` prints the odds ratio and relative risk with 95% CIs. Check which column's risk it is reporting (the "Cohort (Col1 Risk)" line is the response risk ratio).

## Interpreting the output

- **Response rates: 60% vs 42%.** The risk difference is **+18.0 percentage points**.
- **95% CI for the difference: 4.4 to 31.6 percentage points.** The data are compatible with a true advantage anywhere from a small gain (about 4 points) to a large one (about 32 points). The interval excludes 0.
- **z = 2.55 (z² = 6.48), p = 0.011.** If the two arms truly had the same response rate, a difference this large would occur in about 1% of trials.
- **Risk ratio 1.43 (95% CI 1.08 to 1.89).** Responders are about 43% more likely on the drug.
- **Odds ratio 2.07 (95% CI 1.18 to 3.64).** The odds of response are roughly doubled. Note that the OR is further from 1 than the RR, because response is common (it is not a rare outcome).
- **Yates correction.** With R's default `correct = TRUE`, p = 0.016. The conclusion does not change here, but the p-value does, which is why the option should be stated.
- **R and Python agree** on z² = X² = 6.4826, p = 0.01089, and on the interval for the risk difference.

**How to report it:**
> Week-12 response occurred in 60/100 patients on drug and 42/100 on placebo (risk difference 18.0 percentage points, 95% CI 4.4 to 31.6; risk ratio 1.43, 95% CI 1.08 to 1.89; two-proportion z-test p = 0.011).

⚠️ **Important for trials:** the Wald CI for a risk difference can behave poorly with small samples or rates near 0% or 100%. Regulatory analyses often use the Miettinen–Nurminen or Newcombe interval, and often adjust for stratification factors (for example with [CMH](#/test/cmh)). State the interval method in the analysis plan.

## Check your understanding

<details class="check"><summary>1. Twelve patients are tested before and after a treatment for a binary outcome (infected / not infected). Is the two-proportion z-test appropriate?</summary>

No. The two measurements come from the same patients, so the groups are not independent. Use [McNemar's test](#/test/mcnemar) on the discordant pairs.
</details>

<details class="check"><summary>2. A colleague says "we can use the z-test because n is above 30 in each arm". What should be checked instead?</summary>

The expected counts in the 2 × 2 table, not the raw n. With 200 patients but a response rate of 1%, expected counts would be tiny and the normal approximation poor. Check that n·p̂ and n·(1 − p̂) are at least 5 in each arm (some authors use 10); otherwise use [Fisher's exact test](#/test/fisher-exact).
</details>

<details class="check"><summary>3. R's default `prop.test(resp, n)` gives p = 0.016, but our write-up says p = 0.011. Which is right?</summary>

Both are valid calculations of different things. R applies the Yates continuity correction by default (p = 0.016); `correct = FALSE` gives the uncorrected z-test (p = 0.011), which matches statsmodels. State which you used.
</details>

<details class="check"><summary>4. Suppose the response rates stayed at 60% and 42% but each arm had only 50 patients. What happens to z?</summary>

The standard error grows by about √2, so z falls from 2.55 to about 1.8 and p rises to about 0.07. The effect (18 points) is the same, but the evidence is weaker: the CI would also be wider. Note that a non-significant result here would not mean the drug has no effect.
</details>

## Common mistakes

- **Letting R apply Yates silently.** `prop.test` and `chisq.test` use continuity correction on 2 × 2 tables by default, while `proportions_ztest` does not. Set `correct = FALSE` (R) or note the setting so results can be reproduced.
- **Reporting only the p-value.** Give the risk difference (absolute effect) and its CI. The OR alone overstates the effect when the outcome is common.
- **Mixing up the direction.** The difference is group 1 − group 2, in the order you supplied. Label it.
- **Confusing the odds ratio with the risk ratio.** "Twice the odds" is not "twice as likely".
- **Using the test with small expected counts.** Switch to [Fisher's exact test](#/test/fisher-exact).
- **Treating a CI and a p-value as always consistent.** The test uses a pooled SE and the CI an unpooled SE, so near the boundary one can cross 0 while the other does not.
- **Analysing paired data this way.** Before/after on the same patients needs [McNemar's test](#/test/mcnemar).

## Related tests

- [Chi-square test of independence](#/test/chi-square-independence): the same test for tables larger than 2 × 2
- [Fisher's exact test](#/test/fisher-exact): small samples
- [McNemar's test](#/test/mcnemar): paired binary data
- [Cochran–Mantel–Haenszel test](#/test/cmh): stratified 2 × 2 tables
- [Logistic regression](#/test/logistic-regression): covariate-adjusted odds ratios
- [One-sample proportion test](#/test/one-proportion): one group against a target rate

## Key takeaways

- The z statistic is the difference in proportions over its standard error; z² equals the 2 × 2 χ² (without Yates).
- Report the **risk difference with a CI** first; add RR and OR as secondary measures.
- The test needs adequate **expected counts** (≥ 5 per cell), not just a large total n.
- R's `prop.test` applies Yates by default; state your continuity-correction choice.
- Paired binary data need McNemar; small counts need Fisher.

## References

- Newcombe RG. Interval estimation for the difference between independent proportions: comparison of eleven methods. *Stat Med.* 1998;17:873–890.
- Agresti A. *Categorical Data Analysis.* 3rd ed. Wiley; 2013. Ch. 2.
- ICH E9: Statistical Principles for Clinical Trials (1998)
- R documentation: `?prop.test`; statsmodels documentation: `proportions_ztest`, `confint_proportions_2indep`
- SAS/STAT User's Guide: The FREQ Procedure — RISKDIFF and RELRISK options
