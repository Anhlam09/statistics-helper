# One-sample proportion (exact binomial) test

> Tests whether the response rate observed in a single group differs from a fixed target rate, using the exact binomial distribution, and gives an exact (Clopper–Pearson) confidence interval for the true rate.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Compare one group's event rate with a fixed reference rate | Binary outcome, one group | Independent patients; same probability each; p₀ fixed in advance | 14 of 40 lymphoma patients respond vs a 20% historical rate |

## Learning objectives

- Explain how the binomial distribution under H₀ turns a count into a p-value
- Run the exact binomial test in R and Python and report the Clopper–Pearson CI
- Explain why the two-sided exact p-value can differ between software packages
- Interpret a single-arm result against a historical rate with appropriate caution

## What it is

You count how many of n patients had the event (a "success", such as tumour response) and ask whether the underlying probability is compatible with a reference value p₀. The **exact binomial test** does not rely on a normal approximation: it adds up binomial probabilities directly, so it stays valid for small samples and for rates near 0% or 100%. The two-sided p-value is the total probability of all outcomes that are no more likely than the one you observed, when the true rate is p₀.

- **Test statistic:** X = number of successes, which follows Binomial(n, p₀) under H₀ (there is no z or χ² to look up)
- **Effect size:** the observed proportion p̂ = X / n with its **exact Clopper–Pearson 95% CI** (primary); the difference p̂ − p₀ is a useful secondary summary

## The intuition

Think of a coin that lands heads 20% of the time, and 40 tosses. You would expect about 8 heads (40 × 0.20), give or take 2.5 or so. If you saw 14, you would start to suspect the coin is not a 20% coin. The exact binomial test measures that suspicion directly: how often does a genuinely 20% process give a result as extreme as 14 of 40?

- **Signal:** the gap between what we saw (14 responders) and what H₀ predicts (8).
- **Noise:** the natural binomial scatter, SD = √(40 × 0.2 × 0.8) ≈ 2.5 responders.

Under H₀ the exact probabilities for the number of responders are:

| Responders | 2 | 3 | 4 | 6 | 8 | 10 | 12 | 13 | 14 | 15 |
|---|---|---|---|---|---|---|---|---|---|---|
| P(X = x) | 0.0065 | 0.0205 | 0.0475 | 0.1246 | 0.1560 | 0.1075 | 0.0443 | 0.0238 | 0.0115 | 0.0050 |

Notice that 8 is by far the most likely result and 14 is rare (about 1%). The two-sided p-value adds up every outcome at most as likely as 14: the upper tail P(X ≥ 14) = 0.0194 plus the lower outcomes 0, 1, and 2 (total 0.0079), which gives p = 0.0273.

As a rough picture, the normal approximation puts the observed rate 35% at z = (0.35 − 0.20) / √(0.20 × 0.80 / 40) ≈ 2.37. The approximate p-value (0.018) is smaller than the exact one (0.027), a reminder of why exact methods matter for small n.

<div class="viz" data-viz="z" data-stat="2.37"></div>

What to notice: **this curve is the large-sample normal approximation, not the exact test.** The observed z sits in the right tail, but the shaded two-sided area is about 0.018, whereas the exact p-value is 0.027. They differ because the binomial distribution is discrete and, with n = 40 and p₀ = 0.2, slightly skewed, while the normal curve is smooth and symmetric (a continuity correction, z ≈ 2.17, would give about 0.030). Report the exact 0.027.

## When to use it

Use the exact binomial test when **all** of these are true:

- The outcome is **binary** (responder / non-responder, infected / not infected)
- There is **one group** and a **pre-specified reference rate** p₀ (from historical data, a guideline, or a clinically meaningful threshold)
- Subjects are independent and each has the same probability of success
- A single-arm comparison with a reference rate is **hypothesis-generating at best**: historical rates come from different patients, so treat a significant result as a signal to proceed to a randomized trial, not as proof of benefit

**Clinical scenario.** A single-arm phase II study treats 40 patients with relapsed lymphoma using a new agent. Historical data suggest a 20% objective response rate with standard salvage therapy, and the sponsor wants to know whether the new agent clearly does better. *Is the response rate different from 20%?* Fourteen of the 40 patients respond.

**Don't use it when:**

- There are **two independent groups** to compare → use the [two-proportion z-test](#/test/two-proportion-z) or [Fisher's exact test](#/test/fisher-exact)
- The same patients are measured twice → use [McNemar's test](#/test/mcnemar)
- You have **three or more categories** and expected proportions → use the [chi-square goodness-of-fit test](#/test/chi-square-gof)
- The outcome is a time to event with censoring → use the [Kaplan–Meier estimator](#/test/kaplan-meier)
- You need to adjust for covariates → use [logistic regression](#/test/logistic-regression)

## Assumptions

1. **Each patient is a single independent Bernoulli trial.** No patient is counted twice, and outcomes do not influence each other. Check this from the study design (for example, no clustering by site that drives response).
2. **The success probability is the same for every patient.** In a heterogeneous population the "true rate" is an average. Check by looking at obvious subgroups, or model with [logistic regression](#/test/logistic-regression).
3. **p₀ was fixed before seeing the data.** Choosing the reference rate after the fact invalidates the p-value. Check the protocol or statistical analysis plan (ICH E9).
4. **n is fixed in advance.** If the sample size depended on interim results (for example a two-stage Simon design), the standard p-value and CI are no longer exactly correct; use the design's own analysis.

## Hypotheses

- **H₀:** p = 0.20 (the true response rate equals the reference rate)
- **H₁:** p ≠ 0.20 (two-sided). A one-sided alternative p > 0.20 is common in phase II screening, but it must be pre-specified.

## Example data

Single-arm phase II study, n = 40 patients.

| Outcome | Patients | Percent |
|---|---|---|
| Responder | 14 | 35.0% |
| Non-responder | 26 | 65.0% |
| **Total** | **40** | 100% |

## How to run it

### R

```r
# Phase II single arm: 14 responders out of 40; reference rate 20%
x <- 14
n <- 40

# Exact binomial test with Clopper-Pearson CI
res <- binom.test(x, n, p = 0.20, alternative = "two.sided", conf.level = 0.95)
res

# Observed proportion
x / n
```

**What the code does**

- `binom.test(x, n, p = 0.20)` takes the **number of successes** and the **number of trials**, not the raw 0/1 data. For a 0/1 vector use `binom.test(sum(v), length(v), ...)`.
- `p = 0.20` is the null value p₀. The default is 0.5, which is rarely what you want.
- `alternative = "two.sided"` is the default. Use `"greater"` only if a one-sided test was pre-specified.
- The confidence interval is always the **exact Clopper–Pearson** interval. R's `prop.test` would give a different (Wilson-type, continuity-corrected) interval, so do not mix them up.

**Output**

```text
	Exact binomial test

data:  x and n
number of successes = 14, number of trials = 40, p-value = 0.02735
alternative hypothesis: true probability of success is not equal to 0.2
95 percent confidence interval:
 0.2062825 0.5168445
sample estimates:
probability of success 
                  0.35 

[1] 0.35
```

### Python

```python
from scipy import stats

# Phase II single arm: 14 responders out of 40; reference rate 20%
x, n = 14, 40

# Exact binomial test
res = stats.binomtest(x, n, p=0.20, alternative="two-sided")

# Clopper-Pearson CI (method="exact")
ci = res.proportion_ci(confidence_level=0.95, method="exact")

print(f"Observed proportion = {res.statistic:.4f} ({x}/{n})")
print(f"Exact two-sided p = {res.pvalue:.4f}")
print(f"95% Clopper-Pearson CI = ({ci.low:.4f}, {ci.high:.4f})")
```

**What the code does**

- `stats.binomtest` (SciPy ≥ 1.7) replaces the old `binom_test`. Its arguments are the successes `k`, trials `n`, and null probability `p`.
- `proportion_ci(method="exact")` gives the Clopper–Pearson interval. The default method is also `"exact"`, but naming it makes the intent clear. Other choices such as `"wilson"` give different intervals.
- `res.statistic` is the observed proportion, not a test statistic in the z or t sense.

**Output**

```text
Observed proportion = 0.3500 (14/40)
Exact two-sided p = 0.0273
95% Clopper-Pearson CI = (0.2063, 0.5168)
```

### SAS

```sas
/* One row per outcome category with its frequency */
data resp;
  input response $ count;
  datalines;
Yes 14
No  26
;
run;

proc freq data=resp order=data;
  weight count;
  tables response / binomial(p=0.20 level="Yes" exact) alpha=0.05;
  exact binomial;
run;
```

**What the code does**

- `weight count;` tells PROC FREQ that each row stands for `count` patients. Without it, each row would count once.
- `binomial(p=0.20 level="Yes" ...)` requests the one-sample binomial analysis for the "Yes" category against the null value 0.20. Without `level=`, SAS uses the first category in sorted or data order.
- The `exact` option inside `binomial(...)` adds the **Clopper–Pearson exact confidence limits** to the default Wald limits.
- `exact binomial;` requests the **exact test p-value**. Without it, SAS reports only the asymptotic z test.
- **Two-sided p-values differ by convention.** PROC FREQ's exact two-sided p-value is twice the one-sided tail probability (here 2 × 0.0194 = 0.0388), whereas R and SciPy sum the probabilities of all outcomes no more likely than the observed one (0.0273). Both are valid exact tests; the Clopper–Pearson limits are identical. Report which convention you used.
- `alpha=0.05` gives 95% limits.

## Interpreting the output

- **Observed response rate = 35.0% (14/40).** This is the point estimate.
- **95% Clopper–Pearson CI: 20.6% to 51.7%.** The data are compatible with a true response rate from roughly 21% to 52%. The interval is wide because n = 40, and its lower end is only just above the 20% reference rate.
- **Exact two-sided p = 0.027.** If the true rate were 20%, an outcome at least this unusual (in either direction) would occur in about 2.7% of studies. This is evidence against H₀ at the 5% level, but not overwhelming.
- **Agreement between R and Python.** Both give the same p-value and the same interval, because both use the exact method. SAS gives the same interval but a larger two-sided p-value (0.039) because it doubles the one-sided tail; the one-sided exact p-value (0.0194) is the same in all three.

**How to report it:**
> Fourteen of 40 patients responded (35.0%; exact 95% CI 20.6% to 51.7%), which differed from the 20% reference rate (exact binomial test, two-sided p = 0.027).

⚠️ **Important for trials:** statistical significance against a historical rate is not proof of benefit. Historical controls differ in patient selection, supportive care, and response criteria. In phase II, the decision usually rests on a pre-specified rule (for example, the lower CI bound exceeding a minimum clinically interesting rate), not on p < 0.05 alone.

## Check your understanding

<details class="check"><summary>1. A trial gives 14 of 40 responders and the p-value against 20% is 0.027. Does that mean the new drug raises the response rate above 20%?</summary>

It says the data are unlikely if the true rate were exactly 20%, and the CI (20.6% to 51.7%) agrees. But the 20% is a historical figure. If current patients are healthier or response criteria differ, the "improvement" may not be due to the drug. Randomized comparison is needed to claim benefit.
</details>

<details class="check"><summary>2. Why do R and SAS give different two-sided exact p-values (0.0273 vs 0.0388) for the same data?</summary>

They define "two-sided" differently. R and SciPy add probabilities of all outcomes no more likely than the observed one (0.0194 + 0.0079). SAS doubles the one-sided tail (2 × 0.0194). Both are valid exact tests, so state which you used. The one-sided p-value and the CI are the same.
</details>

<details class="check"><summary>3. What if only 11 of 40 had responded (27.5%)? Would you reject p = 0.20 at the 5% level?</summary>

Almost certainly not. The expected count under H₀ is 8, and P(X ≥ 11) = 0.16, so the two-sided p is about 0.24. The CI would also include 20%. Observing 27.5% against 20% is within what chance produces with n = 40.
</details>

<details class="check"><summary>4. A study of 10 patients has 0 responders. Can you say the true response rate is 0%?</summary>

No. The exact 95% CI for 0/10 runs from 0% to about 31%. Zero events in a small sample is compatible with a real rate as high as roughly 3 in 10. This is the origin of the "rule of three": the upper 95% limit for 0 events in n trials is about 3/n.
</details>

## Common mistakes

- **Using the normal-approximation (Wald) CI with small n or extreme rates.** It can fall below 0 or above 1 and has poor coverage. Use the exact or Wilson interval.
- **Passing proportions or percentages instead of counts.** `binom.test(0.35, 40)` is an error or nonsense; give the number of successes.
- **Leaving the null at 0.5 by default.** Always set p₀ explicitly.
- **Choosing p₀ or the sidedness after seeing the data.**
- **Reading "exact" as "unbiased in coverage".** Clopper–Pearson is conservative: the actual coverage is at least 95%, so the interval is slightly wider than necessary.
- **Ignoring the two-stage design.** If the trial used a Simon two-stage rule, report the design-based analysis.

## Related tests

- [Two-proportion z-test](#/test/two-proportion-z): compares rates in two independent groups
- [Fisher's exact test](#/test/fisher-exact): exact test for two groups with small counts
- [Chi-square goodness-of-fit test](#/test/chi-square-gof): one sample, more than two categories
- [Logistic regression](#/test/logistic-regression): response rate adjusted for covariates

## Key takeaways

- The exact binomial test compares X successes in n trials with Binomial(n, p₀); no approximation is needed.
- Always report the observed proportion with an exact (Clopper–Pearson) CI, not just p.
- Exact two-sided p-values differ by convention across software; the one-sided p and the CI do not.
- p₀ must be fixed before the data are seen, and a win against a historical rate is not proof of benefit.
- Zero events in a small sample does not mean a zero rate.

## References

- Clopper CJ, Pearson ES. The use of confidence or fiducial limits illustrated in the case of the binomial. *Biometrika.* 1934;26:404–413.
- Agresti A. *Categorical Data Analysis.* 3rd ed. Wiley; 2013. Ch. 1.
- R documentation: `?binom.test`; SciPy documentation: `scipy.stats.binomtest`
- SAS/STAT User's Guide: The FREQ Procedure — BINOMIAL option and EXACT statement
