# McNemar's test

> Compares two paired binary outcomes, such as the same patients before and after treatment, by testing whether the patients who switch from "yes" to "no" outnumber those who switch from "no" to "yes".

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Binary outcome measured twice on the same patients (or matched pairs) | Paired binary (yes/no) | Independent pairs; enough discordant pairs for the χ² version, else exact | Nocturnal asthma symptoms at baseline vs week 12 in 100 patients |

## Learning objectives

- Recognize paired binary data and set them out as a 2 × 2 table of pairs
- Explain why only the discordant pairs carry information
- Run the exact and asymptotic versions in R and Python and know their defaults
- Report counts of improvers and worseners with a paired difference and CI

## What it is

When the same patients are classified twice, the data form a 2 × 2 table of **pairs**: yes/yes, yes/no, no/yes, no/no. Patients who give the same answer both times carry no information about a change. Only the **discordant pairs** matter. Let b = patients who go yes → no, and c = patients who go no → yes. Under H₀ (no change in the marginal proportions) each discordant pair is equally likely to fall in b or c, so b ~ Binomial(b + c, 0.5).

- **Test statistic (asymptotic):** χ² = (b − c)² / (b + c), df = 1. The **exact version** uses the binomial distribution directly. With continuity correction, χ² = (|b − c| − 1)² / (b + c).
- **Effect size:** the **difference in paired proportions** (b − c)/n with a 95% CI (primary), and the **odds ratio for discordant pairs** b/c (secondary; it is the conditional odds ratio of change)

## The intuition

McNemar's test boils the question down to a coin-flip. Patients who answered the same both times (25 symptomatic both times, 49 symptom-free both times) tell us nothing about change. Only the 26 patients who **switched** matter. If the programme did nothing, a switcher is equally likely to have improved or worsened, like a fair coin. The signal is how lopsided the switches are (20 improved vs 6 worsened); the noise is how lopsided 26 fair coin flips can be by luck.

Analogy: two doctors independently judge the same 100 X-rays as "abnormal" or "normal". Where they agree, there is nothing to compare. To learn whether one doctor is stricter, count the films where they disagree and ask who was the stricter one more often.

<div class="viz" data-viz="chisq" data-stat="7.54" data-df="1"></div>

**What to notice:** the asymptotic statistic χ² = (20 − 6)² / 26 = 7.54 lies far in the upper tail of the χ² curve with 1 df. This curve is the large-sample approximation (p = 0.006); the exact binomial test (6 or fewer, or 20 or more, out of 26 switchers), which is the result to report here, gives p = 0.0094 and the same message. Only the upper tail is shaded because squaring turns both "too many improvers" and "too many worseners" into a large value.

| | Switched to better | Switched to worse |
|---|---|---|
| Patients | 20 | 6 |

The concordant cells (25 and 49) would not change the p-value even if they were 2500 and 4900. They do matter for the CI, because they set the denominator n.

## When to use it

Use McNemar's test when **all** of these are true:

- The outcome is **binary** and measured **twice on the same subjects** (before/after, test A vs test B, matched case–control pairs)
- Each subject (or matched pair) contributes **one 2 × 2 cell**
- You want to compare the **marginal proportions** (for example prevalence at baseline vs follow-up)

**Clinical scenario.** In 100 patients with asthma, the presence of nocturnal symptoms (yes/no) is recorded at baseline and again after a 12-week self-management programme. At baseline 45 patients have symptoms; at week 12, 31 do. *Did the proportion with symptoms change?* Because the same patients are measured twice, the two proportions are not independent.

**Don't use it when:**

- The two groups are **different people** → use the [two-proportion z-test](#/test/two-proportion-z) or [Fisher's exact test](#/test/fisher-exact)
- There are **three or more** time points → use [GEE](#/test/gee) or Cochran's Q
- You want agreement between raters rather than a change in prevalence → use [Cohen's kappa](#/test/cohens-kappa)
- You need to adjust for covariates → use [GEE](#/test/gee) or conditional logistic regression
- The outcome is continuous → use the [paired t-test](#/test/paired-t)

## Assumptions

1. **Pairs are independent of each other.** One pair per patient; no patient enters twice. Check from the design.
2. **The outcome is binary at both time points**, defined identically.
3. **The pairing is real.** Rows of the table are the first measurement and columns the second for the same patient (or matched pair).
4. **For the asymptotic test, b + c should be at least about 25** (some use 10). With fewer discordant pairs, use the **exact** test. Check b + c from the table; here it is 26.

## Hypotheses

- **H₀:** the probability of switching yes → no equals that of switching no → yes (marginal homogeneity: P(yes at baseline) = P(yes at week 12))
- **H₁:** the two switching probabilities differ (the proportion with symptoms has changed)

## Example data

Paired counts (n = 100 patients).

| Baseline \ Week 12 | Symptoms | No symptoms | Total |
|---|---|---|---|
| Symptoms | 25 | **20** (b) | 45 |
| No symptoms | **6** (c) | 49 | 55 |
| **Total** | 31 | 69 | 100 |

Discordant pairs: b + c = 26.

## How to run it

### R

```r
# Rows: baseline (Yes/No); columns: week 12 (Yes/No); same 100 patients
tab <- matrix(c(25, 20,
                 6, 49), nrow = 2, byrow = TRUE,
              dimnames = list(Baseline = c("Yes", "No"), Week12 = c("Yes", "No")))

mcnemar.test(tab, correct = FALSE)   # asymptotic, no continuity correction
mcnemar.test(tab, correct = TRUE)    # asymptotic with continuity correction (R default)
binom.test(6, 26, p = 0.5)           # exact: c out of b + c discordant pairs

# Difference in paired proportions (baseline - week 12) with Wald CI
b <- 20; c <- 6; n <- 100
d  <- (b - c) / n
se <- sqrt((b + c) - (b - c)^2 / n) / n
cat("Paired difference:", d, " Wald CI", round(d + c(-1, 1) * 1.96 * se, 3), "\n")
```

**What the code does**

- `mcnemar.test(tab)` expects the paired 2 × 2 table. **`correct = TRUE` is the default** and applies the continuity correction, which gives a larger p-value. Set it explicitly.
- R has **no built-in exact McNemar**. The exact test is `binom.test(c, b + c, 0.5)` on the discordant counts. Its p-value is two-sided.
- The concordant cells (25 and 49) do not affect any of these tests.
- The confidence interval for the paired difference uses the standard large-sample variance for paired proportions, (b + c) − (b − c)²/n, divided by n.

**Output**

```text
	McNemar's Chi-squared test

data:  tab
McNemar's chi-squared = 7.5385, df = 1, p-value = 0.00604

	McNemar's Chi-squared test with continuity correction

data:  tab
McNemar's chi-squared = 6.5, df = 1, p-value = 0.01079

	Exact binomial test

data:  6 and 26
number of successes = 6, number of trials = 26, p-value = 0.009355
alternative hypothesis: true probability of success is not equal to 0.5
95 percent confidence interval:
 0.08974011 0.43647510
sample estimates:
probability of success 
             0.2307692 

Paired difference: 0.14  Wald CI 0.044 0.236 
```

### Python

```python
import numpy as np
from statsmodels.stats.contingency_tables import mcnemar

tab = np.array([[25, 20],
                [ 6, 49]])      # rows: baseline Yes/No; columns: week 12 Yes/No

# Exact, asymptotic, and asymptotic with continuity correction
for label, kw in [("Exact (binomial)", dict(exact=True)),
                  ("Asymptotic, no correction", dict(exact=False, correction=False)),
                  ("Asymptotic, continuity corr.", dict(exact=False, correction=True))]:
    r = mcnemar(tab, **kw)
    print(f"{label:30s} statistic = {r.statistic:.4f}, p = {r.pvalue:.5f}")

# Difference in paired proportions with Wald CI
b, c, n = 20, 6, 100
d  = (b - c) / n
se = np.sqrt((b + c) - (b - c)**2 / n) / n
print(f"Paired difference = {d:.3f}, Wald 95% CI ({d-1.96*se:.3f}, {d+1.96*se:.3f})")
```

**What the code does**

- `mcnemar(table, exact=..., correction=...)` from statsmodels. The defaults are **`exact=True`** and **`correction=True`**, so a bare `mcnemar(tab)` runs the exact test, whereas R's bare `mcnemar.test(tab)` runs the corrected asymptotic test. Always state the settings.
- With `exact=True`, the "statistic" is min(b, c) = 6 and the p-value is the two-sided binomial p-value. With `exact=False`, the statistic is χ².
- `correction` only matters when `exact=False`.
- statsmodels does not provide a CI for the paired difference, so it is computed from the formula.

**Output**

```text
Exact (binomial)               statistic = 6.0000, p = 0.00936
Asymptotic, no correction      statistic = 7.5385, p = 0.00604
Asymptotic, continuity corr.   statistic = 6.5000, p = 0.01079
Paired difference = 0.140, Wald 95% CI (0.044, 0.236)
```

### SAS

```sas
data asthma;
  input baseline $ week12 $ count;
  datalines;
Yes Yes 25
Yes No  20
No  Yes 6
No  No  49
;
run;

proc freq data=asthma order=data;
  weight count;
  tables baseline*week12 / agree;
  exact mcnem;
run;
```

**What the code does**

- `weight count;` stands for the number of patients in each cell, and `order=data` puts Yes before No in both dimensions.
- `agree` prints **McNemar's test** (the asymptotic χ², **without continuity correction**) and the simple and weighted kappa. It is the only way to get McNemar's test from PROC FREQ.
- `exact mcnem;` adds the **exact McNemar test** (the binomial test on the discordant pairs).
- SAS reports the asymptotic version without the correction, so compare with R's `correct = FALSE`.
- PROC FREQ does not give the paired difference in proportions directly here; compute (b − c)/n and its CI by hand (as in the R and Python code), or check your SAS version's documentation for a paired risk-difference option.

## Interpreting the output

- **Discordant pairs: 20 patients improved (yes → no), 6 worsened (no → yes).** The 74 concordant patients do not enter the test.
- **Paired difference: 14.0 percentage points (95% CI 4.4 to 23.6).** The proportion with symptoms fell from 45% to 31%. The CI excludes 0.
- **Exact p = 0.0094.** The recommended result here. If improvement and worsening were equally likely, 20 or more of 26 discordant patients (or 6 or fewer) improving would be rare.
- **Asymptotic χ² = 7.54, p = 0.006**, and **with continuity correction χ² = 6.50, p = 0.011**. The three p-values differ, but all lead to the same conclusion. The uncorrected statistic is the smallest p-value, and the correction brings it closer to the exact result.
- **Discordant odds ratio = 20/6 = 3.3.** Patients who changed were more than three times as likely to change in the favourable direction.
- **R and Python agree** for all three versions.

**How to report it:**
> Nocturnal symptoms were present in 45% of patients at baseline and 31% at week 12 (difference 14.0 percentage points, 95% CI 4.4 to 23.6; 20 patients improved and 6 worsened; exact McNemar test p = 0.009).

⚠️ **Important for trials:** a before/after change in a single group does not show that the programme caused it. Regression to the mean, seasonal variation, and natural history can all produce a change. You need a concurrent control group, usually analysed with [logistic regression](#/test/logistic-regression) or [GEE](#/test/gee).

## Check your understanding

<details class="check"><summary>1. A study compares the cure rate in 80 patients given drug A and 80 different patients given drug B. Can you use McNemar's test?</summary>

No. The groups are different patients, so the data are not paired. Use the [two-proportion z-test](#/test/two-proportion-z) or [Fisher's exact test](#/test/fisher-exact).
</details>

<details class="check"><summary>2. A colleague says the 74 patients with no change show the programme did not work for most people. What does McNemar's test say about them?</summary>

Nothing. Concordant pairs do not affect the test, because it asks whether switches go one way more than the other. Whether most patients stay the same is a descriptive point, not part of this hypothesis.
</details>

<details class="check"><summary>3. Suppose only 3 patients improved and 1 worsened (b + c = 4). Which version should you use?</summary>

The exact test. With so few discordant pairs the χ² approximation is poor; the asymptotic statistic would be (3 − 1)²/4 = 1, p = 0.32, but the exact binomial p-value of 0.62 is the valid one. Either way, the data are far too sparse to tell anything.
</details>

<details class="check"><summary>4. What would change if 20 improvers and 6 worseners became 40 and 12, with the same proportion switching favourably?</summary>

The ratio b/c stays about 3.3, but the evidence strengthens: χ² = (40 − 12)²/52 = 15.1. More discordant pairs means less noise, so p falls to about 0.0001.
</details>

## Common mistakes

- **Using the chi-square test of independence on paired data.** It ignores the pairing and tests the wrong hypothesis.
- **Not stating whether the exact or asymptotic version, with or without correction, was used.** Software defaults differ: R corrects, statsmodels is exact, and SAS reports the uncorrected χ².
- **Using the asymptotic test with few discordant pairs.** With b + c below about 25, use the exact test.
- **Putting the table the wrong way round.** Rows and columns must be the two measurements, not groups and outcomes.
- **Treating the concordant cells as informative.** A large "no/no" cell does not add evidence of a change.
- **Reporting only p.** Report the two proportions, the paired difference with its CI, and the discordant counts.
- **Applying it to matched case–control data with several controls per case.** Use conditional logistic regression.

## Related tests

- [Paired t-test](#/test/paired-t): continuous paired outcomes
- [Cohen's kappa](#/test/cohens-kappa): agreement between two paired classifications
- [Two-proportion z-test](#/test/two-proportion-z): independent groups
- [Fisher's exact test](#/test/fisher-exact): independent groups, small counts
- [GEE](#/test/gee): correlated binary outcomes with covariates

## Key takeaways

- Only the **discordant pairs** carry information; concordant pairs do not enter the test.
- Under H₀, each switcher is equally likely to go either way: b ~ Binomial(b + c, 0.5).
- Use the **exact** test when b + c is small (below about 25); state the version used because software defaults differ.
- Report the discordant counts, the paired difference with its CI, and both proportions.
- A before/after change in one group does not prove the intervention caused it.

## References

- McNemar Q. Note on the sampling error of the difference between correlated proportions or percentages. *Psychometrika.* 1947;12:153–157.
- Agresti A. *Categorical Data Analysis.* 3rd ed. Wiley; 2013. Ch. 10.
- R documentation: `?mcnemar.test`; statsmodels documentation: `statsmodels.stats.contingency_tables.mcnemar`
- SAS/STAT User's Guide: The FREQ Procedure — AGREE option and EXACT MCNEM
