# Fisher's exact test

> Tests for association in a 2 × 2 table using the exact hypergeometric distribution, so it stays valid when counts are small and the chi-square approximation is not.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Test association between two binary variables, exactly, especially with small counts | Two binary variables (2 × 2 table); r × c tables too | Independent patients, each in one cell; margins treated as fixed | Surgical-site infection with vs without prophylaxis (2/15 vs 8/14) |

## Learning objectives

- Explain how the hypergeometric distribution gives an exact p-value without an approximation
- Run Fisher's test in R and Python and distinguish the conditional MLE from the sample odds ratio
- Recognize that Fisher's test is exact for any table size, and not only a fix for tiny 2 × 2 tables
- Report the two proportions, odds ratio, and exact CI together

## What it is

Fix the row and column totals of a table (a 2 × 2 table here, though the idea extends to r × c). Under the null hypothesis of no association, the count in the top-left cell follows a **hypergeometric distribution**, like drawing patients without replacement. Fisher's test lists every table with the same totals, computes its probability, and adds up the probabilities of all tables that are as unlikely as, or less likely than, the observed one. No large-sample approximation is involved.

- **Test statistic:** none in the usual sense. The p-value is computed directly from the hypergeometric probabilities.
- **Effect size:** the **odds ratio** with an exact 95% CI (primary). The odds ratio that R's `fisher.test` reports is the **conditional maximum-likelihood estimate (conditional MLE)**, which differs slightly from the simple sample odds ratio (a·d)/(b·c). Also give the risk difference or the two proportions, since they are easier for clinicians to read.

## The intuition

Imagine 29 patients, 10 of whom get an infection, and 15 of whom are handed prophylaxis. Now suppose prophylaxis does nothing, and the 10 infections fall at random among the 29 patients, like dealing 10 "infection" cards from a deck of 29. On average, 15/29 of the infections (5.2 of 10) would land in the prophylaxis group. We *saw* only 2. Fisher's test asks how often random dealing gives 2 or something even more lopsided.

- **Signal:** how far the observed count (2) is from the expected count (5.2) in the top-left cell.
- **Noise:** the natural scatter of random dealing. It is a hypergeometric distribution with these margins.

The exact probabilities for the number of infections in the prophylaxis arm (given 10 infections, 15 prophylaxis, 29 patients) are:

| Infections in prophylaxis arm | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Probability | 0.0000 | 0.0015 | 0.0157 | 0.0780 | 0.2046 | 0.3001 | 0.2501 | 0.1169 | 0.0292 | 0.0035 | 0.0001 |

Notice that the most likely table has 5 infections in the prophylaxis arm and that 2 is in the thin lower tail. The two-sided p-value adds all tables no more probable than ours (0, 1, 2 on the low side, and 9, 10 on the high side): 0.0157 + 0.0015 + 0.0000 + 0.0035 + 0.0001 = 0.0209 (matching R and Python's 0.02094 after rounding).

For comparison, the chi-square approximation (Pearson, no continuity correction) gives χ² = 6.15 with 1 degree of freedom and an approximate p of 0.013, which is **too small**. That is the approximation overstating the evidence when counts are small.

<div class="viz" data-viz="chisq" data-stat="6.15" data-df="1"></div>

What to notice: **this curve is the large-sample approximation, not Fisher's test.** It shades an upper-tail area of about 0.013, whereas the exact p-value is 0.021. They differ because the χ² curve is smooth and continuous, while the true distribution is a lumpy hypergeometric one with only a few possible tables, so the smooth curve misjudges the tail when counts are small. Both point the same way, but trust the exact 0.021. Fisher's test is also not limited to 2 × 2 tables with tiny counts: it is exact for any counts, and software extends it to larger tables, but with large counts it matches the chi-square test closely and can be slow to compute.

## When to use it

Use Fisher's exact test when **all** of these are true:

- Both variables are **binary** (a 2 × 2 table). For larger r × c tables with sparse cells, the Fisher–Freeman–Halton extension (`fisher.test` in R, `exact` in SAS PROC FREQ) is available
- The groups are **independent**
- Counts are **small or sparse**: a typical trigger is an expected count below 5, or a total n below about 20 to 40. This is a rule of thumb, not a law: Fisher's test is valid for any counts, it simply matters most when the chi-square approximation is shaky

**Clinical scenario.** A surgical unit pilots antibiotic prophylaxis for a high-risk procedure. Of 15 patients given prophylaxis, 2 develop a surgical-site infection. Of 14 patients managed without prophylaxis, 8 do. *Is infection associated with prophylaxis?* The counts are too small to trust a chi-square approximation.

**Don't use it when:**

- Counts are large (all expected at least 5) → the [two-proportion z-test](#/test/two-proportion-z) or [chi-square test](#/test/chi-square-independence) is fine and gives a risk-difference CI
- You want a less conservative exact test → consider the **mid-p** version or Boschloo's test, which use the same exact logic with a smaller penalty
- The data are paired → use [McNemar's test](#/test/mcnemar)
- The data are stratified (by site, say) → use the [Cochran–Mantel–Haenszel test](#/test/cmh)
- You need covariate adjustment → use [logistic regression](#/test/logistic-regression) (with Firth correction if sparse)

## Assumptions

1. **Observations are independent**, and each patient is in exactly one cell. Check from the design.
2. **Both variables are binary** and defined before looking at the data.
3. **Margins are treated as fixed** (the "conditional" inference). This is a mathematical convenience; the test is valid even if only the row totals were fixed by design, as in a randomised trial. A consequence is that the test is **conservative**: its true type I error is typically below 5%.
4. **No assumption about expected counts**, which is the reason to use it. Still, check them (`chisq.test(tab)$expected`) to confirm whether the chi-square test would have been adequate.

## Hypotheses

- **H₀:** the odds ratio is 1 (no association between treatment and infection)
- **H₁:** the odds ratio is not 1 (two-sided)

A one-sided test is possible but must be pre-specified. For two-sided tests, the p-value sums the probabilities of tables no more probable than the observed table (the default in R and SciPy).

## Example data

| Group | Infection | No infection | n | Infection rate |
|---|---|---|---|---|
| Prophylaxis | 2 | 13 | 15 | 13.3% |
| Control | 8 | 6 | 14 | 57.1% |

The smallest expected count is 4.83, just below 5, so the chi-square approximation is questionable.

## How to run it

### R

```r
# Rows: prophylaxis, control; columns: infection, no infection
tab <- matrix(c(2, 13,
                8,  6), nrow = 2, byrow = TRUE,
              dimnames = list(Group = c("Prophylaxis", "Control"),
                              Outcome = c("Infection", "No infection")))

# Expected counts: is the chi-square approximation OK?
chisq.test(tab)$expected

# Fisher's exact test (two-sided), conditional MLE odds ratio and exact CI
fisher.test(tab)

# Risk difference (prophylaxis - control)
2/15 - 8/14
```

**What the code does**

- `fisher.test(tab)` takes the 2 × 2 matrix, with rows as groups and columns as outcomes.
- The **odds ratio it prints is the conditional MLE** (0.1255), not the sample odds ratio (2 × 6)/(13 × 8) = 0.1154. The CI is the exact conditional interval, based on the noncentral hypergeometric distribution.
- `alternative = "two.sided"` is the default. Use `"less"` or `"greater"` only for pre-specified one-sided tests. `conf.int = TRUE` and `conf.level = 0.95` are the defaults for 2 × 2 tables.
- `chisq.test(tab)$expected` shows why Fisher's test is needed; R also warns that the approximation may be incorrect (the warning is expected here).
- The risk difference is computed by hand because `fisher.test` does not report it.

**Output**

```text
             Outcome
Group         Infection No infection
  Prophylaxis  5.172414     9.827586
  Control      4.827586     9.172414

	Fisher's Exact Test for Count Data

data:  tab
p-value = 0.02094
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 0.01003686 0.89024721
sample estimates:
odds ratio 
  0.125534 

[1] -0.4380952
```

### Python

```python
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio

# Rows: prophylaxis, control; columns: infection, no infection
tab = np.array([[2, 13],
                [8,  6]])

# Exact p-value; the "statistic" is the SAMPLE odds ratio
res = fisher_exact(tab, alternative="two-sided")
print(f"p = {res.pvalue:.5f}, sample OR = {res.statistic:.4f}")

# Conditional MLE odds ratio with exact CI (matches R)
cond = odds_ratio(tab, kind="conditional")
ci = cond.confidence_interval(confidence_level=0.95)
print(f"Conditional MLE OR = {cond.statistic:.4f}, 95% CI ({ci.low:.4f}, {ci.high:.4f})")
print(f"Risk difference = {2/15 - 8/14:.3f}")
```

**What the code does**

- `fisher_exact` returns the p-value and the **sample (unconditional) odds ratio** (a·d)/(b·c). This differs from the OR printed by R's `fisher.test`, which is the conditional MLE.
- `odds_ratio(tab, kind="conditional")` (SciPy ≥ 1.10) returns the **conditional MLE** and its exact CI, which match R. `kind="sample"` gives the simple ratio.
- `fisher_exact` does not return a CI, so use `odds_ratio(...).confidence_interval()` for that.
- The p-value is the same in both languages because both sum the probabilities of tables no more likely than the observed one.

**Output**

```text
p = 0.02094, sample OR = 0.1154
Conditional MLE OR = 0.1255, 95% CI (0.0100, 0.8903)
Risk difference = -0.438
```

### SAS

```sas
data ssi;
  input group $ outcome $ count;
  datalines;
Prophylaxis Infection   2
Prophylaxis NoInfection 13
Control     Infection   8
Control     NoInfection 6
;
run;

proc freq data=ssi order=data;
  weight count;
  tables group*outcome / fisher relrisk riskdiff expected;
  exact or;
run;
```

**What the code does**

- `weight count;` makes each row stand for `count` patients; `order=data` keeps Infection as column 1.
- `fisher` prints Fisher's exact test for 2 × 2 tables, with the left-sided, right-sided and **two-sided** p-values. PROC FREQ also prints it automatically for 2 × 2 tables whenever `chisq` is requested; naming it makes the intent clear.
- `exact or;` requests **exact confidence limits for the odds ratio** (the conditional inference).
- `relrisk` gives the sample odds ratio and relative risks, and `riskdiff` the risk difference. The odds ratio line in the standard output is the sample estimate, so expect 0.1154 there and exact limits from the `exact or` table.
- `expected` shows the expected counts so you can document why the exact test was chosen.

## Interpreting the output

- **Infection rates: 13.3% with prophylaxis vs 57.1% with no prophylaxis.** The risk difference is **−43.8 percentage points**.
- **Exact two-sided p = 0.021.** If prophylaxis and infection were unrelated, such an imbalance (or a more extreme one) would occur in about 2% of tables with these margins.
- **Conditional MLE odds ratio = 0.126 (95% CI 0.010 to 0.890).** The odds of infection with prophylaxis are about one-eighth of those without. The CI is very wide, from roughly a 99% reduction in odds to an 11% reduction, because there are only 10 infections.
- **Sample odds ratio = 0.115.** Slightly further from 1 than the conditional MLE. Both are legitimate; say which you report.
- **Smallest expected count = 4.83**, justifying the exact test. R and Python give the same p-value and the same conditional estimate and CI.
- Fisher's test is **conservative**; the p-value is slightly larger than the true probability of such a result, so some analysts prefer the mid-p version for exploratory work.

**How to report it:**
> Surgical-site infection occurred in 2/15 patients (13.3%) with prophylaxis and 8/14 (57.1%) without (risk difference −43.8 percentage points; conditional odds ratio 0.13, exact 95% CI 0.01 to 0.89; Fisher's exact test, two-sided p = 0.021).

⚠️ **Important for trials:** a pilot of 29 patients can show a large, significant effect that does not replicate. The very wide CI is the honest summary. If allocation was not randomised, confounding (for example surgeon, case mix) may explain the difference.

## Check your understanding

<details class="check"><summary>1. A trial has 400 patients per arm, and one cell expected count is 12. Must you use Fisher's test?</summary>

No. With all expected counts comfortably above 5, the chi-square or two-proportion z-test is fine and gives a risk-difference CI directly. Fisher's test would also be valid and give nearly the same p-value, so the choice is not critical, but you must pre-specify it.
</details>

<details class="check"><summary>2. In the infection example R reports an odds ratio of 0.1255, but the hand calculation (2 × 6)/(13 × 8) gives 0.1154. Which is right?</summary>

Both are legitimate estimates. 0.1154 is the sample odds ratio. R's `fisher.test` reports the **conditional MLE** (0.1255), computed given the margins. State which one you report, and give the exact CI that goes with it.
</details>

<details class="check"><summary>3. Fisher's test gives p = 0.12 in a pilot of 20 patients. A colleague writes "prophylaxis has no effect on infection". What is wrong?</summary>

p = 0.12 means the data are not clearly incompatible with no effect. With 20 patients, power is low and the CI for the odds ratio would be very wide, so a large benefit (or harm) cannot be excluded. Report the CI and call the result inconclusive.
</details>

<details class="check"><summary>4. What if the table were 4/30 vs 16/28 instead of 2/15 vs 8/14 (same rates, twice the sample)?</summary>

The rates are the same (13.3% and 57.1%) but the evidence is stronger: the p-value would be much smaller and the odds ratio CI narrower. The same proportion difference with more patients gives more precision.
</details>

## Common mistakes

- **Reading the sample odds ratio and the conditional MLE as the same number.** R and SciPy's two functions report different ones.
- **Choosing Fisher's test only after seeing that chi-square is non-significant.** Decide from the expected counts before analysis.
- **Using it for stratified data** without a stratified method such as [CMH](#/test/cmh).
- **Reporting only p.** Give the two proportions and an effect size with a CI.
- **Interpreting a non-significant result as "no effect".** With small n, power is low and the CI is wide.
- **Ignoring that it is conservative.** The nominal 5% level is rarely reached exactly.
- **Flipping the table** (swapping rows or columns) without noticing that the odds ratio inverts, so 0.13 becomes 7.97.

## Related tests

- [Chi-square test of independence](#/test/chi-square-independence): the large-sample alternative
- [Two-proportion z-test](#/test/two-proportion-z): large-sample risk difference
- [McNemar's test](#/test/mcnemar): paired binary data
- [Cochran–Mantel–Haenszel test](#/test/cmh): stratified 2 × 2 tables
- [Logistic regression](#/test/logistic-regression): adjusted odds ratios

## Key takeaways

- Fisher's exact test conditions on the margins and sums hypergeometric probabilities. It is exact for any counts.
- It matters most when counts are small, but it is not only a 2 × 2 fix: it extends to r × c tables.
- R reports the conditional MLE odds ratio; SciPy's `fisher_exact` reports the sample odds ratio.
- It is conservative: the true type I error is usually below 5%. The mid-p and Boschloo tests are less conservative.
- With small n, a non-significant result is inconclusive, not evidence of no effect.

## References

- Fisher RA. *Statistical Methods for Research Workers.* Oliver & Boyd; 1925.
- Agresti A. *Categorical Data Analysis.* 3rd ed. Wiley; 2013. Ch. 3.
- R documentation: `?fisher.test`; SciPy documentation: `scipy.stats.fisher_exact`, `scipy.stats.contingency.odds_ratio`
- SAS/STAT User's Guide: The FREQ Procedure — FISHER option and EXACT statement
