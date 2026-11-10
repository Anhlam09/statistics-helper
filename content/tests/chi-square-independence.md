# Chi-square test of independence

> Tests whether two categorical variables are associated, by comparing the counts in an r × c table with the counts expected if the variables were unrelated.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Test association between two categorical variables in an r × c table | Counts in categories | Independent patients; each in one cell; expected counts not too small | Response grade (none/partial/complete) by dose arm |

## Learning objectives

- Explain observed vs expected counts and why large gaps argue against independence
- Compute and interpret χ², its degrees of freedom, Cramér's V, and standardised residuals
- Run the test in R and Python and check the expected counts sensibly
- Report row percentages alongside the test, not only the p-value

## What it is

Cross-tabulate two categorical variables (for example treatment arm and response category). If they were independent, the count in each cell would equal (row total × column total) / grand total. The Pearson χ² statistic sums the squared gaps between observed and expected counts, each scaled by the expected count. A large value means the pattern of counts is far from what independence predicts.

- **Test statistic:** χ² = Σ (O − E)² / E, summed over all cells
- **Degrees of freedom:** (r − 1)(c − 1)
- **Effect size:** **Cramér's V** = √[ χ² / (n · min(r − 1, c − 1)) ], from 0 (no association) to 1. Also report the row percentages (the clinically readable effect), and look at **standardised residuals** to see which cells drive the result. For a 2 × 2 table, report the risk difference or odds ratio instead (see [two-proportion z-test](#/test/two-proportion-z)).

## The intuition

If the dose had nothing to do with response, the response mix would be the same in every arm. Overall, 70 of 210 patients (33%) had no response, so you would *expect* 33% non-responders in each arm: 20 of 60 on placebo, 23.3 of 70 on low dose, and 26.7 of 80 on high dose. The χ² statistic asks, cell by cell, "how far is what we saw from what independence predicts, compared with how much that count would bounce around by chance?" The expected count E in the denominator is the noise scale (a count of E bounces by roughly √E), and (O − E) is the signal.

| Arm | None: observed / expected | Complete: observed / expected |
|---|---|---|
| Placebo | 30 / 20.0 | 10 / 17.1 |
| Low dose | 25 / 23.3 | 17 / 20.0 |
| High dose | 15 / 26.7 | 33 / 22.9 |

What to notice: the big gaps (30 vs 20 on placebo, 15 vs 26.7 on high dose, 33 vs 22.9 on high dose) are exactly the cells with large standardised residuals (+3.2, −3.5 and +3.2), while the low dose row is close to expected. Squaring and summing those scaled gaps across all 9 cells gives χ² = 18.65 on (3 − 1)(3 − 1) = 4 degrees of freedom. If independence held, χ² would typically be around its degrees of freedom (4), and the curve below shows how unlikely 18.65 is.

<div class="viz" data-viz="chisq" data-stat="18.65" data-df="4"></div>

Because we square the gaps, only the upper tail counts: a value near zero means the data agree suspiciously well with independence, not that anything is wrong.

## When to use it

Use the chi-square test of independence when **all** of these are true:

- Both variables are **categorical** (nominal, or ordinal if you accept ignoring the order)
- Each patient falls into **exactly one cell**, and observations are independent
- Expected counts are adequate (a common guide is no expected count below 1 and at most about 20% of cells below 5; a rough guide, not a hard rule)

**Clinical scenario.** A dose-finding study randomises 210 patients with rheumatoid arthritis to placebo (n = 60), low dose (n = 70) or high dose (n = 80). Response at week 24 is graded as none, partial, or complete. *Is the distribution of response grades different across the three arms?*

**Don't use it when:**

- The table is 2 × 2 with small expected counts → use [Fisher's exact test](#/test/fisher-exact)
- You only want to compare two proportions and report a risk difference → use the [two-proportion z-test](#/test/two-proportion-z)
- The same patients are classified twice (paired data) → use [McNemar's test](#/test/mcnemar)
- You have one categorical variable and a hypothesised distribution → use the [chi-square goodness-of-fit test](#/test/chi-square-gof)
- The data are stratified (for example by site) and you want a pooled test → use the [Cochran–Mantel–Haenszel test](#/test/cmh)
- You need to adjust for covariates → use [logistic regression](#/test/logistic-regression)

## Assumptions

1. **Observations are independent.** One row per patient, no repeated measures or matched pairs. Check from the design.
2. **Categories are mutually exclusive and exhaustive.** Each patient is in one and only one cell.
3. **Expected counts are large enough.** Print the expected table (`$expected` in R, `chi2_contingency(...)[3]` in Python, `expected` in SAS) and check that none are very small. The classic Cochran guide (all at least 1 and at least 80% at least 5) is a rough guide and is often conservative, not a hard threshold. If many cells are sparse, merge categories sensibly, or use an exact or Monte Carlo p-value.
4. **Counts, not percentages.** The test needs raw frequencies.
5. **The total sample is fixed or random, but each patient is classified once.** The test is valid for both a fixed-row-total design (a randomised trial) and a cross-sectional sample.

## Hypotheses

- **H₀:** the two variables are independent (the response distribution is the same in every treatment arm)
- **H₁:** the variables are associated (the response distribution differs in at least one arm)

The test is omnibus: a small p-value says the arms differ somewhere, but not where.

## Example data

Counts of patients (n = 210).

| Arm | None | Partial | Complete | Total |
|---|---|---|---|---|
| Placebo | 30 | 20 | 10 | 60 |
| Low dose | 25 | 28 | 17 | 70 |
| High dose | 15 | 32 | 33 | 80 |
| **Total** | 70 | 80 | 60 | 210 |

## How to run it

### R

```r
# Rows: arm; columns: response category
tab <- matrix(c(30, 20, 10,
                25, 28, 17,
                15, 32, 33), nrow = 3, byrow = TRUE,
              dimnames = list(Arm = c("Placebo", "Low dose", "High dose"),
                              Response = c("None", "Partial", "Complete")))

res <- chisq.test(tab)       # Pearson test; no Yates correction for a 3 x 3 table
res
res$expected                 # check the expected counts first
res$stdres                   # which cells contribute most

# Cramer's V (effect size)
sqrt(unname(res$statistic) / (sum(tab) * (min(dim(tab)) - 1)))
```

**What the code does**

- `chisq.test(tab)` on a matrix runs the independence test. If you pass two factor vectors, it cross-tabulates them first.
- **Yates' continuity correction** is applied by default **only to 2 × 2 tables** (`correct = TRUE`). It is not used here because the table is 3 × 3, but it is the main reason R and other software can disagree for 2 × 2 data. Set `correct = FALSE` explicitly for a 2 × 2 table when you want the uncorrected statistic.
- `res$expected` holds the expected counts. R warns "Chi-squared approximation may be incorrect" if any is below 5.
- `res$stdres` are **adjusted standardised residuals**; values beyond about ±2 mark cells that depart from independence (roughly a z score).
- Cramér's V is computed by hand because base R does not provide it (the `DescTools` and `rcompanion` packages do).

**Output**

```text
	Pearson's Chi-squared test

data:  tab
X-squared = 18.65, df = 4, p-value = 0.0009205

           Response
Arm             None  Partial Complete
  Placebo   20.00000 22.85714 17.14286
  Low dose  23.33333 26.66667 20.00000
  High dose 26.66667 30.47619 22.85714

           Response
Arm               None    Partial   Complete
  Placebo    3.2403703 -0.8987170 -2.4152295
  Low dose   0.5175492  0.4019185 -0.9721111
  High dose -3.5167867  0.4458885  3.1904455

[1] 0.2107261
```

### Python

```python
import numpy as np
from scipy.stats import chi2_contingency
from scipy.stats.contingency import association

# Rows: Placebo, Low dose, High dose; columns: None, Partial, Complete
tab = np.array([[30, 20, 10],
                [25, 28, 17],
                [15, 32, 33]])

# Pearson chi-square test of independence
chi2, p, dof, expected = chi2_contingency(tab, correction=False)
print(f"chi2 = {chi2:.4f}, df = {dof}, p = {p:.5f}")
print("Expected counts:\n", expected.round(2))

# Cramer's V
print(f"Cramer's V = {association(tab, method='cramer'):.3f}")
```

**What the code does**

- `chi2_contingency(tab, correction=False)` returns the statistic, p-value, degrees of freedom, and the table of expected counts, in that order.
- `correction=True` is the SciPy **default** and applies the Yates correction **only when df = 1** (a 2 × 2 table). Setting it explicitly to `False` makes the code behave the same whatever the table size, and makes it match R's `correct = FALSE`.
- `association(..., method="cramer")` returns Cramér's V (SciPy ≥ 1.7).
- The expected counts are checked by printing them; there is no automatic warning.

**Output**

```text
chi2 = 18.6503, df = 4, p = 0.00092
Expected counts:
 [[20.   22.86 17.14]
 [23.33 26.67 20.  ]
 [26.67 30.48 22.86]]
Cramer's V = 0.211
```

### SAS

```sas
data ra;
  input arm $ response $ count;
  datalines;
Placebo   None     30
Placebo   Partial  20
Placebo   Complete 10
Low       None     25
Low       Partial  28
Low       Complete 17
High      None     15
High      Partial  32
High      Complete 33
;
run;

proc freq data=ra order=data;
  weight count;
  tables arm*response / chisq expected cellchi2 nopercent norow;
run;
```

**What the code does**

- `weight count;` makes each row stand for `count` patients.
- `order=data` keeps Placebo, Low, High and None, Partial, Complete in the order entered, instead of alphabetical order.
- `chisq` prints the Pearson χ² (the one to report), the likelihood-ratio χ², the Mantel–Haenszel (linear-by-linear) χ² which is relevant for ordered categories, and **Phi and Cramér's V**.
- `expected` prints the expected counts, and `cellchi2` each cell's contribution to χ².
- `nopercent norow` removes cell and row percentages to keep the table short; the column percentages remain. Remove `norow` if you want row percentages across arms.

## Interpreting the output

- **χ²(4) = 18.65, p = 0.0009.** The distribution of response grades differs across the three arms. Both R and Python give identical values.
- **Expected counts range from 17.1 to 30.5**, all well above 5, so the approximation is reliable.
- **Cramér's V = 0.21.** This is a modest-to-moderate strength of association.
- **Look at the percentages.** "None" occurs in 50% of the placebo arm (30/60) but 18.8% of the high-dose arm (15/80). "Complete" occurs in 16.7% of placebo patients and 41.3% of high-dose patients.
- **Standardised residuals.** Placebo has more non-responders than expected (+3.2) and fewer complete responders (−2.4). High dose has fewer non-responders (−3.5) and more complete responders (+3.2). Low dose is close to independence in every cell (all within ±1).
- The test does not use the **order** of the categories. Here a dose–response trend is plausible, and a trend test (the Mantel–Haenszel line in SAS, or an ordinal model) is more powerful if it was pre-specified.

**How to report it:**
> Response grade differed by treatment arm (χ²(4) = 18.65, p < 0.001; Cramér's V = 0.21). Complete response occurred in 16.7% of placebo, 24.3% of low-dose and 41.3% of high-dose patients.

## Check your understanding

<details class="check"><summary>1. A study classifies 200 patients by smoking status (never / former / current) and by whether they develop a postoperative infection (yes / no). Which test is appropriate?</summary>

The chi-square test of independence on a 3 × 2 table (df = 2), provided the patients are independent and the expected counts are not tiny. If the infection data were 2 × 2 with very few events, [Fisher's exact test](#/test/fisher-exact) would be better. Because smoking is ordered, a trend test could also be more powerful.
</details>

<details class="check"><summary>2. One cell of a 3 × 3 table has an expected count of 4.2. Does the "expected count must be at least 5" rule mean the test is invalid?</summary>

No. The "5" is a rough guide, not a law. The rule of thumb often cited is that no expected count is below 1 and at most about 20% of cells are below 5, and even that is conservative: simulations show the Pearson test holds its error rate in many tables with smaller expected counts, especially when n is moderate. If you are worried, compare with an exact or Monte Carlo p-value (`chisq.test(..., simulate.p.value = TRUE)`) or merge sparse categories sensibly. The check protects against a poor approximation, nothing more.
</details>

<details class="check"><summary>3. Someone says "χ² = 18.65 and Cramér's V = 0.21, so the high dose works best." What does the test actually tell you, and what does it not tell you?</summary>

It says the response distribution differs across arms somewhere (omnibus). It does not say which arm is best; the row percentages and standardised residuals do that (complete response 41.3% on high dose vs 16.7% on placebo). It also does not prove the dose is causal beyond what randomisation gives, and V = 0.21 is a modest association, not a measure of how many patients benefit.
</details>

<details class="check"><summary>4. What would happen to χ² if every cell count were doubled (n = 420) with the same percentages?</summary>

χ² would double (to about 37.3), and the p-value would become far smaller. The percentages, and therefore Cramér's V = 0.21, do not change. This is why a p-value alone cannot measure how strong an association is: significance depends on n, V does not.
</details>

## Common mistakes

- **Running the test on percentages or proportions** instead of counts.
- **Ignoring small expected counts.** Check the expected table, not the observed counts.
- **Forgetting the continuity correction on 2 × 2 tables.** R and SciPy apply it by default; SAS PROC FREQ prints both the uncorrected Pearson χ² and the continuity-adjusted χ² for 2 × 2 tables. State which one you report.
- **Treating a significant omnibus test as "all arms differ".** Use standardised residuals or targeted comparisons (with [multiplicity control](#/test/multiple-comparisons)).
- **Applying it to paired or repeated data.** Use [McNemar's test](#/test/mcnemar) or [GEE](#/test/gee).
- **Reading p-value as strength of association.** With large n, tiny differences are significant; report Cramér's V and the percentages.
- **Pooling strata blindly.** If site or another factor drives both variables, Simpson's paradox is possible; see [CMH](#/test/cmh).

## Related tests

- [Two-proportion z-test](#/test/two-proportion-z): the 2 × 2 special case, with a risk difference
- [Fisher's exact test](#/test/fisher-exact): exact alternative for sparse tables
- [Chi-square goodness-of-fit test](#/test/chi-square-gof): one variable against expected proportions
- [Cochran–Mantel–Haenszel test](#/test/cmh): association controlling for a stratifying factor
- [Logistic regression](#/test/logistic-regression): adjusts for covariates

## Key takeaways

- χ² = Σ (O − E)² / E compares observed counts with those expected under independence. df = (r − 1)(c − 1).
- It is an omnibus test. Use row percentages and standardised residuals (beyond about ±2) to see where the association lies.
- Small expected counts weaken the approximation, but "all E ≥ 5" is a rough guide, not a law. Use an exact or simulated p-value when in doubt.
- Doubling the sample doubles χ² but leaves Cramér's V unchanged: report V and percentages, not just p.
- Use counts, not percentages, and only for independent (unpaired) classifications.

## References

- Agresti A. *Categorical Data Analysis.* 3rd ed. Wiley; 2013. Ch. 3.
- Altman DG. *Practical Statistics for Medical Research.* Chapman & Hall; 1991. Ch. 10.
- R documentation: `?chisq.test`; SciPy documentation: `scipy.stats.chi2_contingency`
- SAS/STAT User's Guide: The FREQ Procedure — CHISQ option
