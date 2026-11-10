# Chi-square goodness-of-fit test

> Tests whether the observed counts across the categories of one variable match a set of expected proportions that were specified in advance.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Compare observed category counts with proportions fixed in advance | One categorical variable (counts) | Independent subjects; expected counts mostly at least 5 | Cohort ABO blood groups vs population frequencies |

## Learning objectives

- Recognize a one-variable question that compares counts with a pre-specified distribution
- Compute expected counts and Pearson residuals by hand
- Run the test in R and Python and read χ², df, and p
- Explain why p > 0.05 does not show that the data "fit" the reference

## What it is

You classify n subjects into k categories and compare the counts you see with the counts you would expect if a reference distribution were true (for example, the known blood-group frequencies of the general population). The statistic adds up (observed − expected)² / expected over the categories. If the observed pattern is close to the reference, χ² is small.

- **Test statistic:** χ² = Σ (Oᵢ − Eᵢ)² / Eᵢ, with Eᵢ = n · pᵢ (pᵢ the hypothesised proportion)
- **Degrees of freedom:** k − 1, minus one more for each parameter estimated from the same data (for example, the allele frequency in a Hardy–Weinberg check)
- **Effect size:** compare **observed vs expected percentages** category by category (primary); Pearson residuals (O − E)/√E show which categories deviate. Cohen's w = √(χ²/n) is a standardised summary.

## The intuition

The test is a **signal-to-noise** check on the whole pattern of counts. For each category, the signal is the gap between observed and expected counts. The noise is the usual random wobble of a count, roughly √expected. Dividing the gap by √E gives the **Pearson residual**, a z-score for that category. Squaring and adding the residuals gives χ².

An everyday picture: you roll a die 60 times and expect about 10 of each face. Getting 8 or 12 of one face is unremarkable. Getting 25 sixes is not. The residual tells you how many "wobbles" away each face is.

In our blood-group example, the largest residual is +1.71 for group B (30 observed vs 22 expected). The others are small. Residuals within about ±2 are ordinary noise; the plot marks 0 as "observed = expected".

<div class="viz" data-viz="dots" data-unit="Pearson residual (O, A, B, AB)" data-groups='{"Residuals":[-0.84,-0.22,1.71,0.71]}' data-zero="1" data-caption="Four residuals for O, A, B, AB. None is beyond 2, so no single blood group stands out."></div>

The χ² statistic is 4.17 with 3 df. If the cohort truly matched the population, χ² values would follow the curve below. Ours sits near the middle, and the shaded upper-tail p-value is large (0.24).

<div class="viz" data-viz="chisq" data-stat="4.17" data-df="3"></div>

Notice that only the **upper** tail is shaded: a χ² near 0 means a very close match, not evidence against H₀.

## When to use it

Use the goodness-of-fit test when **all** of these are true:

- One **categorical** variable with k ≥ 2 categories
- A **reference distribution fixed in advance** (from the population, a genetic model, or a protocol assumption)
- Independent subjects, and **expected counts of at least 5** in nearly all categories

**Clinical scenario.** A transfusion-medicine centre enrols 200 patients in a registry and wants to check that the cohort is representative of the local population, where ABO blood-group frequencies are O 45%, A 40%, B 11%, and AB 4%. *Does the cohort's blood-group mix differ from the population?*

**Don't use it when:**

- You have **two categorical variables** and want to test for association → use the [chi-square test of independence](#/test/chi-square-independence)
- There are **two categories only** → use the [one-sample proportion test](#/test/one-proportion), which is exact
- Expected counts are small in several categories → combine categories, or use an exact multinomial test
- The variable is continuous and you want to test its distribution → use the [Shapiro–Wilk test](#/test/shapiro-wilk) for normality

## Assumptions

1. **Observations are independent**, and each subject falls in exactly one category. Check from the design.
2. **The hypothesised proportions are specified before looking at the data** and sum to 1. Check the protocol; R's `chisq.test` stops with an error if they do not sum to 1 (unless `rescale.p = TRUE`).
3. **Expected counts are large enough.** Aim for every Eᵢ to be at least 5 (a common minimum is that none falls below 1 and no more than about 20% fall below 5). Compute them as n · pᵢ. Here the smallest is 8 (AB).
4. **The reference proportions were not estimated from these data.** If they were, reduce the degrees of freedom accordingly (see "Common mistakes").

## Hypotheses

- **H₀:** the true category proportions equal the specified ones (p_O = 0.45, p_A = 0.40, p_B = 0.11, p_AB = 0.04)
- **H₁:** at least one category proportion differs from its specified value

## Example data

| Blood group | O | A | B | AB | Total |
|---|---|---|---|---|---|
| Observed | 82 | 78 | 30 | 10 | 200 |
| Population proportion | 0.45 | 0.40 | 0.11 | 0.04 | 1.00 |
| Expected (n × p) | 90 | 80 | 22 | 8 | 200 |

## How to run it

### R

```r
# Observed counts (O, A, B, AB) and hypothesised population proportions
obs <- c(O = 82, A = 78, B = 30, AB = 10)
p0  <- c(0.45, 0.40, 0.11, 0.04)

res <- chisq.test(obs, p = p0)
res
res$expected        # expected counts
res$residuals       # Pearson residuals (O - E) / sqrt(E)
```

**What the code does**

- `chisq.test(x, p = ...)` with a **vector** of counts runs the goodness-of-fit test. The default `p` is equal proportions, which is almost never what you want, so always supply it.
- `p` must sum to 1 (a check R performs). Pass percentages divided by 100.
- The degrees of freedom are k − 1 = 3. If you estimated parameters from the data, this df is wrong and you must compute the p-value yourself with `pchisq(stat, df, lower.tail = FALSE)`.
- The Yates continuity correction is **not** applied in the goodness-of-fit setting, so there is no hidden default to worry about here.
- `res$residuals` are the Pearson residuals; `res$stdres` are adjusted and slightly larger.

**Output**

```text
	Chi-squared test for given probabilities

data:  obs
X-squared = 4.1702, df = 3, p-value = 0.2437

 O  A  B AB 
90 80 22  8

         O          A          B         AB 
-0.8432740 -0.2236068  1.7056057  0.7071068
```

### Python

```python
import numpy as np
from scipy.stats import chisquare

obs = np.array([82, 78, 30, 10])             # O, A, B, AB
p0  = np.array([0.45, 0.40, 0.11, 0.04])

# SciPy needs expected COUNTS, not proportions
expected = p0 * obs.sum()
res = chisquare(obs, f_exp=expected)         # df = k - 1 by default (ddof=0)

print(f"chi2 = {res.statistic:.4f}, df = {len(obs)-1}, p = {res.pvalue:.4f}")
print("Expected:", expected)
print("Pearson residuals:", ((obs - expected) / np.sqrt(expected)).round(2))
```

**What the code does**

- `chisquare(f_obs, f_exp)` expects **expected counts** in `f_exp`, not proportions. If you pass proportions, the sums of observed and expected will not match and SciPy will raise an error (or, in older versions, return a wrong answer).
- Without `f_exp`, SciPy assumes **equal** expected counts in every category.
- The `ddof` argument adjusts the degrees of freedom: use `ddof=1` if one parameter was estimated from the data, and so on. The default `ddof=0` gives df = k − 1.
- The residuals are computed by hand to see which categories contribute most.

**Output**

```text
chi2 = 4.1702, df = 3, p = 0.2437
Expected: [90. 80. 22.  8.]
Pearson residuals: [-0.84 -0.22  1.71  0.71]
```

### SAS

```sas
data abo;
  input bloodgroup $ count;
  datalines;
O  82
A  78
B  30
AB 10
;
run;

proc freq data=abo order=data;
  weight count;
  tables bloodgroup / chisq testp=(0.45 0.40 0.11 0.04);
run;
```

**What the code does**

- `weight count;` treats each row as `count` subjects.
- `order=data` keeps O, A, B, AB in the order entered. **The order of `testp=` values must match the order of the categories**, so this option is essential.
- `chisq testp=(...)` requests the goodness-of-fit test against the given proportions (use `testf=(...)` for expected frequencies instead). Without `testp=`, SAS tests equal proportions.
- The proportions must be nonnegative and sum to 1 (percentages summing to 100 are also accepted); check your SAS version's documentation for how it treats values that do not sum correctly.
- The output table shows observed, expected (test) percentages and the cumulative frequencies, followed by the χ² statistic, DF, and p-value.

## Interpreting the output

- **χ²(3) = 4.17, p = 0.24.** There is no evidence that the cohort's blood-group distribution differs from the population. R and Python give identical results.
- **Observed vs expected counts:** O 82 vs 90, A 78 vs 80, B 30 vs 22, AB 10 vs 8. In percentages: O 41.0% vs 45%, A 39.0% vs 40%, B 15.0% vs 11%, AB 5.0% vs 4%.
- **Residuals** are all within ±1.71, with the largest for group B (+1.71, more than expected). None is beyond ±2, so no single category stands out.
- **A non-significant result is not proof of agreement.** With n = 200 the test has limited power to detect modest departures such as a few percentage points.

**How to report it:**
> The blood-group distribution of the cohort (O 41.0%, A 39.0%, B 15.0%, AB 5.0%; n = 200) was not significantly different from the population distribution (χ² goodness-of-fit, df = 3, χ² = 4.17, p = 0.24).

## Check your understanding

<details class="check"><summary>1. You compare the sex ratio of 300 trial participants with the 50:50 split in the target population. Is a goodness-of-fit test needed, or something simpler?</summary>

With only two categories, the exact [one-sample proportion test](#/test/one-proportion) is simpler and more precise. The χ² test would give the same conclusion in large samples.
</details>

<details class="check"><summary>2. The p-value is 0.24. The registry team writes "the cohort matches the population." What is wrong?</summary>

A large p-value only means no departure was detected. With n = 200 the test cannot detect modest departures, such as blood group B being 15% instead of 11%. Report the observed vs expected percentages and judge whether differences like these matter for the study.
</details>

<details class="check"><summary>3. What if the counts were doubled (164, 156, 60, 20; n = 400) with the same percentages? What happens to the test?</summary>

Each gap doubles and each expected count doubles, so every (O − E)²/E term doubles: χ² = 2 × 4.17 = 8.34, and the p-value falls to about 0.039. The same percentage departure becomes significant with more data, which is why effect sizes (percentages, Cohen's w) matter alongside p.
</details>

<details class="check"><summary>4. A Hardy–Weinberg check uses an allele frequency estimated from the same patients. How many degrees of freedom do you use for 3 genotypes?</summary>

df = 3 − 1 − 1 = 1. One degree of freedom is lost for the estimated parameter. R's `chisq.test` does not do this for you, so compute the p-value from the statistic yourself, or use `ddof=1` in SciPy.
</details>

## Common mistakes

- **Forgetting to give the expected proportions.** The default is equal proportions, which will almost always be rejected for the wrong reason.
- **Passing proportions where counts are needed** (SciPy `f_exp`) or counts where proportions are needed (R `p`).
- **Not adjusting df when parameters are estimated.** For a Hardy–Weinberg check with the allele frequency estimated from the data, df = 3 genotypes − 1 − 1 = 1, not 2. In SciPy use `ddof=1`; in R compute the p-value from the statistic.
- **Low expected counts.** With expected counts below 5 in several categories, merge categories or use an exact test.
- **Using percentages as the data.**
- **Mixing up the order** of categories in the data and in the expected proportions (especially in SAS `testp=`).
- **Treating "p > 0.05" as proof the data fit the model.**

## Related tests

- [Chi-square test of independence](#/test/chi-square-independence): two categorical variables
- [One-sample proportion test](#/test/one-proportion): the exact two-category case
- [Shapiro–Wilk test](#/test/shapiro-wilk): checking normality of a continuous variable
- [Fisher's exact test](#/test/fisher-exact): exact methods for sparse tables

## Key takeaways

- Goodness-of-fit compares observed counts with expected counts n × p from a reference fixed **in advance**.
- χ² = Σ (O − E)² / E, with df = k − 1 (minus one per estimated parameter).
- Pearson residuals (O − E)/√E show *which* categories deviate; about ±2 or beyond is notable.
- Expected counts should be at least 5 in nearly all categories.
- A large p-value means no detected departure, not proof of a good fit.

## References

- Agresti A. *Categorical Data Analysis.* 3rd ed. Wiley; 2013. Ch. 1.
- Pearson K. On the criterion that a given system of deviations from the probable in the case of a correlated system of variables is such that it can be reasonably supposed to have arisen from random sampling. *Philos Mag.* 1900;50:157–175.
- R documentation: `?chisq.test`; SciPy documentation: `scipy.stats.chisquare`
- SAS/STAT User's Guide: The FREQ Procedure — TESTP= and TESTF= options
