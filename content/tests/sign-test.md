# Sign test

> The simplest paired test: ignore how big each change is, count how many subjects went up and how many went down, and ask whether that split differs from 50:50.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Paired data where only the direction of change is trustworthy or differences are skewed | Paired continuous or ordinal | Pairs independent; each pair can be classified up, down, or tied | Pain before vs after physiotherapy in 14 patients |

## Learning objectives

- Explain the sign test as a binomial test on the number of positive vs negative differences
- Handle zero differences (ties) correctly and report how many pairs were used
- Run it in R and Python and compare it with the Wilcoxon signed-rank test
- Decide when the sign test's weaker power is an acceptable price for fewer assumptions

## What it is

For each pair, compute the difference (after − before) and record only its **sign**: positive, negative, or zero. Drop the zeros. If the intervention had no effect, each remaining pair is equally likely to be positive or negative, so the number of negatives follows a Binomial(n, 0.5) distribution, where n is the number of non-zero pairs.

- **Test statistic:** the count of negative (or positive) differences. No t or z is needed; the p-value is exact.
- **Effect size:** the proportion of pairs that improved, with an exact (Clopper–Pearson) 95% CI. If you need a magnitude, report the median difference.
- **Hypothesis:** the median of the differences is 0.

## The intuition

Imagine flipping a fair coin for every patient: heads, they improve; tails, they get worse. If the physiotherapy did nothing, you would see roughly half heads. Seeing 12 improvements out of 13 is like getting 12 heads in 13 flips, which is very unlikely with a fair coin.

The sign test is robust because it uses only the direction. A patient who improved by 25 points counts the same as one who improved by 2. That makes it immune to outliers and to skewed or oddly scaled data, but also wasteful: it throws away information about magnitude, which the [Wilcoxon signed-rank test](#/test/wilcoxon-signed-rank) keeps by ranking the sizes.

Look at the changes below. All but two are below the zero line (one is exactly zero, one is slightly above).

<div class="viz" data-viz="dots" data-unit="VAS points, after − before" data-groups='{"Change":[-12,-6,-11,0,-14,2,-25,-7,-14,-12,-9,-15,-9,-14]}' data-zero="1" data-caption="12 of 14 changes are negative, one is zero, one is positive. The sign test only counts these."></div>

## When to use it

- You have **paired** observations (before/after, left/right, matched pairs)
- The data are ordinal, heavily skewed, or have outliers, or only the **direction** of change is reliable (e.g., "better / worse / same" ratings)
- You want the minimum assumptions (no normality, **no symmetry** of the differences)

**Clinical scenario.** Fourteen patients with chronic low back pain complete a physiotherapy programme. Pain is rated on a 0–100 visual analogue scale before and after. One patient improved by a record 25 points, and the VAS differences look skewed. *Did pain generally decrease?*

**Don't use it when:**

- The differences are roughly normal and you want the most power → use the [paired t-test](#/test/paired-t)
- The differences are roughly symmetric and magnitude matters → use the [Wilcoxon signed-rank test](#/test/wilcoxon-signed-rank), which is more powerful
- The outcome is binary → use [McNemar's test](#/test/mcnemar)
- There are three or more conditions → use the [Friedman test](#/test/friedman)
- The groups are independent → use the [Mann–Whitney test](#/test/mann-whitney)

## Assumptions

1. **Pairs are independent of each other.** Check the design.
2. **Each difference can be classified** as positive, negative, or zero. The two measurements must be on at least an ordinal scale, so that "higher" is meaningful.
3. **Every pair has the same chance of a positive difference** (all differences share the same median, which is 0 under H₀). The differences need **not** be symmetric or normal.
4. **Few ties.** Many zero differences reduce the effective n sharply. Report how many pairs were used.

## Hypotheses

- **H₀:** the median difference is 0 (equivalently, P(difference > 0) = P(difference < 0) = 0.5)
- **H₁:** the median difference is not 0 (two-sided)

## Example data

| Patient | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Before | 62 | 55 | 71 | 48 | 66 | 59 | 80 | 52 | 64 | 70 | 58 | 75 | 61 | 67 |
| After | 50 | 49 | 60 | 48 | 52 | 61 | 55 | 45 | 50 | 58 | 49 | 60 | 52 | 53 |
| Difference | −12 | −6 | −11 | 0 | −14 | +2 | −25 | −7 | −14 | −12 | −9 | −15 | −9 | −14 |

Patient 4 has a zero difference and is dropped, leaving n = 13.

## How to run it

### R

```r
# Pain (0-100 VAS) before and after a physiotherapy programme, 14 patients
before <- c(62, 55, 71, 48, 66, 59, 80, 52, 64, 70, 58, 75, 61, 67)
after  <- c(50, 49, 60, 48, 52, 61, 55, 45, 50, 58, 49, 60, 52, 53)
d <- after - before
d
n_pos <- sum(d > 0); n_neg <- sum(d < 0); n_zero <- sum(d == 0)
c(positive = n_pos, negative = n_neg, zero = n_zero)

# Sign test: drop the zero, test the number of decreases against p = 0.5
binom.test(n_neg, n_pos + n_neg, p = 0.5)

# Wilcoxon signed-rank test on the same differences, for comparison
wilcox.test(after, before, paired = TRUE, exact = FALSE)
```

**What the code does**

- `d > 0`, `d < 0`, `d == 0` count the three kinds of pair. Zeros are excluded from the trials.
- `binom.test(x, n, p = 0.5)` is the exact sign test, with `x` = number of decreases and `n` = non-zero pairs. Its CI is for the proportion of pairs that decreased.
- `wilcox.test(..., paired = TRUE, exact = FALSE)` uses the normal approximation because the data contain ties and a zero (R cannot give an exact p-value then).

**Output**

```text
 [1] -12  -6 -11   0 -14   2 -25  -7 -14 -12  -9 -15  -9 -14
positive negative     zero 
       1       12        1 

	Exact binomial test

data:  n_neg and n_pos + n_neg
number of successes = 12, number of trials = 13, p-value = 0.003418
alternative hypothesis: true probability of success is not equal to 0.5
95 percent confidence interval:
 0.6397026 0.9980544
sample estimates:
probability of success 
             0.9230769 

	Wilcoxon signed rank test with continuity correction

data:  after and before
V = 1, p-value = 0.002066
alternative hypothesis: true location shift is not equal to 0
```

### Python

```python
import numpy as np
from scipy import stats

# Pain (0-100 VAS) before and after a physiotherapy programme, 14 patients
before = np.array([62, 55, 71, 48, 66, 59, 80, 52, 64, 70, 58, 75, 61, 67])
after  = np.array([50, 49, 60, 48, 52, 61, 55, 45, 50, 58, 49, 60, 52, 53])
d = after - before
print("d =", d)
n_pos, n_neg, n_zero = (d > 0).sum(), (d < 0).sum(), (d == 0).sum()
print(f"positive = {n_pos}, negative = {n_neg}, zero = {n_zero}")

# Sign test: drop zeros, binomial test on the number of decreases
res = stats.binomtest(int(n_neg), int(n_pos + n_neg), p=0.5)
ci = res.proportion_ci(confidence_level=0.95, method="exact")
print(f"Sign test: {n_neg}/{n_pos + n_neg} decreases, p = {res.pvalue:.5f}, "
      f"95% CI for proportion ({ci.low:.3f}, {ci.high:.3f})")

# Wilcoxon signed-rank (normal approximation, zeros dropped) for comparison
w = stats.wilcoxon(after, before, zero_method="wilcox", method="approx", correction=True)
print(f"Wilcoxon: W = {w.statistic:.1f}, p = {w.pvalue:.5f}")
```

**What the code does**

- `stats.binomtest(k, n, p=0.5)` is SciPy's exact binomial test (older versions: `binom_test`). SciPy has no dedicated sign-test function, and statsmodels' `sign_test` gives the same result.
- `proportion_ci(method="exact")` gives the Clopper–Pearson interval, identical to R's.
- `stats.wilcoxon(..., zero_method="wilcox", method="approx")` drops zeros and uses the normal approximation with a continuity correction, matching `wilcox.test(exact = FALSE)`. W is the smaller of the two rank sums.

**Output**

```text
d = [-12  -6 -11   0 -14   2 -25  -7 -14 -12  -9 -15  -9 -14]
positive = 1, negative = 12, zero = 1
Sign test: 12/13 decreases, p = 0.00342, 95% CI for proportion (0.640, 0.998)
Wilcoxon: W = 1.0, p = 0.00207
```

## Interpreting the output

- **12 decreases, 1 increase, 1 zero.** The zero is dropped, so n = 13.
- **Sign test p = 0.0034.** If pain were equally likely to go up or down, 12 or more of the same direction out of 13 would be seen about 3 times in 1,000 (two-sided).
- **Proportion improved = 92.3% (95% CI 64% to 99.8%).** The CI is wide because n is small, but it excludes 50%.
- **Wilcoxon p = 0.0021.** Slightly smaller than the sign test because it also uses the *sizes* of the changes (the single increase was small, and the big decreases get high ranks). In general the signed-rank test has more power when the differences are roughly symmetric; the sign test needs no symmetry but has less power. The conclusion is the same here.
- **The sign test does not estimate how much pain fell.** Add the median change (−11.5 VAS points) and its interval (for example, from a bootstrap).

**How to report it:**
> Pain decreased in 12 of 13 patients with a non-zero change (92%, 95% CI 64% to 99.8%; one patient had no change); exact sign test p = 0.003. The median change was −11.5 points on the 0–100 VAS.

## Check your understanding

<details class="check"><summary>1. In a study of 20 pairs, 6 differences are exactly zero. How many pairs go into the sign test?</summary>

14. The sign test drops zeros because they give no information about direction. Report the number dropped, since many zeros weaken the test and can make the conclusion fragile.
</details>

<details class="check"><summary>2. A colleague says "the sign test tests whether the mean difference is zero." Is that right?</summary>

No. It tests whether the **median** difference is zero (equivalently, whether positive and negative differences are equally likely). The mean can be dragged by outliers, which the sign test ignores.
</details>

<details class="check"><summary>3. Your differences are clearly skewed (one huge outlier) and symmetry is doubtful. Should you use Wilcoxon or the sign test?</summary>

The signed-rank test assumes the differences are symmetric about the median; with doubtful symmetry its p-value may not be valid. The sign test needs no such assumption, so it is the safer choice, at the cost of some power. You could run both and see if they agree.
</details>

<details class="check"><summary>4. Using our data, what would the sign-test p-value be if only 9 of the 13 pairs decreased (4 increased)?</summary>

Two-sided p = 2 × P(X ≥ 9 | n = 13, p = 0.5) ≈ 0.267. A 9:4 split is easily seen from a fair coin; the evidence is far weaker than the 12:1 observed. (In R: `binom.test(9, 13)`.)
</details>

## Common mistakes

- **Counting zeros as "no change" successes or failures.** They must be dropped (or handled with a specified tie rule) and reported.
- **Calling it a test of means.** It tests the median.
- **Using a sign test when a paired t-test or Wilcoxon would do.** You lose power for no reason if the data meet the assumptions.
- **Using independent groups.** The sign test needs pairs.
- **Reporting only the p-value.** Give the counts, the proportion improved, and the median change.

## Related tests

- [Wilcoxon signed-rank test](#/test/wilcoxon-signed-rank): more powerful if differences are symmetric
- [Paired t-test](#/test/paired-t): parametric alternative
- [McNemar's test](#/test/mcnemar): paired binary outcomes
- [One-sample proportion (exact binomial) test](#/test/one-proportion): the engine behind the sign test
- [Friedman test](#/test/friedman): three or more paired conditions

## Key takeaways

- The sign test counts positive vs negative differences and runs an exact binomial test against 50:50.
- Zero differences are dropped; report how many.
- It needs no symmetry or normality, but it discards magnitude, so it has less power than Wilcoxon or the paired t-test.
- It tests the median difference; add the median change as an effect size.
- Use it when only direction is trustworthy, or when symmetry is doubtful.

## References

- Dixon WJ, Mood AM. The statistical sign test. *J Am Stat Assoc.* 1946;41:557–66.
- Conover WJ. *Practical Nonparametric Statistics.* 3rd ed. Wiley; 1999. Ch. 3.
- R documentation: `?binom.test`
- SciPy documentation: `scipy.stats.binomtest`
