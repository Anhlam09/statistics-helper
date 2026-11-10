# Bootstrap and permutation tests

> Resampling methods that get a p-value or confidence interval by recomputing the statistic on shuffled or resampled copies of your own data, instead of relying on a formula or a normal distribution.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Small or odd samples, or statistics without a simple formula (median, ratio) | Any; here continuous and skewed | Observations independent; groups exchangeable under H₀ (permutation); sample representative (bootstrap) | Length of stay in two small groups; CI for a difference in medians |

## Learning objectives

- Explain how a **permutation test** builds the null distribution by shuffling group labels
- Explain how the **bootstrap** estimates sampling variability by resampling with replacement
- Run a permutation test and a bootstrap CI (percentile and BCa) in R and Python
- Choose between them, and recognize what each can and cannot do

## What it is

**Permutation test (a test of H₀).** If the treatment has no effect, the group labels are arbitrary: any patient could just as well have been in the other group. So shuffle the labels many times, recompute the statistic (here the difference in means) each time, and see how often a shuffled difference is as extreme as the observed one. That proportion is the p-value. With few observations you can list all possible shuffles (**exact**); otherwise you sample random shuffles (**Monte Carlo**). The Monte Carlo p-value is (count of shuffles as extreme + 1) / (number of shuffles + 1).

**Bootstrap (an estimate of uncertainty).** Treat your sample as a stand-in for the population. Draw new samples **with replacement**, of the same size, recompute the statistic each time, and use the spread of those values as the sampling distribution. For a CI:

- **Percentile CI:** the 2.5th and 97.5th percentiles of the bootstrap statistics.
- **BCa (bias-corrected and accelerated):** adjusts the percentiles for bias and skewness in the bootstrap distribution. It is usually more accurate than the plain percentile CI, but needs enough data.

**Effect size to report:** the difference in means or medians with a CI (bootstrap), together with the permutation p-value if you test.

## The intuition

A permutation test asks one question: "If group membership meant nothing, how often would shuffling produce a gap this big?" Picture the 30 observed stays written on 30 cards. Deal them at random into a pile of 15 labelled "new" and a pile of 15 labelled "standard", compute the gap in pile means, and repeat 9999 times. The histogram of those gaps is centred at 0 (here with SD about 1.5 days, and 95% of the shuffles fall between −3.0 and +2.9 days). Our real gap of −3.3 days lies just outside that central band, so shuffled gaps this large are rare (about 3%).

The bootstrap asks a different question: "How much would my statistic wobble if I could repeat the study?" You cannot, so you reuse the data as if they were the population. Each resample re-draws patients with replacement, and the spread of the statistic over 9999 resamples is the standard error.

Here are the raw data. Notice the right skew in the standard-care group (one stay of 21 days and one of 18), which is why we use medians for the bootstrap.

<div class="viz" data-viz="dots" data-unit="days" data-groups='{"New pathway":[5,7,4,9,6,8,5,15,6,4,7,8,12,5,6],"Standard care":[6,9,5,12,8,14,7,21,10,6,11,9,18,7,13]}' data-caption="Hospital stay in days (15 patients per group). Group means are 7.1 vs 10.4 days. The distributions are right-skewed with a few long stays, so means and normal-theory tests can be shaky."></div>

## When to use it

- The sample is **small** or the data are **skewed or odd**, and you do not trust normal-theory formulas
- The statistic has **no simple standard error** (median difference, ratio of means, a custom index)
- You want an answer that **makes few assumptions** and is easy to explain

**Clinical scenario.** A surgical unit compares hospital stay (days) for 15 patients on a new recovery pathway versus 15 on standard care. Stays are right-skewed. *Is the difference in mean stay larger than chance, and what is a plausible range for the difference in median stay?*

**Don't use it when:**

- The data are clearly normal and the sample is reasonable: a [t-test](#/test/independent-t) is simpler and just as good
- The sample is **very small** (n below about 8 per group) for the bootstrap: it cannot invent information you do not have
- Observations are **not independent** (repeated measures, clustered patients): resample whole clusters or use a [mixed model](#/test/mixed-models)
- You want a rank-based test: use [Mann–Whitney](#/test/mann-whitney)
- Data come from a biased sample: resampling reproduces the bias

## Assumptions

1. **Independent observations.** Each patient contributes once. Resampling individuals from clustered data gives wrong answers.
2. **Exchangeability under H₀ (permutation test).** If H₀ is true, labels can be swapped freely. This holds in a randomized trial; in observational data, shuffling does not remove confounding. Note that it also assumes the groups have the same *distribution*, not just the same mean, so with very different variances the permutation test of means can be misleading.
3. **The sample is representative (bootstrap).** The bootstrap treats the sample as the population. It works poorly for extreme statistics (the maximum), tiny samples, and heavy tails.
4. **Enough resamples.** 9999 is typical. The Monte Carlo error of a p-value near 0.03 with 9999 shuffles is about ±0.002.

## Hypotheses

- **Permutation test:** H₀: the group labels are exchangeable, so the two groups have the same distribution (in particular the same mean). H₁: the mean differs (two-sided).
- **Bootstrap:** no hypothesis; the target is the difference in medians (new − standard) with a 95% CI. If the CI excludes 0, that indicates a median difference at the 5% level; it is a different statistic from the permutation test of means.

## Example data

Hospital stay (days) after surgery, 15 patients per arm.

| Group | Stay (days) | Mean | Median |
|---|---|---|---|
| New pathway | 5, 7, 4, 9, 6, 8, 5, 15, 6, 4, 7, 8, 12, 5, 6 | 7.1 | 6 |
| Standard care | 6, 9, 5, 12, 8, 14, 7, 21, 10, 6, 11, 9, 18, 7, 13 | 10.4 | 9 |

## How to run it

### R

```r
library(boot)

# Length of hospital stay (days) after surgery: standard care vs new pathway
std <- c(6, 9, 5, 12, 8, 14, 7, 21, 10, 6, 11, 9, 18, 7, 13)
new <- c(5, 7, 4, 9, 6, 8, 5, 15, 6, 4, 7, 8, 12, 5, 6)
y <- c(new, std); g <- rep(c("new", "std"), each = 15)
obs <- mean(new) - mean(std)
cat(sprintf("Observed mean difference (new - std): %.2f days\n", obs))

# 1. Permutation test: shuffle the group labels 9999 times
set.seed(2024)
perm <- replicate(9999, { s <- sample(g); mean(y[s == "new"]) - mean(y[s == "std"]) })
# Monte Carlo p-value, counting the observed arrangement itself (+1)
p <- (sum(abs(perm) >= abs(obs)) + 1) / (9999 + 1)
cat(sprintf("Permutation p-value (two-sided) = %.4f\n", p))

# 2. Bootstrap CI for the DIFFERENCE IN MEDIANS, resampling within each group
stat <- function(data, i) {
  d <- data[i, ]
  median(d$y[d$g == "new"]) - median(d$y[d$g == "std"])
}
dat <- data.frame(y, g)
b <- boot(dat, stat, R = 9999, strata = factor(dat$g))
cat(sprintf("Observed median difference: %.1f days\n", b$t0))
ci <- boot.ci(b, type = c("perc", "bca"))
cat(sprintf("Percentile 95%% CI: %.1f to %.1f\n", ci$percent[4], ci$percent[5]))
cat(sprintf("BCa 95%% CI:        %.1f to %.1f\n", ci$bca[4], ci$bca[5]))
```

**What the code does**

- `sample(g)` shuffles the labels; `replicate()` repeats it 9999 times. The +1 in numerator and denominator counts the observed arrangement itself, so the p-value can never be exactly 0.
- `set.seed(2024)` fixes the shuffles so you can reproduce the result. (The `coin` package's `oneway_test(..., distribution = approximate(9999))` does the same job in one line.)
- `boot()` with `strata =` resamples **within each group**, which keeps the group sizes fixed; `stat` recomputes the median difference on each resample.
- `boot.ci(type = c("perc", "bca"))` returns both CIs. Entries 4 and 5 are the lower and upper limits.

**Output**

```text
Observed mean difference (new - std): -3.27 days
Permutation p-value (two-sided) = 0.0281
Observed median difference: -3.0 days
Percentile 95% CI: -6.0 to 0.0
BCa 95% CI:        -8.0 to -1.0
```

### Python

```python
import numpy as np
from scipy import stats

# Length of hospital stay (days) after surgery: standard care vs new pathway
std = np.array([6, 9, 5, 12, 8, 14, 7, 21, 10, 6, 11, 9, 18, 7, 13])
new = np.array([5, 7, 4, 9, 6, 8, 5, 15, 6, 4, 7, 8, 12, 5, 6])
obs = new.mean() - std.mean()
print(f"Observed mean difference (new - std): {obs:.2f} days")

# 1. Permutation test: shuffle group labels (Monte Carlo, 9999 resamples)
def mean_diff(x, y, axis=-1):
    return x.mean(axis=axis) - y.mean(axis=axis)
pt = stats.permutation_test((new, std), mean_diff, permutation_type="independent",
                            n_resamples=9999, alternative="two-sided",
                            vectorized=True, random_state=2024)
print(f"Permutation p-value (two-sided) = {pt.pvalue:.4f}")

# 2. Bootstrap CI for the DIFFERENCE IN MEDIANS (each group resampled separately)
def med_diff(x, y, axis=-1):
    return np.median(x, axis=axis) - np.median(y, axis=axis)
print(f"Observed median difference: {np.median(new) - np.median(std):.1f} days")
for method in ["percentile", "BCa"]:
    r = stats.bootstrap((new, std), med_diff, vectorized=True, n_resamples=9999,
                        confidence_level=0.95, method=method, random_state=2024)
    ci = r.confidence_interval
    print(f"{method:10s} 95% CI: {ci.low:.1f} to {ci.high:.1f}")
```

**What the code does**

- `scipy.stats.permutation_test` with `permutation_type="independent"` shuffles the pooled data between the two groups. `vectorized=True` means the statistic accepts an `axis` argument, which makes it fast.
- With `n_resamples=9999` smaller than the number of possible arrangements, it is a Monte Carlo test. If `n_resamples` is at least the total number of arrangements, SciPy enumerates them all and the p-value is **exact**.
- `scipy.stats.bootstrap` resamples each sample separately; `method="percentile"` and `method="BCa"` give the two CIs.
- `random_state` fixes the stream.

**Output**

```text
Observed mean difference (new - std): -3.27 days
Permutation p-value (two-sided) = 0.0320
Observed median difference: -3.0 days
percentile 95% CI: -6.0 to 0.0
BCa        95% CI: -7.0 to 0.0
```

**Why the digits differ.** Both programs use fixed seeds, but R and Python use different random-number generators, so the shuffles and resamples differ. The permutation p-value is 0.028 in R and 0.032 in Python, a gap of about the Monte Carlo error. The percentile CIs agree; the BCa CIs differ (−8 to −1 vs −7 to 0) because medians of whole-day data are very discrete, so the bootstrap distribution is lumpy and the BCa corrections are sensitive to implementation details. Neither is "wrong". More resamples shrink Monte Carlo noise, but not the lumpiness. For comparison, the Welch t-test gives p = 0.030, in line with the permutation result.

## Interpreting the output

- **Observed mean difference: −3.27 days.** Stays on the new pathway were about 3 days shorter on average (7.1 vs 10.4).
- **Permutation p = 0.028 (R) or 0.032 (Python).** Under shuffled labels, a difference at least this large in either direction occurs in about 3% of shuffles. That is evidence against "no difference in means". The two values agree within Monte Carlo error, so report "p ≈ 0.03".
- **Median difference: −3 days (6 vs 9).** The bootstrap CI measures how precisely we know it.
- **Percentile 95% CI: −6.0 to 0.0 days.** The upper limit touches 0, so the evidence for a median difference is borderline.
- **BCa 95% CI: −8 to −1 (R), −7 to 0 (Python).** With whole-day data and only 15 per group, the interval is rough. The two BCa limits disagree across programs, so prefer the percentile interval here and say "about −3 days, with a plausible range from about −6 to 0".
- **Why the CI and the p-value look different.** The p-value is for means; the CI is for medians. Different statistics, different precision; both are legitimate. Prespecify which one is primary.

**How to report it:**
> Mean hospital stay was 7.1 days on the new pathway and 10.4 days on standard care (difference −3.3 days; permutation test with 9999 resamples, two-sided p = 0.03). The median difference was −3 days (bootstrap percentile 95% CI −6 to 0, 9999 resamples, seed 2024).

## Check your understanding

<details class="check"><summary>1. "A permutation test makes no assumptions, so it is valid for any data." What is wrong?</summary>

It still assumes **exchangeability** under H₀: that the labels could be swapped without changing the distribution. That holds in a randomized experiment. In observational data, group differences may be due to confounding, and shuffling labels does not address that. It also assumes independent observations, and tests of means can be misleading when the groups have very different variances.
</details>

<details class="check"><summary>2. A bootstrap resample has the same size as the original sample. Why does it contain only about 63% of the distinct original patients?</summary>

Sampling **with replacement**: some patients are drawn several times, and others not at all. The chance that a given patient is missed in n draws is (1 − 1/n)ⁿ ≈ 0.37, so about 63% of distinct patients appear. This is normal and not a problem; it is how the bootstrap mimics repeated sampling. It does mean that a few extreme patients can drive the results when n is small.
</details>

<details class="check"><summary>3. You rerun the R code with a different seed and the permutation p-value changes from 0.028 to 0.031. Is the test unreliable?</summary>

No. This is **Monte Carlo error**: with 9999 shuffles and p near 0.03, the standard error is roughly √(0.03 × 0.97 / 9999) ≈ 0.0017. Changes in the third decimal are expected. If the conclusion depends on it, increase the number of shuffles (for example 99,999) or compute the exact p-value. Report the number of resamples and the seed.
</details>

<details class="check"><summary>4. What would change if the true gap between the groups were much larger, say 8 days instead of 3.3?</summary>

The observed difference would likely lie far in the tail of the permutation distribution (whose SD is about 1.5 days). Almost no shuffle would beat it, so the Monte Carlo p-value would be at or near the smallest possible, (0 + 1)/(9999 + 1) = 0.0001. Never report "p = 0": report p ≤ 0.0001 (the smallest value 9999 resamples can resolve), not p < 0.0001.
</details>

## Common mistakes

- **Resampling the pooled data for a two-sample bootstrap CI.** Resample **within each group** so group sizes and the group structure are preserved.
- **Reporting p = 0** from a Monte Carlo test. The smallest possible value is 1/(B + 1).
- **Too few resamples** (such as 100). Use at least 9999 for tests and CIs.
- **Forgetting the seed.** Fix and report it so results are reproducible.
- **Treating the bootstrap as a cure for small samples.** With a very small n, it cannot add information that is not there, and percentile CIs tend to be too narrow.
- **Bootstrapping dependent data** (repeated measures, clusters) by resampling individual rows. Resample whole subjects or clusters.

## Related tests

- [Independent two-sample t-test](#/test/independent-t): the parametric counterpart of the permutation test of means
- [Mann–Whitney U test](#/test/mann-whitney): rank-based alternative, an exact permutation test on ranks
- [Confidence intervals](#/test/confidence-intervals): the bootstrap gives CIs without a formula
- [Hypothesis testing and p-values](#/test/hypothesis-testing): the permutation test is a very literal null distribution
- [Fisher's exact test](#/test/fisher-exact): an exact permutation test for 2×2 tables

## Key takeaways

- A permutation test shuffles labels to build the null distribution; the p-value is the share of shuffles at least as extreme as the data (+1 for Monte Carlo).
- The bootstrap resamples with replacement to estimate the variability of any statistic; use percentile or BCa intervals.
- Neither method removes bias or confounding, and neither works for dependent data unless you resample the right units.
- Fix the seed and report the number of resamples; results differ slightly between languages and seeds.
- With small, skewed samples, they are a transparent alternative to formula-based tests.

## References

- Efron B, Tibshirani RJ. *An Introduction to the Bootstrap.* Chapman & Hall; 1993.
- Davison AC, Hinkley DV. *Bootstrap Methods and their Application.* Cambridge University Press; 1997.
- Phipson B, Smyth GK. Permutation p-values should never be zero. *Stat Appl Genet Mol Biol.* 2010;9:Article 39.
- SciPy documentation: `scipy.stats.permutation_test`, `scipy.stats.bootstrap`
