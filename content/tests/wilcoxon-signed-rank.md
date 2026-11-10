# Wilcoxon signed-rank test

> Tests whether the **paired differences** from the same subjects are centred on zero, using the ranks of the differences instead of their raw values, so it tolerates skewed data and outliers.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Two paired measurements when the differences are skewed or n is small | Continuous or ordinal paired outcome | Independent pairs; differences roughly symmetric about their centre | Back-pain VAS before vs after physiotherapy |

## Learning objectives

- Explain how ranking the absolute differences makes the test robust to outliers
- Run the exact test in R and Python and read W⁺, W⁻, the Hodges–Lehmann estimate, and its CI
- State what the test actually tests, and why the symmetry assumption matters
- Choose between the paired t-test, signed-rank test, and sign test

## What it is

The Wilcoxon signed-rank test is the rank-based counterpart of the [paired t-test](#/test/paired-t). For each pair you compute the difference (after − before), drop any zero differences, rank the **absolute** differences from smallest to largest, and then add up the ranks separately for the positive differences (W⁺) and the negative differences (W⁻). If the treatment did nothing, positive and negative differences should be about equally large and W⁺ ≈ W⁻.

Because it uses ranks, a single extreme difference (say −60 mm instead of −29 mm) counts only as "the largest", not as a huge number that drags a mean around.

- **Test statistic:** W⁺ (sum of ranks of positive differences), as reported by R and SAS (SAS prints S = W⁺ − n(n+1)/4). SciPy reports min(W⁺, W⁻) for a two-sided test. W⁺ + W⁻ = n(n+1)/2.
- **Null distribution:** exact (by enumerating all 2ⁿ sign patterns) for small n without ties; a normal approximation, with tie correction, otherwise.
- **Effect size:** the **Hodges–Lehmann (HL) estimate** (median of all Walsh averages (dᵢ + dⱼ)/2 for i ≤ j) with its confidence interval; the matched-pairs rank-biserial correlation r = (W⁺ − W⁻)/(W⁺ + W⁻) as a standardised measure.

## The intuition

Imagine each patient casting a vote: negative if pain fell, positive if it rose. The paired t-test counts how *big* each vote is. The signed-rank test counts big votes more than small ones, but only by **rank**: the largest change is worth 12 points whether it is 28 mm or 280 mm. Then it asks whether the votes for "pain fell" carry far more total weight than the votes for "pain rose". If physiotherapy did nothing, small and large changes would be split about evenly between the two directions.

- **Signal:** the imbalance between W⁺ and W⁻. Here W⁻ = 76 and W⁺ = 2 out of 78 possible points.
- **Noise:** how lopsided that split could be by chance alone if every difference were equally likely to have either sign.

Look at the differences below. Eleven of the twelve dots sit left of zero, and the lone positive one (+3) is small, so it earns only rank 2.

<div class="viz" data-viz="dots" data-unit="mm VAS, after − before" data-groups='{"Change":[-16,-12,-4,-9,-5,-26,-2,-28,-7,-6,3,-19]}' data-zero="1" data-caption="Within-patient pain changes. Almost all fall below zero; the single increase is small."></div>

The exact null distribution of W⁺ for n = 12 has mean 39 and SD about 12.7. The diagram below is the **normal approximation**, not the exact test: it puts our W⁺ = 2 at z ≈ −2.90, which shades a two-sided area of about 0.004. The exact p-value is **0.0015**, less than half of that. They differ because with only 12 pairs the true distribution of W⁺ is discrete and has lighter tails than the normal curve, so the approximation overstates how likely such an extreme split is.

<div class="viz" data-viz="z" data-stat="-2.90"></div>

Notice how far out in the tail the observed value lies: a rank split this one-sided would be rare if the treatment did nothing, and the exact test says it is rarer still.

## When to use it

Use the Wilcoxon signed-rank test when **all** of these are true:

- You have **paired** data: two measurements on the same subjects, or deliberately matched pairs
- The outcome is at least **ordinal** with meaningful differences (VAS score, lab value, count)
- The differences are **roughly symmetric** around their centre, but not necessarily normal
- The sample is small, skewed, or has outliers, so the paired t-test is doubtful

**Clinical scenario.** Twelve patients with chronic low-back pain rate their pain on a 0–100 mm visual analogue scale (VAS) before and after 4 weeks of supervised physiotherapy. A few patients improve dramatically, most modestly, and one gets slightly worse, so the differences are skewed. *Did pain change, and by how much?*

**Don't use it when:**

- The differences are approximately normal and n is reasonable → the [paired t-test](#/test/paired-t) is more efficient and gives a mean difference
- The two groups contain **different people** → use the [Mann–Whitney U test](#/test/mann-whitney)
- There are **three or more** repeated measurements → use the [Friedman test](#/test/friedman)
- The outcome is binary → use [McNemar's test](#/test/mcnemar)
- The differences are **clearly asymmetric**: the test then no longer targets a median shift, and a significant result could come from skewness alone. The **sign test** (which needs no symmetry), or a transformation such as the log for ratios, is a safer, if less powerful, choice

## Assumptions

1. **Pairs are independent of each other.** One patient's difference says nothing about another's. Check from the study design (no clustering, e.g. two eyes of one patient counted as two pairs).
2. **The outcome is measured on at least an ordinal scale, and differences can be ranked.** Check that "a 10 mm change" means roughly the same thing wherever it occurs on the scale.
3. **The distribution of the differences is symmetric about its median under H₀.** This is a real assumption, not a technicality, and it is often forgotten: the test is *not* assumption-free and does not compare medians in general. Check with a histogram or box plot of the differences. Without symmetry, a significant result only shows that the differences are not balanced around zero, not a clean statement about the median. In small samples such as ours (n = 12) you cannot reliably confirm symmetry, so say plainly that it is assumed.
4. **Few ties and few zero differences.** Ties among |d| force a normal approximation in R, and zeros are dropped by default. Count them and report them.

## Hypotheses

- **H₀:** the distribution of the paired differences is symmetric about 0 (the median, or pseudo-median, difference is 0)
- **H₁:** the differences are shifted away from 0 (two-sided; the default for regulators)

## Example data

Pain VAS (mm), 12 patients. Differences = after − before, so negative means pain relief.

| Patient | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Before | 68 | 72 | 55 | 80 | 63 | 75 | 59 | 90 | 66 | 71 | 62 | 84 |
| After | 52 | 60 | 51 | 71 | 58 | 49 | 57 | 62 | 59 | 65 | 65 | 65 |
| Difference | −16 | −12 | −4 | −9 | −5 | −26 | −2 | −28 | −7 | −6 | +3 | −19 |

There are no zero differences and no tied absolute differences, so an exact p-value is possible.

## How to run it

### R

```r
# Pain VAS (mm, 0-100) before and after 4 weeks of physiotherapy, 12 patients
before <- c(68, 72, 55, 80, 63, 75, 59, 90, 66, 71, 62, 84)
after  <- c(52, 60, 51, 71, 58, 49, 57, 62, 59, 65, 65, 65)
d <- after - before

# Exact Wilcoxon signed-rank test with Hodges-Lehmann estimate and 95% CI
wilcox.test(after, before, paired = TRUE, exact = TRUE,
            conf.int = TRUE, conf.level = 0.95)

# Matched-pairs rank-biserial correlation (effect size)
W_plus  <- sum(rank(abs(d))[d > 0])
W_minus <- sum(rank(abs(d))[d < 0])
c(W_plus = W_plus, W_minus = W_minus,
  r_rb = (W_plus - W_minus) / (W_plus + W_minus))
```

**What the code does**

- `wilcox.test(after, before, paired = TRUE)` runs the signed-rank test on `after − before`. `paired = TRUE` is essential: without it R runs the two-sample rank-sum test, which is a different test.
- `exact = TRUE` requests the exact null distribution. R's default is exact only when n < 50 **and** there are no ties or zeros; with ties it silently switches to a normal approximation (with a warning), so set the option explicitly and read the first line of the output ("exact test" vs "with continuity correction").
- `correct = TRUE` (the default) applies a continuity correction, but only in the normal-approximation case; it has no effect on exact tests.
- `conf.int = TRUE` adds the **Hodges–Lehmann estimate** (labelled "(pseudo)median") and its confidence interval.
- The statistic `V` is **W⁺**, the sum of ranks of the positive differences. Here V = 2 means that almost all the rank mass is on the negative (pain-reducing) side.
- The last lines compute W⁺, W⁻ and the rank-biserial correlation by hand.

**Output**

```text

	Wilcoxon signed rank exact test

data:  after and before
V = 2, p-value = 0.001465
alternative hypothesis: true location shift is not equal to 0
95 percent confidence interval:
 -17.5  -4.5
sample estimates:
(pseudo)median 
         -10.5 

    W_plus    W_minus       r_rb 
 2.0000000 76.0000000 -0.9487179 
```

### Python

```python
import itertools
import numpy as np
from scipy import stats

# Pain VAS (mm, 0-100) before and after 4 weeks of physiotherapy, 12 patients
before = np.array([68, 72, 55, 80, 63, 75, 59, 90, 66, 71, 62, 84])
after  = np.array([52, 60, 51, 71, 58, 49, 57, 62, 59, 65, 65, 65])
d = after - before
n = len(d)

# 1. Exact two-sided signed-rank test (no ties or zero differences here)
res = stats.wilcoxon(d, method="exact", zero_method="wilcox")
ranks = stats.rankdata(np.abs(d))
w_plus, w_minus = ranks[d > 0].sum(), ranks[d < 0].sum()
print(f"SciPy statistic = min(W+, W-) = {res.statistic:.0f}, p = {res.pvalue:.5f}")
print(f"W+ = {w_plus:.0f}, W- = {w_minus:.0f}")

# 2. Hodges-Lehmann estimate = median of the Walsh averages (all pairs i <= j)
i, j = np.triu_indices(n)
walsh = np.sort((d[i] + d[j]) / 2)
print(f"Hodges-Lehmann estimate = {np.median(walsh):.1f} mm")

# 3. Exact 95% CI: enumerate the null distribution of W+ (2^12 sign patterns)
null = np.array([(np.array(s) * np.arange(1, n + 1)).sum()
                 for s in itertools.product([0, 1], repeat=n)])
cdf = np.array([(null <= c).mean() for c in range(n * (n + 1) // 2 + 1)])
c = int(np.max(np.where(cdf <= 0.025)[0]))          # largest c with P(W+ <= c) <= 0.025
print(f"95% CI = ({walsh[c]:.1f}, {walsh[-1 - c]:.1f}) mm")

# 4. Matched-pairs rank-biserial correlation
print(f"Rank-biserial r = {(w_plus - w_minus) / (w_plus + w_minus):.2f}")
```

**What the code does**

- `stats.wilcoxon(d, ...)` takes the **differences** (you can also pass two paired arrays). For a two-sided test SciPy reports the statistic as **min(W⁺, W⁻)**, so it gives 2 here, which matches R's V only because W⁺ happens to be the smaller sum; compute W⁺ and W⁻ yourself when comparing packages.
- `method="exact"` uses the exact null distribution. SciPy's default `method="auto"` is exact for n ≤ 50 without zeros or ties and otherwise uses the normal approximation, much like R.
- `zero_method="wilcox"` (the default) drops zero differences, as R does. `"pratt"` keeps them in the ranking, and `"zsplit"` splits them between W⁺ and W⁻. Pick one and state it.
- `correction=False` is SciPy's default for the continuity correction; R's default is `correct = TRUE`. This matters only when the normal approximation is used, so with ties **set both explicitly**.
- SciPy does not return a Hodges–Lehmann estimate here, so the code builds the Walsh averages with `np.triu_indices` and takes their median. The exact CI limits are the c-th smallest and c-th largest Walsh averages, where c is the largest value with P(W⁺ ≤ c) ≤ 0.025 under the null, obtained by enumerating the 2¹² sign patterns. This is the same construction R uses.

**Output**

```text
SciPy statistic = min(W+, W-) = 2, p = 0.00146
W+ = 2, W- = 76
Hodges-Lehmann estimate = -10.5 mm
95% CI = (-17.5, -4.5) mm
Rank-biserial r = -0.95
```

### SAS

```sas
/* One row per patient: pain VAS (mm) before and after physiotherapy */
data pain;
  input patient before after;
  diff = after - before;            /* negative = pain went down */
  datalines;
1 68 52
2 72 60
3 55 51
4 80 71
5 63 58
6 75 49
7 59 57
8 90 62
9 66 59
10 71 65
11 62 65
12 84 65
;
run;

/* Signed-rank test of H0: median difference = 0 (exact p-value for n <= 20) */
proc univariate data=pain mu0=0;
  var diff;
  ods select TestsForLocation;
run;

/* Hodges-Lehmann estimate: median of all Walsh averages (d_i + d_j)/2, i <= j */
data walsh;
  set pain(keep=diff rename=(diff=di)) nobs=n;
  do j = _n_ to n;
    set pain(keep=diff rename=(diff=dj)) point=j;
    walsh = (di + dj) / 2;
    output;
  end;
  keep walsh;
run;

proc means data=walsh median;
  var walsh;
run;
```

**What the code does**

- The data are in **wide format**; `diff = after - before` is computed in the DATA step, so negative means pain relief.
- `proc univariate ... mu0=0` prints "Tests for Location", which includes Student's t, the Sign test and the **Signed Rank** test. SAS reports **S = W⁺ − n(n+1)/4** (here 2 − 39 = −37), not W⁺ itself, together with `Pr >= |S|`. For n ≤ 20 SAS computes the **exact** p-value by default; for larger n it uses the normal approximation. It drops zero differences. The exact p-value should agree with R's 0.0015.
- `PROC UNIVARIATE` does not print a Hodges–Lehmann estimate for paired data (`PROC NPAR1WAY ... HL` works for two independent samples only). The `walsh` DATA step therefore creates all Walsh averages with `set ... point=j`, and `proc means median` returns the HL estimate (−10.5 mm).
- For a confidence interval, use a macro or `PROC IML` to pick the order statistics of the Walsh averages, or use `PROC UNIVARIATE` with `cipctldf` on the differences as a rough, distribution-free interval for the median.

## Interpreting the output

- **W⁺ = 2, W⁻ = 76.** Only the single patient whose pain rose (+3 mm, rank 2) contributes positive ranks; the other eleven contribute 76 of the 78 possible rank points. The two sums are far from equal, which is evidence of a shift.
- **p = 0.0015 (exact, two-sided).** If physiotherapy had no effect on pain, a rank imbalance this extreme would occur in about 1 in 700 samples. R (0.001465) and SciPy (0.00146) agree.
- **Hodges–Lehmann estimate = −10.5 mm, 95% CI −17.5 to −4.5 mm.** The best estimate of the typical within-patient change is a reduction of about 10 mm, and the data are compatible with reductions between roughly 5 and 18 mm. The CI excludes 0, matching the p-value. Compare it with the minimal clinically important difference for VAS pain (often taken as about 10–15 mm in the literature; check your population).
- **Rank-biserial r = −0.95.** Almost all of the rank-weighted evidence points the same way (pain went down). Values near ±1 mean near-complete consistency; 0 means balance.

**How to report it:**
> Pain VAS fell after 4 weeks of physiotherapy (Hodges–Lehmann median change −10.5 mm, 95% CI −17.5 to −4.5; Wilcoxon signed-rank exact p = 0.0015, n = 12; 11 of 12 patients improved).

⚠️ **Important for trials:** a significant before/after change in a single group does not show that the treatment caused it (regression to the mean, natural history, placebo effects). Prefer a controlled design, analysed for example with [ANCOVA](#/test/ancova).

## Check your understanding

<details class="check"><summary>1. A crossover trial measures each patient's HbA1c on drug A and drug B. The differences look strongly right-skewed with n = 15. Which test, and what caveat?</summary>

The data are paired, so a paired method is required. With strong skew the paired t-test is doubtful; the signed-rank test is a candidate, but skew also violates its symmetry assumption. A log-transform of a ratio-type outcome followed by a paired t-test, or the sign test, may be more defensible.
</details>

<details class="check"><summary>2. A colleague says "the signed-rank test makes no assumptions, so I do not need to check anything." Correct?</summary>

No. It needs independent pairs, rankable differences, and symmetry of the differences about their centre. It is robust to outliers and non-normality, but not assumption-free.
</details>

<details class="check"><summary>3. Suppose the +3 mm patient had instead changed by +28 mm, the same size as the biggest improvement (−28 mm). What happens to W⁺ and W⁻?</summary>

The two 28s now tie for the top ranks and share the average rank 11.5, so W⁺ = 11.5 and W⁻ = 78 − 11.5 = 66.5. The imbalance shrinks and p grows (and ties force a normal approximation). A paired t-test would react more strongly, because 28 mm enters its mean and SD at full size, whereas here it only adds a rank.
</details>

<details class="check"><summary>4. The test gives p = 0.30 in a small study. Can you conclude the treatment does not change pain?</summary>

No. A large p-value reflects limited evidence, not absence of an effect. Look at the Hodges–Lehmann estimate and its CI. If the CI includes clinically important reductions, the study is inconclusive.
</details>

## Common mistakes

- **Forgetting `paired = TRUE` in R.** `wilcox.test(after, before)` is the Mann–Whitney test and ignores the pairing.
- **Treating it as a test of means or of the raw medians.** It tests the distribution of *differences*; the difference of the two group medians is not the same as the median of the differences. Report the HL estimate.
- **Comparing the statistic across packages.** R gives W⁺, SciPy gives min(W⁺, W⁻), SAS gives S = W⁺ − n(n+1)/4. They are different numbers for the same test.
- **Ignoring ties and zeros.** Zero differences are dropped (by default) and ties require an approximation. Report how many pairs were used, and which option you chose.
- **Using it on very asymmetric differences.** The symmetry assumption is real; check a plot of the differences.
- **Reporting only the p-value.** Always add an effect estimate with a CI.

## Related tests

- [Paired t-test](#/test/paired-t): parametric alternative, estimates the mean difference
- [Mann–Whitney U test](#/test/mann-whitney): the two-independent-groups analogue
- [Friedman test](#/test/friedman): three or more repeated measurements
- [McNemar's test](#/test/mcnemar): paired binary outcomes
- [Shapiro–Wilk test](#/test/shapiro-wilk): check normality of the differences before choosing

## Key takeaways

- The signed-rank test ranks the absolute paired differences and compares the rank sums of positive (W⁺) and negative (W⁻) differences.
- It is robust to outliers and non-normality but assumes the differences are roughly **symmetric**. If not, consider the sign test.
- Report the Hodges–Lehmann estimate with its CI, not just p.
- R reports W⁺, SciPy reports min(W⁺, W⁻), and SAS reports S; they are the same test.
- Always set `paired = TRUE` in R, or you run the wrong test.

## References

- Hollander M, Wolfe DA, Chicken E. *Nonparametric Statistical Methods.* 3rd ed. Wiley; 2014. Ch. 3.
- Conover WJ. *Practical Nonparametric Statistics.* 3rd ed. Wiley; 1999.
- Documentation: R `?wilcox.test`; SciPy `scipy.stats.wilcoxon`; SAS/STAT User's Guide, The UNIVARIATE Procedure — Tests for Location
