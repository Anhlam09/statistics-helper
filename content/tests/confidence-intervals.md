# Confidence intervals

> A 95% confidence interval gives a range of effect sizes compatible with your data; it shows both how big the effect looks and how precisely you know it, which a p-value alone cannot.

## Learning objectives

- Give the correct frequentist interpretation of a 95% CI, and avoid the common wrong one
- Explain how CI width depends on sample size and variability
- Use the duality between a CI and a two-sided test
- Read a CI against a minimal clinically important difference (MCID) to judge clinical relevance

## The intuition

Suppose you want the true mean SBP reduction of a new drug in all patients like yours. You can only study a sample, and a different sample would give a different average. The **sampling variability** of the mean is the noise. Look at 20 repeated samples of 10 patients versus 20 repeated samples of 100 patients, drawn from the same population (true mean 140 mmHg, SD 12):

<div class="viz" data-viz="dots" data-unit="mmHg" data-groups='{"20 samples of n=10":[139.2,136.8,135.6,145.7,138.9,136.4,142.7,143.8,143.6,138.7,141.9,139.6,137.6,136.2,137,145.1,148.8,137.6,140.2,138.4],"20 samples of n=100":[141.3,139.7,141.7,138.2,138.6,138.9,141.6,140.4,139.8,137.5,139.2,141.1,140.4,140.4,138.5,139.6,140.7,140.3,140.3,138.7]}' data-caption="Each dot is one sample mean. With n = 100 the means cluster tightly around the truth (140); with n = 10 they scatter widely."></div>

Notice how much tighter the n = 100 means are. The typical wobble of a mean is its **standard error**, SD/√n. A 95% CI wraps a margin of about two standard errors around your estimate:

**CI = estimate ± t × SE**, where t is about 2 for moderate samples (2.05 for 28 degrees of freedom).

### The correct interpretation

The true mean is a fixed number. The interval is what is random: it changes with every sample. "95%" describes the **procedure**: if you repeated the study many times and built an interval each time, about 95% of those intervals would contain the true value.

So you may *not* say "there is a 95% probability that the true effect lies in this particular interval". After the data are in, the interval either contains the truth or it does not; we just do not know which. The practical reading is gentler and still useful: the interval is the set of effect sizes that are reasonably **compatible with the data**.

Let's watch the procedure work. The simulation below runs 100 trials of n = 25 patients where the true mean SBP change is −8 mmHg. Each trial gives a CI. About 95 should contain −8:

<div class="viz" data-viz="dots" data-unit="mmHg" data-groups='{"First 20 sample means (n=25)":[-2.4,-7.8,-6,-9.1,-5.2,-6.4,-7.9,-6.3,-8.6,-8.9,-9.5,-6.7,-7.1,-9.9,-7.7,-9.7,-7.5,-8.4,-9.2,-5]}' data-caption="Sample means from the first 20 simulated trials (true value -8). Every one is an honest estimate; none is exactly right, and the first (-2.4) is far off, yet its wide interval (-8.3 to 3.5) still covers the truth."></div>

The real coverage numbers come from the code below.

## How width depends on n and SD

Half-width = t × SD / √n. That formula says:

- **Bigger n narrows the interval, but with diminishing returns.** To halve the width you need **four times** the patients.
- **Bigger SD (noisier outcome) widens it.** Tighter measurement or a more homogeneous population gives a narrower interval.
- **Higher confidence widens it.** A 99% CI is wider than a 95% CI from the same data.

The table in the output below shows the width for several n and SD values.

## The CI–test duality

A 95% CI and a two-sided test at α = 0.05 are two views of the same calculation:

- If the 95% CI for a difference **excludes 0**, then p < 0.05.
- If it **includes 0**, then p > 0.05.
- Each value inside the interval is one you could not reject at the 5% level; each outside is rejected.

The t distribution below shows why. A statistic of ±2.05 (28 df) leaves exactly 5% in the two tails, and that is the multiplier used to build the 95% interval.

<div class="viz" data-viz="t" data-stat="2.05" data-df="28"></div>

The interval has more to offer than the p-value, though: it also shows the *size* and the *precision*, and it lets you compare the whole range with what matters clinically.

## Reading a CI against the MCID

Statistical significance asks "is the effect different from 0?". The clinical question is "is the effect big enough to matter?". Say the **minimal clinically important difference (MCID)** for a blood-pressure drug is 5 mmHg. Draw the MCID on the picture and read three things: where the estimate is, whether the interval excludes 0, and whether it reaches the MCID.

| Scenario | Estimate (95% CI) | p | Reading |
|---|---|---|---|
| **A.** Huge trial, 5,000 per arm | −1.0 mmHg (−1.5 to −0.5) | 0.00003 | Significant but **trivial**: the whole interval is below the MCID. |
| **B.** Our small trial, 15 per arm | −4.6 mmHg (−11.1 to +1.9) | 0.16 | Not significant, but **compatible with benefit**, even a 10 mmHg drop. Inconclusive, not negative. |
| **C.** Large trial, ideal case (illustrative numbers) | −9 mmHg (−11 to −7) | < 0.001 | Significant **and** clinically important: the whole interval is beyond the MCID. |

You may conclude "no meaningful effect" only when the whole CI sits inside a pre-specified zone of equivalence ([equivalence testing](#/test/tost-equivalence)). Scenario B is the usual way a "negative" small trial fails: the interval is wide and p > 0.05 is read as "no effect". Scenario A is the opposite error: a tiny p is read as "important".

## How to run it

The simulation builds many CIs where we *know* the truth (true mean change −8 mmHg, SD 12), counts how many cover it, tabulates width by n and SD, and recomputes scenarios A and B. R and Python use different random streams, so the 100-interval counts differ slightly between them, as real replications would.

### R

```r
# 1. CI for the drug-vs-placebo trial (same 30 patients as the p-value lesson)
drug    <- c(133, 130, 144, 126, 125, 125, 133, 128, 137, 134, 133, 131, 124, 130, 121)
placebo <- c(133, 140, 122, 124, 154, 129, 156, 147, 139, 128, 130, 136, 122, 137, 126)
tt <- t.test(drug, placebo, var.equal = TRUE)
cat("diff =", round(diff(rev(tt$estimate)), 2), " 95% CI:", round(tt$conf.int, 2), " p =", round(tt$p.value, 3), "\n")

# 2. 100 trials: true mean SBP change = -8 mmHg, SD 12, n = 25 per trial
set.seed(7)
mu <- -8; sd <- 12; n <- 25
ci_one <- function() { x <- rnorm(n, mu, sd); t.test(x)$conf.int }
ci100 <- t(replicate(100, ci_one()))
covers <- ci100[, 1] <= mu & mu <= ci100[, 2]
cat("Of 100 intervals,", sum(covers), "cover the true mean of -8\n")
cat("Misses:", which(!covers), "\n")
ci1e4 <- t(replicate(10000, ci_one()))
cat("Coverage over 10,000 intervals:", mean(ci1e4[, 1] <= mu & mu <= ci1e4[, 2]), "\n")

# 3. Width of a 95% CI for one mean: 2 * t * SD / sqrt(n)
width <- function(n, sd) 2 * qt(0.975, n - 1) * sd / sqrt(n)
w <- outer(c(10, 25, 100, 400), c(8, 12, 16), Vectorize(width))
dimnames(w) <- list(paste("n =", c(10, 25, 100, 400)), paste("SD =", c(8, 12, 16)))
print(round(w, 1))

# 4. Large trial: -1 mmHg difference, SD 12, 5000 per arm (summary statistics)
se <- 12 * sqrt(2 / 5000); d <- -1
cat("Large trial: diff = -1, 95% CI:", round(d + c(-1, 1) * qnorm(0.975) * se, 2),
    " p =", format.pval(2 * pnorm(-abs(d / se)), digits = 2), "\n")
```

**What the code does**

- `t.test(...)$conf.int` returns the 95% CI; the duality means it excludes 0 exactly when p < 0.05.
- `ci_one()` draws a fresh sample and returns its CI; `replicate` repeats this to build many independent "studies".
- `covers` checks whether each interval contains the true value of −8.
- The width table applies 2 × t × SD / √n for each n and SD.
- Scenario A is computed from summary statistics: SE of a difference = SD × √(2/n per arm).

**Output**

```text
diff = -4.6  95% CI: -11.11 1.91  p = 0.159 
Of 100 intervals, 98 cover the true mean of -8
Misses: 36 75 
Coverage over 10,000 intervals: 0.9523 
        SD = 8 SD = 12 SD = 16
n = 10    11.4    17.2    22.9
n = 25     6.6     9.9    13.2
n = 100    3.2     4.8     6.3
n = 400    1.6     2.4     3.1
Large trial: diff = -1, 95% CI: -1.47 -0.53  p = 3.1e-05 
```

### Python

```python
import numpy as np
from scipy import stats

# 1. CI for the drug-vs-placebo trial (same 30 patients as the p-value lesson)
drug    = np.array([133, 130, 144, 126, 125, 125, 133, 128, 137, 134, 133, 131, 124, 130, 121])
placebo = np.array([133, 140, 122, 124, 154, 129, 156, 147, 139, 128, 130, 136, 122, 137, 126])
r = stats.ttest_ind(drug, placebo)
ci = r.confidence_interval(0.95)
print(f"diff = {drug.mean() - placebo.mean():.2f}  95% CI: {ci.low:.2f} {ci.high:.2f}  p = {r.pvalue:.3f}")

# 2. 100 trials: true mean SBP change = -8 mmHg, SD 12, n = 25 per trial
rng = np.random.default_rng(7)
mu, sd, n = -8, 12, 25
def cis(k):                                     # k simulated trials, all at once
    x = rng.normal(mu, sd, (k, n))
    m, se = x.mean(1), x.std(1, ddof=1) / np.sqrt(n)
    h = stats.t.ppf(0.975, n - 1) * se          # half-width of each interval
    return m - h, m + h
lo, hi = cis(100)
covers = (lo <= mu) & (mu <= hi)
print("Of 100 intervals,", covers.sum(), "cover the true mean of -8")
print("Misses:", np.where(~covers)[0] + 1)
lo, hi = cis(10000)
print("Coverage over 10,000 intervals:", np.mean((lo <= mu) & (mu <= hi)))

# 3. Width of a 95% CI for one mean: 2 * t * SD / sqrt(n)
print("Width (mmHg)   SD=8  SD=12  SD=16")
for n_ in (10, 25, 100, 400):
    w = [2 * stats.t.ppf(0.975, n_ - 1) * s / np.sqrt(n_) for s in (8, 12, 16)]
    print(f"n = {n_:<4}      " + "  ".join(f"{v:5.1f}" for v in w))

# 4. Large trial: -1 mmHg difference, SD 12, 5000 per arm (summary statistics)
se, d = 12 * np.sqrt(2 / 5000), -1
print("Large trial: diff = -1, 95% CI:", np.round(d + np.array([-1, 1]) * stats.norm.ppf(0.975) * se, 2),
      " p =", f"{2 * stats.norm.sf(abs(d / se)):.1e}")
```

**What the code does**

- `r.confidence_interval(0.95)` gives the CI for the difference in means (SciPy 1.10 or later).
- `cis(k)` builds k intervals at once from a matrix of simulated data, one row per trial; `ddof=1` gives the sample SD.
- `covers.sum()` counts intervals that contain −8; with only 100 intervals, the count bounces around 95 by chance.
- The 10,000-interval run pins coverage down very close to 95%.

**Output**

```text
diff = -4.60  95% CI: -11.11 1.91  p = 0.159
Of 100 intervals, 94 cover the true mean of -8
Misses: [ 1 19 22 65 68 75]
Coverage over 10,000 intervals: 0.9485
Width (mmHg)   SD=8  SD=12  SD=16
n = 10         11.4   17.2   22.9
n = 25          6.6    9.9   13.2
n = 100         3.2    4.8    6.3
n = 400         1.6    2.4    3.1
Large trial: diff = -1, 95% CI: [-1.47 -0.53]  p = 3.1e-05
```

## Interpreting the output

- **Our trial: −4.6 mmHg, 95% CI −11.1 to +1.9.** The interval contains 0, so p = 0.159 (duality). It is also wide: it covers a trivial harm up to an 11 mmHg benefit.
- **100 intervals: 98 covered in R, 94 in Python.** Both are normal for a nominal 95%; the count of 100 is itself a random number (roughly 95 ± 2).
- **10,000 intervals: 95.2% (R) and 94.9% (Python).** The long-run coverage is 95%, as promised. That is the entire meaning of "95%".
- **Width table:** going from n = 25 to n = 100 (four times the patients) halves the width from 9.9 to 4.8 mmHg at SD 12. Raising SD from 8 to 16 doubles the width at any n.
- **Scenario A:** −1.0 mmHg (−1.47 to −0.53), p = 3×10⁻⁵. Highly significant, but even the best end of the interval, −1.5 mmHg, is far below a 5 mmHg MCID.

**How to report it:**
> Mean SBP was 4.6 mmHg lower on drug than placebo (95% CI −11.1 to +1.9 mmHg). The interval includes both no effect and reductions beyond the 5 mmHg MCID, so the trial is inconclusive.

## Check your understanding

<details class="check"><summary>1. A paper reports a mean reduction of 6 mmHg, 95% CI 2 to 10. Is it correct to say "there is a 95% chance the true effect is between 2 and 10"?</summary>

Strictly, no. The 95% describes how often the *method* captures the truth over repeated studies. For this one interval the truth is either inside or outside. In practice, readers treat it as a range of values compatible with the data, which is fine, but the probability statement belongs to Bayesian credible intervals.
</details>

<details class="check"><summary>2. Two trials of the same drug: Trial 1 gives −3 mmHg (95% CI −6 to 0), Trial 2 gives −3 mmHg (95% CI −12 to +6). Which has the same estimate but stronger evidence, and why?</summary>

Trial 1: its interval is much narrower (more patients or less variability), so the effect is estimated precisely. Trial 2 has the same point estimate but a very wide interval that is compatible with large benefit or large harm. Trial 2 mainly tells you it was too small.
</details>

<details class="check"><summary>3. A 95% CI for a mean difference is (−0.4, 8.2) mmHg. What can you say about the two-sided p-value for H₀: difference = 0, and about a hypothesised difference of 10?</summary>

The interval contains 0, so the two-sided p is above 0.05. The value 10 lies outside the interval, so a test of H₀: difference = 10 would give p < 0.05. This is the CI–test duality: a 95% CI is the set of null values you cannot reject at α = 0.05.
</details>

<details class="check"><summary>4. You want a CI that is half as wide. By what factor must n increase? What if instead SD fell by half?</summary>

Width is proportional to SD/√n. Halving the width needs √n to double, so **n must be four times larger**. Halving the SD achieves the same thing without extra patients: better measurement or a more homogeneous sample is often cheaper than recruiting.
</details>

## Key takeaways

- A 95% CI is built by a procedure that captures the true value in 95% of repeated studies; any single interval either covers it or does not.
- Width shrinks with √n and grows with SD: four times the patients halves the interval.
- A 95% CI excludes 0 exactly when the two-sided p < 0.05, but the CI also shows size and precision.
- Compare the whole interval with the MCID: significant can still be trivial; non-significant can still be compatible with benefit.
- Report the estimate with its CI in every result.

## References

- Gardner MJ, Altman DG. Confidence intervals rather than P values: estimation rather than hypothesis testing. *BMJ.* 1986;292:746-750.
- Greenland S, et al. Statistical tests, P values, confidence intervals, and power: a guide to misinterpretations. *European Journal of Epidemiology.* 2016;31:337-350.
- Spiegelhalter D. *The Art of Statistics.* Pelican; 2019. Ch. 8.
