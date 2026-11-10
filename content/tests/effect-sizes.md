# Effect sizes and clinical significance

> A p-value says whether an effect is detectable; an effect size says how big it is, and only the effect size (with its confidence interval) tells you whether it matters to patients.

## Learning objectives

- Distinguish statistical significance from clinical significance, using the MCID
- Report an effect as a raw difference in natural units and as a standardized measure (Cohen's d), knowing the limits of each
- Recognize RR, OR, RD, and NNT as effect measures for binary outcomes
- Always pair an effect size with its confidence interval

## The intuition

A p-value answers "is the effect distinguishable from zero?" It does not answer "how large is it?". A huge trial can make a trivial effect "highly significant", and a tiny trial can make an important effect "non-significant". Cohen (1994) put the problem sharply: the null hypothesis of exactly zero effect is almost never literally true, so rejecting it with enough data is guaranteed and uninformative.

Here are two hypothetical trials of an antihypertensive, both with patient SD 12 mmHg and nearly the same p-value:

- **Trial A:** 1,360 patients per arm, SBP difference −1.0 mmHg, p = 0.030.
- **Trial B:** 17 patients per arm, SBP difference −9.0 mmHg, p = 0.036.

Where do their test statistics fall on the null curve?

<div class="viz" data-viz="z" data-stat="-2.17"></div>

<div class="viz" data-viz="t" data-stat="-2.19" data-df="32"></div>

Notice that both statistics land in about the same place and the shaded tail area is almost identical, so the two p-values are about the same. Yet one trial found a 1 mmHg change, which no patient would notice and which sits far below any plausible clinical threshold, while the other found a 9 mmHg drop, which would matter. The p-value cannot tell them apart. The effect size and its confidence interval can: A has a 95% CI of −1.9 to −0.1 mmHg (a precisely-measured trivial effect); B has −17.4 to −0.6 mmHg (a large but imprecisely-measured effect, compatible with anything from negligible to enormous).

## Statistical vs clinical significance

- **Statistical significance:** the data are hard to explain by chance alone (p < α).
- **Clinical significance:** the effect is large enough to change what a patient feels, does, or survives, or what a clinician decides.

They are separate questions. A result can be significant and trivial (Trial A), or non-significant and potentially important (a small trial whose CI includes a big benefit).

### The MCID

The **minimal clinically important difference (MCID)** is the smallest change that patients or clinicians would regard as worthwhile. Examples: about 5 mmHg for a blood-pressure drug in many settings; about 1.0 percentage point of HbA1c; about 10 mm on a 100 mm pain VAS. The MCID comes from clinical judgment, patient-anchored studies, or guideline thresholds, not from statistics. It should be **fixed before the trial** and used in the sample size calculation. At the end, compare the whole confidence interval with it, as in [Confidence intervals](#/test/confidence-intervals).

## Raw vs standardized effects

**Raw (unstandardized) effects** stay in natural units: a mean difference of −9 mmHg, a difference of 0.5% HbA1c. Prefer them whenever the unit is meaningful. A clinician understands "9 mmHg" immediately.

**Standardized effects** divide by the SD so that results from different scales can be compared or combined. The most common is **Cohen's d**:

d = (mean₁ − mean₂) / pooled SD

For Trial A d = −1.0/12 = −0.08; for Trial B d = −9/12 = −0.75.

Caveats:

- **Benchmarks are not universal.** Cohen's rough labels (0.2 small, 0.5 medium, 0.8 large) were meant as a last resort when nothing else was known. Whether d = 0.3 matters depends on the outcome: for mortality it can be enormous; for a symptom score it may be trivial.
- **d depends on the SD of the sample.** The same 5 mmHg change gives a different d in a narrow, homogeneous sample than in a broad one, so d is not a pure property of the treatment.
- **d is hard to explain to patients.** "0.5 SD" means little to a patient; "5 mmHg lower" does.
- **Report the CI for d** as well. In small samples use Hedges' g, which corrects the slight upward bias of d.

## Binary outcomes: RR, OR, RD, NNT

For yes/no outcomes the effect is a comparison of risks. Say 30 of 200 patients on a drug have a stroke, versus 45 of 200 on placebo (risks 15% vs 22.5%).

| Measure | Formula | Value |
|---|---|---|
| Risk difference (RD) | 0.15 − 0.225 | −0.075 (7.5 percentage points fewer strokes) |
| Risk ratio (RR) | 0.15 / 0.225 | 0.67 (a one-third relative reduction) |
| Odds ratio (OR) | (30/170) / (45/155) | 0.61 |
| Number needed to treat (NNT) | 1 / \|RD\| | 13.3 (treat about 14 patients to prevent one stroke) |

Relative measures (RR, OR) travel well across populations but hide the baseline risk: a 33% relative reduction of a 0.1% risk is tiny in absolute terms. Absolute measures (RD, NNT) show the real benefit to a patient. Report both. The details, including how to choose and how to compute CIs, are in [Effect measures](#/test/effect-measures).

## Always report the confidence interval

An effect size without a CI hides how precisely it was measured. In the stroke example the RD is −7.5 percentage points, but its 95% CI runs from −15.1 to +0.1 points: compatible with a large benefit and with essentially none (p = 0.055). The point estimate alone would mislead. Note that the NNT of 13.3 has no ordinary CI here: because the RD interval crosses 0, the NNT interval would run from a benefit of about 7 patients through infinity to a harm, so report the RD interval instead (see [Effect measures](#/test/effect-measures)). Give the estimate, its 95% CI, and a sentence on whether the CI includes effects that are clinically important.

## How to run it

These programs compute the effect size, 95% CI, and p-value for the two hypothetical trials from summary statistics (same p, different effect), then the RD, RR, OR, and NNT for the stroke example. Everything is deterministic, so R and Python give identical numbers.

### R

```r
# Two trials, SD 12 mmHg, same p-value but very different effects (from summary statistics)
trial <- function(diff, n, sd = 12) {              # n = patients per arm, difference in means
  se <- sd * sqrt(2 / n); df <- 2 * n - 2
  ci <- diff + c(-1, 1) * qt(0.975, df) * se
  c(n = n, diff = diff, lo = ci[1], hi = ci[2], d = diff / sd, p = 2 * pt(-abs(diff / se), df))
}
res <- rbind("Big trial, tiny effect" = trial(-1, 1360), "Small trial, big effect" = trial(-9, 17))
print(round(res, 3))

# Binary outcome: stroke in 30/200 on drug vs 45/200 on placebo
e <- c(30, 45); n <- c(200, 200); r <- e / n
rd <- r[1] - r[2]; rr <- r[1] / r[2]
or <- (e[1] / (n[1] - e[1])) / (e[2] / (n[2] - e[2]))
se_rd <- sqrt(sum(r * (1 - r) / n))
cat("Risks:", r, " RD =", round(rd, 3), " RR =", round(rr, 2), " OR =", round(or, 2), " NNT =", round(1 / abs(rd), 1), "\n")
cat("RD 95% CI:", round(rd + c(-1, 1) * 1.96 * se_rd, 3), " p =", round(prop.test(e, n, correct = FALSE)$p.value, 3), "\n")
```

**What the code does**

- `trial()` builds the standard error of a difference of two means, SD × √(2/n), then the t statistic, CI, Cohen's d, and two-sided p.
- Cohen's d is just the raw difference divided by the SD.
- The binary block computes RD, RR, and OR directly from the counts, and the NNT as 1/|RD|.
- The RD CI uses the Wald standard error, and `prop.test(correct = FALSE)` gives the uncorrected test for two proportions.

**Output**

```text
                           n diff      lo     hi      d     p
Big trial, tiny effect  1360   -1  -1.902 -0.098 -0.083 0.030
Small trial, big effect   17   -9 -17.384 -0.616 -0.750 0.036
Risks: 0.15 0.225  RD = -0.075  RR = 0.67  OR = 0.61  NNT = 13.3 
RD 95% CI: -0.151 0.001  p = 0.055 
```

### Python

```python
import numpy as np
from scipy import stats

# Two trials, SD 12 mmHg, same p-value but very different effects (from summary statistics)
def trial(diff, n, sd=12):                         # n = patients per arm, difference in means
    se, df = sd * np.sqrt(2 / n), 2 * n - 2
    h = stats.t.ppf(0.975, df) * se
    p = 2 * stats.t.sf(abs(diff / se), df)
    return [n, diff, diff - h, diff + h, diff / sd, p]
print("                          n  diff    lo     hi      d      p")
for name, args in (("Big trial, tiny effect", (-1, 1360)), ("Small trial, big effect", (-9, 17))):
    v = trial(*args)
    print(f"{name:<24}{v[0]:>5}{v[1]:>6.1f}{v[2]:>7.2f}{v[3]:>7.2f}{v[4]:>7.3f}{v[5]:>7.3f}")

# Binary outcome: stroke in 30/200 on drug vs 45/200 on placebo
e, n = np.array([30, 45]), np.array([200, 200]); r = e / n
rd, rr = r[0] - r[1], r[0] / r[1]
odds = (e[0] / (n[0] - e[0])) / (e[1] / (n[1] - e[1]))
se_rd = np.sqrt(np.sum(r * (1 - r) / n))
print("Risks:", r, " RD =", round(rd, 3), " RR =", round(rr, 2), " OR =", round(odds, 2), " NNT =", round(1 / abs(rd), 1))
pbar = e.sum() / n.sum(); z = rd / np.sqrt(pbar * (1 - pbar) * (1 / n[0] + 1 / n[1]))   # uncorrected z test
print("RD 95% CI:", np.round(rd + np.array([-1, 1]) * 1.96 * se_rd, 3), " p =", round(2 * stats.norm.sf(abs(z)), 3))
```

**What the code does**

- `trial()` returns the same quantities as the R version: n, difference, CI, d, and p, using the t distribution with 2n − 2 degrees of freedom.
- The loop prints each trial as one aligned row.
- The binary block computes RD, RR, OR, and NNT from the counts; the p-value comes from a pooled-variance z test, which equals R's uncorrected `prop.test`.

**Output**

```text
                          n  diff    lo     hi      d      p
Big trial, tiny effect   1360  -1.0  -1.90  -0.10 -0.083  0.030
Small trial, big effect    17  -9.0 -17.38  -0.62 -0.750  0.036
Risks: [0.15  0.225]  RD = -0.075  RR = 0.67  OR = 0.61  NNT = 13.3
RD 95% CI: [-0.151  0.001]  p = 0.055
```

## Interpreting the output

- **Both trials are "significant" with nearly the same p (0.030 and 0.036).** If you stopped at the p-value you would treat them as equivalent.
- **Trial A:** −1.0 mmHg (95% CI −1.9 to −0.1), d = −0.08. Even the most optimistic end of the CI is under 2 mmHg, well below a 5 mmHg MCID. Statistically significant, clinically trivial.
- **Trial B:** −9.0 mmHg (95% CI −17.4 to −0.6), d = −0.75. The estimate exceeds the MCID, but the interval is huge: the true effect might be anything from negligible to very large. Promising, but too imprecise to be confident.
- **Stroke example:** RD = −7.5 percentage points, RR = 0.67, OR = 0.61, NNT = 13.3. The OR looks more impressive than the RR because the outcome is not rare; do not read an OR as a risk ratio when events are common.
- **RD 95% CI −15.1 to +0.1 points, p = 0.055.** The upper limit is just above 0, so the result is borderline and the data are compatible with a large benefit.

**How to report it:**
> Stroke occurred in 15.0% on drug and 22.5% on placebo (risk difference −7.5 percentage points, 95% CI −15.1 to +0.1; RR 0.67; NNT 14, but its CI is not defined because the RD interval includes 0).

## Check your understanding

<details class="check"><summary>1. A 3,000-patient trial reports a 1.2 mmHg reduction in SBP, p < 0.001. The authors call it "highly clinically significant". Comment.</summary>

Statistical significance is not clinical significance. A p < 0.001 only says the effect is unlikely to be exactly zero. A 1.2 mmHg reduction is well below typical MCIDs for blood pressure (around 5 mmHg), though across a whole population even small shifts can matter for public health. The authors should compare the CI to the MCID, not use p to judge importance.
</details>

<details class="check"><summary>2. A small pilot reports a 10 mmHg reduction, p = 0.09, and the authors conclude "the drug is ineffective because the result was not significant". What is the more appropriate reading?</summary>

The point estimate is a large, clinically important effect, and p = 0.09 only reflects low precision (small n). The CI will be wide and almost certainly include both 0 and large benefit. The right conclusion is "inconclusive; a larger trial is warranted", not "ineffective".
</details>

<details class="check"><summary>3. In the stroke example the RR is 0.67. If the baseline stroke risk were only 1% (1.0% vs 0.67%), what would the RD and NNT be?</summary>

RD = 0.67% − 1.0% = −0.33 percentage points, so NNT = 1/0.0033 ≈ 303. The same relative reduction gives a much smaller absolute benefit when the baseline risk is low: you would treat about 300 people to prevent one stroke. Always report an absolute measure alongside the relative one.
</details>

<details class="check"><summary>4. A paper says "the effect was large (d = 0.9)" in a study where the sample SD of the outcome was unusually small because of strict eligibility criteria. Why might d be misleading?</summary>

Cohen's d divides the raw difference by the sample SD. A narrow, homogeneous sample has a small SD, which inflates d for the same raw difference. In a broader clinical population the same raw effect would produce a smaller d. Check the raw difference in natural units and compare it with the MCID instead of relying on benchmark labels.
</details>

## Key takeaways

- A p-value shows detectability, not importance. The same p can come from a trivial effect in a huge trial or a big effect in a tiny one.
- Report raw effects in natural units (mmHg) first; add standardized effects (Cohen's d) for comparison across scales, but do not trust benchmarks blindly.
- Define the MCID before the study and compare the whole CI with it.
- For binary outcomes report both relative (RR, OR) and absolute (RD, NNT) measures.
- Always give the 95% CI alongside the effect size.

## References

- Cohen J. The earth is round (p < .05). *American Psychologist.* 1994;49(12):997-1003.
- Cohen J. *Statistical Power Analysis for the Behavioral Sciences.* 2nd ed. Lawrence Erlbaum; 1988.
- Sullivan GM, Feinn R. Using effect size, or why the P value is not enough. *Journal of Graduate Medical Education.* 2012;4(3):279-282.
- Jaeschke R, Singer J, Guyatt GH. Measurement of health status: ascertaining the minimal clinically important difference. *Controlled Clinical Trials.* 1989;10(4):407-415.
