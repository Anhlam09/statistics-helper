# Equivalence (TOST) & non-inferiority tests

> Shows that two treatments differ by **less than a pre-specified margin** (equivalence) or that a new treatment is **not worse than a comparator by more than a margin** (non-inferiority), by reversing the usual null and alternative hypotheses.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Show two treatments are about the same, or a new one is not meaningfully worse | Continuous outcome (extends to other types) | Pre-specified clinical margin; independent groups; roughly normal data | Once-weekly vs once-daily diabetes drug on HbA1c |

## Learning objectives

- Explain why a non-significant t-test does **not** show that treatments are equivalent
- Define an equivalence or non-inferiority margin and why it must be set before the trial
- Run TOST in R and Python and compare the 90% CI with the margin
- Report equivalence and non-inferiority results with the CI and the margin

## What it is

An ordinary t-test with p > 0.05 does **not** show that two treatments are the same: it only says the data cannot rule out zero difference. Equivalence and non-inferiority tests turn this around. You first state, **before the trial**, a margin Δ: the largest difference that you would still call clinically unimportant. Then the burden of proof is on showing that the true difference lies inside that margin.

- **Equivalence (two-sided), TOST.** The *two one-sided tests* procedure tests H₀₁: μ_new − μ_ref ≤ −Δ and H₀₂: μ_new − μ_ref ≥ +Δ. Equivalence is concluded only if **both** are rejected at level α, so the TOST p-value is the **larger** of the two one-sided p-values.
- **Non-inferiority (one-sided).** For an outcome where higher is worse, test H₀: μ_new − μ_ref ≥ Δ against H₁: μ_new − μ_ref < Δ. Only the upper margin matters. If higher is better, flip the sign and test H₀: μ_new − μ_ref ≤ −Δ.
- **Confidence-interval approach (equivalent to TOST).** Equivalence at level α is concluded if the **100(1 − 2α)% CI**, i.e. the **90% CI** for α = 0.05, lies entirely inside (−Δ, +Δ). Non-inferiority at one-sided α = 0.025 is concluded if the upper limit of the usual **95% two-sided CI** is below Δ. The CI view is also the best way to report results, because it shows the margin and the uncertainty together.
- **Test statistic:** two t statistics (one per margin), each tested one-sided, here with Welch's unequal-variance degrees of freedom
- **Effect size:** the difference in means with its 90% CI (equivalence) or 95% CI (non-inferiority), always shown against the margin Δ. For bioequivalence trials, the analysis is typically done on the log scale and the margin is the ratio 0.80–1.25.

The choice and justification of the margin is a **clinical and regulatory** decision, not a statistical one. For a non-inferiority margin, the guidance asks that it be smaller than the effect the active comparator has been shown to have over placebo in earlier trials, and that enough of that effect is preserved (the M₁ and M₂ margins in FDA's terminology).

## The intuition

A normal t-test has a built-in bias: it starts by assuming "no difference" and asks whether the data **disprove** that. A big p-value then only means "we could not disprove it", which is also what happens when you have too little data. **Absence of evidence is not evidence of absence.**

Everyday analogy: a defendant is "not guilty" when the prosecution fails, which is different from "proven innocent". Equivalence testing flips the burden of proof: it assumes the treatments **differ by at least Δ** (the margin) and asks the data to **disprove that**. The margin Δ is the largest difference you would still call clinically unimportant, and it must be fixed before the trial.

Signal vs noise: the signal is how close to zero the estimated difference is. The noise is its uncertainty. Equivalence needs the whole confidence interval to fit **inside** the margin, so it needs a small difference *and* a narrow interval.

The plot shows the 30 HbA1c changes per arm (more negative is better). The two clouds overlap almost completely, and the group means (−1.22 and −1.20) are almost the same. Our estimated difference is −0.02 points, with 90% CI −0.27 to +0.22, which sits inside ±0.4.

<div class="viz" data-viz="dots" data-unit="% change in HbA1c" data-groups='{"New drug":[-0.56,-0.84,-1.16,-1.22,-0.46,-0.39,-0.81,-1.17,-1.77,-1.72,-2.02,-0.84,-0.64,-0.94,-2.9,-1.18,-1.61,-0.85,-1.98,-2.21,-0.47,-0.93,-0.53,-1.22,-0.88,-1.95,-1.17,-0.72,-2.23,-1.16],"Comparator":[-1.97,-1.36,-1.46,-0.52,-1.45,-0.97,-0.99,-1.45,-1.49,-0.94,-1.66,-1.18,-0.91,-0.09,-1.84,-1.11,-1.0,-0.81,-0.83,-0.11,-1.48,-1.6,-1.54,-1.75,-1.55,-1.97,-0.69,-1.31,-0.63,-1.2]}' data-caption="Almost complete overlap and nearly identical means. The question is whether the 90% CI for the difference stays within the margin of plus or minus 0.4."></div>

| Result | Difference | Interval | Margin | Conclusion |
|---|---|---|---|---|
| Equivalence (TOST) | −0.02 | 90% CI −0.27 to +0.22 | ±0.4 | Inside, so equivalent |
| Non-inferiority | −0.02 | 95% CI upper limit 0.27 | 0.3 | Below, so non-inferior |
| Ordinary t-test | −0.02 | p = 0.88 | none | Says nothing about equivalence |

Notice the last row. A p-value of 0.88 sounds like "the same", but the t-test would give a similar p-value for a much noisier study whose CI ran from −1.5 to +1.5. Only the CI against the margin tells you whether "the same" is justified.

## When to use it

Use TOST or a non-inferiority test when **all** of these are true:

- The research question is "**no meaningfully worse**" or "**about the same**", not "different"
- A clinically justified, **pre-specified margin** exists (set in the protocol and statistical analysis plan)
- The outcome is continuous and the groups are independent (the same logic extends to proportions, hazard ratios and paired or crossover designs)

**Clinical scenario.** A new once-weekly drug is compared with an established once-daily comparator for type 2 diabetes in 30 patients per arm. The primary endpoint is the change in HbA1c (%) from baseline to week 26 (more negative is better). The sponsor wants to show that the new drug is **equivalent** within ±0.4 percentage points and, as a separate (non-inferiority) claim, that it is **not worse than the comparator by 0.3 percentage points or more**. *Does the evidence support these claims?*

**Don't use it when:**

- You only want to detect a **difference** (superiority) → use the [independent t-test](#/test/independent-t) or [ANCOVA](#/test/ancova)
- The margin was chosen after seeing the data: the test is then not interpretable, and regulators will not accept it
- The outcome is a proportion or a time-to-event; use the corresponding equivalence/non-inferiority versions (for example a CI for a risk difference against the margin, or a [Cox model](#/test/cox-ph) hazard-ratio CI against the margin)
- The data are clearly non-normal with small n → use a rank-based or bootstrap approach, e.g. a Hodges–Lehmann CI from the [Mann–Whitney test](#/test/mann-whitney) compared with the margin

## Assumptions

1. **Independent observations** within and between groups (for paired or crossover data use the paired TOST). Check the design.
2. **Approximately normal outcomes or enough observations** for the t-test to be reliable (about 30 per group is usually adequate). Check dot plots or Q–Q plots; a [Shapiro–Wilk test](#/test/shapiro-wilk) can support this but should not be used to choose between methods.
3. **Variances:** the Welch version used here does not assume equal variances; the pooled version does (check with [Levene's test](#/test/levene)).
4. **A margin Δ justified clinically and fixed in advance**, including for the non-inferiority claim the historical evidence on the comparator's effect (ICH E10). Review the protocol and the statistical analysis plan.
5. **Trial conduct and analysis population:** poor adherence, missing data and cross-over between arms bias results towards "no difference", which favours equivalence. For equivalence and non-inferiority trials, analyse both the intention-to-treat and per-protocol populations and check that they agree (ICH E9).
6. **Assay sensitivity:** the trial must be capable of detecting a difference if one existed. This rests on the comparator having a well-established effect and on the study being run like the earlier trials.

## Hypotheses

For **equivalence** with margin Δ = 0.4, where δ = μ_new − μ_ref:

- **H₀:** δ ≤ −0.4 **or** δ ≥ +0.4 (the treatments differ by at least the margin)
- **H₁:** −0.4 < δ < +0.4 (the treatments are equivalent)

For **non-inferiority** with margin Δ = 0.3 (higher HbA1c = worse):

- **H₀:** δ ≥ 0.3 (new is inferior by at least the margin)
- **H₁:** δ < 0.3 (new is non-inferior)

## Example data

Change in HbA1c (%) from baseline to week 26 (lower is better), 30 patients per arm. Both arms improved by roughly 1.2 percentage points on average.

| Arm | Values (first 10 of 30 shown) | Mean |
|---|---|---|
| New drug | −0.56, −0.84, −1.16, −1.22, −0.46, −0.39, −0.81, −1.17, −1.77, −1.72, … | −1.218 |
| Comparator | −1.97, −1.36, −1.46, −0.52, −1.45, −0.97, −0.99, −1.45, −1.49, −0.94, … | −1.195 |

The full 60 values appear in the code blocks below. The data were simulated once in R (seed 2024) and are identical in all three languages.

## How to run it

### R

```r
# Change in HbA1c (%) from baseline to week 26, n = 30 per arm (lower = better)
new <- c(-0.56, -0.84, -1.16, -1.22, -0.46, -0.39, -0.81, -1.17, -1.77, -1.72,
         -2.02, -0.84, -0.64, -0.94, -2.90, -1.18, -1.61, -0.85, -1.98, -2.21,
         -0.47, -0.93, -0.53, -1.22, -0.88, -1.95, -1.17, -0.72, -2.23, -1.16)
ref <- c(-1.97, -1.36, -1.46, -0.52, -1.45, -0.97, -0.99, -1.45, -1.49, -0.94,
         -1.66, -1.18, -0.91, -0.09, -1.84, -1.11, -1.00, -0.81, -0.83, -0.11,
         -1.48, -1.60, -1.54, -1.75, -1.55, -1.97, -0.69, -1.31, -0.63, -1.20)
delta <- 0.4   # equivalence margin (+/- 0.4 %); non-inferiority margin below is 0.3

# 1. CI approach: 90% CI for the difference (new - ref), Welch
t.test(new, ref, conf.level = 0.90)$conf.int

# 2. TOST: two one-sided Welch tests against -delta and +delta
p_low <- t.test(new, ref, mu = -delta, alternative = "greater")$p.value  # H0: diff <= -delta
p_upp <- t.test(new, ref, mu =  delta, alternative = "less")$p.value     # H0: diff >= +delta
c(p_lower = p_low, p_upper = p_upp, p_TOST = max(p_low, p_upp))

# 3. Non-inferiority (higher HbA1c = worse), margin 0.3: one-sided 97.5% upper bound
t.test(new, ref, alternative = "less", mu = 0.3, conf.level = 0.975)

# 4. Usual two-sided Welch test of "no difference" (NOT a test of equivalence)
t.test(new, ref)$p.value
```

**What the code does**

- `t.test(new, ref, conf.level = 0.90)$conf.int` gives the **90% CI** for new − ref. R's `t.test` uses Welch's test by default (`var.equal = FALSE`). The CI-inclusion rule says: equivalence if the CI lies within (−0.4, +0.4).
- The two `t.test` calls with `mu = −delta`, `alternative = "greater"` and `mu = +delta`, `alternative = "less"` are the **two one-sided tests**. `mu` is the null difference, so each call tests one margin. The TOST p-value is `max()` of the two. Both tests use α = 0.05; you do **not** split the α.
- Step 3 is the non-inferiority test: `alternative = "less", mu = 0.3` tests H₀: δ ≥ 0.3. Using `conf.level = 0.975` gives the one-sided 97.5% upper bound (the CI line says `-Inf 0.2695`). The usual regulatory convention for non-inferiority is a one-sided α of 0.025, equivalent to a two-sided 95% CI.
- The last line shows the ordinary two-sided Welch p-value, which is large and says nothing about equivalence.
- The **TOSTER** package (`t_TOST`) reproduces this with one call and a plot, if you prefer.

**Output**

```text
[1] -0.2659536  0.2212869
attr(,"conf.level")
[1] 0.9
    p_lower     p_upper      p_TOST 
0.006077746 0.002677393 0.006077746 

	Welch Two Sample t-test

data:  new and ref
t = -2.2138, df = 54.691, p-value = 0.01552
alternative hypothesis: true difference in means is less than 0.3
97.5 percent confidence interval:
      -Inf 0.2694961
sample estimates:
mean of x mean of y 
-1.217667 -1.195333 

[1] 0.8786581
```

### Python

```python
import numpy as np
from scipy import stats
from statsmodels.stats.weightstats import ttost_ind

# Change in HbA1c (%) from baseline to week 26, n = 30 per arm (lower = better)
new = np.array([-0.56, -0.84, -1.16, -1.22, -0.46, -0.39, -0.81, -1.17, -1.77, -1.72,
                -2.02, -0.84, -0.64, -0.94, -2.90, -1.18, -1.61, -0.85, -1.98, -2.21,
                -0.47, -0.93, -0.53, -1.22, -0.88, -1.95, -1.17, -0.72, -2.23, -1.16])
ref = np.array([-1.97, -1.36, -1.46, -0.52, -1.45, -0.97, -0.99, -1.45, -1.49, -0.94,
                -1.66, -1.18, -0.91, -0.09, -1.84, -1.11, -1.00, -0.81, -0.83, -0.11,
                -1.48, -1.60, -1.54, -1.75, -1.55, -1.97, -0.69, -1.31, -0.63, -1.20])
delta = 0.4   # equivalence margin (+/- 0.4 %); non-inferiority margin below is 0.3

# 1. CI approach: 90% CI for the difference (new - ref), Welch
res = stats.ttest_ind(new, ref, equal_var=False)
ci = res.confidence_interval(confidence_level=0.90)
print(f"Difference = {new.mean() - ref.mean():.3f}, 90% CI ({ci.low:.3f}, {ci.high:.3f})")

# 2. TOST with unequal variances (Welch)
p, lower, upper = ttost_ind(new, ref, -delta, delta, usevar="unequal")
print(f"Lower test: t = {lower[0]:.3f}, p = {lower[1]:.4f}, df = {lower[2]:.1f}")
print(f"Upper test: t = {upper[0]:.3f}, p = {upper[1]:.4f}, df = {upper[2]:.1f}")
print(f"TOST p = {p:.4f}")

# 3. Non-inferiority, margin 0.3 (higher = worse): one-sided 97.5% upper bound and p
ni = stats.ttest_ind(new, ref, equal_var=False, alternative="less")
ci1 = ni.confidence_interval(confidence_level=0.975)
p_ni = stats.ttest_ind(new - 0.3, ref, equal_var=False, alternative="less").pvalue  # H0: diff >= 0.3
print(f"Upper 97.5% bound = {ci1.high:.3f} (margin 0.3), one-sided p = {p_ni:.4f}")

# 4. Usual two-sided test of "no difference" (NOT a test of equivalence)
print(f"Two-sided Welch p = {res.pvalue:.4f}")
```

**What the code does**

- `ttost_ind(x1, x2, low, upp, usevar="unequal")` from `statsmodels.stats.weightstats` runs TOST for two independent samples. It returns the overall TOST p-value (the maximum of the two) and the `(t, p, df)` triple for each one-sided test: the first tests against the lower margin (H₀: δ ≤ low), the second against the upper margin (H₀: δ ≥ upp). `usevar="unequal"` is Welch's version, matching R. **The default is `"pooled"`**, which assumes equal variances and would give slightly different degrees of freedom.
- The margins are given as `low = −delta` and `upp = +delta` on the scale of the **difference** x1 − x2. Check the sign convention: swapping the arguments reverses the difference.
- `stats.ttest_ind(..., equal_var=False).confidence_interval(0.90)` gives the 90% CI (SciPy ≥ 1.11). `equal_var=False` is Welch's test; the default is the pooled Student test.
- `alternative="less"` with `confidence_level=0.975` yields the one-sided upper bound for the non-inferiority comparison (margin 0.3). For paired or crossover data use `ttost_paired`.

**Output**

```text
Difference = -0.022, 90% CI (-0.266, 0.221)
Lower test: t = 2.594, p = 0.0061, df = 54.7
Upper test: t = -2.901, p = 0.0027, df = 54.7
TOST p = 0.0061
Upper 97.5% bound = 0.269 (margin 0.3), one-sided p = 0.0155
Two-sided Welch p = 0.8787
```

### SAS

```sas
/* HbA1c change (%) from baseline, 30 patients per arm.
   Arms are named so that PROC TTEST computes the difference New - Ref
   (it subtracts the second class level from the first, in alphabetical order). */
data hba1c;
  input arm $ chg @@;
  datalines;
  New -0.56  New -0.84  New -1.16  New -1.22  New -0.46
  New -0.39  New -0.81  New -1.17  New -1.77  New -1.72
  New -2.02  New -0.84  New -0.64  New -0.94  New -2.90
  New -1.18  New -1.61  New -0.85  New -1.98  New -2.21
  New -0.47  New -0.93  New -0.53  New -1.22  New -0.88
  New -1.95  New -1.17  New -0.72  New -2.23  New -1.16
  Ref -1.97  Ref -1.36  Ref -1.46  Ref -0.52  Ref -1.45
  Ref -0.97  Ref -0.99  Ref -1.45  Ref -1.49  Ref -0.94
  Ref -1.66  Ref -1.18  Ref -0.91  Ref -0.09  Ref -1.84
  Ref -1.11  Ref -1.00  Ref -0.81  Ref -0.83  Ref -0.11
  Ref -1.48  Ref -1.60  Ref -1.54  Ref -1.75  Ref -1.55
  Ref -1.97  Ref -0.69  Ref -1.31  Ref -0.63  Ref -1.20
;
run;

/* 1. TOST for equivalence, margins -0.4 to +0.4, alpha = 0.05 per one-sided test */
proc ttest data=hba1c tost(-0.4, 0.4) alpha=0.05;
  class arm;
  var chg;
run;

/* 2. 90% two-sided CI for the difference (compare it with the margins) */
proc ttest data=hba1c alpha=0.10 plots=none;
  class arm;
  var chg;
run;

/* 3. Non-inferiority, margin 0.3 (higher HbA1c = worse):
      H0: New - Ref >= 0.3  vs  H1: New - Ref < 0.3, one-sided alpha = 0.025 */
proc ttest data=hba1c sides=L h0=0.3 alpha=0.025 plots=none;
  class arm;
  var chg;
run;
```

**What the code does**

- The data are in **long format** with `arm` and `chg`. PROC TTEST forms the difference as the **first class level minus the second** (alphabetical order), so naming the arms `New` and `Ref` gives New − Ref, the same direction as R and Python.
- `tost(-0.4, 0.4)` on the PROC TTEST statement (SAS/STAT 14.1 or later; on older releases use the CI approach below) requests the **two one-sided tests** for equivalence with the given lower and upper margins. SAS reports the one-sided t tests for the lower and upper margins and an overall TOST p-value, for both the pooled and Satterthwaite (Welch) methods. Read the **Satterthwaite** rows to compare with R and Python.
- The second `proc ttest ... alpha=0.10` prints the 90% two-sided CI for the difference. This is the CI-inclusion approach: check that both limits lie within (−0.4, +0.4).
- The third step performs the **non-inferiority** test: `h0=0.3` sets the null difference and `sides=L` makes the alternative "difference less than 0.3". With `alpha=0.025`, the printed one-sided CI has an upper limit that you compare with 0.3. A one-sided p below 0.025 means non-inferiority.
- Check that `class arm` ordering matches your intended contrast; in the output, the line "Diff (1-2)" should be New − Ref.

## Interpreting the output

- **Difference (new − ref) = −0.022 percentage points** (−1.218 vs −1.195). The point estimate is essentially zero.
- **90% CI: −0.266 to +0.221.** Both limits lie inside the equivalence margin of ±0.4, so equivalence is shown. The data are compatible with a true difference of about a quarter of a percentage point in either direction, but not 0.4.
- **TOST: lower test t = 2.59, p = 0.0061; upper test t = −2.90, p = 0.0027; overall TOST p = 0.0061 (the larger).** Both one-sided nulls are rejected at α = 0.05. R, SciPy/statsmodels agree.
- **Non-inferiority (margin 0.3):** the upper limit of the 95% two-sided CI is 0.269 (R's one-sided 97.5% bound, 0.2695), below 0.3, and the one-sided p-value is 0.0155 (< 0.025), so non-inferiority is shown. Note that the margin 0.3 is tighter than 0.4; with this CI width, a margin of 0.25 would **not** have been met. The margin must be fixed in advance, not after seeing the CI.
- **Ordinary Welch p = 0.88.** This only says there is no evidence of a difference. It does not show equivalence: with an imprecise estimate, a non-significant p-value would also arise if the groups differed by a lot. The equivalence conclusion comes from the CI against the margin.

**How to report it:**
> The mean change in HbA1c was −1.22% with the new drug and −1.20% with the comparator (difference −0.02 percentage points, 90% CI −0.27 to 0.22). Because the 90% CI lay within the pre-specified equivalence margin of ±0.4 (TOST p = 0.006), the treatments were considered equivalent. The upper limit of the 95% CI (0.27) was below the non-inferiority margin of 0.3 (one-sided p = 0.016).

⚠️ **Important for trials:** equivalence and non-inferiority conclusions depend on a margin that is clinically justified and pre-specified in the protocol, on assay sensitivity, and on high-quality trial conduct. Report both ITT and per-protocol results. Follow ICH E9/E10 and the relevant FDA guidance; a sample size must be calculated for the equivalence or non-inferiority hypothesis, not for superiority.

## Check your understanding

<details class="check"><summary>1. A trial of 12 patients per arm gives p = 0.60 on a t-test and a 95% CI for the difference of −1.1 to +0.8 mmHg (margin ±0.5). Are the treatments equivalent?</summary>

No. The CI extends well beyond ±0.5 on both sides, so the data cannot rule out a clinically important difference. A non-significant t-test with a wide CI shows the study was uninformative, not that the treatments are the same.
</details>

<details class="check"><summary>2. Why does TOST at α = 0.05 use a 90% CI instead of a 95% CI?</summary>

TOST runs two one-sided tests, each at the full α = 0.05. Together these correspond to the 100(1 − 2α)% = 90% two-sided CI lying inside (−Δ, +Δ). Using 95% would be too conservative, since it implies each one-sided test is run at 2.5%.
</details>

<details class="check"><summary>3. After seeing a 90% CI of −0.27 to +0.22, the team decides the margin should have been ±0.3 and writes that into the report. Is this acceptable?</summary>

No. Choosing or changing the margin after seeing the data makes the conclusion uninterpretable and regulators will not accept it. The margin must be justified clinically and fixed in the protocol and statistical analysis plan.
</details>

<details class="check"><summary>4. What if the true difference were still −0.02 but the standard deviation were twice as large? Would equivalence within ±0.4 still be shown with 30 patients per arm?</summary>

Probably not. Doubling the SD doubles the standard error, so the 90% CI would stretch to roughly −0.5 to +0.5, going outside ±0.4. Equivalence needs a small difference **and** enough precision, which is why equivalence trials often need larger samples than superiority trials.
</details>

## Common mistakes

- **Concluding equivalence from a non-significant t-test.** "No significant difference" is not "equivalent".
- **Using a 95% CI for TOST at α = 0.05.** The correct interval for two one-sided 5% tests is the **90% CI**. Using 95% is too conservative.
- **Choosing the margin after seeing the data.** The margin must be pre-specified and justified clinically; changing it afterwards invalidates the conclusion.
- **Splitting α.** In TOST each one-sided test is done at the full α (0.05); the overall type I error remains ≤ 0.05.
- **Reading the contrast direction backwards.** In SAS the contrast is the first class minus the second (alphabetical); in R `t.test(x, y)` is x − y; in statsmodels it is x1 − x2. For non-inferiority, be sure which direction is "worse".
- **Using only the intention-to-treat population.** ITT tends to bias toward equivalence (dilution), so analyse per-protocol as well.
- **Sizing the trial for superiority.** Equivalence and non-inferiority designs need their own sample size calculations, often larger.

## Related tests

- [Independent two-sample t-test](#/test/independent-t): the superiority test
- [Paired t-test](#/test/paired-t): the paired data version (use `ttost_paired` in Python)
- [Mann–Whitney U test](#/test/mann-whitney): a rank-based Hodges–Lehmann CI to compare with the margin
- [Two-proportion z-test](#/test/two-proportion-z): for binary endpoints, compare the CI of the risk difference with the margin
- [Bland–Altman analysis](#/test/bland-altman): agreement between measurement methods, a different question

## Key takeaways

- A non-significant t-test (p > 0.05) does **not** show equivalence; it may simply reflect low power.
- TOST flips the hypotheses: H₀ is "differs by at least Δ", and both one-sided tests must reject.
- Equivalence is concluded when the **90% CI** lies entirely within (−Δ, +Δ) for α = 0.05.
- The margin Δ is a clinical and regulatory decision, fixed **before** seeing the data.
- Equivalence and non-inferiority trials need their own sample-size calculation and a check of both ITT and per-protocol results.

## References

- ICH E9: Statistical Principles for Clinical Trials (1998); ICH E10: Choice of Control Group and Related Issues in Clinical Trials (2000).
- US FDA. *Non-Inferiority Clinical Trials to Establish Effectiveness: Guidance for Industry.* November 2016.
- Schuirmann DJ. A comparison of the two one-sided tests procedure and the power approach for assessing the equivalence of average bioavailability. *J Pharmacokinet Biopharm.* 1987;15:657–680.
- Documentation: statsmodels `weightstats.ttost_ind`; R `?t.test`; SAS/STAT User's Guide, The TTEST Procedure — TOST option.
