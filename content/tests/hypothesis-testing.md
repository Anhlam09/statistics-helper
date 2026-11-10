# Hypothesis testing and p-values

> A hypothesis test asks how surprising your data would be if there were truly no effect; the p-value is the answer, and it is one of the most misread numbers in medicine.

## Learning objectives

- State a null and an alternative hypothesis for a clinical question
- Explain a test statistic as signal divided by noise, and a p-value as a tail area under H₀
- Recognize what a p-value is **not** (not P(H₀ is true), not the effect size, not proof of "no effect")
- Run a simulation that shows p-values are uniform when H₀ is true and that about 5% fall below 0.05

## The intuition

Imagine a court case. The defendant is presumed innocent (the **null hypothesis**, H₀) until the evidence is strong enough to doubt it. A hypothesis test does the same: it starts by assuming "nothing is going on", then asks how well the data fit that assumption.

Our running example. A trial randomizes 30 patients with hypertension to a new drug or placebo (15 per arm). After 8 weeks we measure systolic blood pressure (SBP, mmHg). The drug group averaged 130.3 mmHg and the placebo group 134.9 mmHg, a difference of −4.6 mmHg. Is that a real effect, or could chance alone produce a gap like that between two groups of 15?

<div class="viz" data-viz="dots" data-unit="mmHg" data-groups='{"Drug":[133,130,144,126,125,125,133,128,137,134,133,131,124,130,121],"Placebo":[133,140,122,124,154,129,156,147,139,128,130,136,122,137,126]}' data-caption="The drug group sits a little lower, but the two clouds overlap a lot. Patients differ from each other far more than the arms differ."></div>

Look at the picture: the means differ, but individual patients vary by 20 mmHg or more within each arm. That within-group spread is the **noise**. The gap between the means is the **signal**. Every test statistic is a version of

**test statistic = signal / noise**, for example t = (difference in means) / (standard error of that difference).

Here t = −4.6 / 3.18 = −1.45. The signal is only about 1.4 times the noise. Now the key move: if H₀ were true (no drug effect), repeated trials would give t values that follow the bell-shaped t curve with 28 degrees of freedom, centered on 0. The p-value is the **tail area**: the probability, assuming H₀, of getting a statistic at least as extreme as ours.

<div class="viz" data-viz="t" data-stat="-1.45" data-df="28"></div>

Notice that the shaded tails (both sides, because the two-sided test counts a surprisingly *large* positive t as equally interesting) cover about 16% of the curve. A result like ours is quite ordinary in a world with no effect. For contrast, a statistic of 1.96 on a z scale marks the edge of the classic 5% region:

<div class="viz" data-viz="z" data-stat="1.96"></div>

The shaded area here is almost exactly 0.05. That is where the famous threshold comes from: a statistic beyond ±1.96 (for large samples) puts you in the most extreme 5% of the null world.

## Hypotheses

- **H₀ (null):** the true mean SBP is the same on drug and placebo, μ_drug − μ_placebo = 0.
- **H₁ (alternative):** the true means differ, μ_drug − μ_placebo ≠ 0 (two-sided).

You choose α (usually 0.05) **before** the trial. If p ≤ α you "reject H₀"; otherwise you "fail to reject H₀".

## The p-value, defined carefully

> The p-value is the probability of getting a test statistic at least as extreme as the one observed, **if H₀ and all the other model assumptions are true**.

In our trial, p = 0.159. In words: if the drug did nothing, about 16 trials in 100 would show a gap at least this large (in either direction) just from random patient-to-patient variation.

### What p is NOT

- **Not the probability that H₀ is true.** p = 0.159 does *not* mean "a 16% chance the drug does nothing". The p-value is calculated *assuming* H₀; it cannot tell you the probability of H₀. That needs prior information (a Bayesian calculation).
- **Not the size or importance of the effect.** A tiny p can come from a trivial effect in a huge trial; a large p can hide a big effect in a tiny trial (see [Effect sizes and clinical significance](#/test/effect-sizes)).
- **Not proof of no effect when p > 0.05.** In the simulated world below where a drug truly lowers SBP by 7 mmHg, a trial this size often gives p > 0.05, just as ours did (p = 0.16). "Not significant" often just means "not enough data".
- **Not a measure of reproducibility.** If p = 0.03, the chance that a replication also gives p < 0.05 is nowhere near 97%.

### The ASA statement

In 2016 the American Statistical Association took the unusual step of issuing a formal statement on p-values (Wasserstein & Lazar). Its six principles, paraphrased: (1) p-values can indicate how incompatible the data are with a specified model; (2) they do not measure the probability that the hypothesis is true, or that the data arose by chance alone; (3) conclusions and policy should not rest only on whether p crosses a threshold; (4) proper inference requires full reporting and transparency (no hidden multiple looks); (5) a p-value does not measure effect size or importance; (6) by itself, a p-value is not a good measure of evidence about a model or hypothesis. The practical lesson: report the estimate and its confidence interval, and treat p as one piece of evidence.

### One-sided vs two-sided

A **two-sided** test counts surprising results in *both* directions (drug much better or much worse). A **one-sided** test counts only one direction. For our data, the two-sided p is 0.159 and the one-sided p (drug lower) is 0.079, exactly half. Because a one-sided test makes it easier to reach p < 0.05, regulators generally expect two-sided tests, or a one-sided test at α = 0.025, and the direction must be fixed **before** seeing the data. Choosing one-sided after peeking at the results is a form of cheating.

### "Fail to reject" is not "accept"

A court says "not guilty", never "innocent". Likewise, a non-significant test says only that the data did not provide enough evidence against H₀. It does not show H₀ is true. The confidence interval tells you more: ours runs from −11.1 to +1.9 mmHg, which includes no effect and also includes a clinically useful 10 mmHg drop. If you really need to show two treatments are the same, design an [equivalence study](#/test/tost-equivalence).

## How to run it

Both programs analyse the same 30 patients, then simulate 5,000 trials in two worlds: one where H₀ is true and one where the drug truly lowers SBP by 7 mmHg (SD 12, 15 per arm). The simulated data differ between R and Python because the random number generators differ, but the conclusions agree.

### R

```r
# Worked example: SBP (mmHg) after 8 weeks, 15 patients per arm
drug    <- c(133, 130, 144, 126, 125, 125, 133, 128, 137, 134, 133, 131, 124, 130, 121)
placebo <- c(133, 140, 122, 124, 154, 129, 156, 147, 139, 128, 130, 136, 122, 137, 126)
t.test(drug, placebo, var.equal = TRUE)          # pooled-variance t-test, two-sided
t.test(drug, placebo, var.equal = TRUE, alternative = "less")$p.value  # one-sided

# Simulate 5,000 trials in which H0 is TRUE (both arms mean 140, SD 12)
set.seed(2024)
nsim <- 5000; n <- 15; sd <- 12
sim_p <- function(delta) replicate(nsim, {
  x <- rnorm(n, 140 - delta, sd); y <- rnorm(n, 140, sd)
  t.test(x, y, var.equal = TRUE)$p.value
})
p0 <- sim_p(0)                                   # no true effect
cat("H0 true: share with p < 0.05 =", mean(p0 < 0.05), "\n")
print(table(cut(p0, seq(0, 1, 0.2), include.lowest = TRUE)) / nsim)  # flat = uniform

p1 <- sim_p(7)                                   # true effect = 7 mmHg
cat("H1 true (7 mmHg): share with p < 0.05 =", mean(p1 < 0.05), "\n")
cat("H1 true: median p =", round(median(p1), 3), "\n")
```

**What the code does**

- `t.test(..., var.equal = TRUE)` runs the classic pooled t-test; `alternative = "less"` asks the one-sided question "is the drug lower?".
- `sim_p()` runs a whole fake trial 5,000 times and keeps only the p-value from each.
- With `delta = 0` every trial is a null trial, so any p < 0.05 is by definition a false positive.
- `table(cut(...))` counts p-values in five equal bins; under H₀ each bin should hold about 20%.
- With `delta = 7`, H₁ is true, so p-values pile up near 0 but many still exceed 0.05.

**Output**

```text
	Two Sample t-test

data:  drug and placebo
t = -1.4475, df = 28, p-value = 0.1588
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -11.109428   1.909428
sample estimates:
mean of x mean of y 
 130.2667  134.8667 

[1] 0.07942496
H0 true: share with p < 0.05 = 0.0496 

  [0,0.2] (0.2,0.4] (0.4,0.6] (0.6,0.8]   (0.8,1] 
   0.2002    0.2018    0.1938    0.1968    0.2074 
H1 true (7 mmHg): share with p < 0.05 = 0.3322 
H1 true: median p = 0.12 
```

### Python

```python
import numpy as np
from scipy import stats

# Worked example: SBP (mmHg) after 8 weeks, 15 patients per arm
drug    = np.array([133, 130, 144, 126, 125, 125, 133, 128, 137, 134, 133, 131, 124, 130, 121])
placebo = np.array([133, 140, 122, 124, 154, 129, 156, 147, 139, 128, 130, 136, 122, 137, 126])
r = stats.ttest_ind(drug, placebo)               # pooled-variance, two-sided
print(f"t = {r.statistic:.3f}, df = {r.df:.0f}, two-sided p = {r.pvalue:.4f}")
print(f"one-sided p (drug lower) = {stats.ttest_ind(drug, placebo, alternative='less').pvalue:.4f}")

# Simulate 5,000 trials in which H0 is TRUE (both arms mean 140, SD 12)
rng = np.random.default_rng(2024)
nsim, n, sd = 5000, 15, 12
def sim_p(delta):
    x = rng.normal(140 - delta, sd, (nsim, n))   # one row per simulated trial
    y = rng.normal(140, sd, (nsim, n))
    return stats.ttest_ind(x, y, axis=1).pvalue  # t-test on every row at once

p0 = sim_p(0)                                    # no true effect
print("H0 true: share with p < 0.05 =", np.mean(p0 < 0.05))
print("Share of p-values per 0.2 bin:", np.histogram(p0, bins=5, range=(0, 1))[0] / nsim)

p1 = sim_p(7)                                    # true effect = 7 mmHg
print("H1 true (7 mmHg): share with p < 0.05 =", np.mean(p1 < 0.05))
print("H1 true: median p =", round(float(np.median(p1)), 3))
```

**What the code does**

- `stats.ttest_ind` is the independent-samples t-test (pooled variance by default); `alternative='less'` gives the one-sided p.
- Passing 2-D arrays with `axis=1` runs 5,000 t-tests in one vectorized call, which is why it is fast.
- `rng = np.random.default_rng(2024)` fixes the seed so your numbers match the output below.
- The histogram counts should be near 0.20 in every bin if p-values are uniform under H₀.

**Output**

```text
t = -1.448, df = 28, two-sided p = 0.1588
one-sided p (drug lower) = 0.0794
H0 true: share with p < 0.05 = 0.0542
Share of p-values per 0.2 bin: [0.201  0.1978 0.193  0.203  0.2052]
H1 true (7 mmHg): share with p < 0.05 = 0.3428
H1 true: median p = 0.115
```

## Interpreting the output

- **t(28) = −1.45, p = 0.159.** The observed −4.6 mmHg gap is well within what chance produces. We fail to reject H₀.
- **One-sided p = 0.079** is exactly half the two-sided p, because all of the tail area is on one side.
- **When H₀ is true, about 5% of trials (0.050 in R, 0.054 in Python) give p < 0.05.** That is what "α = 0.05" means: a 5% false-positive rate per test when nothing is going on.
- **The five bins hold about 20% each.** Under H₀ the p-value is **uniform** between 0 and 1: p = 0.9 is just as likely as p = 0.02. That is why a single p = 0.04 is weak evidence.
- **When the drug truly works (7 mmHg), only 33% to 34% of trials reach p < 0.05.** The trial is badly underpowered. The median p is 0.12, so the *typical* trial of this size would miss a real, clinically relevant effect. Our non-significant result is exactly what we would expect to see from a drug that works. See [Type I and II errors, and power](#/test/errors-and-power).

**How to report it:**
> Mean SBP was 4.6 mmHg lower on drug than placebo (95% CI −11.1 to 1.9; t(28) = −1.45, two-sided p = 0.16). The study could not exclude either no effect or a clinically useful reduction.

## Check your understanding

<details class="check"><summary>1. A trial reports p = 0.03. Does that mean there is a 3% chance the drug has no effect?</summary>

No. The p-value is the probability of data this extreme **given** H₀ is true, not the probability that H₀ is true given the data. Turning it around needs prior probability of an effect, which a p-value does not use. A drug with a very implausible mechanism and p = 0.03 may still be probably ineffective.
</details>

<details class="check"><summary>2. Our trial gave p = 0.16. A colleague concludes "the drug does not lower blood pressure". Is that justified?</summary>

No. Failing to reject H₀ is not the same as showing H₀ is true. The 95% CI (−11.1 to +1.9 mmHg) is compatible with no effect and with a 10 mmHg drop. In the simulation where the drug truly lowers SBP by 7 mmHg, a trial this size misses significance about two times in three.
</details>

<details class="check"><summary>3. In a world where H₀ is true, you run 20 independent trials. Roughly how many will give p < 0.05, and what does the distribution of their p-values look like?</summary>

About 1 of 20 (5%) will be "significant" by chance. The p-values from null trials are uniformly spread between 0 and 1, so each tenth of the range is equally likely. This is why running many tests and reporting only the significant ones is misleading.
</details>

<details class="check"><summary>4. Your two-sided p is 0.08. You then decide the drug "obviously can only lower BP" and switch to a one-sided test, getting p = 0.04. Is that acceptable?</summary>

No. The direction of the alternative must be fixed before the data are seen, and justified by the design (for example, a harm in the other direction is impossible or irrelevant). Switching after seeing the data halves the p-value to cross the threshold and inflates the false-positive rate. Report the pre-specified analysis.
</details>

## Key takeaways

- A test statistic is **signal / noise**; the p-value is the tail area of its distribution **assuming H₀**.
- p is **not** P(H₀ true), not the effect size, and p > 0.05 is **not** evidence of no effect.
- Under H₀, p-values are uniform, so about 5% of true-null tests give p < 0.05; with low power, many real effects also give p > 0.05.
- Decide one- vs two-sided, and α, before the data. Say "fail to reject", never "accept H₀".
- Always report the estimate with its 95% CI; follow the ASA's advice not to rest conclusions on a threshold alone.

## References

- Wasserstein RL, Lazar NA. The ASA statement on p-values: context, process, and purpose. *The American Statistician.* 2016;70(2):129-133.
- Greenland S, et al. Statistical tests, P values, confidence intervals, and power: a guide to misinterpretations. *European Journal of Epidemiology.* 2016;31:337-350.
- Spiegelhalter D. *The Art of Statistics.* Pelican; 2019. Ch. 9.
