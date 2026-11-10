# Type I and II errors, and power

> Every test can be wrong in two ways, a false alarm or a missed signal; power is the chance of catching a real effect, and multiple testing quietly multiplies the false alarms.

## Learning objectives

- Fill in the 2×2 table of truth versus decision and name each cell
- Define α, β, and power, and say what drives power: effect size, SD, n, and α
- Calculate how multiple testing inflates the chance of at least one false positive
- Explain why "observed power" computed after the study is meaningless

## The intuition

Think of a smoke alarm. It can fail in two ways. It can ring when there is no fire (a **false alarm**), or stay silent during a real fire (a **miss**). You can make it more sensitive and get fewer misses, but more false alarms; the way to improve both is a better detector. A hypothesis test works the same way, and the "better detector" is a larger sample or a less noisy measurement.

Take our running trial: a drug that truly lowers SBP by 7 mmHg versus placebo, patient SD 12 mmHg, 15 patients per arm. Whether the trial detects the effect depends on how far apart the two groups' averages are relative to the noise. The picture below shows one such trial (the data from [Hypothesis testing and p-values](#/test/hypothesis-testing)); the drug group is lower, but with this much overlap the evidence is weak.

<div class="viz" data-viz="dots" data-unit="mmHg" data-groups='{"Drug":[133,130,144,126,125,125,133,128,137,134,133,131,124,130,121],"Placebo":[133,140,122,124,154,129,156,147,139,128,130,136,122,137,126]}' data-caption="A real 7 mmHg effect, but the noise is large with only 15 per arm."></div>

For this trial t = −1.45 on 28 degrees of freedom, and the test does not reach significance, a **miss** (type II error) given that the drug really works:

<div class="viz" data-viz="t" data-stat="-1.45" data-df="28"></div>

The two tails you see shaded are 16% of the null curve. To be "significant" a trial needs |t| ≥ 2.05, the point where only 5% of the curve lies beyond. That 5% is the false-alarm rate α.

## The 2×2 table of truth vs decision

|  | **Test says "significant" (reject H₀)** | **Test says "not significant" (fail to reject)** |
|---|---|---|
| **H₀ true** (no real effect) | **Type I error** (false positive), probability **α** | Correct, probability 1 − α |
| **H₁ true** (real effect) | Correct, probability **1 − β = power** | **Type II error** (false negative), probability **β** |

- **α** (significance level) is set by you, usually 0.05. It is the false-positive rate *when H₀ is true*.
- **β** is the chance of missing a real effect of a stated size.
- **Power = 1 − β.** Trials are usually designed for 80% or 90% power.

Neither α nor power is the probability that a particular published "positive" result is true. That depends on how many tested hypotheses were real in the first place (see the Ioannidis paper below).

## What drives power

Power increases when:

1. **The true effect is larger.** Big effects are easy to see.
2. **The SD is smaller.** Less noise means a clearer signal.
3. **The sample size is larger.** The standard error shrinks like 1/√n.
4. **α is larger.** A looser threshold catches more real effects, at the cost of more false alarms.

For our drug (7 mmHg, SD 12, 15 per arm) power is only about 34%: two trials out of three would miss it. About 47 patients per arm are needed for 80% power. How to plan this properly is covered in [Sample size & power: means and proportions](#/test/sample-size-means-proportions) and, for time-to-event outcomes, in [Sample size for time-to-event studies](#/test/sample-size-survival).

## Multiple testing inflates false positives

If you run k independent tests, each at α = 0.05, and all the null hypotheses are true, the chance of at least one false positive is

**P(≥ 1 false positive) = 1 − 0.95^k**

| k tests | 1 | 5 | 10 | 20 |
|---|---|---|---|---|
| P(≥ 1 false positive) | 5.0% | 22.6% | 40.1% | 64.2% |

A trial analysing 20 outcomes or 20 subgroups has a better-than-even chance of finding "something significant" even when the treatment does nothing. Remedies: pre-specify one primary outcome, adjust the threshold (Bonferroni, Holm, or the methods in [Multiple comparisons](#/test/multiple-comparisons)), and label unplanned analyses as exploratory. Ioannidis (2005) argued that, because of low power, multiple testing, flexible analysis, and publication bias, a large share of published positive findings may be false.

## Why "observed power" is meaningless

Some journals or reviewers ask for the "post-hoc power" of a finished study, computed by plugging the *observed* effect into a power formula. This is circular. Observed power is just a re-expression of the p-value: it is about 50% when p = 0.05 and falls as p grows, so a non-significant result always shows "low power". It gives no information beyond p, and it does not tell you whether the study could have detected a *clinically relevant* effect. The right approach is to compute power **before** the study for the smallest effect worth detecting, and **after** the study to look at the **confidence interval**, which shows which effects are compatible with the data.

## How to run it

The simulations repeat the drug-versus-placebo trial thousands of times at several sample sizes and compare the share of significant trials with the exact formula (`power.t.test` in R; `TTestIndPower` in statsmodels). They also check the 1 − 0.95^k formula and demonstrate observed power. R and Python use different random streams, so simulated numbers differ slightly.

### R

```r
# Setting: drug lowers SBP by 7 mmHg (true), SD 12, alpha = 0.05, two-sided
power.t.test(n = 15, delta = 7, sd = 12, sig.level = 0.05)$power     # theory, n = 15/arm
power.t.test(power = 0.80, delta = 7, sd = 12)$n                     # n/arm for 80% power

# Estimate power by simulation: share of trials with p < 0.05 when H1 is true
set.seed(42); nsim <- 3000
sim_power <- function(n, delta, sd = 12)
  mean(replicate(nsim, t.test(rnorm(n, delta, sd), rnorm(n, 0, sd), var.equal = TRUE)$p.value < 0.05))
for (n in c(15, 30, 60)) cat("n =", n, "per arm: simulated power =", sim_power(n, 7),
    " theory =", round(power.t.test(n = n, delta = 7, sd = 12)$power, 3), "\n")
cat("Type I error (delta = 0, n = 15):", sim_power(15, 0), "\n")

# Multiple testing: k independent true-null tests, P(at least one p < 0.05)
k <- c(1, 5, 10, 20)
cat("Theory 1 - 0.95^k:", round(1 - 0.95^k, 3), "\n")
fwer <- sapply(k, function(kk) mean(replicate(1000, any(replicate(kk, t.test(rnorm(15), rnorm(15))$p.value) < 0.05))))
cat("Simulated:         ", fwer, "\n")

# "Observed power": a trial with exactly p = 0.05 (n = 15/arm), plugging in the observed effect
t_obs <- qt(0.975, 28); d_obs <- t_obs * 12 * sqrt(2 / 15)
cat("Observed power when p = 0.05:", round(power.t.test(n = 15, delta = d_obs, sd = 12)$power, 3), "\n")
t_obs <- qt(1 - 0.30 / 2, 28); d_obs <- t_obs * 12 * sqrt(2 / 15)
cat("Observed power when p = 0.30:", round(power.t.test(n = 15, delta = d_obs, sd = 12)$power, 3), "\n")
```

**What the code does**

- `power.t.test` gives the exact (noncentral t) power for a two-sample, two-sided t-test; supplying `power` instead of `n` solves for the sample size per arm.
- `sim_power()` runs a fake trial `nsim` times and reports the share with p < 0.05; with a true effect this estimates power, with `delta = 0` it estimates the Type I error.
- The multiple-testing block runs k true-null tests per "study" and counts studies with at least one p < 0.05.
- The last block converts a p-value into the observed effect and then computes power at that effect, which is the circular "observed power".

**Output**

```text
[1] 0.3382965
[1] 47.11207
n = 15 per arm: simulated power = 0.3343333  theory = 0.338 
n = 30 per arm: simulated power = 0.5963333  theory = 0.603 
n = 60 per arm: simulated power = 0.8883333  theory = 0.887 
Type I error (delta = 0, n = 15): 0.04266667 
Theory 1 - 0.95^k: 0.05 0.226 0.401 0.642 
Simulated:          0.05 0.211 0.377 0.633 
Observed power when p = 0.05: 0.507 
Observed power when p = 0.30: 0.174 
```

### Python

```python
import numpy as np
from scipy import stats
from statsmodels.stats.power import TTestIndPower

# Setting: drug lowers SBP by 7 mmHg (true), SD 12, alpha = 0.05, two-sided
pw = TTestIndPower()
d = 7 / 12                                             # Cohen's d = delta / SD
print("Theory power, n = 15/arm:", round(pw.power(d, nobs1=15, alpha=0.05), 4))
print("n/arm for 80% power:", round(pw.solve_power(d, power=0.80, alpha=0.05), 2))

# Estimate power by simulation: share of trials with p < 0.05 when H1 is true
rng = np.random.default_rng(42); nsim = 5000
def sim_power(n, delta, sd=12):
    x = rng.normal(delta, sd, (nsim, n)); y = rng.normal(0, sd, (nsim, n))
    return np.mean(stats.ttest_ind(x, y, axis=1).pvalue < 0.05)
for n in (15, 30, 60):
    print(f"n = {n} per arm: simulated power = {sim_power(n, 7)}  theory = {pw.power(d, nobs1=n, alpha=0.05):.3f}")
print("Type I error (delta = 0, n = 15):", sim_power(15, 0))

# Multiple testing: k independent true-null tests, P(at least one p < 0.05)
ks = [1, 5, 10, 20]
print("Theory 1 - 0.95^k:", [round(1 - 0.95**k, 3) for k in ks])
def fwer(k, reps=5000):
    p = stats.ttest_ind(rng.normal(size=(reps, k, 15)), rng.normal(size=(reps, k, 15)), axis=2).pvalue
    return np.mean((p < 0.05).any(axis=1))
print("Simulated:        ", [round(float(fwer(k)), 3) for k in ks])

# "Observed power": a trial with exactly p = 0.05 (n = 15/arm), plugging in the observed effect
for p_obs in (0.05, 0.30):
    d_obs = stats.t.ppf(1 - p_obs / 2, 28) * np.sqrt(2 / 15)
    print(f"Observed power when p = {p_obs}:", round(pw.power(d_obs, nobs1=15, alpha=0.05), 3))
```

**What the code does**

- `TTestIndPower().power()` takes the standardized effect d = Δ/SD and `nobs1` (patients in the first arm; the ratio defaults to 1).
- `solve_power()` leaves one quantity unknown (here the sample size) and solves for it.
- The simulations generate all trials as one matrix, one row per trial, so 5,000 t-tests take a fraction of a second.
- The 3-D array in `fwer()` holds `reps` studies × `k` tests × 15 patients; the `axis=2` t-test runs all of them at once.

**Output**

```text
Theory power, n = 15/arm: 0.3385
n/arm for 80% power: 47.11
n = 15 per arm: simulated power = 0.345  theory = 0.339
n = 30 per arm: simulated power = 0.6172  theory = 0.603
n = 60 per arm: simulated power = 0.8856  theory = 0.887
Type I error (delta = 0, n = 15): 0.052
Theory 1 - 0.95^k: [0.05, 0.226, 0.401, 0.642]
Simulated:         [0.05, 0.217, 0.405, 0.634]
Observed power when p = 0.05: 0.507
Observed power when p = 0.3: 0.175
```

## Interpreting the output

- **Power at n = 15 per arm is 0.34 (theory), and simulation agrees (0.33 in R, 0.35 in Python).** The trial would miss a real 7 mmHg effect about two times in three.
- **Power grows with n:** 0.34 at 15, 0.60 at 30, 0.89 at 60 per arm. About 47 per arm are needed for 80%.
- **Type I error when H₀ is true: 0.043 (R) and 0.052 (Python).** Both are consistent with α = 0.05 given simulation noise.
- **Multiple testing:** with 20 independent tests the simulated chance of at least one false positive is 0.63 in both languages, matching 1 − 0.95²⁰ = 0.64.
- **Observed power is 0.51 at p = 0.05 and 0.17 at p = 0.30.** It is a one-to-one function of p, so it tells you nothing new.

**How to report it:**
> The trial was designed with 80% power to detect a 7 mmHg difference in SBP (SD 12 mmHg, two-sided α = 0.05), requiring 48 patients per arm.

## Check your understanding

<details class="check"><summary>1. A study of 20 outcomes finds one significant at p = 0.04, with no adjustment. The authors say "the drug works for outcome 14". What is the problem?</summary>

With 20 independent tests at α = 0.05, the chance of at least one false positive is 1 − 0.95²⁰ = 64%, even if the drug does nothing. One p = 0.04 among 20 tests is what chance alone often produces. Unless outcome 14 was the pre-specified primary outcome, treat it as hypothesis-generating, or adjust for multiplicity.
</details>

<details class="check"><summary>2. A trial has 80% power. The result is non-significant (p = 0.20). Does that mean the drug has no effect?</summary>

No. 80% power means that **if** the true effect equals the planned effect, the trial would be significant 80% of the time. A p = 0.20 does not show the effect is zero. Check the CI: if it is compatible with the planned effect, the result is inconclusive; if the CI excludes any effect of interest, you have stronger evidence of no important effect.
</details>

<details class="check"><summary>3. Using this lesson's numbers (true effect 7 mmHg, SD 12, 15 per arm, power about 34%), what happens to power if the SD is reduced to 8 mmHg through tighter measurement, with everything else fixed?</summary>

Power rises sharply. The standardized effect goes from 7/12 = 0.58 to 7/8 = 0.88, so the signal-to-noise ratio is larger and the test statistic is further from 0 on average. With 15 per arm, power for d = 0.88 is about 64%, compared with 34%. Reducing noise is as valuable as adding patients.
</details>

<details class="check"><summary>4. A reviewer asks you to compute "post-hoc power" using your observed effect. Your result had p = 0.30. What will it show, and should you report it?</summary>

It will show low power (about 17% in our example), but that figure is only a restatement of p = 0.30, so it adds nothing and reasoning from it is circular. Do not report it. Instead report the 95% CI and, if relevant, state the power the study had for a clinically meaningful effect, computed with pre-specified assumptions.
</details>

## Key takeaways

- Type I error (α) is a false positive when H₀ is true; Type II error (β) is a miss when H₁ is true; power = 1 − β.
- Power rises with larger effect, smaller SD, larger n, and larger α; sample size must be planned before the study.
- k independent true-null tests give P(≥ 1 false positive) = 1 − 0.95^k: 64% at k = 20.
- Observed (post-hoc) power is a function of p and is uninformative; interpret the CI instead.
- Low power makes non-significant results hard to interpret and makes significant findings overestimate effects.

## References

- Ioannidis JPA. Why most published research findings are false. *PLoS Medicine.* 2005;2(8):e124.
- Hoenig JM, Heisey DM. The abuse of power: the pervasive fallacy of power calculations for data analysis. *The American Statistician.* 2001;55(1):19-24.
- Altman DG, Bland JM. Absence of evidence is not evidence of absence. *BMJ.* 1995;311:485.
