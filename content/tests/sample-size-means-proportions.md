# Sample size & power: means and proportions

> Calculates how many patients per arm a two-group trial needs to detect a clinically important difference in a mean (continuous outcome) or in a response rate (binary outcome) with a chosen type I error and power.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Planning a two-arm trial before data collection | Continuous or binary primary outcome | Realistic δ and σ (or p₁, p₂); independent patients; analysis matches the formula | 92 per arm to detect a 5 mmHg SBP difference (SD 12) with 80% power |

## Learning objectives

- Name the four linked quantities (α, power, effect, variability) and say which one you solve for
- Choose δ as the minimal clinically important difference, not the hoped-for effect
- Run the calculation in R and Python, round up, and inflate correctly for dropout
- Explain why post-hoc "observed power" is uninformative and how sensitive power is to the inputs

## What it is

Sample-size calculation answers one planning question: *how many participants do we need so that, if the true treatment effect is at least as large as the smallest effect we care about, our trial will detect it with high probability?* It is done **before** data collection and fixed in the protocol and statistical analysis plan (ICH E9).

Four quantities are linked, and fixing any three determines the fourth:

- **Significance level α:** the type I error rate, usually 0.05 two-sided
- **Power 1 − β:** the probability of a significant result when the true effect equals the planned effect, usually 80% or 90%
- **Effect to detect:** δ = μ₁ − μ₂ for means, or p₁ vs p₂ for proportions
- **Variability:** the standard deviation σ for means (for proportions it follows from p₁ and p₂)
- **Formula for means (per arm, equal allocation):** n = 2 (z₁₋α/₂ + z₁₋β)² σ² / δ², plus a small t-correction that `power.t.test` and statsmodels apply
- **Standardized effect size:** Cohen's d = δ / σ (means) or Cohen's h = 2·asin(√p₁) − 2·asin(√p₂) (proportions). Cohen's "small/medium/large" labels are rough conventions, not universal standards: what counts as large depends on the field and the clinical stakes, so anchor δ in the MCID.

**Where do δ and σ come from?**

- **δ should be the minimal clinically important difference (MCID):** the smallest change patients or clinicians would consider worthwhile. It is a clinical judgement, not a statistical one. Do not pick the effect you *hope* for, because trials sized for an optimistic effect are underpowered for the realistic one.
- **σ comes from pilot data, earlier trials in the same population, or registries and the literature.** A pilot SD is itself an estimate and is surprisingly imprecise. With only 15 patients, the 95% CI for a true SD of 12 mmHg runs from about 8.8 to 18.9 mmHg, and because n scales with σ², the required sample size could be off by a factor of 0.5 to 2.5. Prefer large published studies, use the upper end of the SD's confidence interval when in doubt, and re-check the SD at a blinded interim review if the protocol allows.

## The intuition

A sample-size calculation is a **signal-to-noise budget**. The signal is the effect you want to detect (δ = 5 mmHg). The noise is how much a difference in group means would bounce around by chance, which shrinks as you enrol more patients: the standard error of the difference is σ√(2/n). You are asking: *how large must n be so that, when the true effect is δ, the observed signal-to-noise ratio will usually clear the significance bar?*

Think of a diagnostic scale in a clinic. To tell whether a patient gained 200 g, a scale that jitters by ±500 g is useless, and you would need to average many weighings. A smaller jitter (or more weighings) makes a 200 g change visible. Sample size is the number of weighings.

With δ = 5, σ = 12 and n = 92 per arm, the expected t statistic if the planned effect is exactly true is about t = (5/12) × √(92/2) = 2.83. That is the typical trial result under H₁. Look at where it falls on the null curve below: it sits just beyond the critical value (about ±1.97), in the shaded tail. Because real trials scatter around this typical value, about 80% of them land past the critical value and 20% fall short. That 80% is the power.

<div class="viz" data-viz="t" data-stat="2.83" data-df="182"></div>

**What to notice:** the shaded area (p ≈ 0.005) is the p-value of the *typical* trial at n = 92. A trial landing at the edge of the shading would be only just significant. Halving the effect (2.5 mmHg) would need about four times the patients, because n scales with 1/δ² and with σ².

The sensitivity table in "Interpreting the output" shows the same idea in numbers: power drops from 80% to 61% if the true effect is 4 mmHg instead of 5.

## When to use it

Use these calculations when **all** of these are true:

- You are designing a **two-arm, parallel-group** study with 1:1 (or a known fixed) allocation
- The primary outcome is **continuous** (compared with a two-sample t-test) or **binary** (compared with a two-proportion test)
- You can name a clinically meaningful effect and a credible SD or control-group event rate

**Clinical scenario.** A trial will randomise adults with stage 1 hypertension to a new drug or standard care and compare the change in systolic blood pressure (SBP) at 12 weeks. The clinical team agrees that a 5 mmHg between-group difference is worth detecting, and earlier trials suggest an SD of about 12 mmHg. A second trial in the same programme uses a binary outcome: the proportion of patients achieving BP control (response), expected to be 40% on standard care and 55% on the new drug. *How many patients per arm does each trial need?*

**Don't use it when:**

- The outcome is **time to event** → use [sample size for time-to-event studies](#/test/sample-size-survival)
- The aim is to show **equivalence or non-inferiority** → use the design for a [TOST equivalence test](#/test/tost-equivalence), which has a different hypothesis and formula
- Outcomes are **repeated measures, clustered, or analysed with covariates** → baseline-adjusted [ANCOVA](#/test/ancova) needs fewer patients (roughly n × (1 − ρ²)); mixed-model and cluster designs need a design effect
- The study is **already finished** → do not compute "observed power" (see Common mistakes); report the confidence interval instead

## Assumptions

1. **The planning inputs are realistic.** The calculation is only as good as δ and σ (or p₁ and p₂). Check them against the literature and run a sensitivity table (see below).
2. **Observations are independent.** Patients are not clustered within clinics or families. If they are, inflate n by the design effect 1 + (m − 1)·ICC.
3. **The analysis will match the calculation.** `power.t.test` assumes an equal-variance two-sample t-test and `power.prop.test` a two-proportion z-test (χ² without continuity correction). If you plan a Welch test, a continuity-corrected test, or Fisher's exact test, expect slightly larger n.
4. **Outcome distributions are approximately normal (means) or the sample is not tiny (proportions).** Both formulas are large-sample approximations. For very rare outcomes (p near 0 or 1) or n below about 10 per cell, use an exact method or simulation.
5. **Equal allocation unless stated.** Equal arms minimise the total n for a given power; unequal allocation is covered below.

## Hypotheses

The design hypotheses are the usual two-sided superiority hypotheses; the sample size is chosen so the test rejects H₀ with probability 1 − β **if H₁ is true at the planned effect**.

- **H₀:** μ_treatment − μ_control = 0 (means) or p_treatment = p_control (proportions)
- **H₁:** the difference is not zero (two-sided), with planned effect δ = 5 mmHg or p₁ = 0.40 vs p₂ = 0.55
- **Design parameters:** α = 0.05 (two-sided), power = 80% and 90%

## Example data

These are design inputs, not patient data.

| Input | Means scenario (SBP) | Proportions scenario (BP control) |
|---|---|---|
| Effect to detect | δ = 5 mmHg | p₁ = 0.40 (standard care), p₂ = 0.55 (new drug) |
| Variability | σ = 12 mmHg | implied by p₁, p₂ |
| Standardized effect | d = 5/12 = 0.417 | h = 0.3015 |
| Two-sided α | 0.05 | 0.05 |
| Power | 80% and 90% | 80% and 90% |
| Allocation | 1:1 | 1:1 |
| Expected dropout | 15% | 15% |

## How to run it

### R

Means first, then proportions. Both use base R `stats`.

```r
delta <- 5; sigma <- 12                       # MCID (mmHg) and common SD

# 1. n per arm for 80% and 90% power (two-sided alpha = 0.05)
power.t.test(delta = delta, sd = sigma, sig.level = 0.05, power = 0.80)$n
power.t.test(delta = delta, sd = sigma, sig.level = 0.05, power = 0.90)$n

# 2. Round UP, then check the achieved power with n = 92 per arm
power.t.test(n = 92, delta = delta, sd = sigma)$power

# 3. Inflate for 15% dropout
ceiling(92 / (1 - 0.15))

# 4. Power vs true difference, n = 92 per arm
sapply(c(3, 4, 5, 6, 7), function(d) power.t.test(n = 92, delta = d, sd = sigma)$power)

# 5. 2:1 allocation (normal approximation), n in control arm
k <- 2; z <- qnorm(0.975) + qnorm(0.80)
ceiling((1 + 1/k) * (sigma * z / delta)^2)
```

**What the code does**

- `power.t.test(delta, sd, ...)` works on the raw scale: give δ and σ **separately**. Leave out exactly one of `n`, `power`, `delta`, `sd`, `sig.level` and R solves for it.
- Defaults to remember: `type = "two.sample"`, `alternative = "two.sided"`, `sig.level = 0.05`, equal-variance t-test. `n` is the number **per group**.
- `$n` is fractional (91.39). **Always round up**, which gives 92 per arm, and re-check the achieved power (80.3%).
- `ceiling(92 / (1 - 0.15))` inflates for 15% dropout: you must enrol n / (1 − dropout), not n × 1.15.
- The `sapply` line gives power at n = 92 for true differences of 3, 4, 5, 6 and 7 mmHg.
- `power.t.test` has no allocation-ratio argument, so for 2:1 allocation the code uses the normal-approximation formula n_control = (1 + 1/k)(σ z / δ)², where k = n_treatment / n_control and z = z₁₋α/₂ + z₁₋β.

**Output**

```text
[1] 91.38944
[1] 122.0139
[1] 0.8026334
[1] 109
[1] 0.3923003 0.6136600 0.8026334 0.9212093 0.9758906
[1] 68
```

```r
p1 <- 0.40; p2 <- 0.55

# 1. n per arm for 80% and 90% power (two-sided alpha = 0.05)
power.prop.test(p1 = p1, p2 = p2, sig.level = 0.05, power = 0.80)$n
power.prop.test(p1 = p1, p2 = p2, sig.level = 0.05, power = 0.90)$n

# 2. Achieved power with 173 per arm
power.prop.test(n = 173, p1 = p1, p2 = p2)$power

# 3. Dropout inflation (15%) applied to the rounded-up n
ceiling(173 / (1 - 0.15))
```

- `power.prop.test(p1, p2, ...)` uses the **pooled-variance under H₀ and unpooled under H₁** formula; the default is two-sided, α = 0.05, n per group.
- The result 172.8 is rounded up to 173 per arm for 80% power, and 230.83 to 231 for 90% power.
- Dropout inflation is applied to the rounded-up n: 173 / 0.85 = 203.5, so 204 per arm.

**Output**

```text
[1] 172.7999
[1] 230.8303
[1] 0.8004574
[1] 204
```

### Python

```python
import numpy as np
from scipy.stats import norm
from statsmodels.stats.power import TTestIndPower

delta, sigma = 5, 12                  # MCID (mmHg) and common SD
d = delta / sigma                     # standardized effect size (Cohen's d)
tp = TTestIndPower()

# 1. n per arm for 80% and 90% power (two-sided alpha = 0.05)
for pw in (0.80, 0.90):
    n = tp.solve_power(effect_size=d, alpha=0.05, power=pw, ratio=1.0)
    print(f"power {pw:.0%}: n per arm = {n:.3f} -> {int(np.ceil(n))}")

# 2. Achieved power with n = 92 per arm
print(f"power at n=92: {tp.power(effect_size=d, nobs1=92, alpha=0.05):.4f}")

# 3. Inflate for 15% dropout
print("with 15% dropout:", int(np.ceil(92 / (1 - 0.15))))

# 4. Power vs true difference, n = 92 per arm
for dl in (3, 4, 5, 6, 7):
    print(f"delta {dl}: power = {tp.power(effect_size=dl/sigma, nobs1=92, alpha=0.05):.3f}")

# 5. 2:1 allocation, normal approximation (same formula as the R tab)
z = norm.ppf(0.975) + norm.ppf(0.80)
n_ctrl = (1 + 1/2) * (sigma * z / delta) ** 2
print(f"2:1 -> control {np.ceil(n_ctrl):.0f}, treatment {np.ceil(2*n_ctrl):.0f}")
```

**What the code does**

- `TTestIndPower` takes the **standardized** effect size d = δ / σ (here 5/12 = 0.4167), not δ and σ separately. Passing δ = 5 directly is a very common error that gives an absurdly small n.
- `solve_power(...)` solves for whichever of `effect_size`, `nobs1`, `alpha`, `power` is left out. `nobs1` is the size of group 1 (per arm when `ratio = 1`).
- `ratio` is n₂ / n₁. `alternative` defaults to `"two-sided"`.
- `tp.power(...)` returns power for a given n; here it matches R (0.8026 at n = 92).
- For the 2:1 split, the code uses the same normal approximation as the R tab so that both languages agree. Using `solve_power(..., ratio=2.0)` includes the t-correction and would give 69 control and 137 treatment (from 68.46 and 136.92).

**Output**

```text
power 80%: n per arm = 91.389 -> 92
power 90%: n per arm = 122.014 -> 123
power at n=92: 0.8026
with 15% dropout: 109
delta 3: power = 0.392
delta 4: power = 0.614
delta 5: power = 0.803
delta 6: power = 0.921
delta 7: power = 0.976
2:1 -> control 68, treatment 136
```

```python
import numpy as np
from scipy.stats import norm
from statsmodels.stats.power import NormalIndPower
from statsmodels.stats.proportion import proportion_effectsize

p1, p2 = 0.40, 0.55

# A. statsmodels: arcsine (Cohen's h) approximation
h = proportion_effectsize(p2, p1)
print(f"Cohen's h = {h:.4f}")
for pw in (0.80, 0.90):
    n = NormalIndPower().solve_power(effect_size=h, alpha=0.05, power=pw, ratio=1.0)
    print(f"statsmodels  power {pw:.0%}: n per arm = {n:.3f}")

# B. Exact reproduction of R's power.prop.test formula
def n_prop_test(p1, p2, alpha=0.05, power=0.80):
    pbar = (p1 + p2) / 2
    za, zb = norm.ppf(1 - alpha / 2), norm.ppf(power)
    num = (za * np.sqrt(2 * pbar * (1 - pbar)) + zb * np.sqrt(p1*(1-p1) + p2*(1-p2))) ** 2
    return num / (p1 - p2) ** 2

for pw in (0.80, 0.90):
    n = n_prop_test(p1, p2, power=pw)
    print(f"R formula    power {pw:.0%}: n per arm = {n:.3f}")

# C. Dropout inflation (15%) on the rounded-up n
print("with 15% dropout:", int(np.ceil(173 / (1 - 0.15))))
```

**What the code does**

- `proportion_effectsize(p2, p1)` computes Cohen's h = 2·asin(√p₂) − 2·asin(√p₁) (0.3015). `NormalIndPower` then uses the normal approximation on that **arcsine-transformed** scale.
- Because R's `power.prop.test` works on the original proportion scale with explicit variances, the two methods give slightly different answers: 172.66 vs 172.80 at 80% power, and 231.14 vs 230.83 at 90% power. Neither is "wrong"; they are different approximations.
- `n_prop_test()` re-implements the `power.prop.test` formula, n = (z₁₋α/₂ √(2p̄(1 − p̄)) + z₁₋β √(p₁(1 − p₁) + p₂(1 − p₂)))² / (p₁ − p₂)², with p̄ = (p₁ + p₂)/2. It matches R to all printed digits (172.800 and 230.830).
- `norm.ppf(1 - alpha/2)` is the two-sided critical value. `np.ceil` rounds up.

**Output**

```text
Cohen's h = 0.3015
statsmodels  power 80%: n per arm = 172.659
statsmodels  power 90%: n per arm = 231.142
R formula    power 80%: n per arm = 172.800
R formula    power 90%: n per arm = 230.830
with 15% dropout: 204
```

SAS users can run the same calculations with PROC POWER (`twosamplemeans` and `twosamplefreq` statements); they are listed under References.

## Interpreting the output

**Means (SBP, δ = 5 mmHg, σ = 12 mmHg)**

- **80% power: 91.39, rounded up to 92 per arm** (184 in total). R and statsmodels agree (91.389). At n = 92 the achieved power is 80.3%.
- **90% power: 122.01, rounded up to 123 per arm** (246 in total). Going from 80% to 90% power costs about 34% more patients.
- **Dropout:** if 15% of patients are expected to be lost before the primary endpoint, enrol 109 per arm (92 / 0.85 = 108.2) for 80% power. This keeps the *analysable* sample at 92 only if dropout is non-informative; the dropout itself can bias the estimate.
- **2:1 allocation (normal approximation):** 68 control and 136 treatment, so 204 in total, versus 184 with 1:1. Unequal allocation always costs total sample size; it is justified only by other reasons (more safety data on the new drug, a more attractive trial for patients).

**Sensitivity: power at n = 92 per arm, SD = 12**

| True difference (mmHg) | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|
| Power | 39.2% | 61.4% | 80.3% | 92.1% | 97.6% |

The table shows how fragile power is: if the true effect is 4 mmHg instead of 5, power falls from 80% to 61%. Present a table like this in the protocol, and consider a sensitivity analysis for σ as well.

**Proportions (40% vs 55%)**

- **80% power: 172.80 (R) → 173 per arm; 90% power: 230.83 → 231 per arm.** With 173 per arm, achieved power is 80.0%.
- **statsmodels (arcsine):** 172.66 and 231.14, which round to 173 and 232 per arm. The 90% result differs by one patient per arm from R only because of the different approximation.
- **Python re-implementation of the R formula** gives 172.800 and 230.830, identical to R. Use it when you need your Python result to match an R-based protocol.
- **Dropout:** 173 / 0.85 = 203.5, so 204 per arm.

**How to report it:**
> We calculated that 92 patients per arm would provide 80% power to detect a between-group difference in systolic blood pressure of 5 mmHg (SD 12 mmHg, two-sided α = 0.05, two-sample t-test); allowing for 15% attrition, we will randomise 109 patients per arm. For the binary outcome, 173 patients per arm give 80% power to detect an increase in response from 40% to 55% (two-sided α = 0.05).

⚠️ **Do not compute "observed" or post-hoc power.** Power calculated from the observed effect of a completed trial is just a one-to-one transformation of the p-value: a non-significant result always has low "observed power", so it adds no information. After the trial, interpret the confidence interval instead: it shows which effects are compatible with the data.

## Check your understanding

<details class="check"><summary>1. The pilot SD was 12 mmHg, but a larger published trial reports 15 mmHg. By roughly what factor does the required n change?</summary>

n scales with σ², so it grows by (15/12)² ≈ 1.56. The 92 per arm becomes about 143 per arm. This is why the SD input deserves as much scrutiny as δ.
</details>

<details class="check"><summary>2. A trial with n = 92 per arm ends with p = 0.30. A colleague computes "observed power" of 18% and concludes the trial was simply underpowered. Is that reasoning useful?</summary>

No. Observed power is a one-to-one transformation of the p-value, so a non-significant result always yields low observed power and tells you nothing new. Report the effect estimate and its 95% CI, which shows which effects are still compatible with the data.
</details>

<details class="check"><summary>3. What happens to power if the true difference is 4 mmHg rather than the planned 5, with the same n = 92 per arm?</summary>

Power falls from 80.3% to about 61.4% (see the sensitivity table). A smaller true effect means less signal against the same noise.
</details>

<details class="check"><summary>4. You need 92 analysable patients per arm and expect 15% dropout. Is 92 × 1.15 = 106 enough?</summary>

No. Dropout removes 15% of those *enrolled*, so enrol n / (1 − 0.15) = 109 per arm. Multiplying by 1.15 gives only 106, and about 90 would remain analysable.
</details>

## Common mistakes

- **Not rounding up.** 91.39 patients means 92. Rounding down leaves you slightly under the target power.
- **Using the effect you hope for rather than the MCID.** Optimistic effects give small trials that miss realistic benefits.
- **Treating a pilot SD as truth.** Small pilots give imprecise SDs, and n scales with σ². Use an upper confidence limit or a larger external source.
- **Passing δ instead of d = δ / σ to statsmodels.** `TTestIndPower` expects the standardized effect size.
- **Confusing per-arm and total n.** Both `power.t.test` and `statsmodels` return n **per group**.
- **Multiplying by 1.15 for 15% dropout.** Use n / (1 − 0.15); 92 × 1.15 = 106 is too few.
- **Ignoring multiplicity.** With several primary comparisons, use the adjusted α in the calculation (see [multiple comparisons](#/test/multiple-comparisons)).
- **Post-hoc (observed) power.** It is not informative; report confidence intervals.
- **Assuming R and Python agree for proportions.** They use different approximations, so state the method in the protocol.

## Related tests

- [Sample size for time-to-event studies](#/test/sample-size-survival): events and patients for a hazard ratio
- [Independent t-test](#/test/independent-t): the analysis that the means calculation assumes
- [Two-proportion z-test](#/test/two-proportion-z): the analysis that the proportions calculation assumes
- [ANCOVA](#/test/ancova): baseline adjustment reduces the required n
- [TOST equivalence test](#/test/tost-equivalence): equivalence and non-inferiority designs need a different calculation
- SAS users can use PROC POWER (`twosamplemeans`, `twosamplefreq`) for the same calculations.

## Key takeaways

- Four quantities are linked (α, power, effect, variability); fix three and solve for the fourth.
- Set δ to the **minimal clinically important difference**, not the effect you hope for.
- n scales with σ²/δ²: halving the effect needs about four times the patients.
- Always round **up**, and inflate for dropout with n / (1 − dropout).
- Never compute observed power after the trial; report the confidence interval instead.

## References

- Chow S-C, Shao J, Wang H, Lokhnygina Y. *Sample Size Calculations in Clinical Research.* 3rd ed. Chapman & Hall/CRC; 2017.
- ICH E9: *Statistical Principles for Clinical Trials.* 1998 (section 3.5, sample size).
- Hoenig JM, Heisey DM. The abuse of power: the pervasive fallacy of power calculations for data analysis. *Am Stat.* 2001;55:19-24.
- R documentation: `?power.t.test`, `?power.prop.test`; statsmodels documentation: `statsmodels.stats.power.TTestIndPower`, `NormalIndPower`
