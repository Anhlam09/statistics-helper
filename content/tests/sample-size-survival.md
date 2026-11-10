# Sample size for time-to-event studies

> Calculates the number of events, and then the number of patients, needed to detect a given hazard ratio in a two-arm trial with a time-to-event primary endpoint such as overall survival.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Planning a two-arm trial with a time-to-event primary endpoint | Time to event with censoring | Proportional hazards; exponential planning model; uniform accrual | Events and patients needed to detect HR 0.70 in overall survival |

## Learning objectives

- Explain why power in a survival trial depends on the number of events, not patients
- Compute required events with Schoenfeld's formula and convert them to patients
- Run the calculation in R and Python and read events, P(event), and N
- Report a survival sample size with all assumptions, and adjust for dropout and interim looks

## What it is

In time-to-event trials, **statistical power depends on the number of events** (deaths, progressions), not directly on the number of patients. The planning therefore has two steps: (1) compute the required number of events, and (2) work out how many patients must be enrolled, and followed for how long, to observe that many events.

- **Step 1, Schoenfeld's formula:** d = (z₁₋α/₂ + z₁₋β)² / (p (1 − p) (log HR)²), where p is the proportion of patients allocated to the treatment arm (p = 0.5 for 1:1). It assumes proportional hazards and applies to the log-rank test and the Cox model.
- **Step 2, patients from events:** N = d / P(event), where P(event) is the overall probability that a randomly chosen patient has the event during the study. Under exponential survival with hazard λ, uniform accrual over A months, and F further months of follow-up after the last patient enters, P(event) = 1 − (1/A) ∫₀ᴬ exp(−λ (F + u)) du. A patient entering at time u has A + F − u months of potential follow-up, so the integral averages survival over the uniform entry times (substituting v = A − u gives the form shown). The overall P(event) is the allocation-weighted average of the two arms.
- **Effect size:** the hazard ratio (HR) with its 95% CI. HR = 0.70 means a 30% lower instantaneous event rate in the treatment arm.

**Where do the inputs come from?** The HR should be the smallest clinically worthwhile benefit, not the optimistic one. The control-arm median (or event rate at a fixed time) comes from historical data or registries; under an exponential model the hazard is λ = ln 2 / median.

## The intuition

Think of a survival trial as a **news poll that only counts people who give an answer**. Enrolling 1,000 patients who all survive past the end of the study tells you almost nothing about survival: the information arrives only when an event happens. A patient who dies contributes information about the hazard ratio; a patient who is still alive at the final analysis contributes only "no event yet". So the unit of information is the **event**.

The signal is the log hazard ratio, log(HR). The noise is how precisely the log HR can be estimated, and that precision grows with the number of events d (its standard error is roughly 2/√d for 1:1 allocation). To detect a smaller signal, you need more events, and the requirement grows with the **square** of 1/log(HR). The dot plot below uses the lesson's own formula (80% power, two-sided α = 0.05, 1:1) to show how fast the events needed climb as the target HR moves toward 1.

<div class="viz" data-viz="dots" data-unit="events needed" data-groups='{"HR 0.50, 0.60, 0.70, 0.80 (left to right)":[65,120,247,631]}' data-caption="Required events for HR 0.50, 0.60, 0.70 and 0.80. Notice that going from HR 0.70 to 0.80 more than doubles the requirement (247 to 631), while the jump from 0.50 to 0.60 adds only about 55 events."></div>

Events are only half of the plan. The second step answers "how many patients do I have to enrol to *see* 247 deaths?". That depends on how long people are followed and how quickly they die: in our example about 67% of patients are expected to have an event during the 36-month study, so 247 events need about 370 patients. A faster-dying population, or longer follow-up, means fewer patients per event.

## When to use it

Use this approach when **all** of these are true:

- The primary endpoint is **time to a single event type** (overall survival, progression-free survival, time to first MACE)
- You plan to compare two arms with a **log-rank test or Cox model**, and a constant hazard ratio is reasonable
- You can specify the **accrual period, follow-up period, and control-arm event rate**

**Clinical scenario.** A phase III trial randomises patients with advanced cancer 1:1 to a new drug or standard therapy, with overall survival as the primary endpoint. Median survival on standard therapy is 12 months. The trial wants 80% power (two-sided α = 0.05) to detect HR = 0.70, which would raise median survival to about 17 months. Patients will be accrued over 24 months and followed for a further 12 months after the last patient is enrolled. *How many events, and how many patients, are needed?*

**Don't use it when:**

- The endpoint is a **binary response** or a **continuous** measure → use [sample size for means and proportions](#/test/sample-size-means-proportions)
- Hazards are **clearly non-proportional** (crossing curves, delayed effect, cure fraction) → use simulation or methods based on the restricted mean survival time
- **Competing risks** dominate → see [competing risks](#/test/competing-risks) for the appropriate estimand and design
- The analysis is only **descriptive** → see [Kaplan–Meier](#/test/kaplan-meier)

## Assumptions

1. **Proportional hazards.** The HR is constant over time. Check plausibility from earlier trials and with [Cox model](#/test/cox-ph) diagnostics after the study.
2. **Exponential survival in the planning model.** Event times in each arm are exponential, so the hazard is constant. This is a planning simplification; with a Weibull or piecewise model, compute P(event) with the actual survival function.
3. **Uniform accrual** over A months, with no early departure other than the planned follow-up. Real accrual is often slower at the start, which delays events.
4. **Non-informative censoring and no loss to follow-up** in the basic calculation. Inflate N for expected dropout: N / (1 − dropout).
5. **Fixed allocation and a single final analysis.** Interim looks need an adjustment (below).

## Hypotheses

- **H₀:** HR = 1 (the hazard of the event is the same in both arms)
- **H₁:** HR ≠ 1 (two-sided), with planned HR = 0.70
- **Design parameters:** α = 0.05 (two-sided), power = 80%, allocation 1:1 (p = 0.5)

## Example data

These are design inputs, not patient data.

| Input | Value |
|---|---|
| Hazard ratio to detect | 0.70 |
| Two-sided α | 0.05 |
| Power | 80% |
| Allocation | 1:1 (p = 0.5) |
| Control median survival | 12 months (λ_c = ln 2 / 12 = 0.0578 per month) |
| Treatment hazard | λ_t = 0.70 × λ_c = 0.0404 per month |
| Accrual period A | 24 months (uniform) |
| Additional follow-up F | 12 months |
| Total study length | 36 months |

## How to run it

### R

```r
hr <- 0.70; alpha <- 0.05; power <- 0.80; p <- 0.5   # p = share in treatment arm
A <- 24; Fu <- 12; med <- 12                           # accrual, follow-up, control median (months)
lam_c <- log(2) / med; lam_t <- lam_c * hr             # exponential hazards

# 1. Schoenfeld: required number of events
z <- qnorm(1 - alpha/2) + qnorm(power)
d <- z^2 / (p * (1 - p) * log(hr)^2)
d

# 2. P(event) = 1 - average of S(t) over entry times, i.e. 1 - (1/A) * integral S(F + u) du, u in [0, A]
p_event <- function(lam) 1 - integrate(function(u) exp(-lam * (Fu + u)), 0, A)$value / A
pe_c <- p_event(lam_c); pe_t <- p_event(lam_t)
pe_bar <- (1 - p) * pe_c + p * pe_t
c(control = pe_c, treatment = pe_t, average = pe_bar)

# 3. Convert events to patients and round up (even total for 1:1)
N <- d / pe_bar
c(events = ceiling(d), N = ceiling(N), per_arm = ceiling(N / 2))

# 4. Freedman's formula (events) for comparison
(1 + hr)^2 * z^2 / (1 - hr)^2
```

**What the code does**

- `z` is z₁₋α/₂ + z₁₋β = 1.960 + 0.842 = 2.802. Schoenfeld's `d` is the number of **events**, not patients.
- `p_event()` integrates the survival function over the entry times. A patient who enters at time u (0 to A) has A + F − u months before the analysis, so the integral over u of exp(−λ (A + F − u)) is the same as the integral of exp(−λ (F + u)) over u in [0, A], which is the form used in the code.
- `integrate(f, 0, A)$value` does the numerical integration; dividing by A averages over the uniform entry times.
- `pe_bar` averages P(event) across arms with weights (1 − p) and p. This matters: the treatment arm has fewer events, so the average is lower than the control value.
- `N <- d / pe_bar`, then round up. For 1:1 allocation, the per-arm number is N / 2, rounded up to a whole patient.
- The last line is Freedman's older formula d = (1 + HR)² z² / (1 − HR)², shown for comparison. It gives a slightly larger number of events than Schoenfeld's because it approximates the variance differently.

**Output**

```text
[1] 246.7871
  control treatment   average 
0.7294947 0.6060268 0.6677607 
 events       N per_arm 
    247     370     185 
[1] 252.0362
```

### Python

```python
import numpy as np
from scipy.stats import norm
from scipy.integrate import quad

hr, alpha, power, p = 0.70, 0.05, 0.80, 0.5     # p = share in treatment arm
A, F, med = 24, 12, 12                           # accrual, follow-up, control median (months)
lam_c = np.log(2) / med; lam_t = lam_c * hr      # exponential hazards

# 1. Schoenfeld: required number of events
z = norm.ppf(1 - alpha / 2) + norm.ppf(power)
d = z**2 / (p * (1 - p) * np.log(hr) ** 2)
print(f"events d = {d:.3f}")

# 2. P(event) = 1 - (1/A) * integral of S(F + u) du over u in [0, A]
def p_event(lam):
    val, _ = quad(lambda u: np.exp(-lam * (F + u)), 0, A)
    return 1 - val / A
pe_c, pe_t = p_event(lam_c), p_event(lam_t)
pe_bar = (1 - p) * pe_c + p * pe_t
print(f"P(event): control {pe_c:.4f}, treatment {pe_t:.4f}, average {pe_bar:.4f}")

# 3. Convert events to patients and round up
N = d / pe_bar
print(f"events {int(np.ceil(d))}, N = {N:.1f} -> {int(np.ceil(N))}, per arm {int(np.ceil(N/2))}")

# 4. Freedman's formula (events) for comparison
print(f"Freedman events = {(1 + hr)**2 * z**2 / (1 - hr)**2:.1f}")
```

**What the code does**

- `norm.ppf` gives the normal quantiles, and the Schoenfeld formula is identical to the R code.
- `scipy.integrate.quad` returns a tuple `(value, abs_error)`; the code keeps only the value. It is the counterpart of R's `integrate`.
- `np.exp(-lam * (F + u))` is the survival probability at the end of the study for a patient whose entry time is A − u (u is measured backwards from the end of accrual, so u + F is that patient's follow-up time), so averaging over u is exact for uniform accrual.
- `np.ceil` rounds events and patients up. Everything matches R to the printed digits.

**Output**

```text
events d = 246.787
P(event): control 0.7295, treatment 0.6060, average 0.6678
events 247, N = 369.6 -> 370, per arm 185
Freedman events = 252.0
```

SAS users can use PROC POWER with the `twosamplesurvival` statement for comparable calculations; it is listed under References.

## Interpreting the output

- **Required events: 246.79, rounded up to 247.** About 247 deaths are needed in total, split across the two arms, for 80% power at HR = 0.70. This number is fixed by α, power and the HR alone; accrual and follow-up change only how long it takes to observe them.
- **P(event) by arm: 73.0% control and 60.6% treatment, average 66.8%.** Under exponential survival with a 12-month median, about 73% of control patients are expected to die within the study window. With HR = 0.70, 61% of treated patients are.
- **Patients: 247 / 0.6678 = 369.6, rounded up to 370 in total, or 185 per arm.** Roughly 1.5 patients must be enrolled for every event.
- **Freedman: 252 events,** about 2% more than Schoenfeld (247) at this HR. The difference grows when the HR moves away from 1. Both are approximations, so state which one you used.
- **Dropout:** if 10% of patients are expected to be lost to follow-up, divide N by 0.90 (370 / 0.90 = 411.1, or 412 patients).
- **Trade-offs:** longer follow-up reduces the patients needed, because more events occur per patient, and a lower HR (a bigger effect) reduces the events needed at all. With the total study length held fixed, a shorter accrual period gives every patient more follow-up and so needs fewer patients; with F held fixed, it shortens the trial instead. Neither changes the number of events needed. Changing the control median changes P(event) and so the patient count, but not the number of events.
- **Group-sequential designs:** if the trial has interim efficacy analyses (for example O'Brien–Fleming boundaries), the final-analysis critical value is slightly stricter, so the maximum number of events rises by a small percentage (often 2–5%, depending on the number of looks and spending function). Use a dedicated tool (for example the R package gsDesign) to compute it.

**How to report it:**
> Assuming exponential survival with a control median of 12 months, a hazard ratio of 0.70, two-sided α = 0.05, and 80% power, 247 events are required (Schoenfeld). With 24 months of uniform accrual and 12 months of additional follow-up, 370 patients (185 per arm) are expected to yield these events.

## Check your understanding

<details class="check"><summary>1. A sponsor says: "We enrolled 400 patients, so we have the power we planned." The protocol called for 247 deaths and only 150 have occurred. Is the sponsor right?</summary>

No. Power in a log-rank or Cox analysis is driven by the number of events. With 150 events instead of 247, the study is underpowered for HR = 0.70, no matter how many patients were enrolled. Either follow patients longer until the planned event count is reached, or accept a lower power.
</details>

<details class="check"><summary>2. The team is tempted to plan for HR = 0.60 because an early phase II study suggested it. The smallest worthwhile benefit is HR = 0.75. Which should drive the sample size, and why?</summary>

Plan for the smallest clinically worthwhile effect (HR = 0.75). Phase II estimates are usually optimistic and imprecise. Powering for HR = 0.60 needs only about 120 events, but if the true HR is closer to 0.75, the trial would be badly underpowered and could miss a benefit that matters to patients.
</details>

<details class="check"><summary>3. Median control survival is lengthened from 12 to 18 months with everything else unchanged. Does the required number of events change? Does the number of patients?</summary>

The number of events does not change: Schoenfeld's formula uses only α, power, allocation, and the HR. The number of patients **does** change, because with a longer median fewer patients die within the 36-month window, so P(event) falls and more patients must be enrolled to produce the same 247 events.
</details>

<details class="check"><summary>4. Calculation: with the same design, a statistician finds that 10% of patients will be lost to follow-up. How many patients should be enrolled?</summary>

Divide by the retention fraction: 370 / 0.90 = 411.1, so enrol 412 patients. Loss to follow-up reduces the events you will observe, so the inflation is needed to still reach the 247 events.
</details>

## Common mistakes

- **Reporting patients without events.** The protocol should state the **number of events** that triggers the final analysis, since that is what delivers the power.
- **Using the control P(event) alone.** The overall P(event) must average both arms.
- **Confusing the HR direction.** HR = 0.70 and HR = 1/0.70 = 1.43 give the same required events, but the control hazard (and so P(event)) must be assigned to the right arm.
- **Treating median ratio and HR as interchangeable.** They agree only under exponential survival. Under non-exponential models, convert the effect carefully.
- **Ignoring dropout and slow accrual.** Real trials usually take longer than the plan; build in a margin and monitor events, not calendar time.
- **Forgetting interim analyses.** Group-sequential boundaries inflate the maximum sample size slightly.
- **Rounding down.** Round events and patients **up**.

## Related tests

- [Sample size & power: means and proportions](#/test/sample-size-means-proportions): sample size for continuous and binary outcomes
- [Log-rank test](#/test/log-rank): the test that Schoenfeld's formula targets
- [Cox proportional hazards](#/test/cox-ph): the model used to estimate the HR
- [Kaplan–Meier](#/test/kaplan-meier): describing survival in each arm
- SAS users can use PROC POWER (`twosamplesurvival`) for comparable calculations.

## Key takeaways

- Power in a survival trial comes from **events**, not patients. State the event count that triggers the final analysis.
- Schoenfeld: d = (z₁₋α/₂ + z₁₋β)² / (p(1 − p)(log HR)²). For HR 0.70, 80% power, α = 0.05, 1:1 this is about 247 events.
- Patients = events / P(event). P(event) depends on the control event rate, the HR, accrual, and follow-up.
- Plan for the smallest worthwhile HR, round up, and then inflate for dropout and for interim analyses.
- The requirement grows with 1 / (log HR)², so modest HRs near 1 need very large trials.

## References

- Schoenfeld DA. The asymptotic properties of nonparametric tests for comparing survival distributions. *Biometrika.* 1981;68:316-319.
- Freedman LS. Tables of the number of patients required in clinical trials using the logrank test. *Stat Med.* 1982;1:121-129.
- ICH E9: *Statistical Principles for Clinical Trials.* 1998; Jennison C, Turnbull BW. *Group Sequential Methods with Applications to Clinical Trials.* Chapman & Hall/CRC; 2000.
- R documentation: `?integrate`; SciPy documentation: `scipy.integrate.quad`
