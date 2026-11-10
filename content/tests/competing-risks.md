# Competing risks (cumulative incidence & Fine–Gray)

> Estimates the probability of an event of interest (for example relapse) when another event (for example death without relapse) can remove the patient from risk, using the **cumulative incidence function** instead of 1 − Kaplan–Meier, and compares groups with cause-specific or Fine–Gray hazard models.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Time to first of several exclusive events; another event can prevent the one of interest | Time to event with event type and censoring | Independent censoring; correct event coding; proportional (sub)hazards | Relapse after transplant, with death without relapse as competitor |

## Learning objectives

- Explain why 1 − Kaplan–Meier overestimates risk when competing events are common
- Distinguish the cumulative incidence function, the cause-specific hazard, and the Fine–Gray subdistribution hazard
- Run these analyses in R and Python and read CIFs, Gray's test, and the two kinds of HR
- Choose and report the right estimand, and always report the competing event too

## What it is

In many studies a patient can experience more than one type of event, and the first one that happens prevents the others. After a stem-cell transplant a patient may **relapse** (the event of interest) or **die without relapse** (a *competing event*). Someone who has died can no longer relapse. The ordinary [Kaplan–Meier estimator](#/test/kaplan-meier) treats death as censoring, which quietly assumes that these patients would have gone on to relapse at the same rate as those still alive. That is not true, and the result is biased.

**Why 1 − KM overestimates.** Kaplan–Meier answers "what would the relapse risk be *in a world where nobody can die of anything else*?". Censoring a death moves its weight onto the patients who remain at risk, so the curve rises faster than reality. The risk of relapse in the real world is the **cumulative incidence function** (CIF), the probability of having had a relapse *as the first event* by time t. The CIFs of all event types plus the probability of being event-free always add up to 1, whereas 1 − KM for each cause taken separately can add up to more than 1. The more common the competing event, the bigger the gap.

- **Cumulative incidence (Aalen–Johansen / Kalbfleisch–Prentice):** CIFₖ(t) = Σ over event times ≤ t of S(t⁻) · dₖ(t) / n(t), where S is the overall (any-event) event-free survival, dₖ the number of events of type k and n(t) the number at risk.
- **Gray's test:** compares CIFs for one cause between groups (a competing-risks analogue of the log-rank test).
- **Cause-specific hazard (CSH):** the instantaneous rate of cause k among patients *still event-free*. Fitted with an ordinary [Cox model](#/test/cox-ph) that treats the other event as censored. Effect size: cause-specific HR.
- **Subdistribution hazard (Fine–Gray):** the instantaneous rate of cause k among patients who are event-free **or have already had a competing event**. Patients who died stay in the risk set. Effect size: subdistribution HR (SHR), which is tied directly to the CIF (it is not a hazard ratio on the usual risk set of patients still event-free).

**Which hazard to report.**

| Question | Use | Why |
|---|---|---|
| Etiology: does the treatment change the biological rate of relapse? | Cause-specific HR | Not influenced by how often patients die; compares like with like among those still at risk |
| Prognosis or absolute risk: what is the chance a patient relapses by 1 year? | CIF, with Fine–Gray SHR for the effect | The SHR is tied to the CIF, so it reflects both the relapse rate and the competing event |

A treatment can lower the relapse hazard and still *raise* the CIF of relapse if it keeps patients alive longer, or the reverse. Report the CSH for **both** events, plus the CIFs, to avoid misleading readers. SAS users can use PROC LIFETEST (CIF and Gray's test) and PROC PHREG with `EVENTCODE=` (cause-specific and Fine–Gray models).

## The intuition

Imagine a bathtub with two drains. Water leaves through the "relapse" drain or the "death without relapse" drain, and each patient leaves through whichever opens first. If you want the volume that left through the relapse drain, you must also count the water that was drained away by the other drain before it could reach the relapse drain. Kaplan–Meier pretends the second drain is plugged and the dead patients were only "lost to follow-up", so it overstates how much would have gone through the relapse drain. The cumulative incidence function (CIF) counts what actually happened.

The signal in a comparison depends on the question. For the biology (does the regimen change the *rate* of relapse among patients still at risk?), the signal is the cause-specific hazard ratio. For absolute risk (what fraction of patients relapse by 12 months?), the signal is the difference in CIFs, and the competing deaths are part of that story.

The two dot plots show the real times (months) in our 50-patient example. Look at the relapses first, then the deaths without relapse.

<div class="viz" data-viz="dots" data-unit="months to relapse" data-groups='{"Standard":[9.93,2.83,2.77,4.04,2.26,24.11,5.09,3.04],"New regimen":[8.61,10.21,6.46,4.18,4.28]}' data-caption="Relapses: fewer in the new regimen arm (5 vs 8), at broadly similar times, with just a handful of patients per arm."></div>

<div class="viz" data-viz="dots" data-unit="months to death without relapse" data-groups='{"Standard":[7.01,16.51,17.34,10.7,9.04,3.75],"New regimen":[13.69,7.12,0.07,0.41,2.95,9.91,20.8,0.11,7.08,1.29,8.56,9.28]}' data-caption="Deaths without relapse: twice as many in the new regimen arm (12 vs 6), many of them very early. Each of these patients can no longer relapse, which is why the new arm's relapse CIF is lower (20.9%) than the naive 1 - KM suggests (30.3%)."></div>

What to notice: the new arm has fewer relapses but many early deaths. A naive analysis makes the two arms look identical for relapse at 12 months (1 − KM 30.2% vs 30.3%), while the CIFs are 28.9% vs 20.9%. Gray's test for death without relapse gives χ² = 4.68 on 1 df. Under H₀ (equal CIFs of death in both arms), it follows the curve below.

<div class="viz" data-viz="chisq" data-stat="4.68" data-df="1"></div>

## When to use it

Use competing-risks methods when **all** of these are true:

- The outcome is a **time to the first of several mutually exclusive event types** (relapse vs death without relapse, cancer death vs non-cancer death, transplant vs death on a waiting list)
- The competing event **prevents** the event of interest (death) or changes its meaning, so censoring it is not harmless
- You want the **probability** of the event of interest over time, or a comparison of groups that respects the competing event

**Clinical scenario.** A single-centre study follows 50 patients after allogeneic transplant for 36 months: 25 receive standard conditioning and 25 a new conditioning regimen. *Does the new regimen change the risk of relapse, given that some patients die of other causes (infection, graft-versus-host disease) before they can relapse?* Status is 0 = censored, 1 = relapse, 2 = death without relapse.

**Don't use it when:**

- There is only **one** event type, or the competing event is negligible → use the [Kaplan–Meier estimator](#/test/kaplan-meier), the [log-rank test](#/test/log-rank) and [Cox regression](#/test/cox-ph)
- The competing event does **not** remove patients (for example you study overall survival, which counts death of any cause as the event) → standard survival methods are appropriate
- The same patient can have **repeated** events → use recurrent-event models (frailty or marginal models)
- The outcome is binary at a fixed time with no censoring → use [logistic regression](#/test/logistic-regression)

## Assumptions

1. **Independent censoring.** Patients censored at a given time have the same future risk of every event type as those still followed (non-informative censoring, as in [Kaplan–Meier](#/test/kaplan-meier)).
2. **Mutually exclusive, correctly classified events.** Each patient contributes the **first** event only. Check the coding: a patient who relapsed and later died is a relapse (1), not a death (2).
3. **Proportional (sub)distribution hazards** for the Fine–Gray model, and proportional cause-specific hazards for the Cox models. Check with Schoenfeld residuals (`cox.zph` for CSH) and by plotting the CIFs; for Fine–Gray, add a covariate × time interaction in `crr` (argument `cov2`, `tf`).
4. **The Fine–Gray model assumes the competing-event pattern in the study stays the same.** Because patients who died remain in the risk set, an SHR depends on the competing-event rate and does not transport to settings with different non-relapse mortality.
5. **Enough events of each type.** A traditional rough guide is about 10 events per covariate for the event of interest (a guide, not a threshold). This 50-patient example (13 relapses, 18 deaths) is a teaching dataset only.
6. **Independent patients,** one row per patient.

## Hypotheses

- **Gray's test, per cause:** H₀: the cumulative incidence of cause k is the same in the two arms at all times; H₁: it differs at some time.
- **Cause-specific Cox, per cause:** H₀: HR = 1 (equal cause-specific hazards); H₁: HR ≠ 1.
- **Fine–Gray:** H₀: SHR = 1 (equal subdistribution hazards, hence no proportional difference in the CIF); H₁: SHR ≠ 1.

## Example data

The data are simulated once in R (`set.seed(2026)`, exponential latent times for relapse and death with independent uniform censoring between 6 and 36 months, times rounded to 0.01 month) and typed in compactly in both languages: 25 patients per arm, in months from transplant. All times are distinct, so **there are no ties** (see the note on ties under Python).

| Status | Standard (n = 25) | New regimen (n = 25) |
|---|---|---|
| 0 = censored | 11 | 8 |
| 1 = relapse | 8 | 5 |
| 2 = death without relapse | 6 | 12 |

## How to run it

### R

```r
library(survival); library(cmprsk)
# Months from transplant; status 0 = censored, 1 = relapse, 2 = death without relapse
time <- c(9.93,2.83,6.81,7.01,2.77,34.03,9.38,4.04,17.36,7.79,16.51,20.34,2.26,30.06,24.11,17.34,10.7,12.71,9.04,32.4,3.75,35.76,21.75,5.09,3.04,
          33.89,8.61,15.52,13.69,31.1,10.29,7.02,10.21,7.12,9,11.06,0.07,0.41,2.95,10.22,9.91,20.8,0.11,7.08,1.29,6.46,8.56,4.18,9.28,4.28)
status <- c(1,1,0,2,1,0,0,1,0,0,2,0,1,0,1,2,2,0,2,0,2,0,0,1,1,
            0,1,0,2,0,0,0,1,2,0,0,2,2,2,0,2,2,2,2,2,1,2,1,2,1)
arm <- factor(rep(c("Standard", "New"), each = 25), levels = c("Standard", "New"))
d <- data.frame(time, status, arm)
table(d$arm, d$status)                # rows: arm; columns: status

# 1. Naive 1 - KM for relapse (death treated as censoring) vs the CIF at 12 months
km <- survfit(Surv(time, status == 1) ~ arm, data = d)
round(1 - summary(km, times = 12)$surv, 4)   # Standard, New
ci <- cuminc(d$time, d$status, group = d$arm)   # CIF for every status code except 0
timepoints(ci, times = c(12, 24))     # rows: arm + cause (1 = relapse, 2 = death)
round(ci$Tests, 4)                    # Gray's test, one row per cause

# 2. Fine-Gray subdistribution hazard for relapse
x <- model.matrix(~ arm, data = d)[, -1, drop = FALSE]
fg <- crr(d$time, d$status, cov1 = x, failcode = 1, cencode = 0)
summary(fg)

# 3. Cause-specific Cox models (the other event treated as censored)
round(summary(coxph(Surv(time, status == 1) ~ arm, data = d))$conf.int, 3)   # relapse
round(summary(coxph(Surv(time, status == 2) ~ arm, data = d))$conf.int, 3)   # death
```

**What the code does**

- `status` has **three** codes. `cuminc(ftime, fstatus, group)` treats 0 as censoring and estimates a separate CIF for every other code, per group. Row labels such as `Standard 1` mean "arm Standard, cause 1 (relapse)".
- `survfit(Surv(time, status == 1) ~ arm)` is the **naive** analysis: it recodes death as censored. `1 - summary(km, times = 12)$surv` gives the biased 1 − KM at 12 months, one value per arm.
- `timepoints(ci, times = ...)` reads the CIF (`$est`) and its variance (`$var`) at chosen times. The CIF is a step function that only moves at relapse times; no relapse occurs between months 12 and 24 in either arm (the last relapses are at 10.21 and 24.11 months), so the relapse values at 12 and 24 months are equal.
- `ci$Tests` is **Gray's test**, one row per cause. Row 1 is relapse, row 2 is death without relapse.
- `crr()` fits the **Fine–Gray** model. It takes a numeric covariate matrix (`cov1`), not a formula, so `model.matrix` builds the 0/1 arm indicator. `failcode = 1` is the event of interest and `cencode = 0` the censoring code. Patients with status 2 are kept in the risk set, which is what defines the subdistribution hazard.
- The two `coxph` calls are **cause-specific** models: `status == 1` makes deaths censored and vice versa. `conf.int` holds the HR and its 95% CI. R uses Efron's method for tied times, but there are no ties here.

**Output**

```text
            0  1  2
  Standard 11  8  6
  New       8  5 12
[1] 0.3022 0.3026
$est
                  12        24
Standard 1 0.2886516 0.2886516
New 1      0.2088636 0.2088636
Standard 2 0.1761810 0.2832145
New 2      0.4202273 0.5933182
    stat     pv df
1 1.1178 0.2904  1
2 4.6750 0.0306  1

         coef exp(coef) se(coef)     z p-value
armNew -0.599      0.55    0.548 -1.09    0.27

       exp(coef) exp(-coef) 2.5% 97.5%
armNew      0.55       1.82 0.188  1.61

       exp(coef) exp(-coef) lower .95 upper .95
armNew     0.735      1.361     0.239      2.26
       exp(coef) exp(-coef) lower .95 upper .95
armNew     2.796      0.358     1.031     7.579
```

### Python

```python
import numpy as np, pandas as pd
from lifelines import AalenJohansenFitter, KaplanMeierFitter, CoxPHFitter

# Months from transplant; status 0 = censored, 1 = relapse, 2 = death without relapse
time = [9.93,2.83,6.81,7.01,2.77,34.03,9.38,4.04,17.36,7.79,16.51,20.34,2.26,30.06,24.11,17.34,10.7,12.71,9.04,32.4,3.75,35.76,21.75,5.09,3.04,
        33.89,8.61,15.52,13.69,31.1,10.29,7.02,10.21,7.12,9,11.06,0.07,0.41,2.95,10.22,9.91,20.8,0.11,7.08,1.29,6.46,8.56,4.18,9.28,4.28]
status = [1,1,0,2,1,0,0,1,0,0,2,0,1,0,1,2,2,0,2,0,2,0,0,1,1,
          0,1,0,2,0,0,0,1,2,0,0,2,2,2,0,2,2,2,2,2,1,2,1,2,1]
d = pd.DataFrame({"time": time, "status": status, "new": np.repeat([0, 1], 25)})

def at(curve, t):                       # value of a step function at time t
    return float(curve[curve.index <= t].iloc[-1])

for arm, g in d.groupby("new"):
    name = ["Standard", "New"][arm]
    km = KaplanMeierFitter().fit(g.time, g.status == 1)          # death = censored
    ajf = AalenJohansenFitter(seed=1).fit(g.time, g.status, event_of_interest=1)
    ajd = AalenJohansenFitter(seed=1).fit(g.time, g.status, event_of_interest=2)
    print(f"{name:9s} 12 mo: 1-KM relapse = {1 - km.predict(12):.4f}, "
          f"CIF relapse = {at(ajf.cumulative_density_.iloc[:, 0], 12):.4f}, "
          f"CIF death = {at(ajd.cumulative_density_.iloc[:, 0], 12):.4f}")

# Cause-specific Cox HRs (the other event is treated as censored)
for ev, label in [(1, "relapse"), (2, "death w/o relapse")]:
    dd = d.assign(event=(d.status == ev).astype(int))[["time", "event", "new"]]
    cph = CoxPHFitter().fit(dd, "time", "event")
    s = cph.summary.loc["new"]
    print(f"Cause-specific HR, {label}: {s['exp(coef)']:.3f} "
          f"(95% CI {s['exp(coef) lower 95%']:.3f} to {s['exp(coef) upper 95%']:.3f}), p = {s['p']:.3f}")
```

**What the code does**

- `AalenJohansenFitter` estimates the CIF for **one** event type at a time (`event_of_interest`); the other codes are treated as competing events and 0 as censoring. Call it once per cause and per arm. `cumulative_density_` is a step-function table indexed by time, and `at()` reads its value at a chosen time.
- `KaplanMeierFitter().fit(g.time, g.status == 1)` gives the **naive** curve with death as censoring, for the 1 − KM comparison.
- **Ties.** The Aalen–Johansen estimator in lifelines does not accept tied event times: it adds a tiny random jitter to break them (the `seed` argument makes this reproducible, and `jitter_level` controls its size) and warns when it does so. Our data have **no tied times**, so no jitter is applied and the results equal R's `cuminc` to four decimals. With your own tied data, expect differences in the later decimals, and always set a seed and say so.
- The cause-specific HRs use `CoxPHFitter` on a 0/1 `event` column, so the other event is censored. lifelines uses Efron's method for ties, the same as R's `coxph` (SAS PHREG defaults to Breslow).
- **No Fine–Gray in mainstream Python.** lifelines, statsmodels and scikit-survival do not provide a Fine–Gray regression, and the few small third-party packages are not widely validated. For the subdistribution hazard (and for Gray's test), use R `cmprsk::crr` / `cuminc` (or the `tidycmprsk` and `riskRegression` packages), or SAS PROC PHREG. A pragmatic Python workflow is to report CIFs plus cause-specific HRs for **both** events.

**Output**

```text
Standard  12 mo: 1-KM relapse = 0.3022, CIF relapse = 0.2887, CIF death = 0.1762
New       12 mo: 1-KM relapse = 0.3026, CIF relapse = 0.2089, CIF death = 0.4202
Cause-specific HR, relapse: 0.735 (95% CI 0.239 to 2.260), p = 0.591
Cause-specific HR, death w/o relapse: 2.796 (95% CI 1.031 to 7.579), p = 0.043
```

## Interpreting the output

- **Counts.** Standard: 8 relapses, 6 deaths without relapse, 11 censored. New regimen: 5 relapses, 12 deaths without relapse, 8 censored. The new arm has fewer relapses but twice as many competing deaths.
- **1 − KM versus the CIF at 12 months (R and Python agree).** Standard: 1 − KM = 30.2%, CIF = 28.9%, a small gap because few patients died first. New regimen: 1 − KM = 30.3%, CIF = **20.9%**, so the naive estimate is about 9 percentage points (45%) too high, because 12 of 25 patients died and could not relapse. The naive analysis would make the arms look identical (30.2% vs 30.3%), whereas the correct CIFs are 28.9% vs 20.9%.
- **CIF of death without relapse.** At 12 months 17.6% (standard) vs 42.0% (new); at 24 months 28.3% vs 59.3%. Gray's test: relapse p = 0.29 (no clear difference in relapse incidence), death without relapse p = 0.031.
- **Cause-specific HR for relapse = 0.735 (95% CI 0.239 to 2.260, p = 0.59).** Among patients still alive and relapse-free, the relapse rate with the new regimen is about 27% lower, but the CI is very wide and includes 1. **Cause-specific HR for death without relapse = 2.80 (1.03 to 7.58, p = 0.043):** the new regimen is associated with a higher non-relapse mortality rate. R and Python give identical values.
- **Fine–Gray SHR for relapse = 0.55 (95% CI 0.19 to 1.61, p = 0.27).** This is *further* from 1 than the cause-specific HR because many patients in the new arm die first, which lowers their chance of ever relapsing even if the relapse rate itself is only modestly reduced. The SHR is consistent with the lower CIF curve, but it is not an effect on the biology of relapse. Python has no equivalent, so this number comes from R only.
- **Both CIs include 1 for relapse.** With 13 relapses the study is far too small to conclude anything about relapse. Treat the example as an illustration of the methods.

**How to report it:**
> The 12-month cumulative incidence of relapse was 28.9% with standard conditioning and 20.9% with the new regimen (Gray's test p = 0.29), with death without relapse as a competing event (cumulative incidence 17.6% vs 42.0%; p = 0.031). The Fine–Gray subdistribution hazard ratio for relapse was 0.55 (95% CI 0.19 to 1.61); the cause-specific hazard ratios were 0.74 (0.24 to 2.26) for relapse and 2.80 (1.03 to 7.58) for death without relapse.

⚠️ **Important for trials:** report the competing event too. A regimen that appears to reduce relapse while increasing non-relapse mortality is not necessarily beneficial, and pre-specify which hazard (cause-specific or subdistribution) is the primary estimand (see ICH E9(R1) on estimands and intercurrent events such as death).

## Check your understanding

<details class="check"><summary>1. In a study of older adults, the outcome is dementia diagnosis and 30% of participants die of other causes before they can be diagnosed. The analyst plots 1 − KM for dementia, censoring deaths. What is wrong?</summary>

1 − KM answers "what would the dementia risk be if nobody could die of anything else?" and overstates the real-world risk, because deaths are treated as if those people would go on to be diagnosed at the same rate as survivors. Use the cumulative incidence function (Aalen–Johansen) for dementia, with death as a competing event.
</details>

<details class="check"><summary>2. A new drug reduces the cause-specific hazard of relapse (HR = 0.7) but increases deaths from infection. Could the cumulative incidence of relapse still go up or stay the same? And which HR is the one for the question "does the drug act on relapse biology?"</summary>

The CIF of relapse can go up, stay the same, or go down depending on how many patients are kept alive (and so at risk) or die first. The cause-specific HR (here 0.7) is the one for the biological question, because it compares patients still at risk. The Fine–Gray SHR and the CIF answer the absolute-risk question and must be read together with the competing event. Report the cause-specific HR for both events plus the CIFs.
</details>

<details class="check"><summary>3. In the example the Fine–Gray SHR for relapse (0.55) is further from 1 than the cause-specific HR (0.735). A colleague says the new regimen "cuts relapse hazard by 45%." Is that a fair reading?</summary>

No. The SHR keeps patients who died in the risk set, so it reflects the competing deaths as well as the relapse rate. Much of its distance from 1 comes from the 12 early deaths in the new arm, not from a biological effect on relapse. The cause-specific HR of 0.735 (95% CI 0.24 to 2.26) is the better estimate of the effect on the relapse rate, and both CIs include 1, so neither supports a claim of benefit.
</details>

<details class="check"><summary>4. What would happen to the gap between 1 − KM and the CIF if the competing death were rare (for example 2% of patients)?</summary>

The gap would shrink to almost nothing. 1 − KM and the CIF differ by the weight of the competing deaths, so with rare competing events the ordinary Kaplan–Meier analysis is a good approximation. The standard arm in the example (6 of 25 deaths, small gap 30.2% vs 28.9%) shows this.
</details>

## Common mistakes

- **Using 1 − KM (or the log-rank test) for the event of interest when death is common.** This overestimates the absolute risk. Use the CIF.
- **Reporting only the Fine–Gray SHR as if it were a biological effect.** It mixes the event rate with the competing-event rate. For etiology use the cause-specific HR.
- **Reporting only one cause.** Always show the competing event as well (its CIF and HR).
- **Treating the sum of the two cause-specific 1 − KM curves as meaningful.** They can add up to more than 100%, while the CIFs plus event-free survival add up to exactly 100%.
- **Miscoding the status.** Only the first event counts; passing a 0/1 status to `cuminc` or `AalenJohansenFitter` silently gives the naive analysis.
- **Ignoring ties in lifelines.** The Aalen–Johansen fitter adds random jitter to tied times, so set `seed` and report it.
- **Interpreting the SHR as a hazard among patients at risk.** Patients who died remain in its risk set, which is why its interpretation is less direct than that of an ordinary HR.
- **Forgetting to check proportional hazards** for the Cox and Fine–Gray models, and fitting too many covariates for 13 relapses.

## Related tests

- [Kaplan–Meier estimator](#/test/kaplan-meier): the single-event estimator that competing risks corrects
- [Log-rank test](#/test/log-rank): the cause-specific test analogue of Gray's test
- [Cox proportional hazards regression](#/test/cox-ph): used for the cause-specific hazard ratio
- [Sample size for time-to-event studies](#/test/sample-size-survival): planning the number of events
- [Logistic regression](#/test/logistic-regression): when only a fixed-time binary outcome is analysed

## Key takeaways

- When another event can prevent the event of interest, 1 − KM overestimates risk. Use the cumulative incidence function (CIF).
- Cause-specific HR answers "what does the treatment do to the rate among patients still at risk?" (etiology). The Fine–Gray SHR is tied to the CIF (absolute risk).
- A treatment can lower the cause-specific hazard of an event yet fail to lower its CIF. Report both events.
- Gray's test compares CIFs. It is the competing-risks analogue of the log-rank test.
- Code every patient with their first event only, and pre-specify the estimand.

## References

- Fine JP, Gray RJ. A proportional hazards model for the subdistribution of a competing risk. *J Am Stat Assoc.* 1999;94:496–509.
- Austin PC, Lee DS, Fine JP. Introduction to the analysis of survival data in the presence of competing risks. *Circulation.* 2016;133:601–609.
- Gray RJ. A class of K-sample tests for comparing the cumulative incidence of a competing risk. *Ann Stat.* 1988;16:1141–1154.
- R documentation: `cmprsk::cuminc`, `cmprsk::crr`, `survival::coxph`; lifelines: `AalenJohansenFitter`, `CoxPHFitter`
- SAS/STAT User's Guide: The LIFETEST Procedure (cumulative incidence) and the PHREG Procedure (`EVENTCODE=` option)
