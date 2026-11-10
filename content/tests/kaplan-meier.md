# Kaplan–Meier estimator

> Estimates the survival curve (the probability of still being event-free over time) from time-to-event data in which some patients are **censored**, and gives medians and landmark survival rates with confidence intervals.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Describe time to an event when some patients are still event-free (censored) | Time-to-event with censoring | Non-informative censoring; common time origin; independent patients | Progression-free survival in two trial arms |

## Learning objectives

- Explain **censoring** and why simply discarding censored patients biases results
- Describe how the KM curve is built as a running product of survival fractions
- Run the estimator in R and Python and read medians and landmark survival with CIs
- Report KM results and avoid mistakes like comparing curves by CI overlap

## What it is

Survival analysis deals with outcomes of the form "how long until something happens": time to death, disease progression, relapse, or hospital discharge. The difficulty is **censoring**. A patient who is still alive and progression-free when the study ends, who moves away, or who withdraws consent has not had the event yet. We only know that their event time is *longer than* the time they were observed. Throwing these patients away wastes information and biases the result, and treating their last visit as an event is simply wrong.

The Kaplan–Meier (KM) estimator uses everyone for as long as they are under observation. At each time an event occurs, it asks: *of the patients still at risk just before this moment, what fraction survived it?* The survival curve is the running product of those fractions.

- **Estimator:** Ŝ(t) = ∏ (1 − dᵢ / nᵢ) over all event times tᵢ ≤ t, where dᵢ is the number of events at tᵢ and nᵢ the number of patients at risk just before tᵢ. Censored patients leave the risk set without causing a drop in the curve.
- **Variance:** Greenwood's formula, which is used to build pointwise confidence intervals.
- **What to report:** the **median survival time** (the first time Ŝ(t) ≤ 0.5) with 95% CI, the **survival probability at a clinically meaningful landmark** (e.g., 12-month progression-free survival) with 95% CI, and the number at risk beneath the plot. The KM curve is descriptive. It does not by itself test or quantify a difference between arms (see [Log-rank test](#/test/log-rank) and [Cox regression](#/test/cox-ph)).

## The intuition

The KM estimator has no test statistic, so there is no "signal vs noise" here. The idea is **a running tally of who is still in the race**. At every moment an event occurs, we ask: of the people still at risk just before, what fraction got through? We multiply those fractions together as time passes.

Everyday analogy: a marathon where some runners stop to leave (censored) and others drop out injured (events). To estimate the chance of finishing a given mile, you only consider the runners still on the course at that point. Runners who leave the course early still count up to the moment they leave, but they are not counted as failures.

Here are the first three steps for the Control arm in our data:

| Month | At risk just before | Events | Fraction surviving | Cumulative S(t) |
|---|---|---|---|---|
| 1 | 18 | 2 | 16/18 | 0.889 |
| 3 | 16 | 1 | 15/16 | 0.833 |
| 4 | 15 | 1 | 14/15 | 0.778 |

The dot plot below shows every patient's follow-up time, split by arm and by whether they had an event or were censored. Notice how many Treatment patients are censored at long follow-up times (right-hand side), which is why the Treatment median is only reached at month 25 and is imprecise.

<div class="viz" data-viz="dots" data-unit="months from randomisation" data-groups='{"Control event":[3,13,4,19,10,12,1,6,1,7],"Control censored":[9,17,12,20,17,16,10,6],"Treatment event":[5,25,20,19,15,16,7],"Treatment censored":[18,10,9,24,26,25,7,15,27,26,16]}' data-caption="Events happen earlier in Control. Censored patients (no event yet) contribute information up to their last visit."></div>

Discarding the censored patients would throw away 19 of 36 patients and badly underestimate how long people stay progression-free. KM uses everyone for as long as they are observed.

## When to use it

Use the Kaplan–Meier estimator when **all** of these are true:

- The outcome is a **time from a clear origin** (randomisation, diagnosis, surgery) **to an event** (death, progression, relapse)
- Some patients have **not had the event** by the end of follow-up (right censoring)
- Censoring is **unrelated to prognosis** (non-informative censoring)
- You want a descriptive summary, overall or by group, rather than adjusted effect estimates

**Clinical scenario.** A two-arm trial randomises 36 patients with advanced cancer to standard care (Control, n = 18) or a new agent (Treatment, n = 18). The outcome is time to disease progression or death, in months. Several patients are still progression-free when the data are cut off. *What is the median progression-free survival in each arm, and what proportion of patients are progression-free at 12 months?* The same data set is used on the [Log-rank test](#/test/log-rank) and [Cox regression](#/test/cox-ph) pages, so you can follow the progression from description to comparison to adjusted effect.

**Don't use it when:**

- You want to **test** whether two curves differ → use the [log-rank test](#/test/log-rank)
- You need a **hazard ratio** or must adjust for covariates such as age → use [Cox proportional hazards regression](#/test/cox-ph)
- The outcome is a binary status at a fixed time with **no censoring** (e.g., 30-day mortality, everyone followed for 30 days) → a [proportion](#/test/one-proportion) or [two-proportion z-test](#/test/two-proportion-z) is simpler
- The outcome is continuous (e.g., change in blood pressure) → use the [t-tests](#/test/independent-t) or [linear regression](#/test/linear-regression)

## Assumptions

1. **Censoring is non-informative.** Patients who are censored must have the same future risk as those who remain under observation. This fails if, for example, patients drop out because they are getting sicker. Check by comparing baseline characteristics of censored and uncensored patients, reporting reasons for censoring, and running sensitivity analyses (e.g., treating dropouts as events).
2. **Event times are measured from a common, well-defined origin.** Use randomisation, not treatment start, in randomised trials. Check the protocol definition of time zero.
3. **Observations are independent.** One row per patient, with no repeated events per patient in the same curve. For recurrent events use other methods.
4. **Event ascertainment is the same for all patients.** If one arm is scanned every 6 weeks and the other every 12 weeks, progression times are shifted (interval censoring). Check the visit schedules.
5. **Follow-up is adequate and similar across groups.** Estimates in the tail rest on very few patients. Check the **number at risk** table and be wary of the curve where fewer than about 10 patients remain.

## Hypotheses

The Kaplan–Meier estimator is **descriptive**, so there is no H₀ or H₁. It estimates the survival function S(t) = P(T > t), the probability that the event time T is later than t, for each arm, together with derived quantities (median survival, S(12 months)) and their 95% confidence intervals. Testing whether S(t) is the same in the two arms is the job of the [log-rank test](#/test/log-rank).

## Example data

Thirty-six patients; time is months from randomisation to progression or death. A **+** marks a **censored** time (progression-free at last follow-up). Of the 36 patients, 17 had an event and 19 were censored.

**Control**

| Patient | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Months | 3 | 13 | 9+ | 17+ | 4 | 19 | 10 | 12+ | 20+ | 1 | 12 | 17+ | 6 | 1 | 7 | 16+ | 10+ | 6+ |
| Age | 70 | 66 | 61 | 60 | 71 | 72 | 66 | 61 | 52 | 53 | 49 | 66 | 69 | 64 | 36 | 61 | 55 | 66 |

**Treatment**

| Patient | 19 | 20 | 21 | 22 | 23 | 24 | 25 | 26 | 27 | 28 | 29 | 30 | 31 | 32 | 33 | 34 | 35 | 36 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Months | 18+ | 10+ | 9+ | 24+ | 26+ | 5 | 25 | 25+ | 7+ | 20 | 19 | 15+ | 27+ | 15 | 16 | 26+ | 16+ | 7 |
| Age | 49 | 46 | 71 | 65 | 70 | 60 | 65 | 50 | 61 | 67 | 46 | 61 | 50 | 59 | 58 | 72 | 58 | 65 |

Patient 3 (Control, 9+) was progression-free when last seen at month 9. They count in the risk set up to month 9 and then drop out without being counted as an event.

## How to run it

### R

```r
library(survival)
# Shared dataset: 36 patients, time in months, event 1 = progression/death, 0 = censored
trial <- data.frame(
  arm   = factor(rep(c("Control", "Treatment"), each = 18),
                 levels = c("Control", "Treatment")),
  time  = c(3,13,9,17,4,19,10,12,20,1,12,17,6,1,7,16,10,6,
            18,10,9,24,26,5,25,25,7,20,19,15,27,15,16,26,16,7),
  event = c(1,1,0,0,1,1,1,0,0,1,1,0,1,1,1,0,0,0,
            0,0,0,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1),
  age   = c(70,66,61,60,71,72,66,61,52,53,49,66,69,64,36,61,55,66,
            49,46,71,65,70,60,65,50,61,67,46,61,50,59,58,72,58,65)
)
fit_ll <- survfit(Surv(time, event) ~ arm, data = trial, conf.type = "log-log")
print(fit_ll)
summary(fit_ll, times = 12)
fit_log <- survfit(Surv(time, event) ~ arm, data = trial)
print(fit_log)
```

**What the code does**

- `Surv(time, event)` builds the survival object. `event` must be 1 for an event and 0 for censored. R prints censored times with a `+`.
- `survfit(... ~ arm)` fits one KM curve per arm. Use `~ 1` for a single overall curve.
- `conf.type = "log-log"` sets the confidence-interval transformation. **R's default is `"log"`, whereas SAS PROC LIFETEST defaults to log-log** (and lifelines also builds its intervals on the log(−log) scale). The log-log interval always stays inside (0, 1) and behaves better in the tails, so I set it explicitly here so that all three languages agree. The last block shows that the default can change the interval: the lower limits of the median CI move from 6 and 16 to 7 and 19.
- `print(fit)` shows the number of events and the **median with 95% CI**. `NA` for the upper limit means the confidence band never fell below 0.5 within follow-up, so that limit is "not reached".
- `summary(fit, times = 12)` gives the survival probability at 12 months with its standard error and CI. Note that `n.event` here is the cumulative number of events up to that time.

**Output**

```text
Call: survfit(formula = Surv(time, event) ~ arm, data = trial, conf.type = "log-log")
               n events median 0.95LCL 0.95UCL
arm=Control   18     10     13       6      NA
arm=Treatment 18      7     25      16      NA
Call: survfit(formula = Surv(time, event) ~ arm, data = trial, conf.type = "log-log")
                arm=Control 
        time       n.risk      n.event     survival      std.err lower 95% CI 
      12.000        8.000        8.000        0.521        0.125        0.261 
upper 95% CI 
       0.729 
                arm=Treatment 
        time       n.risk      n.event     survival      std.err lower 95% CI 
     12.0000      13.0000       2.0000       0.8889       0.0741       0.6242 
upper 95% CI 
      0.9710 
Call: survfit(formula = Surv(time, event) ~ arm, data = trial)
               n events median 0.95LCL 0.95UCL
arm=Control   18     10     13       7      NA
arm=Treatment 18      7     25      19      NA
```

### Python

```python
import pandas as pd

# Shared dataset: 36 patients, time in months, event 1 = progression/death, 0 = censored
trial = pd.DataFrame({
    "arm": ["Control"] * 18 + ["Treatment"] * 18,
    "time": [3,13,9,17,4,19,10,12,20,1,12,17,6,1,7,16,10,6,
             18,10,9,24,26,5,25,25,7,20,19,15,27,15,16,26,16,7],
    "event": [1,1,0,0,1,1,1,0,0,1,1,0,1,1,1,0,0,0,
              0,0,0,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1],
    "age": [70,66,61,60,71,72,66,61,52,53,49,66,69,64,36,61,55,66,
            49,46,71,65,70,60,65,50,61,67,46,61,50,59,58,72,58,65],
})
from lifelines import KaplanMeierFitter
from lifelines.utils import median_survival_times

for arm, g in trial.groupby("arm"):
    kmf = KaplanMeierFitter(alpha=0.05)
    kmf.fit(g["time"], event_observed=g["event"], label=arm)
    lo, hi = median_survival_times(kmf.confidence_interval_).iloc[0]
    band = kmf.confidence_interval_survival_function_
    s12 = kmf.survival_function_at_times(12).iloc[0]
    l12, u12 = band[band.index <= 12].iloc[-1]
    print(f"{arm}: median = {kmf.median_survival_time_}, 95% CI ({lo}, {hi})")
    print(f"   S(12) = {s12:.3f}, 95% CI ({l12:.3f}, {u12:.3f})")
```

**What the code does**

- `KaplanMeierFitter` is from the `lifelines` package. `fit(durations, event_observed)` takes the event indicator as 1 = event, 0 = censored.
- lifelines computes confidence intervals with the **exponential Greenwood formula**, which is a log(−log) interval. This is the same interval as `conf.type = "log-log"` in R and the SAS default, so no extra argument is needed. `alpha=0.05` gives a 95% CI.
- `median_survival_times(kmf.confidence_interval_)` finds where the lower and upper confidence bands cross 0.5. `inf` means "not reached", the counterpart of `NA` in R.
- `survival_function_at_times(12)` returns Ŝ(12). The CI at 12 months is read from the last row of the confidence band at or before month 12, because KM curves are step functions.

**Output**

```text
Control: median = 13.0, 95% CI (6.0, inf)
   S(12) = 0.521, 95% CI (0.261, 0.729)
Treatment: median = 25.0, 95% CI (16.0, inf)
   S(12) = 0.889, 95% CI (0.624, 0.971)
```

### SAS

```sas
data trial;
  input id arm $ months event age;   /* event: 1 = event, 0 = censored */
  datalines;
1 Control 3 1 70
2 Control 13 1 66
3 Control 9 0 61
4 Control 17 0 60
5 Control 4 1 71
6 Control 19 1 72
7 Control 10 1 66
8 Control 12 0 61
9 Control 20 0 52
10 Control 1 1 53
11 Control 12 1 49
12 Control 17 0 66
13 Control 6 1 69
14 Control 1 1 64
15 Control 7 1 36
16 Control 16 0 61
17 Control 10 0 55
18 Control 6 0 66
19 Treatment 18 0 49
20 Treatment 10 0 46
21 Treatment 9 0 71
22 Treatment 24 0 65
23 Treatment 26 0 70
24 Treatment 5 1 60
25 Treatment 25 1 65
26 Treatment 25 0 50
27 Treatment 7 0 61
28 Treatment 20 1 67
29 Treatment 19 1 46
30 Treatment 15 0 61
31 Treatment 27 0 50
32 Treatment 15 1 59
33 Treatment 16 1 58
34 Treatment 26 0 72
35 Treatment 16 0 58
36 Treatment 7 1 65
;
run;

ods graphics on;
proc lifetest data=trial method=km conftype=loglog alpha=0.05
              timelist=6 12 18 24 reduceout plots=survival(atrisk cl);
  time months*event(0);              /* (0) = value that means censored */
  strata arm;
run;
```

**What the code does**

- `time months*event(0);` names the follow-up time and the event indicator. The number in parentheses is the value that marks a **censored** observation, which is easy to get backwards. I named the variable `months` instead of `time` to avoid confusion with the TIME statement.
- `strata arm;` produces one KM curve per arm. Without a STRATA statement you get one overall curve.
- `method=km` requests the Kaplan–Meier (product-limit) estimator, which is the default. `method=lt` would give an actuarial life table.
- `conftype=loglog` is the **SAS default** and is written out for clarity. Other choices are `log`, `linear`, `asinsqrt`, and `logit`. Use `conftype=log` to reproduce R's default.
- `timelist=6 12 18 24` prints survival estimates with 95% limits at those landmark times, and `reduceout` limits the printed product-limit table to those times.
- The "Quartile Estimates" table gives the median (50th percentile) with its 95% CI for each arm.
- `plots=survival(atrisk cl)` draws the KM curves with a number-at-risk table and confidence bands. The output is not shown here because SAS cannot be run in this environment. With these options SAS should report the same estimates as R and Python.

## Interpreting the output

- **Events and censoring.** Control: 10 events among 18 patients (8 censored). Treatment: 7 events among 18 patients (11 censored). Both languages agree.
- **Median survival, Control = 13 months (95% CI 6 to not reached).** Half of Control patients are estimated to have progressed or died by month 13. The interval is wide (the upper limit is "not reached") because there are only 10 events.
- **Median survival, Treatment = 25 months (95% CI 16 to not reached).** The curve only drops below 50% at month 25, when few patients remain at risk, so this median is imprecise. Do not over-read it.
- **Survival at 12 months, Control = 52.1% (95% CI 26.1% to 72.9%)**, with 8 patients still at risk.
- **Survival at 12 months, Treatment = 88.9% (95% CI 62.4% to 97.1%)**, with 13 patients at risk.
- **Overlapping intervals are not a test.** The CIs at 12 months only partly overlap, but whether the arms differ is answered by the [log-rank test](#/test/log-rank), and the size of the difference by the hazard ratio from [Cox regression](#/test/cox-ph).
- **Same results in R and Python.** Medians, survival at 12 months and all CI limits match to three decimals once the CI type is the same (log-log).

**How to report it:**
> Median progression-free survival was 13 months (95% CI 6 to not reached) in the control arm and 25 months (95% CI 16 to not reached) in the treatment arm. Estimated 12-month progression-free survival was 52.1% (95% CI 26.1 to 72.9) versus 88.9% (95% CI 62.4 to 97.1). Kaplan–Meier estimates with log-log confidence intervals; 19 of 36 patients were censored.

Also report the **median follow-up** (preferably computed by the reverse Kaplan–Meier method, flipping the event indicator) and show the number at risk under the plot.

## Check your understanding

<details class="check"><summary>1. A patient is alive and progression-free at their last visit at month 9, then moves away. How is this patient handled?</summary>

As **censored** at month 9. They count in the risk set up to month 9 and then leave it without causing a drop in the curve. They are neither an event nor ignored.
</details>

<details class="check"><summary>2. The median survival in the treatment arm is 25 months with 95% CI 16 to "not reached." Is this a sign something went wrong?</summary>

No. "Not reached" means the confidence band did not fall below 50% during follow-up. With only 7 events and few patients at risk late in follow-up, the upper limit is unbounded. It signals imprecision, not an error.
</details>

<details class="check"><summary>3. The 12-month survival CIs for the two arms overlap slightly (26% to 73% vs 62% to 97%). Can you conclude there is no difference?</summary>

No. Overlapping CIs are not a test, and survival at one time point does not capture the whole curve. Use the pre-specified [log-rank test](#/test/log-rank) and report the hazard ratio from [Cox regression](#/test/cox-ph).
</details>

<details class="check"><summary>4. What if the censored Control patients had actually dropped out because they felt much worse (informative censoring)? How would the KM curve be affected?</summary>

The curve would look too optimistic. The patients who left were at higher risk, so the remaining risk set is healthier than the full group. Report reasons for censoring and do sensitivity analyses, for example treating dropouts as events.
</details>

## Common mistakes

- **Reporting the mean of observed times or a simple proportion of events.** Censoring makes both biased. Use the KM median or survival at a landmark time.
- **Coding the event indicator backwards.** If `event` = 1 means censored, the curve will show the wrong thing. In SAS, the value in `event(0)` must be the *censoring* value.
- **Assuming the CI type is the same across software.** R defaults to `log`, SAS to `loglog`, and lifelines to an exponential-Greenwood (log-log) interval. Results can differ in the lower limit, as seen above. State the CI type in your methods and set it explicitly.
- **Trusting the tail of the curve.** With only a few patients at risk, one event creates a big drop. Show the number at risk and truncate the plot where it is small.
- **Reading a median that is "not reached" as "no events".** It means the curve did not fall to 50% (or the CI did not cross 50%) during follow-up.
- **Informative censoring.** If sicker patients drop out, the KM curve is too optimistic. Report reasons for censoring and do a sensitivity analysis.
- **Comparing curves visually or by CI overlap alone.** Always pair KM curves with a pre-specified test and effect estimate (ICH E9).

## Related tests

- [Log-rank test](#/test/log-rank): tests whether two or more KM curves differ
- [Cox proportional hazards regression](#/test/cox-ph): hazard ratio and covariate adjustment
- [Cochran–Mantel–Haenszel test](#/test/cmh): the stratified analysis idea for binary outcomes
- [Two-proportion z-test](#/test/two-proportion-z): compare event proportions at a fixed time when there is no censoring
- [Poisson & negative binomial regression](#/test/poisson-regression): event rates per person-time

## Key takeaways

- Kaplan–Meier estimates S(t) as a running product of (1 − events/at-risk) at each event time.
- Censored patients count while observed and then leave the risk set. Never drop them or call them events.
- Report median survival and landmark survival (for example 12-month) with 95% CIs, plus the number at risk.
- "Median not reached" means the curve or its CI never crossed 50%.
- KM is descriptive. Compare arms with the log-rank test or Cox model, not by CI overlap.

## References

- Kaplan EL, Meier P. Nonparametric estimation from incomplete observations. *J Am Stat Assoc.* 1958;53:457–481.
- Bland JM, Altman DG. Survival probabilities (the Kaplan–Meier method). *BMJ.* 1998;317:1572.
- Collett D. *Modelling Survival Data in Medical Research.* 3rd ed. CRC Press; 2015. Ch. 2.
- R documentation: `?survfit`, `?survfit.formula` (argument `conf.type`); lifelines documentation: `KaplanMeierFitter`; SAS/STAT User's Guide: The LIFETEST Procedure (CONFTYPE= option)
