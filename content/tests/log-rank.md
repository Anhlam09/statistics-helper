# Log-rank test

> Tests whether two or more survival curves are the same by comparing, at every event time, the observed number of events in each group with the number expected if the groups had identical hazards.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Compare survival curves between two or more groups | Time-to-event with right censoring | Non-informative censoring; roughly proportional hazards; independent patients | Progression-free survival, treatment vs control |

## Learning objectives

- Explain the log-rank test as observed vs expected events at every event time
- Run it in R and Python and read observed, expected, χ², and p
- Recognize when crossing curves make the standard log-rank test a poor choice
- Pair the test with a hazard ratio and CI instead of reporting p alone

## What it is

The log-rank test is the standard way to compare time-to-event outcomes between groups when some patients are censored (see [Kaplan–Meier estimator](#/test/kaplan-meier) for what censoring means). At every time an event occurs, the test looks at the patients still at risk and asks how many events *each group* would be expected to have if the event was equally likely in all of them. It then adds up the observed minus expected events over all event times. If one arm consistently has more events than expected, the groups differ.

- **Test statistic:** χ² = (O₁ − E₁)² / V, where O₁ and E₁ are the observed and expected events in group 1 and V is the hypergeometric variance summed over the event times. For two groups it follows a χ² distribution with 1 df. For k groups the df is k − 1.
- **Effect size:** the log-rank test gives a p-value, not an effect size. Pair it with the **hazard ratio** and 95% CI from [Cox regression](#/test/cox-ph) and with the KM medians or landmark survival rates. A rough hazard ratio can be read off as (O₁/E₁) / (O₂/E₂).
- **Weighting:** the standard log-rank test gives every event time the same weight, which is most powerful when the hazard ratio is constant over time. Other weights (Wilcoxon/Gehan, Peto–Peto, Fleming–Harrington) emphasise early events and should be pre-specified.

## The intuition

Picture the trial as a series of small races. Each time a patient progresses, look at everyone still being followed at that moment, from both arms. If the arms were identical, the person who progressed would be a random draw from that group: an arm making up 60% of those at risk should account for about 60% of the events. The log-rank test tallies, over all event times, how many events each arm actually had versus how many that "fair draw" predicts.

- **Signal:** the total gap O − E. Control had 10 events where 5.87 were expected, a gap of about +4.1.
- **Noise:** the random variation in who gets drawn, measured by the hypergeometric variance V (summed over event times).
- **Censored patients** do not count as events, but they stay in the "at risk" pool until they are censored, which is how the test uses their information.

Look at the follow-up times below (events and censored times mixed, as the dot plot cannot show censoring). Treatment patients tend to be followed longer, but the overlap shows why we need the censoring-aware test instead of comparing average times.

<div class="viz" data-viz="dots" data-unit="months of follow-up" data-groups='{"Control":[3,13,9,17,4,19,10,12,20,1,12,17,6,1,7,16,10,6],"Treatment":[18,10,9,24,26,5,25,25,7,20,19,15,27,15,16,26,16,7]}' data-caption="Time on study per patient. Many of the longest times are censored (patients still event-free at last visit), so averaging these times would be misleading."></div>

The χ² statistic is (O₁ − E₁)² / V = 4.94 with 1 degree of freedom. Under H₀, χ² values follow the curve below, and large values are in the shaded upper tail.

<div class="viz" data-viz="chisq" data-stat="4.94" data-df="1"></div>

What to notice: 4.94 is out in the right tail, but not dramatically so (p = 0.026). With only 17 events in total, the evidence is modest, which is why the hazard-ratio CI later is so wide.

## When to use it

Use the log-rank test when **all** of these are true:

- The outcome is **time to an event** with **right censoring** (death, progression, relapse)
- You are comparing **two or more independent groups** (treatment arms, biomarker status)
- You expect the **hazards to be roughly proportional** (the curves separate steadily and do not cross)
- Censoring is non-informative and similar in the groups

**Clinical scenario.** In the 36-patient oncology trial introduced on the [Kaplan–Meier page](#/test/kaplan-meier), progression-free survival is compared between Control (n = 18) and Treatment (n = 18). The KM curves suggest the new agent delays progression. *Is the difference larger than we would expect by chance?* The same data are analysed with [Cox regression](#/test/cox-ph) to get a hazard ratio.

**Don't use it when:**

- You need an **effect size or adjustment** for covariates such as age → use [Cox proportional hazards regression](#/test/cox-ph)
- You only want a description of one group → use the [Kaplan–Meier estimator](#/test/kaplan-meier)
- The curves **cross** (hazards are clearly non-proportional), because the early and late differences cancel and the test loses power. Use a weighted log-rank test, a pre-specified max-combo test, or compare the **restricted mean survival time** (RMST), described on the [Cox page](#/test/cox-ph)
- The outcome is binary at a fixed time with no censoring → use the [two-proportion z-test](#/test/two-proportion-z) or [Fisher's exact test](#/test/fisher-exact)
- You want to test for a trend across ordered groups → use Cox regression with an ordinal covariate (or the TREND option of the SAS STRATA statement)

## Assumptions

1. **Censoring is non-informative and unrelated to group.** If patients in one arm drop out for reasons linked to their prognosis, the test is biased. Compare the censoring patterns and reasons between arms.
2. **Hazards are proportional (the curves do not cross).** The test is most powerful for this situation. Check by plotting the KM curves and the log(−log) survival curves, or with a Schoenfeld-residual test (see [Cox page](#/test/cox-ph)). If curves cross, the p-value may be misleading.
3. **Patients are independent**, with one row per patient.
4. **Event times are measured the same way in every group** (same assessment schedule and definition of the event).
5. **Group membership is fixed at time zero** (for example, randomised arm). Comparing "responders" and "non-responders" defined after baseline creates immortal-time bias. Use a landmark or time-dependent analysis instead.

## Hypotheses

- **H₀:** S₁(t) = S₂(t) for all t (the survival curves, and therefore the hazard functions, are identical in the two groups)
- **H₁:** S₁(t) ≠ S₂(t) for some t (two-sided; the default). The test has the most power when the hazards differ by a constant ratio.

## Example data

The same 36 patients as the [Kaplan–Meier page](#/test/kaplan-meier). Months from randomisation to progression or death; **+** marks a censored time. Age is used only for the stratified analysis (stratum: under 62 versus 62 and over).

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
# Standard log-rank test
lr <- survdiff(Surv(time, event) ~ arm, data = trial)
lr
1 - pchisq(lr$chisq, df = 1)          # p-value with more digits

# Stratified log-rank: compare arms within age strata (<62 vs 62+)
trial$agegrp <- factor(ifelse(trial$age >= 62, "62+", "<62"))
survdiff(Surv(time, event) ~ arm + strata(agegrp), data = trial)
```

**What the code does**

- `Surv(time, event)` uses 1 = event and 0 = censored.
- `survdiff(Surv(time, event) ~ arm)` is the log-rank test. By default `rho = 0`, which is the standard log-rank test. `rho = 1` would give a Peto–Peto-type weighted test that emphasises early differences.
- The printed table gives observed and expected events per arm. The `(O-E)^2/V` column is the χ² statistic. `p= 0.03` is rounded, so I recompute the p-value with `1 - pchisq(...)`.
- `strata(agegrp)` makes the **stratified log-rank test**. The comparison between arms is done within each age stratum, and the (O − E) and V values are then summed over strata. Use it for a stratification factor used in randomisation.
- `survdiff` does not report a hazard ratio or CI. Use `coxph` for that.

**Output**

```text
Call:
survdiff(formula = Surv(time, event) ~ arm, data = trial)
               N Observed Expected (O-E)^2/E (O-E)^2/V
arm=Control   18       10     5.87      2.90      4.94
arm=Treatment 18        7    11.13      1.53      4.94
 Chisq= 4.9  on 1 degrees of freedom, p= 0.03 
[1] 0.02616886
Call:
survdiff(formula = Surv(time, event) ~ arm + strata(agegrp), 
    data = trial)
               N Observed Expected (O-E)^2/E (O-E)^2/V
arm=Control   18       10     5.96      2.73      4.75
arm=Treatment 18        7    11.04      1.48      4.75
 Chisq= 4.8  on 1 degrees of freedom, p= 0.03 
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
import numpy as np
from scipy import stats
from lifelines.statistics import logrank_test

ctl = trial[trial.arm == "Control"]
trt = trial[trial.arm == "Treatment"]

# Standard log-rank test
res = logrank_test(ctl["time"], trt["time"],
                   event_observed_A=ctl["event"], event_observed_B=trt["event"])
print(f"Log-rank: chi2 = {res.test_statistic:.2f}, df = 1, p = {res.p_value:.4f}")

# Stratified log-rank (lifelines has no built-in): sum (O - E) and V over strata
trial["agegrp"] = np.where(trial["age"] >= 62, "62+", "<62")
O_minus_E, V = 0.0, 0.0
for _, s in trial.groupby("agegrp"):
    for t in np.sort(s.loc[s.event == 1, "time"].unique()):
        n  = (s.time >= t).sum()                       # at risk, both arms
        n1 = ((s.time >= t) & (s.arm == "Control")).sum()
        d  = ((s.time == t) & (s.event == 1)).sum()    # events (both arms)
        d1 = ((s.time == t) & (s.event == 1) & (s.arm == "Control")).sum()
        O_minus_E += d1 - d * n1 / n
        if n > 1:
            V += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1)
chi2 = O_minus_E**2 / V
print(f"Stratified log-rank: chi2 = {chi2:.2f}, df = 1, p = {stats.chi2.sf(chi2, 1):.4f}")
```

**What the code does**

- `logrank_test(durations_A, durations_B, event_observed_A, event_observed_B)` is the standard two-group log-rank test from `lifelines`. The result has `test_statistic` (χ² with 1 df) and `p_value`.
- `lifelines` does not provide a stratified log-rank test, so the second part computes it directly. At each event time within each age stratum it accumulates the observed minus expected Control events and the hypergeometric variance, then forms the χ² from the sums. This is the same calculation R's `strata()` does.
- For more than two groups use `multivariate_logrank_test(durations, groups, events)`.

**Output**

```text
Log-rank: chi2 = 4.94, df = 1, p = 0.0262
Stratified log-rank: chi2 = 4.75, df = 1, p = 0.0292
```

### SAS

```sas
data trial;
  input id arm $ months event age;   /* event: 1 = event, 0 = censored */
  agegrp = (age >= 62);              /* stratification factor: 0 = <62, 1 = 62+ */
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

/* Standard log-rank test (and Wilcoxon for comparison) */
proc lifetest data=trial notable;
  time months*event(0);
  strata arm / test=(logrank wilcoxon);
run;

/* Stratified log-rank test: arms compared within age strata */
proc lifetest data=trial notable;
  time months*event(0);
  strata agegrp / group=arm test=logrank;
run;
```

**What the code does**

- `time months*event(0);` identifies the time variable and the event indicator. The value in parentheses (0) means **censored**.
- `strata arm / test=(logrank wilcoxon);` compares the arms with the log-rank test and also the Gehan–Wilcoxon test, which weights early events more heavily. Choose the primary test in the statistical analysis plan, not after seeing the results. `test=peto` is also available.
- `notable` suppresses the long product-limit table so that only the tests are printed.
- In the second call the roles are different. `strata agegrp` defines the **strata**, and `group=arm` names the variable whose levels are **compared**. The result is the stratified log-rank test, whose χ² should agree with R and Python.
- The first call prints a "Test of Equality over Strata" table and the second a "Test of Equality over Group" table (each with χ², DF, Pr > Chi-Square). It was not run here, but it should show χ² of about 4.94 and p = 0.026 for the first call and about 4.75 and p = 0.029 for the stratified call.

## Interpreting the output

- **Observed vs expected events.** Control had 10 events but only 5.87 were expected under the null. Treatment had 7 events against 11.13 expected. The Control arm has more events than expected, and the Treatment arm fewer.
- **χ² = 4.94, df = 1, p = 0.026.** If the survival curves were truly identical, a difference this large would occur in about 2.6% of trials. R and Python give the same χ² and p-value.
- **Stratified log-rank: χ² = 4.75, p = 0.029.** Adjusting the comparison for the age group changes little. The conclusion is unchanged. R and the hand-coded Python version agree.
- **A rough effect size.** (O/E) in Control = 10/5.87 = 1.70 and in Treatment = 7/11.13 = 0.63, giving an approximate hazard ratio for Treatment versus Control of 0.63/1.70 ≈ 0.37. The proper estimate from [Cox regression](#/test/cox-ph) is 0.32 (95% CI 0.11 to 0.91).
- **Statistically significant is not the same as clinically important.** With 17 events the test is only moderately powered. Always show the KM curves, the medians, and a hazard ratio with its CI.

**How to report it:**
> Progression-free survival differed between arms (log-rank χ²(1) = 4.94, p = 0.026): 7 events occurred in the treatment arm (11.1 expected) and 10 in the control arm (5.9 expected). The age-stratified log-rank test gave a similar result (χ²(1) = 4.75, p = 0.029). The hazard ratio is 0.32 (95% CI 0.11 to 0.91; Cox model).

## Check your understanding

<details class="check"><summary>1. Two arms have the same median survival (14 months), but one KM curve is above the other for the first year and below it afterwards. What happens to the log-rank test?</summary>

The curves cross, so early excess events in one arm are offset by late excess events in the other. O − E partly cancels and the test can miss a real difference. Look at the curves first, and use a pre-specified alternative such as a weighted test, max-combo, or a restricted mean survival time comparison.
</details>

<details class="check"><summary>2. The log-rank test gives p = 0.026. A colleague says "so the treatment reduces progression by 74%". What is missing?</summary>

The log-rank test provides no effect size. The 74% (a hazard ratio of 0.26) is a number the test never produced. Report the hazard ratio from a Cox model (here 0.32, 95% CI 0.11 to 0.91), the KM medians or landmark survival, and note the wide CI from only 17 events.
</details>

<details class="check"><summary>3. Patients who responded within 3 months are compared with those who did not, using a log-rank test from randomisation. What is wrong?</summary>

Response status is only known after follow-up begins, so responders had to survive long enough to respond. This **immortal-time bias** makes responders look better. Use a landmark analysis (start the clock at 3 months) or a time-dependent covariate.
</details>

<details class="check"><summary>4. What if the same 17 events had occurred in a trial of 360 patients rather than 36?</summary>

The statistic depends on events, not patients. With the same number of events and the same O/E ratios, the χ² and p would be similar. The test's power comes from the **number of events**, so a large sample with few events is still underpowered.
</details>

## Common mistakes

- **Reporting only the log-rank p-value.** It says nothing about how big the difference is. Give the hazard ratio with CI and the KM medians or landmark rates.
- **Using it when the curves cross.** Early benefit and late harm (or the reverse) cancel out and the test has little power. Inspect the curves first, and prefer a pre-specified alternative such as a weighted test or RMST difference.
- **Choosing the weighting after looking at the data.** Switching from log-rank to Wilcoxon because it gives a smaller p-value is p-hacking. Pre-specify it (ICH E9).
- **Ignoring stratification factors.** If randomisation was stratified, the primary analysis should usually be the stratified log-rank test and a stratified Cox model.
- **Comparing the median survival times and calling it a test.** Medians can be equal while curves differ, and CIs for medians are wide.
- **Defining groups by something that happens after baseline** (for example, "responders"), which gives immortal-time bias.
- **Backwards censoring code.** In SAS the number in `event(0)` is the censoring value. In R and Python, 1 = event.

## Related tests

- [Kaplan–Meier estimator](#/test/kaplan-meier): the survival curves being compared
- [Cox proportional hazards regression](#/test/cox-ph): hazard ratio, covariate adjustment, and the PH check
- [Cochran–Mantel–Haenszel test](#/test/cmh): the analogous stratified test for binary outcomes
- [Fisher's exact test](#/test/fisher-exact): binary outcome at a fixed time with no censoring
- [Mann–Whitney U test](#/test/mann-whitney): compare a continuous outcome with no censoring

## Key takeaways

- At each event time, the log-rank test compares observed events in an arm with the number expected if the arms had identical hazards, then sums O − E.
- It is most powerful when the hazard ratio is constant. If the curves cross, it can miss real differences.
- It gives a p-value only. Always add a hazard ratio with CI and KM medians.
- Power depends on the number of events, not the number of patients.
- Groups must be defined at time zero, or you risk immortal-time bias.

## References

- Mantel N. Evaluation of survival data and two new rank order statistics arising in its consideration. *Cancer Chemother Rep.* 1966;50:163–170.
- Bland JM, Altman DG. The logrank test. *BMJ.* 2004;328:1073.
- ICH E9: Statistical Principles for Clinical Trials (1998), and the E9(R1) addendum on estimands (2019).
- R documentation: `?survdiff`; lifelines documentation: `lifelines.statistics.logrank_test`; SAS/STAT User's Guide: The LIFETEST Procedure (STRATA statement, TEST= and GROUP= options)
