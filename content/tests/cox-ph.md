# Cox proportional hazards regression

> Models the hazard (instantaneous event rate) of a time-to-event outcome as a function of covariates, giving **hazard ratios** with confidence intervals while handling censoring and without assuming a shape for the baseline hazard.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Time-to-event outcome with censoring; hazard ratio or covariate adjustment | Time to event plus event/censor indicator | Proportional hazards; non-informative censoring; log-linear covariates | Progression-free survival, new agent vs control, adjusted for age (17 events) |

## Learning objectives

- Explain the hazard and the hazard ratio in plain language
- Describe how the Cox model uses risk sets and censoring without specifying a baseline hazard
- Fit a Cox model in R and Python and read HR, 95% CI, and the PH check
- Judge when the hazard ratio is a poor summary and consider RMST instead

## What it is

Where the [log-rank test](#/test/log-rank) only tells you *whether* two survival curves differ, the Cox model tells you *by how much*, and lets you adjust for other factors such as age. It models the **hazard** h(t), the instantaneous risk of the event at time t among patients who are still event-free:

h(t | x) = h₀(t) · exp(β₁x₁ + β₂x₂ + …)

The baseline hazard h₀(t) can have any shape and is never estimated explicitly (the model is *semiparametric*). Parameters are estimated by maximising a **partial likelihood** that compares, at each event time, the patient who had the event with everyone else still at risk. Censored patients contribute to the risk sets until they leave.

- **Effect size:** the **hazard ratio** HR = exp(β) with a 95% CI. HR = 0.5 means that at any time the event rate in the group is half that of the reference. For a continuous covariate the HR is per one-unit increase.
- **Test statistics:** Wald z = β̂ / SE(β̂) for each coefficient, and the likelihood-ratio, Wald, and score tests for the whole model. With one binary covariate the score test equals the log-rank test.
- **Key assumption:** the hazard ratio is **constant over time** (proportional hazards, PH).
- **Ties:** when several patients have the event at the same time, the partial likelihood needs an approximation. **R and lifelines use Efron's method; SAS PHREG uses Breslow's by default.**

## The intuition

The Cox model compares **rates**, not totals. The hazard is the instantaneous event rate among patients still event-free. The **signal** is the log hazard ratio (here −1.145, i.e. HR 0.32): how much lower the event rate is on treatment at any moment. The **noise** is its standard error (0.538), which depends mainly on the number of **events**, not the number of patients. The Wald z is −1.145 / 0.538 = −2.13.

Analogy: two roads with traffic cameras. The hazard is the accident rate per hour among cars currently on the road. If road B has a third of the accident rate of road A at every hour of the day, the hazard ratio is 0.33, whether you watch the roads for one hour or ten. The model assumes this constancy (proportional hazards) and does not care what the rate looks like over the day.

The plot shows the times at which events happened. Notice that control events are earlier and more concentrated (many in the first 13 months), while treatment events are fewer and later. What the plot cannot show is the censored patients, who contribute follow-up until they leave and are the reason a simple comparison of average times would mislead. The Cox model handles them through the risk sets.

<div class="viz" data-viz="dots" data-unit="months to progression or death (events only)" data-groups='{"Control (10 events)":[3,13,4,19,10,1,12,6,1,7],"Treatment (7 events)":[5,25,20,19,15,16,7]}' data-caption="Event times only. Censored patients (8 in control, 11 in treatment) are not shown but they still count in the risk sets."></div>

<div class="viz" data-viz="z" data-stat="-2.13"></div>

**What to notice:** the treatment coefficient lies just beyond the usual cut-off, with p = 0.033, and the CI (0.11 to 0.91) is wide because only 17 events inform it. It is a hint of benefit, not a precise estimate.

## When to use it

Use Cox regression when **all** of these are true:

- The outcome is a **time to event** with **right censoring**
- You need a **hazard ratio** for a treatment or exposure, or want to **adjust** for covariates (age, stage, stratification factors)
- The hazards are **approximately proportional** over follow-up
- There are enough events: a common guideline is at least 10 events per covariate in the model

**Clinical scenario.** The 36-patient oncology trial from the [Kaplan–Meier](#/test/kaplan-meier) and [log-rank](#/test/log-rank) pages. *What is the hazard ratio for progression or death with the new agent versus control, adjusted for age?* The model includes treatment arm and age at randomisation.

**Don't use it when:**

- The curves **cross** or the hazard ratio clearly changes over time → compare the **restricted mean survival time** (RMST), use a time-varying effect or stratify by the offending variable
- You only need a description of survival in each group → use the [Kaplan–Meier estimator](#/test/kaplan-meier)
- You only need a significance test of two groups with no adjustment → the [log-rank test](#/test/log-rank) is enough
- The outcome is binary at a fixed time with no censoring → use [logistic regression](#/test/logistic-regression)
- The outcome is a count of events per person-time → use [Poisson regression](#/test/poisson-regression)
- Events are recurrent or patients are clustered → use extensions (frailty or marginal models), or [GEE](#/test/gee) for correlated data

## Assumptions

1. **Proportional hazards.** The HR for each covariate is constant over time. Check with a test of the **scaled Schoenfeld residuals** against time (`cox.zph` in R, `proportional_hazard_test` in lifelines, `ASSESS PH` in SAS), with plots of the residuals, and by looking at whether the KM curves cross. A non-significant test is reassuring but is weak with few events.
2. **Non-informative censoring.** Censoring is unrelated to prognosis after accounting for covariates.
3. **Log-linearity for continuous covariates.** The log hazard changes linearly with age. Check with martingale residuals, or add splines or categories.
4. **Independent observations**, one row per patient.
5. **Enough events relative to covariates.** Around 10 events per variable. Here there are 17 events and 2 covariates, which is on the low side, and this is a teaching example rather than a real analysis.
6. **No influential outliers.** Check dfbeta residuals if estimates look unstable.

## Hypotheses

For each coefficient β (log hazard ratio):

- **H₀:** β = 0, that is HR = 1 (the covariate has no association with the hazard, adjusting for the others)
- **H₁:** β ≠ 0, that is HR ≠ 1 (two-sided)

For the proportional-hazards check the null hypothesis is that the **residuals are unrelated to time**, which means the PH assumption holds, so here a *large* p-value is the reassuring one.

## Example data

The same 36 patients used on the [Kaplan–Meier](#/test/kaplan-meier) and [log-rank](#/test/log-rank) pages. Months from randomisation to progression or death; **+** marks a censored time. Age (years) is the covariate for adjustment.

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
# Cox model: treatment effect adjusted for age; Efron tie handling (R default, stated explicitly)
fit <- coxph(Surv(time, event) ~ arm + age, data = trial, ties = "efron")
summary(fit)

# Proportional-hazards check (scaled Schoenfeld residuals vs. transformed time)
cox.zph(fit)

# Restricted mean survival time up to 24 months, by arm (alternative if PH fails)
print(survfit(Surv(time, event) ~ arm, data = trial), rmean = 24)
```

**What the code does**

- `coxph(Surv(time, event) ~ arm + age)` fits the model. `arm` is a factor with Control as the reference level, so the coefficient `armTreatment` is Treatment versus Control. Age enters as a linear term, so its HR is per year.
- `ties = "efron"` is the **R default** (written out for clarity). Efron is more accurate than Breslow when there are many ties. SAS defaults to Breslow, so set `TIES=EFRON` in SAS to match.
- `summary(fit)` gives `coef` (log HR), `exp(coef)` (the HR), the standard error, Wald z and p, and a second table with the 95% CI for the HR. It also reports concordance and the overall likelihood-ratio, Wald, and score tests.
- `cox.zph(fit)` tests the PH assumption for each covariate and globally, using scaled Schoenfeld residuals against a transformed time (Kaplan–Meier scale by default). Use `plot(cox.zph(fit))` to see the residual trend.
- `survfit(..., rmean = 24)` prints the **restricted mean survival time** up to 24 months for each arm, the area under the KM curve up to 24 months. Use it as an alternative effect measure when PH is doubtful. The `median` and CI columns in that table use R's default `log` CI type.

**Output**

```text
Call:
coxph(formula = Surv(time, event) ~ arm + age, data = trial, 
    ties = "efron")
  n= 36, number of events= 17 
                  coef exp(coef)  se(coef)      z Pr(>|z|)  
armTreatment -1.145089  0.318195  0.538388 -2.127   0.0334 *
age           0.007232  1.007258  0.028788  0.251   0.8016  
             exp(coef) exp(-coef) lower .95 upper .95
armTreatment    0.3182     3.1427    0.1108    0.9141
age             1.0073     0.9928    0.9520    1.0657
Concordance= 0.674  (se = 0.064 )
Likelihood ratio test= 4.86  on 2 df,   p=0.09
Wald test            = 4.62  on 2 df,   p=0.1
Score (logrank) test = 5.05  on 2 df,   p=0.08
         chisq df    p
arm    0.98570  1 0.32
age    0.00155  1 0.97
GLOBAL 0.98684  2 0.61
Call: survfit(formula = Surv(time, event) ~ arm, data = trial)
               n events rmean* se(rmean) median 0.95LCL 0.95UCL
arm=Control   18     10   13.3      2.10     13       7      NA
arm=Treatment 18      7   19.9      1.45     25      19      NA
    * restricted mean with upper limit =  24 
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
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import proportional_hazard_test
from lifelines.utils import restricted_mean_survival_time

# lifelines needs numeric covariates: 0 = Control, 1 = Treatment
df = trial.assign(treat=(trial["arm"] == "Treatment").astype(int)).drop(columns="arm")

cph = CoxPHFitter()                      # Efron ties (the only method lifelines offers)
cph.fit(df, duration_col="time", event_col="event")
cols = ["coef", "exp(coef)", "se(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]
print(cph.summary[cols].round(3).to_string())

# Proportional-hazards test (Schoenfeld residuals, Kaplan-Meier time scale)
proportional_hazard_test(cph, df, time_transform="km").print_summary(decimals=3)

# RMST up to 24 months, by arm
for arm, g in trial.groupby("arm"):
    kmf = KaplanMeierFitter().fit(g["time"], g["event"])
    print(f"{arm}: RMST(24) = {restricted_mean_survival_time(kmf, t=24):.2f} months")
```

**What the code does**

- `CoxPHFitter` needs a data frame with numeric covariates, so the arm is recoded as `treat` (0 = Control, 1 = Treatment). The reference group is the one coded 0.
- `fit(df, duration_col="time", event_col="event")` uses all other columns as covariates. Be careful not to leave an ID column in the data frame.
- lifelines handles tied event times with **Efron's method**, the same as R's default. (The `baseline_estimation_method = breslow` shown in `print_summary` refers only to how the baseline hazard is estimated and does not change the coefficients.)
- `cph.summary` holds the coefficient, HR (`exp(coef)`), SE, 95% CI for the HR, and p-value.
- `proportional_hazard_test(..., time_transform="km")` is the Schoenfeld-residual test of proportional hazards on the Kaplan–Meier time scale, the same time transformation R uses. The statistics are not identical to R's `cox.zph` because the variance calculation differs, but the conclusion is the same. `cph.check_assumptions(df)` also reports this test and can draw plots.
- `restricted_mean_survival_time(kmf, t=24)` gives the RMST up to 24 months from a fitted KM curve.

**Output**

```text
            coef  exp(coef)  se(coef)  exp(coef) lower 95%  exp(coef) upper 95%      p
covariate                                                                             
age        0.007      1.007     0.029                0.952                1.066  0.802
treat     -1.145      0.318     0.538                0.111                0.914  0.033
Proportional hazards test (time_transform = km)
       test_statistic     p  -log2(p)
age             0.006 0.940     0.090
treat           0.602 0.438     1.191
Control: RMST(24) = 13.28 months
Treatment: RMST(24) = 19.95 months
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
proc phreg data=trial;
  class arm(ref='Control') / param=ref;           /* Control = reference group */
  model months*event(0) = arm age / ties=efron risklimits;
  assess ph / resample seed=2024;                 /* PH check: Schoenfeld-type resampling */
run;
```

**What the code does**

- `class arm(ref='Control') / param=ref;` declares `arm` as categorical with Control as the reference level, so the hazard ratio is Treatment versus Control. `param=ref` (reference coding) is the PHREG default and is written out for clarity; `ref='Control'` matters because PHREG's default reference level is the *last* level, which here would be Treatment and would flip the sign of the coefficient.
- `model months*event(0) = arm age;` fits the Cox model. The value in parentheses (0) is the **censoring** value.
- `ties=efron` is needed to match R and lifelines. **The SAS default is `TIES=BRESLOW`**, which gives slightly different estimates when there are tied event times, as in this data set.
- `risklimits` prints the hazard ratio with its 95% confidence limits in the "Analysis of Maximum Likelihood Estimates" table. Use `alpha=` in the MODEL statement to change the confidence level. The `Hazard Ratio` column is exp(β̂).
- `assess ph / resample;` checks the PH assumption using cumulative sums of Schoenfeld residuals and gives a supremum test p-value for each covariate. The `seed=` option makes the simulated p-values reproducible. The graphs need `ods graphics on`.
- For a classical Schoenfeld test, use `output out=resid ressch=sch_arm sch_age;` and correlate the residuals with time or a rank of event time.
- The output is not shown because SAS cannot be run here. With `ties=efron`, the log HR, SE and HR CI should agree with R and Python.

## Interpreting the output

- **Sample.** 36 patients, 17 events (R: `n= 36, number of events= 17`).
- **Treatment: HR = 0.318 (95% CI 0.111 to 0.914), p = 0.033.** At any time during follow-up, patients on the new agent had about 68% lower hazard of progression or death than controls, adjusted for age. The CI excludes 1, consistent with the p-value, but it is wide (roughly 9% to 89% lower), as expected with 17 events. The coefficient is −1.145 (SE 0.538).
- **Age: HR = 1.007 per year (95% CI 0.952 to 1.066), p = 0.80.** There is no evidence that age is associated with progression in this small sample. The wide CI means a clinically relevant effect is not excluded.
- **Concordance = 0.674.** In 67% of comparable pairs of patients, the one with the higher predicted risk had the event first (0.5 would be chance).
- **Overall model tests (2 df).** Likelihood ratio χ² = 4.86, p = 0.09. The overall test is weaker than the treatment test because the age term adds a degree of freedom without adding information.
- **Proportional hazards check.** R `cox.zph`: arm p = 0.32, age p = 0.97, global p = 0.61. lifelines: arm p = 0.44, age p = 0.94. No evidence against PH, though with 17 events the test has low power. Here the KM curves also do not cross.
- **R and Python match.** Coefficients, SEs, HRs, CIs and p-values agree. Only the PH test statistics differ slightly because of different variance estimates.
- **Ties matter slightly.** In R, `ties = "breslow"` gives a treatment coefficient of −1.124 (HR 0.325) against −1.145 (HR 0.318) for Efron. SAS with its default would report the Breslow value, so set `TIES=EFRON` when comparing software.
- **RMST (restricted to 24 months).** Control 13.3 months and Treatment 19.9 months (R prints 13.3 and 19.9; lifelines gives 13.28 and 19.95). Patients on treatment live progression-free on average 6.7 months longer within the first 24 months. This difference is easy to interpret and does not depend on proportional hazards.

**How to report it:**
> In a Cox proportional hazards model adjusted for age, treatment was associated with a lower hazard of progression or death compared with control (HR 0.32, 95% CI 0.11 to 0.91; p = 0.033; 17 events in 36 patients). Age was not associated with the outcome (HR 1.01 per year, 95% CI 0.95 to 1.07). There was no evidence against proportional hazards (Schoenfeld residual test, global p = 0.61). Ties were handled with the Efron method.

## Check your understanding

<details class="check"><summary>1. A HR of 0.32 for treatment means patients on treatment are 68% less likely to have had the event by 24 months. True?</summary>

Not exactly. The HR compares the instantaneous event **rates** among patients still at risk. It is not a ratio of cumulative risks at a fixed time. It describes the relative speed of events, and cumulative differences also depend on the baseline hazard.
</details>

<details class="check"><summary>2. The Schoenfeld PH test gives p = 0.61 in this 17-event trial. Is proportional hazards confirmed?</summary>

No. With few events the test has low power, so a high p-value is only mildly reassuring. Combine it with residual plots and check whether the KM curves cross. If PH is doubtful, report RMST (13.3 vs 19.9 months here).
</details>

<details class="check"><summary>3. You want to add 8 more covariates to this model. What is the problem?</summary>

With 17 events, a common rule is about 10 events per covariate, so one or two covariates are the practical limit. More covariates make coefficients unstable and the model overfit. Count events, not patients.
</details>

<details class="check"><summary>4. The same trial is re-analysed in SAS and gives a treatment coefficient of −1.124 instead of −1.145. Is something wrong?</summary>

Probably not. SAS PHREG defaults to Breslow's method for tied event times; R and lifelines use Efron's. Set `TIES=EFRON` in SAS to match.
</details>

## Common mistakes

- **Not checking proportional hazards.** An HR is only a single summary if the effect is constant over time. If the PH test is significant or the curves cross, report time-specific effects, stratify, or use RMST.
- **Reading the HR as a risk ratio or "x times faster".** The HR compares instantaneous rates. It is not the ratio of the proportions with the event by a fixed time.
- **Different ties methods across software.** R and lifelines default to Efron, SAS to Breslow. Set `TIES=EFRON` in SAS.
- **Wrong reference level or coding.** In R the first factor level is the reference, in lifelines it is the group coded 0, and in SAS it is set in the CLASS statement. Check the sign of the coefficient.
- **Too many covariates for the number of events.** With 17 events, adding more than one or two variables leads to unstable, overfit models. Count *events*, not patients.
- **Leaving an ID column or a post-baseline variable in the model.** In lifelines every extra column is a covariate. Post-baseline variables (for example, response status) cause immortal-time bias.
- **Treating a non-significant PH test as proof.** Small trials have low power to detect non-proportionality. Look at the plots as well.
- **Reporting only p-values.** Always give the HR, its 95% CI, and the number of events.

## Related tests

- [Kaplan–Meier estimator](#/test/kaplan-meier): descriptive survival curves, medians, and landmark rates
- [Log-rank test](#/test/log-rank): unadjusted comparison, equal to the score test of a Cox model with one binary covariate
- [Logistic regression](#/test/logistic-regression): binary outcome at a fixed time
- [Poisson & negative binomial regression](#/test/poisson-regression): event counts and rates per person-time
- [Linear regression](#/test/linear-regression): continuous outcomes

## Key takeaways

- The Cox model gives a **hazard ratio**: the relative event rate between groups, assumed constant over time.
- It does not specify the baseline hazard and handles censoring through risk sets.
- Precision depends on the **number of events** (about 10 per covariate), not the number of patients.
- Check proportional hazards; if it fails or is doubtful, report RMST or time-specific effects.
- Software defaults for ties differ (Efron in R and lifelines, Breslow in SAS).

## References

- Cox DR. Regression models and life-tables. *J R Stat Soc B.* 1972;34:187–220.
- Therneau TM, Grambsch PM. *Modeling Survival Data: Extending the Cox Model.* Springer; 2000.
- Royston P, Parmar MKB. Restricted mean survival time: an alternative to the hazard ratio for the design and analysis of randomized trials with a time-to-event outcome. *BMC Med Res Methodol.* 2013;13:152.
- R documentation: `?coxph`, `?cox.zph`; lifelines documentation: `CoxPHFitter`; SAS/STAT User's Guide: The PHREG Procedure (MODEL statement TIES= option, ASSESS statement)
