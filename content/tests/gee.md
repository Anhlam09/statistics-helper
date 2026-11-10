# Generalized estimating equations (GEE)

> Fits a regression model to correlated outcomes, such as a binary response measured at several visits on the same patient, and estimates population-averaged effects with standard errors that stay valid even if the correlation structure is misspecified.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Repeated or clustered outcomes, population-averaged effect wanted | Binary, count or continuous, correlated within patient | Independent clusters; correct mean model; enough clusters; MCAR missingness | Asthma control (yes/no) at 3 visits, drug vs placebo |

## Learning objectives

- Explain why repeated measures on the same patient are correlated and why that matters for SEs
- Describe the three GEE choices: mean model, distribution, and working correlation
- Fit a GEE in R and Python and read the robust SEs and the population-averaged OR
- Distinguish population-averaged (GEE) from subject-specific (mixed model) odds ratios

## What it is

GEE extends generalised linear models (logistic, Poisson, linear) to **clustered or repeated measurements**. You choose three things: the mean model (for example a logit link with treatment and time), the outcome distribution (binomial), and a **working correlation** structure for the repeated measurements within a patient (independence, exchangeable, AR(1), unstructured). The coefficients are estimated by solving estimating equations that use this working correlation, and the standard errors come from the **robust (sandwich) variance estimator**, which remains consistent for large numbers of clusters even if the working correlation is wrong.

The key interpretation point is that GEE targets **population-averaged (marginal)** effects. For a logit link, exp(β) is the odds ratio comparing the *populations* of treated and untreated patients. A subject-specific odds ratio from a random-effects logistic model (GLMM) compares two patients with the same random effect, and is numerically larger in magnitude.

- **Test statistic:** Wald z = β̂ / SE_robust (some programs report z²)
- **Effect size:** the population-averaged odds ratio exp(β̂) with a 95% CI exp(β̂ ± 1.96 × SE_robust); for count outcomes, the rate ratio

## The intuition

Forty patients each contribute three yes/no answers, so the 120 rows contain far less than 120 pieces of independent information. If a patient is well controlled at week 4, they are more likely than a random patient to be well controlled at week 8. It is like polling 40 households and asking each person three times: three answers from one person are not three independent voters. A naive logistic regression would treat the 120 rows as independent and report SEs that are too small.

GEE keeps the usual regression for the **average** response, adds a **working guess** of how visits within a patient are correlated (here "exchangeable": one common correlation, estimated as 0.22), and then computes the SE with a **sandwich (robust) estimator** that checks the guess against the data's actual scatter. Signal: the log-OR for drug (0.678). Noise: its robust SE (0.474). Ratio: Wald z = 1.43.

The dot plot counts, for each patient, the number of the three visits at which asthma was well controlled (0 to 3).

<div class="viz" data-viz="dots" data-unit="visits with good control (of 3)" data-groups='{"Placebo":[0,0,1,0,2,0,3,3,2,0,1,1,0,1,0,0,0,1,0,1],"Drug":[2,3,0,2,1,1,2,0,2,1,0,1,2,1,0,1,2,3,0,1]}' data-caption="Notice that patients cluster at the extremes: several patients are controlled at 0 or 3 visits, which is what within-patient correlation looks like. The drug group has a higher average (1.25 vs 0.80 visits out of 3), but the overlap is large with only 20 patients per arm."></div>

Under H₀ (no drug effect on the population-averaged odds), the Wald z follows the standard normal curve. Our z = 1.43 is not far enough into the tail to be convincing (p = 0.153).

<div class="viz" data-viz="z" data-stat="1.43"></div>

## When to use it

Use GEE when **all** of these are true:

- The outcome is **binary, count or continuous**, measured **repeatedly** on each patient (or on patients clustered in sites, families or clinics)
- You want an effect that describes the **average across the population**, which is usually what a trial's estimand is
- There are enough clusters for the sandwich variance to be reliable (a rough guide is at least 30–40; fewer is risky, but there is no sharp cutoff)
- Missing visits, if any, are plausibly **missing completely at random**, or you will weight the equations (weighted GEE)

**Clinical scenario.** In a trial, 40 patients with asthma are randomised to placebo or a new inhaled therapy (20 per arm). At weeks 4, 8 and 12 each patient is classified as having well-controlled asthma (yes/no). The three assessments on one patient are correlated. *Does the therapy raise the odds of good control, averaged over the study period?*

**Don't use it when:**

- You have one measurement per patient → use [logistic regression](#/test/logistic-regression)
- The outcome is continuous with missing visits that depend on earlier outcomes → use a [mixed model (MMRM)](#/test/mixed-models), which is valid under MAR; standard GEE needs MCAR
- You need subject-specific effects or predictions for an individual patient → fit a random-effects model (GLMM)
- Only two paired binary measurements per patient are compared → use [McNemar's test](#/test/mcnemar)
- There are very few clusters (below roughly 20) → the robust SEs are too optimistic; consider a mixed model or small-sample corrections
- The outcome is a count of events over follow-up time → see [Poisson & negative binomial regression](#/test/poisson-regression) (which GEE can extend to repeated counts)

## Assumptions

1. **Clusters (patients) are independent**; observations within a cluster may be correlated. Check the study design.
2. **The mean model is correct.** The link function and the way covariates enter (for example linear in `week`) must be adequate, because GEE is not protected against a wrong mean model. Try visit as a factor or add a treatment × time interaction and compare.
3. **The number of clusters is large enough** for the sandwich variance estimator to be accurate. With few clusters, SEs and CIs are too narrow.
4. **Missing data are missing completely at random** (or covered by inverse-probability weights). Unlike a likelihood-based mixed model, standard GEE is not valid under MAR.
5. **The working correlation is a reasonable approximation.** It need not be exactly right for valid inference, but a good choice improves efficiency. Compare the estimated correlations with the choice made.

## Hypotheses

For the treatment coefficient β_trt in a logit model:

- **H₀:** β_trt = 0, that is, population-averaged OR = 1
- **H₁:** β_trt ≠ 0, that is, OR ≠ 1 (two-sided)

## Example data

40 patients, response at weeks 4, 8 and 12 (1 = asthma well controlled, 0 = not). Placebo: patients 1–20; drug: patients 21–40. Observed control rates are 20%, 30%, 30% on placebo and 40%, 35%, 50% on the drug at weeks 4, 8, 12.

| Placebo id | Wk 4, 8, 12 | Drug id | Wk 4, 8, 12 |
|---|---|---|---|
| 1 | 0 0 0 | 21 | 1 0 1 |
| 2 | 0 0 0 | 22 | 1 1 1 |
| 3 | 0 1 0 | 23 | 0 0 0 |
| 4 | 0 0 0 | 24 | 1 0 1 |
| 5 | 0 1 1 | 25 | 0 1 0 |
| 6 | 0 0 0 | 26 | 0 0 1 |
| 7 | 1 1 1 | 27 | 1 1 0 |
| 8 | 1 1 1 | 28 | 0 0 0 |
| 9 | 1 0 1 | 29 | 0 1 1 |
| 10 | 0 0 0 | 30 | 0 0 1 |
| 11 | 0 1 0 | 31 | 0 0 0 |
| 12 | 0 0 1 | 32 | 0 1 0 |
| 13 | 0 0 0 | 33 | 1 1 0 |
| 14 | 0 0 1 | 34 | 0 0 1 |
| 15 | 0 0 0 | 35 | 0 0 0 |
| 16 | 0 0 0 | 36 | 1 0 0 |
| 17 | 0 0 0 | 37 | 1 0 1 |
| 18 | 1 0 0 | 38 | 1 1 1 |
| 19 | 0 0 0 | 39 | 0 0 0 |
| 20 | 0 1 0 | 40 | 0 0 1 |

## How to run it

### R

```r
library(geepack)

# 40 patients (1-20 placebo, 21-40 drug); response = asthma well controlled (1) or not (0)
wide <- read.table(text = "
id trt w4 w8 w12
1 0 0 0 0
2 0 0 0 0
3 0 0 1 0
4 0 0 0 0
5 0 0 1 1
6 0 0 0 0
7 0 1 1 1
8 0 1 1 1
9 0 1 0 1
10 0 0 0 0
11 0 0 1 0
12 0 0 0 1
13 0 0 0 0
14 0 0 0 1
15 0 0 0 0
16 0 0 0 0
17 0 0 0 0
18 0 1 0 0
19 0 0 0 0
20 0 0 1 0
21 1 1 0 1
22 1 1 1 1
23 1 0 0 0
24 1 1 0 1
25 1 0 1 0
26 1 0 0 1
27 1 1 1 0
28 1 0 0 0
29 1 0 1 1
30 1 0 0 1
31 1 0 0 0
32 1 0 1 0
33 1 1 1 0
34 1 0 0 1
35 1 0 0 0
36 1 1 0 0
37 1 1 0 1
38 1 1 1 1
39 1 0 0 0
40 1 0 0 1
", header = TRUE)

# Wide -> long, ordered by patient (geeglm needs rows of a cluster to be contiguous)
long <- reshape(wide, direction = "long", varying = c("w4", "w8", "w12"), v.names = "resp",
                timevar = "vi", times = 1:3, idvar = "id")
long$week <- c(4, 8, 12)[long$vi]
long <- long[order(long$id, long$vi), ]

# GEE: logit link, exchangeable working correlation, robust (sandwich) SEs
fit <- geeglm(resp ~ trt + week, id = id, data = long, family = binomial,
              corstr = "exchangeable")
s <- summary(fit)$coefficients
print(round(s, 3))
est <- s["trt", "Estimate"]; se <- s["trt", "Std.err"]
cat("Population-averaged OR (drug vs placebo):", round(exp(est), 2),
    " 95% CI", round(exp(est + c(-1.96, 1.96) * se), 2), "\n")
cat("Working correlation alpha:", round(summary(fit)$corr$Estimate, 3), "\n")
```

**What the code does**

- The data are typed inline in wide format, then reshaped to **long** format (one row per patient-visit) and ordered by patient. `geeglm` requires the rows of a cluster to be contiguous and assumes the data are sorted by `id`; unsorted data silently give wrong answers.
- `id = id` names the cluster variable. `family = binomial` uses the logit link by default.
- `corstr = "exchangeable"` sets the working correlation to one common correlation between any two visits of the same patient.
- `summary(fit)$coefficients` reports **robust** (sandwich) standard errors; `geeglm` uses them by default. The "Wald" column is z², not z.
- `exp(est ± 1.96 * se)` turns the log-odds coefficient for `trt` into the population-averaged OR and its 95% CI.
- `week` is entered as a numeric linear trend; the coding is 4, 8, 12.

**Output**

```text
            Estimate Std.err  Wald Pr(>|W|)
(Intercept)   -1.477   0.555 7.076    0.008
trt            0.678   0.474 2.046    0.153
week           0.057   0.049 1.354    0.245
Population-averaged OR (drug vs placebo): 1.97  95% CI 0.78 4.99 
Working correlation alpha: 0.22 
```

### Python

```python
import io
import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf

# Same 40 patients (1-20 placebo, 21-40 drug); response = well controlled (1) or not (0)
wide = pd.read_csv(io.StringIO("""
id trt w4 w8 w12
1 0 0 0 0
2 0 0 0 0
3 0 0 1 0
4 0 0 0 0
5 0 0 1 1
6 0 0 0 0
7 0 1 1 1
8 0 1 1 1
9 0 1 0 1
10 0 0 0 0
11 0 0 1 0
12 0 0 0 1
13 0 0 0 0
14 0 0 0 1
15 0 0 0 0
16 0 0 0 0
17 0 0 0 0
18 0 1 0 0
19 0 0 0 0
20 0 0 1 0
21 1 1 0 1
22 1 1 1 1
23 1 0 0 0
24 1 1 0 1
25 1 0 1 0
26 1 0 0 1
27 1 1 1 0
28 1 0 0 0
29 1 0 1 1
30 1 0 0 1
31 1 0 0 0
32 1 0 1 0
33 1 1 1 0
34 1 0 0 1
35 1 0 0 0
36 1 1 0 0
37 1 1 0 1
38 1 1 1 1
39 1 0 0 0
40 1 0 0 1
"""), sep=r"\s+")

# Wide -> long, ordered by patient and visit
long = wide.melt(id_vars=["id", "trt"], var_name="vis", value_name="resp")
long["week"] = long["vis"].str[1:].astype(int)
long = long.sort_values(["id", "week"]).reset_index(drop=True)

# GEE: logit link, exchangeable working correlation; robust SEs are the default
model = smf.gee("resp ~ trt + week", groups="id", data=long,
                family=sm.families.Binomial(), cov_struct=sm.cov_struct.Exchangeable())
res = model.fit()
print(res.summary().tables[1])
b, se = res.params["trt"], res.bse["trt"]
print(f"Population-averaged OR (drug vs placebo): {np.exp(b):.2f}  "
      f"95% CI ({np.exp(b-1.96*se):.2f}, {np.exp(b+1.96*se):.2f})")
print(f"Working correlation alpha: {model.cov_struct.dep_params:.3f}")
```

**What the code does**

- `smf.gee(..., groups="id", ...)` defines the clusters. Rows are assumed to be sorted by cluster, so the code sorts by `id` and `week` first.
- `family=sm.families.Binomial()` uses the logit link by default.
- `cov_struct=sm.cov_struct.Exchangeable()` sets the exchangeable working correlation.
- `.fit()` returns **robust** covariance estimates by default (`cov_type="robust"`), so the printed standard errors are sandwich SEs. The p-values are z-based.
- `model.cov_struct.dep_params` is the estimated common within-patient correlation.
- The R and Python estimates, SEs and correlation agree to the digits shown; small differences can appear in the last digit because the two programs estimate the correlation parameter slightly differently.

**Output**

```text
==============================================================================
                 coef    std err          z      P>|z|      [0.025      0.975]
------------------------------------------------------------------------------
Intercept     -1.4767      0.555     -2.660      0.008      -2.565      -0.389
trt            0.6784      0.474      1.430      0.153      -0.251       1.608
week           0.0573      0.049      1.163      0.245      -0.039       0.154
==============================================================================
Population-averaged OR (drug vs placebo): 1.97  95% CI (0.78, 4.99)
Working correlation alpha: 0.220
```

### SAS

```sas
data asthma;
  input id trt w4 w8 w12;
  datalines;
1 0 0 0 0
2 0 0 0 0
3 0 0 1 0
4 0 0 0 0
5 0 0 1 1
6 0 0 0 0
7 0 1 1 1
8 0 1 1 1
9 0 1 0 1
10 0 0 0 0
11 0 0 1 0
12 0 0 0 1
13 0 0 0 0
14 0 0 0 1
15 0 0 0 0
16 0 0 0 0
17 0 0 0 0
18 0 1 0 0
19 0 0 0 0
20 0 0 1 0
21 1 1 0 1
22 1 1 1 1
23 1 0 0 0
24 1 1 0 1
25 1 0 1 0
26 1 0 0 1
27 1 1 1 0
28 1 0 0 0
29 1 0 1 1
30 1 0 0 1
31 1 0 0 0
32 1 0 1 0
33 1 1 1 0
34 1 0 0 1
35 1 0 0 0
36 1 1 0 0
37 1 1 0 1
38 1 1 1 1
39 1 0 0 0
40 1 0 0 1
;
run;

/* Long format: one row per patient-visit */
data long;
  set asthma;
  array w{3} w4 w8 w12;
  do vi = 1 to 3;
    resp = w{vi};
    week = 4 * vi;
    output;
  end;
  keep id trt week resp;
run;

proc sort data=long; by id week; run;

/* GEE: exchangeable working correlation, empirical (robust) SEs */
proc genmod data=long;
  class id;
  model resp(event='1') = trt week / dist=binomial link=logit;
  repeated subject=id / type=exch corrw;
  estimate 'Drug vs placebo (log OR)' trt 1 / exp;
run;
```

**What the code does**

- The `array` loop reshapes to long format. `week = 4 * vi` gives 4, 8, 12, matching the numeric coding used in R and Python.
- `class id` declares the cluster variable. `trt` and `week` are numeric, so they are not in `CLASS`.
- `resp(event='1')` makes SAS model the probability of `resp = 1`. Without it, PROC GENMOD models the *lower* ordered level (0) and every coefficient flips sign.
- `repeated subject=id / type=exch corrw` requests GEE with an exchangeable working correlation and prints the working correlation matrix. Other options are `type=ind`, `type=ar(1)` and `type=un`.
- The parameter table is labelled "Analysis Of GEE Parameter Estimates" with **empirical standard error estimates**, which are the robust SEs. Model-based SEs are shown only if you add `modelse` to the `REPEATED` statement; report the empirical ones.
- `estimate ... trt 1 / exp` prints the log-OR contrast and, with `exp`, the odds ratio with its confidence limits. Exponentiating the `trt` row of the parameter table gives the same OR.
- The estimates should agree closely with R and Python, because all three fit the same model with the same working correlation (exact digits can differ slightly in how the correlation parameter is estimated).

## Interpreting the output

- **Treatment coefficient.** β̂_trt = 0.678 on the log-odds scale (robust SE 0.474, z = 1.43, p = 0.153).
- **Population-averaged OR = 1.97** (95% CI 0.78 to 4.99). Across the population, the odds of well-controlled asthma are about twice as high on the drug, but the interval is wide and includes 1, so with 40 patients the data are compatible with no effect and with a large benefit.
- **Time.** The coefficient for `week` is 0.057 per week (SE 0.049, p = 0.245), an OR of about 1.06 per week, or 1.26 per 4 weeks, with no clear trend.
- **Intercept = −1.477.** On the logit scale this is the placebo log odds at week 0 (extrapolated), corresponding to a probability of about 0.19.
- **Working correlation α = 0.22.** The estimated correlation between two visits of the same patient is modest and positive. The exchangeable structure affects efficiency; the robust SEs protect inference even if the true structure differs.
- **Why not ignore the clustering?** Treating the 120 rows as independent would pretend there are 120 independent observations rather than 40 patients, which usually makes SEs too small and p-values too optimistic.

**How to report it:**
> Asthma control (yes/no) at weeks 4, 8 and 12 was analysed with a GEE logistic model (exchangeable working correlation, robust standard errors) including treatment and week. The population-averaged odds ratio for good control with the drug versus placebo was 1.97 (95% CI 0.78 to 4.99; p = 0.15).

## Check your understanding

<details class="check"><summary>1. A trial measures a binary response in 200 patients at 4 visits. Analyst A uses ordinary logistic regression on all 800 rows; Analyst B uses GEE with patient as the cluster. Whose SEs are more trustworthy, and in which direction would A be wrong?</summary>

B's. Ignoring the within-patient correlation treats 800 rows as independent, so A's SEs for between-patient effects (such as treatment) are typically too small and the p-values too optimistic. GEE's robust SEs account for the clustering.
</details>

<details class="check"><summary>2. The population-averaged OR for the drug is 1.97 (95% CI 0.78 to 4.99, p = 0.153). A colleague writes "the drug does not improve asthma control." What would you change?</summary>

A non-significant result is not evidence of no effect. The interval ranges from a 22% reduction in the odds to a five-fold increase, so the study (40 patients) is simply inconclusive. Write that the estimate was an OR of 1.97 with a wide CI that includes 1.
</details>

<details class="check"><summary>3. The same data are analysed with a random-intercept logistic model (GLMM) and the drug OR comes out at 2.5 instead of 1.97. Is one of the analyses wrong?</summary>

Not necessarily. They answer different questions. GEE gives the population-averaged OR (comparing populations of treated and untreated patients), while the GLMM gives a subject-specific OR (comparing two patients who share the same random effect). With a logit link the subject-specific OR is further from 1 whenever there is between-patient variation.
</details>

<details class="check"><summary>4. What would happen if patients who do poorly early are more likely to drop out before week 12, and you use a standard GEE?</summary>

Standard GEE is valid only if the missing data are missing completely at random. Dropout that depends on earlier responses (missing at random, MAR) can bias the estimates. Use weighted GEE (inverse-probability weights) or a likelihood-based method such as a mixed model, which is valid under MAR.
</details>

## Common mistakes

- **Interpreting a GEE odds ratio as subject-specific.** It is a marginal (population-averaged) effect. A GLMM odds ratio answers a different, "same patient" question and is usually further from 1.
- **Unsorted data.** `geeglm` requires the rows of a cluster to be adjacent, and time-dependent structures such as AR(1) need visits in time order in every program. Sort by id and time first.
- **Reporting the model-based SEs.** Report the robust (empirical) SEs; they are the point of GEE.
- **Using GEE with few clusters.** With fewer than about 20–30 clusters the sandwich estimator underestimates variance.
- **Assuming it handles dropout like a mixed model.** Standard GEE requires MCAR; if dropout depends on earlier responses, use weighted GEE or a likelihood-based approach (see [mixed models](#/test/mixed-models)).
- **Choosing the working correlation from the p-value.** Pre-specify it; exchangeable is a sensible default for a few visits, AR(1) for many equally spaced visits.
- **Comparing AIC across working correlations.** GEE is not likelihood-based. If you need a criterion, use QIC.

## Related tests

- [Logistic regression](#/test/logistic-regression): the independent-observations version
- [Linear mixed-effects models (MMRM)](#/test/mixed-models): subject-specific / likelihood-based alternative, valid under MAR
- [McNemar's test](#/test/mcnemar): two paired binary measurements
- [Poisson & negative binomial regression](#/test/poisson-regression): count outcomes
- [Cochran–Mantel–Haenszel test](#/test/cmh): stratified binary data without a regression model

## Key takeaways

- GEE fits regression models to correlated data (repeated measures, clusters) and estimates **population-averaged** effects.
- You choose the mean model, the distribution, and a working correlation. The robust (sandwich) SEs protect inference if the working correlation is wrong.
- Report the robust SEs. Ignoring correlation makes SEs too small.
- A GEE OR is marginal. A GLMM OR is subject-specific and usually larger in magnitude.
- Standard GEE needs MCAR missingness and a reasonable number of clusters (a rough guide is dozens, not just a few).

## References

- Liang K-Y, Zeger SL. Longitudinal data analysis using generalized linear models. *Biometrika.* 1986;73:13–22.
- Hardin JW, Hilbe JM. *Generalized Estimating Equations.* 2nd ed. Chapman & Hall/CRC; 2012.
- R documentation: `?geepack::geeglm`; statsmodels documentation: `GEE`
- SAS/STAT User's Guide: The GENMOD Procedure (REPEATED statement)
