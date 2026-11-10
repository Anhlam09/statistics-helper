# Linear mixed-effects models (MMRM)

> Analyses repeated or clustered continuous measurements while using every observed visit, so patients with missing visits stay in the analysis; the mixed model for repeated measures (MMRM) is the standard primary analysis for many longitudinal trial endpoints.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Repeated or clustered continuous measurements, some visits missing | Continuous outcome, several visits per patient | Correct mean and covariance model; normal residuals; missing at random | Change in SBP at weeks 4, 8, 12, with 9 of 24 patients missing week 12 |

## Learning objectives

- Explain why repeated measures on one patient are correlated and why ordinary regression gets the SEs wrong
- Distinguish a random-intercept model from an MMRM with unstructured covariance
- Fit both models in R (and the random-intercept one in Python) and read the week-12 treatment contrast
- Explain what "missing at random" buys you, and why it cannot be tested

## What it is

Repeated measurements on the same patient are correlated, so an ordinary regression that treats every row as independent gives wrong standard errors. A **linear mixed-effects model** (LMM) fits the mean structure (treatment, visit, treatment × visit, baseline) and models the correlation explicitly. This page covers two versions that are often confused:

- **Random-intercept LMM.** Each patient gets a random intercept uᵢ ~ N(0, σ²_u), so all visits of one patient share a patient-specific shift. The implied covariance is *compound symmetry*: equal variance at every visit and equal correlation between any two visits, ρ = σ²_u / (σ²_u + σ²_e) (this ρ is also the ICC).
- **MMRM.** Visit is treated as a categorical factor and the residuals within a patient follow an **unstructured (UN)** covariance matrix, with a free variance for each visit and a free covariance for each pair of visits. There is no random effect; the correlation is modelled directly in the residual covariance ("R-side"). This is the usual regulatory-friendly choice because it makes the fewest assumptions about the correlation.

Both are fitted by (restricted) maximum likelihood, REML, which uses all observed data under the **missing at random (MAR)** assumption.

- **Test statistic:** Wald t (or F) for a linear contrast of the fixed effects, L̂β̂ / SE, with Satterthwaite or Kenward–Roger degrees of freedom
- **Effect size:** the model-based treatment difference at the visit of interest (for example, week 12) with 95% CI; variance components or ICC as secondary summaries

## The intuition

Think of **signal vs noise** with two kinds of noise. Between patients, people differ in their starting level and in how they respond. Within a patient, visits are scattered around that patient's own trajectory. A mixed model separates the two, so that the treatment signal is judged against the right noise. Ignoring the correlation is like counting three visits from one patient as three independent patients: it overstates how much information you have.

Analogy: classroom test scores across the school year. A strong student scores high on every test, so a student's three scores are not three independent pieces of evidence about the class average. The model accounts for that "same student" effect.

The second job of the model is missing data. In this trial, 9 patients have no week-12 value. Look at their week-4 changes in the plot: the eventual dropouts improved less early on (mean about −4.0 mmHg) than the patients who stayed in (mean about −8.6 mmHg). If we analyse only completers, we are comparing the best responders. That is why the complete-case analysis is risky here.

<div class="viz" data-viz="dots" data-unit="week-4 change in SBP, mmHg" data-groups='{"Later dropped out":[3.2,-6.3,-4.9,-4.6,-3.6,-9.0,-4.3,0.9,-7.7],"Completed week 12":[-12.1,-13.4,-8.9,-5.0,-6.1,-15.6,-9.3,0.0,-11.2,-4.6,-4.7,-12.7,-9.2,-12.1,-3.7]}' data-zero="1" data-caption="Patients who later dropped out had improved less at week 4. Dropout is related to earlier response."></div>

A mixed model uses the week-4 and week-8 data from the dropouts, and the correlation between visits, to adjust for this. That works if missingness depends only on *observed* values (MAR).

<div class="viz" data-viz="dots" data-unit="week-12 change in SBP, mmHg (observed completers)" data-groups='{"Placebo":[-21.1,-11.8,-19.0,-1.6,-7.5,-11.6],"Drug":[-18.1,-7.3,-19.5,2.1,-15.3,-2.5,-23.5,-24.3,-10.3]}' data-zero="1" data-caption="Only 15 patients reached week 12. The groups overlap heavily, so the model-based difference of about -4 to -5 mmHg has a wide CI."></div>

**What to notice:** both groups improve, the drug only a little more on average, and the scatter is wide. With so few observed week-12 values the confidence intervals (about −12.7 to +3.0 mmHg) comfortably include zero. This is an example of mechanics, not evidence for or against the drug.

## When to use it

Use a mixed model when **all** of these are true:

- The outcome is **continuous** and approximately normal given the model (change in SBP, HbA1c, FEV₁, a symptom score)
- Patients are measured at **two or more visits**, or are clustered (patients within sites)
- Some visits may be **missing** and you want to keep partially observed patients in the analysis

**Clinical scenario.** A trial randomises 24 patients with hypertension to placebo or a new drug (12 per arm) and records the change in systolic blood pressure (mmHg) from baseline at weeks 4, 8 and 12. Nine patients have no week-12 value because they dropped out, and the dropouts were those who had responded least in the preceding visits. *What is the difference in SBP change between drug and placebo at week 12?*

**Don't use it when:**

- There are only two measurements per person and none are missing → the [paired t-test](#/test/paired-t) is simpler
- There are exactly balanced complete data and a single within-subject factor → [repeated-measures ANOVA](#/test/repeated-measures-anova) gives the same answer, but it drops incomplete patients
- Only the week-12 value matters and nobody has missing data → [ANCOVA](#/test/ancova)
- The outcome is binary or a count and a population-averaged effect is wanted → [GEE](#/test/gee) or [logistic regression](#/test/logistic-regression) with a GLMM
- Dropout depends on the unobserved outcome itself (missing not at random) → use sensitivity analyses such as pattern-mixture or tipping-point models, not a plain MMRM

## Assumptions

1. **Linearity and correct mean model.** The fixed effects (baseline, treatment, visit, treatment × visit) are adequately specified. Check residual-versus-fitted plots.
2. **Residuals (and random effects) are approximately normal.** Check Q–Q plots of the residuals. REML estimates are fairly robust with moderate samples.
3. **The covariance structure is reasonable.** Compound symmetry assumes constant variance and correlation; UN assumes nothing but needs enough patients per visit pair. Compare structures with AIC/BIC (same fixed effects, REML) and inspect the estimated covariance matrix.
4. **Missing data are MAR.** Missingness may depend on *observed* values (earlier visits, baseline, treatment) but not on the unobserved outcome after conditioning on them. This cannot be tested from the data; justify it clinically and run sensitivity analyses.
5. **Patients are independent of each other.** Only the repeated measures within a patient are correlated.

## Hypotheses

For the treatment contrast at week 12 (θ = μ_Drug,12 − μ_Placebo,12, the adjusted mean difference in change from baseline):

- **H₀:** θ = 0
- **H₁:** θ ≠ 0 (two-sided)

## Example data

24 patients, change from baseline in SBP (mmHg; negative = improvement), one row per patient. `.`/`NA` marks a missing visit. Every patient has week 4; 5 patients are missing week 8 and 9 are missing week 12 (monotone dropout). `trt` is 0 = placebo, 1 = drug.

| id | trt | base | w4 | w8 | w12 |
|---|---|---|---|---|---|
| 1 | 0 | 154 | 3.2 | NA | NA |
| 2 | 0 | 141 | −12.1 | −13.6 | −21.1 |
| 3 | 0 | 151 | −13.4 | −17.5 | −11.8 |
| 4 | 0 | 149 | −6.3 | −1.8 | NA |
| 5 | 0 | 145 | −4.9 | NA | NA |
| 6 | 0 | 130 | −8.9 | −18.1 | −19.0 |
| 7 | 0 | 144 | −4.6 | −1.2 | NA |
| 8 | 0 | 142 | −3.6 | NA | NA |
| 9 | 0 | 151 | −9.0 | −0.5 | NA |
| 10 | 0 | 146 | −5.0 | −12.2 | −1.6 |
| 11 | 0 | 147 | −6.1 | −1.8 | −7.5 |
| 12 | 0 | 144 | −15.6 | −8.5 | −11.6 |
| 13 | 1 | 148 | −9.3 | −13.9 | −18.1 |
| 14 | 1 | 148 | 0.0 | −7.0 | −7.3 |
| 15 | 1 | 130 | −11.2 | −15.8 | −19.5 |
| 16 | 1 | 161 | −4.3 | −1.9 | NA |
| 17 | 1 | 155 | −4.6 | −8.2 | 2.1 |
| 18 | 1 | 152 | −4.7 | −22.2 | −15.3 |
| 19 | 1 | 144 | −12.7 | 0.3 | −2.5 |
| 20 | 1 | 156 | 0.9 | NA | NA |
| 21 | 1 | 147 | −7.7 | NA | NA |
| 22 | 1 | 149 | −9.2 | −21.6 | −23.5 |
| 23 | 1 | 139 | −12.1 | −9.4 | −24.3 |
| 24 | 1 | 162 | −3.7 | −5.5 | −10.3 |

## How to run it

**Match the models before you compare the software.** R, Python and SAS give the same numbers only if they fit the *same model*: same fixed effects, same covariance structure, REML (not ML), and the same reference level for the factors. In this page:

- The **random-intercept LMM** is available in all three languages and matches (`lmer`, `statsmodels.MixedLM`, `PROC MIXED ... RANDOM INTERCEPT`). Estimates, variance components and the REML fit agree. Standard errors differ slightly because the denominator degrees of freedom are handled differently (Satterthwaite in lmerTest, Kenward–Roger in the SAS code below, large-sample z in statsmodels).
- The **MMRM with unstructured covariance** is run in R (`nlme::gls`) and SAS (`PROC MIXED ... REPEATED / TYPE=UN`). `statsmodels` has no unstructured residual covariance (its `MixedLM` supports random effects only), so there is no equivalent Python block. The dedicated `mmrm` R package (`mmrm(chg ~ base + trt * visit + us(visit | id), data)`) fits the same model with a faster, more robust optimiser and Kenward–Roger df, and is worth considering for regulatory work.

### R

```r
suppressPackageStartupMessages({library(lmerTest); library(nlme)})
wide <- read.table(text = "
id trt base w4 w8 w12
1 0 154 3.2 NA NA
2 0 141 -12.1 -13.6 -21.1
3 0 151 -13.4 -17.5 -11.8
4 0 149 -6.3 -1.8 NA
5 0 145 -4.9 NA NA
6 0 130 -8.9 -18.1 -19.0
7 0 144 -4.6 -1.2 NA
8 0 142 -3.6 NA NA
9 0 151 -9.0 -0.5 NA
10 0 146 -5.0 -12.2 -1.6
11 0 147 -6.1 -1.8 -7.5
12 0 144 -15.6 -8.5 -11.6
13 1 148 -9.3 -13.9 -18.1
14 1 148 0.0 -7.0 -7.3
15 1 130 -11.2 -15.8 -19.5
16 1 161 -4.3 -1.9 NA
17 1 155 -4.6 -8.2 2.1
18 1 152 -4.7 -22.2 -15.3
19 1 144 -12.7 0.3 -2.5
20 1 156 0.9 NA NA
21 1 147 -7.7 NA NA
22 1 149 -9.2 -21.6 -23.5
23 1 139 -12.1 -9.4 -24.3
24 1 162 -3.7 -5.5 -10.3", header = TRUE)

# Wide -> long; week 12 is the reference visit, so the trt coefficient = week-12 difference
long <- reshape(wide, direction = "long", varying = c("w4","w8","w12"), v.names = "chg",
                timevar = "vnum", times = 1:3, idvar = "id")
long$visit <- relevel(factor(c("w4","w8","w12")[long$vnum]), ref = "w12")
long$trt   <- factor(long$trt, labels = c("Placebo", "Drug"))
long <- long[order(long$id, long$vnum), ]
long <- long[!is.na(long$chg), ]          # drop missing visits (MAR handled by the model)

# A. Random-intercept LMM, REML (lmerTest adds Satterthwaite df and p-values)
m1 <- lmer(chg ~ base + trt * visit + (1 | id), data = long, REML = TRUE)
a <- summary(m1)$coefficients["trtDrug", ]
cat("LMM  Drug-Placebo at wk12:", round(a["Estimate"], 2), " SE", round(a["Std. Error"], 2),
    " 95% CI", round(a["Estimate"] + c(-1, 1) * qt(.975, a["df"]) * a["Std. Error"], 2),
    " p =", round(a["Pr(>|t|)"], 3), "\n")
print(VarCorr(m1))

# B. MMRM: unstructured covariance across visits, REML (nlme::gls)
m2 <- gls(chg ~ base + trt * visit, data = long, method = "REML",
          correlation = corSymm(form = ~ vnum | id),
          weights = varIdent(form = ~ 1 | visit))
b <- summary(m2)$tTable["trtDrug", ]
ci <- intervals(m2)$coef["trtDrug", ]
cat("MMRM Drug-Placebo at wk12:", round(b["Value"], 2), " SE", round(b["Std.Error"], 2),
    " 95% CI", round(ci[c("lower", "upper")], 2), " p =", round(b["p-value"], 3), "\n")
print(round(getVarCov(m2, individual = "2"), 1))   # estimated 3 x 3 covariance of the residuals

# C. Complete-case ANCOVA at week 12 (completers only)
cc <- lm(w12 ~ base + factor(trt), data = wide)
cat("Complete-case n =", nobs(cc), " diff =", round(coef(cc)[3], 2),
    " 95% CI", round(confint(cc)[3, ], 2), "\n")
```

**What the code does**

- `read.table(... NA ...)` reads the wide data; `reshape()` converts to **long** format (one row per patient-visit), and `chg` rows that are `NA` are removed. The model uses every remaining row, so a patient with only week 4 still contributes.
- `relevel(..., ref = "w12")` makes week 12 the reference visit. With a `trt * visit` interaction, the `trtDrug` coefficient is then the drug-minus-placebo difference *at week 12*, not an average over visits.
- `lmer(... + (1 | id))` fits the random-intercept model; `REML = TRUE` is the default. Loading **lmerTest** (rather than plain lme4) adds Satterthwaite df and p-values.
- `gls(...)` fits the MMRM. `corSymm(form = ~ vnum | id)` gives a free correlation for each pair of visits (`vnum` is the visit's position, which also handles missing visits), and `varIdent(form = ~ 1 | visit)` gives each visit its own variance. Together they form an unstructured covariance.
- `getVarCov(m2, individual = "2")` prints the estimated 3 × 3 covariance matrix for a patient with all three visits (patient 2).
- The `lm()` block is a **complete-case** ANCOVA at week 12, kept for comparison only.

**Output**

```text
LMM  Drug-Placebo at wk12: -3.86  SE 3.24  95% CI -10.37 2.64  p = 0.238 
 Groups   Name        Std.Dev.
 id       (Intercept) 3.4875  
 Residual             5.1482  
MMRM Drug-Placebo at wk12: -4.87  SE 3.91  95% CI -12.72 2.98  p = 0.218 
Marginal variance covariance matrix
     [,1] [,2] [,3]
[1,] 18.1  5.1 10.8
[2,]  5.1 50.3 30.5
[3,] 10.8 30.5 60.6
  Standard Deviations: 4.2544 7.0922 7.7846 
Complete-case n = 15  diff = -3.18  95% CI -12.41 6.06 
```

### Python

```python
import io
import pandas as pd
import statsmodels.formula.api as smf

wide = pd.read_csv(io.StringIO("""id trt base w4 w8 w12
1 0 154 3.2 NA NA
2 0 141 -12.1 -13.6 -21.1
3 0 151 -13.4 -17.5 -11.8
4 0 149 -6.3 -1.8 NA
5 0 145 -4.9 NA NA
6 0 130 -8.9 -18.1 -19.0
7 0 144 -4.6 -1.2 NA
8 0 142 -3.6 NA NA
9 0 151 -9.0 -0.5 NA
10 0 146 -5.0 -12.2 -1.6
11 0 147 -6.1 -1.8 -7.5
12 0 144 -15.6 -8.5 -11.6
13 1 148 -9.3 -13.9 -18.1
14 1 148 0.0 -7.0 -7.3
15 1 130 -11.2 -15.8 -19.5
16 1 161 -4.3 -1.9 NA
17 1 155 -4.6 -8.2 2.1
18 1 152 -4.7 -22.2 -15.3
19 1 144 -12.7 0.3 -2.5
20 1 156 0.9 NA NA
21 1 147 -7.7 NA NA
22 1 149 -9.2 -21.6 -23.5
23 1 139 -12.1 -9.4 -24.3
24 1 162 -3.7 -5.5 -10.3"""), sep=r"\s+")

# Wide -> long; drop missing visits (the model uses all observed data)
long = wide.melt(id_vars=["id", "trt", "base"], var_name="visit", value_name="chg").dropna()

# Random-intercept LMM, REML; week 12 is the reference visit,
# so the trt coefficient is the Drug - Placebo difference at week 12
m = smf.mixedlm("chg ~ base + C(trt) * C(visit, Treatment('w12'))",
                data=long, groups=long["id"]).fit(reml=True)
est, se = m.params["C(trt)[T.1]"], m.bse["C(trt)[T.1]"]
print(f"N obs = {len(long)}, N patients = {long['id'].nunique()}")
print(f"LMM Drug-Placebo at wk12: {est:.2f}  SE {se:.2f}  "
      f"95% CI ({est-1.96*se:.2f}, {est+1.96*se:.2f})  p = {m.pvalues['C(trt)[T.1]']:.3f}")
print(f"Random-intercept SD = {m.cov_re.iloc[0,0]**0.5:.3f}, residual SD = {m.scale**0.5:.3f}")

# Complete-case ANCOVA at week 12 for comparison
cc = smf.ols("w12 ~ base + C(trt)", data=wide.dropna(subset=["w12"])).fit()
lo, hi = cc.conf_int().loc["C(trt)[T.1]"]
print(f"Complete-case n = {int(cc.nobs)}  diff = {cc.params['C(trt)[T.1]']:.2f}  95% CI ({lo:.2f}, {hi:.2f})")
```

**What the code does**

- `melt(...).dropna()` converts to long format and removes missing visits; 58 of the 72 possible patient-visits are observed.
- `C(visit, Treatment('w12'))` sets week 12 as the reference level, exactly as in the R code, so `C(trt)[T.1]` is the week-12 treatment difference.
- `smf.mixedlm(..., groups=long["id"])` with no `re_formula` fits a **random intercept** per patient. `.fit(reml=True)` is the default but is written out so the intent is clear.
- statsmodels reports **z-based** Wald p-values and CIs (no small-sample df), so its SE (3.28) differs slightly from lmerTest's Satterthwaite SE (3.24).
- `m.cov_re` is σ²_u and `m.scale` is σ²_e, so the square roots are the SDs that R prints in `VarCorr`.
- There is **no unstructured MMRM** in statsmodels. Do not substitute a random-intercept model and call it an MMRM.

**Output**

```text
N obs = 58, N patients = 24
LMM Drug-Placebo at wk12: -3.86  SE 3.28  95% CI (-10.29, 2.56)  p = 0.239
Random-intercept SD = 3.488, residual SD = 5.148
Complete-case n = 15  diff = -3.18  95% CI (-12.41, 6.06)
```

### SAS

```sas
data bp;
  input id trt base w4 w8 w12;
  datalines;
1 0 154 3.2 . .
2 0 141 -12.1 -13.6 -21.1
3 0 151 -13.4 -17.5 -11.8
4 0 149 -6.3 -1.8 .
5 0 145 -4.9 . .
6 0 130 -8.9 -18.1 -19.0
7 0 144 -4.6 -1.2 .
8 0 142 -3.6 . .
9 0 151 -9.0 -0.5 .
10 0 146 -5.0 -12.2 -1.6
11 0 147 -6.1 -1.8 -7.5
12 0 144 -15.6 -8.5 -11.6
13 1 148 -9.3 -13.9 -18.1
14 1 148 0.0 -7.0 -7.3
15 1 130 -11.2 -15.8 -19.5
16 1 161 -4.3 -1.9 .
17 1 155 -4.6 -8.2 2.1
18 1 152 -4.7 -22.2 -15.3
19 1 144 -12.7 0.3 -2.5
20 1 156 0.9 . .
21 1 147 -7.7 . .
22 1 149 -9.2 -21.6 -23.5
23 1 139 -12.1 -9.4 -24.3
24 1 162 -3.7 -5.5 -10.3
;
run;

/* Wide -> long, keeping only observed visits (visit = 1, 2, 3 = weeks 4, 8, 12) */
data long;
  set bp;
  array w{3} w4 w8 w12;
  do visit = 1 to 3;
    chg = w{visit};
    if not missing(chg) then output;
  end;
  keep id trt base visit chg;
run;

/* A. Random-intercept LMM (REML is the PROC MIXED default) */
proc mixed data=long method=reml;
  class id trt visit;
  model chg = base trt|visit / solution ddfm=kr;
  random intercept / subject=id;
  estimate 'Drug - Placebo, week 12' trt -1 1 trt*visit 0 0 -1 0 0 1 / cl;
run;

/* B. MMRM: unstructured residual covariance across visits */
proc mixed data=long method=reml;
  class id trt visit;
  model chg = base trt|visit / solution ddfm=kr;
  repeated visit / subject=id type=un r rcorr;
  lsmeans trt*visit / diff cl;
  estimate 'Drug - Placebo, week 12' trt -1 1 trt*visit 0 0 -1 0 0 1 / cl;
run;
```

**What the code does**

- The `array` loop reshapes the data to **long** format. `PROC MIXED` needs one row per patient-visit, and rows with a missing `chg` are simply not output, which is what lets the model use partially observed patients.
- `class id trt visit` declares the factors. `visit` must be in `CLASS` for an MMRM so that each visit gets its own mean and covariance parameters. `trt` is numeric 0/1, so its levels sort as placebo, drug.
- `model chg = base trt|visit` expands to `trt visit trt*visit` plus baseline. `ddfm=kr` requests Kenward–Roger df, the usual choice for MMRM in trials; the default methods (containment with `RANDOM`, residual otherwise) ignore the uncertainty in the estimated covariance parameters, which matters with an unstructured covariance and missing data.
- `random intercept / subject=id` makes block A the random-intercept model. `method=reml` is the default and is written out for clarity.
- `repeated visit / subject=id type=un` makes block B the MMRM: a free variance for each visit and a free covariance for each pair. `r` and `rcorr` print the estimated covariance and correlation matrices. Use `type=cs` here to get compound symmetry as a residual structure instead.
- `estimate ... trt -1 1 trt*visit 0 0 -1 0 0 1` is the week-12 drug-minus-placebo contrast. The `trt*visit` cells are ordered (trt, visit) = (0,1), (0,2), (0,3), (1,1), (1,2), (1,3), so the last cell of each arm is week 12. The `/ cl` option gives the 95% CI.
- Expect block A to agree with `lmer` to about two decimals (the estimate and variance components are identical; SE and df differ slightly because of Kenward–Roger), and block B to agree closely with `gls`.

## Interpreting the output

- **Sample.** 58 of 72 possible patient-visits are observed (24 patients; 15 have week 12). All 24 patients contribute to the mixed models.
- **Random-intercept LMM.** The week-12 difference (drug − placebo) is **−3.86 mmHg** (SE 3.24, 95% CI −10.37 to 2.64, p = 0.238). The patient SD is 3.49 and the residual SD is 5.15, so the implied within-patient correlation is 3.49² / (3.49² + 5.15²) = 0.31. Python gives the same estimate and variance components, with z-based SE 3.28 and CI −10.29 to 2.56.
- **MMRM (unstructured).** The week-12 difference is **−4.87 mmHg** (SE 3.91, 95% CI −12.72 to 2.98, p = 0.218). The estimated covariance matrix shows variances that **grow over time** (SD 4.3, 7.1, 7.8 mmHg at weeks 4, 8, 12) and correlations that differ by pair (about 0.17 between weeks 4 and 8, 0.33 between weeks 4 and 12, and 0.55 between weeks 8 and 12). That pattern is why the unstructured model is preferred here: compound symmetry would force equal variances and equal correlations.
- **Complete-case ANCOVA.** Using only the 15 completers, the difference is −3.18 mmHg with a wider CI (−12.41 to 6.06). It discards 9 patients (6 placebo, 3 drug), and it assumes that dropout is unrelated to outcome (completely at random), which is not the case here.
- **Why MMRM copes better with missing visits.** Under MAR, a patient who dropped out after week 4 still provides information on the week-4 response and, through the estimated correlations, on what is likely to have happened later. The MMRM uses that information; complete-case analysis throws it away and can be biased when dropout is related to earlier response. In this example none of the intervals exclude 0, so the study is inconclusive. It is an illustration of the mechanics, not evidence for or against the drug.

**How to report it:**
> Change from baseline in systolic blood pressure was analysed with an MMRM (unstructured covariance; fixed effects for baseline SBP, treatment, visit and treatment-by-visit interaction; REML; Kenward–Roger df). All 24 randomised patients were included under a missing-at-random assumption. At week 12 the adjusted mean difference (drug − placebo) was −4.9 mmHg (95% CI −12.7 to 3.0; p = 0.22).

## Check your understanding

<details class="check"><summary>1. A trial analyses three visits per patient with an ordinary linear regression and treats each row as an independent observation. What goes wrong?</summary>

Rows from the same patient are correlated, so the model thinks it has more independent information than it does. Standard errors are too small, CIs too narrow, and p-values too optimistic. Use a mixed model (or GEE) that accounts for within-patient correlation.
</details>

<details class="check"><summary>2. What is the difference between the random-intercept model and the MMRM on this page?</summary>

The random-intercept model gives each patient a random shift, which implies equal variance at every visit and equal correlation between any two visits (compound symmetry). The MMRM treats visit as a factor and lets every visit have its own variance and every pair its own correlation (unstructured). Here the variances grow over time (SD 4.3, 7.1, 7.8), so the MMRM is the better fit.
</details>

<details class="check"><summary>3. Patients drop out because their blood pressure is rising, in ways not captured by earlier visits. Does the MMRM fix this?</summary>

No. That is missing not at random (MNAR). MMRM assumes the dropout depends only on observed data. MAR cannot be verified from the data, so run sensitivity analyses such as tipping-point or pattern-mixture models.
</details>

<details class="check"><summary>4. The week-12 estimate is −4.9 mmHg with a CI of −12.7 to 3.0. Write one sentence on what you can and cannot conclude.</summary>

The drug may lower SBP by up to about 13 mmHg more than placebo at week 12, or may raise it by up to about 3 mmHg; the study is too small to tell. The result is inconclusive, not evidence that the drug does not work.
</details>

## Common mistakes

- **Calling a random-intercept model an MMRM.** An MMRM uses an unstructured (or other chosen) residual covariance with visit as a categorical factor. They give different answers when variances change over time.
- **Treating visit as a continuous covariate.** This forces a straight-line trend in the mean. Use a categorical visit for the primary analysis.
- **Comparing software output without matching the model.** Differences in REML vs ML, covariance structure, reference level, or df method are the usual reasons that R, SAS and Python disagree.
- **Using ML to compare covariance structures with different fixed effects, or REML to compare different fixed effects.** Compare covariance structures with REML and the same fixed effects; compare fixed effects with ML.
- **Forgetting to pre-specify the covariance structure.** ICH E9 expects the primary model, including its covariance structure and a fallback if it fails to converge, to be fixed in the statistical analysis plan.
- **Assuming MAR is proven.** It is an untestable assumption; add sensitivity analyses (see ICH E9(R1) on estimands and intercurrent events).
- **Reading the `trt` main effect as an overall effect.** With a `trt * visit` interaction, it is the difference at the reference visit only.

## Related tests

- [Repeated-measures ANOVA](#/test/repeated-measures-anova): the older, complete-case version for balanced designs
- [ANCOVA](#/test/ancova): single follow-up time point, adjusted for baseline
- [Paired t-test](#/test/paired-t): two time points on the same patients
- [GEE](#/test/gee): population-averaged alternative, especially for binary and count outcomes
- [Linear regression](#/test/linear-regression): the independent-observations special case

## Key takeaways

- Repeated measures on one patient are correlated; a mixed model accounts for this and gives honest SEs.
- MMRM = visit as a factor + unstructured residual covariance + REML; a random-intercept model is *not* an MMRM.
- It uses all observed visits, so patients with missing visits stay in, if data are **missing at random**.
- MAR is untestable. Justify it clinically and add sensitivity analyses for dropout.
- Report the model-based treatment difference at the key visit with a 95% CI, and state the covariance structure and df method.

## References

- Mallinckrodt CH, Lane PW, Schnell D, et al. Recommendations for the primary analysis of continuous endpoints in longitudinal clinical trials. *Drug Information Journal.* 2008;42:303–319.
- ICH E9(R1). *Addendum on estimands and sensitivity analysis in clinical trials.* 2019.
- Bates D, Mächler M, Bolker B, Walker S. Fitting linear mixed-effects models using lme4. *Journal of Statistical Software.* 2015;67(1). See also R documentation `?nlme::gls`, `?nlme::corSymm`.
- SAS/STAT User's Guide: The MIXED Procedure (REPEATED and RANDOM statements); statsmodels documentation: `MixedLM`.
