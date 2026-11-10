# How to choose a statistical test

> Choosing a test is mostly about your study design, not your data. Answer the five questions below **before** you look at results, and write the answers into your statistical analysis plan (SAP).

## The five questions

### 1. What is the outcome, and what type is it?

| Outcome type | Examples | Typical summary |
|---|---|---|
| **Continuous (≈ normal)** | Systolic BP, HbA1c, LDL-C, FEV₁ | Mean (SD), mean difference |
| **Continuous (skewed) / ordinal** | Length of stay, CRP, pain 0–10, mRS | Median (IQR), Hodges–Lehmann shift |
| **Binary** | Responder yes/no, 30-day mortality, AE occurred | Proportion, risk difference / RR / OR |
| **Categorical (> 2 levels)** | Tumour grade, NYHA class | Counts and percentages |
| **Count / rate** | Exacerbations, hypoglycaemic events | Rate per patient-year, rate ratio |
| **Time-to-event** | Overall survival, time to relapse | KM median, hazard ratio |

### 2. How many groups or conditions are compared?

One group against a reference value, two groups, or three or more.

### 3. Are the observations independent or paired?

**Paired/repeated** means the same subject contributes more than one value: before/after, crossover periods, both eyes, or several visits. Matched case–control pairs also count. Ignoring the pairing is one of the most common errors in the medical literature.

### 4. Do you need to adjust for anything?

Do you need to account for baseline values, stratification factors, or confounders? Then you need a **model** (ANCOVA, regression, Cox, a mixed model) rather than a simple two-group test. In randomized trials, ICH E9 and EMA guidance recommend adjusting for the stratification factors and for important baseline covariates that you pre-specify.

### 5. Are the assumptions reasonable?

Check them with plots (histograms, Q–Q plots, residuals), not just with formal tests. With large samples, the t-test and ANOVA are robust to non-normality because of the central limit theorem. With small samples and clear skew, use a rank-based test or transform the data.

## Quick reference table

| Outcome | 1 group vs reference | 2 independent groups | 2 paired | 3+ independent | 3+ repeated | Adjust for covariates |
|---|---|---|---|---|---|---|
| Continuous, normal | [One-sample t](#/test/one-sample-t) | [Welch t-test](#/test/independent-t) | [Paired t](#/test/paired-t) | [One-way ANOVA](#/test/one-way-anova) | [Mixed model](#/test/mixed-models) / [RM-ANOVA](#/test/repeated-measures-anova) | [ANCOVA](#/test/ancova), [linear regression](#/test/linear-regression) |
| Skewed / ordinal | [Signed-rank](#/test/wilcoxon-signed-rank) | [Mann–Whitney](#/test/mann-whitney) | [Wilcoxon signed-rank](#/test/wilcoxon-signed-rank) | [Kruskal–Wallis](#/test/kruskal-wallis) | [Friedman](#/test/friedman) | [Proportional-odds regression](#/test/ordinal-regression) |
| Binary | [Exact binomial](#/test/one-proportion) | [χ²](#/test/chi-square-independence) / [two-proportion z](#/test/two-proportion-z) / [Fisher](#/test/fisher-exact) | [McNemar](#/test/mcnemar) | [χ² r × c](#/test/chi-square-independence) | [GEE](#/test/gee) | [Logistic regression](#/test/logistic-regression), [CMH](#/test/cmh) |
| Count / rate | [Poisson](#/test/poisson-regression) | [Poisson / NB](#/test/poisson-regression) | [GEE](#/test/gee) | [Poisson / NB](#/test/poisson-regression) | [GEE](#/test/gee) | [Poisson / NB regression](#/test/poisson-regression) |
| Time-to-event | [Kaplan–Meier](#/test/kaplan-meier) | [Log-rank](#/test/log-rank) | Stratified Cox | [Log-rank](#/test/log-rank) | — | [Cox PH](#/test/cox-ph); [competing risks](#/test/competing-risks) |

**Other questions:**

- Agreement between raters → [Cohen's kappa](#/test/cohens-kappa) (categorical) · [ICC](#/test/icc) (continuous)
- Agreement between two measurement methods → [Bland–Altman](#/test/bland-altman)
- Diagnostic accuracy → [ROC / AUC](#/test/roc-auc)
- Equivalence or non-inferiority → [TOST](#/test/tost-equivalence)
- Association between two variables → [Pearson](#/test/pearson) · [Spearman](#/test/spearman)
- Many endpoints or many arms → [multiplicity adjustment](#/test/multiple-comparisons)
- Another event (e.g., death) prevents the event of interest → [competing risks](#/test/competing-risks)
- Planning a study → [sample size: means & proportions](#/test/sample-size-means-proportions) · [sample size: time-to-event](#/test/sample-size-survival) · [calculators](#/calculators)

## Parametric and nonparametric equivalents

Parametric tests compare means and assume a distribution, usually normal. Rank-based (nonparametric) tests make weaker assumptions but answer a slightly different question: they compare whole distributions or a location shift rather than means.

| Design | Parametric | Rank-based alternative | Robust parametric option |
|---|---|---|---|
| 1 group vs a reference value | [One-sample t](#/test/one-sample-t) | [Wilcoxon signed-rank](#/test/wilcoxon-signed-rank) (one-sample) | — |
| 2 independent groups | [Student t](#/test/independent-t) | [Mann–Whitney U](#/test/mann-whitney) | [Welch t](#/test/independent-t) |
| 2 paired measurements | [Paired t](#/test/paired-t) | [Wilcoxon signed-rank](#/test/wilcoxon-signed-rank) | — |
| 3+ independent groups | [One-way ANOVA](#/test/one-way-anova) | [Kruskal–Wallis](#/test/kruskal-wallis) | [Welch's ANOVA](#/test/welch-anova) |
| 3+ repeated measurements | [Repeated-measures ANOVA](#/test/repeated-measures-anova) | [Friedman](#/test/friedman) | [Mixed model](#/test/mixed-models) |
| Association of two variables | [Pearson r](#/test/pearson) | [Spearman ρ](#/test/spearman) | — |
| Ordinal outcome with covariates | — | — | [Proportional-odds regression](#/test/ordinal-regression) |

**Switching to a rank-based test is reasonable when:**

- the outcome is ordinal (a pain score, a Likert item, the mRS)
- the sample is small *and* the plots show clear skew or outliers
- the outcome has floor or ceiling effects (many zeros, values piled at the maximum)

**It is not needed just because:**

- **A Shapiro–Wilk test is "significant" in a large sample.** With n in the hundreds, trivial departures from normality reach p < 0.05, and the t-test is still fine.
- **The raw data are skewed when the test actually assumes normal residuals or differences.** The paired t-test, for example, only needs the *differences* to be roughly normal.

A third option is often better than either: log-transform a right-skewed lab value (CRP, triglycerides, viral load), then report the ratio of geometric means.

## Assumptions checklist

| Assumption | Applies to | How to check | If it fails |
|---|---|---|---|
| **Independent observations** | Every test | Study design: one row per subject? Repeated visits, both eyes, or clustered sites? | Paired tests, [mixed models](#/test/mixed-models), [GEE](#/test/gee) |
| **Approximate normality** of residuals or differences | t-tests, ANOVA, linear regression, Pearson | Q–Q plot and histogram first; [Shapiro–Wilk](#/test/shapiro-wilk) as a supplement | Transform, use a rank-based test, or rely on the CLT with a large n |
| **Equal variances** | Student t, classic ANOVA | Compare SDs (a ratio > 2 is a warning sign); [Levene / Brown–Forsythe](#/test/levene) | [Welch t](#/test/independent-t), [Welch's ANOVA](#/test/welch-anova) |
| **Expected counts ≥ 5** | χ² tests | Look at the *expected* table, not the observed one | [Fisher's exact](#/test/fisher-exact), exact binomial |
| **Linearity** | Pearson, linear and logistic regression (continuous predictors) | Scatterplot; residuals vs fitted; splines | Transform, add splines, or use [Spearman](#/test/spearman) |
| **No severe multicollinearity** | Multiple regression | VIF (> 5 is a concern, > 10 is serious) | Drop or combine predictors |
| **Enough events per predictor** | Logistic and Cox regression | Count the events, not the patients (rough guide: ≥ 10 per parameter) | Fewer predictors, penalized regression |
| **Proportional hazards** | Cox, log-rank | Schoenfeld residuals; log(−log S) plots | Stratify, time-varying effects, RMST |
| **Proportional odds** | Ordinal regression | Score/Brant test; compare separate binary logits | Partial proportional odds |
| **Non-informative censoring** | Kaplan–Meier, log-rank, Cox | Design: why do patients drop out? | Sensitivity analyses; [competing-risks](#/test/competing-risks) methods if the "censoring" is really another event |

**A note on "n ≥ 30":** the rule that n ≥ 30 makes everything normal is a rough heuristic, not a law. Heavily skewed data, such as costs or length of stay, can need far more than 30, while symmetric data need fewer. Likewise, "use a z-test when n ≥ 30" is outdated: if the SD is estimated from the data, the [t-test](#/test/one-sample-t) is the correct test at any sample size (see the [z-test page](#/test/one-sample-z)).

## How to use each test page

Every test page follows the same structure, so you can jump straight to what you need:

1. **What it is**: the idea in plain words and the effect size you should report.
2. **When to use it**: a checklist, a realistic clinical scenario, and when *not* to use it, with links to the right alternative.
3. **Assumptions**: what must hold and how to check it.
4. **Hypotheses**: H₀ and H₁ in symbols and words.
5. **Example data**: a small dataset that is used in all three languages.
6. **How to run it**: tabs for **R** and **Python** (plus **SAS** on many pages). Your language choice is remembered. Each tab has copy-ready code, a line-by-line explanation, and the real output.
7. **Interpreting the output**: what every number means, plus a sentence you can adapt for a manuscript.
8. **Common mistakes** and **related tests**.

**Tip:** run the code yourself. The R and Python examples are self-contained, so you can paste them into RStudio, Jupyter, or VS Code and get the same numbers shown on the page.

## Software setup

**R** (≥ 4.2). Base R covers most tests. Add these packages as needed:

```r
install.packages(c("survival", "lme4", "lmerTest", "car", "emmeans",
                   "irr", "pROC", "geepack", "DescTools", "TOSTER"))
```

**Python** (≥ 3.10):

```bash
pip install numpy pandas scipy statsmodels pingouin lifelines scikit-learn matplotlib
```

**SAS** (optional): where SAS is shown, it uses SAS/STAT procedures (TTEST, NPAR1WAY, FREQ, GLM, MIXED, GENMOD, LOGISTIC, LIFETEST, PHREG). SAS OnDemand for Academics is free for learners.

## Good practice in clinical research

- **Pre-specify.** Write down the primary endpoint, the test, the covariates, and how you'll handle missing data in the protocol and SAP *before* unblinding (ICH E9).
- **Estimates over p-values.** Report the effect size with its 95% CI. A p-value alone says nothing about clinical relevance.
- **Don't choose the test from the results.** Switching from a t-test to Mann–Whitney because it "gives a smaller p" inflates the false-positive rate.
- **"Not significant" ≠ "no difference."** To claim similarity, use an [equivalence or non-inferiority](#/test/tost-equivalence) design.
- **Mind multiplicity.** Multiple endpoints, arms, or interim looks need a [pre-planned adjustment](#/test/multiple-comparisons).
- **Define the estimand.** ICH E9(R1) asks you to define how intercurrent events, such as treatment discontinuation or rescue medication, are handled. That choice can matter more than which test you pick.

## Glossary

- **α (significance level)**: The false-positive rate you accept *before* the study, usually 0.05 (two-sided).
- **Censoring**: Follow-up ended before the event was observed (the study ended, or the patient was lost to follow-up). The patient's time still counts up to that point.
- **Confidence interval (95% CI)**: A range of effect sizes compatible with the data. If the study were repeated many times, 95% of such intervals would contain the true effect.
- **Effect size**: The magnitude of a difference or association in meaningful units: a mean difference in mmHg, a risk difference, an OR, RR, or HR, or a standardized measure such as Cohen's d. Report it with its CI.
- **Estimand**: A precise description of the treatment effect a trial targets, including how intercurrent events such as rescue medication are handled (ICH E9(R1)).
- **Hazard ratio (HR)**: The ratio of instantaneous event rates between groups. HR = 0.70 means a 30% lower event rate at any given time, assuming proportional hazards.
- **Independent vs paired**: Independent groups contain different subjects. Paired (or repeated) data come from the same subject, or from deliberately matched subjects.
- **MCID**: The minimal clinically important difference: the smallest change patients would notice. It is often used as δ in [sample-size calculations](#/test/sample-size-means-proportions).
- **Null hypothesis (H₀)**: The "no effect" statement the test tries to reject, such as μ₁ = μ₂ or HR = 1.
- **Odds ratio vs risk ratio**: The RR compares probabilities and the OR compares odds. They are similar only when the outcome is rare (< ~10%). Logistic regression reports ORs.
- **p-value**: How surprising the data would be *if H₀ were true*. It is not the probability that H₀ is true, and it says nothing about the size of the effect.
- **Parametric vs nonparametric**: Parametric tests assume a distribution and estimate its parameters (such as means). Nonparametric (rank-based) tests make fewer assumptions and usually compare distributions.
- **Power (1 − β)**: The probability of detecting a true effect of a given size. It is planned before the study, usually at 80–90%. Never compute "observed power" after the fact.
- **Type I and Type II error**: A false positive (rejecting a true H₀; rate α) and a false negative (missing a real effect; rate β).
