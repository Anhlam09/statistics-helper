# Cochran–Mantel–Haenszel test

> Tests the association between a binary treatment and a binary outcome while controlling for a stratifying factor such as study site, and gives a single combined (Mantel–Haenszel) odds ratio.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Binary treatment vs binary outcome, stratified by site or another factor | Binary outcome in K 2 × 2 tables | Independent patients; similar odds ratio across strata | Psoriasis PASI 75 response, drug vs placebo, 4 sites |

## Learning objectives

- Explain why pooling strata can mislead (Simpson's paradox) and how CMH avoids it
- Run the CMH test and read the Mantel–Haenszel common odds ratio and its CI
- Use the Breslow–Day test and the site-specific ORs to judge whether one common OR makes sense
- Report an odds ratio as an odds ratio, not as a risk ratio

## What it is

In a multi-centre trial, response rates and the patient mix differ between sites. Pooling all patients into a single 2 × 2 table can mislead (Simpson's paradox), and analysing each site separately wastes information. The CMH approach keeps one 2 × 2 table per site (stratum), compares the observed number of responders in the treatment arm with the number expected under no association **within each site**, and combines the differences across sites. The **Mantel–Haenszel common odds ratio** is a weighted average of the site-specific odds ratios, which assumes that the treatment effect is about the same in every site. The **Breslow–Day test** checks that assumption.

- **Test statistic:** CMH χ² with 1 df, built from (Σ[observed − expected] over strata)², scaled by the sum of the hypergeometric variances
- **Effect size:** the **Mantel–Haenszel common odds ratio** with a 95% CI (Robins–Breslow–Greenland variance), plus the site-specific odds ratios
- **Homogeneity check:** **Breslow–Day** χ² with K − 1 df for K strata

## The intuition

A multi-centre trial is like several small trials run in parallel. Sites differ: some enrol sicker patients, so their baseline response rate is lower. CMH compares drug and placebo **within** each site (like for like), then combines the within-site evidence. The signal is the excess of responders on the drug over what would be expected if the drug did nothing, summed across sites; the noise is the combined sampling variance of those counts. Because site is held fixed, differences between sites cannot masquerade as a treatment effect.

Think of comparing two runners' times: you only compare times set on the same course, because a hilly course would otherwise confuse the comparison. Here each site is a course.

The dot plot shows the drug-minus-placebo response difference at each site, in percentage points. Each site is one dot, and the dashed line at 0 means "no benefit".

<div class="viz" data-viz="dots" data-unit="percentage points (drug minus placebo)" data-groups='{"Site differences":[26.7,29.2,23.7,23.0]}' data-zero="1" data-caption="Notice that all four sites lie well above zero and close to each other (23 to 29 points). This consistency is what lets CMH combine the sites into one common odds ratio. If one site pointed below zero, a single common OR would hide the story."></div>

Combining the four sites gives CMH χ² = 14.85 with 1 degree of freedom. Under H₀ (no association within sites), the statistic follows the curve below, and 14.85 sits far into the upper tail.

<div class="viz" data-viz="chisq" data-stat="14.85" data-df="1"></div>

## When to use it

Use the CMH test when **all** of these are true:

- The outcome and the exposure (treatment) are both **binary**
- The data fall into **K strata** defined by a categorical factor (site, sex, disease stage, centre), each giving a 2 × 2 table
- You want a **stratified** test of association and a **common** effect estimate

**Clinical scenario.** A multi-centre trial randomises 231 patients with moderate psoriasis to a new drug or placebo at 4 sites. The outcome is PASI 75 response at week 16 (yes/no). Randomisation was stratified by site, so the pre-specified primary analysis is stratified by site. *Does the drug increase the response rate after accounting for site, and by how much?*

**Don't use it when:**

- There is only one stratum (or you do not care about site) → use the [two-proportion z-test](#/test/two-proportion-z) or [Fisher's exact test](#/test/fisher-exact)
- Counts within strata are tiny and sparse → use the exact conditional test, or [logistic regression](#/test/logistic-regression) with caution
- You need to adjust for continuous covariates or several factors at once → use [logistic regression](#/test/logistic-regression)
- Treatment effects clearly go in different directions across strata (strong qualitative interaction) → a single common OR is misleading; report strata separately
- The data are paired → use [McNemar's test](#/test/mcnemar)

## Assumptions

1. **Observations are independent**, and each patient belongs to one stratum. Check from the design.
2. **The odds ratio is homogeneous across strata** (the treatment effect is similar at each site). Check with the **Breslow–Day test** and by inspecting the site-specific ORs. The CMH test itself is still valid when ORs differ in magnitude but go in the same direction, but the common OR then describes only an average.
3. **Strata are defined before analysis**, ideally by the randomisation scheme (ICH E9).
4. **Sample size is adequate in aggregate.** The test remains reliable with small strata as long as the totals are large (that is its advantage over per-stratum chi-square tests). A common rough guide asks for the summed expected counts to be at least 5, but this is not a strict threshold; with very sparse tables use an exact conditional test.

## Hypotheses

- **H₀:** the odds ratio is 1 in every stratum (conditional independence of treatment and response given site)
- **H₁:** the common odds ratio is not 1 (the drug is associated with response, consistently across sites)

For the Breslow–Day test: **H₀:** all site-specific ORs are equal; **H₁:** at least one differs.

## Example data

Responders / patients in each arm, by site.

| Site | Drug | Placebo | Site OR |
|---|---|---|---|
| A | 20 / 30 (66.7%) | 12 / 30 (40.0%) | 3.00 |
| B | 15 / 24 (62.5%) | 8 / 24 (33.3%) | 3.33 |
| C | 24 / 38 (63.2%) | 15 / 38 (39.5%) | 2.63 |
| D | 12 / 23 (52.2%) | 7 / 24 (29.2%) | 2.65 |
| **Total** | 71 / 115 (61.7%) | 42 / 116 (36.2%) | |

## How to run it

### R

```r
# 2 x 2 x 4 array: Treatment x Response x Site
# fill order within a site: Drug-Yes, Placebo-Yes, Drug-No, Placebo-No
dat <- array(c(20, 12, 10, 18,    # Site A
               15,  8,  9, 16,    # Site B
               24, 15, 14, 23,    # Site C
               12,  7, 11, 17),   # Site D
             dim = c(2, 2, 4),
             dimnames = list(Treatment = c("Drug", "Placebo"),
                             Response  = c("Yes", "No"),
                             Site      = c("A", "B", "C", "D")))

# CMH test and Mantel-Haenszel common odds ratio with 95% CI
mantelhaen.test(dat, correct = FALSE)

# Breslow-Day test of homogeneity of the odds ratios
library(DescTools)
BreslowDayTest(dat, correct = FALSE)

# Site-specific odds ratios
apply(dat, 3, function(t) (t[1, 1] * t[2, 2]) / (t[1, 2] * t[2, 1]))
```

**What the code does**

- `mantelhaen.test()` takes a 2 × 2 × K array (or three factors). The **third dimension is the stratum**.
- **`correct = TRUE` is the default** and applies a continuity correction to the CMH χ². SAS and statsmodels do not correct by default, so set `correct = FALSE` to match them.
- It also prints the **Mantel–Haenszel common odds ratio** and its 95% CI. The odds ratio is for the first row and first column relative to the second, here Drug vs Placebo for Yes.
- `DescTools::BreslowDayTest()` is not in base R. `correct = FALSE` gives the standard Breslow–Day test; `TRUE` adds Tarone's correction.
- The array fills column by column, so each block of four numbers is read as Drug-Yes, Placebo-Yes, Drug-No, Placebo-No. Check the layout by printing `dat`.

**Output**

```text
	Mantel-Haenszel chi-squared test without continuity correction

data:  dat
Mantel-Haenszel X-squared = 14.851, df = 1, p-value = 0.0001164
alternative hypothesis: true common odds ratio is not equal to 1
95 percent confidence interval:
 1.673116 4.891188
sample estimates:
common odds ratio 
         2.860686 

	Breslow-Day test on Homogeneity of Odds Ratios

data:  dat
X-squared = 0.11962, df = 3, p-value = 0.9894

       A        B        C        D 
3.000000 3.333333 2.628571 2.649351 
```

### Python

```python
import numpy as np
from statsmodels.stats.contingency_tables import StratifiedTable

# One 2x2 table per site: [[Drug Yes, Drug No], [Placebo Yes, Placebo No]]
tables = [np.array([[20, 10], [12, 18]]),   # Site A
          np.array([[15,  9], [ 8, 16]]),   # Site B
          np.array([[24, 14], [15, 23]]),   # Site C
          np.array([[12, 11], [ 7, 17]])]   # Site D
st = StratifiedTable(tables)

# CMH test (no continuity correction)
t = st.test_null_odds(correction=False)
print(f"CMH chi2 = {t.statistic:.4f}, df = 1, p = {t.pvalue:.5f}")

# Mantel-Haenszel common odds ratio with 95% CI
lo, hi = st.oddsratio_pooled_confint()
print(f"MH common OR = {st.oddsratio_pooled:.4f}, 95% CI ({lo:.4f}, {hi:.4f})")

# Breslow-Day homogeneity test (adjust=False; True gives Tarone's version)
bd = st.test_equal_odds(adjust=False)
print(f"Breslow-Day chi2 = {bd.statistic:.4f}, df = 3, p = {bd.pvalue:.4f}")

# Site-specific odds ratios
print("Site ORs:", [round(float(a[0, 0] * a[1, 1] / (a[0, 1] * a[1, 0])), 3) for a in tables])
```

**What the code does**

- `StratifiedTable` takes a list of 2 × 2 arrays, each laid out with the exposure in rows and the outcome in columns.
- `test_null_odds(correction=False)` is the CMH test. The default is also no correction, so it matches `mantelhaen.test(correct = FALSE)`.
- `oddsratio_pooled` is the Mantel–Haenszel OR, and `oddsratio_pooled_confint()` gives the RBG confidence interval, the same one R reports.
- `test_equal_odds(adjust=False)` is the standard Breslow–Day test; `adjust=True` applies Tarone's adjustment.
- `StratifiedTable.from_data(exposure, outcome, strata, data)` builds the tables directly from patient-level data.

**Output**

```text
CMH chi2 = 14.8506, df = 1, p = 0.00012
MH common OR = 2.8607, 95% CI (1.6731, 4.8912)
Breslow-Day chi2 = 0.1196, df = 3, p = 0.9894
Site ORs: [3.0, 3.333, 2.629, 2.649]
```

### SAS

```sas
data psor;
  input site $ trt $ response $ count;
  datalines;
A Drug    Yes 20
A Drug    No  10
A Placebo Yes 12
A Placebo No  18
B Drug    Yes 15
B Drug    No   9
B Placebo Yes  8
B Placebo No  16
C Drug    Yes 24
C Drug    No  14
C Placebo Yes 15
C Placebo No  23
D Drug    Yes 12
D Drug    No  11
D Placebo Yes  7
D Placebo No  17
;
run;

proc freq data=psor order=data;
  weight count;
  tables site*trt*response / cmh nopercent norow nocol;
run;
```

**What the code does**

- The `tables` statement with **three variables** treats the first as the stratification factor (`site`), the second as the row variable (`trt`), and the third as the column variable (`response`).
- `cmh` produces the **Cochran–Mantel–Haenszel statistics** (general association, row mean scores, and nonzero correlation). For a 2 × 2 × K table, the **nonzero correlation** line, with 1 df, is the usual CMH test for a binary outcome.
- For 2 × 2 × K tables, the output also includes the **"Estimates of the Common Relative Risk"** table: the **Mantel–Haenszel and logit common odds ratios and relative risks with 95% CIs**, followed by the **Breslow–Day test for homogeneity of the odds ratios**. Add the `bdt` option to the `tables` statement to also get Tarone's adjusted version.
- `order=data` puts Drug first and Yes first, so the odds ratio is for Drug vs Placebo with Yes as the event. If SAS orders alphabetically (No before Yes), the OR would be inverted.
- `nopercent norow nocol` keeps the per-site tables short.
- The SAS CMH statistic does not apply a continuity correction.

## Interpreting the output

- **Site-specific ORs: 3.00, 3.33, 2.63, 2.65.** All are well above 1 and similar. The response rate is higher on the drug at every site (by 23 to 29 percentage points).
- **CMH χ²(1) = 14.85, p = 0.0001.** After controlling for site, the drug is associated with a higher response rate.
- **Mantel–Haenszel common OR = 2.86 (95% CI 1.67 to 4.89).** The odds of response are about 2.9 times higher on the drug, with data compatible with a true common OR between about 1.7 and 4.9. This is an odds ratio, not a risk ratio: because response is common, the crude risk ratio is smaller (61.7% / 36.2% = 1.7).
- **Breslow–Day χ²(3) = 0.12, p = 0.99.** There is no sign that the odds ratio differs by site, which supports using a single common OR. The Breslow–Day test has low power with few strata or small tables, so also look at the site-specific ORs and their consistency of direction.
- **R and Python agree** to four decimals on the CMH χ², the common OR, its CI, and the Breslow–Day test. (Without `correct = FALSE`, R's default would give a slightly smaller χ².)

**How to report it:**
> After stratification by site, response at week 16 was more frequent with drug than placebo (61.7% vs 36.2%; Mantel–Haenszel common odds ratio 2.86, 95% CI 1.67 to 4.89; CMH χ²(1) = 14.85, p < 0.001). The odds ratio did not differ significantly between sites (Breslow–Day p = 0.99).

⚠️ **Important for trials:** in a stratified randomised trial, the pre-specified analysis should respect the stratification factor. The common OR is an *odds* ratio; if the sponsor wants a risk difference, use a stratified (Mantel–Haenszel or Miettinen–Nurminen) risk-difference estimator instead. Do not choose between a stratified and an unstratified analysis based on which gives the smaller p-value.

## Check your understanding

<details class="check"><summary>1. Randomisation in a trial was stratified by site and baseline disease severity (mild / severe). Which test is the natural primary analysis for a binary response?</summary>

A stratified analysis: CMH over the 8 site-by-severity strata, or logistic regression with those factors as covariates. The analysis should reflect the stratification used in the randomisation, as ICH E9 recommends. Be careful if some strata have very few patients (sparse tables).
</details>

<details class="check"><summary>2. The Breslow–Day test gives p = 0.99. A reviewer says "this proves the effect is identical at every site." Is that right?</summary>

No. It shows no evidence of heterogeneity, but the test has low power with few strata and small tables. The better check is to look at the site-specific ORs (3.00, 3.33, 2.63, 2.65) and their consistency of direction and magnitude. Here they are similar, which is what supports a common OR.
</details>

<details class="check"><summary>3. The Mantel–Haenszel common OR is 2.86. A press release says "patients on the drug were 2.9 times as likely to respond." What is wrong?</summary>

An odds ratio is not a risk ratio. Response was common (61.7% on the drug vs 36.2% on placebo), so the risk ratio is only about 1.7 (0.617 / 0.362). The OR of 2.86 overstates it. The accurate statement is "the odds of response were about 2.9 times higher", or better, report the response rates or the risk difference (about 25 percentage points).
</details>

<details class="check"><summary>4. What would happen to the CMH χ² if the drug had worked at sites A and B but had no effect at C and D (ORs near 1 there)?</summary>

The CMH χ² would shrink (the combined excess of responders is smaller) and the Breslow–Day test would be more likely to flag heterogeneity. A common OR near 1.9 would then describe no site well. When effects differ materially between strata, report them separately and treat the pooled estimate with caution.
</details>

## Common mistakes

- **Pooling strata into one 2 × 2 table** when the response rate or allocation varies by site. This can hide, exaggerate, or even reverse the effect.
- **Ignoring the Breslow–Day test, or over-reading it.** A non-significant result does not prove homogeneity, and a significant one does not by itself invalidate the CMH test; look at the direction and size of each stratum's OR.
- **Reporting the common OR when the effects clearly go in opposite directions.**
- **Letting R apply a continuity correction by default** while SAS and statsmodels do not. Set `correct = FALSE` to compare.
- **Getting the array dimension order wrong.** The stratum must be the third dimension in R, and the first variable in a SAS three-way `tables` statement.
- **Interpreting the common OR as a risk ratio.** The OR overstates the RR when the outcome is common (here response is 40–60%).
- **Using too many sparse strata.** The MH estimator copes with sparse strata better than most, but with very sparse data consider exact methods.

## Related tests

- [Two-proportion z-test](#/test/two-proportion-z): unstratified comparison of two proportions
- [Chi-square test of independence](#/test/chi-square-independence): unstratified association
- [Fisher's exact test](#/test/fisher-exact): small single table
- [Logistic regression](#/test/logistic-regression): adjusts for site plus other covariates
- [McNemar's test](#/test/mcnemar): paired binary data

## Key takeaways

- CMH tests drug vs placebo **within strata** (sites) and combines the evidence, avoiding Simpson's paradox.
- The Mantel–Haenszel common OR is a weighted average of the stratum ORs. It assumes they are similar.
- Check homogeneity with the site-specific ORs first, and with Breslow–Day second (it has low power).
- A common odds ratio is not a risk ratio: with common outcomes the OR overstates the RR.
- Respect the randomisation strata in the pre-specified analysis, and do not pick stratified vs unstratified by the smaller p.

## References

- Mantel N, Haenszel W. Statistical aspects of the analysis of data from retrospective studies of disease. *J Natl Cancer Inst.* 1959;22:719–748.
- Agresti A. *Categorical Data Analysis.* 3rd ed. Wiley; 2013. Ch. 6.
- ICH E9: Statistical Principles for Clinical Trials (1998)
- R documentation: `?mantelhaen.test`; statsmodels documentation: `StratifiedTable`
- SAS/STAT User's Guide: The FREQ Procedure — CMH option
