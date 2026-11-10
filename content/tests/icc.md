# Intraclass correlation coefficient (ICC)

> Quantifies the reliability of continuous measurements made by several raters or repeated on the same subjects, as the share of total variance that is due to real differences between patients.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Quantify reliability of repeated or multi-rater continuous measurements | Continuous, same subjects measured 2+ times or by 2+ raters | Random sample of patients; independent raters; roughly normal; balanced design | Tumour diameter measured by 3 radiologists on 20 CT scans |

## Learning objectives

- Explain the ICC as the share of total variance due to real differences between patients
- Choose the right ICC form: one-way, two-way, agreement vs consistency, single vs average
- Run the ICC in R and Python and report it with the form and 95% CI
- Explain why the ICC depends on the spread of the sample, and why p-values for it are rarely useful

## What it is

When several raters (or repeated runs of a device) measure the same patients, the total variance of the readings can be split into **true between-patient variance**, **rater (systematic) variance**, and **random error**. The ICC is the proportion of total variance attributable to the patients:

ICC = σ²_patients / (σ²_patients + σ²_raters + σ²_error)

(for the forms that count rater differences as error). It ranges from 0 (the measurements say nothing about which patient is which) to 1 (raters would give identical values). Unlike Pearson correlation, the ICC can be pulled down by raters who are consistently higher or lower than the others, if you choose a form that counts that as error.

There are many ICCs, and choosing the right one is the main task. The Shrout–Fleiss labels (ICC(1,1), ICC(2,1), ICC(3,1), and the "k" versions) and the McGraw–Wong labels (one-way, two-way random/mixed; consistency vs absolute agreement; single vs average) describe the same six quantities:

| Shrout–Fleiss | McGraw–Wong name | Raters are… | Counts rater bias as error? | Pingouin label |
|---|---|---|---|---|
| ICC(1,1) | One-way random, single | a different random set for each patient | yes | ICC(1,1) |
| ICC(2,1) | Two-way random, **absolute agreement**, single | a random sample of possible raters | **yes** | ICC(A,1) |
| ICC(3,1) | Two-way mixed, **consistency**, single | the only raters of interest (fixed) | **no** | ICC(C,1) |
| ICC(1,k), ICC(2,k), ICC(3,k) | the same models, **average** of k raters | as above | as above | ICC(1,k), ICC(A,k), ICC(C,k) |

**Which one should I choose?**

1. **Will you use a single rater's reading or the average of k raters in practice?** Single → the "1" forms. Mean of k raters → the "k" forms. Most clinical applications use a single reading, so ICC(2,1) or ICC(3,1) is the usual answer.
2. **Are the same raters rating every patient?** Yes → a two-way model (2 or 3). No (different raters for different patients) → one-way, ICC(1,1).
3. **Do you want raters to be interchangeable in the readings' actual values (absolute agreement) or only to rank patients the same way (consistency)?** If a reading from rater A must be directly comparable to one from rater B (the usual clinical case), use **absolute agreement**, ICC(2,1). If you will only compare patients rated by the same rater, or can correct for rater offset, **consistency**, ICC(3,1), is enough.
4. **Do the raters represent a larger pool you want to generalise to?** Yes → random (2). Only these raters (for example the central readers of a trial) → mixed/fixed (3).

When in doubt for an inter-rater study, report **ICC(2,1) (two-way random, absolute agreement, single rater)** with its 95% CI, and state the form explicitly.

- **Test statistic:** F = MS_patients / MS_error (two-way) tests H₀: ICC = 0, rarely informative
- **Effect size:** the ICC with its 95% CI. Koo & Li's conventional bands for the *lower* CI limit are < 0.50 poor, 0.50–0.75 moderate, 0.75–0.90 good, > 0.90 excellent.

## The intuition

Reliability is a signal-to-noise idea. The **signal** is how much patients truly differ from one another (a tumour of 15 mm versus 85 mm). The **noise** is everything that makes repeat readings of the *same* patient disagree: random measurement error and, depending on the form chosen, systematic differences between raters. The ICC is signal ÷ (signal + noise). If patients differ hugely and raters barely disagree, the ICC is near 1. If raters disagree as much as patients differ, it is near 0.

An everyday analogy is a bathroom scale. A scale that reads 70.1, 70.0, and 70.2 kg on repeated weighings of the same person is reliable *for telling that person from someone who weighs 90 kg*. A scale that wanders by ±10 kg is not, even if it averages right.

Look at the three raters' readings below. Within each rater the 20 tumour sizes spread from about 15 to 90 mm (a wide between-patient spread, SD about 16 mm for rater 1). Between raters the dot clouds are almost on top of each other, with just a small shift (means 63.1, 64.9, and 68.5 mm).

<div class="viz" data-viz="dots" data-unit="mm, tumour diameter" data-groups='{"Rater 1":[52.4,79.2,74.8,69,55.9,74.2,45.5,52.8,69.1,82.5,15.3,54.9,57.1,56.6,63.9,64.9,72.9,84.9,62.2,74.6],"Rater 2":[55.4,73,66.8,73.3,67.6,75,46.5,63.1,65.3,81.5,21.5,62.6,69.7,63.6,56.3,61,73.7,83.5,70.2,68],"Rater 3":[60,78.2,68.7,77.5,66.4,80.3,51.7,61.5,67.4,89.7,23.2,60.8,65.1,65.7,67.3,81.5,76.4,90.5,57.9,79.6]}' data-caption="Between-patient spread (signal) is large compared with the rater-to-rater shift (noise), which is why the ICC is high."></div>

Rater 3 reads about 5 mm larger than rater 1 on average. That small, consistent offset is what lowers ICC(2,1) = 0.884 relative to ICC(3,1) = 0.911: absolute agreement counts the offset as noise; consistency ignores it.

The F test of H₀: ICC = 0 compares between-patient and error mean squares. Here F = 31.6 with 19 and 38 degrees of freedom. Under H₀ the F curve has its bulk near 1, so ours is far out in the tail.

<div class="viz" data-viz="f" data-stat="31.59" data-df1="19" data-df2="38"></div>

What to notice: the F statistic is enormous, so the test rejects ICC = 0 easily. But "ICC is not zero" is a very low bar. What matters is the size of the ICC and how well the CI pins it down.

## When to use it

Use the ICC when **all** of these are true:

- The measurement is **continuous** (tumour diameter in mm, bone density, a lab value, a 0–100 score)
- The same patients are measured by **two or more** raters, devices, or repeated occasions (test–retest)
- The question is **reliability**: can the instrument or raters distinguish between patients consistently?

**Clinical scenario.** Three radiologists independently measure the longest diameter (mm) of a lung tumour on the same 20 CT scans. A trial intends to use a **single radiologist's** reading for each patient. *How reliable is a single reading, and does it matter that one radiologist tends to measure larger?*

**Don't use it when:**

- The rating is categorical → use [Cohen's kappa](#/test/cohens-kappa)
- You want the size and direction of the difference between two methods in clinical units → use [Bland–Altman analysis](#/test/bland-altman)
- You only want to know whether two variables are linearly related → use [Pearson correlation](#/test/pearson)
- The sample is very small (n < 15–20 subjects) or has very little spread between patients, since the ICC is then unstable and low
- Different samples have very different spreads and you plan to compare ICCs across studies: the ICC depends on the between-patient variance of the sample

## Assumptions

1. **Measurements are continuous and roughly normally distributed**, with the same error variance for all raters. Check histograms of the readings or the residuals.
2. **Patients are a random sample** of the population you want to describe; the ICC is specific to the spread of patients measured.
3. **Raters work independently** and are blinded to one another's measurements.
4. **A complete, balanced design** (every rater measures every patient) for the standard ANOVA estimators. With missing readings, use a mixed-model-based ICC.
5. **No strong interaction between patients and raters** (raters do not rank patients differently). Plot each rater's values against the patient mean to check.

## Hypotheses

- **H₀:** ICC = 0 (readings carry no information about which patient is measured)
- **H₁:** ICC > 0

The hypothesis test is rarely useful. Report the ICC with its CI and judge the lower limit against the reliability you need.

## Example data

Tumour diameter (mm) for 20 patients measured by the same 3 radiologists. Rater means are 63.1, 64.9 and 68.5 mm, so rater 3 reads about 5 mm larger than rater 1.

| Patient | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| Rater 1 | 52.4 | 79.2 | 74.8 | 69 | 55.9 | 74.2 | 45.5 | 52.8 | 69.1 | 82.5 |
| Rater 2 | 55.4 | 73 | 66.8 | 73.3 | 67.6 | 75 | 46.5 | 63.1 | 65.3 | 81.5 |
| Rater 3 | 60 | 78.2 | 68.7 | 77.5 | 66.4 | 80.3 | 51.7 | 61.5 | 67.4 | 89.7 |

| Patient | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 |
|---|---|---|---|---|---|---|---|---|---|---|
| Rater 1 | 15.3 | 54.9 | 57.1 | 56.6 | 63.9 | 64.9 | 72.9 | 84.9 | 62.2 | 74.6 |
| Rater 2 | 21.5 | 62.6 | 69.7 | 63.6 | 56.3 | 61 | 73.7 | 83.5 | 70.2 | 68 |
| Rater 3 | 23.2 | 60.8 | 65.1 | 65.7 | 67.3 | 81.5 | 76.4 | 90.5 | 57.9 | 79.6 |

## How to run it

### R

```r
library(irr)

# Tumour diameter (mm) measured on 20 patients by the same 3 radiologists
dat <- read.table(text = "
patient r1 r2 r3
1 52.4 55.4 60
2 79.2 73 78.2
3 74.8 66.8 68.7
4 69 73.3 77.5
5 55.9 67.6 66.4
6 74.2 75 80.3
7 45.5 46.5 51.7
8 52.8 63.1 61.5
9 69.1 65.3 67.4
10 82.5 81.5 89.7
11 15.3 21.5 23.2
12 54.9 62.6 60.8
13 57.1 69.7 65.1
14 56.6 63.6 65.7
15 63.9 56.3 67.3
16 64.9 61 81.5
17 72.9 73.7 76.4
18 84.9 83.5 90.5
19 62.2 70.2 57.9
20 74.6 68 79.6
", header = TRUE)
x <- dat[, c("r1", "r2", "r3")]          # patients in rows, raters in columns

# Shrout-Fleiss / McGraw-Wong forms
run <- function(model, type, unit) {
  r <- icc(x, model = model, type = type, unit = unit)
  c(ICC = r$value, lower = r$lbound, upper = r$ubound, F = r$Fvalue)
}
res <- rbind(
  "ICC(1,1)  one-way, single"       = run("oneway",  "agreement",   "single"),
  "ICC(2,1)  absolute agreement"    = run("twoway",  "agreement",   "single"),
  "ICC(3,1)  consistency"           = run("twoway",  "consistency", "single"),
  "ICC(2,k)  abs. agreement, mean"  = run("twoway",  "agreement",   "average"),
  "ICC(3,k)  consistency, mean"     = run("twoway",  "consistency", "average"))
print(round(res, 3))
```

**What the code does**

- `x` is a patients × raters matrix. `irr::icc()` expects **rows = subjects, columns = raters**.
- `model = "oneway"` or `"twoway"` chooses the one-way or two-way ANOVA model. In a two-way model, `type = "agreement"` counts rater offsets as error (Shrout–Fleiss ICC(2,·)), while `type = "consistency"` ignores them (ICC(3,·)).
- `unit = "single"` is the reliability of one rater's reading; `"average"` is the reliability of the mean of the k raters.
- `type` is ignored in a one-way model, which always counts everything except patient variance as error.
- `r$lbound` and `r$ubound` are the 95% confidence limits (`conf.level = 0.95` is the default).
- The extra `F` column is MS_patients / MS_error for the corresponding model.

**Output**

```text
                                 ICC lower upper      F
ICC(1,1)  one-way, single      0.883 0.776 0.948 23.715
ICC(2,1)  absolute agreement   0.884 0.743 0.952 31.585
ICC(3,1)  consistency          0.911 0.823 0.961 31.585
ICC(2,k)  abs. agreement, mean 0.958 0.895 0.983 31.585
ICC(3,k)  consistency, mean    0.968 0.933 0.986 31.585
```

### Python

```python
import io
import pandas as pd
import pingouin as pg

# Tumour diameter (mm) measured on 20 patients by the same 3 radiologists
dat = pd.read_csv(io.StringIO("""
patient r1 r2 r3
1 52.4 55.4 60
2 79.2 73 78.2
3 74.8 66.8 68.7
4 69 73.3 77.5
5 55.9 67.6 66.4
6 74.2 75 80.3
7 45.5 46.5 51.7
8 52.8 63.1 61.5
9 69.1 65.3 67.4
10 82.5 81.5 89.7
11 15.3 21.5 23.2
12 54.9 62.6 60.8
13 57.1 69.7 65.1
14 56.6 63.6 65.7
15 63.9 56.3 67.3
16 64.9 61 81.5
17 72.9 73.7 76.4
18 84.9 83.5 90.5
19 62.2 70.2 57.9
20 74.6 68 79.6
"""), sep=r"\s+")

# pingouin wants long format: one row per (patient, rater) measurement
long = dat.melt(id_vars="patient", var_name="rater", value_name="diameter")

icc = pg.intraclass_corr(data=long, targets="patient", raters="rater", ratings="diameter")
print(icc.round(3).to_string(index=False))
```

**What the code does**

- `pingouin.intraclass_corr` needs **long** format, one row per patient-rater measurement, with columns for `targets` (patients), `raters`, and `ratings`.
- It returns all six forms at once. `ICC(A,1)` is ICC(2,1), **absolute agreement**; `ICC(C,1)` is ICC(3,1), **consistency**; the "k" rows are the average-of-k versions. (Older pingouin versions call these ICC2 and ICC3.)
- `CI95` is the 95% confidence interval, rounded to two decimals.
- The values are the same as in R (the R table leaves out ICC(1,k), which is rarely used).

**Output**

```text
    Type   ICC      F  df1  df2  pval         CI95
ICC(1,1) 0.883 23.715   19   40   0.0 [0.78, 0.95]
ICC(A,1) 0.884 31.585   19   38   0.0 [0.74, 0.95]
ICC(C,1) 0.911 31.585   19   38   0.0 [0.82, 0.96]
ICC(1,k) 0.958 23.715   19   40   0.0 [0.91, 0.98]
ICC(A,k) 0.958 31.585   19   38   0.0  [0.9, 0.98]
ICC(C,k) 0.968 31.585   19   38   0.0 [0.93, 0.99]
```

### SAS

```sas
data tumour;
  input patient r1 r2 r3;
  datalines;
1 52.4 55.4 60
2 79.2 73 78.2
3 74.8 66.8 68.7
4 69 73.3 77.5
5 55.9 67.6 66.4
6 74.2 75 80.3
7 45.5 46.5 51.7
8 52.8 63.1 61.5
9 69.1 65.3 67.4
10 82.5 81.5 89.7
11 15.3 21.5 23.2
12 54.9 62.6 60.8
13 57.1 69.7 65.1
14 56.6 63.6 65.7
15 63.9 56.3 67.3
16 64.9 61 81.5
17 72.9 73.7 76.4
18 84.9 83.5 90.5
19 62.2 70.2 57.9
20 74.6 68 79.6
;
run;

/* Long format: one row per patient-rater measurement */
data long;
  set tumour;
  array r{3} r1-r3;
  do rater = 1 to 3;
    y = r{rater};
    output;
  end;
  keep patient rater y;
run;

/* ICC(2,1): patients and raters both random; absolute agreement */
proc mixed data=long method=reml;
  class patient rater;
  model y = ;
  random patient rater;
  ods output covparms=cp_abs;
run;

/* ICC(3,1): raters fixed; consistency */
proc mixed data=long method=reml;
  class patient rater;
  model y = rater;
  random patient;
  ods output covparms=cp_con;
run;

/* ICC = patient variance / (sum of the variance components in the model) */
data icc;
  set cp_abs(in=a) cp_con(in=b);
  if estimate > 0;
  length form $12;
  form = ifc(a, 'ICC(2,1)', 'ICC(3,1)');
run;

proc sql;
  select form, sum(case when lowcase(covparm)='patient' then estimate else 0 end)
              / sum(estimate) as icc format=6.3
  from icc group by form;
quit;
```

**What the code does**

- The first `PROC MIXED` fits `y = μ + patient + rater + error` with both effects random (REML). Its variance components give ICC(2,1) = σ²_patient / (σ²_patient + σ²_rater + σ²_error).
- The second declares `rater` as a **fixed** effect (`model y = rater`) so that systematic rater differences are removed from the error and only patient and residual variance remain. ICC(3,1) = σ²_patient / (σ²_patient + σ²_error).
- `ods output covparms=` saves the variance-component table. The final `DATA` step and `PROC SQL` divide the patient variance by the sum of the model's components, giving one ICC per form.
- For a **balanced** design with positive variance components, REML estimates equal the ANOVA estimates, so these ICCs equal the Shrout–Fleiss values from R and Python (about 0.884 and 0.911 for these data). Estimates truncated at zero or unbalanced data make the two approaches differ.
- This approach gives point estimates only. For confidence intervals you can derive F-based limits from the ANOVA mean squares (as `irr` does), bootstrap, or use a published SAS macro such as `%INTRACC`; if you use a macro, read its documentation to confirm which Shrout–Fleiss form it computes, and compare its output to R or Python before relying on it.

## Interpreting the output

- **ICC(2,1) = 0.884** (95% CI 0.743 to 0.952). The point estimate falls in the "good" band (0.75–0.90), but the lower limit of 0.74 is just short of the 0.75 cut-off, so judged by the lower limit the reliability of a single reading is **moderate to good**.
- **ICC(3,1) = 0.911** (95% CI 0.823 to 0.961). It is higher than ICC(2,1) because it ignores the systematic offset between raters (rater 3 reads about 5 mm larger than rater 1). The gap between the two (0.911 vs 0.884) is the part of the unreliability caused by rater bias, and shows that the radiologists are not fully interchangeable in absolute terms.
- **ICC(1,1) = 0.883** (0.776 to 0.948) is similar to ICC(2,1) here because rater offsets are small compared with patient-to-patient spread.
- **Average-of-3 forms** ICC(2,3) = 0.958 and ICC(3,3) = 0.968: averaging the three readings reduces random error and raises reliability. Use these only if the trial will actually average three readers.
- **F = 31.6** (df 19, 38) is the test of ICC = 0 for the two-way models, which is clearly rejected; this is not informative on its own.
- Python reproduces all of these (for example ICC(A,1) = 0.884, CI 0.74 to 0.95).
- For a trial using one reader per patient, the relevant number is ICC(2,1) = 0.88.

**How to report it:**
> Inter-rater reliability of tumour diameter was assessed in 20 patients measured by 3 radiologists. Using a two-way random-effects, absolute-agreement, single-rater model, ICC(2,1) was 0.88 (95% CI 0.74 to 0.95), indicating moderate-to-good reliability (lower confidence limit 0.74).

## Check your understanding

<details class="check"><summary>1. Two thermometers always differ by exactly 1.0 degrees (one reads higher). A study reports Pearson r = 1.00 and concludes they agree. Is that right?</summary>

No. Pearson r ignores systematic offset: perfectly parallel readings give r = 1. An absolute-agreement ICC would be pulled down by the constant 1.0-degree difference, which is why it is the better reliability measure when readings must be interchangeable.
</details>

<details class="check"><summary>2. A trial will use a single radiologist per patient, but the reliability study reports only ICC(2,3) = 0.958. Is that the right number?</summary>

No. ICC(2,3) is the reliability of the *average of three* raters. For a single-reader trial you need the single-rater form, ICC(2,1) = 0.884 (CI 0.74 to 0.95). Reporting the average-of-k value overstates reliability.
</details>

<details class="check"><summary>3. In a cohort restricted to patients with tumours of 55 to 65 mm, the same radiologists give ICC(2,1) = 0.35. Did the radiologists get worse?</summary>

No. Same raters, same error, but the between-patient variance (the signal) shrank because the patients are now almost identical, so the ratio falls. The ICC describes a measurement *in a population with a given spread*. Report the SD or SEM of the measurements too.
</details>

<details class="check"><summary>4. What if rater 3's readings were all 20 mm larger than they are now (a bigger constant offset)?</summary>

ICC(3,1), the consistency form, would be unchanged, because it ignores constant offsets. ICC(2,1), absolute agreement, would drop sharply because the offset becomes a large part of the total variance. This is a useful way to see what each form measures.
</details>

## Common mistakes

- **Not saying which ICC was used.** "ICC = 0.88" is meaningless without the model, type and unit. State, for example, "two-way random, absolute agreement, single rater".
- **Using consistency when absolute agreement is needed.** If raters differ systematically, ICC(3,1) hides it. Use ICC(2,1) when readings from different raters must be directly comparable.
- **Reporting the average-of-k ICC for a study that will use single readings.** It is always higher and will overstate the reliability of one reader.
- **Interpreting the ICC without its CI** or using the conventional bands on the point estimate; use the lower CI limit.
- **Using Pearson correlation as a reliability measure.** It ignores rater offset and scale differences.
- **Comparing ICCs across populations** with different between-patient spread. A homogeneous sample gives a low ICC even for a precise instrument; report the SD of the measurements or the SEM too.
- **Using the ICC for agreement between two methods in clinical units.** Combine it with [Bland–Altman](#/test/bland-altman) limits of agreement, which show how far apart the readings may be.

## Related tests

- [Cohen's kappa](#/test/cohens-kappa): agreement for categorical ratings
- [Bland–Altman analysis](#/test/bland-altman): bias and limits of agreement between two methods
- [Pearson correlation](#/test/pearson): linear association, not agreement
- [One-way ANOVA](#/test/one-way-anova): the ANOVA calculations that underlie the ICC
- [Linear mixed-effects models](#/test/mixed-models): variance-components estimation, and ICCs for unbalanced or clustered data

## Key takeaways

- ICC = between-patient variance ÷ total variance: the share of variation in readings that reflects real patient differences.
- Always state the form (for example, two-way random, absolute agreement, single rater) and report the 95% CI. Judge reliability by the lower limit.
- Absolute agreement counts rater bias as error. Consistency ignores it. Single-rater and average-of-k ICCs answer different questions.
- The ICC depends on the spread of patients in the sample, so it cannot be compared across populations without context.
- The p-value for ICC = 0 is rarely informative. Do not use Pearson correlation as a reliability measure.

## References

- Shrout PE, Fleiss JL. Intraclass correlations: uses in assessing rater reliability. *Psychological Bulletin.* 1979;86:420–428.
- McGraw KO, Wong SP. Forming inferences about some intraclass correlation coefficients. *Psychological Methods.* 1996;1:30–46.
- Koo TK, Li MY. A guideline of selecting and reporting intraclass correlation coefficients for reliability research. *Journal of Chiropractic Medicine.* 2016;15:155–163.
- R documentation: `?irr::icc`; pingouin documentation: `intraclass_corr`; SAS/STAT User's Guide: The MIXED Procedure (COVPARMS table)
