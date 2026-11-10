# Multiple-comparison corrections (Bonferroni, Holm, BH)

> Adjusts p-values when many hypotheses are tested at once, so that the overall chance of false positives (family-wise error rate) or the expected share of false discoveries (false discovery rate) stays under control.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Many hypotheses tested at once (endpoints, subgroups, pairwise comparisons) | A set of raw p-values from any tests | Valid raw p-values; family defined in advance | 10 secondary endpoints in one phase 3 trial |

## Learning objectives

- Explain why testing many hypotheses inflates the chance of a false positive
- Distinguish **family-wise error rate (FWER)** from **false discovery rate (FDR)**
- Run Bonferroni, Holm, and Benjamini–Hochberg adjustments in R and Python
- Report adjusted p-values and pre-specify the multiplicity strategy

## What it is

If you test one true null hypothesis at α = 0.05, you have a 5% chance of a false positive. Test 10 independent true nulls and the chance of **at least one** false positive is 1 − 0.95¹⁰ ≈ 40%. Multiplicity corrections adjust the p-values (or the significance thresholds) to account for the number of tests, m. Three methods cover most clinical needs:

- **Bonferroni:** adjusted p = min(1, m · p). Controls the **family-wise error rate (FWER)**, the probability of any false positive. Simple and valid under any dependence, but conservative.
- **Holm (step-down Bonferroni):** sort the p-values p₍₁₎ ≤ … ≤ p₍ₘ₎ and compute (m − i + 1) · p₍ᵢ₎, enforcing monotonicity. Also controls the FWER under any dependence, is **never less powerful than Bonferroni**, and so should replace it.
- **Benjamini–Hochberg (BH):** adjusted p₍ᵢ₎ = p₍ᵢ₎ · m / i (made monotone). Controls the **false discovery rate (FDR)**: among the findings declared significant, the expected proportion that are false. It is more powerful than the FWER methods, and appropriate for exploratory work with many endpoints (biomarkers, genomics).
- **Test statistic and effect size:** these methods take the raw p-values as input; there is no new statistic. Report the **estimates and CIs** (adjusted where possible) along with the adjusted p-values.

## The intuition

The core idea is **more looks, more flukes**. If you buy one lottery ticket, winning is unlikely. If you buy ten, the chance that at least one wins is much higher. Likewise, each endpoint tested at α = 0.05 is one ticket in the "false positive" lottery. With ten true nulls, P(at least one false positive) = 1 − 0.95¹⁰ ≈ 40%.

Corrections respond by making each p-value harder to beat. Bonferroni multiplies every p by m = 10. Holm is a smarter step-down version of the same idea. Benjamini–Hochberg is more lenient because it only promises to control the *share* of false findings among the discoveries.

The plot below shows the ten raw p-values and their adjusted versions on the same 0 to 1 scale. Seven raw p-values lie below 0.05. After Holm, only one does. After BH, three do. Notice how the adjusted values (especially Holm, the stricter of the two shown) are pushed to the right, past the 0.05 line. Bonferroni is not drawn; its values (10 × raw, capped at 1) would sit even further right.

<div class="viz" data-viz="dots" data-unit="p-value" data-groups='{"Raw":[0.001,0.008,0.013,0.021,0.034,0.041,0.049,0.18,0.42,0.77],"Holm":[0.01,0.072,0.104,0.147,0.204,0.205,0.205,0.54,0.84,0.84],"BH":[0.01,0.04,0.0433,0.0525,0.068,0.0683,0.07,0.225,0.4667,0.77]}' data-caption="Same ten endpoints, three versions of the p-value. Count how many points fall left of 0.05 in each row: 7, 1, and 3."></div>

The cost of protection is **power**. Stricter control (FWER) means fewer false alarms but also more missed real effects. FDR control (BH) trades a few expected false leads for more discoveries, which suits exploratory screening but usually not a regulatory claim.

## When to use it

Use a multiplicity correction when **all** of these are true:

- You are testing **several hypotheses** that would each be interpreted on their own (multiple endpoints, several subgroups, several pairwise comparisons, many biomarkers)
- You have a set of **raw p-values** from any test
- A claim of "at least one significant result" (FWER) or a list of "discoveries" (FDR) is the goal

**Clinical scenario.** A phase 3 trial reports 10 pre-specified secondary endpoints (SBP, DBP, LDL-C, HbA1c, weight, eGFR, CRP, 6-minute walk distance, quality of life, NT-proBNP), each tested with its own two-sample test. Seven have raw p < 0.05. *How many remain significant after accounting for ten looks?*

**Don't use it when:**

- You want all pairwise group comparisons after ANOVA → use [Tukey's HSD](#/test/tukey-hsd), which exploits the correlation structure
- You are comparing several arms with one control → use [Dunnett's test](#/test/dunnett), which is more powerful than Bonferroni for this family
- There is **one** pre-specified primary endpoint → no correction is needed for it; correct only the secondary family, preferably with a pre-specified hierarchical (gatekeeping) strategy
- The endpoints are tested in a fixed order (hierarchical testing) → no p-value adjustment is needed, but the order must be pre-specified

## Assumptions

1. **The raw p-values are valid**, that is, uniformly distributed under their true null hypotheses. Check that each underlying test's assumptions hold.
2. **The family of hypotheses is defined in advance** (the number m is the number of tests you *performed*, including those you chose not to report).
3. **Bonferroni and Holm:** no assumption on dependence between tests.
4. **BH:** valid for independent or positively dependent tests (common for correlated clinical endpoints); the Benjamini–Yekutieli variant (`"BY"`) is valid under any dependence but more conservative.
5. **Which error rate you control matches the purpose.** FWER for confirmatory claims (regulatory); FDR for exploratory screening.

## Hypotheses

For each of the m tests i:

- **H₀ᵢ:** no effect for endpoint i
- **H₁ᵢ:** an effect for endpoint i

The **global** goal of Bonferroni and Holm is P(reject at least one true H₀ᵢ) ≤ α. The goal of BH is E[false discoveries / total discoveries] ≤ q.

## Example data

Raw p-values for 10 secondary endpoints:

| Endpoint | SBP | DBP | LDL | HbA1c | Weight | eGFR | CRP | 6MWD | QoL | NT-proBNP |
|---|---|---|---|---|---|---|---|---|---|---|
| Raw p | 0.001 | 0.008 | 0.013 | 0.021 | 0.034 | 0.041 | 0.049 | 0.180 | 0.420 | 0.770 |

## How to run it

### R

```r
# Raw p-values for 10 secondary endpoints in one trial (each from its own test)
endpoint <- c("SBP","DBP","LDL","HbA1c","Weight","eGFR","CRP","6MWD","QoL","NT-proBNP")
p <- c(0.001, 0.008, 0.013, 0.021, 0.034, 0.041, 0.049, 0.180, 0.420, 0.770)

res <- data.frame(endpoint, raw = p,
  bonferroni = p.adjust(p, method = "bonferroni"),  # min(1, m * p)
  holm       = p.adjust(p, method = "holm"),        # step-down, FWER
  BH         = p.adjust(p, method = "BH"))          # step-up, FDR
print(res, digits = 3)

# How many endpoints are "significant" at 0.05 after each adjustment?
colSums(res[, -1] < 0.05)
```

**What the code does**

- `p.adjust(p, method = ...)` returns **adjusted p-values** in the original order. You compare these directly with 0.05; there is no separate threshold to track.
- `"bonferroni"` multiplies each p by m = 10 and caps at 1. `"holm"` is the default method of `p.adjust` and is the one recommended for FWER control. `"BH"` (alias `"fdr"`) controls the FDR; `"hochberg"` and `"BY"` are also available.
- The vector must contain **all** p-values in the family: `p.adjust` uses `n = length(p)` by default, and dropping the non-significant ones would understate m.
- `colSums(res[, -1] < 0.05)` counts how many endpoints pass each rule.

**Output**

```text
    endpoint   raw bonferroni  holm     BH
1        SBP 0.001       0.01 0.010 0.0100
2        DBP 0.008       0.08 0.072 0.0400
3        LDL 0.013       0.13 0.104 0.0433
4      HbA1c 0.021       0.21 0.147 0.0525
5     Weight 0.034       0.34 0.204 0.0680
6       eGFR 0.041       0.41 0.205 0.0683
7        CRP 0.049       0.49 0.205 0.0700
8       6MWD 0.180       1.00 0.540 0.2250
9        QoL 0.420       1.00 0.840 0.4667
10 NT-proBNP 0.770       1.00 0.840 0.7700
       raw bonferroni       holm         BH 
         7          1          1          3 
```

### Python

```python
import pandas as pd
from statsmodels.stats.multitest import multipletests

# Raw p-values for 10 secondary endpoints in one trial (each from its own test)
endpoint = ["SBP","DBP","LDL","HbA1c","Weight","eGFR","CRP","6MWD","QoL","NT-proBNP"]
p = [0.001, 0.008, 0.013, 0.021, 0.034, 0.041, 0.049, 0.180, 0.420, 0.770]

res = pd.DataFrame({"endpoint": endpoint, "raw": p})
for label, method in [("bonferroni", "bonferroni"), ("holm", "holm"), ("BH", "fdr_bh")]:
    reject, p_adj, _, _ = multipletests(p, alpha=0.05, method=method)
    res[label] = p_adj
print(res.round(3).to_string(index=False))

# Number significant at 0.05 after each adjustment
print((res[["raw", "bonferroni", "holm", "BH"]] < 0.05).sum())
```

**What the code does**

- `multipletests(pvals, alpha, method)` returns four things: a boolean `reject` array, the **adjusted p-values** `p_adj`, and the corrected α for Šidák and Bonferroni. Keep the second.
- Method names differ from R: `"bonferroni"`, `"holm"`, and `"fdr_bh"` (Benjamini–Hochberg; `"fdr_by"` is the arbitrary-dependence version). `"holm"` here is Holm's step-down; `"hs"` is Holm–Šidák.
- The adjusted values are identical to R's `p.adjust`; they are only rounded differently in the print (R shows 0.0525, pandas shows 0.052 after `round(3)`).
- `alpha` only affects the `reject` flags, not the adjusted p-values.

**Output**

```text
 endpoint   raw  bonferroni  holm    BH
      SBP 0.001        0.01 0.010 0.010
      DBP 0.008        0.08 0.072 0.040
      LDL 0.013        0.13 0.104 0.043
    HbA1c 0.021        0.21 0.147 0.052
   Weight 0.034        0.34 0.204 0.068
     eGFR 0.041        0.41 0.205 0.068
      CRP 0.049        0.49 0.205 0.070
     6MWD 0.180        1.00 0.540 0.225
      QoL 0.420        1.00 0.840 0.467
NT-proBNP 0.770        1.00 0.840 0.770
raw           7
bonferroni    1
holm          1
BH            3
dtype: int64
```

### SAS

```sas
/* One row per endpoint; the p-value variable MUST be named Raw_P */
data pvals;
  length endpoint $ 9;                 /* default $ length is 8: would cut NT-proBNP */
  input endpoint $ Raw_P;
  datalines;
SBP       0.001
DBP       0.008
LDL       0.013
HbA1c     0.021
Weight    0.034
eGFR      0.041
CRP       0.049
6MWD      0.180
QoL       0.420
NT-proBNP 0.770
;
run;

/* Adjust the p-values: Bonferroni, Holm (stepdown Bonferroni), Benjamini-Hochberg */
proc multtest inpvalues=pvals bonferroni stepbon fdr out=adjusted;
run;

proc print data=adjusted noobs;
  var endpoint Raw_P bon_p stpbon_p fdr_p;
run;
```

**What the code does**

- `inpvalues=pvals` tells PROC MULTTEST to read p-values from a data set instead of running tests itself. The data set must contain a variable called **`Raw_P`**; other variables (here `endpoint`) are carried through to the output.
- `bonferroni` requests Bonferroni adjustment; `stepbon` requests **Holm's step-down Bonferroni**; `fdr` requests **Benjamini–Hochberg FDR** control. (`hochberg` is also available for the step-up FWER method.)
- `out=adjusted` saves a data set with the adjusted p-values in `bon_p`, `stpbon_p`, and `fdr_p`. These correspond to R's `"bonferroni"`, `"holm"`, and `"BH"`.
- PROC MULTTEST can also compute the raw p-values itself (for example from `test mean(...)` or `ca`/`trend` tests with a `class` variable) and adjust them, or use permutation or bootstrap resampling; the `inpvalues=` route is the one to use when p-values come from elsewhere.
- The `length endpoint $ 9;` statement is needed because list input with a bare `$` gives character variables a default length of 8, which would cut `NT-proBNP`.

## Interpreting the output

- **Raw p < 0.05:** 7 of 10 endpoints (SBP through CRP). If these were reported unadjusted, the chance of at least one false positive could be very high (up to 40% with ten independent null endpoints).
- **Bonferroni:** adjusted p = 10 × raw, capped at 1. Only **SBP** survives (0.001 → 0.010). DBP becomes 0.08.
- **Holm:** only **SBP** survives (0.010). The Holm adjusted p-values are smaller than or equal to Bonferroni's for every endpoint (DBP 0.072 vs 0.080, LDL 0.104 vs 0.130), but in this example not by enough to change the decisions. The rule is "never worse", not "always different".
- **Benjamini–Hochberg:** SBP (0.010), DBP (0.040), and LDL (0.043) are significant at an FDR of 5%; HbA1c (0.0525) just misses. The interpretation is different: of the endpoints declared significant, we expect about 5% to be false discoveries, not that the chance of any false positive is 5%.
- **Choosing between them:** the FWER methods fit confirmatory claims (labelling, regulatory submissions); BH fits exploratory screening where some false leads are acceptable.

**How to report it:**
> Of the 10 secondary endpoints, systolic blood pressure remained statistically significant after Holm adjustment for multiplicity (adjusted p = 0.010); the other nine did not (adjusted p ≥ 0.072). Under Benjamini–Hochberg control of the false discovery rate at 5%, diastolic blood pressure (adjusted p = 0.040) and LDL-C (adjusted p = 0.043) were also flagged as exploratory findings.

⚠️ **Important for trials:** the multiplicity strategy must be pre-specified in the protocol or statistical analysis plan (ICH E9; see also FDA's guidance on multiple endpoints). Choosing the method after seeing which one gives more significant results defeats its purpose.

## Check your understanding

<details class="check"><summary>1. A trial has one pre-specified primary endpoint and ten secondary endpoints. Do you need to adjust the primary endpoint's p-value for the ten secondaries?</summary>

No. The primary endpoint is tested on its own at α = 0.05. Multiplicity adjustment (or a pre-specified hierarchy) applies to the secondary family so that claims about them stay credible.
</details>

<details class="check"><summary>2. After BH adjustment at 5% FDR, three endpoints are significant. A colleague says "there is only a 5% chance any of these three is a false positive." Correct?</summary>

No. FDR 5% means that, on average, about 5% of the declared discoveries are expected to be false, which is about 0.15 of three findings on average. It is not a 5% chance that any single finding, or any finding at all, is false. The probability of *any* false positive is what FWER methods (Bonferroni, Holm) control.
</details>

<details class="check"><summary>3. You ran 10 tests but only wish to report the 3 that have p &lt; 0.05. What m should you use?</summary>

m = 10. The family is every test you performed, whether or not you report it. Using m = 3 understates the multiplicity and defeats the correction.
</details>

<details class="check"><summary>4. What if SBP's raw p-value were 0.006 instead of 0.001? Would it survive Holm adjustment at 0.05?</summary>

The smallest p-value is multiplied by m = 10 under Holm, giving 0.060, which is above 0.05. So it would no longer be significant. That shows how much a single threshold gets tightened: with 10 tests, the smallest p-value must be below 0.005 to survive.
</details>

## Common mistakes

- **Adjusting only the "significant" p-values or dropping tests from the family.** m must include every test in the family, including the ones not reported.
- **Using Bonferroni when Holm is available.** Holm has the same guarantee and is never less powerful.
- **Mixing up FWER and FDR.** An FDR of 5% does not mean a 5% chance of any false positive.
- **Correcting for the wrong family.** The primary endpoint, a pre-specified hierarchy, and exploratory analyses are usually handled differently; applying one blanket correction to all of them can be needlessly conservative.
- **Comparing adjusted p-values with an adjusted α.** Adjusted p-values are compared with the usual 0.05; adjusting both double-counts.
- **Adjusting the p-value but reporting an unadjusted CI.** Where possible, give multiplicity-adjusted CIs (for example Dunnett or Tukey simultaneous intervals).
- **Ignoring multiplicity in subgroup analyses and interim looks.** These are part of the family too (or need alpha-spending approaches).

## Related tests

- [Tukey's HSD](#/test/tukey-hsd): all pairwise comparisons after ANOVA
- [Dunnett's test](#/test/dunnett): several treatments vs one control
- [One-way ANOVA](#/test/one-way-anova): the omnibus test that often precedes post-hoc comparisons
- [Repeated-measures ANOVA](#/test/repeated-measures-anova): pairwise visit comparisons need adjustment
- [Kruskal–Wallis test](#/test/kruskal-wallis): nonparametric omnibus test, followed by adjusted pairwise tests
- [Independent t-test](#/test/independent-t): the source of many of the raw p-values

## Key takeaways

- Each extra test adds a chance of a false positive: 10 true nulls give about a 40% chance of at least one.
- **FWER** (Bonferroni, Holm) controls the chance of *any* false positive; use it for confirmatory claims.
- **FDR** (Benjamini–Hochberg) controls the *share* of false discoveries; use it for exploratory screening.
- Holm is never less powerful than Bonferroni, so prefer it.
- The family (m) must be defined in advance and include every test performed.

## References

- Holm S. A simple sequentially rejective multiple test procedure. *Scand J Stat.* 1979;6:65–70.
- Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc B.* 1995;57:289–300.
- ICH E9 *Statistical Principles for Clinical Trials* (1998); FDA *Multiple Endpoints in Clinical Trials: Guidance for Industry* (2022).
- R documentation: `?p.adjust`; statsmodels: `statsmodels.stats.multitest.multipletests`; SAS/STAT User's Guide: The MULTTEST Procedure
