# Cochran's Q test

> Tests whether the proportion of "yes" responses differs across three or more conditions measured on the **same subjects**; the binary-outcome analogue of repeated-measures ANOVA and an extension of McNemar's test.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Same subjects, binary outcome under 3+ conditions | Binary (0/1) repeated measures | Subjects independent; conditions applied to all; enough discordant subjects | Pain relief (yes/no) with drugs A, B, C in a crossover |

## Learning objectives

- Recognize a repeated binary outcome across 3+ conditions
- Explain how Cochran's Q uses subjects as their own controls
- Run it in R and Python, then follow up with adjusted pairwise McNemar tests
- Report the result and avoid over-reading a small-sample follow-up

## What it is

Cochran's Q generalises [McNemar's test](#/test/mcnemar) from 2 to k conditions. Each subject has a row of k 0/1 values. Subjects whose responses are the same in every condition (all 0 or all 1) carry no information about differences between conditions; the test is driven entirely by the subjects whose responses *differ*.

- **Test statistic:** Q = (k − 1)·[k·ΣTⱼ² − N²] / [k·N − ΣLᵢ²], where Tⱼ is the number of "yes" in condition j, Lᵢ the number of "yes" for subject i, and N the total number of "yes".
- **Distribution:** approximately χ² with k − 1 degrees of freedom under H₀.
- **Effect size:** the response proportion in each condition, and the pairwise differences in proportions (with CIs).

## The intuition

Signal: how different the column totals (response rates per condition) are. Noise: how much that would wobble by chance. Because the same people appear in every column, their stable tendencies (some patients always respond, some never do) cancel out. Only people who respond to one drug but not another help the comparison.

In our example, 12 patients each tried three analgesics. Drug A relieved pain in 9 (75%), B in 5 (42%), and C in 3 (25%). Cochran's Q is 9.33 on 2 df.

<div class="viz" data-viz="chisq" data-stat="9.33" data-df="2"></div>

Notice that Q sits well into the right tail of the χ²(2) curve, so the shaded area (the p-value) is small, about 0.009. Three patients (000) never respond and three (111) always respond; both kinds give the same result under every drug, so they contribute nothing. That is how the test automatically ignores uninformative subjects. Only the six patients with mixed rows drive Q.

## When to use it

- The outcome is **binary** and measured under **three or more** conditions or time points
- **Every subject** is measured under **every** condition (crossover, repeated visits, three raters on the same cases)
- Subjects are independent of each other

**Clinical scenario.** In a three-period crossover, 12 patients with migraine each take drug A, B, and C in randomised order (with washout). For each, "pain relief at 2 hours" is recorded as yes/no. *Do the three drugs differ in the proportion of patients who get relief?*

**Don't use it when:**

- There are only **two** conditions → use [McNemar's test](#/test/mcnemar)
- Different patients receive each condition → use the [chi-square test](#/test/chi-square-independence)
- The outcome is continuous or ordinal → use [repeated-measures ANOVA](#/test/repeated-measures-anova) or the [Friedman test](#/test/friedman)
- You have missing values for some conditions, covariates, or carry-over effects to model → use [GEE](#/test/gee) or [mixed models](#/test/mixed-models) for a binary outcome

## Assumptions

1. **Subjects are independent.** One row per subject, no repeated subjects across rows.
2. **Every subject is measured under every condition.** Complete data only; a subject missing a condition is dropped.
3. **Binary outcome with the same coding in every condition.**
4. **Condition order is randomised or balanced** (so period or carry-over effects do not mimic a drug effect).
5. **Enough informative subjects.** The χ² approximation is poor with few informative subjects (a common rule of thumb: at least 4 subjects with mixed responses and N × k ≥ 24); use an exact or permutation version for tiny samples. Our example has 6 mixed subjects and 12 × 3 = 36, so it is acceptable but small.

## Hypotheses

- **H₀:** the probability of "yes" is the same in all k conditions
- **H₁:** at least one condition has a different probability of "yes"

## Example data

| Patient | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Drug A | 1 | 1 | 0 | 1 | 1 | 0 | 1 | 1 | 1 | 0 | 1 | 1 |
| Drug B | 1 | 0 | 0 | 1 | 0 | 0 | 1 | 0 | 1 | 0 | 1 | 0 |
| Drug C | 0 | 0 | 0 | 1 | 0 | 0 | 1 | 0 | 0 | 0 | 1 | 0 |

(1 = pain relief at 2 h)

## How to run it

### R

```r
library(DescTools)
# Pain relief at 2 h (1 = yes, 0 = no): the same 12 patients tried drugs A, B, C
d <- data.frame(
  A = c(1,1,0,1,1,0,1,1,1,0,1,1),
  B = c(1,0,0,1,0,0,1,0,1,0,1,0),
  C = c(0,0,0,1,0,0,1,0,0,0,1,0))
y <- as.matrix(d)
colMeans(y)                      # response rate per drug

# Cochran's Q
CochranQTest(y)

# By formula: Q = (k-1) * (k*sum(T_j^2) - N^2) / (k*N - sum(L_i^2))
k <- ncol(y); Tj <- colSums(y); Li <- rowSums(y); N <- sum(y)
Q <- (k - 1) * (k * sum(Tj^2) - N^2) / (k * N - sum(Li^2))
c(Q = Q, df = k - 1, p = pchisq(Q, k - 1, lower.tail = FALSE))

# Pairwise McNemar tests (exact), Holm-adjusted
prs <- combn(names(d), 2, simplify = FALSE)
res <- t(sapply(prs, function(p) {
  b <- sum(d[[p[1]]] == 1 & d[[p[2]]] == 0)   # A yes, B no
  c <- sum(d[[p[1]]] == 0 & d[[p[2]]] == 1)
  c(discordant_1 = b, discordant_2 = c,
    p_raw = if (b + c > 0) binom.test(b, b + c, 0.5)$p.value else 1)
}))
rownames(res) <- sapply(prs, paste, collapse = " vs ")
cbind(res, p_holm = p.adjust(res[, "p_raw"], "holm"))
```

**What the code does**

- `CochranQTest(y)` takes a subjects × conditions matrix of 0/1. Rows are subjects.
- The formula block reproduces the same Q from column totals (`Tj`) and row totals (`Li`).
- The pairwise step counts discordant pairs per pair of drugs. The exact McNemar test is a binomial test on those counts, with H₀ of 50:50.
- `p.adjust(..., "holm")` controls the family-wise error rate across the 3 comparisons.

**Output**

```text
        A         B         C 
0.7500000 0.4166667 0.2500000 

	Cochran's Q test

data:  y
Q = 9.3333, df = 2, p-value = 0.009404

          Q          df           p 
9.333333333 2.000000000 0.009403563 
       discordant_1 discordant_2   p_raw  p_holm
A vs B            4            0 0.12500 0.25000
A vs C            6            0 0.03125 0.09375
B vs C            2            0 0.50000 0.50000
```

### Python

```python
import numpy as np
from itertools import combinations
from scipy.stats import binomtest
from statsmodels.stats.contingency_tables import cochrans_q
from statsmodels.stats.multitest import multipletests

# Pain relief at 2 h (1 = yes, 0 = no): the same 12 patients tried drugs A, B, C
y = np.array([[1,1,0],[1,0,0],[0,0,0],[1,1,1],[1,0,0],[0,0,0],
              [1,1,1],[1,0,0],[1,1,0],[0,0,0],[1,1,1],[1,0,0]])
print("Response rates:", y.mean(axis=0).round(3))

# Cochran's Q
q = cochrans_q(y)
print(f"Cochran Q = {q.statistic:.4f}, df = {q.df}, p = {q.pvalue:.6f}")

# Pairwise exact McNemar (binomial on discordant pairs), Holm-adjusted
names = ["A", "B", "C"]
rows, raw = [], []
for i, j in combinations(range(3), 2):
    b = int(np.sum((y[:, i] == 1) & (y[:, j] == 0)))
    c = int(np.sum((y[:, i] == 0) & (y[:, j] == 1)))
    p = binomtest(b, b + c, 0.5).pvalue if b + c else 1.0
    rows.append((f"{names[i]} vs {names[j]}", b, c)); raw.append(p)
holm = multipletests(raw, method="holm")[1]
for (lab, b, c), p, ph in zip(rows, raw, holm):
    print(f"{lab}: {b} vs {c} discordant, raw p = {p:.4f}, Holm p = {ph:.4f}")
```

**What the code does**

- `cochrans_q(y)` takes the same subjects × conditions array and returns the statistic, df, and p-value.
- Each pair is tested with an exact binomial test on the discordant counts (`binomtest`), which is the exact McNemar test.
- `multipletests(..., method="holm")` returns adjusted p-values as the second element.

**Output**

```text
Response rates: [0.75  0.417 0.25 ]
Cochran Q = 9.3333, df = 2, p = 0.009404
A vs B: 4 vs 0 discordant, raw p = 0.1250, Holm p = 0.2500
A vs C: 6 vs 0 discordant, raw p = 0.0312, Holm p = 0.0938
B vs C: 2 vs 0 discordant, raw p = 0.5000, Holm p = 0.5000
```

## Interpreting the output

- **Response rates: A 75%, B 42%, C 25%.** Drug A looks best and C worst.
- **Q = 9.33, df = 2, p = 0.0094.** If the three drugs truly had the same response probability, a Q this large would be seen in about 1 study in 100. Both software packages agree.
- **Pairwise McNemar (Holm-adjusted).** A vs C: 6 patients responded to A but not C, none the other way (raw p = 0.031, Holm p = 0.094). A vs B: 4 vs 0 (Holm p = 0.25). B vs C: 2 vs 0 (Holm p = 0.50). **No single pair is significant after adjustment**, even though the omnibus test is.
- **Why the mismatch?** With only 12 patients, each pairwise test has very few discordant pairs and little power. The omnibus test pools all the information. Report all three differences with CIs rather than only stating which pairs "are significant."

**How to report it:**
> Pain relief at 2 hours occurred in 75% (A), 42% (B), and 25% (C) of 12 patients. Response rates differed across drugs (Cochran's Q = 9.33, df = 2, p = 0.009). In Holm-adjusted exact McNemar comparisons, the largest difference was A vs C (6 vs 0 discordant patients, adjusted p = 0.094); no pairwise difference reached significance in this small sample.

## Check your understanding

<details class="check"><summary>1. The same 40 patients are tested for a rash after each of 4 sunscreens (rash yes/no). Which test, and is it paired?</summary>

Cochran's Q. The outcome is binary, there are four conditions, and every patient is measured under all four, so the data are paired (repeated).
</details>

<details class="check"><summary>2. Patients 3, 6, and 10 never respond to any drug. A reviewer says to remove them as "non-responders." What happens to Q if you do?</summary>

Nothing. Subjects with the same response under every condition contribute zero to the numerator and the denominator adjustment of Q, so Q is unchanged. (Removing them would, however, distort the response rates you report.)
</details>

<details class="check"><summary>3. Q is significant (p = 0.009). Can you say "A is better than B and C"?</summary>

No. Q says only that *at least one* condition differs. You need adjusted pairwise comparisons, and here none was significant after Holm. Describe the pattern, report effect sizes with CIs, and acknowledge the low power.
</details>

<details class="check"><summary>4. Calculation: with 3 pairwise tests, what is the Bonferroni-adjusted p for A vs C (raw p = 0.031)? How does Holm compare?</summary>

Bonferroni: 3 × 0.03125 = 0.094. Holm gives the same here (0.094) because A vs C has the smallest raw p, so it is multiplied by 3. Holm is never less powerful than Bonferroni and is slightly better for the other comparisons.
</details>

## Common mistakes

- **Using a chi-square test on repeated data.** It treats the 36 observations as independent and overstates the evidence.
- **Using Q with incomplete rows.** Subjects with a missing condition must be dropped or handled with GEE or mixed models.
- **Skipping the multiplicity adjustment** in follow-up pairwise tests.
- **Interpreting a non-significant pairwise test as "no difference"** in a small sample.
- **Ignoring period and carry-over effects** in crossover studies.

## Related tests

- [McNemar's test](#/test/mcnemar): the two-condition case
- [Friedman test](#/test/friedman): the same design with an ordinal or continuous outcome
- [Multiple-comparison corrections](#/test/multiple-comparisons): Holm, Bonferroni, BH
- [GEE](#/test/gee) and [mixed models](#/test/mixed-models): repeated binary outcomes with missing data or covariates
- [Chi-square test of independence](#/test/chi-square-independence): independent groups

## Key takeaways

- Cochran's Q compares proportions across 3+ conditions in the **same subjects**; df = k − 1.
- Only subjects with **mixed** responses contribute information.
- A significant Q must be followed by adjusted pairwise McNemar tests (Holm or Bonferroni).
- A significant omnibus test with no significant pair is common in small samples; report all effect sizes.
- Complete data are required; otherwise use GEE or mixed models.

## References

- Cochran WG. The comparison of percentages in matched samples. *Biometrika.* 1950;37:256–66.
- Agresti A. *Categorical Data Analysis.* 3rd ed. Wiley; 2013. Ch. 11.
- statsmodels documentation: `cochrans_q`
- DescTools documentation: `CochranQTest`
