# Kendall's tau

> A rank correlation that measures association between two ordinal or non-normal variables by comparing every pair of subjects: how many pairs are ordered the same way (concordant) versus opposite ways (discordant).

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Monotonic association between two ordinal or non-normal variables, especially with many ties or small n | Ordinal or continuous, paired measurements | Independent subjects; monotonic relationship | Disease severity grade vs 6-minute walk distance |

## Learning objectives

- Define concordant and discordant pairs and compute tau by hand from them
- Explain how tau-b corrects for ties
- Run `cor.test(method = "kendall")` and `scipy.stats.kendalltau`, and read tau and p
- Contrast Kendall's tau with Spearman's rho and choose between them

## What it is

Take any two patients, i and j. If the patient with the higher value of X also has the higher value of Y, the pair is **concordant**; if the order is reversed, it is **discordant**; if either variable is tied, the pair is a tie. Kendall's tau summarises the balance.

- **Tau-a** = (C − D) / [n(n − 1)/2], with C concordant and D discordant pairs. It is only right when there are no ties.
- **Tau-b** = (C − D) / √[(n₀ − n₁)(n₀ − n₂)], where n₀ = n(n − 1)/2 and n₁, n₂ are the numbers of pairs tied on X and on Y. It adjusts the denominator for ties and is the default in R and SciPy.
- **Range:** −1 (all pairs discordant) to +1 (all concordant); 0 means no monotonic association.
- **Effect size:** tau itself. With no ties, (1 + τ)/2 is the probability that a random pair is concordant. With ties, tau-b is slightly smaller in size than the share (C − D)/(C + D) computed on untied pairs, so quote the pair counts as well.

## The intuition

Tau answers a simple question: if I pick two patients at random, how much more likely is it that they are ordered the *same* way on both measures than the *opposite* way (ignoring pairs that tie)? That makes it easy to say aloud. Spearman's rho, in contrast, is a Pearson correlation on ranks, which is harder to translate into a statement about patients.

In our data, 15 patients have a severity grade (1 to 5) and a walk distance. Higher grades go with shorter walks. The dot plot shows walk distance by grade group: the mean falls at each step.

<div class="viz" data-viz="dots" data-unit="6MWD, m" data-groups='{"Grade 1-2":[520,480,500,450,470],"Grade 3":[410,430,390,440],"Grade 4-5":[350,380,400,300,340,310]}' data-caption="Walk distance falls as severity rises. Almost every pair of patients is discordant (higher grade, shorter walk)."></div>

Standardised, tau-b gives z = −4.32. The picture below shows where that falls under the null (no association).

<div class="viz" data-viz="z" data-stat="-4.32"></div>

The observed z is far out in the left tail, so the two-sided p-value area is minuscule.

## When to use it

- Two variables, each **ordinal or continuous**, measured on the same subjects
- You expect a **monotonic** (not necessarily linear) relationship
- The data have **many ties** (e.g., a 5-point grade) or a **small sample**, where tau's p-values behave better than Spearman's
- You want an interpretation in terms of **pairs**

**Clinical scenario.** In 15 patients with chronic obstructive pulmonary disease, disease severity is graded 1 to 5 (an ordinal scale with many ties) and the 6-minute walk distance (6MWD, metres) is measured. *Is greater severity associated with a shorter walk distance, and how strongly?*

**Don't use it when:**

- Both variables are continuous, roughly bivariate normal, and the relationship is linear → use [Pearson correlation](#/test/pearson)
- You prefer a rank correlation with an R²-like reading and have little tie → [Spearman correlation](#/test/spearman) works
- You need to adjust for other variables or predict one from the other → use [linear regression](#/test/linear-regression) or [ordinal regression](#/test/ordinal-regression)
- You want to test agreement, not association → use [Cohen's kappa](#/test/cohens-kappa), the [ICC](#/test/icc), or [Bland–Altman](#/test/bland-altman)

## Assumptions

1. **Independent subjects.** One pair of values per subject; no repeated measures per subject.
2. **Both variables at least ordinal.** Larger must mean more, in both.
3. **Monotonic relationship.** Check with a scatter plot. A U-shape can give tau near 0 despite a strong relationship.
4. **Ties are handled.** Use tau-b (or tau-c for rectangular tables); the plain tau-a is wrong with ties.

## Hypotheses

- **H₀:** τ = 0 (no monotonic association; concordant and discordant pairs equally likely)
- **H₁:** τ ≠ 0 (two-sided)

## Example data

| Patient | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Severity grade | 1 | 1 | 2 | 2 | 2 | 3 | 3 | 3 | 3 | 4 | 4 | 4 | 5 | 5 | 5 |
| 6MWD (m) | 520 | 480 | 500 | 450 | 470 | 410 | 430 | 390 | 440 | 350 | 380 | 400 | 300 | 340 | 310 |

## How to run it

### R

```r
# Disease severity grade (ordinal, 1-5, many ties) vs 6-minute walk distance (m), 15 patients
grade <- c(1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 5)
walk  <- c(520, 480, 500, 450, 470, 410, 430, 390, 440, 350, 380, 400, 300, 340, 310)

# Kendall's tau-b (handles ties); exact = FALSE because of ties
cor.test(grade, walk, method = "kendall", exact = FALSE)

# Spearman for comparison
cor.test(grade, walk, method = "spearman", exact = FALSE)

# Count concordant and discordant pairs by hand
n <- length(grade); C <- D <- 0
for (i in 1:(n-1)) for (j in (i+1):n) {
  s <- sign(grade[i] - grade[j]) * sign(walk[i] - walk[j])
  if (s > 0) C <- C + 1 else if (s < 0) D <- D + 1
}
c(concordant = C, discordant = D, total_pairs = n*(n-1)/2)
```

**What the code does**

- `cor.test(..., method = "kendall")` reports tau-b and a z statistic. R switches to the normal approximation (with a tie-corrected variance) whenever ties are present; `exact = FALSE` makes this explicit.
- The Spearman call is for comparison only.
- The loop counts pairs. The product of the two signs is +1 for a concordant pair, −1 for a discordant pair, and 0 for a tie on either variable.

**Output**

```text
	Kendall's rank correlation tau

data:  grade and walk
z = -4.3171, p-value = 1.581e-05
alternative hypothesis: true tau is not equal to 0
sample estimates:
       tau 
-0.8792842 

	Spearman's rank correlation rho

data:  grade and walk
S = 1096.1, p-value = 2.199e-08
alternative hypothesis: true rho is not equal to 0
sample estimates:
       rho 
-0.9573663 

 concordant  discordant total_pairs 
          2          87         105 
```

### Python

```python
import numpy as np
from scipy import stats

# Disease severity grade (ordinal, 1-5, many ties) vs 6-minute walk distance (m), 15 patients
grade = np.array([1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 5])
walk  = np.array([520, 480, 500, 450, 470, 410, 430, 390, 440, 350, 380, 400, 300, 340, 310])

# Kendall's tau-b (SciPy default variant is "b"); p-value uses a normal approximation with ties
k = stats.kendalltau(grade, walk, variant="b", method="asymptotic")
print(f"Kendall tau-b = {k.statistic:.4f}, p = {k.pvalue:.3e}")

# Spearman for comparison
s = stats.spearmanr(grade, walk)
print(f"Spearman rho = {s.statistic:.4f}, p = {s.pvalue:.3e}")

# Count concordant and discordant pairs by hand
C = D = 0
for i in range(len(grade) - 1):
    for j in range(i + 1, len(grade)):
        sgn = np.sign(grade[i] - grade[j]) * np.sign(walk[i] - walk[j])
        C += sgn > 0; D += sgn < 0
n = len(grade)
print(f"concordant = {C}, discordant = {D}, total pairs = {n*(n-1)//2}")
```

**What the code does**

- `stats.kendalltau` returns tau-b by default (`variant="c"` gives tau-c). `method="asymptotic"` uses the normal approximation, matching R when there are ties; `method="exact"` is only valid without ties.
- `stats.spearmanr` gives rho for comparison.
- The loop counts concordant and discordant pairs; the remaining 16 pairs are tied on grade.

**Output**

```text
Kendall tau-b = -0.8793, p = 1.581e-05
Spearman rho = -0.9574, p = 2.199e-08
concordant = 2, discordant = 87, total pairs = 105
```

## Interpreting the output

- **Pairs.** Of 105 possible pairs of patients, 16 are tied on severity grade and are not informative, 2 are concordant, and 87 are discordant.
- **Tau-b = −0.88.** (2 − 87) / √[(105 − 16)(105)] = −0.879. A very strong negative association: among pairs with different grades, about 98% (87/89) are discordant. In words: if you pick two patients with different severity, the sicker one almost always walks less far.
- **z = −4.32, p < 0.001.** Strong evidence against no monotonic association.
- **Spearman's rho = −0.96.** Same direction and same conclusion, but numerically larger. Rho is usually larger than tau for the same data (about 1.5 times for moderate associations; closer to 1 for a very strong one like this), so never compare them directly or apply the same "strong/weak" cut-offs.
- **No CI printed.** Tau itself should be reported with a CI where possible (bootstrap, or the Fisher-type interval from add-on packages).

**How to report it:**
> Disease severity grade was strongly and inversely associated with 6-minute walk distance (Kendall's τ_b = −0.88, z = −4.32, p < 0.001, n = 15); among patient pairs with different severity grades, 98% were discordant.

⚠️ Association is not causation, and a strong tau does not mean agreement, nor does it say how much walk distance changes per grade; that is a regression question.

## Check your understanding

<details class="check"><summary>1. Which would you use to relate a 5-grade pain scale to a 4-grade satisfaction rating in 25 patients: Pearson or Kendall?</summary>

Kendall's tau (or Spearman). Both variables are ordinal with heavy ties, so means and linearity are not meaningful. Use tau-b.
</details>

<details class="check"><summary>2. A paper reports Spearman ρ = 0.60 and a second paper reports Kendall τ = 0.45 for similar variables. Is the first association stronger?</summary>

Not necessarily. For the same data, |ρ| is usually larger than |τ| (roughly 1.5 times for moderate values), so the two numbers are on different scales. Compare tau with tau and rho with rho.
</details>

<details class="check"><summary>3. Tau-b = 0.05, p = 0.70. A researcher concludes the two variables are unrelated. What else could be going on?</summary>

Tau only detects **monotonic** association. A U-shaped or inverted-U relationship would give tau near 0 despite a strong dependency. Always plot the data. The small sample may also be underpowered.
</details>

<details class="check"><summary>4. Calculation: what would tau-a (no tie correction) be for our data, and why is it smaller in magnitude than tau-b?</summary>

Tau-a = (2 − 87) / 105 = −0.81. Tau-b divides by √[(105 − 16)(105)] ≈ 96.7 instead of 105, because 16 pairs are tied on grade and cannot be ordered. Tau-a counts those tied pairs as "no association", which understates the relationship, so |τ_a| < |τ_b| (0.81 vs 0.88).
</details>

## Common mistakes

- **Using tau-a (no tie correction) with tied data.** Use tau-b; R and SciPy do this by default.
- **Using `exact = TRUE` with ties.** An exact p-value is only valid for untied data.
- **Comparing tau with rho or Pearson r.** They are on different scales.
- **Ignoring non-monotonic patterns.** Plot the data.
- **Reporting only the p-value.** Give tau, n, and the direction.
- **Treating tau as a slope or as agreement.** It measures rank association only.

## Related tests

- [Spearman rank correlation](#/test/spearman): the other rank correlation
- [Pearson correlation](#/test/pearson): linear association for continuous data
- [Linear regression](#/test/linear-regression): the size of the relationship, with adjustment
- [Tests for trend](#/test/trend-tests): ordered groups and a binary or continuous outcome
- [Cohen's kappa](#/test/cohens-kappa): agreement rather than association

## Key takeaways

- Tau = (concordant − discordant) / pairs, ranging from −1 to +1.
- Use tau-b when there are ties; it is the default in R and SciPy.
- Tau is often easier to interpret than rho (concordant vs discordant pairs) and behaves better with small n and many ties.
- Tau and rho are on different scales; |ρ| is usually larger.
- Tau detects monotonic, not any, association: plot the data.

## References

- Kendall MG. A new measure of rank correlation. *Biometrika.* 1938;30:81–93.
- Newson R. Parameters behind "nonparametric" statistics: Kendall's tau, Somers' D and median differences. *Stata J.* 2002;2:45–64.
- R documentation: `?cor.test`
- SciPy documentation: `scipy.stats.kendalltau`
