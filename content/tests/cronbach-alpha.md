# Cronbach's alpha

> Measures the **internal consistency** of a multi-item scale, i.e. how strongly the items covary, so that a total score can be treated as a reasonably reliable measurement of one underlying trait such as quality of life.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Check that items of a questionnaire hang together before summing them | Ordinal or continuous item scores (e.g., 1–5 Likert) | Items measure one construct; similar scaling; uncorrelated errors; all keyed in the same direction | Six-item quality-of-life scale in 120 patients |

## Learning objectives

- Explain what internal consistency is and what Cronbach's alpha does and does not tell you
- Run `psych::alpha` in R and `pingouin.cronbach_alpha` in Python, including alpha-if-item-dropped
- Interpret alpha, item-rest correlations and the average inter-item correlation
- Describe the key caveats (alpha is not unidimensionality; it rises with item count) and when McDonald's omega is a better choice

## What it is

Cronbach's alpha (α) summarizes how much the items of a scale vary **together** compared with how much they vary **individually**. If the items all measure the same underlying trait, people who score high on one item will tend to score high on the others, the items' variances will add up to far less than the variance of the total score, and α will be high.

α = [k / (k − 1)] · [1 − Σ σ²ᵢ / σ²_total]

where k is the number of items, σ²ᵢ is the variance of item i and σ²_total is the variance of the summed score. The **standardized alpha** is based on the average inter-item correlation r̄: α_std = k·r̄ / (1 + (k − 1)·r̄).

- **Statistic:** α, from 0 to 1 in practice (it can be negative if items are inversely related, usually because of un-reversed items)
- **Benchmarks:** about 0.70 is the usual minimum for group-level research, 0.80 is good, 0.90 or above is needed for individual decisions. Values above about 0.95 often mean redundant items.
- **What to report:** α with a **95% CI**, the average inter-item correlation, and item-level diagnostics (item-rest correlation, alpha-if-item-dropped).

## The intuition

Imagine six thermometers placed in the same room. If they all track the room temperature, their readings move together, and the average of the six is a stable number. If one thermometer ignores the room, it adds only noise: the readings no longer move together so tightly. Alpha is a measure of that "moving together" for the items of a questionnaire.

The dot plot shows total scores (6 to 30) for the first 40 of 120 patients. The scores spread widely (the SD is 5.2 across all patients) because the items agree about who is doing better or worse. If the items were unrelated, each patient's high and low item scores would cancel out and the totals would bunch near the middle.

<div class="viz" data-viz="dots" data-unit="total score (6 to 30)" data-groups='{"Total QoL score":[17,14,22,17,15,19,6,17,28,15,23,21,19,11,17,20,17,14,16,15,15,20,16,19,21,24,13,18,24,19,10,16,26,13,20,16,23,25,21,23]}' data-caption="Total quality-of-life score (sum of six 1-5 items) for the first 40 of 120 patients. The wide spread comes from items that agree about who is better or worse off."></div>

Alpha also grows with the **number of items**, even if the average inter-item correlation stays the same. Our six items have average r̄ = 0.411 and α = 0.807 (standardized). The Spearman–Brown formula shows what the same r̄ would give for other scale lengths:

| Number of items | α (with the same r̄ = 0.411) |
|---|---|
| 2 | 0.583 |
| 3 | 0.677 |
| 4 | 0.736 |
| 6 | 0.807 |
| 12 | 0.893 |
| 20 | 0.933 |

So a long scale can reach a "good" alpha with only mediocre inter-item correlations. Alpha is a statement about **length plus correlation**, not about item quality alone.

## When to use it

Use Cronbach's alpha when:

- You have a **multi-item scale or index** meant to measure **one construct** and plan to sum or average the items
- Items are **ordinal or continuous** and all scored in the **same direction** (reverse-code negatively worded items first)
- You are developing, validating or **reporting the reliability** of a patient-reported outcome measure in your own sample
- You have at least about 5 items and 100+ respondents (smaller samples give wide CIs)

**Clinical scenario.** A rehabilitation clinic develops a six-item quality-of-life screening scale (physical function, pain, energy, mood, social participation and sleep), each scored 1 (worst) to 5 (best). It is given to 120 patients. *Do the items hang together well enough to report a single total score?*

**Don't use it when:**

- You want agreement between **raters** or repeated measurements of the same quantity → use the [intraclass correlation (ICC)](#/test/icc) or [Cohen's kappa](#/test/cohens-kappa)
- You want agreement between **two measurement methods** → use [Bland–Altman analysis](#/test/bland-altman)
- The scale has **several dimensions** (subscales) → compute alpha or omega for each subscale after checking the structure (factor analysis)
- Items are **formative** (they cause the construct, like a list of symptoms or comorbidities) rather than reflective → internal consistency is not meaningful
- You want **criterion validity** or discrimination → use [ROC analysis](#/test/roc-auc) or correlations with an external criterion

## Assumptions

1. **Items measure a single underlying construct** (unidimensionality). Alpha does **not** test this; check with a factor analysis or scree plot.
2. **Items are scored in the same direction.** Reverse-code first; otherwise alpha is artificially low or negative.
3. **Tau-equivalence:** items have about equal loadings on the construct (equal "importance"). When loadings differ a lot, alpha underestimates reliability; McDonald's omega does not need this.
4. **Errors are uncorrelated** across items. Near-duplicate wording or items sharing a stem inflate alpha.
5. **Reasonably continuous items.** 5-point Likert items work acceptably; with very few categories or skewed responses, polychoric-based alpha or ordinal omega can be better.
6. **Adequate sample size** (at least about 100 respondents, ideally 5–10 per item) and no floor or ceiling effects.

## Hypotheses

Alpha is an **estimation** method, not a hypothesis test. It estimates the reliability coefficient of the total score. A formal test is rarely used; the useful questions are "is α at least about 0.70?" and "how wide is the CI?", the latter via Feldt's CI or a bootstrap.

## Example data

Simulated responses from 120 patients on six 1–5 items: physical function, pain, energy, mood, social participation, sleep. All items are scored so that higher means better. Item means range 2.8 to 3.0 and SDs 1.2 to 1.3.

| Item | Mean | SD |
|---|---|---|
| Physical function | 2.92 | 1.18 |
| Pain | 3.02 | 1.22 |
| Energy | 2.78 | 1.18 |
| Mood | 2.89 | 1.21 |
| Social participation | 2.93 | 1.28 |
| Sleep | 3.03 | 1.28 |

The 720 responses are typed inline in each code block below (one string of digits per item, one digit per patient).

## How to run it

### R

```r
library(psych)
# 120 patients answer six quality-of-life items scored 1 (worst) to 5 (best)
bits <- function(s) as.integer(strsplit(gsub("\\s", "", s), "")[[1]])
items <- data.frame(
  physical = bits("33432313533322242231244245355414423245432532121321544351342453215341
    3323223352342413324422323333242453213351124135332322"),
  pain     = bits("23543313522341424334232453225422534354342333422132343442251453113251
    2424232251542324445224234453212453322351323325312514"),
  energy   = bits("43322314534432232223334343333412534343334552122312533541154331112432
    2232434241342323245114324133122451213341121323431414"),
  mood     = bits("42443412534322422331421333333213322245341543142333344334143353311354
    1245222521442523234325423121322552212451332235342521"),
  social   = bits("22412313425442244223152225133223414344542544132121444341152452313452
    2434223341452213325322145254141354224322343435353113"),
  sleep    = bits("21233312425442353233333535125332523324352542334531334131232541115422
    2235424452523324535221345241142355433431134445233511"))

# 1. Cronbach's alpha, with alpha-if-item-dropped and item-rest correlations
a <- alpha(items)
round(a$total[, c("raw_alpha", "std.alpha", "average_r")], 3)
c(lower = a$feldt$lower.ci[[1]], upper = a$feldt$upper.ci[[1]])   # Feldt 95% CI
round(cbind(a$alpha.drop[, "raw_alpha", drop = FALSE],
            r.drop = a$item.stats$r.drop), 3)               # one row per item

# 2. McDonald's omega (total) from a one-factor model
ld <- fa(items, nfactors = 1)$loadings[, 1]
round(sum(ld)^2 / (sum(ld)^2 + sum(1 - ld^2)), 3)
```

**What the code does**

- `psych::alpha()` computes raw alpha (from covariances), standardized alpha (from correlations), the average inter-item correlation, and item-level statistics. It warns about reverse-keyed items; here all items are scored in the same direction.
- `a$total` has the scale-level summaries; `a$feldt` has Feldt's 95% CI for alpha.
- `a$alpha.drop` shows **alpha if each item is dropped**, and `r.drop` is the **corrected item-rest correlation** (each item against the sum of the others).
- `fa(items, nfactors = 1)` fits a one-factor model; **McDonald's omega** (total) is then (Σλ)² / [(Σλ)² + Σ(1 − λ²)], using the loadings λ.

**Output**

```text
 raw_alpha std.alpha average_r
     0.806     0.807     0.411
    lower     upper 
0.7473142 0.8554912 
         raw_alpha r.drop
physical     0.754  0.666
pain         0.772  0.585
energy       0.773  0.579
mood         0.799  0.461
social       0.766  0.609
sleep        0.792  0.497
[1] 0.81
```

### Python

```python
import numpy as np, pandas as pd
import pingouin as pg

# 120 patients answer six quality-of-life items scored 1 (worst) to 5 (best)
bits = lambda s: np.array(list(s.replace(" ", "")), dtype=int)
items = pd.DataFrame({
    "physical": bits("334323135333222422312442453554144232454325321213215443513424532153 "
    "413323223352342413324422323333242453213351124135332322"),
    "pain":     bits("235433135223414243342324532254225343543423334221323434422514531132 "
    "512424232251542324445224234453212453322351323325312514"),
    "energy":   bits("433223145344322322233343433334125343433345521223125335411543311124 "
    "322232434241342323245114324133122451213341121323431414"),
    "mood":     bits("424434125343224223314213333332133222453415431423333443341433533113 "
    "541245222521442523234325423121322552212451332235342521"),
    "social":   bits("224123134254422442231522251332234143445425441321214443411524523134 "
    "522434223341452213325322145254141354224322343435353113"),
    "sleep":    bits("212333124254423532333335351253325233243525423345313341312325411154 "
    "222235424452523324535221345241142355433431134445233511")})

# 1. Cronbach's alpha with its 95% CI
alpha, ci = pg.cronbach_alpha(data=items)
print(f"alpha = {alpha:.3f}, 95% CI ({ci[0]:.3f}, {ci[1]:.3f})")

# 2. Alpha if each item is dropped, and corrected item-rest correlations
for c in items.columns:
    rest = items.drop(columns=c)
    print(f"{c:<9} alpha if dropped = {pg.cronbach_alpha(data=rest)[0]:.3f}, "
          f"item-rest r = {np.corrcoef(items[c], rest.sum(axis=1))[0, 1]:.3f}")

# 3. Average inter-item correlation and standardized alpha
R = items.corr().values
rbar = R[np.triu_indices(6, 1)].mean()
print(f"mean r = {rbar:.3f}, standardized alpha = {6 * rbar / (1 + 5 * rbar):.3f}")
```

**What the code does**

- `pg.cronbach_alpha(data=items)` returns alpha and a 95% CI (the Feldt interval, same as R). It expects one **column per item**, one row per respondent.
- The loop recomputes alpha without each item (alpha if item dropped) and each item's item-rest correlation.
- The last lines compute the average inter-item correlation and standardized alpha by hand from the correlation matrix.
- `pingouin` has no omega; compute it from a one-factor model (R `psych::omega`, or the `factor_analyzer` package in Python).

**Output**

```text
alpha = 0.806, 95% CI (0.747, 0.855)
physical  alpha if dropped = 0.754, item-rest r = 0.666
pain      alpha if dropped = 0.772, item-rest r = 0.585
energy    alpha if dropped = 0.773, item-rest r = 0.579
mood      alpha if dropped = 0.799, item-rest r = 0.461
social    alpha if dropped = 0.766, item-rest r = 0.609
sleep     alpha if dropped = 0.792, item-rest r = 0.497
mean r = 0.411, standardized alpha = 0.807
```

## Interpreting the output

- **α = 0.806 (95% CI 0.747 to 0.855).** Good internal consistency by the usual benchmarks (≥ 0.70 acceptable, ≥ 0.80 good). The CI shows that the true value could plausibly be anywhere from acceptable to very good. R and Python agree (0.806; 0.747, 0.855).
- **Standardized α = 0.807; average inter-item correlation r̄ = 0.411.** Raw and standardized alpha are almost identical because the items have similar SDs. An average inter-item correlation of 0.2–0.5 is typical for a coherent scale; much higher suggests redundancy.
- **Alpha if item dropped:** every value is below 0.806 (range 0.754 to 0.799), so removing any item would **lower** alpha. No item is hurting the scale. Physical function is the most valuable item (dropping it gives 0.754); mood (0.799) and sleep (0.792) are the least.
- **Corrected item-rest correlations:** from 0.461 (mood) to 0.666 (physical function). All are comfortably above the common 0.30 rule of thumb; mood is the weakest link, though still acceptable.
- **McDonald's omega (R) = 0.81**, essentially equal to alpha here because the loadings are fairly similar (tau-equivalence roughly holds). When loadings differ more, omega will be higher than alpha.
- **What alpha does not say:** it does not show that the six items form one dimension, that the scale is valid, or that a patient's score is precise. It does not say the scale measures quality of life, only that its items are consistent with each other.

**How to report it:**
> The six-item quality-of-life scale showed good internal consistency in this sample (n = 120): Cronbach's α = 0.81 (95% CI 0.75–0.86; McDonald's ω = 0.81), with a mean inter-item correlation of 0.41 and corrected item-total correlations of 0.46–0.67. Deleting any single item reduced α.

## Caveats and alternatives

- **Alpha is not unidimensionality.** A scale made of two distinct subscales (for example, physical items and psychological items) can still yield a high alpha, especially with many items. Check the structure with an exploratory factor analysis or a scree plot, and compute alpha separately for each subscale.
- **Alpha rises with the number of items.** Doubling a scale's length raises alpha even if each added item is no better, as the table above showed. Compare alphas only across scales of similar length, and report r̄ too.
- **Very high alpha (> 0.95) can signal redundancy** (near-identical items) rather than quality.
- **Alpha is not a fixed property of a scale;** it describes a scale **in a sample**. Report it each time, with a CI.
- **Alpha is a lower bound** of reliability that is exact only under tau-equivalence. **McDonald's omega (ω)** uses the factor loadings and is generally preferred when loadings differ; it is easy to compute (`psych::omega`) and should be reported alongside or instead of alpha. For ordinal items with few categories or skewed responses, use polychoric correlations (ordinal alpha/omega).
- **Alpha is not test–retest reliability or inter-rater agreement.** For stability over time or between raters, use the [ICC](#/test/icc) or [kappa](#/test/cohens-kappa).

## Check your understanding

<details class="check"><summary>1. A 20-item symptom checklist has α = 0.92. The authors conclude it measures a single construct. Is that justified?</summary>

No. A high alpha is compatible with several dimensions, especially with 20 items; alpha rises with the number of items. A factor analysis (scree plot, loadings) is needed to support unidimensionality, and omega or alpha per subscale would then be reported.
</details>

<details class="check"><summary>2. In our output, dropping the "mood" item gives α = 0.799, compared with 0.806 for all six. Should we remove mood to improve the scale?</summary>

No. Dropping it gives a **lower** alpha (0.799 < 0.806), so mood is helping, even if only slightly. An item is a candidate for removal only if alpha **increases** when it is dropped (or its item-rest correlation is very low). Content validity matters too: mood may be essential to the construct even if it correlates less with the others.
</details>

<details class="check"><summary>3. A researcher adds 6 more items with the same average inter-item correlation (r̄ = 0.411) to make a 12-item scale. What happens to alpha, and is the scale "better"?</summary>

Standardized alpha rises from 0.807 to 12 × 0.411 / (1 + 11 × 0.411) = 0.893, with no improvement in how strongly the items relate to each other. The scale is longer and more reliable in the sense of a lower error in the total score, but the patient burden is also doubled. Alpha alone cannot tell whether the extra items are worthwhile.
</details>

<details class="check"><summary>4. A scale has an alpha of −0.2 in your data. What is the most likely reason?</summary>

One or more items are scored in the opposite direction (for example, a negatively worded item that has not been reverse-coded), so the item covariances are negative and cancel the positive ones. Reverse-code those items and recompute. Another possible reason is that the items measure unrelated things.
</details>

## Common mistakes

- **Reading alpha as proof of one dimension.** It is not; check the factor structure.
- **Forgetting to reverse-code** negatively worded items, which gives a falsely low or negative alpha.
- **Chasing a high alpha by adding near-duplicate items.** It inflates alpha and adds burden but no information.
- **Applying alpha to formative indices** (symptom counts, comorbidity lists), where items need not correlate.
- **Reporting alpha without a CI** or without the sample size, especially with fewer than 100 respondents.
- **Using alpha for raters or repeated measurements.** Use the ICC or kappa.
- **Dropping items solely on statistical grounds,** ignoring content validity.
- **Quoting "α ≥ 0.70" as a universal pass mark.** The needed level depends on the use, group-level research versus individual decisions.

## Related tests

- [Intraclass correlation (ICC)](#/test/icc): reliability for raters, repeated measures and agreement
- [Cohen's kappa](#/test/cohens-kappa): agreement for categorical ratings
- [Bland–Altman analysis](#/test/bland-altman): agreement between two measurement methods
- [Pearson correlation](#/test/pearson): the inter-item correlations behind alpha
- [Spearman correlation](#/test/spearman): rank-based item correlations for ordinal items

## Key takeaways

- Alpha measures **internal consistency**: how strongly the items of a scale covary. It is estimated in your sample and needs a CI.
- About 0.70 is acceptable for group research, 0.80 good, above 0.95 often redundant.
- Alpha is **not** evidence of unidimensionality, and it **rises with the number of items** even when item quality is unchanged.
- Use **alpha-if-item-dropped** and item-rest correlations to spot weak items, but weigh content validity too.
- **McDonald's omega** relaxes the equal-loadings assumption and is a better modern default; report it alongside alpha.

## References

- Cronbach LJ. Coefficient alpha and the internal structure of tests. *Psychometrika.* 1951;16(3):297–334.
- Tavakol M, Dennick R. Making sense of Cronbach's alpha. *Int J Med Educ.* 2011;2:53–55.
- McNeish D. Thanks coefficient alpha, we'll take it from here. *Psychol Methods.* 2018;23(3):412–433.
- R documentation: `?psych::alpha`; pingouin documentation: `pingouin.cronbach_alpha`
