# Meta-analysis (fixed and random effects)

> Pools the effect estimates from several studies into one overall estimate, weighting each study by its precision, and asks how much the studies disagree beyond chance.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Combining effect estimates (log-OR, mean difference) from several similar studies | Effect estimate + standard error per study | Studies independent; estimates approximately normal; comparable question | Seven trials of a drug on mortality, pooled as an odds ratio |

## Learning objectives

- Explain how inverse-variance weighting turns many trial results into one pooled estimate
- Distinguish the fixed-effect model from the random-effects model and say when each answers the question
- Run both in R and Python and read Q, I², and τ²
- Interpret a forest plot and recognize the limits of meta-analysis (publication bias, heterogeneity)

## What it is

A meta-analysis treats each study as one data point: an **effect estimate** yᵢ (for example a log odds ratio or a mean difference) and its **standard error** SEᵢ. The pooled estimate is a weighted average, with weights wᵢ = 1/SEᵢ², so precise (usually large) studies count more.

- **Fixed-effect model:** assumes every study estimates the **same** true effect θ. Differences between studies are sampling error only. Pooled θ̂ = Σ wᵢyᵢ / Σ wᵢ, with SE = 1/√Σ wᵢ.
- **Random-effects model:** assumes each study has its **own** true effect, drawn from a distribution with mean μ and between-study variance **τ²**. Weights become 1/(SEᵢ² + τ²), which pulls the weights closer together and widens the CI. The result estimates the **average** effect across the population of studies.
- **DerSimonian–Laird (DL)** is the classic moment estimator of τ². R's `metafor::rma()` defaults to **REML**, which is generally preferred today. DL can underestimate τ² (and give too-narrow CIs) when there are few studies, and is not recommended as the only analysis; REML, often with the Hartung–Knapp CI adjustment, is safer. This lesson uses `method = "DL"` in both languages so the numbers match.
- **Heterogeneity statistics:** **Q** = Σ wᵢ(yᵢ − θ̂)² (compare with χ² on k − 1 df); **τ²** (between-study variance, on the scale of yᵢ²); **I²** = (Q − df)/Q, the share of total variability due to real heterogeneity rather than chance.
- **Effect size to report:** pooled estimate with 95% CI (for ratios, back-transform from the log scale), plus τ² and a **95% prediction interval** (where the true effect in a new study is expected to lie).

## The intuition

Think of each trial as a noisy measurement of the treatment effect. A big trial is a precise instrument, a small trial a shaky one. Averaging while trusting the precise instruments more gives a better estimate than any single trial.

The model question is the only real difference between fixed and random effects. If the trials differ in dose, population, or follow-up, there is no single "true effect", only a spread of true effects. The random-effects model admits that spread by adding τ² to every study's variance. When the studies agree closely, τ² is near 0 and the two models give nearly the same answer. When they disagree, random effects gives a wider CI and weights small studies relatively more.

Here are the seven study log-ORs from our example (negative = treatment lowers the odds of death). Notice that most points lie left of zero, but the spread is wider than the individual standard errors would suggest, and two trials sit on the wrong side of zero.

<div class="viz" data-viz="dots" data-unit="log odds ratio" data-groups='{"Study log-ORs":[-0.85,-0.10,-0.70,0.25,-0.40,0.05,-0.60]}' data-zero="1" data-caption="Seven trials, all plotted on the log-OR scale. The line at 0 is no effect. Most trials favor treatment, two do not."></div>

## When to use it

- You have **two or more** studies that ask a similar clinical question with a comparable outcome and comparison
- Each study provides an **effect estimate and a standard error** (or the counts or means and SDs to compute them)
- The studies were found with a **systematic, reproducible search** and not just the ones you happened to know

**Clinical scenario.** Seven randomized trials compare a new anticoagulant with standard care for 90-day mortality. Individually, most are too small to be conclusive. *What is the overall odds ratio, and is it consistent across trials?*

**Reporting standard.** A meta-analysis should follow a protocol and be reported with the **PRISMA** checklist: a flow diagram of studies identified, screened, and included; the search strategy; risk-of-bias assessment; and the statistical methods. Register the protocol (for example in PROSPERO) before looking at results.

**Don't use it when:**

- There are too few or very different studies (clinically unlike populations or outcomes): a narrative synthesis is more honest
- You only have **one** trial: use the single-study test, such as a [two-proportion z-test](#/test/two-proportion-z) or [logistic regression](#/test/logistic-regression)
- You have **patient-level** data from several centers: use a [mixed model](#/test/mixed-models) or stratified analysis such as [Cochran–Mantel–Haenszel](#/test/cmh)
- You want to compare ORs, RRs, and risk differences before choosing the effect measure: see [Effect measures](#/test/effect-measures)

## Assumptions

1. **Studies are independent.** The same patients must not appear in two studies. Check for overlapping publications.
2. **Study estimates are approximately normal with known SEs.** Reasonable for log-ORs and mean differences in moderately sized studies; with sparse events (zero cells) use special methods.
3. **The studies are similar enough to pool.** Judge this clinically first (population, dose, outcome definition), then look at Q, I², and τ².
4. **The set of studies is not selectively missing.** Check a funnel plot (see below) and consider a search of trial registries.
5. **(Random effects)** True effects are roughly normally distributed around a mean. With few studies (k < 10) τ² is poorly estimated, so interpret with care.

## Hypotheses

Meta-analysis is mainly estimation, but two tests are usually reported.

- **Pooled effect:** H₀: θ = 0 (for ratio measures, OR = 1) vs H₁: θ ≠ 0
- **Heterogeneity (Q test):** H₀: all studies share one true effect (τ² = 0) vs H₁: the true effects differ. This test has **low power with few studies**, so a non-significant Q does not prove homogeneity.
- **Estimated:** the pooled log-OR (and OR), τ², and I².

## Example data

Seven hypothetical trials; the effect is the log odds ratio of 90-day death, treatment vs standard care.

| Trial | Aalto 2008 | Brandt 2010 | Chen 2012 | Diaz 2013 | Evans 2015 | Fujita 2017 | Gomez 2019 |
|---|---|---|---|---|---|---|---|
| log OR | −0.85 | −0.10 | −0.70 | 0.25 | −0.40 | 0.05 | −0.60 |
| SE | 0.30 | 0.15 | 0.35 | 0.25 | 0.12 | 0.28 | 0.20 |

## How to run it

### R

```r
library(metafor)

# Seven trials: log odds ratio (treatment vs control, outcome = death) and its SE
dat <- data.frame(
  trial = c("Aalto 2008", "Brandt 2010", "Chen 2012", "Diaz 2013",
            "Evans 2015", "Fujita 2017", "Gomez 2019"),
  yi  = c(-0.85, -0.10, -0.70, 0.25, -0.40, 0.05, -0.60),
  sei = c( 0.30,  0.15,  0.35, 0.25,  0.12, 0.28,  0.20))

# Fixed effect (inverse variance) and random effects (DerSimonian-Laird).
# metafor's default is REML; we ask for DL so Python matches.
fe <- rma(yi, sei = sei, data = dat, method = "FE")
dl <- rma(yi, sei = sei, data = dat, method = "DL")

show <- function(lab, m) cat(sprintf(
  "%-14s logOR = %.3f (%.3f, %.3f)   OR = %.3f (%.3f, %.3f)\n", lab,
  coef(m), m$ci.lb, m$ci.ub, exp(coef(m)), exp(m$ci.lb), exp(m$ci.ub)))
show("fixed effect", fe); show("random effect", dl)

# Heterogeneity: Q, p, tau^2, I^2
cat(sprintf("Q = %.3f (df = %d), p = %.3f\n", dl$QE, dl$k - 1, dl$QEp))
cat(sprintf("tau2 = %.4f, tau = %.3f, I2 = %.1f%%\n", dl$tau2, sqrt(dl$tau2), dl$I2))

# 95% prediction interval for the true OR in a new trial (t with k - 2 df)
pi <- coef(dl) + c(-1, 1) * qt(0.975, dl$k - 2) * sqrt(dl$tau2 + dl$se^2)
cat(sprintf("Prediction interval: OR %.3f to %.3f\n", exp(pi[1]), exp(pi[2])))

# Forest plot (saved to a file; atransf = exp shows the axis as odds ratios)
png("forest.png", width = 700, height = 400)
forest(dl, slab = dat$trial, atransf = exp, xlab = "Odds ratio (log scale)")
invisible(dev.off())
```

**What the code does**

- `rma(yi, sei = sei, method = "FE")` is the fixed-effect model; `method = "DL"` is DerSimonian–Laird random effects. Leave the method out and you get REML, which will differ slightly.
- `dl$QE`, `dl$QEp`, `dl$tau2`, and `dl$I2` hold Q, its p-value, τ², and I².
- `forest()` draws the forest plot (one row per study, a diamond for the pooled estimate). `atransf = exp` relabels the log axis as odds ratios.

**Output**

```text
fixed effect   logOR = -0.309 (-0.452, -0.166)   OR = 0.734 (0.636, 0.847)
random effect  logOR = -0.313 (-0.569, -0.057)   OR = 0.731 (0.566, 0.945)
Q = 15.777 (df = 6), p = 0.015
tau2 = 0.0679, tau = 0.261, I2 = 62.0%
Prediction interval: OR 0.346 to 1.547
```

### Python

```python
import numpy as np
import matplotlib
matplotlib.use("Agg")
from scipy.stats import chi2
from statsmodels.stats.meta_analysis import combine_effects

trial = ["Aalto 2008", "Brandt 2010", "Chen 2012", "Diaz 2013",
         "Evans 2015", "Fujita 2017", "Gomez 2019"]
yi  = np.array([-0.85, -0.10, -0.70, 0.25, -0.40, 0.05, -0.60])  # log OR
sei = np.array([ 0.30,  0.15,  0.35, 0.25,  0.12, 0.28,  0.20])

# combine_effects takes VARIANCES (sei**2). method_re="dl" = DerSimonian-Laird.
res = combine_effects(yi, sei**2, method_re="dl", row_names=trial)
sf = res.summary_frame()
for row in ["fixed effect", "random effect"]:
    e, lo, hi = sf.loc[row, ["eff", "ci_low", "ci_upp"]]
    print(f"{row:14s} logOR = {e:.3f} ({lo:.3f}, {hi:.3f})   "
          f"OR = {np.exp(e):.3f} ({np.exp(lo):.3f}, {np.exp(hi):.3f})")

# Heterogeneity: Q, p, tau^2, I^2
k = len(yi)
print(f"Q = {res.q:.3f} (df = {k-1}), p = {chi2.sf(res.q, k-1):.3f}")
print(f"tau2 = {res.tau2:.4f}, tau = {np.sqrt(res.tau2):.3f}, "
      f"I2 = {max(0, (res.q-(k-1))/res.q)*100:.1f}%")

# 95% prediction interval for the true OR in a new trial (t with k - 2 df)
from scipy.stats import t
mu, se = sf.loc["random effect", ["eff", "sd_eff"]]
half = t.ppf(0.975, k - 2) * np.sqrt(res.tau2 + se**2)
print(f"Prediction interval: OR {np.exp(mu - half):.3f} to {np.exp(mu + half):.3f}")

# Forest plot: studies, fixed and random pooled estimates (OR scale)
fig = res.plot_forest(use_exp=True)
fig.savefig("forest.png", dpi=100)
```

**What the code does**

- `combine_effects` expects **variances** (`sei**2`), not standard errors, a common slip.
- `method_re="dl"` selects DerSimonian–Laird. The summary frame also has rows for the scale-adjusted ("wls") variants, which we ignore here.
- The function reports `q` and `tau2`; the p-value and I² are computed by hand from Q.
- `plot_forest(use_exp=True)` draws the forest plot on the odds-ratio scale.

**Output**

```text
fixed effect   logOR = -0.309 (-0.452, -0.166)   OR = 0.734 (0.636, 0.847)
random effect  logOR = -0.313 (-0.569, -0.057)   OR = 0.731 (0.566, 0.945)
Q = 15.777 (df = 6), p = 0.015
tau2 = 0.0679, tau = 0.261, I2 = 62.0%
Prediction interval: OR 0.346 to 1.547
```

R and Python agree to the third decimal because both use DL.

## Interpreting the output

- **Fixed effect:** OR = 0.734 (95% CI 0.636 to 0.847). If all trials estimate one effect, treatment lowers the odds of death by about 27%.
- **Random effects (DL):** OR = 0.731 (95% CI 0.566 to 0.945). The point estimate is almost identical, but the CI is much **wider** because it now includes between-trial variation. This is the more cautious summary here.
- **Q = 15.8 on 6 df, p = 0.015.** Evidence that the trials differ by more than sampling error.
- **τ² = 0.068 (τ = 0.26).** The true log-OR varies between trials with an SD of about 0.26. That is an odds ratio ranging roughly from 0.56 to 0.95 for trials one SD either side of the mean (e^(−0.313 ± 0.26)). This is a spread of true trial effects, not the CI of the pooled mean; it only happens to look similar to that CI here.
- **Prediction interval: OR 0.35 to 1.55.** This is the range in which the true effect of a *new* trial is expected to fall (it combines τ² with the uncertainty in the pooled mean, using t with k − 2 df). It is much wider than the CI and includes 1, so although the average effect is beneficial, a future setting could show no benefit or harm.
- **I² = 62%.** About 62% of the variability in the estimates reflects real differences between trials rather than chance: moderate to substantial heterogeneity. I² describes **proportion**, not the size of the disagreement. τ² (or a prediction interval) tells you how large it is.
- **Forest plot.** Each square is a trial (size proportional to weight), each line its 95% CI, and the diamond the pooled estimate. Look for CIs that do not overlap the diamond. Here Diaz 2013 and Fujita 2017 sit on the other side of 1.

**Funnel plots and publication bias.** A funnel plot shows each study's effect against its precision. Without bias, points form a symmetric funnel around the pooled value. If small, unimpressive studies are missing, the plot is lopsided. In R, `funnel(dl)` and `regtest(dl)` (Egger-type test) do this. Use them only with about 10 or more studies, and remember that asymmetry can also come from real heterogeneity.

**How to report it:**
> Seven randomized trials were pooled with a DerSimonian–Laird random-effects model. Treatment was associated with lower 90-day mortality (OR 0.73, 95% CI 0.57 to 0.95), with moderate heterogeneity (Q = 15.8, df = 6, p = 0.015; I² = 62%; τ² = 0.068; 95% prediction interval 0.35 to 1.55). The fixed-effect estimate was similar (OR 0.73, 95% CI 0.64 to 0.85).

## Check your understanding

<details class="check"><summary>1. A meta-analysis reports I² = 0%. A reviewer writes "the studies agree, so we can treat them as identical." What is wrong?</summary>

I² = 0% means the observed variation is **no more than** expected by chance. It does not show the studies are identical, especially with few small studies where Q has little power. The CIs of the individual studies may be wide enough to hide clinically important differences. Judge similarity from the populations, interventions, and outcomes too, and look at τ² or a prediction interval.
</details>

<details class="check"><summary>2. Which model, fixed or random, gives a wider CI when the trials are heterogeneous, and why?</summary>

Random effects. It adds τ² to each study's variance, so the weights are less extreme and the pooled SE is larger. Here the CI widened from 0.636–0.847 to 0.566–0.945. When τ² = 0 the two models coincide.
</details>

<details class="check"><summary>3. A large new trial with a very small SE and a result opposite to the others is added. How do the fixed-effect and random-effects pooled estimates respond?</summary>

The fixed-effect estimate is dominated by the big trial, because its weight 1/SE² is huge. In the random-effects model, the new heterogeneity raises τ², which flattens the weights, so the big trial gets **relatively less** influence and the CI widens. Neither result is "right"; the right response is to ask why the trials disagree.
</details>

<details class="check"><summary>4. Someone says "the pooled result is statistically significant, so the evidence is solid." Name two reasons that may be false.</summary>

First, the pooled result is only as good as the included trials: bias in each trial (poor randomization, selective reporting) carries through. Second, **publication bias** or selective searching can make the included set unrepresentative. Add to that heterogeneity: a significant average effect may hide trials with no benefit. A significant pooled p-value says nothing about these.
</details>

## Common mistakes

- **Pooling raw counts across trials** (summing events and patients) instead of pooling effect estimates. This can reverse the direction of the effect (Simpson's paradox).
- **Choosing the model by whether Q is significant.** Decide on the question and on clinical grounds, and report τ² and a prediction interval.
- **Passing SEs where variances are expected** (or the reverse). Check which the function wants.
- **Mixing scales.** Pool log-ORs, then back-transform; do not average ORs directly.
- **Reading I² as the size of the heterogeneity.** It is a proportion that depends on study precision.
- **Searching only for the studies you know.** Report a systematic search, follow PRISMA, and assess publication bias.

## Related tests

- [Effect measures](#/test/effect-measures): which effect measure to pool (OR, RR, risk difference)
- [Cochran–Mantel–Haenszel test](#/test/cmh): stratified pooling of 2×2 tables, a fixed-effect cousin
- [Cochran's Q test](#/test/cochrans-q): a different test with a similar name, for repeated binary outcomes
- [Confidence intervals](#/test/confidence-intervals): how to read the pooled CI
- [Linear mixed-effects models](#/test/mixed-models): random effects for patient-level clustered data

## Key takeaways

- Meta-analysis pools **effect estimates**, weighting each by 1/variance, never raw counts.
- Fixed effect assumes one true effect; random effects assumes a distribution of effects and gives a wider, more cautious CI when studies differ.
- Q tests heterogeneity, τ² measures its size, and I² gives the share of variability that is not chance. I² = 0 does not mean the studies agree.
- Always look at the forest plot, and check a funnel plot for publication bias when there are enough studies.
- Report the method (DL vs REML), the search, and PRISMA compliance, not just the pooled number.

## References

- DerSimonian R, Laird N. Meta-analysis in clinical trials. *Controlled Clinical Trials.* 1986;7:177–188.
- Higgins JPT, Thompson SG. Quantifying heterogeneity in a meta-analysis. *Statistics in Medicine.* 2002;21:1539–1558.
- Page MJ et al. The PRISMA 2020 statement. *BMJ.* 2021;372:n71.
- Viechtbauer W. Conducting meta-analyses in R with the metafor package. *J Stat Softw.* 2010;36(3).
