# Welch's ANOVA

> Compares the means of three or more independent groups **without assuming equal variances**, using a weighted F-test with adjusted degrees of freedom; a safer default than the classic one-way ANOVA.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Compare means of 3+ independent groups when SDs or group sizes differ | Continuous outcome | Independent observations; roughly normal within groups; unequal variances allowed | CRP across four illness-severity groups with SDs from 0.8 to 4.0 mg/L |

## Learning objectives

- Explain why the classic one-way ANOVA can mislead when variances are unequal
- Describe how Welch's ANOVA down-weights noisy groups and uses a fractional denominator df
- Run Welch's ANOVA and Games–Howell in R and Python and read F, df, and adjusted p
- Report an omnibus plus post-hoc result without over-reading a significant F

## What it is

Classic [one-way ANOVA](#/test/one-way-anova) pools the within-group variance into a single number, which is only sensible if every group has the same spread. Welch's ANOVA drops that requirement. It weights each group mean by n / s² (so noisy groups count less) and corrects the denominator degrees of freedom, just as Welch's t-test does for two groups. With two groups it reduces exactly to the Welch t-test.

- **Test statistic:** F_W, a weighted between-group F; it is compared with an F distribution with **df₁ = k − 1** and a **fractional df₂** that is estimated from the group sizes and variances
- **Effect size:** group means with differences and CIs (primary); ω² or η² (secondary, see below)
- **Post-hoc:** [Games–Howell](#/test/multiple-comparisons) pairwise comparisons, which also do not assume equal variances

**Why make it the default?** When group variances are unequal, the classic F-test can have a wrong type I error rate, especially if the **smaller groups have the larger variances**. When variances happen to be equal, Welch's ANOVA loses very little power. Because pre-testing for equal variances (with [Levene's test](#/test/levene)) and then choosing a method distorts the error rates and is unreliable in small samples, several methodological papers recommend simply using the Welch version routinely, for instance Delacre, Lakens and Leys (2017) for two groups and Delacre et al. (2019) for the multi-group case.

## The intuition

ANOVA is still a **signal-to-noise ratio**. The signal is how far apart the group means are. The noise is how much each group mean would wobble by chance, which depends on that group's own spread and size (s²/n). The classic F-test averages all the within-group variances into one pooled number. That is fine if every group is equally noisy, but misleading if one group is far noisier than the rest.

Analogy: four clinics each report an average waiting time. One clinic's times are very consistent; another's range from 5 minutes to 3 hours. You would trust the consistent clinic's average more than the erratic one's. Welch's ANOVA does exactly this: it weights each group mean by n/s², so noisy groups count for less, and it lowers the denominator df to reflect that the noise estimate is itself uncertain (24.1 instead of 67 here).

The plot shows the CRP values. Notice how the Healthy group is tightly packed while Moderate and Severe are spread out, with a single Severe patient at 17.1 mg/L. The Moderate mean is clearly highest, but with only 9 patients and a wide spread, its mean is itself uncertain.

<div class="viz" data-viz="dots" data-unit="CRP, mg/L" data-groups='{"Healthy":[1.4,1.8,0.9,0.9,3.1,1.2,3.3,2.4,1.8,1.1,1.2,1.5],"Mild":[1.1,2.3,1.4,2.7,2.4,4.5,1.9,1.9,1.8,2.7,2.1,3.3,2.8,2.7,2.4,1.7,2.4,1.5,1.4,1.6,4,1.1,1.6,3.6,2.4,6.6,1.9,2.2,1.1,2.7],"Moderate":[8.2,1.7,8.2,2.7,0.5,8.1,7.1,4.5,7.5],"Severe":[1.6,0.9,3.1,2.3,3.2,1.7,2.9,7.8,2.1,0.9,17.1,1.7,0.3,3.2,9.4,0.9,0.7,1.1,3.9,3]}' data-caption="Spread grows with severity. Group means differ, but the wide groups carry more uncertainty."></div>

The Welch statistic is F = 5.22 with 3 and 24.1 degrees of freedom. If all four true means were equal, F values this large would be rare. ANOVA only looks at the upper tail.

<div class="viz" data-viz="f" data-stat="5.22" data-df1="3" data-df2="24.1"></div>

**What to notice:** the observed F sits well into the right tail, and the shaded area (p = 0.006) is small. The F test says "the means are not all equal", but not *which* ones differ. That is the job of the post-hoc step.

## When to use it

Use Welch's ANOVA when **all** of these are true:

- The outcome is **continuous** (CRP, creatinine, tumour volume, pain score)
- There are **three or more independent groups**
- The group SDs are **unequal** or you do not want to assume that they are equal, which is especially likely when group sizes also differ
- Within each group the data are approximately normal (or n is moderately large)

**Clinical scenario.** C-reactive protein (CRP, mg/L), an inflammatory marker, is measured in 71 patients in four groups: healthy controls (n = 12), mild illness (n = 30), moderate illness (n = 9) and severe illness (n = 20). CRP spread rises with disease severity, with SDs from 0.8 to 4.0 mg/L, and group sizes are unequal. *Do mean CRP levels differ across the four groups?*

**Don't use it when:**

- Variances are plainly equal and the design is balanced → the classic [one-way ANOVA](#/test/one-way-anova) is equivalent and simpler to explain
- The data are strongly skewed or ordinal, or have extreme outliers → use the [Kruskal–Wallis test](#/test/kruskal-wallis), or analyse log-CRP
- There are only **two groups** → use the [Welch independent t-test](#/test/independent-t)
- You need to adjust for covariates (baseline values, age) → use [ANCOVA](#/test/ancova)
- The same patients are measured repeatedly → use [repeated-measures ANOVA](#/test/repeated-measures-anova) or a [mixed model](#/test/mixed-models)

## Assumptions

1. **Observations are independent**, within and between groups.
2. **The outcome is roughly normal within each group** (equivalently, the residuals are roughly normal). Check Q–Q plots of the residuals or per group; formal tests like [Shapiro–Wilk](#/test/shapiro-wilk) have little power in small groups. Welch's ANOVA fixes unequal variances, **not** skewness. CRP is typically right-skewed, so a log transformation is often worth considering in real data (this page keeps the raw scale for clarity).
3. **The outcome is on an interval or ratio scale.**
4. **Variances may differ**: this is what makes Welch's ANOVA different from the classic test. You do not need to test for equal variances first, but looking at the SDs by group is good practice.

## Hypotheses

- **H₀:** μ₁ = μ₂ = … = μ_k (all group means are equal)
- **H₁:** at least one group mean differs

## Example data

CRP (mg/L) by illness severity (simulated once in R with `set.seed(11)` from skewed distributions, then rounded to 0.1 mg/L and hard-coded below):

| Group | n | Mean | SD |
|---|---|---|---|
| Healthy | 12 | 1.72 | 0.81 |
| Mild | 30 | 2.39 | 1.15 |
| Moderate | 9 | 5.39 | 3.08 |
| Severe | 20 | 3.39 | 3.95 |

The largest SD is almost five times the smallest, and group sizes range from 9 to 30. That combination of unequal n and unequal variances is where the classic F-test can be unreliable. The raw values are in the code blocks.

## How to run it

### R

```r
# CRP (mg/L) by illness severity: unequal n and unequal SDs
healthy  <- c(1.4, 1.8, 0.9, 0.9, 3.1, 1.2, 3.3, 2.4, 1.8, 1.1, 1.2, 1.5)
mild     <- c(1.1, 2.3, 1.4, 2.7, 2.4, 4.5, 1.9, 1.9, 1.8, 2.7, 2.1, 3.3, 2.8, 2.7,
              2.4, 1.7, 2.4, 1.5, 1.4, 1.6, 4, 1.1, 1.6, 3.6, 2.4, 6.6, 1.9, 2.2, 1.1, 2.7)
moderate <- c(8.2, 1.7, 8.2, 2.7, 0.5, 8.1, 7.1, 4.5, 7.5)
severe   <- c(1.6, 0.9, 3.1, 2.3, 3.2, 1.7, 2.9, 7.8, 2.1, 0.9, 17.1, 1.7, 0.3, 3.2,
              9.4, 0.9, 0.7, 1.1, 3.9, 3)
d <- data.frame(
  crp   = c(healthy, mild, moderate, severe),
  group = factor(rep(c("Healthy", "Mild", "Moderate", "Severe"),
                     c(12, 30, 9, 20)),
                 levels = c("Healthy", "Mild", "Moderate", "Severe")))

# 1. Look at n, mean and SD per group, then test equal variances
aggregate(crp ~ group, d, function(x) c(n = length(x), mean = mean(x), sd = sd(x)))
car::leveneTest(crp ~ group, d, center = median)

# 2. Welch's ANOVA vs the classic equal-variance ANOVA
oneway.test(crp ~ group, data = d, var.equal = FALSE)
oneway.test(crp ~ group, data = d, var.equal = TRUE)

# 3. Effect size from the classic ANOVA table: eta-squared and omega-squared
a  <- anova(lm(crp ~ group, d))
ss <- a$"Sum Sq"; df <- a$Df; ms_w <- ss[2] / df[2]
eta2   <- ss[1] / sum(ss)
omega2 <- (ss[1] - df[1] * ms_w) / (sum(ss) + ms_w)
round(c(eta2 = eta2, omega2 = omega2), 3)

# 4. Games-Howell post-hoc (Tukey-type test with Welch df for each pair)
m <- tapply(d$crp, d$group, mean); v <- tapply(d$crp, d$group, var)
n <- tapply(d$crp, d$group, length); k <- length(m)
gh <- do.call(rbind, lapply(combn(k, 2, simplify = FALSE), function(ij) {
  i <- ij[1]; j <- ij[2]
  se  <- sqrt(v[i] / n[i] + v[j] / n[j])
  dfw <- se^4 / ((v[i] / n[i])^2 / (n[i] - 1) + (v[j] / n[j])^2 / (n[j] - 1))
  q   <- abs(m[i] - m[j]) / se * sqrt(2)          # studentized range statistic
  data.frame(contrast = paste(names(m)[j], "-", names(m)[i]),
             diff = m[j] - m[i], se = se, df = dfw,
             p_adj = ptukey(q, k, dfw, lower.tail = FALSE))
}))
print(gh, digits = 3, row.names = FALSE)
```

**What the code does**

- The four vectors are stacked into one long data frame with a `group` factor; `levels = ` fixes the order of the groups.
- `aggregate(...)` shows n, mean and SD per group; always look at these first.
- `car::leveneTest(..., center = median)` is the Brown–Forsythe version of Levene's test. It is shown here for information; as explained above, the decision to use Welch's ANOVA need not depend on it.
- `oneway.test(..., var.equal = FALSE)` **is** Welch's ANOVA; `FALSE` is the default, so `oneway.test(crp ~ group, d)` alone also gives Welch. With `var.equal = TRUE` it reproduces the classic one-way ANOVA, as in `aov()`.
- The effect size is taken from the classic ANOVA table: η² = SS_between / SS_total, and ω² = (SS_between − df_between × MS_within) / (SS_total + MS_within), which is less biased in small samples.
- Games–Howell is written out here because it is not in base R. For each pair it computes a Welch standard error and Welch–Satterthwaite df, converts the t value to a studentized range value (`q = |diff| / se × √2`), and gets the family-wise adjusted p-value from `ptukey()`. (Packages such as `rstatix::games_howell_test` or `PMCMRplus::gamesHowellTest` do the same.)

**Output**

```text
     group      crp.n   crp.mean     crp.sd
1  Healthy 12.0000000  1.7166667  0.8144528
2     Mild 30.0000000  2.3933333  1.1488925
3 Moderate  9.0000000  5.3888889  3.0836036
4   Severe 20.0000000  3.3900000  3.9519349

Levene's Test for Homogeneity of Variance (center = median)
      Df F value  Pr(>F)  
group  3  3.2367 0.02756 *
      67                  

	One-way analysis of means (not assuming equal variances)

data:  crp and group
F = 5.216, num df = 3.000, denom df = 24.114, p-value = 0.006437


	One-way analysis of means

data:  crp and group
F = 4.5345, num df = 3, denom df = 67, p-value = 0.00592

  eta2 omega2 
 0.169  0.130 
           contrast   diff    se    df  p_adj
     Mild - Healthy  0.677 0.315 28.61 0.1626
 Moderate - Healthy  3.672 1.054  8.84 0.0297
   Severe - Healthy  1.673 0.914 21.60 0.2871
    Moderate - Mild  2.996 1.049  8.68 0.0767
      Severe - Mild  0.997 0.908 21.16 0.6949
  Severe - Moderate -1.999 1.356 19.67 0.4707
```

### Python

```python
import numpy as np, pandas as pd
import pingouin as pg
from scipy import stats

# CRP (mg/L) by illness severity: unequal n and unequal SDs
healthy  = [1.4, 1.8, 0.9, 0.9, 3.1, 1.2, 3.3, 2.4, 1.8, 1.1, 1.2, 1.5]
mild     = [1.1, 2.3, 1.4, 2.7, 2.4, 4.5, 1.9, 1.9, 1.8, 2.7, 2.1, 3.3, 2.8, 2.7,
            2.4, 1.7, 2.4, 1.5, 1.4, 1.6, 4, 1.1, 1.6, 3.6, 2.4, 6.6, 1.9, 2.2, 1.1, 2.7]
moderate = [8.2, 1.7, 8.2, 2.7, 0.5, 8.1, 7.1, 4.5, 7.5]
severe   = [1.6, 0.9, 3.1, 2.3, 3.2, 1.7, 2.9, 7.8, 2.1, 0.9, 17.1, 1.7, 0.3, 3.2,
            9.4, 0.9, 0.7, 1.1, 3.9, 3]
names = ["Healthy", "Mild", "Moderate", "Severe"]
groups = [healthy, mild, moderate, severe]
d = pd.DataFrame({"crp": np.concatenate(groups),
                  "group": np.repeat(names, [len(g) for g in groups])})
d["group"] = pd.Categorical(d["group"], categories=names)

# 1. n, mean, SD per group; Levene (median-centred = Brown-Forsythe)
print(d.groupby("group", observed=True)["crp"].agg(["count", "mean", "std"]).round(3))
print(stats.levene(*groups, center="median"))

# 2. Welch's ANOVA (pingouin) vs classic ANOVA (scipy)
print(pg.welch_anova(data=d, dv="crp", between="group").round(4))
print(stats.f_oneway(*groups))

# 3. Effect size from the classic ANOVA: eta- and omega-squared
grand = d["crp"].mean()
ssb = sum(len(g) * (np.mean(g) - grand) ** 2 for g in groups)
ssw = sum(((np.array(g) - np.mean(g)) ** 2).sum() for g in groups)
dfb, dfw = len(groups) - 1, len(d) - len(groups)
msw = ssw / dfw
print(f"eta2 = {ssb / (ssb + ssw):.3f}, "
      f"omega2 = {(ssb - dfb * msw) / (ssb + ssw + msw):.3f}")

# 4. Games-Howell post-hoc
gh = pg.pairwise_gameshowell(data=d, dv="crp", between="group")
gh["diff"] = gh["mean_B"] - gh["mean_A"]          # B minus A, as in the R code
print(gh[["A", "B", "diff", "se", "df", "pval"]].round(4).to_string(index=False))
```

**What the code does**

- `pg.welch_anova` is Welch's ANOVA from the `pingouin` package; SciPy and statsmodels have no one-line version. It returns F, df₁ (`ddof1`), df₂ (`ddof2`), the p-value (`p_unc`) and `np2`.
- `np2` here is the partial η², which in a one-way design equals η² of the classic ANOVA (0.1688). It is **not** computed from Welch's weights.
- `stats.f_oneway` is the classic equal-variance ANOVA, for comparison with R's `var.equal = TRUE`.
- `stats.levene(..., center="median")` is the Brown–Forsythe-type Levene test; SciPy's default is already `"median"`.
- ω² is computed by hand from the ANOVA sums of squares, using the same formula as R.
- `pg.pairwise_gameshowell` gives Games–Howell. pingouin reports `diff` as A − B, so the code recomputes it as B − A to match R's contrasts. The df and adjusted p-values match R's.

**Output**

```text
          count   mean    std
group                        
Healthy      12  1.717  0.814
Mild         30  2.393  1.149
Moderate      9  5.389  3.084
Severe       20  3.390  3.952
LeveneResult(statistic=np.float64(3.236688282068302), pvalue=np.float64(0.027557698250252716))
  Source  ddof1    ddof2      F   p_unc     np2
0  group      3  24.1142  5.216  0.0064  0.1688
F_onewayResult(statistic=np.float64(4.534508159406568), pvalue=np.float64(0.0059201211247353475))
eta2 = 0.169, omega2 = 0.130
       A        B    diff     se      df   pval
 Healthy     Mild  0.6767 0.3151 28.6057 0.1626
 Healthy Moderate  3.6722 1.0544  8.8414 0.0297
 Healthy   Severe  1.6733 0.9144 21.5982 0.2871
    Mild Moderate  2.9956 1.0491  8.6760 0.0767
    Mild   Severe  0.9967 0.9082 21.1574 0.6949
Moderate   Severe -1.9989 1.3555 19.6715 0.4707
```

## Interpreting the output

- **Group means:** 1.72 (healthy), 2.39 (mild), 5.39 (moderate), 3.39 (severe) mg/L. SDs rise from 0.81 to 3.95 mg/L, so equal variances are not tenable; Levene's test agrees (p = 0.028).
- **Welch's ANOVA: F(3, 24.1) = 5.22, p = 0.0064.** There is evidence that mean CRP is not the same in all four groups. R and Python give identical F, df and p. Note the **fractional denominator df of 24.1**, far below the 67 of the classic test, which is Welch's price for unequal variances.
- **Classic ANOVA: F(3, 67) = 4.53, p = 0.0059.** In this dataset the two tests reach the same conclusion and the p-values are close. That is not guaranteed: with a different pattern of n and variances (for instance, the smallest group with the largest variance) the classic test can be badly wrong. The point of using Welch's ANOVA is that you do not have to hope.
- **Effect size:** η² = 0.17 and ω² = 0.13, i.e., group membership accounts for roughly 13–17% of the variance in CRP. These are computed from the classic sums of squares, and are a descriptive approximation under unequal variances. Prefer to also report the group means and the pairwise differences with CIs.
- **Games–Howell:** only Moderate vs Healthy is significant (difference 3.67 mg/L, adjusted p = 0.030, Welch df 8.8). Moderate vs Mild is borderline (3.00 mg/L, p = 0.077). The other four comparisons are not significant. Severe has the largest SD (3.95), which is why the Severe comparisons have wide uncertainty despite n = 20.
- **Overall vs pairwise:** a significant omnibus test does not guarantee that every pair differs; here it is mainly driven by the moderate group.

**How to report it:**
> CRP differed across the four severity groups (Welch's ANOVA F(3, 24.1) = 5.22, p = 0.006; ω² = 0.13). Games–Howell post-hoc tests showed that CRP was higher in moderate illness than in healthy controls (mean difference 3.7 mg/L, adjusted p = 0.030); no other pairwise comparison was significant at α = 0.05.

**Brown–Forsythe as an alternative.** The Brown–Forsythe *ANOVA* (not to be confused with the Brown–Forsythe *variance* test, i.e., median-centred Levene) is another heteroscedasticity-robust F-test. It adjusts the F numerator and denominator df with a different approximation. It tends to be a bit more conservative for normal data and may hold up better with skewed data or extreme variance patterns, but Welch's ANOVA is more widely implemented and is the better-studied default. If you pre-specify one, report which.

## Check your understanding

<details class="check"><summary>1. Four groups have SDs of 1, 1, 1 and 6, and the smallest group (n = 8) has the SD of 6. Which ANOVA is safer, and why?</summary>

Welch's ANOVA. The classic F-test pools variances, which understates the noise of the small, highly variable group and can inflate the type I error rate. Welch's version weights by n/s² and corrects the df.
</details>

<details class="check"><summary>2. A colleague runs Levene's test (p = 0.20) and, finding it non-significant, switches to the classic ANOVA. What is wrong with this workflow?</summary>

Pre-testing for equal variances and then choosing the method distorts error rates, and Levene's test has low power in small samples, so p = 0.20 does not show the variances are equal. Using Welch's ANOVA routinely, or pre-specifying the method, is safer.
</details>

<details class="check"><summary>3. The omnibus Welch test has p = 0.006. Does that mean all four groups differ from each other?</summary>

No. It means the four means are not all equal. In this example only Moderate vs Healthy is significant after Games–Howell adjustment (p = 0.030). Always follow with pairwise comparisons that suit the variance pattern.
</details>

<details class="check"><summary>4. If the 17.1 mg/L Severe value were removed, what would you expect to happen to the Severe SD and to the evidence for a Severe effect?</summary>

The Severe SD would fall noticeably (it is driven by that point), and the Severe mean would drop. Welch's ANOVA corrects for unequal variances but not for skewness or outliers, which is why a log transformation of CRP is often worth considering.
</details>

## Common mistakes

- **Choosing the method after a Levene test.** A two-stage procedure (test for equal variances, then pick ANOVA or Welch) changes the type I error rate and Levene's test is underpowered when n is small. Use Welch's ANOVA throughout, or pre-specify the choice.
- **Following Welch's ANOVA with Tukey's HSD.** [Tukey's HSD](#/test/tukey-hsd) assumes equal variances; Games–Howell is the matching post-hoc. See [multiple comparisons](#/test/multiple-comparisons).
- **Believing Welch's ANOVA fixes non-normality.** It handles unequal variances, not skewness or outliers. Consider a log transformation (very often appropriate for CRP) or the [Kruskal–Wallis test](#/test/kruskal-wallis).
- **Reporting the wrong df.** The denominator df is fractional (24.1); report it as given.
- **Mixing up the two Brown–Forsythe procedures.** One tests variances, the other compares means.
- **Reporting only the p-value.** Give group means (or medians), SDs, pairwise differences with CIs and an effect size.

## Related tests

- [One-way ANOVA](#/test/one-way-anova): the equal-variance version
- [Levene's test](#/test/levene): checks equal variances
- [Independent t-test (Student & Welch)](#/test/independent-t): the two-group case
- [Tukey's HSD](#/test/tukey-hsd): post-hoc for equal variances
- [Multiple comparisons](#/test/multiple-comparisons): Games–Howell and other post-hoc adjustments
- [Kruskal–Wallis test](#/test/kruskal-wallis): rank-based alternative

## Key takeaways

- Welch's ANOVA compares 3+ means **without assuming equal variances**; it weights groups by n/s².
- Its denominator df is fractional (here 24.1), which is the price of not pooling variances.
- Use it routinely rather than choosing after a Levene test.
- It fixes unequal variances, **not** skewness or outliers.
- Follow a significant omnibus test with Games–Howell, not Tukey's HSD, when variances differ.

## References

- Welch BL. On the comparison of several mean values: an alternative approach. *Biometrika.* 1951;38(3/4):330–336.
- Delacre M, Leys C, Mora YL, Lakens D. Taking parametric assumptions seriously: arguments for the use of Welch's F-test instead of the classical F-test in one-way ANOVA. *International Review of Social Psychology.* 2019;32(1):13.
- Games PA, Howell JF. Pairwise multiple comparison procedures with unequal n's and/or variances. *Journal of Educational Statistics.* 1976;1(2):113–125.
- R documentation: `?oneway.test`; pingouin documentation: `pingouin.welch_anova`, `pingouin.pairwise_gameshowell`
