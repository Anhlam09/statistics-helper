# Two-way (factorial) ANOVA

> Tests the effects of two categorical factors, and whether they interact, on a continuous outcome, such as whether a drug's effect differs between men and women.

| When to use | Data type | Key assumptions | Example |
|---|---|---|---|
| Effects of two crossed categorical factors, and their interaction | Continuous outcome | Independent observations; roughly normal residuals; similar variances | LDL-C reduction by statin dose and sex |

## Learning objectives

- Distinguish **main effects** from an **interaction** in a factorial design
- Explain why the interaction is tested and interpreted first
- Run a Type III two-way ANOVA in R and Python and read the F tests
- Describe effects within levels of the other factor when an interaction is present

## What it is

A two-way ANOVA splits the variation in an outcome into four parts: the **main effect of factor A** (e.g., drug dose), the **main effect of factor B** (e.g., sex), the **A × B interaction**, and residual error. The interaction asks the question clinicians care about most: *does the effect of one factor depend on the level of the other?* If the answer is yes, the main effects cannot be interpreted on their own.

- **Test statistic:** one F = MS_effect / MS_error for each of A, B, and A × B
- **Degrees of freedom:** (levels − 1) for each main effect, (a − 1)(b − 1) for the interaction, and N − (number of cells) for error
- **Effect size:** cell means and contrasts with 95% CI (primary); partial η² for each term (secondary)
- **Type I vs Type III sums of squares:** with **equal** cell sizes (balanced data) all types agree. With **unequal** cell sizes, which is the norm in clinical data, they do not:
  - **Type I (sequential):** each term is tested after the terms listed before it, so the result depends on the **order** of the terms. This is what R's base `anova()` gives.
  - **Type III (marginal):** each term is tested after all the others, including the interaction, and does not depend on order. This is the SAS default in PROC GLM and what most regulatory submissions use. It requires **sum-to-zero contrasts** in R and Python.

## The intuition

Every F statistic is a **signal-to-noise ratio**: the variation explained by a term, divided by the leftover variation among patients within cells. Two-way ANOVA computes three of them: one for dose, one for sex, and one for the interaction.

An interaction is the idea that "it depends." Imagine two sedatives that both work in young adults, but one loses effect in the elderly. The effect of the drug *depends on* age. Graphically, interaction means the lines connecting cell means are **not parallel**.

Our data show this. In women, the LDL-C reduction rises modestly from placebo to high dose. In men, it climbs much more steeply. Compare the two plots below, drawn on the same LDL-C reduction scale: in men the high-dose group sits far above the others, while in women the three groups overlap more.

<div class="viz" data-viz="dots" data-unit="mg/dL, LDL-C reduction (women)" data-groups='{"Placebo":[6,17,4,10,19,7,8,7],"Low":[15,17,22,12,11],"High":[21,17,21,19]}' data-caption="Women: the dose effect is real but moderate (means 9.75, 15.40, 19.50)."></div>

<div class="viz" data-viz="dots" data-unit="mg/dL, LDL-C reduction (men)" data-groups='{"Placebo":[-5,7,5,11],"Low":[17,19,16,16,11,19],"High":[31,21,20,26,30,35,30,31]}' data-caption="Men: the dose effect is much steeper (means 4.50, 16.33, 28.00). Non-parallel slopes are the interaction."></div>

The interaction test gives F(2, 29) = 5.51. Under H₀ (parallel lines) F follows the curve below; ours falls in the upper tail.

<div class="viz" data-viz="f" data-stat="5.51" data-df1="2" data-df2="29"></div>

Because of the interaction, a single "dose effect" number would be misleading: it would be an average of a modest effect in women and a large effect in men.

## When to use it

Use a two-way ANOVA when **all** of these are true:

- The outcome is **continuous**
- There are **two categorical factors** (crossed, so every combination can occur), each with two or more levels
- Observations are **independent** (different patients in every cell)

**Clinical scenario.** A dose-finding study randomizes 35 adults with hypercholesterolaemia to placebo, low-dose, or high-dose statin. The outcome is the reduction in LDL-C (mg/dL) at 12 weeks. Sex is recorded, and recruitment left the cells unbalanced (4 to 8 patients each). *Does the dose effect differ between women and men?*

**Don't use it when:**

- There is only **one** factor → use [one-way ANOVA](#/test/one-way-anova)
- The same patients are measured repeatedly → use [repeated-measures ANOVA](#/test/repeated-measures-anova) or a [mixed model](#/test/mixed-models)
- A **continuous** covariate such as baseline LDL-C must be adjusted for → use [ANCOVA](#/test/ancova)
- The residuals are clearly non-normal with small n → consider a transformation, or the [Kruskal–Wallis test](#/test/kruskal-wallis) for single-factor questions
- The outcome is binary or a count → use [logistic regression](#/test/logistic-regression) or [Poisson regression](#/test/poisson-regression)

## Assumptions

1. **Independent observations.** Each patient appears in exactly one cell. Check from the design.
2. **Residuals are approximately normal** within cells. Check a Q–Q plot of the model residuals or run [Shapiro–Wilk](#/test/shapiro-wilk) on them.
3. **Equal variances across cells (homoscedasticity).** Check with [Levene's test](#/test/levene) or by comparing cell SDs; a ratio of the largest to smallest SD above about 2 is a warning, especially with unequal cell sizes.
4. **Correct model specification.** If you expect an interaction, include it. Dropping a real interaction biases the main effects.
5. **Adequate cell sizes.** Empty cells make Type III tests of some terms impossible.

## Hypotheses

Three separate tests:

- **Drug:** H₀: the mean LDL-C reduction is the same for all three doses (averaged over sex); H₁: at least one differs
- **Sex:** H₀: the mean reduction is the same for women and men (averaged over dose); H₁: it differs
- **Drug × sex:** H₀: the dose effect is the same in both sexes (parallel lines); H₁: the dose effect differs by sex

Test the interaction **first**. If it is significant, describe the effect of each factor within levels of the other.

## Example data

LDL-C reduction (mg/dL), 35 patients; cell sizes are unequal (n in brackets):

| | Women | Men |
|---|---|---|
| Placebo | 6, 17, 4, 10, 19, 7, 8, 7 (8) | −5, 7, 5, 11 (4) |
| Low dose | 15, 17, 22, 12, 11 (5) | 17, 19, 16, 16, 11, 19 (6) |
| High dose | 21, 17, 21, 19 (4) | 31, 21, 20, 26, 30, 35, 30, 31 (8) |

Cell means: placebo 9.75 (F) and 4.50 (M); low dose 15.40 and 16.33; high dose 19.50 and 28.00.

## How to run it

### R

```r
library(car)
# LDL-C reduction (mg/dL): 3 drug doses x 2 sexes, UNBALANCED cell sizes
cells <- list(
  Placebo.F = c(6,17,4,10,19,7,8,7),   Placebo.M = c(-5,7,5,11),
  Low.F     = c(15,17,22,12,11),       Low.M     = c(17,19,16,16,11,19),
  High.F    = c(21,17,21,19),          High.M    = c(31,21,20,26,30,35,30,31))
d <- do.call(rbind, lapply(names(cells), function(n) {
  p <- strsplit(n, "\\.")[[1]]
  data.frame(drug = p[1], sex = p[2], y = cells[[n]])
}))
d$drug <- factor(d$drug, levels = c("Placebo", "Low", "High"))
d$sex  <- factor(d$sex)
table(d$drug, d$sex)

# Type I (sequential) SS: what base-R anova() gives; depends on term order
anova(lm(y ~ drug * sex, data = d))
anova(lm(y ~ sex * drug, data = d))   # reversed order: Type I SS change

# Type III SS (SAS default): needs sum-to-zero contrasts to be meaningful
options(contrasts = c("contr.sum", "contr.poly"))
fit <- lm(y ~ drug * sex, data = d)
Anova(fit, type = 3)

# Cell means for the interaction
with(d, tapply(y, list(drug, sex), mean))
```

**What the code does**

- The list-to-data-frame step just builds a long data set (one row per patient) from the cell values; `factor(..., levels = ...)` puts the doses in a logical order instead of alphabetical.
- `anova(lm(...))` computes **Type I** (sequential) sums of squares. Compare the two orders: the drug and sex rows change, the interaction row does not (it is always last or adjusted for both).
- `options(contrasts = c("contr.sum", "contr.poly"))` switches from R's default treatment (dummy) coding to **sum-to-zero** coding. This is **essential** for correct Type III tests; with the default coding, `car::Anova(type = 3)` gives misleading main-effect tests when an interaction is present. The option must be set **before** fitting the model.
- `car::Anova(fit, type = 3)` gives Type III SS, matching SAS PROC GLM and Python's `typ=3`. (`type = 2` is also available; it tests main effects without the interaction.)
- `tapply` shows the cell means that the interaction describes.

**Output**

```text
         
          F M
  Placebo 8 4
  Low     5 6
  High    4 8

Response: y     (drug first, then sex: Type I)
          Df  Sum Sq Mean Sq F value    Pr(>F)    
drug       2 1771.60  885.80 38.5688 6.756e-09 ***
sex        1   15.60   15.60  0.6791  0.416636    
drug:sex   2  252.95  126.47  5.5068  0.009392 ** 
Residuals 29  666.03   22.97                      

Response: y     (sex first, then drug: Type I)
          Df  Sum Sq Mean Sq F value   Pr(>F)    
sex        1  234.86  234.86 10.2263 0.003337 ** 
drug       2 1552.33  776.16 33.7952 2.65e-08 ***
sex:drug   2  252.95  126.47  5.5068 0.009392 ** 
Residuals 29  666.03   22.97                     

Anova Table (Type III tests)

Response: y
            Sum Sq Df  F value    Pr(>F)    
(Intercept) 7826.1  1 340.7586 < 2.2e-16 ***
drug        1475.4  2  32.1207  4.42e-08 ***
sex           15.7  1   0.6824  0.415515    
drug:sex     252.9  2   5.5068  0.009392 ** 
Residuals    666.0 29                       

           F        M
Placebo 9.75  4.50000
Low     15.40 16.33333
High    19.50 28.00000
```

### Python

```python
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm

# LDL-C reduction (mg/dL): 3 drug doses x 2 sexes, UNBALANCED cell sizes
cells = {
    ("Placebo","F"): [6,17,4,10,19,7,8,7],  ("Placebo","M"): [-5,7,5,11],
    ("Low","F"):     [15,17,22,12,11],      ("Low","M"):     [17,19,16,16,11,19],
    ("High","F"):    [21,17,21,19],         ("High","M"):    [31,21,20,26,30,35,30,31]}
d = pd.DataFrame([(dr, sx, y) for (dr, sx), ys in cells.items() for y in ys],
                 columns=["drug", "sex", "y"])
print(pd.crosstab(d.drug, d.sex).loc[["Placebo","Low","High"]])

# Type I (sequential) SS with default treatment coding
lv = ['Placebo','Low','High']
m1 = smf.ols("y ~ C(drug, levels=lv) * C(sex)", data=d).fit()
print(anova_lm(m1, typ=1).to_string(float_format=lambda x: f"{x:.4g}"))

# Type III SS: sum-to-zero (Sum) contrasts are essential, otherwise results differ
m3 = smf.ols("y ~ C(drug, Sum) * C(sex, Sum)", data=d).fit()
print(anova_lm(m3, typ=3).to_string(float_format=lambda x: f"{x:.4g}"))

print(d.pivot_table(index="drug", columns="sex", values="y", aggfunc="mean")
        .loc[["Placebo","Low","High"]].round(2))
```

**What the code does**

- The data frame is built from the same cell values as in R. `pd.crosstab` confirms the unbalanced design.
- `anova_lm(model, typ=1)` gives **Type I** (sequential) SS in the order of the terms in the formula, equivalent to R's `anova()`. Writing `C(sex) * C(drug)` instead would reproduce R's second table.
- `C(drug, Sum)` requests **sum-to-zero** (effect) coding. With the default treatment coding, `typ=3` tests for main effects are **not** the same as SAS or `car::Anova` when there is an interaction, which is a very common source of "Python gives a different answer" confusion.
- `anova_lm(m3, typ=3)` gives Type III SS that match R and SAS. The `Intercept` row can be ignored.
- The `levels=lv` argument sets the reference ordering of drug; it does not affect the Type III results.

**Output**

```text
sex      F  M
drug         
Placebo  8  4
Low      5  6
High     4  8
                           df  sum_sq  mean_sq      F    PR(>F)
C(drug, levels=lv)          2    1772    885.8  38.57 6.756e-09
C(sex)                      1    15.6    15.6 0.6791    0.4166
C(drug, levels=lv):C(sex)   2    252.9    126.5  5.507  0.009392
Residual                   29     666    22.97    NaN       NaN
                          sum_sq  df      F    PR(>F)
Intercept                   7826   1  340.8 1.424e-17
C(drug, Sum)                1475   2  32.12  4.42e-08
C(sex, Sum)                15.67   1 0.6824    0.4155
C(drug, Sum):C(sex, Sum)   252.9   2  5.507  0.009392
Residual                     666  29    NaN       NaN
sex          F      M
drug                 
Placebo   9.75   4.50
Low      15.40  16.33
High     19.50  28.00
```

### SAS

```sas
/* One row per patient: dose, sex, LDL-C reduction (mg/dL) */
data ldl;
  input drug $ sex $ y;
  datalines;
Placebo F 6
Placebo F 17
Placebo F 4
Placebo F 10
Placebo F 19
Placebo F 7
Placebo F 8
Placebo F 7
Placebo M -5
Placebo M 7
Placebo M 5
Placebo M 11
Low F 15
Low F 17
Low F 22
Low F 12
Low F 11
Low M 17
Low M 19
Low M 16
Low M 16
Low M 11
Low M 19
High F 21
High F 17
High F 21
High F 19
High M 31
High M 21
High M 20
High M 26
High M 30
High M 35
High M 30
High M 31
;
run;

/* Two-way ANOVA with interaction; request both Type I and Type III SS */
proc glm data=ldl;
  class drug sex;
  model y = drug|sex / ss1 ss3;
  lsmeans drug*sex;
run;
quit;
```

**What the code does**

- `class drug sex;` declares both factors categorical. `drug|sex` is shorthand for `drug sex drug*sex` (both main effects and the interaction).
- `/ ss1 ss3` prints **Type I** and **Type III** tables. PROC GLM prints Type I and Type III by default; the options just make the choice explicit. **Type III is the one most analyses and submissions use.** Its values match `car::Anova(type = 3)` in R and `typ=3` with sum contrasts in Python.
- The Type I table follows the order in the `MODEL` statement (`drug` before `sex`), matching R's first `anova()` table.
- `lsmeans drug*sex;` prints the cell means (equal to the raw cell means here because the model is saturated). For marginal means of one factor, `lsmeans drug` gives **unweighted** (least-squares) means, which is what the Type III test compares.
- SAS parameterizes classification effects internally with dummy variables, but the Type III *tests* do not depend on that choice. This is why SAS needs no `contrasts` setting.

## Interpreting the output

- **Interaction (drug × sex):** F(2, 29) = 5.51, p = 0.009, and identical in every table. The dose effect depends on sex. Women go from 9.75 (placebo) to 19.50 mg/dL (high dose), a gain of about 9.8 mg/dL; men go from 4.50 to 28.00, a gain of 23.5 mg/dL.
- **Drug main effect:** Type III F(2, 29) = 32.12, p < 0.001 (Type I: 38.57; the differing value reflects the order and imbalance, not an error). Because of the interaction, describe the dose effect within each sex rather than as one number.
- **Sex main effect:** Type III F(1, 29) = 0.68, p = 0.42. Averaged over doses, there is no clear difference between women and men. **Do not read this as "sex does not matter":** the interaction shows women respond better on placebo and men much better on high dose, and these opposite differences cancel in the average. Also, the Type I table with sex entered **first** gives F = 10.23, p = 0.003, a spurious "significant" sex effect that arises because men are over-represented in the high-dose arm (8 of 12 high-dose patients). Type I credits to sex any variation that is shared with dose.
- **Residual:** MS_error = 22.97 (SD ≈ 4.8 mg/dL) on 29 df.
- **Which to report:** the Type III table, with the interaction interpreted first, and the unbalanced cell sizes stated.

**How to report it:**
> LDL-C reduction depended on dose (Type III F(2, 29) = 32.12, p < 0.001) and the dose effect differed by sex (dose × sex interaction F(2, 29) = 5.51, p = 0.009): from placebo to high dose, the mean reduction increased by 9.8 mg/dL in women (9.75 to 19.50) and by 23.5 mg/dL in men (4.50 to 28.00).

## Check your understanding

<details class="check"><summary>1. The drug × sex interaction is significant, and the sex main effect has p = 0.42. A colleague concludes "sex does not matter." What is wrong?</summary>

The main effect is an average over doses. Women do better on placebo (9.75 vs 4.50) but men do far better on high dose (28.00 vs 19.50). These differences of opposite sign cancel in the average, which hides a real sex difference that depends on dose. When an interaction is present, report effects within each level of the other factor.
</details>

<details class="check"><summary>2. Why can the same data give a significant sex effect (p = 0.003) from R's `anova()` and a non-significant one (p = 0.42) from `car::Anova(type = 3)`?</summary>

With unbalanced cells, Type I (sequential) sums of squares credit to sex whatever variation it shares with dose, because it was entered first. Men are over-represented in the high-dose arm, so sex picks up part of the dose effect. Type III tests each term after all the others, removing that overlap. Pre-specify Type III (with sum-to-zero contrasts in R and Python) in the analysis plan.
</details>

<details class="check"><summary>3. Three doses are tested in men and women, and the interaction p is 0.30. Can you now say the drug works the same in both sexes?</summary>

No. A non-significant interaction test is often low-powered, especially in trials sized for the main effect. It means the data do not clearly show a difference in the dose effect, not that none exists. Look at the estimated cell means and CIs, and decide whether to keep the interaction in the model in the analysis plan, not after seeing the results.
</details>

<details class="check"><summary>4. What if the men's line were exactly parallel to the women's, with the same placebo-to-high-dose gain of 9.75 mg/dL? What would happen to the interaction F?</summary>

The interaction F would drop towards 0 because the interaction sum of squares measures how far the cell means depart from parallel lines. The dose and sex main effects could then be interpreted more directly as average effects.
</details>

## Common mistakes

- **Using R's `anova()` on unbalanced data and calling it "the ANOVA."** It gives Type I SS; the main-effect results depend on the order of the terms. Use `car::Anova(type = 3)` with sum-to-zero contrasts to match SAS.
- **Forgetting `contr.sum` (R) or `Sum` coding (Python)** when requesting Type III. The table is then silently wrong for main effects.
- **Interpreting main effects when the interaction is significant.** Average effects can hide opposite or very different effects in subgroups; show the cell means or an interaction plot.
- **Dropping a non-significant interaction automatically.** A non-significant interaction test is low-powered; decide the model in the analysis plan.
- **Ignoring unequal variances.** If cell SDs differ widely, consider a variance-stabilizing transformation or a heteroscedasticity-robust approach.
- **Treating a post-hoc subgroup finding as confirmatory.** Sex-by-treatment interactions should be pre-specified, and interaction tests have low power in trials sized for the main effect.

## Related tests

- [One-way ANOVA](#/test/one-way-anova): the single-factor version
- [ANCOVA](#/test/ancova): adds a continuous covariate such as baseline
- [Tukey's HSD](#/test/tukey-hsd): all pairwise comparisons between cell or level means
- [Levene's test](#/test/levene): checking equal variances
- [Linear regression](#/test/linear-regression): the same model written with dummy variables
- [Linear mixed-effects models (MMRM)](#/test/mixed-models): when the factors include repeated measures

## Key takeaways

- Two-way ANOVA tests **two main effects** and their **interaction** with three F ratios.
- An interaction means the effect of one factor **depends on** the level of the other (non-parallel lines).
- Test and interpret the interaction first. If present, describe effects within levels, not as one average.
- With unbalanced cells, Type I SS depend on term order. Report Type III with sum-to-zero contrasts.
- A non-significant interaction is not proof that the effect is the same in every subgroup.

## References

- Searle SR, Gruber MJ. *Linear Models.* 2nd ed. Wiley; 2017. On Type I–III sums of squares.
- Fox J, Weisberg S. *An R Companion to Applied Regression.* 3rd ed. Sage; 2019. Section on Type II and III tests and contrasts.
- R documentation: `?car::Anova`, `?stats::anova.lm`
- statsmodels documentation: `statsmodels.stats.anova.anova_lm`
- SAS/STAT User's Guide: The GLM Procedure — Specification of Effects, and the MODEL statement SS1 and SS3 options
