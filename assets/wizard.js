/*
 * "Which test?" decision tree.
 * Each node is a question with options. An option either points to the next
 * node (`next`) or ends in a recommendation (`result`): a list of
 * { slug, why } where the first entry is the primary recommendation.
 */
(function () {
  const r = (slug, why) => ({ slug, why });

  // Leaves for "compare groups", indexed by outcome → groups → design
  const COMPARE = {
    normal: {
      "1":  { ind: [r("one-sample-t", "One group's mean vs a known reference value."), r("one-sample-z", "Only if the population SD is truly known (e.g., standardized z-scores); otherwise use the t-test.")] },
      "2":  { ind: [r("independent-t", "Two independent groups; use Welch's version by default."), r("ancova", "If you have a baseline measurement (e.g., an RCT), ANCOVA adjusting for baseline is more powerful.")],
              pair: [r("paired-t", "Two measurements on the same subjects (or matched pairs)."), r("bland-altman", "If the goal is agreement between two methods rather than a change, use Bland–Altman.")] },
      "3+": { ind: [r("one-way-anova", "Three or more independent groups."), r("welch-anova", "If group SDs or sizes differ noticeably, use Welch's ANOVA (with Games–Howell post-hoc)."), r("tukey-hsd", "Follow up a significant ANOVA with all pairwise comparisons."), r("dunnett", "If every arm is compared only with a control, use Dunnett."), r("two-way-anova", "If there is a second factor (e.g., sex or site), use a two-way ANOVA.")],
              pair: [r("mixed-models", "Repeated visits per subject; handles missing visits under MAR. Preferred in trials."), r("repeated-measures-anova", "Classic option when every subject has complete data at every time point.")] },
    },
    skewed: {
      "1":  { ind: [r("wilcoxon-signed-rank", "The one-sample signed-rank test compares the distribution's center to a reference value."), r("one-proportion", "If you dichotomize the outcome, use an exact binomial test.")] },
      "2":  { ind: [r("mann-whitney", "Rank-based comparison of two independent groups."), r("ordinal-regression", "For an ordinal scale such as mRS, a proportional-odds model gives a common odds ratio and can adjust for covariates.")],
              pair: [r("wilcoxon-signed-rank", "Rank-based test on within-subject differences."), r("sign-test", "If the differences are not symmetric, the sign test only counts ups vs downs.")] },
      "3+": { ind: [r("kruskal-wallis", "Rank-based comparison of three or more independent groups (follow with Dunn's test)."), r("ordinal-regression", "For an ordinal outcome, a proportional-odds model with covariates.")],
              pair: [r("friedman", "Rank-based comparison of three or more repeated conditions.")] },
    },
    binary: {
      "1":  { ind: [r("one-proportion", "Observed rate vs a target rate, with an exact Clopper–Pearson CI.")] },
      "2":  { ind: [r("two-proportion-z", "Large samples: risk difference with a CI."), r("fisher-exact", "Small samples (any expected cell count < 5)."), r("chi-square-independence", "The equivalent chi-square approach for a 2 × 2 table."), r("cmh", "If randomization is stratified (e.g., by site), use CMH."), r("logistic-regression", "To adjust for covariates, use logistic regression.")],
              pair: [r("mcnemar", "Paired yes/no outcomes on the same subjects.")] },
      "3+": { ind: [r("chi-square-independence", "An r × c contingency table."), r("trend-tests", "If the groups are ordered (e.g., doses), a trend test is more powerful."), r("fisher-exact", "Use the exact version (Fisher–Freeman–Halton) when counts are small."), r("logistic-regression", "To adjust for covariates.")],
              pair: [r("cochrans-q", "Same subjects under 3+ conditions with a yes/no outcome."), r("gee", "Repeated binary outcomes with covariates: population-averaged odds ratios.")] },
    },
    // 3+ unordered categories (ordered categories such as disease stage go to "skewed")
    nominal: {
      "1":  { ind: [r("chi-square-gof", "Do the observed category proportions match expected ones (e.g., a reference population)?")] },
      "2":  { ind: [r("chi-square-independence", "A 2 × k table: compare the distribution across categories between the two groups."), r("fisher-exact", "Small expected counts: the Fisher–Freeman–Halton extension (fisher.test in R handles r × c tables)."), r("multinomial-logistic", "To adjust for covariates: one odds ratio per category vs a reference category.")],
              pair: [r("mcnemar", "The same subjects classified twice into k categories: R's mcnemar.test on a k × k table runs Bowker's test of symmetry (the McNemar extension).")] },
      "3+": { ind: [r("chi-square-independence", "An r × c table: groups by categories."), r("fisher-exact", "Small expected counts: the Fisher–Freeman–Halton extension."), r("multinomial-logistic", "To adjust for covariates.")],
              pair: [r("gee", "Repeated unordered categorical outcomes need a multinomial GEE or mixed model (e.g., R package multgee); worth a statistician's input.")] },
    },
    count: {
      "1":  { ind: [r("poisson-regression", "Estimate an event rate with a CI (an intercept-only model with an exposure offset).")] },
      "2":  { ind: [r("poisson-regression", "Rate ratio between groups; switch to negative binomial if overdispersed.")],
              pair: [r("gee", "Repeated counts per subject.")] },
      "3+": { ind: [r("poisson-regression", "Rate ratios across groups, with an offset for follow-up time.")],
              pair: [r("gee", "Repeated counts per subject."), r("mixed-models", "A subject-specific alternative (a generalized mixed model).")] },
    },
  };

  const outcomeNode = {
    id: "outcome",
    q: "What type of outcome are you comparing?",
    options: [
      { key: "normal", label: "Continuous, roughly normal", hint: "e.g., blood pressure, HbA1c, cholesterol", next: "groups" },
      { key: "skewed", label: "Continuous but skewed, or ordinal (ordered categories)", hint: "e.g., length of stay, pain score 0–10, Likert, disease stage", next: "groups" },
      { key: "binary", label: "Binary (yes/no)", hint: "e.g., responder yes/no, adverse event, 30-day mortality", next: "groups" },
      { key: "nominal", label: "Categorical: 3+ unordered categories", hint: "e.g., discharge destination, blood type, tumour subtype", next: "groups" },
      { key: "count",  label: "Count or rate", hint: "e.g., number of exacerbations, events per patient-year", next: "groups" },
      { key: "tte",    label: "Time to an event", hint: "e.g., overall survival, time to relapse (with censoring)", next: "tteGoal" },
    ],
  };

  const NODES = {
    start: {
      q: "What is your research question?",
      options: [
        { label: "Compare groups or time points", hint: "Treatment vs placebo, before vs after, 3 dose arms…", next: "outcome" },
        { label: "Measure an association or adjust for covariates", hint: "Correlation, regression, prediction", next: "assoc" },
        { label: "Check agreement or reliability", hint: "Two raters, two devices, test–retest", next: "agree" },
        { label: "Evaluate a diagnostic test or biomarker", hint: "Sensitivity, specificity, ROC/AUC", result: [r("roc-auc", "ROC analysis summarizes discrimination across all cutoffs; report sensitivity and specificity at your chosen cutoff.")] },
        { label: "Show two treatments are equivalent or not worse", hint: "Generic vs brand, non-inferiority trial", result: [r("tost-equivalence", "Equivalence and non-inferiority require a pre-specified margin; a non-significant difference test does NOT show equivalence.")] },
        { label: "Check assumptions before testing", hint: "Normality, equal variances", next: "assume" },
        { label: "Plan a study: how many participants do I need?", hint: "Sample size and power", next: "plan" },
        { label: "Something else", hint: "Combine studies, missing data, observational comparisons, resampling", next: "special" },
        { label: "I'm new: teach me the basics first", hint: "p-values, confidence intervals, power, effect sizes", result: [r("hypothesis-testing", "Start here: what a p-value really means."), r("confidence-intervals", "Then: reading a 95% CI."), r("errors-and-power", "Then: false positives, false negatives, and power."), r("effect-sizes", "Finally: statistical vs clinical significance.")] },
      ],
    },
    outcome: outcomeNode,
    groups: {
      q: "How many groups or conditions are you comparing?",
      options: [
        { key: "1",  label: "One group vs a reference value", hint: "e.g., is mean HbA1c different from 7%?", next: "design" },
        { key: "2",  label: "Two groups or conditions", hint: "e.g., drug vs placebo, baseline vs week 12", next: "design" },
        { key: "3+", label: "Three or more", hint: "e.g., 3 dose arms, 4 visits", next: "design" },
      ],
    },
    design: {
      q: "Are the observations independent or paired/repeated?",
      skipIf: (a) => a.groups === "1",
      options: [
        { key: "ind",  label: "Independent", hint: "Different subjects in each group", leaf: true },
        { key: "pair", label: "Paired or repeated", hint: "Same subjects measured more than once, or matched pairs", leaf: true },
      ],
    },
    tteGoal: {
      q: "What do you want to do with the time-to-event data?",
      options: [
        { label: "Describe survival (curves, median, 1-year survival)", result: [r("kaplan-meier", "The non-parametric survival curve, which handles censoring.")] },
        { label: "Compare survival between groups", result: [r("log-rank", "Standard unadjusted comparison of survival curves."), r("kaplan-meier", "Always plot KM curves alongside the test.")] },
        { label: "Estimate an adjusted hazard ratio", result: [r("cox-ph", "Hazard ratios adjusted for covariates; check proportional hazards.")] },
        { label: "Another event can prevent the event of interest", hint: "e.g., death before relapse", result: [r("competing-risks", "Use cumulative incidence functions, not 1 − KM; Fine–Gray for absolute risk, cause-specific Cox for etiology.")] },
      ],
    },
    assoc: {
      q: "What does your outcome look like?",
      options: [
        { label: "Two continuous variables, no outcome/predictor distinction", hint: "e.g., BMI and systolic BP", next: "corr" },
        { label: "Continuous outcome with predictors", result: [r("linear-regression", "Estimate adjusted mean differences/slopes with CIs.")] },
        { label: "Binary outcome with predictors", result: [r("logistic-regression", "Adjusted odds ratios. Watch events-per-variable."), r("firth-logistic", "Few events, a small sample, or separation (huge SEs): use Firth."), r("conditional-logistic", "Matched case–control design: use conditional logistic regression."), r("effect-measures", "Report RD, RR, and NNT alongside the OR.")] },
        { label: "Unordered categorical outcome (3+ levels)", hint: "e.g., discharge destination", result: [r("multinomial-logistic", "One set of coefficients per category vs a reference.")] },
        { label: "Several continuous outcomes together", result: [r("manova", "Tests group differences on a set of correlated outcomes."), r("mixed-models", "If the outcomes are repeated measures of the same variable, use a mixed model instead.")] },
        { label: "Ordinal outcome with predictors", hint: "e.g., modified Rankin Scale, NYHA class", result: [r("ordinal-regression", "A proportional-odds model uses the full ordered scale instead of dichotomizing it.")] },
        { label: "Count or rate outcome", result: [r("poisson-regression", "Rate ratios with a person-time offset; check overdispersion.")] },
        { label: "Repeated measures or clustered data", hint: "Visits within patients, patients within sites", result: [r("mixed-models", "Subject-specific effects; handles unbalanced/missing visits."), r("gee", "Population-averaged effects with robust SEs.")] },
        { label: "Time-to-event outcome", result: [r("cox-ph", "Adjusted hazard ratios.")] },
      ],
    },
    corr: {
      q: "Is the relationship linear, and are both variables roughly normal?",
      options: [
        { label: "Yes, linear and roughly normal", result: [r("pearson", "Strength of linear association, with a CI for r.")] },
        { label: "No, monotonic, skewed, ordinal, or with outliers", result: [r("spearman", "Rank-based; robust to outliers and non-linearity (if monotonic)."), r("kendall-tau", "Many ties or a small sample: Kendall's tau.")] },
      ],
    },
    agree: {
      q: "What kind of measurements are being compared?",
      options: [
        { label: "Categorical ratings by two raters", hint: "e.g., two radiologists: malignant / benign", result: [r("cohens-kappa", "Agreement beyond chance; use weighted kappa for ordinal ratings.")] },
        { label: "Continuous measurements: reliability across raters or repeats", hint: "e.g., 3 nurses measure wound size", result: [r("icc", "Proportion of variance due to true differences between subjects.")] },
        { label: "Items of a questionnaire or scale", hint: "e.g., a 6-item quality-of-life score", result: [r("cronbach-alpha", "Internal consistency of the items.")] },
        { label: "Continuous measurements: two methods or devices", hint: "e.g., new cuff vs reference sphygmomanometer", result: [r("bland-altman", "Bias and 95% limits of agreement; correlation is not agreement.")] },
      ],
    },
    plan: {
      q: "What is the primary endpoint of the study you're planning?",
      options: [
        { label: "Continuous or binary endpoint", hint: "e.g., change in SBP, response rate", result: [r("sample-size-means-proportions", "Sample size for comparing two means or two proportions; inflate for dropout.")] },
        { label: "Time-to-event endpoint", hint: "e.g., overall survival, progression-free survival", result: [r("sample-size-survival", "The number of events drives power (Schoenfeld); then convert events to patients.")] },
      ],
    },
    special: {
      q: "Which of these describes your situation?",
      options: [
        { label: "Combine results from several studies", result: [r("meta-analysis", "Fixed or random effects, with heterogeneity and a forest plot.")] },
        { label: "Some outcome or covariate values are missing", result: [r("multiple-imputation", "Complete-case analysis can be biased; multiple imputation under MAR.")] },
        { label: "Compare treatments in observational (non-randomized) data", result: [r("propensity-scores", "Balance measured confounders by matching or weighting.")] },
        { label: "Assumptions doubtful, or no formula for the CI", result: [r("bootstrap-permutation", "Resampling-based tests and confidence intervals.")] },
      ],
    },
    assume: {
      q: "Which assumption do you want to check?",
      options: [
        { label: "Normality", result: [r("shapiro-wilk", "A formal test; always pair it with a Q–Q plot, especially for large n.")] },
        { label: "Equal variances across groups", result: [r("levene", "Robust to non-normality (the Brown–Forsythe median version); or just use Welch's t-test.")] },
        { label: "I'm running many tests", result: [r("multiple-comparisons", "Control the family-wise error rate or the false discovery rate.")] },
      ],
    },
  };

  function resolveCompare(a) {
    const design = a.groups === "1" ? "ind" : a.design;
    return COMPARE[a.outcome]?.[a.groups]?.[design] || [];
  }

  window.WIZARD = { NODES, resolveCompare };
})();
