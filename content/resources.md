# Further learning

> Hand-picked resources to go deeper: start with the intuition-builders, then move to textbooks and the clinical-research classics.

## Start here: build intuition

- **[Khan Academy: Statistics and probability](https://www.khanacademy.org/math/statistics-probability)**: short visual lessons on sampling distributions, significance tests, and confidence intervals, with practice questions after each one. Best for a first pass at hypothesis testing.
- **[StatQuest with Josh Starmer](https://statquest.org)**: clear, picture-driven videos on p-values, t-tests, linear models, logistic regression, and more. Each topic is explained intuitively before any formula appears.
- **[Penn State STAT 415 and STAT 500 course notes](https://online.stat.psu.edu)**: free university course notes that are thorough and well organized, with learning objectives, derivations, and software examples.
- **[Seeing Theory (Brown University)](https://seeing-theory.brown.edu)**: interactive visualizations of probability, sampling, and inference. Drag the sliders and watch confidence intervals and p-values behave.
- **[OpenIntro Statistics](https://www.openintro.org/book/os/)**: a free, peer-reviewed textbook with exercises and real datasets that covers all the major tests.
- **[UCLA OARC: What statistical analysis should I use?](https://stats.oarc.ucla.edu/other/mult-pkg/whatstat/)**: a classic table that maps study designs to analyses, with worked examples in R, SAS, Stata, and SPSS.

## Textbooks

| Book | Level | Why read it |
|---|---|---|
| *The Art of Statistics*, David Spiegelhalter | Beginner | A modern, very readable account of how statistics answers real questions, with many health examples |
| *Introduction to the Practice of Statistics*, Moore, McCabe & Craig | Introductory | The classic first course: data, sampling, inference, and the core tests |
| *OpenIntro Statistics* (4th ed.), Diez, Çetinkaya-Rundel & Barr | Introductory | Free PDF, ideal for self-study, with exercises and real data |
| *Discovering Statistics Using R* / *…Using IBM SPSS*, Andy Field | Introductory–applied | Very practical and friendly; great for applied researchers |
| *Practical Statistics for Medical Research*, Douglas Altman | Intermediate | The standard medical-statistics text, focused on choosing and interpreting analyses |
| *Essential Medical Statistics*, Kirkwood & Sterne | Intermediate | A concise reference for epidemiology and clinical research methods |
| *Regression Modeling Strategies*, Frank Harrell | Advanced | How to build, validate, and report clinical prediction and regression models |
| *Data Analysis Using Regression and Multilevel/Hierarchical Models*, Gelman & Hill | Advanced | Regression and multilevel modeling with an emphasis on model checking; bridges classical and Bayesian approaches |
| *All of Statistics*, Larry Wasserman | Advanced | A compact, mathematically precise tour of inference |
| *Statistical Inference*, Casella & Berger | Advanced | The rigorous graduate-level theory behind estimation and testing |
| *Testing Statistical Hypotheses*, Lehmann & Romano | Advanced | The definitive theoretical treatment of hypothesis testing |

## Foundational papers

| Paper | Why it matters |
|---|---|
| Fisher (1925). *Statistical Methods for Research Workers* | Introduced significance testing and the p-value as we use them today |
| Neyman & Pearson (1933). *On the problem of the most efficient tests of statistical hypotheses.* Phil Trans R Soc A 231:289–337 | Formalized Type I and II errors, power, and critical regions. See [errors and power](#/test/errors-and-power) |
| Cohen (1994). *The earth is round (p < .05).* American Psychologist 49:997–1003 | A classic critique of mechanical p-value use; argues for effect sizes and power. See [effect sizes](#/test/effect-sizes) |
| Wasserstein & Lazar (2016). *The ASA statement on p-values.* The American Statistician 70:129–133 | The American Statistical Association's official guidance on interpreting p-values. See [p-values](#/test/hypothesis-testing) |
| Ioannidis (2005). *Why most published research findings are false.* PLOS Medicine 2:e124 | Shows how low power, bias, and multiple testing inflate false discoveries |

## Clinical research essentials

- **[BMJ Statistics Notes (Altman & Bland)](https://www.bmj.com/specialties/statistics-notes)**: short, authoritative explanations of single concepts, such as why "absence of evidence is not evidence of absence".
- **[ICH E9: Statistical Principles for Clinical Trials](https://database.ich.org/sites/default/files/E9_Guideline.pdf)** and **[ICH E9(R1) on estimands](https://database.ich.org/sites/default/files/E9-R1_Step4_Guideline_2019_1203.pdf)**: the regulatory foundation for trial analysis plans.
- **[CONSORT](https://www.consort-spirit.org/)**: how to report randomized trials.
- **[SAMPL guidelines](https://www.equator-network.org/reporting-guidelines/sampl/)**: how to report statistical methods and results in biomedical journals.
- **[EQUATOR Network](https://www.equator-network.org/)**: reporting guidelines for every study design (STROBE, STARD, TRIPOD, PRISMA, and more).

## Software

- **R:** *[R for Data Science](https://r4ds.hadley.nz/)* (free online) for data handling, and the CRAN Task View for **[Clinical Trials](https://cran.r-project.org/web/views/ClinicalTrials.html)**.
- **Python:** the **[SciPy stats](https://docs.scipy.org/doc/scipy/reference/stats.html)**, **[statsmodels](https://www.statsmodels.org/)**, and **[lifelines](https://lifelines.readthedocs.io/)** documentation. Each includes worked examples.

## Free tools

- **[R](https://www.r-project.org)**: `t.test()`, `chisq.test()`, `aov()`, and thousands of packages. Every lesson here has R code.
- **[Python: SciPy](https://scipy.org) and [statsmodels](https://www.statsmodels.org)**: `scipy.stats` for classic tests, and statsmodels for regression and ANOVA.
- **[G*Power](https://www.psychologie.hhu.de/arbeitsgruppen/allgemeine-psychologie-und-arbeitspsychologie/gpower)**: a free desktop tool for power and sample-size calculations.
- **[JASP](https://jasp-stats.org)**: free point-and-click software for both frequentist and Bayesian tests. Good for checking your R or Python results.

## Research and papers

- **[Google Scholar](https://scholar.google.com)**: search any method by name; the "Cited by" links help you find tutorials and applications.
- **[PubMed](https://pubmed.ncbi.nlm.nih.gov/)**: find tutorial papers. Try searching "<test name> tutorial" with the filter *Review*.
- **[JSTOR](https://www.jstor.org)**: classic statistics papers, such as Cox (1972) and Kaplan & Meier (1958).
- **[arXiv: Statistics](https://arxiv.org/list/stat/recent)**: free preprints of new methods research.
