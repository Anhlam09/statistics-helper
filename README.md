# Biostat Test Guide

A free, friendly guide to choosing and running statistical tests in clinical research. It has 62 lessons, each with a plain-language explanation, a clinical example, the assumptions to check, and runnable **R** and **Python** code. Some lessons also include **SAS** code, which is marked as not yet tested.

**Live site:** https://anhlam09.github.io/biostat-test-guide/

> **Educational resource.** Always pre-specify your analysis and consult a statistician for trial work.

## Run it locally

The site is static (HTML, CSS and JavaScript, with no build step). Pages are loaded with `fetch`, so serve the folder instead of opening `index.html` directly:

```bash
python3 -m http.server 8765
```

Then open http://localhost:8765.

## Structure

```
index.html              App shell
assets/                 Styles, router, plush mascots, café, claw machine, calculators, diagrams
content/manifest.json   Categories and tests (drives navigation, search and routing)
content/tests/*.md      One lesson per test (paired-t.md is the reference lesson)
content/about.md        About page
CONTENT_GUIDE.md        Rules for writing a lesson
tools/                  Content, calculator and code-example checks
```

## Checks

```bash
python3 tools/check_content.py     # lesson structure
node tools/test_calculators.js     # calculators vs R reference values
node tools/viz_pvalues.js          # p-values shown in diagrams
python3 tools/run_examples.py      # runs every R and Python example (needs R and Python packages)
```

## Author

Made by Anh Lam.
