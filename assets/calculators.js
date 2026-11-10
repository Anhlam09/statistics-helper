/* Interactive calculators. Pure functions live in window.CALC_COMPUTE (testable in node);
   DOM code lives in the mount() functions below and is skipped when there is no document. */
(function (root) {
  "use strict";
  const SM = root.StatMath || (typeof require !== "undefined" ? require("./statmath.js") : null);
  if (!SM) throw new Error("calculators.js needs statmath.js");

  /* =====================================================================
     Pure computation
     ===================================================================== */
  const zq = (p) => SM.normQuantile(p);
  const isInt = (x) => Number.isFinite(x) && Math.floor(x) === x;
  const fin = (...xs) => xs.every(Number.isFinite);

  function ciBounds(est, se, q) { return [est - q * se, est + q * se]; }

  /* a. Two-sample t-test from summary statistics */
  function ttestInd(i) {
    const { m1, s1, n1, m2, s2, n2 } = i;
    const welch = i.welch !== false;
    const conf = i.conf == null ? 0.95 : i.conf;
    const err = [];
    if (!fin(m1, s1, n1, m2, s2, n2)) err.push("Enter mean, SD and n for both groups.");
    else {
      if (s1 < 0 || s2 < 0) err.push("SD cannot be negative.");
      if (!isInt(n1) || !isInt(n2) || n1 < 2 || n2 < 2) err.push("Each group needs a whole-number n of at least 2.");
      if (!err.length && s1 === 0 && s2 === 0) err.push("At least one SD must be greater than 0.");
    }
    if (err.length) return { errors: err };
    const diff = m1 - m2;
    const v1 = s1 * s1, v2 = s2 * s2;
    const sp2 = ((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2);
    let se, df;
    if (welch) {
      const a = v1 / n1, b = v2 / n2;
      se = Math.sqrt(a + b);
      df = ((a + b) * (a + b)) / ((a * a) / (n1 - 1) + (b * b) / (n2 - 1));
    } else {
      se = Math.sqrt(sp2 * (1 / n1 + 1 / n2));
      df = n1 + n2 - 2;
    }
    const t = diff / se;
    const p = SM.tTwoTail(t, df);
    const q = SM.tQuantile(1 - (1 - conf) / 2, df);
    const [ci_lo, ci_hi] = ciBounds(diff, se, q);
    return { welch, diff, se, t, df, p, ci_lo, ci_hi, d: diff / Math.sqrt(sp2), sp: Math.sqrt(sp2) };
  }

  /* b. Paired / one-sample t from summary statistics */
  function ttestPaired(i) {
    const { mean, sd, n } = i;
    const mu0 = i.mu0 == null || Number.isNaN(i.mu0) ? 0 : i.mu0;
    const err = [];
    if (!fin(mean, sd, n)) err.push("Enter the mean, SD and n.");
    else {
      if (!(sd > 0)) err.push("SD must be greater than 0.");
      if (!isInt(n) || n < 2) err.push("n must be a whole number of at least 2.");
    }
    if (err.length) return { errors: err };
    const se = sd / Math.sqrt(n), df = n - 1, t = (mean - mu0) / se;
    const q = SM.tQuantile(0.975, df);
    const [ci_lo, ci_hi] = ciBounds(mean, se, q);
    return { mean, mu0, est: mean, se, t, df, p: SM.tTwoTail(t, df), ci_lo, ci_hi, dz: (mean - mu0) / sd };
  }

  /* Fisher exact two-sided p, R's rule: sum of P(x) <= P(obs) * (1 + 1e-7) */
  function fisherP(a, b, c, d) {
    const m = a + b, k = a + c, N = a + b + c + d;
    const lo = Math.max(0, k - (N - m)), hi = Math.min(m, k);
    const dobs = SM.hyperPmf(a, k, N, m);
    const lim = dobs * (1 + 1e-7);
    let sum = 0;
    for (let x = lo; x <= hi; x++) {
      const px = SM.hyperPmf(x, k, N, m);
      if (px <= lim) sum += px;
    }
    return Math.min(1, sum);
  }

  /* c. 2x2 table */
  function twoByTwo(i) {
    const { a, b, c, d } = i;
    if (![a, b, c, d].every((x) => Number.isFinite(x))) return { errors: ["Enter all four cell counts."] };
    if (![a, b, c, d].every((x) => isInt(x) && x >= 0)) return { errors: ["Cell counts must be whole numbers, 0 or more."] };
    const N = a + b + c + d, r1 = a + b, r2 = c + d, c1 = a + c, c2 = b + d;
    if (N === 0) return { errors: ["The table is empty."] };
    const out = { N, fisher_p: fisherP(a, b, c, d) };
    if (r1 > 0 && r2 > 0 && c1 > 0 && c2 > 0) {
      const ad = Math.abs(a * d - b * c);
      out.chi = (N * ad * ad) / (r1 * r2 * c1 * c2);
      out.chi_p = SM.chisqSf(out.chi, 1);
      const dev = ad / N, y = Math.min(0.5, dev);
      const E = [r1 * c1 / N, r1 * c2 / N, r2 * c1 / N, r2 * c2 / N];
      out.yates = E.reduce((s, e) => s + ((dev - y) * (dev - y)) / e, 0);
      out.yates_p = SM.chisqSf(out.yates, 1);
      out.minE = Math.min(...E);
    } else out.chi_note = "Chi-square is undefined because a row or column total is 0.";
    const z = zq(0.975);
    if (r1 > 0 && r2 > 0) {
      const p1 = a / r1, p2 = c / r2;
      out.p1 = p1; out.p2 = p2;
      out.rd = p1 - p2;
      const se = Math.sqrt((p1 * (1 - p1)) / r1 + (p2 * (1 - p2)) / r2);
      out.rd_lo = out.rd - z * se; out.rd_hi = out.rd + z * se;
    }
    const zero = a === 0 || b === 0 || c === 0 || d === 0;
    out.adjusted = zero;
    const k = zero ? 0.5 : 0;
    const A = a + k, B = b + k, C = c + k, D = d + k;
    const R1 = A + B, R2 = C + D;
    out.rr = (A / R1) / (C / R2);
    const serr = Math.sqrt(1 / A - 1 / R1 + 1 / C - 1 / R2);
    out.rr_lo = Math.exp(Math.log(out.rr) - z * serr); out.rr_hi = Math.exp(Math.log(out.rr) + z * serr);
    out.or = (A * D) / (B * C);
    const seo = Math.sqrt(1 / A + 1 / B + 1 / C + 1 / D);
    out.or_lo = Math.exp(Math.log(out.or) - z * seo); out.or_hi = Math.exp(Math.log(out.or) + z * seo);
    return out;
  }

  /* R's binom.test two-sided p-value */
  function binomTwoSided(x, n, p) {
    if (p === 0) return x === 0 ? 1 : 0;
    if (p === 1) return x === n ? 1 : 0;
    if (p === 0.5) return Math.min(1, 2 * Math.min(SM.binomCdf(x, n, p), SM.binomSf(x - 1, n, p)));
    const relErr = 1 + 1e-7;
    const dd = SM.binomPmf(x, n, p), m = n * p;
    if (x === m) return 1;
    let y = 0;
    if (x < m) {
      for (let i = Math.ceil(m); i <= n; i++) if (SM.binomPmf(i, n, p) <= dd * relErr) y++;
      return Math.min(1, SM.binomCdf(x, n, p) + SM.binomSf(n - y, n, p));
    }
    for (let i = 0; i <= Math.floor(m); i++) if (SM.binomPmf(i, n, p) <= dd * relErr) y++;
    return Math.min(1, SM.binomCdf(y - 1, n, p) + SM.binomSf(x - 1, n, p));
  }

  function clopperPearson(x, n, conf) {
    const al = 1 - (conf == null ? 0.95 : conf);
    return [
      x === 0 ? 0 : SM.betaQuantile(al / 2, x, n - x + 1),
      x === n ? 1 : SM.betaQuantile(1 - al / 2, x + 1, n - x),
    ];
  }
  function wilson(x, n, conf) {
    const z = zq(1 - (1 - (conf == null ? 0.95 : conf)) / 2);
    const p = x / n, z2 = z * z;
    const centre = (p + z2 / (2 * n)) / (1 + z2 / n);
    const half = (z / (1 + z2 / n)) * Math.sqrt((p * (1 - p)) / n + z2 / (4 * n * n));
    return [Math.max(0, centre - half), Math.min(1, centre + half)];
  }

  /* d. McNemar */
  function mcnemar(i) {
    const { b, c } = i;
    if (!fin(b, c) || !isInt(b) || !isInt(c) || b < 0 || c < 0) return { errors: ["b and c must be whole numbers, 0 or more."] };
    const n = b + c;
    const out = { b, c, n };
    if (n === 0) { out.exact_p = 1; out.note = "No discordant pairs: the test is undefined (exact p reported as 1)."; return out; }
    out.exact_p = binomTwoSided(b, n, 0.5);
    out.chi = ((b - c) * (b - c)) / n;
    out.chi_p = SM.chisqSf(out.chi, 1);
    const yates = b !== c ? 1 : 0;
    out.chi_cc = ((Math.abs(b - c) - yates) ** 2) / n;
    out.chi_cc_p = SM.chisqSf(out.chi_cc, 1);
    return out;
  }

  /* e. One proportion */
  function oneProportion(i) {
    const { x, n, p0 } = i;
    if (!fin(x, n, p0)) return { errors: ["Enter x, n and the null proportion."] };
    if (!isInt(n) || n < 1) return { errors: ["n must be a whole number of at least 1."] };
    if (!isInt(x) || x < 0 || x > n) return { errors: ["x must be a whole number between 0 and n."] };
    if (!(p0 > 0 && p0 < 1)) return { errors: ["The null proportion must be strictly between 0 and 1."] };
    const [cp_lo, cp_hi] = clopperPearson(x, n);
    const [wilson_lo, wilson_hi] = wilson(x, n);
    return { phat: x / n, p_exact: binomTwoSided(x, n, p0), cp_lo, cp_hi, wilson_lo, wilson_hi };
  }

  /* f. Diagnostic accuracy */
  function diagnostic(i) {
    const { tp, fp, fn, tn } = i;
    if (![tp, fp, fn, tn].every(Number.isFinite)) return { errors: ["Enter TP, FP, FN and TN."] };
    if (![tp, fp, fn, tn].every((v) => isInt(v) && v >= 0)) return { errors: ["Counts must be whole numbers, 0 or more."] };
    const prev = i.prev == null || Number.isNaN(i.prev) ? null : i.prev;
    if (prev !== null && !(prev > 0 && prev < 1)) return { errors: ["Prevalence must be between 0 and 100 (exclusive)."] };
    if (tp + fp + fn + tn === 0) return { errors: ["The table is empty."] };
    const prop = (x, n) => {
      if (n === 0) return null;
      const [lo, hi] = clopperPearson(x, n);
      return { est: x / n, lo, hi, x, n };
    };
    const sens = prop(tp, tp + fn), spec = prop(tn, tn + fp);
    const out = {
      sens, spec, ppv: prop(tp, tp + fp), npv: prop(tn, tn + fn), acc: prop(tp + tn, tp + fp + fn + tn),
    };
    if (sens && spec) {
      out.lr_pos = spec.est === 1 ? Infinity : sens.est / (1 - spec.est);
      out.lr_neg = spec.est === 0 ? Infinity : (1 - sens.est) / spec.est;
      if (prev !== null) {
        const dp = sens.est * prev + (1 - spec.est) * (1 - prev);
        const dn = (1 - sens.est) * prev + spec.est * (1 - prev);
        out.ppv_prev = dp > 0 ? (sens.est * prev) / dp : null;
        out.npv_prev = dn > 0 ? (spec.est * (1 - prev)) / dn : null;
      }
    }
    return out;
  }

  /* g. Sample size, two means: R power.t.test (two-sample, two-sided, strict = FALSE) */
  function powerT(n, delta, sd, alpha) {
    const nu = (n - 1) * 2;
    const qu = SM.tQuantile(1 - alpha / 2, nu);
    return 1 - SM.ntCdf(qu, nu, (Math.sqrt(n / 2) * Math.abs(delta)) / sd);
  }
  function ssMeans(i) {
    const { delta, sd, alpha, power } = i;
    const dropout = i.dropout == null || Number.isNaN(i.dropout) ? 0 : i.dropout;
    const err = [];
    if (!fin(delta, sd, alpha, power)) err.push("Enter delta, SD, alpha and power.");
    else {
      if (delta === 0) err.push("The difference to detect (delta) cannot be 0.");
      if (!(sd > 0)) err.push("SD must be greater than 0.");
      if (!(alpha > 0 && alpha < 1)) err.push("Alpha must be between 0 and 1.");
      if (!(power > alpha && power < 1)) err.push("Power must be between alpha and 1.");
      if (!(dropout >= 0 && dropout < 100)) err.push("Dropout must be at least 0 and below 100.");
    }
    if (err.length) return { errors: err };
    let lo = 2, hi = 4;
    if (powerT(lo, delta, sd, alpha) >= power) hi = lo;
    else {
      while (powerT(hi, delta, sd, alpha) < power && hi < 1e9) { lo = hi; hi *= 2; }
      for (let k = 0; k < 200 && hi - lo > 1e-11 * hi; k++) {
        const mid = (lo + hi) / 2;
        if (powerT(mid, delta, sd, alpha) < power) lo = mid; else hi = mid;
      }
    }
    const nRaw = hi === 2 ? 2 : (lo + hi) / 2;
    const n = Math.ceil(nRaw);
    const out = { n_raw: nRaw, n, total: 2 * n, power_achieved: powerT(n, delta, sd, alpha) };
    if (dropout > 0) {
      out.n_dropout = Math.ceil(n / (1 - dropout / 100));
      out.total_dropout = 2 * out.n_dropout;
    }
    return out;
  }

  /* h. Sample size, two proportions: R power.prop.test (two-sided, strict = FALSE) */
  function ssProps(i) {
    const { p1, p2, alpha, power } = i;
    const dropout = i.dropout == null || Number.isNaN(i.dropout) ? 0 : i.dropout;
    const err = [];
    if (!fin(p1, p2, alpha, power)) err.push("Enter p1, p2, alpha and power.");
    else {
      if (!(p1 > 0 && p1 < 1) || !(p2 > 0 && p2 < 1)) err.push("p1 and p2 must be strictly between 0 and 1.");
      else if (p1 === p2) err.push("p1 and p2 must differ.");
      if (!(alpha > 0 && alpha < 1)) err.push("Alpha must be between 0 and 1.");
      if (!(power > alpha && power < 1)) err.push("Power must be between alpha and 1.");
      if (!(dropout >= 0 && dropout < 100)) err.push("Dropout must be at least 0 and below 100.");
    }
    if (err.length) return { errors: err };
    const za = zq(1 - alpha / 2), zb = zq(power);
    const pb = (p1 + p2) / 2;
    const num = za * Math.sqrt(2 * pb * (1 - pb)) + zb * Math.sqrt(p1 * (1 - p1) + p2 * (1 - p2));
    const nRaw = (num * num) / ((p1 - p2) * (p1 - p2));
    const n = Math.ceil(nRaw);
    const out = { n_raw: nRaw, n, total: 2 * n };
    if (dropout > 0) { out.n_dropout = Math.ceil(n / (1 - dropout / 100)); out.total_dropout = 2 * out.n_dropout; }
    return out;
  }

  /* i. Survival: Schoenfeld */
  function ssSurvival(i) {
    const { hr, alpha, power } = i;
    const p = i.p == null || Number.isNaN(i.p) ? 0.5 : i.p;
    const ep = i.eventProb == null || Number.isNaN(i.eventProb) ? null : i.eventProb;
    const err = [];
    if (!fin(hr, alpha, power)) err.push("Enter HR, alpha and power.");
    else {
      if (!(hr > 0) || hr === 1) err.push("The hazard ratio must be positive and not equal to 1.");
      if (!(alpha > 0 && alpha < 1)) err.push("Alpha must be between 0 and 1.");
      if (!(power > 0 && power < 1)) err.push("Power must be between 0 and 1.");
      if (!(p > 0 && p < 1)) err.push("Allocation proportion must be strictly between 0 and 1.");
      if (ep !== null && !(ep > 0 && ep <= 1)) err.push("Event probability must be above 0 and at most 1 (or leave blank).");
    }
    if (err.length) return { errors: err };
    const za = zq(1 - alpha / 2), zb = zq(power);
    const lh = Math.log(hr);
    const eventsRaw = ((za + zb) * (za + zb)) / (p * (1 - p) * lh * lh);
    const out = { events_raw: eventsRaw, events: Math.ceil(eventsRaw), za, zb };
    if (ep !== null) {
      out.total_n = Math.ceil(out.events / ep);
      out.n_group1 = Math.ceil(out.total_n * p);
      out.n_group2 = out.total_n - out.n_group1;
    }
    return out;
  }

  const COMPUTE = {
    ttestInd, ttestPaired, twoByTwo, mcnemar, oneProportion, diagnostic, ssMeans, ssProps, ssSurvival,
    fisherP, binomTwoSided, clopperPearson, wilson, powerT,
  };
  root.CALC_COMPUTE = COMPUTE;

  /* =====================================================================
     Formatting helpers
     ===================================================================== */
  const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  function num(x, sig) {
    if (x === null || x === undefined || Number.isNaN(x)) return "n/a";
    if (x === Infinity) return "∞";
    if (x === -Infinity) return "−∞";
    sig = sig || 4;
    if (x === 0) return "0";
    const ax = Math.abs(x);
    if (ax >= 1e7 || ax < 1e-4) return x.toExponential(sig - 1).replace("-", "−");
    return Number(x.toPrecision(sig)).toString().replace("-", "−");
  }
  function pFmt(p) {
    if (p === null || p === undefined || Number.isNaN(p)) return { main: "n/a", exact: "" };
    const main = p < 0.001 ? "p < 0.001" : "p = " + Number(p.toPrecision(3)).toString();
    return { main, exact: "exact p = " + p.toExponential(4) };
  }
  const pct = (x) => num(x * 100, 4) + "%";
  const fmtCI = (lo, hi, f) => "[" + (f || num)(lo) + ", " + (f || num)(hi) + "]";
  const rnum = (x) => (Number.isFinite(x) ? String(Number(x.toPrecision(10))) : String(x));

  /* =====================================================================
     Calculator specs (DOM side)
     ===================================================================== */
  // Each spec: fields, optional selects, run(v) -> view model or {errors}.
  // View model: { key: [{label, value, small}], rows: [[label, value, small]], method, r: [lines], notes: [] }
  const SPECS = [];

  SPECS.push({
    id: "ttest-ind",
    title: "Two-sample t-test (summary data)",
    summary: "Compare two group means from mean, SD and n. Welch or Student, with 95% CI and Cohen's d.",
    related: ["independent-t"],
    selects: [{ k: "welch", label: "Variance assumption", def: "1", options: [["1", "Welch (unequal variances, default)"], ["0", "Student (equal variances)"]] }],
    fields: [
      { k: "m1", label: "Group 1 mean", def: 52.1, step: "any", group: "Group 1" },
      { k: "s1", label: "Group 1 SD", def: 9.8, step: "any", min: 0, group: "Group 1" },
      { k: "n1", label: "Group 1 n", def: 30, step: 1, min: 2, group: "Group 1" },
      { k: "m2", label: "Group 2 mean", def: 47.3, step: "any", group: "Group 2" },
      { k: "s2", label: "Group 2 SD", def: 11.2, step: "any", min: 0, group: "Group 2" },
      { k: "n2", label: "Group 2 n", def: 34, step: 1, min: 2, group: "Group 2" },
    ],
    run(v) {
      const r = ttestInd(Object.assign({}, v, { welch: v.welch === "1" }));
      if (r.errors) return r;
      const pf = pFmt(r.p);
      return {
        key: [
          { label: "Difference (group 1 − group 2)", value: num(r.diff) },
          { label: "95% CI", value: fmtCI(r.ci_lo, r.ci_hi) },
          { label: "p-value (two-sided)", value: pf.main, small: pf.exact },
        ],
        rows: [["t statistic", num(r.t)], ["Degrees of freedom", num(r.df)], ["Standard error", num(r.se)],
          ["Cohen's d (pooled SD)", num(r.d)], ["Pooled SD", num(r.sp)]],
        method: r.welch ? "Welch two-sample t-test (unequal variances)" : "Student two-sample t-test (pooled variance)",
        r: [
          "t.test(x1, x2, var.equal = " + (r.welch ? "FALSE" : "TRUE") + ")   # with raw data",
          "2 * pt(-abs(" + rnum(r.t) + "), df = " + rnum(r.df) + ")   # p-value from the summary numbers",
        ],
      };
    },
  });

  SPECS.push({
    id: "ttest-paired",
    title: "Paired / one-sample t-test (summary data)",
    summary: "Test a mean difference (or a single mean) against a hypothesised value, from mean, SD and n.",
    related: ["paired-t", "one-sample-t"],
    fields: [
      { k: "mean", label: "Mean difference (or mean)", def: 2.4, step: "any" },
      { k: "sd", label: "SD of differences (or SD)", def: 5.1, step: "any", min: 0 },
      { k: "n", label: "n (pairs or subjects)", def: 24, step: 1, min: 2 },
      { k: "mu0", label: "Hypothesised mean μ₀", def: 0, step: "any", optional: true },
    ],
    run(v) {
      const r = ttestPaired(v);
      if (r.errors) return r;
      const pf = pFmt(r.p);
      return {
        key: [
          { label: "Mean − μ₀", value: num(r.mean - r.mu0) },
          { label: "95% CI for the mean", value: fmtCI(r.ci_lo, r.ci_hi) },
          { label: "p-value (two-sided)", value: pf.main, small: pf.exact },
        ],
        rows: [["t statistic", num(r.t)], ["Degrees of freedom", num(r.df)], ["Standard error", num(r.se)],
          ["Standardised effect (d_z)", num(r.dz)]],
        method: "One-sample t-test on the differences (paired t-test)",
        r: [
          "t.test(d, mu = " + rnum(r.mu0) + ")   # d = vector of differences (or t.test(x, y, paired = TRUE))",
          "2 * pt(-abs(" + rnum(r.t) + "), df = " + rnum(r.df) + ")   # p-value from the summary numbers",
        ],
      };
    },
  });

  SPECS.push({
    id: "twobytwo",
    title: "2×2 table: chi-square, Fisher, RR and OR",
    summary: "Chi-square (with and without Yates), Fisher exact, risk difference, relative risk and odds ratio with CIs.",
    related: ["chi-square-independence", "fisher-exact", "two-proportion-z"],
    layout: "table2x2",
    fields: [
      { k: "a", label: "a: Group 1, outcome yes", def: 18, step: 1, min: 0 },
      { k: "b", label: "b: Group 1, outcome no", def: 32, step: 1, min: 0 },
      { k: "c", label: "c: Group 2, outcome yes", def: 8, step: 1, min: 0 },
      { k: "d", label: "d: Group 2, outcome no", def: 42, step: 1, min: 0 },
    ],
    run(v) {
      const r = twoByTwo(v);
      if (r.errors) return r;
      const pf = pFmt(r.fisher_p);
      const rows = [];
      if (r.chi !== undefined) {
        const p1 = pFmt(r.chi_p), p2 = pFmt(r.yates_p);
        rows.push(["Pearson χ² (no correction)", num(r.chi) + ", df 1, " + p1.main, p1.exact]);
        rows.push(["Pearson χ² (Yates)", num(r.yates) + ", df 1, " + p2.main, p2.exact]);
        rows.push(["Smallest expected count", num(r.minE) + (r.minE < 5 ? "  (below 5: prefer Fisher)" : "")]);
      } else rows.push(["Chi-square", r.chi_note]);
      if (r.rd !== undefined) {
        rows.push(["Risk group 1 / group 2", pct(r.p1) + " / " + pct(r.p2)]);
        rows.push(["Risk difference (Wald 95% CI)", num(r.rd) + " " + fmtCI(r.rd_lo, r.rd_hi)]);
      }
      rows.push(["Relative risk (95% CI)", num(r.rr) + " " + fmtCI(r.rr_lo, r.rr_hi)]);
      rows.push(["Odds ratio (95% CI)", num(r.or) + " " + fmtCI(r.or_lo, r.or_hi)]);
      const m = [`matrix(c(${v.a}, ${v.c}, ${v.b}, ${v.d}), nrow = 2)`];
      return {
        key: [
          { label: "Fisher exact p (two-sided)", value: pf.main, small: pf.exact },
          { label: "Odds ratio", value: num(r.or), small: "95% CI " + fmtCI(r.or_lo, r.or_hi) },
          { label: "Relative risk", value: num(r.rr), small: "95% CI " + fmtCI(r.rr_lo, r.rr_hi) },
        ],
        rows,
        method: "Pearson chi-square, Fisher exact test, and log-Wald intervals for RR and OR",
        notes: r.adjusted ? ["A cell is zero, so 0.5 was added to every cell for the RR and OR (and their CIs). Fisher and chi-square use the raw counts."] : [],
        r: [
          "tab <- " + m[0],
          "chisq.test(tab, correct = FALSE); chisq.test(tab, correct = TRUE)",
          "fisher.test(tab)   # two-sided p; its OR is the conditional MLE, not a*d/(b*c)",
        ],
      };
    },
  });

  SPECS.push({
    id: "mcnemar",
    title: "McNemar test (paired binary data)",
    summary: "Exact and chi-square McNemar tests from the two discordant counts.",
    related: ["mcnemar"],
    fields: [
      { k: "b", label: "b: pairs positive then negative", def: 12, step: 1, min: 0 },
      { k: "c", label: "c: pairs negative then positive", def: 4, step: 1, min: 0 },
    ],
    run(v) {
      const r = mcnemar(v);
      if (r.errors) return r;
      const pe = pFmt(r.exact_p);
      const rows = [["Discordant pairs (b + c)", String(r.n)]];
      if (r.chi !== undefined) {
        const p1 = pFmt(r.chi_p), p2 = pFmt(r.chi_cc_p);
        rows.push(["χ² (no correction)", num(r.chi) + ", df 1, " + p1.main, p1.exact]);
        rows.push(["χ² (continuity corrected)", num(r.chi_cc) + ", df 1, " + p2.main, p2.exact]);
      }
      return {
        key: [{ label: "Exact McNemar p (two-sided)", value: pe.main, small: pe.exact }],
        rows,
        method: "McNemar's test: exact binomial on the discordant pairs, and chi-square",
        notes: r.note ? [r.note] : [],
        r: [
          "binom.test(" + v.b + ", " + (v.b + v.c) + ", p = 0.5)   # exact",
          "mcnemar.test(matrix(c(0, " + v.c + ", " + v.b + ", 0), 2), correct = TRUE)   # chi-square (use any concordant counts)",
        ],
      };
    },
  });

  SPECS.push({
    id: "one-proportion",
    title: "One-proportion exact test and CIs",
    summary: "Exact binomial test against p₀, with Clopper–Pearson and Wilson 95% intervals.",
    related: ["one-proportion"],
    fields: [
      { k: "x", label: "Successes x", def: 14, step: 1, min: 0 },
      { k: "n", label: "Trials n", def: 40, step: 1, min: 1 },
      { k: "p0", label: "Null proportion p₀", def: 0.25, step: "any", min: 0, max: 1 },
    ],
    run(v) {
      const r = oneProportion(v);
      if (r.errors) return r;
      const pf = pFmt(r.p_exact);
      return {
        key: [
          { label: "Observed proportion", value: pct(r.phat) },
          { label: "Exact p (two-sided)", value: pf.main, small: pf.exact },
          { label: "Clopper–Pearson 95% CI", value: fmtCI(r.cp_lo, r.cp_hi, (x) => pct(x)) },
        ],
        rows: [["Wilson 95% CI", fmtCI(r.wilson_lo, r.wilson_hi, (x) => pct(x))],
          ["Clopper–Pearson (decimal)", fmtCI(r.cp_lo, r.cp_hi)], ["Wilson (decimal)", fmtCI(r.wilson_lo, r.wilson_hi)]],
        method: "Exact binomial test (R binom.test), Clopper–Pearson and Wilson score intervals",
        r: [
          "binom.test(" + v.x + ", " + v.n + ", p = " + rnum(v.p0) + ")   # exact p and Clopper–Pearson CI",
          "binom::binom.confint(" + v.x + ", " + v.n + ", methods = \"wilson\")   # Wilson CI",
        ],
      };
    },
  });

  SPECS.push({
    id: "diagnostic",
    title: "Diagnostic test accuracy",
    summary: "Sensitivity, specificity, predictive values and likelihood ratios, with an optional prevalence adjustment.",
    related: ["roc-auc"],
    layout: "table2x2",
    fields: [
      { k: "tp", label: "TP: test +, disease +", def: 90, step: 1, min: 0 },
      { k: "fp", label: "FP: test +, disease −", def: 30, step: 1, min: 0 },
      { k: "fn", label: "FN: test −, disease +", def: 10, step: 1, min: 0 },
      { k: "tn", label: "TN: test −, disease −", def: 170, step: 1, min: 0 },
      { k: "prev", label: "Prevalence % (optional)", def: "", step: "any", min: 0, max: 100, optional: true, wide: true, scale: 0.01 },
    ],
    run(v) {
      const r = diagnostic(v);
      if (r.errors) return r;
      const show = (o) => (o ? pct(o.est) + " " + fmtCI(o.lo, o.hi, (x) => pct(x)) : "n/a");
      const rows = [
        ["Sensitivity", show(r.sens)], ["Specificity", show(r.spec)], ["PPV (this sample)", show(r.ppv)],
        ["NPV (this sample)", show(r.npv)], ["Accuracy", show(r.acc)],
        ["LR+", num(r.lr_pos)], ["LR−", num(r.lr_neg)],
      ];
      if (r.ppv_prev !== undefined) {
        rows.push(["PPV at stated prevalence", r.ppv_prev === null ? "n/a" : pct(r.ppv_prev)]);
        rows.push(["NPV at stated prevalence", r.npv_prev === null ? "n/a" : pct(r.npv_prev)]);
      }
      return {
        key: [
          { label: "Sensitivity", value: r.sens ? pct(r.sens.est) : "n/a", small: r.sens ? "95% CI " + fmtCI(r.sens.lo, r.sens.hi, pct) : "" },
          { label: "Specificity", value: r.spec ? pct(r.spec.est) : "n/a", small: r.spec ? "95% CI " + fmtCI(r.spec.lo, r.spec.hi, pct) : "" },
          { label: "LR+ / LR−", value: num(r.lr_pos) + " / " + num(r.lr_neg) },
        ],
        rows,
        method: "Diagnostic accuracy measures with Clopper–Pearson 95% CIs; predictive values by Bayes' rule when prevalence is given",
        notes: r.ppv_prev === undefined ? ["Predictive values above reflect the prevalence in this sample. Enter a prevalence to recompute them for another population."] : [],
        r: [
          "binom.test(" + v.tp + ", " + (v.tp + v.fn) + ")$conf.int   # sensitivity CI",
          "binom.test(" + v.tn + ", " + (v.tn + v.fp) + ")$conf.int   # specificity CI",
        ],
      };
    },
  });

  const SS_COMMON = [
    { k: "alpha", label: "Significance level α (two-sided)", def: 0.05, step: "any", min: 0, max: 1 },
    { k: "power", label: "Power (1 − β)", def: 0.8, step: "any", min: 0, max: 1 },
    { k: "dropout", label: "Expected dropout % (optional)", def: "", step: "any", min: 0, max: 100, optional: true },
  ];

  SPECS.push({
    id: "ss-means",
    title: "Sample size: two means",
    summary: "Participants per group for a two-sample t-test, matching R power.t.test.",
    related: ["sample-size-means-proportions"],
    fields: [
      { k: "delta", label: "Difference to detect δ", def: 5, step: "any" },
      { k: "sd", label: "Standard deviation σ", def: 10, step: "any", min: 0 },
    ].concat(SS_COMMON),
    run(v) {
      const r = ssMeans(v);
      if (r.errors) return r;
      const key = [{ label: "n per group", value: String(r.n), small: "total " + r.total }];
      if (r.n_dropout) key.push({ label: "With " + num(v.dropout) + "% dropout", value: String(r.n_dropout) + " per group", small: "total " + r.total_dropout });
      return {
        key,
        rows: [["Unrounded n per group", num(r.n_raw, 6)], ["Power at rounded n", pct(r.power_achieved)]],
        method: "Two-sample t-test power (noncentral t), two-sided, equal group sizes",
        notes: r.n_dropout ? ["Dropout inflation: ceiling of n ÷ (1 − dropout), applied to the rounded-up n."] : [],
        r: ["power.t.test(delta = " + rnum(v.delta) + ", sd = " + rnum(v.sd) + ", sig.level = " + rnum(v.alpha) + ", power = " + rnum(v.power) +
          ", type = \"two.sample\", alternative = \"two.sided\")"],
      };
    },
  });

  SPECS.push({
    id: "ss-props",
    title: "Sample size: two proportions",
    summary: "Participants per group to compare two proportions, matching R power.prop.test.",
    related: ["sample-size-means-proportions"],
    fields: [
      { k: "p1", label: "Proportion in group 1, p₁", def: 0.3, step: "any", min: 0, max: 1 },
      { k: "p2", label: "Proportion in group 2, p₂", def: 0.45, step: "any", min: 0, max: 1 },
    ].concat(SS_COMMON),
    run(v) {
      const r = ssProps(v);
      if (r.errors) return r;
      const key = [{ label: "n per group", value: String(r.n), small: "total " + r.total }];
      if (r.n_dropout) key.push({ label: "With " + num(v.dropout) + "% dropout", value: String(r.n_dropout) + " per group", small: "total " + r.total_dropout });
      return {
        key,
        rows: [["Unrounded n per group", num(r.n_raw, 6)]],
        method: "Normal-approximation test for two proportions (as in power.prop.test), two-sided",
        notes: r.n_dropout ? ["Dropout inflation: ceiling of n ÷ (1 − dropout), applied to the rounded-up n."] : [],
        r: ["power.prop.test(p1 = " + rnum(v.p1) + ", p2 = " + rnum(v.p2) + ", sig.level = " + rnum(v.alpha) + ", power = " + rnum(v.power) + ")"],
      };
    },
  });

  SPECS.push({
    id: "ss-survival",
    title: "Sample size: survival (events needed)",
    summary: "Schoenfeld formula for the number of events in a log-rank / Cox comparison, and optionally total N.",
    related: ["sample-size-survival"],
    fields: [
      { k: "hr", label: "Hazard ratio to detect", def: 0.7, step: "any", min: 0 },
      { k: "alpha", label: "Significance level α (two-sided)", def: 0.05, step: "any", min: 0, max: 1 },
      { k: "power", label: "Power (1 − β)", def: 0.8, step: "any", min: 0, max: 1 },
      { k: "p", label: "Proportion allocated to group 1", def: 0.5, step: "any", min: 0, max: 1, optional: true },
      { k: "eventProb", label: "Overall event probability (optional)", def: "", step: "any", min: 0, max: 1, optional: true },
    ],
    run(v) {
      const r = ssSurvival(v);
      if (r.errors) return r;
      const p = Number.isFinite(v.p) ? v.p : 0.5;
      const key = [{ label: "Events required", value: String(r.events) }];
      if (r.total_n) key.push({ label: "Total N", value: String(r.total_n), small: r.n_group1 + " + " + r.n_group2 });
      return {
        key,
        rows: [["Unrounded events", num(r.events_raw, 6)], ["z(1 − α/2), z(power)", num(r.za, 5) + ", " + num(r.zb, 5)]],
        method: "Schoenfeld (1983) required number of events",
        notes: r.total_n ? [] : ["Enter the overall probability of observing an event during follow-up to convert events into a total sample size."],
        r: [
          "d <- (qnorm(1 - " + rnum(v.alpha / 2) + ") + qnorm(" + rnum(v.power) + "))^2 / (" + rnum(p) + " * (1 - " + rnum(p) + ") * log(" + rnum(v.hr) + ")^2)",
          "ceiling(d)" + (r.total_n ? "; ceiling(ceiling(d) / " + rnum(v.eventProb) + ")   # total N" : ""),
        ],
      };
    },
  });

  /* =====================================================================
     DOM
     ===================================================================== */
  function parseField(raw, f) {
    const t = String(raw).trim();
    if (t === "") return f.optional ? NaN : NaN;
    return Number(t);
  }

  function validateFields(spec, vals) {
    const errs = [];
    for (const f of spec.fields) {
      const x = vals[f.k];
      if (Number.isNaN(x)) {
        if (!f.optional) errs.push(`Enter a value for "${f.label}".`);
        continue;
      }
      if (f.min !== undefined && x < f.min) errs.push(`"${f.label}" must be at least ${f.min}.`);
      if (f.max !== undefined && x > f.max) errs.push(`"${f.label}" must be at most ${f.max}.`);
    }
    return errs;
  }

  function mountSpec(spec, el, preset) {
    preset = preset || {};
    const init = (k, def) => (preset[k] !== undefined ? preset[k] : def);
    if (typeof document === "undefined") return;
    const uid = "calc-" + spec.id + "-" + Math.random().toString(36).slice(2, 7);
    const selects = (spec.selects || []).map((s) => `
      <div class="calc-field calc-field-wide"><label for="${uid}-${s.k}">${esc(s.label)}</label>
        <select id="${uid}-${s.k}" data-k="${s.k}">${s.options.map((o) => `<option value="${o[0]}"${o[0] === String(init(s.k, s.def)) ? " selected" : ""}>${esc(o[1])}</option>`).join("")}</select></div>`).join("");
    let lastGroup = null;
    const fieldHtml = spec.fields.map((f) => {
      let head = "";
      if (f.group && f.group !== lastGroup) head = `<div class="calc-group-title">${esc(f.group)}</div>`;
      lastGroup = f.group || lastGroup;
      return head + `<div class="calc-field${f.wide ? " calc-field-wide" : ""}"><label for="${uid}-${f.k}">${esc(f.label)}</label>
        <input id="${uid}-${f.k}" data-k="${f.k}" type="number" inputmode="decimal" step="${f.step || "any"}" value="${init(f.k, f.def)}"${f.min !== undefined ? ` min="${f.min}"` : ""}${f.max !== undefined ? ` max="${f.max}"` : ""}></div>`;
    }).join("");
    el.innerHTML = `<div class="calc">
      <div class="calc-inputs"><div class="calc-grid${spec.layout === "table2x2" ? " calc-grid-2" : ""}">${selects}${fieldHtml}</div></div>
      <div class="calc-results" aria-live="polite"></div></div>`;
    const out = el.querySelector(".calc-results");

    function update() {
      const vals = {};
      el.querySelectorAll("select[data-k]").forEach((s) => { vals[s.dataset.k] = s.value; });
      spec.fields.forEach((f) => {
        const x = parseField(el.querySelector(`input[data-k="${f.k}"]`).value, f);
        vals[f.k] = f.scale && Number.isFinite(x) ? x * f.scale : x;
      });
      // Range checks against the displayed (unscaled) values
      const checkVals = {};
      spec.fields.forEach((f) => { checkVals[f.k] = f.scale && Number.isFinite(vals[f.k]) ? vals[f.k] / f.scale : vals[f.k]; });
      let errs = validateFields(spec, checkVals);
      let res = null;
      if (!errs.length) {
        try { res = spec.run(vals); } catch (e) { errs = ["Could not compute: " + e.message]; }
        if (res && res.errors) errs = res.errors;
      }
      if (errs.length) {
        out.innerHTML = `<div class="calc-error" role="alert"><strong>Check your inputs</strong><ul>${errs.map((e) => `<li>${esc(e)}</li>`).join("")}</ul></div>`;
        return;
      }
      out.innerHTML = `
        <div class="calc-key">${res.key.map((k) => `<div class="calc-stat"><div class="calc-stat-label">${esc(k.label)}</div>
          <div class="calc-stat-value">${esc(k.value)}</div>${k.small ? `<div class="calc-small">${esc(k.small)}</div>` : ""}</div>`).join("")}</div>
        ${res.rows && res.rows.length ? `<table class="calc-table"><tbody>${res.rows.map((r) => `<tr><th scope="row">${esc(r[0])}</th><td>${esc(r[1])}${r[2] ? `<div class="calc-small">${esc(r[2])}</div>` : ""}</td></tr>`).join("")}</tbody></table>` : ""}
        ${(res.notes || []).map((n) => `<p class="calc-note">${esc(n)}</p>`).join("")}
        <div class="calc-method"><span class="calc-method-label">Method</span> ${esc(res.method)}</div>
        <div class="calc-r"><div class="calc-r-label">Equivalent R</div><pre class="calc-code"><code>${res.r.map(esc).join("\n")}</code></pre></div>`;
    }
    el.addEventListener("input", update);
    el.addEventListener("change", update);
    update();
  }

  root.CALCULATORS = SPECS.map((s) => ({
    id: s.id, title: s.title, summary: s.summary, related: s.related,
    mount(el, preset) { mountSpec(s, el, preset); },
  }));
})(typeof window !== "undefined" ? window : globalThis);
