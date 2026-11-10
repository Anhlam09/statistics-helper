/* StatMath: small dependency-free numerical library for the calculators. */
(function (root) {
  "use strict";
  const LN_SQRT_PI = 0.5723649429247001; // log(sqrt(pi))
  const SQRT_2DPI = 0.7978845608028654; // sqrt(2/pi)

  /* ---------- log-gamma (shift + Stirling series) ---------- */
  function lgamma(x) {
    if (!(x > 0)) return x === 0 ? Infinity : NaN;
    let shift = 0;
    let prod = 1;
    while (x < 15) { prod *= x; x += 1; shift++; if (prod > 1e250) { shift = -1; break; } }
    // recompute safely if product would overflow (never for x>0 start, <=15 factors)
    const x2 = x * x;
    const series = 1 / (12 * x) - 1 / (360 * x * x2) + 1 / (1260 * x2 * x2 * x) -
      1 / (1680 * x2 * x2 * x2 * x) + 1 / (1188 * x2 * x2 * x2 * x2 * x);
    const st = (x - 0.5) * Math.log(x) - x + 0.9189385332046727 + series;
    return st - Math.log(prod);
  }

  /* ---------- incomplete beta ---------- */
  function betacf(a, b, x) {
    const FPMIN = 1e-300;
    const qab = a + b, qap = a + 1, qam = a - 1;
    let c = 1, d = 1 - (qab * x) / qap;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    d = 1 / d;
    let h = d;
    for (let m = 1; m <= 100000; m++) {
      const m2 = 2 * m;
      let aa = (m * (b - m) * x) / ((qam + m2) * (a + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
      c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
      d = 1 / d; h *= d * c;
      aa = (-(a + m) * (qab + m) * x) / ((a + m2) * (qap + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
      c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
      d = 1 / d;
      const del = d * c; h *= del;
      if (Math.abs(del - 1) < 1e-16) break;
    }
    return h;
  }
  /** Regularized incomplete beta I_x(a, b). */
  function betainc(x, a, b) {
    if (x <= 0) return 0;
    if (x >= 1) return 1;
    const bt = Math.exp(lgamma(a + b) - lgamma(a) - lgamma(b) + a * Math.log(x) + b * Math.log1p(-x));
    if (x < (a + 1) / (a + b + 2)) return (bt * betacf(a, b, x)) / a;
    return 1 - (bt * betacf(b, a, 1 - x)) / b;
  }
  /** Upper-tail version 1 - I_x(a,b), computed without cancellation. */
  function betaincUpper(x, a, b) { return betainc(1 - x, b, a); }

  /** Quantile of Beta(a, b) by bisection. */
  function betaQuantile(p, a, b) {
    if (p <= 0) return 0;
    if (p >= 1) return 1;
    let lo = 0, hi = 1;
    for (let i = 0; i < 2000; i++) {
      const mid = (lo + hi) / 2;
      if (mid === lo || mid === hi) break;
      if (betainc(mid, a, b) < p) lo = mid; else hi = mid;
    }
    return (lo + hi) / 2;
  }

  /* ---------- incomplete gamma ---------- */
  function gammaSeries(a, x) {
    let ap = a, sum = 1 / a, del = sum;
    for (let n = 0; n < 1000000; n++) {
      ap += 1; del *= x / ap; sum += del;
      if (Math.abs(del) < Math.abs(sum) * 1e-17) break;
    }
    return sum * Math.exp(-x + a * Math.log(x) - lgamma(a));
  }
  function gammaCF(a, x) {
    const FPMIN = 1e-300;
    let b = x + 1 - a, c = 1 / FPMIN, d = 1 / b, h = d;
    for (let i = 1; i < 1000000; i++) {
      const an = -i * (i - a);
      b += 2;
      d = an * d + b; if (Math.abs(d) < FPMIN) d = FPMIN;
      c = b + an / c; if (Math.abs(c) < FPMIN) c = FPMIN;
      d = 1 / d;
      const del = d * c; h *= del;
      if (Math.abs(del - 1) < 1e-16) break;
    }
    return Math.exp(-x + a * Math.log(x) - lgamma(a)) * h;
  }
  /** Regularized lower incomplete gamma P(a, x). */
  function gammaP(a, x) {
    if (x <= 0) return 0;
    return x < a + 1 ? gammaSeries(a, x) : 1 - gammaCF(a, x);
  }
  /** Regularized upper incomplete gamma Q(a, x). */
  function gammaQ(a, x) {
    if (x <= 0) return 1;
    return x < a + 1 ? 1 - gammaSeries(a, x) : gammaCF(a, x);
  }

  /* ---------- normal ---------- */
  function normCdf(x) {
    if (x === 0) return 0.5;
    const q = 0.5 * gammaQ(0.5, (x * x) / 2);
    return x < 0 ? q : 1 - q;
  }
  function normSf(x) { return normCdf(-x); }
  /** Wichura AS241 (PPND16). */
  function normQuantile(p) {
    if (!(p > 0 && p < 1)) return p === 0 ? -Infinity : p === 1 ? Infinity : NaN;
    const q = p - 0.5;
    let r, val;
    if (Math.abs(q) <= 0.425) {
      r = 0.180625 - q * q;
      return q * (((((((r * 2509.0809287301226727 + 33430.575583588128105) * r + 67265.770927008700853) * r +
        45921.953931549871457) * r + 13731.693765509461125) * r + 1971.5909503065514427) * r +
        133.14166789178437745) * r + 3.387132872796366608) /
        (((((((r * 5226.495278852545925 + 28729.085735721942674) * r + 39307.89580009271061) * r +
        21213.794301586595867) * r + 5394.1960214247511077) * r + 687.1870074920579083) * r +
        42.313330701600911252) * r + 1);
    }
    r = q < 0 ? p : 1 - p;
    r = Math.sqrt(-Math.log(r));
    if (r <= 5) {
      r -= 1.6;
      val = (((((((r * 7.7454501427834140764e-4 + 0.0227238449892691845833) * r + 0.24178072517745061177) * r +
        1.27045825245236838258) * r + 3.64784832476320460504) * r + 5.7694972214606914055) * r +
        4.6303378461565452959) * r + 1.42343711074968357734) /
        (((((((r * 1.05075007164441684324e-9 + 5.475938084995344946e-4) * r + 0.0151986665636164571966) * r +
        0.14810397642748007459) * r + 0.68976733498510000455) * r + 1.6763848301838038494) * r +
        2.05319162663775882187) * r + 1);
    } else {
      r -= 5;
      val = (((((((r * 2.01033439929228813265e-7 + 2.71155556874348757815e-5) * r + 0.0012426609473880784386) * r +
        0.026532189526576123093) * r + 0.29656057182850489123) * r + 1.7848265399172913358) * r +
        5.4637849111641143699) * r + 6.6579046435011037772) /
        (((((((r * 2.04426310338993978564e-15 + 1.4215117583164458887e-7) * r + 1.8463183175100546818e-5) * r +
        7.868691311456132591e-4) * r + 0.0148753612908506148525) * r + 0.13692988092273580531) * r +
        0.59983220655588793769) * r + 1);
    }
    return q < 0 ? -val : val;
  }

  /* ---------- Student t ---------- */
  /** Two-sided tail P(|T| > |t|). */
  function tTwoTail(t, df) {
    if (!isFinite(t)) return 0;
    return betainc(df / (df + t * t), df / 2, 0.5);
  }
  function tSf(t, df) {
    const half = 0.5 * tTwoTail(t, df);
    return t > 0 ? half : 1 - half;
  }
  function tCdf(t, df) {
    const half = 0.5 * tTwoTail(t, df);
    return t > 0 ? 1 - half : half;
  }
  function tQuantile(p, df) {
    if (!(p > 0 && p < 1)) return p === 0 ? -Infinity : p === 1 ? Infinity : NaN;
    if (p === 0.5) return 0;
    const q = Math.min(p, 1 - p); // upper tail prob
    const target = 2 * q;
    let lo = 0, hi = 1;
    while (tTwoTail(hi, df) > target && hi < 1e300) { lo = hi; hi *= 2; }
    for (let i = 0; i < 300; i++) {
      const mid = (lo + hi) / 2;
      if (mid === lo || mid === hi) break;
      if (tTwoTail(mid, df) > target) lo = mid; else hi = mid;
    }
    const x = (lo + hi) / 2;
    return p < 0.5 ? -x : x;
  }

  /* ---------- chi-square ---------- */
  function chisqCdf(x, df) { return gammaP(df / 2, x / 2); }
  function chisqSf(x, df) { return gammaQ(df / 2, x / 2); }

  /* ---------- binomial / hypergeometric ---------- */
  function lchoose(n, k) { return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1); }
  function binomPmf(x, n, p) {
    if (x < 0 || x > n || x !== Math.floor(x)) return 0;
    if (p === 0) return x === 0 ? 1 : 0;
    if (p === 1) return x === n ? 1 : 0;
    return Math.exp(lchoose(n, x) + x * Math.log(p) + (n - x) * Math.log1p(-p));
  }
  /** P(X <= x) */
  function binomCdf(x, n, p) {
    if (x < 0) return 0;
    if (x >= n) return 1;
    x = Math.floor(x);
    if (p <= 0) return 1;
    if (p >= 1) return 0;
    return betainc(1 - p, n - x, x + 1);
  }
  /** P(X > x) */
  function binomSf(x, n, p) {
    if (x < 0) return 1;
    if (x >= n) return 0;
    x = Math.floor(x);
    if (p <= 0) return 0;
    if (p >= 1) return 1;
    return betainc(p, x + 1, n - x);
  }
  /** P(X = x) with K successes among N, n draws. */
  function hyperPmf(x, K, N, n) {
    if (x < Math.max(0, n - (N - K)) || x > Math.min(n, K)) return 0;
    return Math.exp(lchoose(K, x) + lchoose(N - K, n - x) - lchoose(N, n));
  }

  /* ---------- noncentral t (AS 243, as in R's pnt) ---------- */
  /** Lower-tail P(T <= t) for noncentral t with df and ncp. */
  function ntCdf(t, df, ncp) {
    if (df <= 0) return NaN;
    if (ncp === 0) return tCdf(t, df);
    if (!isFinite(t)) return t < 0 ? 0 : 1;
    let negdel, tt, del;
    if (t >= 0) { negdel = false; tt = t; del = ncp; }
    else {
      if (ncp > 40) return 0;
      negdel = true; tt = -t; del = -ncp;
    }
    if (df > 4e5 || del * del > 2 * Math.LN2 * 1021) {
      const s = 1 / (4 * df);
      const z = (tt * (1 - s) - del) / Math.sqrt(1 + tt * tt * 2 * s);
      const lower = normCdf(z);
      return negdel ? 1 - lower : lower;
    }
    let x = t * t;
    const rxb0 = df / (x + df);
    x = x / (x + df);
    let tnc;
    if (x > 0) {
      const lambda = del * del;
      let p = 0.5 * Math.exp(-0.5 * lambda);
      if (p === 0) return negdel ? 1 : 0;
      let q = SQRT_2DPI * p * del;
      let s = 0.5 - p;
      if (s < 1e-7) s = -0.5 * Math.expm1(-0.5 * lambda);
      let a = 0.5;
      const b = 0.5 * df;
      const rxb = Math.pow(rxb0, b);
      const albeta = LN_SQRT_PI + lgamma(b) - lgamma(0.5 + b);
      let xodd = betainc(x, a, b);
      let godd = 2 * rxb * Math.exp(a * Math.log(x) - albeta);
      tnc = b * x;
      let xeven = tnc < 2.220446049250313e-16 ? tnc : 1 - rxb;
      let geven = tnc * rxb;
      tnc = p * xodd + q * xeven;
      for (let it = 1; it <= 1000; it++) {
        a += 1;
        xodd -= godd;
        xeven -= geven;
        godd *= (x * (a + b - 1)) / a;
        geven *= (x * (a + b - 0.5)) / (a + 0.5);
        p *= lambda / (2 * it);
        q *= lambda / (2 * it + 1);
        tnc += p * xodd + q * xeven;
        s -= p;
        if (s < -1e-10) break;
        if (s >= 0) {
          const errbd = 2 * s * (xodd - godd);
          if (Math.abs(errbd) < 1e-12) break;
        }
      }
    } else tnc = 0;
    tnc += normCdf(-del);
    tnc = Math.min(tnc, 1);
    return negdel ? 1 - tnc : tnc;
  }

  const API = {
    lgamma, lchoose, betainc, betaincUpper, betaQuantile, gammaP, gammaQ,
    normCdf, normSf, normQuantile, tCdf, tSf, tTwoTail, tQuantile,
    chisqCdf, chisqSf, binomPmf, binomCdf, binomSf, hyperPmf, ntCdf,
  };
  root.StatMath = API;
  if (typeof module !== "undefined" && module.exports) module.exports = API;
})(typeof window !== "undefined" ? window : globalThis);
