/*
 * Teaching diagrams rendered from markdown placeholders:
 *   <div class="viz" data-viz="t" data-stat="-7.34" data-df="9"></div>
 *   <div class="viz" data-viz="z" data-stat="1.8"></div>
 *   <div class="viz" data-viz="chisq" data-stat="6.48" data-df="1"></div>
 *   <div class="viz" data-viz="f" data-stat="5.36" data-df1="2" data-df2="21"></div>
 *   <div class="viz" data-viz="dots" data-groups='{"Before":[152,148],"After":[144,145]}' data-unit="mmHg"></div>
 * Null-distribution plots shade the p-value area beyond the observed statistic
 * (two-sided for t and z, upper tail for chi-square and F).
 */
(function () {
  "use strict";
  const SM = () => window.StatMath;
  const NS = "http://www.w3.org/2000/svg";
  const W = 640, H = 240, PAD = { l: 36, r: 16, t: 18, b: 40 };

  const fmt = (x, d = 2) => (Math.abs(x) >= 100 ? x.toFixed(0) : x.toFixed(d));
  const fmtP = (p) => (p < 0.001 ? "p < 0.001" : `p = ${p < 0.01 ? p.toFixed(4) : p.toFixed(3)}`);

  // log-density functions (via lgamma) so large df stay finite
  function densities() {
    const lg = SM().lgamma;
    return {
      z: (x) => Math.exp(-0.5 * x * x) / Math.sqrt(2 * Math.PI),
      t: (x, df) => Math.exp(lg((df + 1) / 2) - lg(df / 2) - 0.5 * Math.log(df * Math.PI) - ((df + 1) / 2) * Math.log(1 + (x * x) / df)),
      chisq: (x, df) => (x <= 0 ? 0 : Math.exp((df / 2 - 1) * Math.log(x) - x / 2 - (df / 2) * Math.log(2) - lg(df / 2))),
      f: (x, d1, d2) => (x <= 0 ? 0 : Math.exp(0.5 * (d1 * Math.log(d1 * x) + d2 * Math.log(d2) - (d1 + d2) * Math.log(d1 * x + d2)) - Math.log(x) - (lg(d1 / 2) + lg(d2 / 2) - lg((d1 + d2) / 2)))),
    };
  }

  function el(tag, attrs, parent) {
    const e = document.createElementNS(NS, tag);
    for (const k in attrs) e.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(e);
    return e;
  }
  function svgRoot(label) {
    const s = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": label, class: "viz-svg" });
    return s;
  }

  function nullPlot(node) {
    const type = node.dataset.viz;
    const stat = parseFloat(node.dataset.stat);
    const df = parseFloat(node.dataset.df), df1 = parseFloat(node.dataset.df1), df2 = parseFloat(node.dataset.df2);
    const D = densities(), S = SM();
    let pdf, p, crit, lo, hi, twoSided, name;
    if (type === "t") {
      pdf = (x) => D.t(x, df); twoSided = true; name = `t distribution, df = ${Number.isInteger(df) ? df : fmt(df, 1)}`;
      p = S.tTwoTail(stat, df); crit = S.tQuantile(0.975, df);
    } else if (type === "z") {
      pdf = D.z; twoSided = true; name = "Standard normal (z) distribution";
      p = 2 * S.normSf(Math.abs(stat)); crit = 1.959964;
    } else if (type === "chisq") {
      pdf = (x) => D.chisq(x, df); twoSided = false; name = `χ² distribution, df = ${df}`;
      p = S.chisqSf(stat, df); crit = qUpper((x) => S.chisqSf(x, df), 0.05);
    } else {
      pdf = (x) => D.f(x, df1, df2); twoSided = false; name = `F distribution, df = ${df1}, ${fmt(df2, 1)}`;
      const fSf = (x) => S.betainc(df2 / (df2 + df1 * x), df2 / 2, df1 / 2);
      p = fSf(stat); crit = qUpper(fSf, 0.05);
    }
    if (twoSided) { const m = Math.max(4, Math.abs(stat) * 1.15, crit * 1.4); lo = -m; hi = m; }
    else { lo = 0; hi = Math.max(crit * 1.8, stat * 1.2, type === "chisq" ? df + 4 * Math.sqrt(2 * df) : 4); }

    const n = 240, xs = [], ys = [];
    for (let i = 0; i <= n; i++) { const x = lo + ((hi - lo) * i) / n; xs.push(x); ys.push(pdf(x)); }
    let ymax = Math.max(...ys.filter(Number.isFinite));
    if (!Number.isFinite(ymax) || ymax <= 0) ymax = 1;
    if (!twoSided && ys[0] > ymax * 0.9 && ys.length > 3) ymax = Math.max(...ys.slice(3)) * 1.1; // df ≤ 2 spikes at 0
    const X = (x) => PAD.l + ((x - lo) / (hi - lo)) * (W - PAD.l - PAD.r);
    const Y = (y) => H - PAD.b - (Math.min(y, ymax) / ymax) * (H - PAD.t - PAD.b);

    const svg = svgRoot(`${name}. Observed statistic ${fmt(stat)}; shaded area ${fmtP(p)}.`);
    const area = (a, b) => {
      const pts = [`${X(a)},${Y(0)}`];
      for (let i = 0; i <= n; i++) { const x = a + ((b - a) * i) / n; pts.push(`${X(x)},${Y(pdf(x))}`); }
      pts.push(`${X(b)},${Y(0)}`);
      return pts.join(" ");
    };
    // rejection region (α = 0.05)
    if (twoSided) {
      el("polygon", { points: area(lo, -crit), class: "viz-reject" }, svg);
      el("polygon", { points: area(crit, hi), class: "viz-reject" }, svg);
    } else el("polygon", { points: area(crit, hi), class: "viz-reject" }, svg);
    // p-value area
    const s = Math.abs(stat);
    if (twoSided) {
      if (s < hi) { el("polygon", { points: area(s, hi), class: "viz-p" }, svg); el("polygon", { points: area(lo, -s), class: "viz-p" }, svg); }
    } else if (stat < hi) el("polygon", { points: area(Math.max(stat, lo), hi), class: "viz-p" }, svg);
    // curve + axis
    el("polyline", { points: xs.map((x, i) => `${X(x)},${Y(ys[i])}`).join(" "), class: "viz-curve" }, svg);
    el("line", { x1: PAD.l, x2: W - PAD.r, y1: Y(0), y2: Y(0), class: "viz-axis" }, svg);
    const ticks = twoSided ? [-crit, 0, crit] : [0, crit];
    ticks.forEach((t) => {
      el("line", { x1: X(t), x2: X(t), y1: Y(0), y2: Y(0) + 5, class: "viz-axis" }, svg);
      const lab = el("text", { x: X(t), y: Y(0) + 18, class: "viz-tick", "text-anchor": "middle" }, svg);
      lab.textContent = fmt(t);
    });
    // observed statistic marker (clamped to the plot)
    const sx = X(Math.max(lo, Math.min(hi, stat)));
    el("line", { x1: sx, x2: sx, y1: PAD.t, y2: Y(0), class: "viz-obs" }, svg);
    const lbl = el("text", { x: sx + (stat > (lo + hi) / 2 ? -6 : 6), y: PAD.t + 12, class: "viz-obs-label", "text-anchor": stat > (lo + hi) / 2 ? "end" : "start" }, svg);
    lbl.textContent = `observed ${type === "chisq" ? "χ²" : type} = ${fmt(stat)}${stat > hi || stat < lo ? " (off the chart)" : ""}`;

    const cap = document.createElement("div");
    cap.className = "viz-caption";
    cap.innerHTML = `<span class="viz-key viz-key-p"></span> p-value area: <strong>${fmtP(p)}</strong>
      <span class="viz-key viz-key-r"></span> rejection region at α = 0.05 (critical value ${fmt(crit)})
      <div class="viz-note">${name}. If H₀ were true, statistics would follow this curve; the shaded tail is how often a result at least this extreme would occur by chance.</div>`;
    node.replaceChildren(svg, cap);
  }

  // bisection for an upper-tail quantile
  function qUpper(sf, alpha) {
    let a = 0, b = 1;
    while (sf(b) > alpha) b *= 2;
    for (let i = 0; i < 80; i++) { const m = (a + b) / 2; if (sf(m) > alpha) a = m; else b = m; }
    return (a + b) / 2;
  }

  function dotPlot(node) {
    let groups;
    try { groups = JSON.parse(node.dataset.groups); } catch { node.textContent = "Invalid diagram data."; return; }
    const names = Object.keys(groups);
    const all = names.flatMap((k) => groups[k]);
    const zero = node.dataset.zero === "1";
    let lo = Math.min(...all, ...(zero ? [0] : [])), hi = Math.max(...all, ...(zero ? [0] : []));
    const padv = (hi - lo || 1) * 0.12; lo -= padv; hi += padv;
    const rowH = Math.min(56, (H - PAD.t - PAD.b) / names.length);
    const h = PAD.t + rowH * names.length + PAD.b;
    const LW = 110;
    const X = (v) => LW + ((v - lo) / (hi - lo)) * (W - LW - PAD.r);
    const svg = el("svg", { viewBox: `0 0 ${W} ${h}`, role: "img", class: "viz-svg", "aria-label": `Dot plot of ${names.join(", ")}` });
    names.forEach((g, gi) => {
      const yc = PAD.t + rowH * gi + rowH / 2;
      const t = el("text", { x: LW - 12, y: yc + 4, class: "viz-group", "text-anchor": "end" }, svg);
      t.textContent = g;
      el("line", { x1: LW, x2: W - PAD.r, y1: yc, y2: yc, class: "viz-grid" }, svg);
      // stack ties vertically so every point is visible
      const seen = {};
      groups[g].forEach((v) => {
        const key = v.toFixed(6); seen[key] = (seen[key] || 0) + 1;
        const k = seen[key] - 1, off = (k % 2 ? 1 : -1) * Math.ceil(k / 2) * 7;
        el("circle", { cx: X(v), cy: yc + off, r: 5, class: `viz-dot viz-dot-${gi % 4}` }, svg);
      });
      const mean = groups[g].reduce((a, b) => a + b, 0) / groups[g].length;
      el("line", { x1: X(mean), x2: X(mean), y1: yc - rowH * 0.38, y2: yc + rowH * 0.38, class: "viz-mean" }, svg);
    });
    if (zero) {
      el("line", { x1: X(0), x2: X(0), y1: PAD.t - 4, y2: PAD.t + rowH * names.length, class: "viz-obs" }, svg);
      const right = X(0) > W * 0.7;
      const zl = el("text", { x: X(0) + (right ? -6 : 6), y: PAD.t + 6, class: "viz-obs-label", "text-anchor": right ? "end" : "start" }, svg);
      zl.textContent = "no change";
    }
    const axisY = PAD.t + rowH * names.length + 6;
    el("line", { x1: LW, x2: W - PAD.r, y1: axisY, y2: axisY, class: "viz-axis" }, svg);
    for (let i = 0; i <= 4; i++) {
      const v = lo + ((hi - lo) * i) / 4;
      const tx = el("text", { x: X(v), y: axisY + 18, class: "viz-tick", "text-anchor": "middle" }, svg);
      tx.textContent = fmt(v, 1);
    }
    const cap = document.createElement("div");
    cap.className = "viz-caption";
    cap.innerHTML = `Each dot is one observation; the vertical bar is the group mean${node.dataset.unit ? ` (${node.dataset.unit})` : ""}. ${node.dataset.caption || ""}`;
    node.replaceChildren(svg, cap);
  }

  // 5×5 dot-matrix loader; each dot's delay comes from its ring distance from the centre (radial ripple)
  window.dotLoader = function (label = "Loading") {
    let dots = "";
    for (let r = 0; r < 5; r++) for (let c = 0; c < 5; c++) {
      const d = Math.max(Math.abs(r - 2), Math.abs(c - 2));
      dots += `<i style="animation-delay:${d * 0.14}s"></i>`;
    }
    return `<div class="dot-loader" role="status"><span class="dot-grid" aria-hidden="true">${dots}</span><span class="dot-label">${label}…</span></div>`;
  };

  function draw(node) {
    try { node.dataset.viz === "dots" ? dotPlot(node) : nullPlot(node); }
    catch (e) { node.textContent = "Diagram unavailable."; console.error(e); }
    node.classList.add("viz-in");
  }

  // Each diagram shows the dot loader until it scrolls into view, ripples briefly, then fades in.
  window.renderViz = function (root) {
    if (!SM()) return;
    const nodes = [...root.querySelectorAll(".viz[data-viz]")];
    const calm = matchMedia("(prefers-reduced-motion: reduce)").matches;
    if (calm || !("IntersectionObserver" in window)) { nodes.forEach(draw); return; }
    const io = new IntersectionObserver((entries) => entries.forEach((e) => {
      if (!e.isIntersecting) return;
      io.unobserve(e.target);
      setTimeout(() => draw(e.target), 650);
    }), { rootMargin: "0px 0px -10% 0px" });
    nodes.forEach((n) => { n.innerHTML = window.dotLoader("Drawing the graph"); io.observe(n); });
  };
})();
