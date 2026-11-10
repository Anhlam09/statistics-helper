(function () {
  "use strict";

  const $ = (sel, el = document) => el.querySelector(sel);
  const main = $("#main");
  const sidebar = $("#sidebar");
  const LANG_LABEL = { r: "R", python: "Python", sas: "SAS" };

  const store = {
    get(k) { try { return localStorage.getItem(k); } catch { return null; } },
    set(k, v) { try { localStorage.setItem(k, v); } catch { /* storage unavailable */ } },
  };

  const esc = (s) => String(s).replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
  const slugify = (s) => s.toLowerCase().replace(/<[^>]+>/g, "").replace(/[^\w]+/g, "-").replace(/^-|-$/g, "");

  let manifest = { categories: [], tests: [] };
  const available = new Set();      // slugs with a written lesson
  const pageCache = new Map();

  /* ---------- Data loading ---------- */
  async function loadManifest() {
    const res = await fetch("content/manifest.json", { cache: "no-cache" });
    manifest = await res.json();
    // The manifest is the source of truth (tools/check_content.py verifies every entry has a
    // lesson file). Mark an unfinished lesson with "soon": true to show it as "coming soon".
    for (const t of manifest.tests) if (!t.soon) available.add(t.slug);
  }

  async function loadPage(path) {
    if (pageCache.has(path)) return pageCache.get(path);
    const res = await fetch(path, { cache: "no-cache" });
    if (!res.ok) throw new Error(res.status);
    const md = await res.text();
    pageCache.set(path, md);
    return md;
  }

  const testBySlug = (slug) => manifest.tests.find((t) => t.slug === slug);
  const catById = (id) => manifest.categories.find((c) => c.id === id);

  /* ---------- Sidebar ---------- */
  function renderSidebar(activeSlug) {
    let html = `<div class="mobile-nav">
      <h4>Navigate</h4>
      <a href="#/choose">Which test?</a><a href="#/tests">All tests</a><a href="#/guide">How to choose</a><a href="#/calculators">Calculators</a><a href="#/resources">Further learning</a><a href="#/about">About</a>
    </div>`;
    for (const c of manifest.categories) {
      html += `<h4><span class="h4-plush">${PLUSH.svg(PLUSH.forCategory(c.id).kind)}</span>${esc(c.name)}</h4>`;
      for (const t of manifest.tests.filter((t) => t.category === c.id)) {
        const cls = [t.slug === activeSlug ? "active" : "", available.has(t.slug) ? "" : "pending"].join(" ");
        html += `<a class="${cls}" href="#/test/${t.slug}">${esc(shortName(t.name))}</a>`;
      }
    }
    sidebar.innerHTML = html;
    sidebar.querySelector("a.active")?.scrollIntoView({ block: "nearest" });
  }
  const shortName = (n) => n.replace(/\s*\(.*\)\s*$/, "");

  function setActiveNav(key) {
    document.querySelectorAll(".topnav a").forEach((a) => a.classList.toggle("active", a.dataset.nav === key));
  }

  /* ---------- Markdown → enhanced HTML ---------- */
  function renderMarkdown(md) {
    const html = marked.parse(md, { gfm: true });
    const wrap = document.createElement("div");
    wrap.className = "article";
    wrap.innerHTML = html;

    // the summary blockquote right after the H1 is the page lead, not a callout
    const h1 = wrap.querySelector("h1");
    if (h1?.nextElementSibling?.tagName === "BLOCKQUOTE") h1.nextElementSibling.classList.add("lead");

    // ids on headings
    wrap.querySelectorAll("h2, h3").forEach((h) => { if (!h.id) h.id = slugify(h.textContent); });

    // tables scroll horizontally on narrow screens
    wrap.querySelectorAll("table").forEach((t) => {
      const d = document.createElement("div"); d.className = "table-scroll";
      t.replaceWith(d); d.appendChild(t);
    });

    // code blocks: highlight + copy buttons
    wrap.querySelectorAll("pre > code").forEach((code) => {
      const pre = code.parentElement;
      const lang = (code.className.match(/language-(\w+)/) || [])[1];
      const box = document.createElement("div");
      box.className = "code-block" + (lang === "text" ? " output" : "");
      pre.replaceWith(box); box.appendChild(pre);
      if (lang && lang !== "text" && window.hljs && hljs.getLanguage(lang)) hljs.highlightElement(code);
      if (lang !== "text") {
        const btn = document.createElement("button");
        btn.className = "copy-btn"; btn.type = "button"; btn.textContent = "Copy";
        btn.addEventListener("click", async () => {
          try { await navigator.clipboard.writeText(code.textContent); btn.textContent = "Copied ✓"; }
          catch { btn.textContent = "Select & copy"; }
          setTimeout(() => (btn.textContent = "Copy"), 1500);
        });
        box.appendChild(btn);
      }
    });

    buildLangTabs(wrap);

    // "At a glance" table: the first table whose header starts with "When to use"
    const glance = [...wrap.querySelectorAll("table")].find((t) => /^when to use/i.test(t.querySelector("th")?.textContent.trim() || ""));
    if (glance) {
      const box = document.createElement("div");
      box.className = "glance";
      const holder = glance.closest(".table-scroll") || glance;
      holder.replaceWith(box); box.appendChild(holder);
    }

    // Learning objectives / Key takeaways become tinted boxes
    [["learning-objectives", "box-objectives"], ["key-takeaways", "box-takeaways"]].forEach(([id, cls]) => {
      const h = wrap.querySelector(`h2#${id}`);
      if (!h) return;
      const box = document.createElement("section");
      box.className = `box ${cls}`;
      h.before(box);
      let n = h.nextElementSibling;
      box.appendChild(h);
      while (n && n.tagName !== "H2") { const next = n.nextElementSibling; box.appendChild(n); n = next; }
    });

    if (window.renderViz) window.renderViz(wrap);
    return wrap;
  }

  // Turn the "How to run it" section's ### R / ### Python / ### SAS into tabs
  function buildLangTabs(wrap) {
    const h2 = [...wrap.querySelectorAll("h2")].find((h) => /how to run/i.test(h.textContent));
    if (!h2) return;
    const panels = [];
    let node = h2.nextElementSibling, current = null;
    const toRemove = [];
    while (node && node.tagName !== "H2") {
      const next = node.nextElementSibling;
      if (node.tagName === "H3") {
        const label = node.textContent.trim();
        const lang = /python/i.test(label) ? "python" : /sas/i.test(label) ? "sas" : /^r\b/i.test(label) ? "r" : null;
        if (lang) { current = { lang, nodes: [] }; panels.push(current); toRemove.push(node); node = next; continue; }
      }
      if (current) { current.nodes.push(node); }
      node = next;
    }
    if (panels.length < 2) return;

    const preferred = store.get("lang") || "r";
    const tabs = document.createElement("div");
    tabs.className = "lang-tabs";
    const list = document.createElement("div");
    list.className = "tab-list"; list.setAttribute("role", "tablist"); list.setAttribute("aria-label", "Language");
    tabs.appendChild(list);

    panels.forEach((p, i) => {
      const btn = document.createElement("button");
      btn.type = "button"; btn.dataset.lang = p.lang; btn.textContent = LANG_LABEL[p.lang];
      btn.setAttribute("role", "tab"); btn.id = `tab-${p.lang}`; btn.setAttribute("aria-controls", `panel-${p.lang}`);
      list.appendChild(btn);
      const panel = document.createElement("div");
      panel.className = "tab-panel"; panel.id = `panel-${p.lang}`; panel.dataset.lang = p.lang;
      panel.setAttribute("role", "tabpanel"); panel.setAttribute("aria-labelledby", btn.id);
      if (p.lang === "sas") {
        const note = document.createElement("p");
        note.className = "callout-warn untested-note";
        note.innerHTML = "<strong>Not yet tested in SAS.</strong> This code follows standard SAS/STAT syntax but hasn't been run, so check the output before relying on it.";
        panel.appendChild(note);
      }
      p.nodes.forEach((n) => panel.appendChild(n));
      tabs.appendChild(panel);
    });
    toRemove.forEach((n) => n.remove());
    h2.after(tabs);

    const select = (lang) => {
      if (!panels.some((p) => p.lang === lang)) lang = panels[0].lang;
      list.querySelectorAll("button").forEach((b) => {
        const on = b.dataset.lang === lang;
        b.setAttribute("aria-selected", on); b.tabIndex = on ? 0 : -1;
      });
      tabs.querySelectorAll(".tab-panel").forEach((p) => (p.hidden = p.dataset.lang !== lang));
    };
    list.addEventListener("click", (e) => {
      const b = e.target.closest("button"); if (!b) return;
      store.set("lang", b.dataset.lang); select(b.dataset.lang);
    });
    list.addEventListener("keydown", (e) => {
      if (!["ArrowLeft", "ArrowRight"].includes(e.key)) return;
      const btns = [...list.querySelectorAll("button")];
      const i = btns.indexOf(document.activeElement);
      const n = btns[(i + (e.key === "ArrowRight" ? 1 : btns.length - 1)) % btns.length];
      n.focus(); n.click();
    });
    select(preferred);
  }

  function buildToc(article) {
    const hs = [...article.querySelectorAll("h2")];
    if (hs.length < 3) return null;
    const toc = document.createElement("nav");
    toc.className = "toc"; toc.setAttribute("aria-label", "On this page");
    toc.innerHTML = `<h5>On this page</h5>` + hs.map((h) => `<a href="#${h.id}" data-target="${h.id}">${esc(h.textContent)}</a>`).join("");
    toc.addEventListener("click", (e) => {
      const a = e.target.closest("a"); if (!a) return;
      e.preventDefault();
      document.getElementById(a.dataset.target)?.scrollIntoView({ behavior: "smooth" });
    });
    const obs = new IntersectionObserver((entries) => {
      entries.forEach((en) => {
        if (en.isIntersecting) {
          toc.querySelectorAll("a").forEach((a) => a.classList.toggle("active", a.dataset.target === en.target.id));
        }
      });
    }, { rootMargin: "-70px 0px -70% 0px" });
    hs.forEach((h) => obs.observe(h));
    return toc;
  }

  /* ---------- Views ---------- */
  // Same "test of the day" for everyone on a given date
  function todaysTest() {
    const list = manifest.tests.filter((t) => available.has(t.slug) && t.category !== "foundations");
    const t = list[Math.floor(Date.now() / 864e5) % list.length];
    return { slug: t.slug, label: (window.SHORTNAME || {})[t.slug] || shortName(t.name) };
  }

  // The p-value café lives in the footer, bottom right
  function drawFooterCafe() {
    if (manifest.tests.length) $("#footerCafe").innerHTML = CAFE.full(todaysTest());
  }

  function viewHome() {
    setActiveNav("");
    const count = available.size;
    main.innerHTML = `
      <section class="home-hero welcome">
        <div class="cloud-puffs" aria-hidden="true"><span></span><span><i class="cloud-face"><b></b><b></b><em></em><em></em><u></u></i></span><span></span><span></span></div>
        <div>
          <h1>Hi, welcome to your learning café!</h1>
          <p>A friendly guide to help you pick the right statistical test, with plain-language explanations and R and Python code for every one.</p>
          <div class="btn-row">
            <a class="btn" href="#claw" data-scroll="claw">Play the claw machine</a>
            <a class="btn ghost" href="#/tests">Browse all ${manifest.tests.length} tests</a>
          </div>
        </div>
      </section>
      <section class="claw-section" id="claw">
        <div id="clawMount"></div>
      </section>
`;
    CLAW.mount($("#clawMount"), { manifest, available, loadPage });
    main.onclick = (e) => {
      const a = e.target.closest('[href="#claw"]');
      if (!a) return;
      e.preventDefault();
      $("#claw").scrollIntoView({ behavior: matchMedia("(prefers-reduced-motion: reduce)").matches ? "auto" : "smooth" });
      $("#clawMount .claw-go")?.focus({ preventScroll: true });
    };
  }

  // Pull the "At a glance" row (When to use | Data type | Key assumptions | Example) out of a lesson
  function glanceRow(md) {
    const lines = md.split("\n");
    const i = lines.findIndex((l) => /^\|\s*when to use/i.test(l));
    if (i < 0 || !lines[i + 2]) return null;
    return lines[i + 2].split("|").slice(1, -1).map((c) => c.trim());
  }

  async function drawSelectorTable(host, list) {
    host.innerHTML = window.dotLoader("Loading the table");
    const rows = await Promise.all(list.map(async (t) => {
      if (!available.has(t.slug)) return { t, cells: null };
      try { return { t, cells: glanceRow(await loadPage(`content/tests/${t.slug}.md`)) }; } catch { return { t, cells: null }; }
    }));
    const cell = (c) => (c ? marked.parseInline(c) : `<span style="color:var(--muted)">—</span>`);
    host.innerHTML = `<div class="table-scroll"><table class="selector-table">
      <thead><tr><th>Test</th><th>When to use</th><th>Data type</th><th>Key assumptions</th><th>Example</th></tr></thead>
      <tbody>${rows.map(({ t, cells }) => `<tr><td><a href="#/test/${t.slug}">${esc(shortName(t.name))}</a></td>
        ${cells ? cells.slice(0, 4).map((c) => `<td>${cell(c)}</td>`).join("") : `<td colspan="4">${esc(t.summary)}</td>`}</tr>`).join("")}</tbody>
    </table></div>`;
  }

  function viewTests(params) {
    setActiveNav("tests");
    let active = params.get("cat") || "all";
    let mode = params.get("view") === "table" ? "table" : "cards";
    const draw = () => {
      const list = manifest.tests.filter((t) => active === "all" || t.category === active);
      main.innerHTML = `
        <h1>All tests</h1>
        <p style="color:var(--text-2)">${manifest.tests.length} lessons across ${manifest.categories.length} topics. Not sure where to start? <a href="#/choose">Use the guided chooser</a>.</p>
        <div class="view-toggle" role="group" aria-label="View">
          <button data-view="cards" aria-pressed="${mode === "cards"}">Cards</button>
          <button data-view="table" aria-pressed="${mode === "table"}">Comparison table</button>
        </div>
        <div class="filters" role="group" aria-label="Filter by topic">
          <button aria-pressed="${active === "all"}" data-cat="all">All</button>
          ${manifest.categories.map((c) => `<button aria-pressed="${active === c.id}" data-cat="${c.id}">${esc(c.name)}</button>`).join("")}
        </div>
        <div id="testList">${mode === "cards" ? `<div class="grid">
          ${list.map((t) => `
            <a class="card ${available.has(t.slug) ? "" : "pending"}" href="#/test/${t.slug}">
              <h3>${esc(t.name)}</h3>
              <p>${esc(t.summary)}</p>
              <div class="count">${esc(catById(t.category)?.name || "")}${t.paired ? " · paired" : ""}${available.has(t.slug) ? "" : " · coming soon"}</div>
            </a>`).join("")}
        </div>` : ""}</div>`;
      if (mode === "table") drawSelectorTable(main.querySelector("#testList"), list.filter((t) => t.category !== "foundations"));
      const sync = () => {
        const q = new URLSearchParams();
        if (active !== "all") q.set("cat", active);
        if (mode === "table") q.set("view", "table");
        history.replaceState(null, "", "#/tests" + (q.toString() ? "?" + q : ""));
      };
      main.querySelector(".filters").addEventListener("click", (e) => {
        const b = e.target.closest("button"); if (!b) return;
        active = b.dataset.cat; sync(); draw();
      });
      main.querySelector(".view-toggle").addEventListener("click", (e) => {
        const b = e.target.closest("button"); if (!b) return;
        mode = b.dataset.view; sync(); draw();
      });
    };
    draw();
  }

  async function viewTest(slug) {
    setActiveNav("");
    const t = testBySlug(slug);
    if (!t) return viewNotFound();
    renderSidebar(slug);
    document.title = `${shortName(t.name)} · Statistics Helper`;
    if (!available.has(slug)) {
      main.innerHTML = `<h1>${esc(t.name)}</h1><p class="lead">${esc(t.summary)}</p>
        <div class="callout-warn">This page hasn't been written yet. Check back soon, or <a href="#/tests">browse the available tests</a>.</div>`;
      return;
    }
    main.innerHTML = window.dotLoader("Loading the lesson");
    try {
      const md = await loadPage(`content/tests/${slug}.md`);
      const article = renderMarkdown(md);
      const cat = catById(t.category);
      const chips = document.createElement("div");
      chips.className = "meta-chips";
      chips.innerHTML = `<span class="chip">Outcome: ${esc(t.outcome)}</span>
        ${t.groups && t.groups !== "n/a" && t.groups !== "any" ? `<span class="chip">Groups: ${esc(t.groups)}</span>` : ""}
        ${t.paired ? `<span class="chip">Paired / repeated</span>` : ""}
        ${t.parametric === true ? `<span class="chip">Parametric</span>` : t.parametric === false ? `<span class="chip">Nonparametric</span>` : ""}`;
      const hero = document.createElement("header");
      hero.className = "lesson-hero";
      hero.innerHTML = `<div class="eyebrow"><a href="#/tests?cat=${t.category}">${esc(cat?.name || "")}</a></div>`;
      const h1 = article.querySelector("h1"), lead = article.querySelector("blockquote.lead");
      if (h1) hero.appendChild(h1);
      if (lead) hero.appendChild(lead);
      hero.appendChild(chips);
      const mascot = PLUSH.forCategory(t.category);
      hero.insertAdjacentHTML("beforeend", `<div class="hero-plush">${PLUSH.svg(mascot.kind)}<span>${esc(mascot.name)}</span></div>`);
      article.prepend(hero);
      embedCalculators(article, slug);
      const wrap = document.createElement("div");
      wrap.className = "article-wrap";
      wrap.appendChild(article);
      const toc = buildToc(article);
      if (toc) wrap.appendChild(toc);
      main.replaceChildren(wrap);
    } catch {
      main.innerHTML = `<h1>${esc(t.name)}</h1><div class="callout-warn">Couldn't load this page. If you opened index.html directly from disk, run a local server instead (see README).</div>`;
    }
  }

  async function viewGuide(file = "content/guide.md", nav = "guide", title = "How to choose") {
    setActiveNav(nav);
    document.title = `${title} · Statistics Helper`;
    try {
      const article = renderMarkdown(await loadPage(file));
      const wrap = document.createElement("div");
      wrap.className = "article-wrap";
      wrap.appendChild(article);
      const toc = buildToc(article);
      if (toc) wrap.appendChild(toc);
      main.replaceChildren(wrap);
    } catch {
      main.innerHTML = `<div class="callout-warn">Couldn't load the guide.</div>`;
    }
  }

  /* ---------- Wizard ---------- */
  function viewChoose() {
    setActiveNav("choose");
    document.title = "Which test? · Statistics Helper";
    const { NODES, resolveCompare } = window.WIZARD;
    const path = [];          // [{ nodeId, option }]
    const answers = () => Object.fromEntries(path.filter((p) => p.option.key).map((p) => [p.nodeId, p.option.key]));

    const draw = (result) => {
      const crumbs = path.map((p, i) => `<button class="chip" data-back="${i}" title="Change this answer">${esc(p.option.label)} ✕</button>`).join("");
      let body;
      if (result) {
        body = `<h2>Recommended</h2>
          ${result.map((x, i) => {
            const t = testBySlug(x.slug);
            if (!t) return "";
            return `<div class="result-card ${i ? "secondary" : ""}">
              <h3><a href="#/test/${t.slug}">${esc(t.name)}</a></h3>
              <p class="why">${esc(x.why)}</p>
              <a class="btn ${i ? "ghost" : ""}" href="#/test/${t.slug}">${available.has(t.slug) ? "Open the guide" : "Coming soon"} →</a>
            </div>`;
          }).join("")}
          <div class="wizard-actions">
            <button class="btn ghost" data-act="back">← Back</button>
            <button class="btn ghost" data-act="restart">Start over</button>
          </div>`;
      } else {
        const nodeId = path.length ? path[path.length - 1].option.next : "start";
        const node = NODES[nodeId];
        body = `<h2>${esc(node.q)}</h2>
          <div class="options">
            ${node.options.map((o, i) => `<button class="option" data-opt="${i}" data-node="${nodeId}">
              <strong>${esc(o.label)}</strong>${o.hint ? `<span>${esc(o.hint)}</span>` : ""}</button>`).join("")}
          </div>
          ${path.length ? `<div class="wizard-actions"><button class="btn ghost" data-act="back">← Back</button></div>` : ""}`;
      }
      main.innerHTML = `<div class="wizard">
        <h1>Which test should I use?</h1>
        <p style="color:var(--text-2)">Answer a few questions about your study design. You can click any answer below to change it.</p>
        ${path.length ? `<div class="crumbs">${crumbs}</div>` : ""}
        ${body}</div>`;
      main.querySelector(".option")?.focus({ preventScroll: true });
    };

    const step = () => {
      const last = path[path.length - 1];
      if (!last) return draw();
      const o = last.option;
      if (o.result) return draw(o.result);
      if (o.leaf) return draw(resolveCompare(answers()));
      const next = NODES[o.next];
      if (next.skipIf && next.skipIf(answers())) return draw(resolveCompare(answers()));
      draw();
    };

    main.onclick = (e) => {
      const opt = e.target.closest("[data-opt]");
      if (opt) {
        const node = NODES[opt.dataset.node];
        path.push({ nodeId: opt.dataset.node, option: node.options[+opt.dataset.opt] });
        return step();
      }
      const back = e.target.closest("[data-back]");
      if (back) { path.length = +back.dataset.back; return step(); }
      const act = e.target.closest("[data-act]")?.dataset.act;
      if (act === "back") { path.pop(); step(); }
      if (act === "restart") { path.length = 0; step(); }
    };
    step();
  }

  /* ---------- Calculators ---------- */
  const calcs = () => window.CALCULATORS || [];

  function mountCalc(calc, el, preset) {
    try { calc.mount(el, preset); }
    catch (err) { el.innerHTML = `<div class="callout-warn">This calculator failed to load.</div>`; console.error(err); }
  }

  // Each lesson's own example data, so "Try it yourself" starts from the numbers the R/Python code printed
  const LESSON_PRESETS = {
    "paired-t": { "ttest-paired": { mean: -6.4, sd: 2.756810, n: 10, mu0: 0 } },
    "one-sample-t": { "ttest-paired": { mean: 7.125, sd: 0.4136863, n: 12, mu0: 7.0 } },
    "independent-t": { "ttest-ind": { m1: -15.75, s1: 13.30840, n1: 12, m2: -3.9, s2: 2.514403, n2: 10, welch: "1" } },
    "two-proportion-z": { twobytwo: { a: 60, b: 40, c: 42, d: 58 } },
    "fisher-exact": { twobytwo: { a: 2, b: 13, c: 8, d: 6 } },
    "mcnemar": { mcnemar: { b: 20, c: 6 } },
    "one-proportion": { "one-proportion": { x: 14, n: 40, p0: 0.2 } },
    "roc-auc": { diagnostic: { tp: 24, fp: 14, fn: 1, tn: 21 } },
    "sample-size-means-proportions": { "ss-means": { delta: 5, sd: 12, alpha: 0.05, power: 0.8 }, "ss-props": { p1: 0.40, p2: 0.55, alpha: 0.05, power: 0.8 } },
    "sample-size-survival": { "ss-survival": { hr: 0.7, alpha: 0.05, power: 0.8, p: 0.5, eventProb: 0.6678 } },
  };

  // Insert a "Try it yourself" section before "Interpreting the output"
  function embedCalculators(article, slug) {
    const related = calcs().filter((c) => (c.related || []).includes(slug));
    if (!related.length) return;
    const anchor = [...article.querySelectorAll("h2")].find((h) => /interpreting/i.test(h.textContent));
    const h2 = document.createElement("h2");
    h2.id = "try-it-yourself"; h2.textContent = "Try it yourself";
    const intro = document.createElement("p");
    intro.textContent = LESSON_PRESETS[slug]
      ? "This starts with the lesson's example data, so the results match the output above. Change any number to see how the answer responds, e.g., halve the sample size or double the SD."
      : "Enter your own summary numbers. Results update as you type and use the same methods as the R code above.";
    const nodes = [h2, intro];
    related.forEach((c) => {
      const box = document.createElement("section");
      box.className = "calc-embed";
      box.setAttribute("aria-label", c.title);
      box.innerHTML = `<h3>${esc(c.title)}</h3>`;
      const body = document.createElement("div");
      box.appendChild(body);
      mountCalc(c, body, (LESSON_PRESETS[slug] || {})[c.id]);
      nodes.push(box);
    });
    if (anchor) nodes.forEach((n) => anchor.before(n));
    else nodes.forEach((n) => article.appendChild(n));
  }

  function viewCalculators() {
    setActiveNav("calculators");
    document.title = "Calculators · Statistics Helper";
    const list = calcs();
    main.innerHTML = `
      <h1>Calculators</h1>
      <p style="color:var(--text-2)">Quick calculations from summary numbers. Every calculator is tested against R (the same functions used on the test pages), runs entirely in your browser, and shows the equivalent R call.</p>
      ${list.length ? `<div class="grid">${list.map((c) => `
        <a class="card" href="#/calc/${c.id}"><h3>${esc(c.title)}</h3><p>${esc(c.summary || "")}</p></a>`).join("")}</div>`
        : `<div class="callout-warn">Calculators are not available.</div>`}`;
  }

  function viewCalc(id) {
    setActiveNav("calculators");
    const c = calcs().find((x) => x.id === id);
    if (!c) return viewNotFound();
    document.title = `${c.title} · Statistics Helper`;
    const links = (c.related || []).map(testBySlug).filter(Boolean)
      .map((t) => `<a class="chip cat" href="#/test/${t.slug}">${esc(shortName(t.name))}</a>`).join(" ");
    main.innerHTML = `<div class="calc-page">
      <p><a href="#/calculators">← All calculators</a></p>
      <h1>${esc(c.title)}</h1>
      <p style="color:var(--text-2)">${esc(c.summary || "")}</p>
      ${links ? `<div class="meta-chips">Learn more: ${links}</div>` : ""}
      <div class="calc-body"></div></div>`;
    mountCalc(c, main.querySelector(".calc-body"));
  }

  function viewNotFound() {
    main.innerHTML = `<h1>Page not found</h1><p><a href="#/">Go home</a></p>`;
  }

  /* ---------- Router ---------- */
  async function route() {
    const raw = location.hash.slice(1) || "/";
    // in-page anchors (from TOC links in rendered markdown) don't start with "/"
    if (!raw.startsWith("/")) return;
    const [path, query] = raw.split("?");
    const params = new URLSearchParams(query || "");
    const parts = path.split("/").filter(Boolean);
    main.onclick = null;
    document.title = "Statistics Helper";
    closeMenu();
    if (parts[0] !== "test") renderSidebar(null);
    if (!parts.length) viewHome();
    else if (parts[0] === "tests") viewTests(params);
    else if (parts[0] === "test" && parts[1]) await viewTest(parts[1]);
    else if (parts[0] === "choose") viewChoose();
    else if (parts[0] === "guide") await viewGuide();
    else if (parts[0] === "resources") await viewGuide("content/resources.md", "resources", "Further learning");
    else if (parts[0] === "about") await viewGuide("content/about.md", "about", "About");
    else if (parts[0] === "calculators") viewCalculators();
    else if (parts[0] === "calc" && parts[1]) viewCalc(parts[1]);
    else viewNotFound();
    window.scrollTo(0, 0);
    main.focus({ preventScroll: true });
  }

  /* ---------- Search ---------- */
  function setupSearch() {
    const input = $("#search"), box = $("#searchResults");
    let sel = -1;
    const score = (t, q) => {
      const hay = `${t.name} ${t.summary} ${catById(t.category)?.name} ${t.slug.replace(/-/g, " ")}`.toLowerCase();
      return q.split(/\s+/).every((w) => hay.includes(w)) ? (t.name.toLowerCase().includes(q) ? 2 : 1) : 0;
    };
    const update = () => {
      const q = input.value.trim().toLowerCase();
      if (!q) { box.hidden = true; return; }
      const hits = manifest.tests.map((t) => [t, score(t, q)]).filter((x) => x[1]).sort((a, b) => b[1] - a[1]).slice(0, 8);
      sel = -1;
      box.innerHTML = hits.length
        ? hits.map(([t]) => `<a href="#/test/${t.slug}">${esc(t.name)}<small>${esc(t.summary)}</small></a>`).join("")
        : `<a href="#/choose">No match. Try the guided chooser →</a>`;
      box.hidden = false;
    };
    input.addEventListener("input", update);
    input.addEventListener("keydown", (e) => {
      const links = [...box.querySelectorAll("a")];
      if (e.key === "ArrowDown" || e.key === "ArrowUp") {
        e.preventDefault();
        sel = (sel + (e.key === "ArrowDown" ? 1 : -1) + links.length) % links.length;
        links.forEach((a, i) => a.classList.toggle("sel", i === sel));
      } else if (e.key === "Enter") {
        const a = links[sel] || links[0];
        if (a) { location.hash = a.getAttribute("href"); input.value = ""; box.hidden = true; input.blur(); }
      } else if (e.key === "Escape") { box.hidden = true; input.blur(); }
    });
    box.addEventListener("click", () => { input.value = ""; box.hidden = true; });
    document.addEventListener("click", (e) => { if (!e.target.closest(".search-wrap")) box.hidden = true; });
    document.addEventListener("keydown", (e) => {
      if (e.key === "/" && !/input|textarea/i.test(document.activeElement.tagName)) { e.preventDefault(); input.focus(); }
    });
  }

  /* ---------- Theme & menu ---------- */
  function applyTheme(theme) {
    if (theme) document.documentElement.dataset.theme = theme;
    const dark = theme ? theme === "dark" : matchMedia("(prefers-color-scheme: dark)").matches;
    $("#hljs-dark").disabled = !dark;
    $("#hljs-light").disabled = dark;
  }
  function setupTheme() {
    applyTheme(store.get("theme"));
    $("#themeBtn").addEventListener("click", () => {
      const cur = document.documentElement.dataset.theme ||
        (matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light");
      const next = cur === "dark" ? "light" : "dark";
      store.set("theme", next); applyTheme(next);
    });
    matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => applyTheme(store.get("theme")));
  }
  function closeMenu() { sidebar.classList.remove("open"); $("#menuBtn").setAttribute("aria-expanded", "false"); }
  function setupMenu() {
    $("#menuBtn").addEventListener("click", (e) => {
      e.stopPropagation();
      const open = sidebar.classList.toggle("open");
      $("#menuBtn").setAttribute("aria-expanded", open);
    });
    document.addEventListener("click", (e) => {
      if (sidebar.classList.contains("open") && !e.target.closest("#sidebar")) closeMenu();
    });
  }

  /* ---------- Boot ---------- */
  (async function init() {
    setupTheme(); setupMenu(); setupSearch();
    try { await loadManifest(); }
    catch {
      main.innerHTML = `<div class="callout-warn"><strong>Couldn't load content.</strong> Browsers block loading local files when you open index.html directly. Start a local server in this folder: <code>python3 -m http.server 8000</code>, then open <a href="http://localhost:8000">http://localhost:8000</a>.</div>`;
      return;
    }
    drawFooterCafe();
    window.addEventListener("hashchange", route);
    route();
  })();
})();
