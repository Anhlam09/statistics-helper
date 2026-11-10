/*
 * Claw machines for the home page.
 *   CLAW.mount(el, { manifest, available, loadPage })
 * "Test" machine grabs a random lesson; "Story" machine reads a study description and
 * grabs the matching test. The prize is a flip card: plush + name on the front,
 * a nutrition-label style "Test facts" panel on the back.
 * CLAW.matchStory(text) is exposed for testing.
 */
(function () {
  "use strict";
  const esc = (s) => String(s).replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
  const reduced = () => matchMedia("(prefers-reduced-motion: reduce)").matches;
  const wait = (ms) => new Promise((r) => setTimeout(r, reduced() ? Math.min(ms, 60) : ms));

  /* ---------- Story -> test ---------- */
  const has = (s, re) => re.test(s);
  const PAIRED = /before.{0,40}after|after.{0,40}before|same (patients|people|subjects|participants|children|women|men)|paired|matched|twice|baseline and|pre.?(and|vs|\/).?post|crossover|each (patient|subject|eye)/;
  const MANY = /\b(three|four|five|six|3|4|5|6)\b.{0,25}(groups|arms|doses|treatments|diets|regimens|hospitals|sites)|more than two|several (groups|doses|arms|treatments)|multiple (groups|doses|arms)/;
  const BINARY = /yes\/no|yes or no|binary|respond(ed|ers)?\b|response rate|proportion|percent|%|cured|infection|infected|readmi|event rate|success|complication|positive or negative|smok(er|ing)\b|vaccinat/;
  const SKEW = /ordinal|likert|pain score|severity (scale|score|grade)|ranks?\b|ranked|skew|not normal|non-?normal|outliers?/;
  const SMALL = /small|few|rare|under 5|less than 5|fewer than 5|tiny/;
  const ADJUST = /adjust|covariat|controlling for|predict|risk factors?|several factors|multiple factors|while accounting/;

  const RULES = [
    [/compet(ing)? risks?|competing event/, "competing-risks", "time to an event where another event can get in the way first"],
    [/surviv|time to|time-to|until (death|relapse|progression)|relapse|progression|censor|hazard|died during follow|follow.?up until/, (s) =>
      has(s, /sample size|how many (patients|participants|subjects)|power/) ? ["sample-size-survival", "planning how many patients a time-to-event study needs"]
      : has(s, ADJUST) || has(s, /hazard ratio/) ? ["cox-ph", "a time-to-event outcome where you want to adjust for other factors"]
      : has(s, /curve|median|estimate|probability of surviving|one group/) ? ["kaplan-meier", "estimating survival over time, with censoring"]
      : ["log-rank", "comparing time-to-event between groups, with censoring"]],
    [/meta.?analy|pool(ed|ing)? (the )?(studies|trials|results)|several (published )?(studies|trials)/, "meta-analysis", "combining results from several studies"],
    [/(repeated|over time|longitudinal|visits|months|weeks).{0,80}(dropouts?|missing)|(dropouts?|missing).{0,80}(repeated|over time|longitudinal)/, "mixed-models", "repeated measurements over time, with dropouts or missing visits"],
    [/missing (data|values)|incomplete data/, "multiple-imputation", "a dataset with missing values"],
    [/propensity|observational.{0,60}(treat|exposure|drug)|confound/, "propensity-scores", "comparing treatments in observational data where groups differ at baseline"],
    [/sample size|how many (patients|participants|subjects|people)|enough power|power calculation/, "sample-size-means-proportions", "planning how many participants you need"],
    [/equivalen|non-?inferior|as good as|no worse than|biosimilar|generic version/, "tost-equivalence", "showing a treatment is not worse, or about the same"],
    [/sensitivity|specificity|diagnos|\broc\b|\bauc\b|biomarker|screening test|detect(s|ing)? (the )?disease|cut.?off/, "roc-auc", "how well a test or marker separates sick from healthy"],
    [/raters?|readers?|radiologists|pathologists|observers|agree(ment)?|two (methods|devices|machines)|reliab|questionnaire items|internal consistency|new device/, (s) =>
      has(s, /internal consistency|questionnaire items|items? (of|on|in) (a|the) (scale|questionnaire)|cronbach/) ? ["cronbach-alpha", "checking that questionnaire items hang together"]
      : has(s, /bland|two (methods|devices|machines)|new device|measur(e|ing) the same|limits of agreement/) ? ["bland-altman", "two methods measuring the same continuous quantity"]
      : has(s, /icc|reliab|continuous|scores?|measurements?/) ? ["icc", "agreement between raters on a continuous measurement"]
      : ["cohens-kappa", "two raters classifying the same cases into categories"]],
    [/shapiro|is (it|my data|the data|this) normal|normally distributed|check(ing)? normality/, "shapiro-wilk", "checking whether data look normally distributed"],
    [/equal variances?|same variance|variances (are|differ)|levene/, "levene", "checking whether groups have similar spread"],
    [/(each|every) (dose|group|arm|treatment) (to|vs\.?|versus|against) (the )?(control|placebo)|compared? (with|to) (a|one|the) (shared |common )?control/, "dunnett", "comparing several groups against one control"],
    [/which (groups|ones|doses|arms) differ|post.?hoc|pairwise|after (an |the )?anova/, "tukey-hsd", "finding which groups differ after an overall test"],
    [/multiple (tests|comparisons|endpoints|outcomes)|many (tests|outcomes|genes)|false discovery|bonferroni/, "multiple-comparisons", "running many tests and controlling false positives"],
    [/bootstrap|permutation|resampl/, "bootstrap-permutation", "resampling when formulas don't fit"],
    [/count(s|ed)?\b|number of (\w+ )?(events|visits|falls|seizures|attacks|admissions|episodes|exacerbations)|rate per|per (person|patient).?year/, "poisson-regression", "a count outcome, like number of events per patient"],
  ];

  function binaryRule(s) {
    if (has(s, PAIRED)) return ["mcnemar", "a yes/no outcome measured twice on the same people"];
    if (has(s, /strat|centers|centres|sites|hospitals/) && !has(s, ADJUST)) return ["cmh", "a yes/no outcome compared across groups, stratified by center"];
    if (has(s, ADJUST) || has(s, /odds/)) return ["logistic-regression", "a yes/no outcome predicted from several factors"];
    if (has(s, /trend|dose.?response|increasing (dose|levels)/)) return ["trend-tests", "a yes/no outcome across ordered dose levels"];
    if (has(s, SMALL)) return ["fisher-exact", "a yes/no outcome in a small sample"];
    if (has(s, /one group|single group|single arm|historical|benchmark|known rate|expected rate/)) return ["one-proportion", "one group's proportion compared with a known value"];
    return ["chi-square-independence", "a yes/no or categorical outcome compared between groups"];
  }

  function matchStory(text) {
    const s = " " + text.toLowerCase() + " ";
    for (const [re, out, why] of RULES) {
      if (!re.test(s)) continue;
      if (typeof out === "function") { const [slug, w] = out(s); return { slug, why: w }; }
      return { slug: out, why };
    }
    if (has(s, /correlat|relationship|associated with|association between|related to|link between/) && !has(s, ADJUST))
      return has(s, SKEW) ? { slug: "spearman", why: "a relationship between two ranked or skewed variables" } : { slug: "pearson", why: "a straight-line relationship between two continuous variables" };
    if (has(s, BINARY)) { const [slug, why] = binaryRule(s); return { slug, why }; }
    if (has(s, SKEW)) {
      if (has(s, PAIRED)) return { slug: "wilcoxon-signed-rank", why: "skewed or ordinal data measured twice on the same people" };
      if (has(s, MANY)) return { slug: "kruskal-wallis", why: "skewed or ordinal data across three or more groups" };
      return { slug: "mann-whitney", why: "skewed or ordinal data in two independent groups" };
    }
    if (has(s, ADJUST)) return has(s, /baseline/) ? { slug: "ancova", why: "comparing groups while adjusting for a baseline value" } : { slug: "linear-regression", why: "a continuous outcome predicted from several factors" };
    if (has(s, /repeated|over time|several visits|multiple visits|longitudinal|every (week|month)|at (weeks|months)/))
      return has(s, /missing|dropout|mmrm|mixed/) ? { slug: "mixed-models", why: "repeated measurements over time with dropouts" } : { slug: "repeated-measures-anova", why: "the same people measured at three or more time points" };
    if (has(s, PAIRED)) return { slug: "paired-t", why: "a continuous measurement taken twice on the same people" };
    if (has(s, /two factors|factorial|and (by )?sex|by (sex|gender) and/)) return { slug: "two-way-anova", why: "a continuous outcome with two factors at once" };
    if (has(s, MANY)) return { slug: "one-way-anova", why: "comparing a continuous outcome across three or more groups" };
    if (has(s, /one group|single group|known value|reference value|target value|national average|normal value/)) return { slug: "one-sample-t", why: "one group's mean compared with a known value" };
    if (has(s, /two groups|placebo|versus|\bvs\.?\b|compare|control group|men and women|treated and untreated|drug a|new drug/)) return { slug: "independent-t", why: "a continuous outcome compared between two separate groups" };
    return null;
  }

  /* ---------- Lesson snippets ---------- */
  const plain = (md) => String(md || "").replace(/\[([^\]]+)\]\([^)]+\)/g, "$1").replace(/[*_`]/g, "").replace(/<[^>]+>/g, "").trim();

  function glance(md) {
    const lines = md.split("\n");
    const i = lines.findIndex((l) => /^\|\s*when to use/i.test(l));
    if (i < 0 || !lines[i + 2]) return null;
    return lines[i + 2].split("|").slice(1, -1).map((c) => plain(c));
  }

  /* ---------- Collection (per-browser) ---------- */
  const COLL = "plushCollection";
  const getColl = () => { try { return new Set(JSON.parse(localStorage.getItem(COLL) || "[]")); } catch { return new Set(); } };
  const addColl = (slug) => { const c = getColl(); c.add(slug); try { localStorage.setItem(COLL, JSON.stringify([...c])); } catch { /* storage unavailable */ } return c; };

  /* ---------- UI ---------- */
  function machineHTML(mode) {
    const story = mode === "story";
    return `<div class="machine machine-${mode}">
      <div class="machine-top">${story ? "Describe your study in a sentence or two, and we'll find the test that fits" : "Grab a test you want to learn"}</div>
      ${story ? "" : `<div class="machine-glass">
        <div class="claw-rope"></div><div class="claw"><span></span><span></span></div>
        <div class="pile"></div>
        <div class="chute">prize</div>
      </div>`}
      <div class="machine-controls">
        ${story ? `<label class="sr-only" for="storyInput">Describe your study</label>
          <textarea id="storyInput" rows="3" placeholder="We measured blood pressure in 20 patients before and after 8 weeks on a new drug."></textarea>
          <div class="claw-err" role="alert"></div>
          <button class="btn claw-go">Find my test</button>`
          : `<p class="machine-hint">Tap any plush to grab its test, or let the claw pick one for you.</p><button class="btn claw-go">Surprise me</button>`}
      </div>
    </div>`;
  }

  function mount(el, { manifest, available, loadPage }) {
    const tests = manifest.tests.filter((t) => available.has(t.slug));
    el.innerHTML = `<div class="claw-machines">${machineHTML("test")}${machineHTML("story")}</div>
      <div class="prize-area" aria-live="polite"></div>`;
    const prize = el.querySelector(".prize-area");
    const PER = 30; // extras stay hidden until the glass is wide enough
    const shown = () => new Set([...el.querySelectorAll(".pile-item")].map((i) => i.dataset.slug));
    const pickNew = () => { const s = shown(), pool = tests.filter((t) => !s.has(t.slug)); return pool[Math.floor(Math.random() * pool.length)] || tests[0]; };
    const fill = (item, t) => {
      item.dataset.slug = t.slug;
      item.setAttribute("aria-label", `Grab ${t.name}`);
      item.innerHTML = `${PLUSH.svg(PLUSH.forCategory(t.category).kind)}<span class="pile-tag">${esc((window.SHORTNAME || {})[t.slug] || t.name)}</span>`;
    };
    el.querySelectorAll(".machine").forEach((m) => {
      const pile = m.querySelector(".pile");
      if (!pile) return;
      for (let i = 0; i < PER; i++) {
        const b = document.createElement("button");
        b.className = "pile-item"; b.type = "button";
        fill(b, pickNew());
        pile.appendChild(b);
      }
      pile.addEventListener("click", (e) => {
        const b = e.target.closest(".pile-item"); if (!b) return;
        const t = tests.find((x) => x.slug === b.dataset.slug);
        if (t) grab(m, t);
      });
    });

    // show only as many plushies as fit in three rows of the glass
    const fit = () => el.querySelectorAll(".pile").forEach((pile) => {
      const cap = Math.max(1, Math.floor(pile.clientWidth / 72)) * 3;
      pile.querySelectorAll(".pile-item").forEach((it, i) => { it.hidden = i >= cap; });
    });
    fit();
    new ResizeObserver(fit).observe(el);

    let busy = false;
    async function grab(machine, test, why) {
      if (busy) return; busy = true;
      prize.innerHTML = "";
      machine.after(prize); // the prize card lands right under the machine you played
      if (!machine.querySelector(".machine-glass")) { await showPrize(test, why); busy = false; return; }
      const glass = machine.querySelector(".machine-glass"), claw = machine.querySelector(".claw"), rope = machine.querySelector(".claw-rope");
      const pile = machine.querySelector(".pile");
      let item = pile.querySelector(`.pile-item[data-slug="${test.slug}"]:not([hidden])`);
      if (!item) { // the story's test isn't in this machine yet: swap it in
        const items = pile.querySelectorAll(".pile-item:not([hidden])");
        item = items[Math.floor(Math.random() * items.length)];
        fill(item, test);
      }
      const gr = glass.getBoundingClientRect(), ir = item.getBoundingClientRect();
      const x = ir.left - gr.left + ir.width / 2, depth = ir.top - gr.top - 6;
      const set = (cx, cy) => { claw.style.left = `${cx}px`; rope.style.left = `${cx}px`; rope.style.height = `${cy}px`; claw.style.top = `${cy}px`; };
      machine.classList.add("running");
      set(x, 14); await wait(700);
      set(x, depth); await wait(700);
      claw.classList.add("closed"); await wait(250);
      item.classList.add("lifted");
      const lift = depth - 18;
      item.style.transform = `translateY(${-lift}px)`;
      set(x, 14 - 0); await wait(700);
      const chuteX = 46; // centre of the prize chute
      item.style.transform = `translate(${chuteX - x}px, ${-lift}px)`;
      set(chuteX, 14); await wait(700);
      claw.classList.remove("closed");
      item.style.transform = `translate(${chuteX - x}px, 40px)`; item.style.opacity = "0";
      await wait(450);
      await showPrize(test, why);
      item.style.transition = "none"; item.style.transform = ""; item.classList.remove("lifted");
      fill(item, pickNew()); // a fresh plush drops into the empty spot
      requestAnimationFrame(() => { item.style.transition = ""; item.style.opacity = ""; });
      machine.classList.remove("running");
      busy = false;
    }

    async function showPrize(t, why) {
      const md = await loadPage(`content/tests/${t.slug}.md`).catch(() => "");
      const g = md ? glance(md) : null;
      const [rLine, pyLine] = (window.KEYCODE || {})[t.slug] || [];
      const p = PLUSH.forCategory(t.category);
      const cat = manifest.categories.find((c) => c.id === t.category);
      const num = manifest.tests.findIndex((x) => x.slug === t.slug) + 1;
      const isNew = !getColl().has(t.slug);
      addColl(t.slug);
      const yes = (b) => (b ? "Yes" : "No");
      const rows = [
        ["Outcome", t.outcome], ["Groups", t.groups], ["Paired", yes(t.paired)],
        ["Type", t.parametric === true ? "Parametric" : t.parametric === false ? "Nonparametric" : "Either"],
        ...(g ? [["Use it for", g[0]], ["Assumes", g[2]], ["Example", g[3]]] : []),
      ].filter(([, v]) => v && v !== "n/a");
      prize.innerHTML = `
        ${why ? `<p class="prize-why">We picked this because your study sounds like <strong>${esc(why)}</strong>.</p>` : ""}
        <div class="prize">
          <button class="flip-card" aria-label="Flip the ${esc(t.name)} card">
            <span class="flip-inner">
              <span class="flip-front" style="--plush-bg:${p.body}">
                <span class="card-top"><span class="dots"><i></i><i></i><i></i></span><span>No. ${String(num).padStart(2, "0")} / ${manifest.tests.length}</span></span>
                <span class="card-art">${PLUSH.svg(p.kind)}<span class="card-plush-name">${esc(p.name)}</span></span>
                <span class="card-name">${esc(t.name.replace(/\s*\(.*\)\s*$/, ""))}</span>
                <span class="card-sub">${esc(cat?.name || "")}<span>tap to flip</span></span>
                ${isNew ? `<span class="card-new">New!</span>` : ""}
              </span>
              <span class="flip-back">
                <span class="facts-h">Test facts</span>
                <span class="facts-serving">Serving: 1 study question</span>
                ${rows.map(([k, v]) => `<span class="facts-row"><span>${esc(k)}</span><span>${esc(v)}</span></span>`).join("")}
                <span class="facts-foot">Daily value of confusion: 0%</span>
              </span>
            </span>
          </button>
          <div class="prize-info">
            <div class="eyebrow">${esc(cat?.name || "")}</div>
            <h3>${esc(t.name)}</h3>
            <p>${esc(t.summary)}</p>
            ${g ? `<p><strong>When to use it:</strong> ${esc(g[0])}</p>` : ""}
            ${rLine ? `<div class="prize-code"><span>R</span><code>${esc(rLine)}</code></div>` : ""}
            ${pyLine ? `<div class="prize-code"><span>Python</span><code>${esc(pyLine)}</code></div>` : ""}
            <div class="btn-row"><a class="btn" href="#/test/${t.slug}">Read the full lesson →</a><button class="btn ghost prize-again">Play again</button></div>
          </div>
        </div>`;
      const card = prize.querySelector(".flip-card");
      card.addEventListener("click", () => card.classList.toggle("flipped"));
      prize.querySelector(".prize-again").addEventListener("click", () => el.querySelector(".machine-test .claw-go").focus());
      prize.scrollIntoView({ behavior: reduced() ? "auto" : "smooth", block: "nearest" });
    }

    const testMachine = el.querySelector(".machine-test"), storyMachine = el.querySelector(".machine-story");
    testMachine.querySelector(".claw-go").addEventListener("click", () => {
      const items = [...testMachine.querySelectorAll(".pile-item:not([hidden])")];
      const pick = items[Math.floor(Math.random() * items.length)];
      const t = pick && tests.find((x) => x.slug === pick.dataset.slug);
      if (t) grab(testMachine, t);
    });
    const input = storyMachine.querySelector("textarea"), err = storyMachine.querySelector(".claw-err");
    input.addEventListener("input", () => { err.textContent = ""; });
    storyMachine.querySelector(".claw-go").addEventListener("click", () => {
      const s = input.value.trim();
      if (!s) { err.textContent = "Type a short description of your study first."; input.focus(); return; }
      const hit = matchStory(s);
      const t = hit && tests.find((x) => x.slug === hit.slug);
      if (!t) { err.innerHTML = `We need a few more clues. Mention your outcome (a measurement, yes/no, time to an event) and your groups, or try the <a href="#/choose">step-by-step chooser</a>.`; return; }
      grab(storyMachine, t, hit.why);
    });
  }

  window.CLAW = { mount, matchStory, glance };
})();
