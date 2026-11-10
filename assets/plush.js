/*
 * Mochi-style plush mascots, one animal per topic.
 *   PLUSH.svg("bear", { label })      -> SVG string (viewBox 0 0 100 112)
 *   PLUSH.forCategory("means")        -> { kind, name, body, line, inner, blurb }
 * Soft egg body, lavender-ish outline, white sticker border, glossy eyes, pink blush.
 */
(function () {
  "use strict";
  const EYE = "#2a2160", IRIS = "#4b3c9a", BLUSH = "#ffa3c4";

  const ANIMALS = {
    bear:    { body: "#f3e6da", line: "#c9a98e", inner: "#f7c6d4" },
    bunny:   { body: "#f6f3ff", line: "#bba9ec", inner: "#ffc4d8" },
    cat:     { body: "#e9e4f7", line: "#9f8ee0", inner: "#ffb3d1" },
    owl:     { body: "#cbbcf3", line: "#8f7bdc", inner: "#f5edff" },
    chick:   { body: "#fff1b8", line: "#e3bf55", inner: "#ffb36b" },
    frog:    { body: "#cdeedb", line: "#7cc39c", inner: "#f7fff9" },
    penguin: { body: "#c9dcf2", line: "#7f9fd0", inner: "#ffffff" },
    hamster: { body: "#ffdcc8", line: "#e3a083", inner: "#fff4ec" },
    seal:    { body: "#e3ecf5", line: "#8fa8c8", inner: "#f7fbff" },
    duck:    { body: "#fff6d6", line: "#e3bf55", inner: "#ffb36b" },
    fox:     { body: "#ffd2b3", line: "#e3946a", inner: "#fff6ef" },
    sloth:   { body: "#eadfd1", line: "#b89c80", inner: "#a88b6e" },
    koala:   { body: "#dcdde6", line: "#9a9cb2", inner: "#ffc4d8" },
  };

  const CATEGORY = {
    foundations:   { kind: "seal",    name: "Null the seal",      blurb: "Your friendly guide to the basics." },
    design:        { kind: "hamster", name: "Power the hamster",  blurb: "Stuffs samples into its cheeks." },
    assumptions:   { kind: "frog",    name: "Shapiro the frog",   blurb: "Hops in to check if you're normal." },
    means:         { kind: "bear",    name: "Tee the bear",       blurb: "Carries two cups, always comparing." },
    posthoc:       { kind: "duck",    name: "Tukey the duck",     blurb: "Lines up the ducklings and compares each pair." },
    nonparametric: { kind: "chick",   name: "Wilcox the chick",   blurb: "Doesn't need normal. Just ranks." },
    categorical:   { kind: "bunny",   name: "Chi the bunny",      blurb: "Counts everything in boxes." },
    regression:    { kind: "owl",     name: "Reggie the owl",     blurb: "Draws the best line through anything." },
    survival:      { kind: "cat",     name: "Kaplan the cat",     blurb: "Has nine lives and tracks every one." },
    agreement:     { kind: "penguin", name: "Kappa the penguin",  blurb: "Waddles in pairs and agrees a lot." },
    diagnostic:    { kind: "fox",     name: "Roc the fox",        blurb: "Sniffs out who's sick and who's not." },
    equivalence:   { kind: "sloth",   name: "Toast the sloth",    blurb: "Close enough is good enough." },
    special:       { kind: "koala",   name: "Boots the koala",    blurb: "Resamples the eucalyptus, again and again." },
  };

  const BODY = "M50 30 C24 30 14 56 14 76 C14 96 30 106 50 106 C70 106 86 96 86 76 C86 56 76 30 50 30Z";
  const FROG_BODY = "M50 36 C22 36 12 58 12 76 C12 96 30 106 50 106 C70 106 88 96 88 76 C88 58 78 36 50 36Z";

  function ears(k, a) {
    const L = a.line, B = a.body, I = a.inner, s = `stroke="${L}" stroke-width="2.5"`;
    switch (k) {
      case "bear": return `<circle cx="25" cy="34" r="12" fill="${B}" ${s}/><circle cx="75" cy="34" r="12" fill="${B}" ${s}/><circle cx="25" cy="34" r="6" fill="${I}"/><circle cx="75" cy="34" r="6" fill="${I}"/>`;
      case "bunny": return `<ellipse cx="37" cy="18" rx="8" ry="22" transform="rotate(-10 37 18)" fill="${B}" ${s}/><ellipse cx="63" cy="18" rx="8" ry="22" transform="rotate(10 63 18)" fill="${B}" ${s}/><ellipse cx="37" cy="20" rx="3.5" ry="15" transform="rotate(-10 37 20)" fill="${I}"/><ellipse cx="63" cy="20" rx="3.5" ry="15" transform="rotate(10 63 20)" fill="${I}"/>`;
      case "cat": return `<path d="M20 44 L24 16 L44 32 Z" fill="${B}" ${s} stroke-linejoin="round"/><path d="M80 44 L76 16 L56 32 Z" fill="${B}" ${s} stroke-linejoin="round"/><path d="M25 38 L27 23 L37 32 Z" fill="${I}"/><path d="M75 38 L73 23 L63 32 Z" fill="${I}"/>`;
      case "fox": return `<path d="M18 46 L22 12 L46 32 Z" fill="${B}" ${s} stroke-linejoin="round"/><path d="M82 46 L78 12 L54 32 Z" fill="${B}" ${s} stroke-linejoin="round"/><path d="M21 22 L22 13 L29 19 Z" fill="#7a4a3a"/><path d="M79 22 L78 13 L71 19 Z" fill="#7a4a3a"/>`;
      case "owl": return `<path d="M22 40 L22 16 L40 30 Z" fill="${L}"/><path d="M78 40 L78 16 L60 30 Z" fill="${L}"/>`;
      case "chick": return `<path d="M48 26 q-2 -10 6 -12 q-4 6 2 10" fill="${L}"/>`;
      case "duck": return `<path d="M50 28 q-6 -10 2 -14 q-1 7 5 9 q-4 1 -7 5" fill="${L}"/>`;
      case "frog": return `<circle cx="33" cy="34" r="14" fill="${B}" ${s}/><circle cx="67" cy="34" r="14" fill="${B}" ${s}/>`;
      case "hamster": return `<circle cx="27" cy="36" r="9" fill="${B}" ${s}/><circle cx="73" cy="36" r="9" fill="${B}" ${s}/><circle cx="27" cy="36" r="4.5" fill="#ffc4d8"/><circle cx="73" cy="36" r="4.5" fill="#ffc4d8"/>`;
      case "koala": return `<circle cx="20" cy="42" r="16" fill="${B}" ${s}/><circle cx="80" cy="42" r="16" fill="${B}" ${s}/><circle cx="20" cy="43" r="9" fill="${I}"/><circle cx="80" cy="43" r="9" fill="${I}"/>`;
      default: return "";
    }
  }

  function extras(k, a) {
    switch (k) {
      case "owl": return `<ellipse cx="50" cy="88" rx="20" ry="17" fill="${a.inner}"/><circle cx="37" cy="62" r="11" fill="none" stroke="#ffd27a" stroke-width="3"/><circle cx="63" cy="62" r="11" fill="none" stroke="#ffd27a" stroke-width="3"/><line x1="48" y1="62" x2="52" y2="62" stroke="#ffd27a" stroke-width="3"/><path d="M47 74 L53 74 L50 79 Z" fill="#ffbf47"/>`;
      case "chick": return `<path d="M46 73 L54 73 L50 78 Z" fill="${a.inner}"/><ellipse cx="18" cy="78" rx="6" ry="9" transform="rotate(25 18 78)" fill="${a.body}" stroke="${a.line}" stroke-width="2"/><ellipse cx="82" cy="78" rx="6" ry="9" transform="rotate(-25 82 78)" fill="${a.body}" stroke="${a.line}" stroke-width="2"/>`;
      case "duck": return `<ellipse cx="50" cy="75" rx="10" ry="5" fill="${a.inner}"/><path d="M42 75 h16" stroke="#e08f45" stroke-width="1.2"/>`;
      case "penguin": return `<path d="M50 40 C30 40 24 62 26 80 C28 98 40 102 50 102 C60 102 72 98 74 80 C76 62 70 40 50 40Z" fill="${a.inner}"/><path d="M46 73 L54 73 L50 78 Z" fill="#ffbf47"/>`;
      case "hamster": return `<ellipse cx="24" cy="78" rx="10" ry="8" fill="${a.inner}"/><ellipse cx="76" cy="78" rx="10" ry="8" fill="${a.inner}"/><ellipse cx="50" cy="94" rx="16" ry="10" fill="${a.inner}"/>`;
      case "bear": return `<ellipse cx="50" cy="77" rx="10" ry="7" fill="#fbf3ec"/><ellipse cx="50" cy="74" rx="3" ry="2.2" fill="${EYE}"/>`;
      case "fox": return `<path d="M50 70 C40 70 24 74 22 86 C30 98 42 100 50 92 C58 100 70 98 78 86 C76 74 60 70 50 70Z" fill="${a.inner}"/><ellipse cx="50" cy="75" rx="3" ry="2.2" fill="${EYE}"/>`;
      case "seal": return `<ellipse cx="50" cy="76" rx="9" ry="6" fill="${a.inner}"/><ellipse cx="50" cy="73" rx="2.6" ry="2" fill="${EYE}"/><g fill="${a.line}"><circle cx="44" cy="77" r=".9"/><circle cx="41" cy="75" r=".9"/><circle cx="56" cy="77" r=".9"/><circle cx="59" cy="75" r=".9"/></g>`;
      case "sloth": return `<ellipse cx="37" cy="63" rx="11" ry="8" transform="rotate(-20 37 63)" fill="${a.inner}"/><ellipse cx="63" cy="63" rx="11" ry="8" transform="rotate(20 63 63)" fill="${a.inner}"/><ellipse cx="50" cy="72" rx="3" ry="2" fill="${EYE}"/>`;
      case "koala": return `<ellipse cx="50" cy="72" rx="6" ry="8" fill="#4a4660"/><ellipse cx="48" cy="69" rx="1.6" ry="2.4" fill="#fff" opacity=".5"/>`;
      default: return "";
    }
  }

  function eyes(k) {
    const ey = k === "frog" ? 34 : 62, dx = k === "frog" ? 17 : 13;
    return [50 - dx, 50 + dx].map((x) => `<ellipse cx="${x}" cy="${ey}" rx="7" ry="8" fill="${EYE}"/><ellipse cx="${x}" cy="${ey + 3}" rx="5" ry="4" fill="${IRIS}"/><circle cx="${x + 2}" cy="${ey - 3}" r="2.8" fill="#fff"/><circle cx="${x - 2.5}" cy="${ey + 3}" r="1.3" fill="#fff"/>`).join("");
  }

  function mouth(k) {
    if (["owl", "chick", "penguin", "duck"].includes(k)) return "";
    const y = { bear: 80, fox: 79, seal: 78, sloth: 76, koala: 83, frog: 60 }[k] || 73;
    return `<path d="M45 ${y} q2.5 3.5 5 0 q2.5 3.5 5 0" fill="none" stroke="${EYE}" stroke-width="1.8" stroke-linecap="round"/>`;
  }

  function svg(kind, opts = {}) {
    const a = ANIMALS[kind] || ANIMALS.bear;
    const body = kind === "frog" ? FROG_BODY : BODY;
    const by = kind === "frog" ? 48 : 72;
    const e = ears(kind, a);
    // white sticker border: the same silhouette, fat white stroke, drawn underneath
    const border = e.replace(/stroke="[^"]*"/g, "").replace(/stroke-width="[^"]*"/g, "").replace(/fill="[^"]*"/g, "");
    const label = opts.label ? `role="img" aria-label="${opts.label}"` : `aria-hidden="true" focusable="false"`;
    return `<svg class="plush ${opts.cls || ""}" viewBox="0 0 100 112" ${label}>
<g fill="#fff" stroke="#fff" stroke-width="9" stroke-linejoin="round">${border}<path d="${body}"/></g>
${e}<path d="${body}" fill="${a.body}" stroke="${a.line}" stroke-width="2.5"/>
<ellipse cx="34" cy="44" rx="9" ry="5" transform="rotate(-25 34 44)" fill="#fff" opacity=".55"/>
${extras(kind, a)}${eyes(kind)}
<ellipse cx="${kind === "frog" ? 28 : 27}" cy="${by}" rx="6" ry="3.5" fill="${BLUSH}" opacity=".85"/><ellipse cx="${kind === "frog" ? 72 : 73}" cy="${by}" rx="6" ry="3.5" fill="${BLUSH}" opacity=".85"/>
${mouth(kind)}</svg>`;
  }

  window.PLUSH = {
    ANIMALS, CATEGORY, svg,
    forCategory: (id) => { const c = CATEGORY[id] || CATEGORY.foundations; return { ...c, ...ANIMALS[c.kind] }; },
  };
})();
