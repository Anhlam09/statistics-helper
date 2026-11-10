/*
 * The p-value café: a cutaway sticker-style coffee shop staffed by the topic plushies.
 *   CAFE.full()   -> SVG string for the home page hero
 * Clickable bits in the full café carry data-href so the home page can route them.
 */
(function () {
  "use strict";
  const P = () => window.PLUSH;
  // nested plush at (x, y) = top-left, w = width
  const pl = (kind, x, y, w, cls = "bob") =>
    `<g class="${cls}"><svg x="${x}" y="${y}" width="${w}" height="${w * 1.12}" viewBox="0 0 100 112">${P().svg(kind).replace(/^<svg[^>]*>|<\/svg>$/g, "")}</svg></g>`;

  // Storefront: big windows into the café, a door that opens as plush customers come and go.
  // window with an arched top
  const archWin = (x, y, w, h) => `<path d="M${x} ${y + h} V${y + 14} a${w / 2} 14 0 0 1 ${w} 0 V${y + h} Z" fill="#fff8ec" stroke="#e08aa4" stroke-width="3"/>`;

  // Two-storey storefront. Upstairs: a study loft. Downstairs: the barista making coffee,
  // customers stepping in and out of the door, runners on the sidewalk.
  function full(today = { slug: "paired-t", label: "Paired t" }) {
    let lights = "";
    for (let i = 0; i < 15; i++) lights += `<circle class="lt" cx="${46 + i * 31}" cy="${84 + (i % 2 ? 3 : 0)}" r="3" fill="${["#ffd27a", "#ffb3d1", "#ffb8cc"][i % 3]}"/>`;
    let awn = "";
    for (let i = 0; i < 10; i++) awn += `<path d="M${30 + i * 46} 90 h46 v12 a23 11 0 0 1 -46 0Z" fill="${i % 2 ? "#ffffff" : "#ffb8cc"}" stroke="#e08aa4" stroke-width="1.2"/>`;
    let flowers = "";
    [[58, "#ffb3d1"], [74, "#ffe5b4"], [90, "#ffb8cc"], [106, "#ffb3d1"], [122, "#c9e4ff"], [138, "#ffe5b4"], [154, "#ffb3d1"], [170, "#ffb8cc"],
     [350, "#ffe5b4"], [366, "#ffb3d1"], [382, "#c9e4ff"], [398, "#ffb8cc"], [414, "#ffb3d1"], [430, "#ffe5b4"], [446, "#ffb8cc"], [462, "#ffb3d1"]]
      .forEach(([x, c], i) => { flowers += `<circle cx="${x}" cy="${240 + (i % 2) * 3}" r="5" fill="#9fd8b5"/><circle cx="${x}" cy="${236 + (i % 2) * 3}" r="3.4" fill="${c}"/>`; });
    let boxes = "";
    [[60, 112], [348, 112]].forEach(([x, w]) => {
      boxes += `<rect x="${x}" y="188" width="${w}" height="7" rx="3" fill="#c9a98e"/>`;
      for (let i = 0; i < 7; i++) boxes += `<circle cx="${x + 8 + i * 16}" cy="186" r="3.4" fill="${["#ffb3d1", "#ffe5b4", "#c9e4ff"][i % 3]}"/>`;
    });
    let rail = "";
    for (let i = 0; i < 9; i++) rail += `<rect x="${202 + i * 14}" y="176" width="2.5" height="14" fill="#6b3a5c"/>`;
    // customers: 10 walk in and 10 walk out on a 20 s loop; each steps up into the doorway and disappears
    // inside (or steps out of it), so you can see people going in and out
    const KINDS = ["hamster", "duck", "fox", "koala", "seal", "sloth", "frog", "chick", "bunny", "cat", "owl", "penguin", "bear"];
    let cust = "";
    for (let i = 0; i < 10; i++) {
      const lane = (i % 3) * 4;
      cust += `<g class="cust cust-in" style="animation-delay:${-i * 2}s">${pl(KINDS[i % 13], 245, 250 + lane, 30, "hop")}</g>`;
      // leaving customers carry a coffee (right) and a pastry (left)
      const y = 252 + lane;
      cust += `<g class="cust cust-out" style="animation-delay:${-i * 2 - 1}s"><g class="hop">${pl(KINDS[(i + 6) % 13], 245, y, 30, "none")}
<rect x="271" y="${y + 18}" width="7" height="9" rx="1.5" fill="#fff" stroke="#c2477a" stroke-width=".8"/><rect x="271" y="${y + 21}" width="7" height="3" fill="#ffb3d1"/><path d="M274 ${y + 16} q-2 -3 0 -6" stroke="#e8a0b8" stroke-width="1" fill="none"/>
<path d="M240 ${y + 25} q5 -9 10 0 q-5 3 -10 0Z" fill="#f2c27b" stroke="#c98a45" stroke-width=".8"/><circle cx="245" cy="${y + 19.5}" r="1.4" fill="#ff8fb0"/></g></g>`;
    }
    const runners = [["chick", 0, 7, "r"], ["frog", 2, 9, "l"], ["bear", 4, 8, "r"], ["bunny", 1, 10, "l"], ["cat", 6, 7, "r"], ["owl", 3, 11, "l"], ["duck", 5, 9, "r"], ["seal", 7, 8, "l"]];
    const run = runners.map(([k, d, dur, dir]) => `<g class="run run-${dir}" style="animation-delay:${-d}s;animation-duration:${dur}s">${pl(k, dir === "r" ? -40 : 520, 272, 24, "hop")}</g>`).join("");
    return `<svg class="cafe-svg" viewBox="24 30 472 382" role="img" aria-label="The p-value café: a two-storey shop with a study loft upstairs, a barista making coffee, and plush customers walking in and out">
<rect x="30" y="34" width="460" height="58" rx="16" fill="#6b3a5c"/>
<rect x="160" y="42" width="200" height="34" rx="10" fill="#fff6e6" stroke="#ffd27a" stroke-width="2"/>
<g transform="translate(176 50)"><path d="M0 6 h14 v7 a7 7 0 0 1 -14 0Z" fill="#c2477a"/><path d="M14 8 q5 0 5 4 t-5 4" fill="none" stroke="#c2477a" stroke-width="2"/><path class="steam" d="M4 4 q-2 -3 0 -6" stroke="#c2477a" stroke-width="1.5" fill="none"/><path class="steam d" d="M10 4 q-2 -3 0 -6" stroke="#c2477a" stroke-width="1.5" fill="none"/></g>
<text x="272" y="66" text-anchor="middle" font-size="19" fill="#c2477a" font-family="Georgia,serif" font-style="italic">p-value café</text>
<!-- second floor: the study loft -->
<rect x="30" y="92" width="460" height="112" fill="#ffd9e5"/>
${lights}
${archWin(52, 104, 128, 84)}${archWin(196, 104, 128, 84)}${archWin(340, 104, 128, 84)}
<!-- left: owl reading at a desk under a lamp -->
<line x1="68" y1="124" x2="68" y2="140" stroke="#6b3a5c" stroke-width="1.5"/><path d="M60 146 a8 6 0 0 1 16 0Z" fill="#ffd27a"/>
${pl("owl", 92, 128, 40, "sway")}
<rect x="60" y="170" width="112" height="5" rx="2" fill="#c9a98e"/>
<path d="M130 168 l10 -4 l10 4 v-9 l-10 -3 l-10 3Z" fill="#fff" stroke="#e08aa4" stroke-width="1"/><line x1="140" y1="156" x2="140" y2="166" stroke="#e08aa4"/><path class="page" d="M140 164 l10 4 v-9 l-10 -3Z" fill="#fff8ec" stroke="#e08aa4" stroke-width="1"/>
<!-- middle: balcony with a cat and a hamster sipping coffee -->
<rect x="232" y="112" width="56" height="13" rx="4" fill="#6b3a5c"/><text x="260" y="122" text-anchor="middle" font-size="8" fill="#ffd27a" font-family="Georgia,serif" font-style="italic">study loft</text>
${pl("cat", 216, 138, 38, "sip")}${pl("hamster", 266, 142, 34, "sip2")}
<rect x="250" y="164" width="9" height="9" rx="2" fill="#fff" stroke="#e08aa4"/><path class="steam" d="M254 160 q-3 -5 0 -9" stroke="#e8a0b8" stroke-width="1.4" fill="none"/>
<rect x="198" y="174" width="124" height="4" rx="2" fill="#6b3a5c"/>${rail}
<!-- right: frog and chick studying stats at a laptop, chalkboard on the wall -->
<rect x="372" y="120" width="60" height="30" rx="3" fill="#2e2550" stroke="#c9a98e" stroke-width="2"/>
<rect x="380" y="138" width="6" height="8" fill="#ffb3d1"/><rect x="390" y="132" width="6" height="14" fill="#ffe5b4"/><rect x="400" y="128" width="6" height="18" fill="#c8f0d8"/>
<text x="418" y="140" text-anchor="middle" font-size="9" fill="#f6f3ff" font-style="italic">x̄</text>
${pl("frog", 350, 146, 30, "type")}${pl("chick", 426, 146, 30, "type2")}
<rect x="350" y="174" width="110" height="5" rx="2" fill="#c9a98e"/>
<rect x="394" y="163" width="22" height="11" rx="2" fill="#e6e8f0" stroke="#e08aa4" stroke-width="1"/>
${boxes}
<rect x="26" y="200" width="468" height="8" rx="3" fill="#e08aa4"/>
<!-- ground floor (drawn 110 below its original coordinates) -->
<g transform="translate(0 110)">
<rect x="30" y="98" width="460" height="152" fill="#ffe6ee"/>
${awn}
<!-- left window: the barista making coffee -->
<rect x="44" y="118" width="140" height="104" rx="10" fill="#fff8ec" stroke="#e08aa4" stroke-width="3"/>
<text x="120" y="134" text-anchor="middle" font-size="9" fill="#e08aa4" font-family="Georgia,serif" font-style="italic">fresh p-values</text>
<g class="barista">${pl("bear", 96, 146, 44)}</g>
<rect x="52" y="190" width="124" height="30" rx="3" fill="#f7c4d4"/><rect x="48" y="186" width="132" height="6" rx="3" fill="#6b3a5c"/>
<rect x="56" y="152" width="32" height="20" rx="4" fill="#e6e8f0" stroke="#e08aa4" stroke-width="1.2"/><rect x="62" y="157" width="20" height="5" rx="1" fill="#ffb8cc"/>
<rect x="70" y="172" width="5" height="4" fill="#8c4a6e"/>
<line class="drip" x1="72.5" y1="176" x2="72.5" y2="180" stroke="#8a5a3c" stroke-width="2" stroke-linecap="round"/>
<rect x="67" y="179" width="11" height="7" rx="1.5" fill="#fff" stroke="#e08aa4" stroke-width="1"/><rect class="fill" x="68" y="181" width="9" height="4" fill="#8a5a3c"/>
<path class="steam" d="M66 150 q-3 -5 0 -10" stroke="#e8a0b8" stroke-width="1.6" fill="none"/>
<g class="slidecup"><rect x="84" y="179" width="9" height="7" rx="1.5" fill="#fff" stroke="#e08aa4" stroke-width="1"/><rect x="84" y="181" width="9" height="2.5" fill="#ffb3d1"/><path d="M88 176 q-2 -3 0 -6" stroke="#e8a0b8" stroke-width="1.2" fill="none"/></g>
${["#ffb3d1", "#ffe5b4", "#c8f0d8"].map((c, i) => `<circle cx="${148 + i * 11}" cy="206" r="4.5" fill="${c}" stroke="#c9a98e" stroke-width=".8"/>`).join("")}
<!-- door, propped open for the crowd -->
<rect x="216" y="116" width="88" height="134" rx="8" fill="#6b3a5c"/>
<rect x="224" y="124" width="72" height="126" rx="5" fill="#8c4a6e"/>
<rect x="232" y="196" width="56" height="4" rx="2" fill="#6b3a5c" opacity=".5"/><circle cx="246" cy="190" r="4" fill="#ffd27a" opacity=".5"/><circle cx="274" cy="190" r="4" fill="#ffd27a" opacity=".5"/>
<g class="door"><rect x="224" y="124" width="72" height="126" rx="5" fill="#ffb8cc" stroke="#e08aa4" stroke-width="1.5"/>
<rect x="232" y="132" width="56" height="58" rx="4" fill="#fff8ec"/><circle cx="284" cy="206" r="3" fill="#ffd27a"/></g>
<line x1="248" y1="132" x2="252" y2="140" stroke="#6b3a5c" stroke-width="1"/><line x1="272" y1="132" x2="268" y2="140" stroke="#6b3a5c" stroke-width="1"/>
<rect x="244" y="140" width="32" height="14" rx="3" fill="#c2477a"/><text x="260" y="150" text-anchor="middle" font-size="8" fill="#fff6e6">OPEN</text>
<!-- right window: customers at a table -->
<rect x="336" y="118" width="140" height="104" rx="10" fill="#fff8ec" stroke="#e08aa4" stroke-width="3"/>
<line x1="350" y1="128" x2="350" y2="138" stroke="#6b3a5c"/><path d="M342 142 a8 5 0 0 1 16 0Z" fill="#ffd27a"/>
<line x1="462" y1="128" x2="462" y2="138" stroke="#6b3a5c"/><path d="M454 142 a8 5 0 0 1 16 0Z" fill="#ffd27a"/>
${pl("bunny", 352, 160, 36, "sip")}${pl("cat", 424, 162, 36, "sip2")}
<ellipse cx="406" cy="204" rx="30" ry="5" fill="#ffd6c9" stroke="#e3a083" stroke-width="1.2"/><rect x="403" y="206" width="6" height="14" fill="#e3a083"/>
<rect x="390" y="194" width="9" height="9" rx="2" fill="#fff" stroke="#e08aa4"/><rect x="413" y="194" width="9" height="9" rx="2" fill="#ffb3d1" stroke="#e08aa4"/>
<path class="steam" d="M394 190 q-3 -5 0 -9" stroke="#e8a0b8" stroke-width="1.4" fill="none"/><path class="steam d" d="M417 190 q-3 -5 0 -9" stroke="#e8a0b8" stroke-width="1.4" fill="none"/>
<rect x="48" y="238" width="132" height="12" rx="4" fill="#c9a98e"/><rect x="340" y="238" width="132" height="12" rx="4" fill="#c9a98e"/>
${flowers}
<rect x="26" y="250" width="468" height="44" rx="6" fill="#fff0f4"/>
<rect x="26" y="250" width="468" height="5" fill="#f7c4d4"/>
<a class="cafe-link" href="#/test/${today.slug}"><title>Test of the day: ${today.label}</title>
<path d="M332 296 L344 262 L356 296" fill="none" stroke="#c9a98e" stroke-width="2.5"/>
<rect x="312" y="258" width="64" height="30" rx="4" fill="#2e2550" stroke="#c9a98e" stroke-width="2"/>
<text x="344" y="270" text-anchor="middle" font-size="6.5" fill="#ffd27a">today's test</text>
<text x="344" y="281" text-anchor="middle" font-size="7" fill="#f6f3ff">${today.label}</text></a>
<ellipse cx="438" cy="276" rx="20" ry="4" fill="#fff" stroke="#e08aa4" stroke-width="1.2"/><rect x="436" y="278" width="4" height="16" fill="#e08aa4"/>
<path d="M416 252 q22 -20 44 0Z" fill="#ffb3d1" stroke="#c2477a" stroke-width="1"/><rect x="437" y="250" width="2" height="26" fill="#e08aa4"/>
${pl("penguin", 450, 258, 30, "bob2")}
<rect x="430" y="270" width="7" height="6" rx="1.5" fill="#fff" stroke="#e08aa4" stroke-width=".8"/>
<clipPath id="cafeWalk"><rect x="26" y="190" width="468" height="112"/></clipPath>
<g clip-path="url(#cafeWalk)">${cust}${run}</g>
</g>
</svg>`;
  }

  window.CAFE = { full };
})();
