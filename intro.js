// Intro: the spots of a section appear, a region is drawn around some of them
// and they turn red, as in the app. It plays once per tab (a reload in the same
// tab skips it), lasts at most 1.5 s, ends as soon as the page is ready after
// the short story, skips on a click or a key, and does not play with reduced
// motion. It is loaded before the rest of the page, so the page never flashes
// first; if it fails, the page shows as usual.
(() => {
  try {
    if (sessionStorage.getItem("piste-intro")) return;
    sessionStorage.setItem("piste-intro", "1");
  } catch { /* storage blocked: play anyway */ }
  if (matchMedia("(prefers-reduced-motion: reduce)").matches || document.visibilityState === "hidden") return;

  const STORY = 1000;      // ms until the region is red
  const LATEST = 1250;     // fade starts by then at the latest
  const FADE = 250;        // so the intro is gone within 1.5 s
  const SKIP_FADE = 150;

  // Hexagonal spot grid, like a Visium section
  const COLS = 12, ROWS = 10, DX = 16, DY = 14, R = 5, PAD = 16;
  const W = 2 * PAD + (COLS - 1) * DX + DX / 2, H = 2 * PAD + (ROWS - 1) * DY;
  const cx = W / 2, cy = H / 2;

  // Tissue outline: a slightly irregular ellipse
  const inTissue = (x, y) => {
    const a = Math.atan2(y - cy, x - cx);
    const k = ((x - cx) / (W * 0.47)) ** 2 + ((y - cy) / (H * 0.46)) ** 2;
    return k <= 1 + 0.1 * Math.sin(3 * a + 0.6);
  };

  // Region drawn with the lasso: a wobbly loop in the lower left
  const rx = cx - W * 0.17, ry = cy + H * 0.13, ra = W * 0.2, rb = H * 0.25;
  const lasso = [];
  for (let i = 0; i < 28; i++) {
    const t = (i / 28) * 2 * Math.PI;
    const w = 1 + 0.08 * Math.sin(3 * t + 1);
    lasso.push([rx + Math.cos(t) * ra * w, ry + Math.sin(t) * rb * w]);
  }
  const inside = (px, py) => {
    let inn = false;
    for (let i = 0, j = lasso.length - 1; i < lasso.length; j = i++) {
      const [xi, yi] = lasso[i], [xj, yj] = lasso[j];
      if ((yi > py) !== (yj > py) && px < ((xj - xi) * (py - yi)) / (yj - yi) + xi) inn = !inn;
    }
    return inn;
  };

  let spots = "";
  for (let row = 0; row < ROWS; row++) {
    for (let col = 0; col < COLS; col++) {
      const x = PAD + col * DX + (row % 2 ? DX / 2 : 0), y = PAD + row * DY;
      if (!inTissue(x, y)) continue;
      // Appear in a wave from the top left, within 0.3 s
      const d = Math.round(300 * (0.6 * row / (ROWS - 1) + 0.4 * col / (COLS - 1)));
      // Spots in the region turn red from its centre outwards, 0.65 to 0.75 s
      const attrs = inside(x, y)
        ? `class="spot in" style="--d:${d}ms;--r:${650 + Math.round(100 * Math.min(1, Math.hypot((x - rx) / ra, (y - ry) / rb)))}ms"`
        : `class="spot" style="--d:${d}ms"`;
      spots += `<circle cx="${x.toFixed(1)}" cy="${y.toFixed(1)}" r="${R}" ${attrs}/>`;
    }
  }
  const path = "M" + lasso.map((p) => p.map((v) => v.toFixed(1)).join(" ")).join("L") + "Z";

  const el = document.createElement("div");
  el.className = "intro";
  el.setAttribute("aria-hidden", "true");
  el.innerHTML = `<svg viewBox="0 0 ${W} ${H}">${spots}<path class="lasso" d="${path}" pathLength="1"/></svg>`;
  document.body.prepend(el);

  const t0 = performance.now();
  let fading = false;
  const end = (ms) => {
    if (fading) return;
    fading = true;
    removeEventListener("keydown", skip, true);
    el.style.transitionDuration = `${ms}ms`;
    el.classList.add("out");
    setTimeout(() => el.remove(), ms + 50);
  };
  const skip = () => end(SKIP_FADE);
  el.addEventListener("pointerdown", skip);
  addEventListener("keydown", skip, true);

  // Start on the next frame, so the transitions run from the hidden state
  requestAnimationFrame(() => requestAnimationFrame(() => el.classList.add("play")));

  // Fade when the page is ready, but not before the story is told, and
  // never later than LATEST
  const loaded = new Promise((res) => (document.readyState === "complete" ? res() : addEventListener("load", res, { once: true })));
  const fonts = document.fonts ? document.fonts.ready : Promise.resolve();
  Promise.all([loaded, fonts]).then(() => {
    setTimeout(() => end(FADE), Math.max(0, STORY - (performance.now() - t0)));
  });
  setTimeout(() => end(FADE), LATEST);
})();
