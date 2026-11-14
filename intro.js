// Intro: the three circles of the logo appear large and covered in red dots,
// like Yayoi Kusama's dotted pumpkins; the dots fade, the circles take their
// colours and shrink into the logo, and the name appears. 1.7 s.
// It plays once per tab (a reload in the same tab skips it), skips on a click
// or any key, and does not play with reduced motion or in a background tab.
// It is loaded before the rest of the page, so the page never flashes first;
// if it fails, the page shows as usual.
(() => {
  try {
    if (sessionStorage.getItem("piste-intro")) return;
    sessionStorage.setItem("piste-intro", "1");
  } catch { /* storage blocked: play anyway */ }
  if (matchMedia("(prefers-reduced-motion: reduce)").matches || document.visibilityState === "hidden") return;

  const BG = "#fbfcfe", RED = "#d23b2b";
  // Logo circles around the logo centre, in the units of logo.svg (viewBox 28)
  const LOGO = [
    { x: -5, y: -5.5, c: "#5368a6", a: 1 },
    { x: 5, y: -2.5, c: RED, a: 0.85 },
    { x: -2, y: 5.5, c: "#f7e0dd", a: 1 },
  ];
  const R = 5;                  // circle radius
  const END = 3.4;              // final scale: logo about 95 px wide
  const TOTAL = 1700, FADE = 250, SKIP_FADE = 200;
  const WORD = 880;             // the name starts to appear

  const clamp = (x) => Math.max(0, Math.min(1, x));
  const seg = (t, a, b) => clamp((t - a) / (b - a));
  const outCubic = (x) => 1 - (1 - x) ** 3;
  const inCubic = (x) => x ** 3;
  const inOutCubic = (x) => (x < 0.5 ? 4 * x ** 3 : 1 - (-2 * x + 2) ** 3 / 2);
  const outBack = (x, s) => (x <= 0 ? 0 : 1 + (s + 1) * (x - 1) ** 3 + s * (x - 1) ** 2);
  const hash = (i) => { const x = Math.sin(i * 12.9898 + 78.233) * 43758.5453; return x - Math.floor(x); };
  const rgb = (h) => [1, 3, 5].map((j) => parseInt(h.slice(j, j + 2), 16));
  const mix = (a, b, k) => { const A = rgb(a), B = rgb(b); return `rgb(${A.map((v, j) => Math.round(v + (B[j] - v) * k)).join(",")})`; };

  const W = innerWidth, H = innerHeight, dpr = Math.min(devicePixelRatio || 1, 2);
  const C = { x: W / 2, y: H / 2 - 26 };        // logo centre; the name sits below
  const S0 = Math.min(W * 0.56, H * 0.52) / 20; // scale of the large circles

  // Dots on a hexagonal grid in each circle, larger in the middle and smaller
  // towards the edge, so each circle reads as a sphere
  const Rb = R * S0, d = Rb / 3.8;
  const discs = LOGO.map((c, i) => {
    const dots = [];
    let k = 0;
    for (let y = -Rb; y <= Rb; y += d * 0.866) {
      const row = Math.round((y + Rb) / (d * 0.866));
      for (let x = -Rb - d; x <= Rb + d; x += d) {
        const px = x + (row % 2 ? d / 2 : 0), q = Math.hypot(px, y) / Rb;
        if (q > 1 - (d * 0.3) / Rb) continue;
        dots.push({
          x: px / S0, y: y / S0,
          size: (d * 0.36 * Math.sqrt(1 - q * q) ** 0.6 + d * 0.06) / S0,
          appear: 65 + i * 85 + q * 265,
          leave: 630 + 150 * hash(i * 1000 + k++),
        });
      }
    }
    return { ...c, i, dots };
  });

  const el = document.createElement("div");
  el.className = "intro";
  el.setAttribute("aria-hidden", "true");
  el.innerHTML = `<canvas></canvas><div class="intro-word">${[..."PISTE"].map((l) => `<span>${l}</span>`).join("")}</div>`;
  document.body.prepend(el);
  const canvas = el.querySelector("canvas"), ctx = canvas.getContext("2d");
  canvas.width = W * dpr; canvas.height = H * dpr;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  const word = el.querySelector(".intro-word"), letters = [...word.children];
  word.style.top = `${C.y + 10.5 * END + 22}px`;

  function render(t) {
    ctx.clearRect(0, 0, W, H);
    const m = inOutCubic(seg(t, 680, 1100));          // large circles shrink into the logo
    for (const c of discs) {
      const p = outCubic(seg(t, 35 + c.i * 85, 435 + c.i * 85));
      if (p <= 0) continue;
      const s = S0 * (0.72 + 0.28 * p) * (1 - m) + END * m;
      const cx = C.x + c.x * s, cy = C.y + c.y * s;
      // Circle: a pale tint of its colour that fills in as it shrinks
      ctx.globalAlpha = p * (1 + (c.a - 1) * m);
      ctx.fillStyle = mix(BG, c.c, 0.14 + 0.86 * m);
      ctx.beginPath();
      ctx.arc(cx, cy, R * s, 0, 2 * Math.PI);
      ctx.fill();
      // Dots
      ctx.globalAlpha = 1;
      ctx.fillStyle = RED;
      ctx.beginPath();
      const breathe = 1 + 0.05 * Math.sin((t / 1000) * 2 * Math.PI * 1.1 + c.i);
      for (const o of c.dots) {
        const a = outBack(seg(t, o.appear, o.appear + 250), 1.8);
        const v = 1 - inCubic(seg(t, o.leave, o.leave + 250));
        const r = o.size * s * Math.min(a, 1.2) * v * breathe;
        if (r <= 0.05) continue;
        ctx.moveTo(cx + o.x * s + r, cy + o.y * s);
        ctx.arc(cx + o.x * s, cy + o.y * s, r, 0, 2 * Math.PI);
      }
      ctx.fill();
    }
    ctx.globalAlpha = 1;
    // The name: letters come into focus one by one while the spacing closes
    const ls = `${0.4 + 0.35 * (1 - outCubic(seg(t, WORD, WORD + 500)))}em`;
    word.style.letterSpacing = ls;
    word.style.paddingLeft = ls;
    letters.forEach((l, i) => {
      const p = outCubic(seg(t, WORD + i * 40, WORD + i * 40 + 300));
      l.style.opacity = p;
      l.style.transform = `translateY(${(1 - p) * 10}px)`;
      l.style.filter = p < 1 ? `blur(${(1 - p) * 5}px)` : "none";
    });
  }

  let t0 = null, raf = 0, fading = false;
  const end = (ms) => {
    if (fading) return;
    fading = true;
    removeEventListener("keydown", skip, true);
    el.style.transitionDuration = `${ms}ms`;
    el.classList.add("out");
    setTimeout(() => { cancelAnimationFrame(raf); el.remove(); }, ms + 50);
  };
  const skip = () => end(SKIP_FADE);
  el.addEventListener("pointerdown", skip);
  addEventListener("keydown", skip, true);
  const frame = (now) => {
    if (t0 === null) t0 = now;
    const t = now - t0;
    render(Math.min(t, TOTAL));
    if (t >= TOTAL - FADE) end(FADE);
    raf = requestAnimationFrame(frame);
  };
  raf = requestAnimationFrame(frame);
})();
