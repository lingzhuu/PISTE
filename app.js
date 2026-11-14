// PISTE regional comparison in the browser.
// R runs in this tab through webR; the PISTE R code is loaded from the
// package source on GitHub, so the site always uses the released functions.

// webR from its own site, with the npm copy on jsDelivr as a fallback
const WEBR_SOURCES = [
  { module: "https://webr.r-wasm.org/v0.6.1/webr.mjs", base: "https://webr.r-wasm.org/v0.6.1/" },
  { module: "https://cdn.jsdelivr.net/npm/webr@0.6.1/dist/webr.js", base: "https://cdn.jsdelivr.net/npm/webr@0.6.1/dist/" },
];

const PISTE_R = "https://raw.githubusercontent.com/lingzhuu/PISTE/main/src/piste/";
const R_FILES = ["R/pi_engine.R", "R/inla_screen.R", "R/observed_fit.R", "R/compare.R", "R/cli_io.R"];

const $ = (id) => document.getElementById(id);
const state = {
  webR: null, ready: false, rFailed: false,
  coords: null,          // {ids, x, y}
  image: null,           // {img, scale}
  expr: null,            // {file, delim, header, rowNames}
  region: null,          // Uint8Array, 1 = region
  tool: "add",
  sets: null,            // {name: genes}
  results: [], filter: "all", sort: { key: "pi", desc: true },
};

// Colours come from the stylesheet tokens, so the plot matches the legend
const css = getComputedStyle(document.documentElement);
const tok = (name, fallback) => css.getPropertyValue(name).trim() || fallback;
const COLOR = {
  hi: tok("--data-hi", "#d23b2b"), lo: tok("--data-lo", "#5368a6"),
  rest: tok("--data-rest", "#c3cbdd"), none: tok("--data-none", "#aab5cf"),
};
const rgba = (hex, a) => {
  const n = parseInt(hex.slice(1), 16);
  return `rgba(${n >> 16},${(n >> 8) & 255},${n & 255},${a})`;
};

// A message the user can act on; anything else is reported as unexpected
class UserError extends Error {}
const userText = (e, fallback) => {
  if (e instanceof UserError) return e.message;
  console.error(e);
  return fallback;
};

// Status line: kind is "", "ok", "warn" or "err"; the stylesheet adds the icon
function setMsg(el, text, kind = "") {
  if (el.dataset.base === undefined) el.dataset.base = el.className.replace(/\b(msg|ok|warn|err)\b/g, "").trim();
  el.className = [el.dataset.base, "msg", kind].filter(Boolean).join(" ");
  el.textContent = text;
}

// ---------------------------------------------------------------------------
// R engine
// ---------------------------------------------------------------------------
function setEngine(kind, text, detail = "") {
  $("engine").className = `engine ${kind}`.trim();
  $("engine-text").textContent = text;
  $("engine-actions").hidden = !(kind === "error" || kind === "warn");
  $("engine-detail").hidden = !detail;
  $("engine-detail-text").textContent = detail;
}
$("btn-reload").addEventListener("click", () => location.reload());

async function startR() {
  try {
    setEngine("", "Downloading R (about 30 MB)…");
    let WebR, base, errors = [];
    for (const src of WEBR_SOURCES) {
      try { ({ WebR } = await import(src.module)); base = src.base; break; }
      catch (e) { errors.push(`${src.module}: ${e.message}`); }
    }
    if (!WebR) throw new Error("could not load webR. " + errors.join(" | "));
    const webR = new WebR({ baseUrl: base });
    const slow = setTimeout(() => setEngine("warn", "R is slow to start. Keep waiting, or reload the page.",
      "R downloads from webr.r-wasm.org. A slow network or a firewall that blocks this site delays it."), 60000);
    await webR.init();
    clearTimeout(slow);
    setEngine("", "Loading PISTE…");
    await webR.FS.mkdir("/piste");
    await webR.FS.mkdir("/piste/lib");
    await webR.FS.mkdir("/piste/lib/PISTE");
    const desc = await fetchText(PISTE_R + "DESCRIPTION");
    await webR.FS.writeFile("/piste/lib/PISTE/DESCRIPTION", new TextEncoder().encode(desc));
    for (const f of R_FILES) {
      const code = await fetchText(PISTE_R + f);
      await webR.FS.writeFile("/piste/" + f.split("/")[1], new TextEncoder().encode(code));
    }
    await webR.evalRVoid(`
      .libPaths(c("/piste/lib", .libPaths()))
      for (f in c(${R_FILES.map((f) => `"/piste/${f.split("/")[1]}"`).join(",")})) source(f)
      piste_one <- function(i, min_genes) {
        set <- gene_sets[i]
        t0 <- proc.time()[["elapsed"]]
        res <- compare(Y, region, set, baseline = "Rest", min_genes = min_genes, verbose = FALSE)
        st <- res$gene_sets
        secs <- proc.time()[["elapsed"]] - t0
        if (st$status != "done") return(c(if (st$status == "skipped") 0 else -1, st$genes_used, 0, 0, 0, 0, secs))
        o <- res$overall[res$overall$region == "Region", ]
        call <- switch(o$call, higher = 1, lower = -1, 0)
        c(1, st$genes_used, o$score, o$evidence_higher, o$evidence_lower, call, secs)
      }
      piste_error <- function(i, min_genes) {
        compare(Y, region, gene_sets[i], baseline = "Rest", min_genes = min_genes, verbose = FALSE)$gene_sets$message
      }
    `);
    const version = await webR.evalRString(`as.character(utils::packageVersion("PISTE"))`);
    state.webR = webR;
    state.ready = true;
    setEngine("ready", `R ready · PISTE ${version}`);
  } catch (e) {
    console.error(e);
    state.rFailed = true;
    setEngine("error", "R did not start. Check your internet connection, then reload the page.",
      e && e.message ? e.message : String(e));
  }
  updateRun();
}

async function fetchText(url) {
  const r = await fetch(url);
  if (!r.ok) throw new Error(`${url} (${r.status})`);
  return r.text();
}

// ---------------------------------------------------------------------------
// File reading (streamed line by line; .gz decompressed in the browser)
// ---------------------------------------------------------------------------
// gzip is recognised by its magic bytes, not by the file name (a server may
// already have decompressed a .gz download)
async function textStream(file) {
  const head = new Uint8Array(await file.slice(0, 2).arrayBuffer());
  let s = file.stream();
  if (head[0] === 0x1f && head[1] === 0x8b) s = s.pipeThrough(new DecompressionStream("gzip"));
  return s.pipeThrough(new TextDecoderStream());
}

async function eachLine(file, fn) {
  const reader = (await textStream(file)).getReader();
  let rest = "", n = 0;
  for (;;) {
    const { value, done } = await reader.read();
    if (done) break;
    const parts = (rest + value).split("\n");
    rest = parts.pop();
    for (const line of parts) if (fn(line.replace(/\r$/, ""), n++) === false) { reader.cancel(); return; }
  }
  if (rest.length) fn(rest.replace(/\r$/, ""), n);
}

const guessDelim = (line) => (line.split("\t").length > line.split(",").length ? "\t" : ",");
const unquote = (s) => s.replace(/^"(.*)"$/, "$1");
const split = (line, d) => line.split(d).map(unquote);

async function readCoords(file) {
  const lines = [];
  await eachLine(file, (l) => { if (l.trim()) lines.push(l); });
  if (lines.length < 2) throw new UserError("This file has no spots. Choose the coordinate file of the section.");
  const d = guessDelim(lines[0]);
  let header = split(lines[0], d);
  const rows = lines.slice(1).map((l) => split(l, d));
  if (header.length === rows[0].length - 1) header = ["spot", ...header];   // R write.csv row names
  const cols = header.map((_, j) => rows.map((r) => r[j]));
  const numeric = cols.map((c) => c.every((v) => v !== "" && v !== undefined && !isNaN(+v)));
  const lower = header.map((h) => h.toLowerCase());
  const pick = (names) => names.map((n) => lower.indexOf(n)).find((j) => j > 0 && numeric[j]);
  // Visium tissue_positions: pixel columns, rows go down the image
  let xi = pick(["x", "pxl_col_in_fullres", "imagecol", "sdimx", "x_coord", "col", "array_col"]);
  let yi = pick(["y", "pxl_row_in_fullres", "imagerow", "sdimy", "y_coord", "row", "array_row"]);
  if (xi === undefined || yi === undefined) {
    const num = numeric.map((v, j) => (v && j > 0 ? j : -1)).filter((j) => j > 0);
    if (num.length < 2) throw new UserError("No x and y columns found. Choose a table with spot name, x and y.");
    [xi, yi] = num.slice(-2);
  }
  const keep = lower.indexOf("in_tissue") > 0 ? cols[lower.indexOf("in_tissue")].map((v) => v === "1") : null;
  const idx = rows.map((_, i) => i).filter((i) => !keep || keep[i]);
  return {
    ids: idx.map((i) => cols[0][i]),
    x: Float64Array.from(idx, (i) => +cols[xi][i]),
    y: Float64Array.from(idx, (i) => +cols[yi][i]),
    xName: header[xi], yName: header[yi],
  };
}

async function scanExpression(file) {
  let delim, header, rowNames = [];
  await eachLine(file, (l, i) => {
    if (i === 0) { delim = guessDelim(l); header = split(l, delim); return; }
    if (!l) return;
    const j = l.indexOf(delim);
    rowNames.push(unquote(l.slice(0, j)));
  });
  if (!header || header.length < 2 || !rowNames.length) {
    throw new UserError("This is not a table with a header row. Choose a CSV or TSV file.");
  }
  return { file, delim, rowNames, fullHeader: header };
}

// Which side of the expression table holds the spots of the coordinate file;
// cached until either file changes
let matchCache = null;
function spotMatch(expr, coords) {
  if (matchCache && matchCache.expr === expr && matchCache.coords === coords) return matchCache.m;
  const ids = new Set(coords.ids);
  const inHeader = expr.fullHeader.slice(1).filter((h) => ids.has(h)).length;
  const inRows = expr.rowNames.filter((r) => ids.has(r)).length;
  const m = inHeader > inRows
    ? { spots: inHeader, genes: expr.rowNames.length, genesInRows: true }
    : { spots: inRows, genes: expr.fullHeader.length - 1, genesInRows: false };
  matchCache = { expr, coords, m };
  return m;
}
const NO_MATCH = "Spot names do not match the coordinate file. Choose files from the same section.";

// Spots x genes matrix of the genes in `genes`, for the spots of the section
async function readMatrix(expr, coords, genes, stopped) {
  const spotIndex = new Map(coords.ids.map((s, i) => [s, i]));
  const head = expr.fullHeader.slice(1);
  const match = spotMatch(expr, coords);
  if (match.spots === 0) throw new UserError(NO_MATCH);
  const genesInRows = match.genesInRows;
  const wanted = new Set(genes);
  const n = coords.ids.length;
  const found = [];          // gene names kept
  const columns = [];        // Float64Array per gene, indexed by coordinate spot
  const present = new Uint8Array(n);

  if (genesInRows) {
    const colSpot = head.map((h) => (spotIndex.has(h) ? spotIndex.get(h) : -1));
    colSpot.forEach((s) => { if (s >= 0) present[s] = 1; });
    await eachLine(expr.file, (l, i) => {
      if (stopped()) return false;
      if (i === 0 || !l) return;
      const f = split(l, expr.delim);
      if (!wanted.has(f[0])) return;
      const v = new Float64Array(n);
      for (let j = 1; j < f.length; j++) if (colSpot[j - 1] >= 0) v[colSpot[j - 1]] = +f[j];
      found.push(f[0]); columns.push(v); wanted.delete(f[0]);
    });
  } else {
    const keepCols = [];
    head.forEach((g, j) => { if (wanted.has(g)) { keepCols.push(j + 1); found.push(g); wanted.delete(g); } });
    keepCols.forEach(() => columns.push(new Float64Array(n)));
    await eachLine(expr.file, (l, i) => {
      if (stopped()) return false;
      if (i === 0 || !l) return;
      const f = split(l, expr.delim);
      const s = spotIndex.get(f[0]);
      if (s === undefined) return;
      present[s] = 1;
      for (let k = 0; k < keepCols.length; k++) columns[k][s] = +f[keepCols[k]];
    });
  }
  const spots = [];
  for (let i = 0; i < n; i++) if (present[i]) spots.push(i);
  const values = new Float64Array(spots.length * found.length);
  columns.forEach((c, k) => spots.forEach((s, r) => { values[k * spots.length + r] = c[s]; }));
  if (values.some((v) => !Number.isFinite(v))) {
    throw new UserError("The expression file has empty or non-numeric values. Fix them and load the file again.");
  }
  return { spots, genes: found, values };
}

function parseGmt(text) {
  const sets = {};
  for (const line of text.split(/\r?\n/)) {
    const f = line.split("\t");
    if (f.length < 3) continue;
    sets[f[0]] = [...new Set(f.slice(2).filter(Boolean))];
  }
  return sets;
}

// ---------------------------------------------------------------------------
// Section plot and lasso
// ---------------------------------------------------------------------------
const canvas = $("canvas");
const ctx = canvas.getContext("2d");
let view = null, lasso = null;

function layout() {
  const c = state.coords;
  if (!c) return;
  const box = $("canvas-box").getBoundingClientRect();
  const dpr = window.devicePixelRatio || 1;
  canvas.width = box.width * dpr; canvas.height = box.height * dpr;
  let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
  for (let i = 0; i < c.x.length; i++) {
    x0 = Math.min(x0, c.x[i]); x1 = Math.max(x1, c.x[i]);
    y0 = Math.min(y0, c.y[i]); y1 = Math.max(y1, c.y[i]);
  }
  const pad = 24;
  const k = Math.min((box.width - 2 * pad) / (x1 - x0 || 1), (box.height - 2 * pad) / (y1 - y0 || 1));
  const ox = (box.width - k * (x1 - x0)) / 2 - k * x0;
  const oy = (box.height - k * (y1 - y0)) / 2 - k * y0;
  const area = k * k * (x1 - x0) * (y1 - y0);
  const r = Math.max(1.2, Math.min(9, 0.5 * Math.sqrt(area / c.x.length)));
  view = { k, ox, oy, r, dpr };
}

function draw() {
  const c = state.coords;
  if (!c || !view) return;
  const { k, ox, oy, r, dpr } = view;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  const im = state.image;
  if (im) ctx.drawImage(im.img, ox, oy, (k * im.img.width) / im.scale, (k * im.img.height) / im.scale);
  for (const pass of [0, 1]) {
    ctx.fillStyle = pass ? (im ? rgba(COLOR.hi, 0.72) : COLOR.hi) : (im ? rgba(COLOR.rest, 0.6) : COLOR.rest);
    ctx.beginPath();
    for (let i = 0; i < c.x.length; i++) {
      if (state.region[i] !== pass) continue;
      const px = ox + k * c.x[i], py = oy + k * c.y[i];
      ctx.moveTo(px + r, py);
      ctx.arc(px, py, r, 0, 2 * Math.PI);
    }
    ctx.fill();
  }
  if (lasso && lasso.length > 1) {
    const col = state.tool === "add" ? COLOR.hi : COLOR.lo;
    ctx.beginPath();
    ctx.moveTo(lasso[0][0], lasso[0][1]);
    for (const p of lasso) ctx.lineTo(p[0], p[1]);
    ctx.closePath();
    ctx.fillStyle = rgba(col, 0.08);
    ctx.fill();
    ctx.setLineDash([5, 4]);
    ctx.lineWidth = 1.5;
    ctx.strokeStyle = col;
    ctx.stroke();
    ctx.setLineDash([]);
  }
}

function inside(px, py, poly) {
  let inn = false;
  for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {
    const [xi, yi] = poly[i], [xj, yj] = poly[j];
    if ((yi > py) !== (yj > py) && px < ((xj - xi) * (py - yi)) / (yj - yi) + xi) inn = !inn;
  }
  return inn;
}

const pos = (e) => { const b = canvas.getBoundingClientRect(); return [e.clientX - b.left, e.clientY - b.top]; };
canvas.addEventListener("pointerdown", (e) => {
  if (!state.coords) return;
  canvas.setPointerCapture(e.pointerId);
  lasso = [pos(e)];
});
canvas.addEventListener("pointermove", (e) => { if (lasso) { lasso.push(pos(e)); draw(); } });
canvas.addEventListener("pointerup", () => {
  if (!lasso) return;
  if (lasso.length > 2) {
    const { k, ox, oy } = view, c = state.coords, v = state.tool === "add" ? 1 : 0;
    for (let i = 0; i < c.x.length; i++) if (inside(ox + k * c.x[i], oy + k * c.y[i], lasso)) state.region[i] = v;
  }
  lasso = null;
  regionChanged();
});

document.querySelectorAll(".toolbar .seg-btn").forEach((b) => b.addEventListener("click", () => {
  document.querySelectorAll(".toolbar .seg-btn").forEach((x) => {
    x.classList.toggle("active", x === b);
    x.setAttribute("aria-pressed", String(x === b));
  });
  state.tool = b.dataset.tool;
}));
$("btn-clear").addEventListener("click", () => { if (state.region) { state.region.fill(0); regionChanged(); } });
function regionChanged() {
  draw();
  const n = state.region.reduce((a, b) => a + b, 0);
  $("region-count").innerHTML = n
    ? `Region <b>${n.toLocaleString()}</b> spots · Rest ${(state.region.length - n).toLocaleString()} spots`
    : "No region yet";
  markSteps();
  updateRun();
}

window.addEventListener("resize", () => { layout(); draw(); });

// ---------------------------------------------------------------------------
// Inputs
// ---------------------------------------------------------------------------
async function loadCoords(file) {
  const st = $("status-coord");
  setMsg(st, "Reading…");
  try {
    const c = await readCoords(file);
    state.coords = c;
    state.region = new Uint8Array(c.ids.length);
    $("drop-coord").classList.add("loaded");
    setMsg(st, `${file.name} · ${c.ids.length.toLocaleString()} spots (x = ${c.xName}, y = ${c.yName})`, "ok");
    $("canvas-empty").hidden = true;
    layout(); regionChanged();
  } catch (e) {
    setMsg(st, userText(e, "Could not read this file. Choose a CSV or TSV file."), "err");
  }
  showExprStatus();
}

async function loadImage(files) {
  const st = $("status-img");
  setMsg(st, "");
  try {
    const imgFile = files.find((f) => /\.(png|jpe?g)$/i.test(f.name));
    const json = files.find((f) => /\.json$/i.test(f.name));
    if (imgFile) {
      const img = new Image();
      img.src = URL.createObjectURL(imgFile);
      try { await img.decode(); }
      catch { throw new UserError("Could not open this image. Choose a PNG or JPG file."); }
      state.image = { img, scale: +$("img-scale").value || 1, name: imgFile.name };
    }
    if (json && state.image) {
      let sf;
      try { sf = JSON.parse(await json.text()); }
      catch { throw new UserError("Could not read the scale factors. Choose scalefactors_json.json from the Space Ranger spatial folder."); }
      const name = state.image.name.toLowerCase();
      const scale = name.includes("lowres") ? sf.tissue_lowres_scalef : sf.tissue_hires_scalef;
      if (scale) { state.image.scale = scale; $("img-scale").value = scale; }
    }
    if (!state.image) throw new UserError("No image selected. Choose a PNG or JPG file.");
    $("drop-img").classList.add("loaded");
    $("scale-row").hidden = false;
    setMsg(st, `${state.image.name} · ${state.image.img.width} × ${state.image.img.height}`, "ok");
    draw();
  } catch (e) {
    setMsg(st, userText(e, "Could not open this image. Choose a PNG or JPG file."), "err");
  }
}
$("file-img").addEventListener("change", (e) => e.target.files.length && loadImage([...e.target.files]));
$("img-scale").addEventListener("change", () => { if (state.image) { state.image.scale = +$("img-scale").value || 1; draw(); } });
$("img-scale").addEventListener("click", (e) => e.preventDefault());

async function loadExpr(file) {
  const st = $("status-expr");
  setMsg(st, "Reading…");
  state.expr = null;
  try {
    state.expr = await scanExpression(file);
    $("drop-expr").classList.add("loaded");
  } catch (e) {
    $("drop-expr").classList.remove("loaded");
    setMsg(st, userText(e, "Could not read this file. Choose a CSV or TSV file."), "err");
  }
  showExprStatus();
  markSteps(); updateRun();
}

// Expression status, checked against the coordinates as soon as both are in
function showExprStatus() {
  const ex = state.expr;
  if (!ex) return;
  const st = $("status-expr"), name = ex.file.name;
  if (!state.coords) {
    setMsg(st, `${name} · ${ex.rowNames.length.toLocaleString()} rows × ${(ex.fullHeader.length - 1).toLocaleString()} columns`, "ok");
    return;
  }
  const m = spotMatch(ex, state.coords);
  if (m.spots === 0) setMsg(st, NO_MATCH, "err");
  else setMsg(st, `${name} · ${m.spots.toLocaleString()} spots, ${m.genes.toLocaleString()} genes`, "ok");
  updateRun();
}

$("file-coord").addEventListener("change", (e) => e.target.files[0] && loadCoords(e.target.files[0]));
$("file-expr").addEventListener("change", (e) => e.target.files[0] && loadExpr(e.target.files[0]));
for (const id of ["drop-expr", "drop-coord", "drop-img"]) {
  const el = $(id);
  el.addEventListener("dragover", (e) => { e.preventDefault(); el.classList.add("over"); });
  el.addEventListener("dragleave", () => el.classList.remove("over"));
  el.addEventListener("drop", (e) => {
    e.preventDefault(); el.classList.remove("over");
    const f = e.dataTransfer.files[0];
    if (!f) return;
    if (id === "drop-img") loadImage([...e.dataTransfer.files]);
    else (id === "drop-expr" ? loadExpr : loadCoords)(f);
  });
}

const EXAMPLE_GENES = "MBP PLP1 MOBP MAG MOG CNP MAL CLDN11 ERMN OPALIN TF CNTN2 ENPP6 UGT8 GJB1 PLLP ASPA FA2H";
$("btn-example").addEventListener("click", async () => {
  // Same site first, then the copy of the website branch on jsDelivr
  const get = async (name) => {
    let last;
    for (const url of ["data/" + name, "https://cdn.jsdelivr.net/gh/lingzhuu/PISTE@website/data/" + name]) {
      try {
        const r = await fetch(url);
        if (r.ok) return new File([await r.blob()], name);
        last = new Error(`${url} (${r.status})`);
      } catch (e) { last = new Error(`${url}: ${e.message}`); }
    }
    throw last;
  };
  $("btn-example").disabled = true;
  try {
    setMsg($("status-coord"), "Downloading the example…");
    await loadCoords(await get("example_coordinates.csv"));
    setMsg($("status-img"), "Downloading the example…");
    await loadImage([await get("example_tissue_hires.jpg"), await get("example_scalefactors.json")]);
    const wm = new Set((await (await get("example_region.txt")).text()).split(/\s+/));
    state.coords.ids.forEach((id, i) => { state.region[i] = wm.has(id) ? 1 : 0; });
    regionChanged();
    $("tab-text").click();
    $("set-name").value = "Myelination";
    $("set-genes").value = EXAMPLE_GENES;
    $("set-genes").dispatchEvent(new Event("input"));
    setMsg($("status-expr"), "Downloading the example expression (6 MB)…");
    await loadExpr(await get("example_expression.csv"));
  } catch (e) {
    console.error(e);
    setMsg($("status-expr"), "Could not download the example. Check your internet connection and try again.", "err");
  }
  $("btn-example").disabled = false;
});

// Gene sets: typed genes (one set) or a list file (plain list = one set, GMT = several)
const splitGenes = (t) => [...new Set(t.split(/[\s,;]+/).map((g) => g.trim().replace(/^"(.*)"$/, "$1")).filter(Boolean))];
let setSource = "text", fileSets = null;

function parseList(text, fileName) {
  const lines = text.split(/\r?\n/).filter((l) => l.trim());
  if (lines.some((l) => l.split("\t").length >= 3)) return parseGmt(text);
  const name = fileName.replace(/\.(txt|csv|tsv|list|gmt)$/i, "");
  const genes = splitGenes(lines.join("\n"));
  return genes.length ? { [name]: genes } : {};
}

$("file-gmt").addEventListener("change", async (e) => {
  const f = e.target.files[0];
  if (!f) return;
  const sets = parseList(await f.text(), f.name);
  const k = Object.keys(sets);
  const st = $("gmt-name");
  if (!k.length) {
    fileSets = null;
    $("drop-gmt").classList.remove("loaded");
    setMsg(st, "No genes found. Choose a text file with one gene per line, or a GMT file.", "err");
  } else {
    fileSets = sets;
    $("drop-gmt").classList.add("loaded");
    setMsg(st, k.length === 1 ? `${f.name} · ${sets[k[0]].length} genes` : `${f.name} · ${k.length} gene sets`, "ok");
  }
  updateRun();
});
$("set-genes").addEventListener("input", () => {
  const n = splitGenes($("set-genes").value).length;
  $("text-count").textContent = n ? `${n} genes` : "";
  updateRun();
});

// Tabs: click, or arrow keys between them
const tabs = [...document.querySelectorAll(".set-tabs .seg-btn")];
function selectTab(b) {
  tabs.forEach((x) => {
    const on = x === b;
    x.classList.toggle("active", on);
    x.setAttribute("aria-selected", String(on));
    x.tabIndex = on ? 0 : -1;
  });
  setSource = b.dataset.src;
  $("panel-text").hidden = setSource !== "text";
  $("panel-file").hidden = setSource !== "file";
  updateRun();
}
tabs.forEach((b, i) => {
  b.addEventListener("click", () => selectTab(b));
  b.addEventListener("keydown", (e) => {
    if (e.key !== "ArrowRight" && e.key !== "ArrowLeft") return;
    const next = tabs[(i + (e.key === "ArrowRight" ? 1 : tabs.length - 1)) % tabs.length];
    selectTab(next); next.focus();
  });
});

function chosenSets() {
  if (setSource === "file") return fileSets;
  const genes = splitGenes($("set-genes").value);
  return genes.length ? { [$("set-name").value.trim() || "Gene set"]: genes } : null;
}

function updateRun() {
  const n = state.region ? state.region.reduce((a, b) => a + b, 0) : 0;
  let why = "", kind = "";
  if (state.rFailed) { why = "R did not start. Reload the page."; kind = "err"; }
  else if (!state.ready) why = "Waiting for R to start…";
  else if (!state.expr) why = "Load the expression file.";
  else if (!state.coords) why = "Load the spot coordinates.";
  else if (spotMatch(state.expr, state.coords).spots === 0) { why = NO_MATCH; kind = "err"; }
  else if (n < 3) why = "Select a region.";
  else if (n > state.region.length - 3) why = "Leave some spots outside the region.";
  else if (!chosenSets()) why = setSource === "text" ? "Enter the genes of a gene set." : "Upload a gene list.";
  $("btn-run").disabled = !!why || running;
  // A message from the last run (done, cancelled, error) stays until there is
  // something new to say
  if (!running) {
    if (why) { setMsg($("run-msg"), why, kind); runNote = false; }
    else if (!runNote) setMsg($("run-msg"), "");
  }
  markSteps();
}

function markSteps() {
  const n = state.region ? state.region.reduce((a, b) => a + b, 0) : 0;
  const done = { data: !!(state.expr && state.coords), region: n > 0, sets: !!chosenSets(), results: state.results.length > 0 };
  document.querySelectorAll(".step").forEach((s) => s.classList.toggle("done", done[s.dataset.step]));
}

// ---------------------------------------------------------------------------
// Comparison
// ---------------------------------------------------------------------------
let running = false;
let runNote = false;         // the run message shows the outcome of the last run
let stopRun = null;         // resolves the race below when Cancel is pressed
let rBusy = null;            // an R fit still running after a cancel
const STOP = Symbol("stop");
$("btn-run").addEventListener("click", run);
$("btn-cancel").addEventListener("click", () => stopRun && stopRun());

// Progress: the finished part is solid; the part being computed pulses, so
// the bar moves from the start even when there is one gene set
function setProgress(done, total, working) {
  const w = 100 / total;
  $("bar").style.width = `${done * w}%`;
  $("cur").style.left = `${done * w}%`;
  $("cur").style.width = working ? `${w}%` : "0";
  $("progress").setAttribute("aria-valuenow", String(Math.round(done * w)));
}

// An R call that Cancel can abandon; R cannot be interrupted here, so the
// call itself finishes in the background before R takes the next one
function rCall(promise, stop) {
  rBusy = promise.finally(() => { if (rBusy === tracked) rBusy = null; });
  const tracked = rBusy;
  rBusy.catch(() => {});
  return Promise.race([promise, stop]);
}

async function run() {
  running = true;
  let stopped = false;
  const stop = new Promise((res) => { stopRun = () => { stopped = true; res(STOP); }; });
  $("btn-cancel").hidden = false;
  updateRun();
  const msg = $("run-msg");
  $("progress").hidden = false;
  setProgress(0, 1, true);
  const webR = state.webR;
  let done = 0, total = 0;
  try {
    if (rBusy) {
      setMsg(msg, "Waiting for the cancelled fit to finish…");
      if (await Promise.race([rBusy.catch(() => {}), stop]) === STOP) throw STOP;
    }
    const sets = chosenSets();
    const names = Object.keys(sets);
    total = names.length;
    const genes = [...new Set(names.flatMap((s) => sets[s]))];
    setMsg(msg, "Reading expression of the gene-set genes…");
    const m = await readMatrix(state.expr, state.coords, genes, () => stopped);
    if (stopped) throw STOP;
    if (!m.genes.length) throw new UserError("None of the genes are in the expression file. Use gene symbols (such as MBP), not Ensembl IDs.");
    const regionLabels = m.spots.map((s) => (state.region[s] ? "Region" : "Rest"));
    const ids = m.spots.map((s) => state.coords.ids[s]);
    const assigned = await rCall(webR.evalRVoid(
      `assign("Y", matrix(v, nrow = length(ids), dimnames = list(ids, genes)), globalenv())
       assign("region", labels, globalenv())
       assign("gene_sets", split(set_genes, factor(set_index, seq_along(set_names), set_names)), globalenv())`,
      { env: { v: Array.from(m.values), ids, genes: m.genes, labels: regionLabels,
               set_names: names, set_genes: names.flatMap((s) => sets[s]),
               set_index: names.flatMap((s, i) => sets[s].map(() => i + 1)) } }), stop);
    if (assigned === STOP) throw STOP;
    const minGenes = Math.max(2, +$("min-genes").value || 5);
    state.results = [];
    state.lastRun = { minGenes };
    const t0 = performance.now();
    for (let i = 0; i < names.length; i++) {
      setProgress(i, total, true);
      setMsg(msg, `Fitting ${i + 1} of ${total}: ${pretty(names[i])}`);
      const r = await rCall(webR.evalRRaw(`piste_one(${i + 1}L, ${minGenes})`, "number[]"), stop);
      if (r === STOP) throw STOP;
      const row = { name: names[i], size: sets[names[i]].length, genes: r[1], status: r[0], secs: r[6] };
      if (r[0] === 1) Object.assign(row, { pi: r[2], hi: r[3], lo: r[4], call: r[5] === 1 ? "higher" : r[5] === -1 ? "lower" : "none" });
      else if (r[0] === 0) row.call = "skipped";
      else {
        row.call = "failed";
        const why = await rCall(webR.evalRString(`piste_error(${i + 1}L, ${minGenes})`), stop);
        if (why === STOP) throw STOP;
        row.message = why;
      }
      state.results.push(row);
      done = i + 1;
      setProgress(done, total, false);
      if (i % 3 === 2 || i === names.length - 1) render();
      await new Promise((res) => setTimeout(res, 0));
    }
    const secs = (performance.now() - t0) / 1000;
    setMsg(msg, `Done in ${secs < 90 ? secs.toFixed(0) + " s" : (secs / 60).toFixed(1) + " min"} · ${m.spots.length.toLocaleString()} spots, ${m.genes.length.toLocaleString()} genes`, "ok");
    render();
    $("step-results").scrollIntoView();
  } catch (e) {
    if (e === STOP) {
      setMsg(msg, total ? `Cancelled. ${done} of ${total} gene sets done.` : "Cancelled.");
      render();
    } else {
      setMsg(msg, userText(e, "The comparison stopped with an error. Reload the page and try again."), "err");
    }
  }
  $("progress").hidden = true;
  $("btn-cancel").hidden = true;
  stopRun = null;
  running = false; runNote = true; updateRun(); markSteps();
}

// ---------------------------------------------------------------------------
// Results
// ---------------------------------------------------------------------------
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
const fmt = (v) => (v === undefined || v === null || Number.isNaN(v) ? "–" : v >= 0.9995 ? ">0.999" : v.toFixed(3));

// MSigDB-style names (HALLMARK_TNFA_SIGNALING_VIA_NFKB) read as
// "TNFA signaling via NFKB"; names typed by the user are left as they are
const SET_PREFIX = /^(HALLMARK|GOBP|GOCC|GOMF|KEGG_MEDICUS|KEGG|REACTOME|WP|BIOCARTA|PID)_/;
const KEEP_CASE = new Set(["DNA", "RNA", "UV", "TNFA", "NFKB", "KRAS", "MTOR", "AKT", "JAK", "STAT", "TGF", "EMT",
  "ROS", "UPR", "IFN", "ATP", "NADH", "GPCR", "MHC", "TCA", "ECM", "NK", "ER", "MAPK", "ERK", "HIF", "VEGF", "WNT", "MYC"]);
function pretty(s) {
  if (!/^[A-Z0-9_]+$/.test(s) || !s.includes("_")) return s;
  const words = s.replace(SET_PREFIX, "").split("_").filter(Boolean);
  return words.map((w, i) => {
    if (/\d/.test(w) || KEEP_CASE.has(w)) return w;
    const low = w.toLowerCase();
    return i === 0 ? low[0].toUpperCase() + low.slice(1) : low;
  }).join(" ");
}

const group = (r) => (r.status === 1 ? r.call : "skipped");
const evidence = (r) => (r.status === 1 ? Math.max(r.hi, r.lo) : undefined);
const CALL_LABEL = { higher: "Higher", lower: "Lower", none: "No call", skipped: "Skipped", failed: "Failed" };
const CALL_ORDER = { higher: 4, none: 3, lower: 2, skipped: 1, failed: 0 };

function render() {
  const rows = state.results;
  $("results-empty").hidden = rows.length > 0;
  $("results").hidden = rows.length === 0;
  $("btn-download").disabled = rows.length === 0;

  // The counts are also the filters
  const n = (g) => rows.filter((r) => group(r) === g).length;
  const filters = [
    ["all", rows.length, "All gene sets", ""], ["higher", n("higher"), "Higher in region", "hi"],
    ["lower", n("lower"), "Lower in region", "lo"], ["none", n("none"), "No call", ""],
  ];
  if (n("skipped")) filters.push(["skipped", n("skipped"), "Skipped or failed", ""]);
  if (!filters.some((f) => f[0] === state.filter)) state.filter = "all";
  $("summary").innerHTML = filters.map(([key, v, label, cls]) =>
    `<button class="stat ${cls}${state.filter === key ? " active" : ""}" data-filter="${key}" aria-pressed="${state.filter === key}">
      <span class="v">${v}</span><span class="l">${label}</span></button>`).join("");

  const key = state.sort.key, sgn = state.sort.desc ? -1 : 1;
  const val = (r) => (key === "name" ? pretty(r.name).toLowerCase() : key === "ev" ? evidence(r) ?? -Infinity
    : key === "call" ? CALL_ORDER[r.call] : r[key] ?? -Infinity);
  const shown = rows
    .filter((r) => state.filter === "all" || group(r) === state.filter)
    .sort((a, b) => (val(a) > val(b) ? sgn : val(a) < val(b) ? -sgn : 0));
  $("table").querySelector("tbody").innerHTML = shown.map((r) => {
    const bar = r.status === 1 ? (() => {
      const lo = Math.min(r.pi, 0.5), hi = Math.max(r.pi, 0.5);
      const col = r.call === "higher" ? COLOR.hi : r.call === "lower" ? COLOR.lo : COLOR.none;
      return `<div class="pibar"><div class="track"><div class="fill" style="left:${lo * 100}%;width:${(hi - lo) * 100}%;background:${col}"></div><div class="mid"></div></div><span class="val">${r.pi.toFixed(3)}</span></div>`;
    })() : r.status === 0
      ? `<span class="row-note">Below the minimum of ${state.lastRun.minGenes} genes</span>`
      : `<span class="row-note err">Fit failed: ${esc(r.message || "unknown error")}</span>`;
    return `<tr>
      <td class="name" title="${esc(r.name)}">${esc(pretty(r.name))}</td>
      <td class="r col-genes">${r.genes}/${r.size}</td>
      <td>${bar}</td>
      <td class="r col-ev">${fmt(evidence(r))}</td>
      <td><span class="badge ${r.call}">${CALL_LABEL[r.call]}</span></td>
    </tr>`;
  }).join("");
  $("filter-empty").hidden = shown.length > 0 || rows.length === 0;
  $("table").querySelectorAll("th").forEach((th) => {
    if (th.dataset.key === key) th.setAttribute("aria-sort", state.sort.desc ? "descending" : "ascending");
    else th.removeAttribute("aria-sort");
  });
  markSteps();
}

$("table").querySelectorAll("th").forEach((th) => th.querySelector("button").addEventListener("click", () => {
  const key = th.dataset.key;
  state.sort = { key, desc: state.sort.key === key ? !state.sort.desc : key !== "name" };
  render();
}));
$("summary").addEventListener("click", (e) => {
  const b = e.target.closest(".stat");
  if (!b) return;
  state.filter = b.dataset.filter;
  render();
  $("summary").querySelector(`[data-filter="${state.filter}"]`).focus();
});

// The CSV keeps the column names of the R package output
$("btn-download").addEventListener("click", () => {
  const head = "gene_set,genes_in_set,genes_used,score,evidence_higher,evidence_lower,call,status";
  const lines = state.results.map((r) => [r.name, r.size, r.genes, r.pi ?? "", r.hi ?? "", r.lo ?? "",
    r.status === 1 ? r.call : "", r.status === 1 ? "done" : r.call].join(","));
  const blob = new Blob([[head, ...lines].join("\n") + "\n"], { type: "text/csv" });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob); a.download = "piste_region_vs_rest.csv"; a.click();
});

// Sidebar: highlight the step in view
const observer = new IntersectionObserver((entries) => {
  for (const e of entries) if (e.isIntersecting) {
    const id = e.target.id.replace("step-", "");
    document.querySelectorAll(".step").forEach((s) => s.classList.toggle("current", s.dataset.step === id));
  }
}, { rootMargin: "-40% 0px -55% 0px" });
document.querySelectorAll(".card").forEach((c) => observer.observe(c));

render(); markSteps(); updateRun();
window.addEventListener("error", (e) => {
  console.error(e.error || e.message);
  if (!state.ready) setEngine("error", "Something went wrong. Reload the page.", e.message);
});
startR();
