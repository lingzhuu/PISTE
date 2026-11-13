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
const CALL_LEVEL = 0.975;

const $ = (id) => document.getElementById(id);
const state = {
  webR: null, ready: false,
  coords: null,          // {ids, x, y}
  image: null,           // {img, scale}
  expr: null,            // {file, delim, header, rowNames}
  region: null,          // Uint8Array, 1 = region
  tool: "add",
  sets: null,            // {name: genes}
  results: [], filter: "all", sort: { key: "pi", desc: true },
};

// ---------------------------------------------------------------------------
// R engine
// ---------------------------------------------------------------------------
async function startR() {
  const box = $("engine");
  try {
    const step = (t) => { $("engine-text").textContent = t; };
    step("Downloading R (about 30 MB)…");
    let WebR, base, errors = [];
    for (const src of WEBR_SOURCES) {
      try { ({ WebR } = await import(src.module)); base = src.base; break; }
      catch (e) { errors.push(`${src.module}: ${e.message}`); }
    }
    if (!WebR) throw new Error("could not load webR. " + errors.join(" | "));
    const webR = new WebR({ baseUrl: base });
    const slow = setTimeout(() => step("Still starting R… (slow network or blocked site webr.r-wasm.org)"), 60000);
    await webR.init();
    clearTimeout(slow);
    $("engine-text").textContent = "Loading PISTE…";
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
    box.classList.add("ready");
    $("engine-text").textContent = `R ready · PISTE ${version}`;
  } catch (e) {
    console.error(e);
    box.classList.add("error");
    $("engine-text").textContent = "R could not start: " + (e && e.message ? e.message : String(e));
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
  const d = guessDelim(lines[0]);
  let header = split(lines[0], d);
  const rows = lines.slice(1).map((l) => split(l, d));
  if (header.length === rows[0].length - 1) header = ["spot", ...header];   // R write.csv row names
  const cols = header.map((_, j) => rows.map((r) => r[j]));
  const numeric = cols.map((c) => c.every((v) => v !== "" && !isNaN(+v)));
  const lower = header.map((h) => h.toLowerCase());
  const pick = (names) => names.map((n) => lower.indexOf(n)).find((j) => j > 0 && numeric[j]);
  // Visium tissue_positions: pixel columns, rows go down the image
  let xi = pick(["x", "pxl_col_in_fullres", "imagecol", "sdimx", "x_coord", "col", "array_col"]);
  let yi = pick(["y", "pxl_row_in_fullres", "imagerow", "sdimy", "y_coord", "row", "array_row"]);
  if (xi === undefined || yi === undefined) {
    const num = numeric.map((v, j) => (v && j > 0 ? j : -1)).filter((j) => j > 0);
    if (num.length < 2) throw new Error("Need two numeric coordinate columns.");
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
  if (!header || header.length < 2) throw new Error("Expected a table with a header row.");
  return { file, delim, rowNames, fullHeader: header };
}

// Spots x genes matrix of the genes in `genes`, for the spots of the section
async function readMatrix(expr, coords, genes) {
  const spotIndex = new Map(coords.ids.map((s, i) => [s, i]));
  const head = expr.fullHeader.slice(1);
  const inHeader = head.filter((h) => spotIndex.has(h)).length;
  const inRows = expr.rowNames.filter((r) => spotIndex.has(r)).length;
  if (Math.max(inHeader, inRows) === 0) throw new Error("No spot names are shared by the expression and coordinate files.");
  const genesInRows = inHeader > inRows;
  const wanted = new Set(genes);
  const n = coords.ids.length;
  const found = [];          // gene names kept
  const columns = [];        // Float64Array per gene, indexed by coordinate spot
  const present = new Uint8Array(n);

  if (genesInRows) {
    const colSpot = head.map((h) => (spotIndex.has(h) ? spotIndex.get(h) : -1));
    colSpot.forEach((s) => { if (s >= 0) present[s] = 1; });
    await eachLine(expr.file, (l, i) => {
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
  if (values.some((v) => !Number.isFinite(v))) throw new Error("The expression file has missing or non-numeric values.");
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
    ctx.fillStyle = pass ? (im ? "rgba(210,59,43,.72)" : "#d23b2b") : (im ? "rgba(83,104,166,.22)" : "#c3cbdd");
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
    ctx.beginPath();
    ctx.moveTo(lasso[0][0], lasso[0][1]);
    for (const p of lasso) ctx.lineTo(p[0], p[1]);
    ctx.closePath();
    ctx.fillStyle = state.tool === "add" ? "rgba(210,59,43,.08)" : "rgba(83,104,166,.08)";
    ctx.fill();
    ctx.setLineDash([5, 4]);
    ctx.lineWidth = 1.5;
    ctx.strokeStyle = state.tool === "add" ? "#d23b2b" : "#5368a6";
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
  document.querySelectorAll(".toolbar .seg-btn").forEach((x) => x.classList.toggle("active", x === b));
  state.tool = b.dataset.tool;
}));
$("btn-clear").addEventListener("click", () => { if (state.region) { state.region.fill(0); regionChanged(); } });
function regionChanged() {
  draw();
  const n = state.region.reduce((a, b) => a + b, 0);
  $("region-count").innerHTML = n ? `<b>${n}</b> region spots · ${state.region.length - n} rest` : "No region yet";
  markSteps();
  updateRun();
}

window.addEventListener("resize", () => { layout(); draw(); });

// ---------------------------------------------------------------------------
// Inputs
// ---------------------------------------------------------------------------
async function loadCoords(file) {
  const st = $("status-coord");
  st.className = "drop-status"; st.textContent = "Reading…";
  try {
    const c = await readCoords(file);
    state.coords = c;
    state.region = new Uint8Array(c.ids.length);
    $("drop-coord").classList.add("loaded");
    st.textContent = `${file.name} · ${c.ids.length.toLocaleString()} spots (x = ${c.xName}, y = ${c.yName})`;
    $("canvas-empty").hidden = true;
    layout(); regionChanged();
  } catch (e) {
    st.className = "drop-status err"; st.textContent = e.message;
  }
}

async function loadImage(files) {
  const st = $("status-img");
  st.className = "drop-status";
  try {
    const imgFile = files.find((f) => /\.(png|jpe?g)$/i.test(f.name));
    const json = files.find((f) => /\.json$/i.test(f.name));
    if (imgFile) {
      const img = new Image();
      img.src = URL.createObjectURL(imgFile);
      await img.decode();
      state.image = { img, scale: +$("img-scale").value || 1, name: imgFile.name };
    }
    if (json && state.image) {
      const sf = JSON.parse(await json.text());
      const name = state.image.name.toLowerCase();
      const scale = name.includes("lowres") ? sf.tissue_lowres_scalef : sf.tissue_hires_scalef;
      if (scale) { state.image.scale = scale; $("img-scale").value = scale; }
    }
    if (!state.image) throw new Error("Choose a PNG or JPG image.");
    $("drop-img").classList.add("loaded");
    $("scale-row").hidden = false;
    st.textContent = `${state.image.name} · ${state.image.img.width} × ${state.image.img.height}`;
    draw();
  } catch (e) {
    st.className = "drop-status err"; st.textContent = e.message;
  }
}
$("file-img").addEventListener("change", (e) => e.target.files.length && loadImage([...e.target.files]));
$("img-scale").addEventListener("change", () => { if (state.image) { state.image.scale = +$("img-scale").value || 1; draw(); } });
$("img-scale").addEventListener("click", (e) => e.preventDefault());

async function loadExpr(file) {
  const st = $("status-expr");
  st.className = "drop-status"; st.textContent = "Reading…";
  try {
    state.expr = await scanExpression(file);
    $("drop-expr").classList.add("loaded");
    st.textContent = `${file.name} · ${state.expr.rowNames.length.toLocaleString()} rows × ${(state.expr.fullHeader.length - 1).toLocaleString()} columns`;
  } catch (e) {
    st.className = "drop-status err"; st.textContent = e.message;
  }
  markSteps(); updateRun();
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
  const get = async (name) => {
    const r = await fetch("data/" + name);
    if (!r.ok) throw new Error(`data/${name} (${r.status})`);
    return new File([await r.blob()], name);
  };
  $("btn-example").disabled = true;
  try {
  $("status-coord").textContent = "Downloading the example…";
  await loadCoords(await get("example_coordinates.csv"));
  $("status-img").textContent = "Downloading the example…";
  await loadImage([await get("example_tissue_hires.jpg"), await get("example_scalefactors.json")]);
  const wm = new Set((await (await fetch("data/example_region.txt")).text()).split(/\s+/));
  state.coords.ids.forEach((id, i) => { state.region[i] = wm.has(id) ? 1 : 0; });
  regionChanged();
  document.querySelector('.set-tabs [data-src="text"]').click();
  $("set-name").value = "Myelination";
  $("set-genes").value = EXAMPLE_GENES;
  $("set-genes").dispatchEvent(new Event("input"));
  $("status-expr").className = "drop-status";
  $("status-expr").textContent = "Downloading the example expression (19 MB)…";
  await loadExpr(await get("example_expression.csv.gz"));
  } catch (e) {
    console.error(e);
    $("status-expr").className = "drop-status err"; $("status-expr").textContent = "Example failed: " + e.message;
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
  return { [name]: splitGenes(lines.join("\n")) };
}

$("file-gmt").addEventListener("change", async (e) => {
  const f = e.target.files[0];
  if (!f) return;
  fileSets = parseList(await f.text(), f.name);
  const k = Object.keys(fileSets);
  $("drop-gmt").classList.add("loaded");
  $("gmt-name").textContent = k.length === 1 ? `${f.name} · ${fileSets[k[0]].length} genes`
    : `${f.name} · ${k.length} gene sets`;
  updateRun();
});
$("set-genes").addEventListener("input", () => {
  const n = splitGenes($("set-genes").value).length;
  $("text-count").textContent = n ? `${n} genes` : "";
  updateRun();
});
document.querySelectorAll(".set-tabs .seg-btn").forEach((b) => b.addEventListener("click", () => {
  document.querySelectorAll(".set-tabs .seg-btn").forEach((x) => x.classList.toggle("active", x === b));
  setSource = b.dataset.src;
  $("panel-text").hidden = setSource !== "text";
  $("panel-file").hidden = setSource !== "file";
  updateRun();
}));

function chosenSets() {
  if (setSource === "file") return fileSets;
  const genes = splitGenes($("set-genes").value);
  return genes.length ? { [$("set-name").value.trim() || "Gene set"]: genes } : null;
}

function updateRun() {
  const n = state.region ? state.region.reduce((a, b) => a + b, 0) : 0;
  let why = "";
  if (!state.ready) why = "Waiting for R to start…";
  else if (!state.expr) why = "Load the expression file.";
  else if (!state.coords) why = "Load the spot coordinates.";
  else if (n < 3) why = "Select a region.";
  else if (n > state.region.length - 3) why = "Leave some spots outside the region.";
  else if (!chosenSets()) why = setSource === "text" ? "Enter the genes of a gene set." : "Upload a gene list.";
  $("btn-run").disabled = !!why || running;
  if (!running && (why || !state.results.length)) { $("run-msg").className = "run-msg"; $("run-msg").textContent = why; }
}

function markSteps() {
  const n = state.region ? state.region.reduce((a, b) => a + b, 0) : 0;
  const done = { data: !!(state.expr && state.coords), region: n > 0, sets: true, results: state.results.length > 0 };
  document.querySelectorAll(".step").forEach((s) => s.classList.toggle("done", done[s.dataset.step]));
}

// ---------------------------------------------------------------------------
// Comparison
// ---------------------------------------------------------------------------
let running = false;
$("btn-run").addEventListener("click", run);

async function run() {
  running = true; updateRun();
  const msg = $("run-msg"), bar = $("bar");
  msg.className = "run-msg";
  $("progress").hidden = false; bar.style.width = "0";
  const webR = state.webR;
  try {
    const sets = chosenSets();
    const names = Object.keys(sets);
    const genes = [...new Set(names.flatMap((s) => sets[s]))];
    msg.textContent = "Reading expression of the gene-set genes…";
    const m = await readMatrix(state.expr, state.coords, genes);
    if (!m.genes.length) throw new Error("None of the gene-set genes are in the expression file. Check that genes are named by symbol.");
    const regionLabels = m.spots.map((s) => (state.region[s] ? "Region" : "Rest"));
    const ids = m.spots.map((s) => state.coords.ids[s]);
    await webR.evalRVoid(
      `assign("Y", matrix(v, nrow = length(ids), dimnames = list(ids, genes)), globalenv())
       assign("region", labels, globalenv())
       assign("gene_sets", split(set_genes, factor(set_index, seq_along(set_names), set_names)), globalenv())`,
      { env: { v: Array.from(m.values), ids, genes: m.genes, labels: regionLabels,
               set_names: names, set_genes: names.flatMap((s) => sets[s]),
               set_index: names.flatMap((s, i) => sets[s].map(() => i + 1)) } });
    const minGenes = Math.max(2, +$("min-genes").value || 5);
    state.results = [];
    const t0 = performance.now();
    for (let i = 0; i < names.length; i++) {
      msg.textContent = `Fitting ${i + 1} of ${names.length}: ${pretty(names[i])}`;
      const r = await webR.evalRRaw(`piste_one(${i + 1}L, ${minGenes})`, "number[]");
      const row = { name: names[i], size: sets[names[i]].length, genes: r[1], status: r[0], secs: r[6] };
      if (r[0] === 1) Object.assign(row, { pi: r[2], hi: r[3], lo: r[4], call: r[5] === 1 ? "higher" : r[5] === -1 ? "lower" : "none" });
      else if (r[0] === 0) row.call = "skipped";
      else { row.call = "failed"; row.message = await webR.evalRString(`piste_error(${i + 1}L, ${minGenes})`); }
      state.results.push(row);
      bar.style.width = `${(100 * (i + 1)) / names.length}%`;
      if (i % 3 === 2 || i === names.length - 1) render();
      await new Promise((res) => setTimeout(res, 0));
    }
    const secs = (performance.now() - t0) / 1000;
    msg.textContent = `Done in ${secs < 90 ? secs.toFixed(0) + " s" : (secs / 60).toFixed(1) + " min"} · ${m.spots.length.toLocaleString()} spots, ${m.genes.length.toLocaleString()} genes`;
    state.lastRun = { spots: m.spots.length, region: regionLabels.filter((l) => l === "Region").length };
    render();
    $("step-results").scrollIntoView();
  } catch (e) {
    console.error(e);
    msg.className = "run-msg err"; msg.textContent = e.message;
  }
  running = false; updateRun(); markSteps();
}

// ---------------------------------------------------------------------------
// Results
// ---------------------------------------------------------------------------
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
const pretty = (s) => s.replace(/^HALLMARK_/, "").replace(/_/g, " ").toLowerCase().replace(/^./, (c) => c.toUpperCase());
const fmt = (v) => (v === undefined || v === null || Number.isNaN(v) ? "–" : v >= 0.9995 ? ">0.999" : v.toFixed(3));

function render() {
  const rows = state.results;
  $("results-empty").hidden = rows.length > 0;
  $("results").hidden = rows.length === 0;
  $("btn-download").disabled = rows.length === 0;
  const count = (c) => rows.filter((r) => r.call === c).length;
  $("summary").innerHTML = rows.length ? `
    <div class="stat"><div class="v">${rows.filter((r) => r.status === 1).length}</div><div class="l">gene sets fitted</div></div>
    <div class="stat hi"><div class="v">${count("higher")}</div><div class="l">higher in region</div></div>
    <div class="stat lo"><div class="v">${count("lower")}</div><div class="l">lower in region</div></div>
    <div class="stat"><div class="v">${count("none")}</div><div class="l">no call</div></div>` : "";

  const key = state.sort.key, sgn = state.sort.desc ? -1 : 1;
  const val = (r) => (key === "name" ? pretty(r.name) : r[key] ?? -Infinity);
  const shown = rows
    .filter((r) => state.filter === "all" || r.call === state.filter || (state.filter === "none" && r.status !== 1))
    .sort((a, b) => (val(a) > val(b) ? sgn : val(a) < val(b) ? -sgn : 0));
  $("table").querySelector("tbody").innerHTML = shown.map((r) => {
    const bar = r.status === 1 ? (() => {
      const lo = Math.min(r.pi, 0.5), hi = Math.max(r.pi, 0.5);
      const col = r.call === "higher" ? "#d23b2b" : r.call === "lower" ? "#5368a6" : "#aab5cf";
      return `<div class="pibar"><div class="track"><div class="fill" style="left:${lo * 100}%;width:${(hi - lo) * 100}%;background:${col}"></div><div class="mid"></div></div><span class="val">${r.pi.toFixed(3)}</span></div>`;
    })() : `<span class="hint">${r.status === 0 ? `fewer than ${$("min-genes").value} genes` : esc(r.message || "failed")}</span>`;
    return `<tr>
      <td class="name">${esc(pretty(r.name))}<small>${esc(r.name)}</small></td>
      <td class="r">${r.genes}/${r.size}</td>
      <td>${bar}</td>
      <td class="r">${fmt(r.hi)}</td><td class="r">${fmt(r.lo)}</td>
      <td><span class="badge ${r.status === 1 ? r.call : "skipped"}">${r.status === 1 ? (r.call === "none" ? "none" : r.call) : r.call}</span></td>
    </tr>`;
  }).join("");
  const keys = ["name", "genes", "pi", "hi", "lo", "call"];
  $("table").querySelectorAll("th").forEach((th, j) => {
    th.classList.toggle("sorted", keys[j] === key);
    th.classList.toggle("asc", keys[j] === key && !state.sort.desc);
  });
  markSteps();
}

$("table").querySelectorAll("th").forEach((th, j) => th.addEventListener("click", () => {
  const key = ["name", "genes", "pi", "hi", "lo", "call"][j];
  state.sort = { key, desc: state.sort.key === key ? !state.sort.desc : key !== "name" };
  render();
}));
document.querySelectorAll(".chip").forEach((c) => c.addEventListener("click", () => {
  document.querySelectorAll(".chip").forEach((x) => x.classList.toggle("active", x === c));
  state.filter = c.dataset.filter; render();
}));

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
window.addEventListener("error", (e) => { $("engine").classList.add("error"); $("engine-text").textContent = "Error: " + e.message; });
startR();
