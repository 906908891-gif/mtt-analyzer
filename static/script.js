

const COLORS = ['#0891b2', '#f59e0b', '#10b981', '#8b5cf6', '#ef4444', '#06b6d4', '#ec4899', '#84cc16', '#6366f1', '#14b8a6'];
const PLATE_KEY_SEP = '::';
function makeWellKey(plateId, wid) { return plateId + PLATE_KEY_SEP + wid; }
function parsePlateFromKey(k) { return k.split(PLATE_KEY_SEP)[0]; }
function parseWellFromKey(k) { return k.split(PLATE_KEY_SEP)[1]; }
function wellIdFromRowCol(row, col) { return String.fromCharCode(64 + row) + col; }
const state = {
  plates: [{ id: 'p1001', name: '板 1', format: 96, wellData: {} }],
  activePlateId: 'p1001',
  groups: [],
  activeGroupId: null,
  selectedWells: new Set(),
  controlGroupId: null,
  lastClickedWell: null,
  chartMode: 'abs',
  statsMode: 'aggregate',
};
function showToast(msg, type) {
  type = type || 'success';
  const t = document.getElementById('toast');
  t.textContent = msg;
  t.className = 'toast show ' + type;
  clearTimeout(t._timer);
  t._timer = setTimeout(function() { t.classList.remove('show'); }, 2200);
}
function escapeHtml(s) {
  return String(s).replace(/[&<>"']/g, function(m) { return ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'})[m]; });
}
function truncate(s, n) { return s.length > n ? s.slice(0, n - 1) + '...' : s; }
function genId(prefix) { return prefix + Date.now().toString(36) + Math.floor(Math.random()*1000).toString(36); }
function getActivePlate() { return state.plates.find(function(p) { return p.id === state.activePlateId; }); }
function plateLayout(p) { p = p || getActivePlate(); return p.format === 96 ? { rows:8, cols:12, well:34, max:96 } : { rows:4, cols:6, well:58, max:24 }; }
function findGroupByWell(plateId, wid) { return state.groups.find(function(g) { return g.wells.has(makeWellKey(plateId, wid)); }); }
function mean(arr) { return arr.reduce(function(a,b){return a+b;}, 0) / arr.length; }
function sd(arr) { if (arr.length < 2) return 0; const m = mean(arr); return Math.sqrt(arr.reduce(function(a,b){return a + (b-m)*(b-m);}, 0) / (arr.length - 1)); }
function sem(arr) { return arr.length < 2 ? 0 : sd(arr) / Math.sqrt(arr.length); }
function cv(arr) { return mean(arr) === 0 ? 0 : (sd(arr) / mean(arr)) * 100; }
function logGamma(x) {
  const c = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231757687645645, 0.1208650973866179e-2, -0.5395239384953e-5];
  let y = x, t = x + 5.5;
  t -= (x + 0.5) * Math.log(t);
  let s = 1.000000000190015;
  for (let j = 0; j < 6; j++) s += c[j] / ++y;
  return -t + Math.log(2.5066282746310005 * s / x);
}
function incompleteBeta(a, b, x) {
  // Regularized incomplete beta function I_x(a, b).
  // Numerical Recipes continued-fraction algorithm (Press et al., 1992).
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  // Symmetry: I_x(a, b) = 1 - I_(1-x)(b, a)
  if (x > (a + 1) / (a + b + 2)) {
    return 1 - incompleteBeta(b, a, 1 - x);
  }
  const bt = Math.exp(
    logGamma(a + b) - logGamma(a) - logGamma(b)
    + a * Math.log(x) + b * Math.log(1 - x)
  );
  return bt * betacf(a, b, x) / a;
}

function betacf(a, b, x) {
  // Continued fraction for I_x(a, b) (NR Press et al., §6.4).
  const MAXIT = 200;
  const EPS = 3e-12;
  const FPMIN = 1e-30;
  const qab = a + b;
  const qap = a + 1.0;
  const qam = a - 1.0;
  let c = 1.0;
  let d = 1.0 - qab * x / qap;
  if (Math.abs(d) < FPMIN) d = FPMIN;
  d = 1.0 / d;
  let h = d;
  for (let m = 1; m <= MAXIT; m++) {
    const m2 = 2 * m;
    let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
    d = 1.0 + aa * d;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1.0 + aa / c;
    if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1.0 / d;
    h *= d * c;
    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
    d = 1.0 + aa * d;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1.0 + aa / c;
    if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1.0 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1.0) < EPS) break;
  }
  return h;
}
function welchTTest(a1, a2) {
  const m1 = mean(a1), m2 = mean(a2);
  const v1 = sd(a1) * sd(a1) / a1.length;
  const v2 = sd(a2) * sd(a2) / a2.length;
  if (v1 + v2 === 0) return { t:0, df:1, p:1 };
  const se = Math.sqrt(v1 + v2);
  const t = (m1 - m2) / se;
  const dDen = (v1 * v1) / (a1.length - 1) + (v2 * v2) / (a2.length - 1);
  const df = dDen === 0 ? 1 : (v1 + v2) * (v1 + v2) / dDen;
  const xx = df / (df + t * t);
  return { t:t, df:df, p:incompleteBeta(df/2, 0.5, xx) };
}
function holmBonferroni(pValues) {
  // Holm-Bonferroni step-down correction.
  // Correct direction: running MAX (not min).
  // Reference: Holm, S. (1979).
  const n = pValues.length;
  if (n === 0) return [];
  const indexed = pValues.map(function(p, i) { return { p: p, i: i }; }).sort(function(a, b) { return a.p - b.p; });
  const adjusted = new Array(n);
  let cumMax = 0;
  for (let k = 0; k < n; k++) {
    const adj = Math.min(1, indexed[k].p * (n - k));
    cumMax = Math.max(cumMax, adj);
    adjusted[indexed[k].i] = cumMax;
  }
  return adjusted;
}
function significanceStars(p) {
  if (p == null) return '';
  if (p < 0.001) return '***';
  if (p < 0.01) return '**';
  if (p < 0.05) return '*';
  return 'ns';
}
function formatP(p) { if (p == null) return '-'; if (p < 0.0001) return '< 0.0001'; return p.toFixed(4); }
function parseConc(s) {
  if (!s) return null;
  s = String(s).trim();
  if (!s) return null;
  const sciMatch = s.match(/^([\d.]+(?:e[-+]?\d+)?)$/i);
  if (sciMatch) { const v = parseFloat(s); if (!isNaN(v) && v > 0) return v; }
  const unitMatch = s.match(/^([\d.]+(?:e[-+]?\d+)?)\s*([pnuKM]?M)$/i);
  if (unitMatch) {
    const value = parseFloat(unitMatch[1]);
    const unit = unitMatch[2].toLowerCase();
    let mult = 1;
    if (unit === 'nm') mult = 1e-9;
    else if (unit === 'um') mult = 1e-6;
    else if (unit === 'mm') mult = 1e-3;
    else if (unit === 'pm') mult = 1e-12;
    return value * mult;
  }
  const num = parseFloat(s);
  return (isNaN(num) || num <= 0) ? null : num;
}
function formatConc(M, digits) {
  digits = digits || 2;
  if (M == null || isNaN(M)) return '-';
  if (M === 0) return '0';
  if (M >= 1) return M.toFixed(digits) + ' M';
  if (M >= 1e-3) return (M * 1e3).toFixed(digits) + ' mM';
  if (M >= 1e-6) return (M * 1e6).toFixed(digits) + ' uM';
  if (M >= 1e-9) return (M * 1e9).toFixed(digits) + ' nM';
  if (M >= 1e-12) return (M * 1e12).toFixed(digits) + ' pM';
  return M.toExponential(digits);
}
function solveLinear4x4(A, b) {
  const n = 4;
  const M = [];
  for (let i = 0; i < n; i++) { M.push(A[i].slice()); M[i].push(b[i]); }
  for (let k = 0; k < n; k++) {
    let maxRow = k;
    for (let i = k + 1; i < n; i++) if (Math.abs(M[i][k]) > Math.abs(M[maxRow][k])) maxRow = i;
    if (maxRow !== k) { const tmp = M[k]; M[k] = M[maxRow]; M[maxRow] = tmp; }
    if (Math.abs(M[k][k]) < 1e-18) return null;
    for (let i = k + 1; i < n; i++) {
      const factor = M[i][k] / M[k][k];
      for (let j = k; j <= n; j++) M[i][j] -= factor * M[k][j];
    }
  }
  const x = [0, 0, 0, 0];
  for (let i = n - 1; i >= 0; i--) {
    let sum = M[i][n];
    for (let j = i + 1; j < n; j++) sum -= M[i][j] * x[j];
    x[i] = sum / M[i][i];
  }
  return x;
}
function invert4x4(M) {
  const n = 4;
  const aug = [];
  for (let i = 0; i < n; i++) {
    aug.push(M[i].slice());
    for (let j = 0; j < n; j++) aug[i].push(i === j ? 1 : 0);
  }
  for (let k = 0; k < n; k++) {
    let maxRow = k;
    for (let i = k + 1; i < n; i++) if (Math.abs(aug[i][k]) > Math.abs(aug[maxRow][k])) maxRow = i;
    if (maxRow !== k) { const tmp = aug[k]; aug[k] = aug[maxRow]; aug[maxRow] = tmp; }
    if (Math.abs(aug[k][k]) < 1e-18) return null;
    const piv = aug[k][k];
    for (let j = 0; j < 8; j++) aug[k][j] /= piv;
    for (let i = 0; i < n; i++) {
      if (i === k) continue;
      const factor = aug[i][k];
      for (let j = 0; j < 8; j++) aug[i][j] -= factor * aug[k][j];
    }
  }
  const inv = [];
  for (let i = 0; i < n; i++) inv.push(aug[i].slice(4));
  return inv;
}
function lm4pl(logX, y) {
  const N = logX.length;
  if (N < 4) return null;
  const yMax = Math.max.apply(null, y);
  const yMin = Math.min.apply(null, y);
  const sortedLX = logX.slice().sort(function(a, b) { return a - b; });
  const lxMed = sortedLX[Math.floor(N / 2)];
  let params = [yMax, yMin, lxMed, 1.0];
  let lambda = 0.001;
  function predict(p) {
    const T = p[0], B = p[1], L = p[2], h = p[3];
    return logX.map(function(lx) {
      const denom = 1 + Math.pow(10, (lx - L) * h);
      return B + (T - B) / denom;
    });
  }
  function residuals(p) {
    const yh = predict(p);
    return y.map(function(yi, i) { return yi - yh[i]; });
  }
  function jacobian(p) {
    const eps = 1e-7;
    const r0 = residuals(p);
    const J = [[], [], [], []];
    for (let j = 0; j < 4; j++) {
      const pp = p.slice();
      pp[j] += eps;
      const re = residuals(pp);
      for (let i = 0; i < N; i++) J[j][i] = (re[i] - r0[i]) / eps;
    }
    return J;
  }
  let converged = false;
  for (let iter = 0; iter < 250; iter++) {
    const r = residuals(params);
    const ss0 = r.reduce(function(s, ri) { return s + ri*ri; }, 0);
    if (ss0 < 1e-14) { converged = true; break; }
    const J = jacobian(params);
    const JtJ = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]];
    const Jtr = [0,0,0,0];
    for (let a = 0; a < 4; a++) {
      for (let i = 0; i < N; i++) {
        Jtr[a] += J[a][i] * r[i];
        for (let b = 0; b < 4; b++) JtJ[a][b] += J[a][i] * J[b][i];
      }
    }
    const JtJd = JtJ.map(function(row, a) { return row.map(function(v, b) { return a === b ? v * (1 + lambda) : v; }); });
    const delta = solveLinear4x4(JtJd, Jtr.map(function(v) { return -v; }));
    if (!delta) { lambda *= 5; if (lambda > 1e10) break; continue; }
    const newParams = params.map(function(p, k) { return p + delta[k]; });
    const rNew = residuals(newParams);
    const ssNew = rNew.reduce(function(s, ri) { return s + ri*ri; }, 0);
    if (ssNew < ss0) {
      params = newParams;
      lambda = Math.max(lambda * 0.7, 1e-12);
      if (iter > 3 && Math.max.apply(null, delta.map(Math.abs)) < 1e-10 && (ss0 - ssNew) / ss0 < 1e-10) { converged = true; break; }
    } else {
      lambda = Math.min(lambda * 5, 1e10);
      if (lambda > 1e10) break;
    }
  }
  const yh = predict(params);
  const yMean = y.reduce(function(s, yi) { return s + yi; }, 0) / N;
  const ssRes = y.reduce(function(s, yi, i) { return s + (yi - yh[i]) * (yi - yh[i]); }, 0);
  const ssTot = y.reduce(function(s, yi) { return s + (yi - yMean) * (yi - yMean); }, 0);
  const R2 = ssTot > 0 ? 1 - ssRes / ssTot : 0;
  const dof = Math.max(1, N - 4);
  const s2 = ssRes / dof;
  const Jf = jacobian(params);
  const JtJf = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]];
  for (let a = 0; a < 4; a++) {
    for (let b = 0; b < 4; b++) {
      let sm = 0;
      for (let i = 0; i < N; i++) sm += Jf[a][i] * Jf[b][i];
      JtJf[a][b] = sm;
    }
  }
  const invJtJ = invert4x4(JtJf);
  const se4 = [NaN, NaN, NaN, NaN];
  if (invJtJ) for (let k = 0; k < 4; k++) se4[k] = Math.sqrt(Math.abs(invJtJ[k][k]) * s2);
  const logEC50 = params[2];
  const EC50 = Math.pow(10, logEC50);
  const seLogEC50 = se4[2];
  return {
    Top: params[0], Bottom: params[1], logEC50: logEC50, HillSlope: params[3],
    EC50: EC50, R2: R2, ssRes: ssRes, ssTot: ssTot, nPoints: N, converged: converged,
    seTop: se4[0], seBottom: se4[1], seLogEC50: seLogEC50, seHillSlope: se4[3],
    ci95Low: isFinite(seLogEC50) ? Math.pow(10, logEC50 - 1.96 * seLogEC50) : NaN,
    ci95High: isFinite(seLogEC50) ? Math.pow(10, logEC50 + 1.96 * seLogEC50) : NaN
  };
}
function getGroupValues(g) { const out=[]; g.wells.forEach(function(k){const pid=parsePlateFromKey(k);const wid=parseWellFromKey(k);const pl=state.plates.find(function(p){return p.id===pid;});if(pl){const v=pl.wellData[wid];if(v!==undefined&&v!==''&&v!==null)out.push(v);}});return out; }
function getGroupValuesForPlate(g, plateId) { const out=[]; g.wells.forEach(function(k){if(parsePlateFromKey(k)!==plateId)return;const wid=parseWellFromKey(k);const pl=state.plates.find(function(p){return p.id===plateId;});if(pl){const v=pl.wellData[wid];if(v!==undefined&&v!==''&&v!==null)out.push(v);}});return out; }
function renderPlatesBar() {
  const bar = document.getElementById('plates-bar');
  let html = state.plates.map(function(p) {
    const isActive = p.id === state.activePlateId;
    const dataCount = Object.keys(p.wellData).length;
    const max = plateLayout(p).max;
    return '<div class="plate-tab ' + (isActive ? 'active' : '') + '" data-plate="' + p.id + '">' +
      '<span class="plate-tab-name">' + escapeHtml(p.name) + '</span>' +
      '<span class="plate-tab-status">' + dataCount + '/' + (max === 96 ? '96' : '24') + '</span>' +
      '</div>';
  }).join('');
  html += '<div class="plate-add-btn" id="add-plate-btn">+ 新建板</div>';
  html += '<div class="plate-actions">' +
    '<button class="icon-btn secondary" id="rename-plate-btn">重命名</button>' +
    '<button class="icon-btn danger" id="delete-plate-btn">删除板</button>' +
    '</div>';
  bar.innerHTML = html;
  document.querySelectorAll('.plate-tab').forEach(function(el) { el.addEventListener('click', function() { if (el.dataset.plate !== state.activePlateId) switchToPlate(el.dataset.plate); }); });
  document.getElementById('add-plate-btn').addEventListener('click', addNewPlate);
  document.getElementById('rename-plate-btn').addEventListener('click', renameActivePlate);
  document.getElementById('delete-plate-btn').addEventListener('click', deleteActivePlate);
  const ap = getActivePlate();
  document.getElementById('plate-format-select').value = ap.format;
  document.getElementById('plate-name-input').value = ap.name;
  document.getElementById('currentPlateLabel').textContent = ap.name + ' (' + (ap.format === 96 ? '96' : '24') + ' 孔)';
  document.getElementById('plateCountBadge').textContent = state.plates.length + ' 板';
  document.getElementById('groupCountBadge').textContent = state.groups.length + ' 组';
}
function switchToPlate(pid) { state.activePlateId = pid; state.selectedWells.clear(); state.lastClickedWell = null; renderPlatesBar(); renderPlate(); renderGroups(); renderGroupSelect(); renderResults(); }
function addNewPlate() { const newFmt = getActivePlate().format; const id = genId('p'); const p = { id:id, name:'板 ' + (state.plates.length + 1), format:newFmt, wellData:{} }; state.plates.push(p); switchToPlate(id); showToast('已新建 ' + p.name); }
function renameActivePlate() { const ap = getActivePlate(); const newName = prompt('给当前板起个名字：', ap.name); if (newName == null || newName.trim() === '') return; ap.name = newName.trim().slice(0, 32); renderPlatesBar(); renderResults(); updateSessionInfo(); showToast('板已重命名', 'info'); }
function deleteActivePlate() {
  if (state.plates.length === 1) { showToast('至少保留一块板', 'error'); return; }
  const ap = getActivePlate();
  if (!confirm('确定要删除「' + ap.name + '」吗？\n所有该板的数据将丢失；已分配该板孔的分组会自动移除这些孔。')) return;
  const pidToDelete = ap.id;
  state.plates = state.plates.filter(function(p) { return p.id !== pidToDelete; });
  state.groups.forEach(function(g) { const toDel = []; g.wells.forEach(function(k) { if (parsePlateFromKey(k) === pidToDelete) toDel.push(k); }); toDel.forEach(function(k) { g.wells.delete(k); }); });
  if (state.controlGroupId) { const cg = state.groups.find(function(g) { return g.id === state.controlGroupId; }); if (!cg || cg.wells.size === 0) state.controlGroupId = null; }
  state.activePlateId = state.plates[0].id;
  state.selectedWells.clear();
  renderPlatesBar(); renderPlate(); renderGroups(); renderGroupSelect(); renderResults();
  updateSessionInfo();
  showToast('已删除板', 'info');
}
function renderPlate() {
  const ap = getActivePlate();
  const layout = plateLayout(ap);
  const rows = layout.rows, cols = layout.cols, wellSize = layout.well;
  const container = document.getElementById('plate-container');
  const template = 'repeat(' + cols + ', ' + wellSize + 'px)';
  let html = '<div class="plate-container"><div style="text-align:center;">';
  html += '<div style="display:grid; grid-template-columns: 22px ' + template + '; gap:4px; margin-bottom:4px;">';
  html += '<div></div>';
  for (let c = 1; c <= cols; c++) html += '<div class="plate-col-header">' + c + '</div>';
  html += '</div>';
  for (let r = 1; r <= rows; r++) {
    html += '<div style="display:grid; grid-template-columns: 22px ' + template + '; gap:4px; margin-bottom:4px;">';
    html += '<div class="plate-row-header">' + String.fromCharCode(64 + r) + '</div>';
    for (let c = 1; c <= cols; c++) {
      const wid = wellIdFromRowCol(r, c);
      const key = makeWellKey(ap.id, wid);
      const data = ap.wellData[wid];
      const isSelected = state.selectedWells.has(key);
      const group = findGroupByWell(ap.id, wid);
      const hasData = data !== undefined && data !== null && data !== '';
      let cls = 'well';
      if (isSelected) cls += ' selected';
      if (hasData) cls += ' has-data';
      let style = 'width:' + wellSize + 'px;height:' + wellSize + 'px;';
      if (group) { cls += ' in-group'; style += ' background:' + group.color + '; border-color:' + group.color + ';'; }
      const title = wid + (hasData ? ': ' + Number(data).toFixed(4) : '') + (group ? ' [' + group.name + ']' : '');
      html += '<div class="' + cls + '" style="' + style + '" data-well="' + wid + '" title="' + escapeHtml(title) + '">';
      if (hasData) html += '<span style="font-size:' + (wellSize < 40 ? '0.48rem' : '0.68rem') + ';line-height:1;">' + Number(data).toFixed(3) + '</span>';
      html += '</div>';
    }
    html += '</div>';
  }
  html += '</div></div>';
  container.innerHTML = html;
  document.querySelectorAll('.well').forEach(function(el) { el.addEventListener('click', function(e) { handleWellClick(e); }); });
}
function handleWellClick(e) {
  const wid = e.currentTarget.dataset.well;
  const pid = state.activePlateId;
  const key = makeWellKey(pid, wid);
  if (e.shiftKey && state.lastClickedWell && parsePlateFromKey(state.lastClickedWell) === pid) {
    const last = parseWellFromKey(state.lastClickedWell);
    const r1 = last.charCodeAt(0) - 64, c1 = parseInt(last.slice(1), 10);
    const r2 = wid.charCodeAt(0) - 64, c2 = parseInt(wid.slice(1), 10);
    const minR = Math.min(r1, r2), maxR = Math.max(r1, r2);
    const minC = Math.min(c1, c2), maxC = Math.max(c1, c2);
    state.selectedWells.clear();
    for (let r = minR; r <= maxR; r++) for (let c = minC; c <= maxC; c++) state.selectedWells.add(makeWellKey(pid, wellIdFromRowCol(r, c)));
  } else if (e.ctrlKey || e.metaKey) {
    if (state.selectedWells.has(key)) state.selectedWells.delete(key); else state.selectedWells.add(key);
    state.lastClickedWell = key;
  } else {
    state.selectedWells.clear();
    state.selectedWells.add(key);
    state.lastClickedWell = key;
  }
  renderPlate();
  renderGroups();
}
function moveGroup(id, direction) {
  const idx = state.groups.findIndex(function(g) { return g.id === id; });
  if (idx < 0) return;
  const newIdx = idx + direction;
  if (newIdx < 0 || newIdx >= state.groups.length) return;
  const moved = state.groups.splice(idx, 1)[0];
  state.groups.splice(newIdx, 0, moved);
  renderGroups(); renderGroupSelect(); renderResults(); renderPlatesBar();
}
function attachDragHandlers() {
  let draggedId = null;
  document.querySelectorAll('.group-item').forEach(function(el) {
    el.addEventListener('dragstart', function(e) {
      if (e.target.tagName === 'INPUT' || e.target.tagName === 'BUTTON') { e.preventDefault(); return; }
      draggedId = el.dataset.id;
      el.classList.add('dragging');
      try { e.dataTransfer.effectAllowed = 'move'; e.dataTransfer.setData('text/plain', draggedId); } catch (ex) {}
    });
    el.addEventListener('dragend', function() {
      el.classList.remove('dragging');
      document.querySelectorAll('.group-item').forEach(function(i) { i.classList.remove('drop-above', 'drop-below'); });
      draggedId = null;
    });
    el.addEventListener('dragover', function(e) {
      e.preventDefault();
      e.dataTransfer.dropEffect = 'move';
      if (el.dataset.id === draggedId) return;
      const rect = el.getBoundingClientRect();
      const midY = rect.top + rect.height / 2;
      if (e.clientY < midY) { el.classList.add('drop-above'); el.classList.remove('drop-below'); }
      else { el.classList.add('drop-below'); el.classList.remove('drop-above'); }
    });
    el.addEventListener('dragleave', function() { el.classList.remove('drop-above', 'drop-below'); });
    el.addEventListener('drop', function(e) {
      e.preventDefault();
      el.classList.remove('drop-above', 'drop-below');
      const sourceId = e.dataTransfer.getData('text/plain') || draggedId;
      const targetId = el.dataset.id;
      if (!sourceId || sourceId === targetId) return;
      const fromIdx = state.groups.findIndex(function(g) { return g.id === sourceId; });
      if (fromIdx < 0) return;
      const rect = el.getBoundingClientRect();
      const midY = rect.top + rect.height / 2;
      const moved = state.groups.splice(fromIdx, 1)[0];
      const newToIdx = state.groups.findIndex(function(g) { return g.id === targetId; });
      const insertIdx = (e.clientY < midY) ? newToIdx : newToIdx + 1;
      state.groups.splice(insertIdx, 0, moved);
      renderGroups(); renderGroupSelect(); renderResults(); renderPlatesBar();
    });
  });
}
function renderGroups() {
  const list = document.getElementById('groups-list');
  if (state.groups.length === 0) {
    list.innerHTML = '<div class="empty"><div class="empty-icon">📋</div>还没有分组<br><span style="font-size:0.75rem;">点击 "+ 新建分组" 创建</span></div>';
    return;
  }
  const ap = getActivePlate();
  list.innerHTML = state.groups.map(function(g, idx) {
    const isActive = state.activeGroupId === g.id;
    const isControl = state.controlGroupId === g.id;
    let totalWells = g.wells.size, wellsWithData = 0, inActivePlate = 0;
    g.wells.forEach(function(k) {
      const pid = parsePlateFromKey(k), wid = parseWellFromKey(k);
      const plate = state.plates.find(function(p) { return p.id === pid; });
      if (plate && plate.wellData[wid] !== undefined) wellsWithData++;
      if (pid === ap.id) inActivePlate++;
    });
    const countBadge = totalWells + (totalWells !== wellsWithData ? ' (' + wellsWithData + '✓)' : '');
    const ctrlTitle = isControl ? '当前对照（点击取消）' : '设为对照';
    const isFirst = idx === 0, isLast = idx === state.groups.length - 1;
    return '<div class="group-item ' + (isActive ? 'active' : '') + '" draggable="true" data-id="' + g.id + '" title="拖动或点击 ↑↓ 重排 · 当前板内 ' + inActivePlate + ' 个孔">' +
      '<div class="group-color-pick" style="background:' + g.color + '"></div>' +
      '<input type="text" class="group-name-input" value="' + escapeHtml(g.name) + '" data-id="' + g.id + '">' +
      '<span class="group-count" title="' + wellsWithData + ' 个有效数据">' + countBadge + '</span>' +
      '<span class="group-control-icon ' + (isControl ? 'is-control' : '') + '" data-id="' + g.id + '" title="' + ctrlTitle + '">●</span>' +
      '<div class="group-order-btns">' +
      '<button class="group-order-btn up" data-id="' + g.id + '" title="上移"' + (isFirst ? ' disabled' : '') + '>▲</button>' +
      '<button class="group-order-btn down" data-id="' + g.id + '" title="下移"' + (isLast ? ' disabled' : '') + '>▼</button>' +
      '</div>' +
      '<button class="group-delete" data-id="' + g.id + '" title="删除分组">✕</button>' +
      '</div>';}).join('');
  attachGroupHandlers();
  attachDragHandlers();
  document.getElementById('assign-btn').textContent = state.selectedWells.size > 0 ? '分配 ' + state.selectedWells.size + ' 个孔' : '分配到选中孔';
  document.getElementById('unassign-btn').textContent = state.selectedWells.size > 0 ? '移除 ' + state.selectedWells.size + ' 个孔' : '取消已分配';
}
function attachGroupHandlers() {
  document.querySelectorAll('.group-item').forEach(function(el) {
    el.addEventListener('click', function(e) {
      if (e.target.closest('input, button, .group-control-icon, .group-color-pick')) return;
      state.activeGroupId = el.dataset.id;
      renderGroups();
    });
  });
  document.querySelectorAll('.group-name-input').forEach(function(el) {
    el.addEventListener('input', function() {
      const g = state.groups.find(function(g2) { return g2.id === el.dataset.id; });
      if (g) g.name = el.value;
      renderResults(); renderGroupSelect(); renderPlatesBar(); autoSaveDebounced();
    });
    el.addEventListener('click', function(e) { e.stopPropagation(); });
  });
  document.querySelectorAll('.group-delete').forEach(function(el) {
    el.addEventListener('click', function(e) {
      e.stopPropagation();
      const id = el.dataset.id;
      state.groups = state.groups.filter(function(g) { return g.id !== id; });
      if (state.activeGroupId === id) state.activeGroupId = null;
      if (state.controlGroupId === id) state.controlGroupId = null;
      renderPlate(); renderGroups(); renderResults(); renderGroupSelect(); renderPlatesBar(); autoSaveDebounced();
    });
  });
  document.querySelectorAll('.group-control-icon').forEach(function(el) {
    el.addEventListener('click', function(e) {
      e.stopPropagation();
      const id = el.dataset.id;
      const g = state.groups.find(function(g2) { return g2.id === id; });
      if (!g || g.wells.size === 0) { showToast('该分组无数据', 'error'); return; }
      const values = getGroupValues(g);
      if (values.length < 2) { showToast('对照至少需要 2 个有效数据点', 'error'); return; }
      state.controlGroupId = state.controlGroupId === id ? null : id;
      renderGroups(); renderGroupSelect(); renderResults(); autoSaveDebounced();
    });
  });
  document.querySelectorAll('.group-color-pick').forEach(function(el) {
    el.addEventListener('click', function(e) {
      e.stopPropagation();
      const rect = el.getBoundingClientRect();
      const id = el.parentElement.dataset.id;
      const popup = document.createElement('div');
      popup.className = 'color-picker-popup';
      popup.style.top = (rect.bottom + 4) + 'px';
      popup.style.left = rect.left + 'px';
      COLORS.forEach(function(color) {
        const d = document.createElement('div');
        d.style.background = color;
        d.addEventListener('click', function(e2) {
          e2.stopPropagation();
          const g = state.groups.find(function(g2) { return g2.id === id; });
          if (g) { g.color = color; renderPlate(); renderGroups(); renderResults(); autoSaveDebounced(); }
          popup.remove();
        });
        popup.appendChild(d);
      });
      document.body.appendChild(popup);
      const closeFn = function(ev) { if (popup.contains && !popup.contains(ev.target)) { popup.remove(); document.removeEventListener('click', closeFn); } };
      setTimeout(function() { document.addEventListener('click', closeFn); }, 0);
    });
  });
  document.querySelectorAll('.group-order-btn.up').forEach(function(el) {
    el.addEventListener('click', function(e) { e.stopPropagation(); moveGroup(el.dataset.id, -1); });
  });
  document.querySelectorAll('.group-order-btn.down').forEach(function(el) {
    el.addEventListener('click', function(e) { e.stopPropagation(); moveGroup(el.dataset.id, 1); });
  });
}
function renderGroupSelect() {
  const sel = document.getElementById('control-select');
  let opts = '<option value="">- 自动（无对照） -</option>';
  state.groups.filter(function(g) { return g.wells.size > 0; }).forEach(function(g) {
    opts += '<option value="' + g.id + '"' + (state.controlGroupId === g.id ? ' selected' : '') + '>' + escapeHtml(g.name) + ' (' + g.wells.size + ')</option>';
  });
  sel.innerHTML = opts;
}
function computeStats() {
  const controlGroup = state.groups.find(function(g) { return g.id === state.controlGroupId; });
  if (!controlGroup) return { groups: [], hasControl: false };
  const filterFn = state.statsMode === 'currentPlate'
    ? function(g) { return getGroupValuesForPlate(g, state.activePlateId); }
    : function(g) { return getGroupValues(g); };
  const controlValues = filterFn(controlGroup);
  if (controlValues.length < 2) return { groups: [], hasControl: false, reason: '对照分组有效数据 < 2（' + (state.statsMode === 'currentPlate' ? '当前板内' : '全部板') + '）' };
  const allStats = [];
  state.groups.forEach(function(g) {
    const values = filterFn(g);
    if (values.length === 0) return;
    const platesSet = new Set();
    g.wells.forEach(function(k) { platesSet.add(parsePlateFromKey(k)); });
    const platesUsed = [];
    platesSet.forEach(function(pid) { const pl = state.plates.find(function(p) { return p.id === pid; }); platesUsed.push(pl ? pl.name : pid); });
    allStats.push({ group: g, n: values.length, mean: mean(values), sd: sd(values), sem: sem(values), cv: cv(values), values: values, platesUsed: platesUsed });
  });
  if (allStats.length === 0) return { groups: [], hasControl: false };
  const controlStats = allStats.find(function(s) { return s.group.id === state.controlGroupId; });
  const compStats = allStats.filter(function(s) { return s.group.id !== state.controlGroupId; });
  const pRaw = compStats.map(function(s) { return welchTTest(s.values, controlStats.values).p; });
  const pAdj = holmBonferroni(pRaw);
  compStats.forEach(function(s, i) {
    s.pValue = pRaw[i]; s.pAdjusted = pAdj[i];
    s.significance = significanceStars(pAdj[i]);
    s.viability = (s.mean / controlStats.mean) * 100;
  });
  controlStats.viability = 100; controlStats.pValue = null; controlStats.pAdjusted = null; controlStats.significance = '';
  return { groups: allStats, control: controlStats, hasControl: true };
}
function renderResults() {
  renderDoseResponse();
  autoSaveDebounced();
  const stats = computeStats();
  const chartContainer = document.getElementById('chart-container');
  const tableContainer = document.getElementById('results-table-container');
  if (state.groups.length === 0) {
    chartContainer.innerHTML = '<div class="empty"><div class="empty-icon">📊</div>创建分组后此处显示图表</div>';
    tableContainer.innerHTML = '';
    return;
  }
  if (!stats.hasControl) {
    chartContainer.innerHTML = '<div class="empty"><div class="empty-icon">📊</div>选择对照分组后显示图表</div>';
    tableContainer.innerHTML = '<div class="empty"><div class="empty-icon">📋</div>' + (stats.reason || '设置对照后显示统计表') + '</div>';
    return;
  }
  if (state.chartMode === 'heatmap') { chartContainer.innerHTML = renderHeatmap(stats); } else { chartContainer.innerHTML = renderBarChart(stats); }
  tableContainer.innerHTML = renderResultsTable(stats);
}
function renderBarChart(stats) {
  const groups = stats.groups;
  const n = groups.length;
  const width = Math.max(420, n * 78 + 80);
  const height = 270;
  const margin = { top: 36, right: 18, bottom: 56, left: 52 };
  const innerWidth = width - margin.left - margin.right;
  const innerHeight = height - margin.top - margin.bottom;
  const isViab = state.chartMode === 'viability';
  const dataArr = groups.map(function(g) {
    return { value: isViab ? g.viability : g.mean, err: isViab ? ((g.sd / (stats.control ? stats.control.mean : g.mean)) * 100) : g.sd };
  });
  const maxVal = Math.max.apply(null, dataArr.map(function(d) { return d.value + d.err; })) * 1.18;
  const slotW = innerWidth / n;
  const barW = Math.min(slotW * 0.62, 64);
  const barGap = (slotW - barW) / 2;
  let bars = '', labels = '', stars = '';
  groups.forEach(function(s, i) {
    const dd = dataArr[i];
    const x = margin.left + i * slotW + barGap;
    const yBase = margin.top + innerHeight;
    const h = innerHeight * (dd.value / maxVal);
    const yErr = innerHeight * (dd.err / maxVal);
    const yTop = yBase - h;
    const cx = x + barW / 2;
    bars += '<rect class="bar" x="' + x + '" y="' + yTop + '" width="' + barW + '" height="' + h + '" fill="' + s.group.color + '" rx="3" />';
    bars += '<line class="err-bar" x1="' + cx + '" y1="' + (yTop - yErr) + '" x2="' + cx + '" y2="' + (yTop + yErr) + '" />';
    bars += '<line class="err-bar" x1="' + (cx - barW * 0.3) + '" y1="' + (yTop - yErr) + '" x2="' + (cx + barW * 0.3) + '" y2="' + (yTop - yErr) + '" />';
    bars += '<line class="err-bar" x1="' + (cx - barW * 0.3) + '" y1="' + (yTop + yErr) + '" x2="' + (cx + barW * 0.3) + '" y2="' + (yTop + yErr) + '" />';
    labels += '<text class="bar-label" x="' + cx + '" y="' + (yTop - yErr - 6) + '" text-anchor="middle">' + dd.value.toFixed(isViab ? 1 : 3) + (isViab ? '%' : '') + '</text>';
    const nm = s.group.id === state.controlGroupId ? s.group.name + ' *' : s.group.name;
    labels += '<text class="axis-label" x="' + cx + '" y="' + (height - 36) + '" text-anchor="middle" transform="rotate(-18 ' + cx + ' ' + (height - 36) + ')">' + escapeHtml(truncate(nm, 14)) + '</text>';
    if (s.significance && s.group.id !== state.controlGroupId) {
      stars += '<text class="star-label" x="' + cx + '" y="' + (yTop - yErr - 22) + '" text-anchor="middle">' + s.significance + '</text>';
    }
  });
  const yTicks = 5;
  let yAxis = '';
  for (let i = 0; i <= yTicks; i++) {
    const val = maxVal * i / yTicks;
    const y = margin.top + innerHeight - innerHeight * (i / yTicks);
    yAxis += '<line class="' + (i === yTicks ? 'axis-line' : 'grid-line') + '" x1="' + margin.left + '" y1="' + y + '" x2="' + (margin.left + innerWidth) + '" y2="' + y + '" />';
    yAxis += '<text class="axis-label" x="' + (margin.left - 6) + '" y="' + (y + 3) + '" text-anchor="end">' + val.toFixed(isViab ? 0 : 2) + '</text>';
  }
  const yLabel = isViab ? 'Cell Viability (%)' : 'Absorbance';
  const modeNote = state.statsMode === 'currentPlate' ? ' [当前板]' : ' [全部板]';
  return '<svg class="bar-chart" viewBox="0 0 ' + width + ' ' + height + '" preserveAspectRatio="xMidYMid meet">' +
    '<text class="title-label" x="' + margin.left + '" y="14">' + yLabel + (isViab ? '' : ' (mean +/- SD)') + modeNote + '</text>' +
    '<text class="axis-label" transform="rotate(-90, 14, ' + (margin.top + innerHeight / 2) + ')" x="14" y="' + (margin.top + innerHeight / 2) + '" text-anchor="middle">' + yLabel + '</text>' +
    yAxis + bars + labels + stars + '</svg>';
}
function renderResultsTable(stats) {
  const modeNote = state.statsMode === 'currentPlate' ? '当前板数据' : '聚合全部板';
  let rows = '';
  stats.groups.forEach(function(s) {
    const isCtrl = s.group.id === state.controlGroupId;
    const pStr = isCtrl ? '-' : (formatP(s.pAdjusted) + ' ' + s.significance);
    const pTitle = isCtrl ? '' : ("Welch's t-test raw p=" + s.pValue.toExponential(2) + ', Holm-adjusted');
    const plates = (s.platesUsed || []).map(function(p) { return '<span class="plate-tag">' + escapeHtml(p) + '</span>'; }).join('');
    rows += '<tr>' +
      '<td><span class="legend-dot" style="background:' + s.group.color + '"></span>' + escapeHtml(s.group.name) + (isCtrl ? ' <span style="color:var(--success);font-size:0.68rem;">*对照</span>' : '') + '</td>' +
      '<td>' + s.n + '</td>' +
      '<td class="plates-cell">' + plates + '</td>' +
      '<td>' + s.mean.toFixed(4) + '</td>' +
      '<td>' + s.sd.toFixed(4) + '</td>' +
      '<td>' + s.sem.toFixed(4) + '</td>' +
      '<td>' + s.cv.toFixed(1) + '%</td>' +
      '<td>' + s.viability.toFixed(1) + '%</td>' +
      '<td class="p-value" title="' + pTitle + '">' + pStr + '</td>' +
      '</tr>';
  });
  return '<table class="results"><thead><tr>' +
    '<th>分组</th><th>n</th><th>板</th><th>Mean</th><th>SD</th><th>SEM</th><th>CV</th><th>Viability</th><th>p (adj.)</th>' +
    '</tr></thead><tbody>' + rows + '</tbody></table>' +
    '<p class="hint" style="margin-top:0.4rem;">统计样本：<strong>' + modeNote + '</strong> · "板" 列显示来源板名 · p 值由 Welch t test 计算，Holm-Bonferroni 校正。<span style="color:var(--danger);font-weight:600;">*</span> p&lt;0.05 <span style="color:var(--danger);font-weight:600;">**</span> p&lt;0.01 <span style="color:var(--danger);font-weight:600;">***</span> p&lt;0.001 ns = 不显著。</p>';
}
function fitDoseResponse() {
  const inputs = document.querySelectorAll('.dr-conc-input');
  const points = [];
  const controlGroup = state.groups.find(function(g) { return g.id === state.controlGroupId; });
  let ctrlMean = null;
  if (controlGroup) {
    const cv = state.statsMode === 'currentPlate' ? getGroupValuesForPlate(controlGroup, state.activePlateId) : getGroupValues(controlGroup);
    if (cv.length >= 1) ctrlMean = mean(cv);
  }
  inputs.forEach(function(inp) {
    const gid = inp.dataset.id;
    const g = state.groups.find(function(g2) { return g2.id === gid; });
    if (!g) return;
    const conc = parseConc(inp.value);
    if (conc == null || conc <= 0) return;
    const absVals = state.statsMode === 'currentPlate' ? getGroupValuesForPlate(g, state.activePlateId) : getGroupValues(g);
    if (absVals.length === 0) return;
    const vVals = ctrlMean ? absVals.map(function(a) { return (a / ctrlMean) * 100; }) : absVals.slice();
    points.push({
      group: g, name: g.name, color: g.color,
      conc: conc, logConc: Math.log10(conc), concRaw: inp.value,
      n: absVals.length,
      meanAbs: mean(absVals), sdAbs: sd(absVals),
      viability: ctrlMean ? mean(vVals) : NaN,
      sdViab: ctrlMean ? sd(vVals) : NaN,
      semViab: ctrlMean ? sem(vVals) : NaN
    });
  });
  if (points.length === 0) return { points: points, fit: null, error: '没有分组同时有浓度和数据', ctrlMean: ctrlMean };
  if (ctrlMean == null) {
    const maxMean = Math.max.apply(null, points.map(function(p) { return p.meanAbs; }));
    ctrlMean = maxMean;
    points.forEach(function(p) {
      p.viability = (p.meanAbs / ctrlMean) * 100;
      p.sdViab = (p.sdAbs / ctrlMean) * 100;
      p.semViab = p.sdViab / Math.sqrt(p.n);
    });
  }
  points.sort(function(a, b) { return a.conc - b.conc; });
  if (points.length < 3) return { points: points, fit: null, error: '至少需要 3 个有浓度数据的分组进行拟合（当前 ' + points.length + ' 个）', ctrlMean: ctrlMean };
  const logX = points.map(function(p) { return p.logConc; });
  const y = points.map(function(p) { return p.viability; });
  const fit = lm4pl(logX, y);
  if (!fit) return { points: points, fit: null, error: '拟合失败（数据可能异常）', ctrlMean: ctrlMean };
  return { points: points, fit: fit, ctrlMean: ctrlMean };
}
function renderDoseResponse() {
  const container = document.getElementById('dose-response-container');
  if (!container) return;
  if (state.groups.length === 0) {
    container.innerHTML = '<div class="empty"><div class="empty-icon">📈</div>创建分组后可拟合 4PL 曲线</div>';
    return;
  }
  let gridHtml = '<div class="dr-grid">';
  gridHtml += '<div class="dr-grid-row dr-header"><div>分组</div><div>n</div><div>Mean Abs</div><div>浓度 (M)</div></div>';
  state.groups.forEach(function(g) {
    let n = 0, m = NaN;
    if (g.wells.size > 0) {
      const vals = state.statsMode === 'currentPlate' ? getGroupValuesForPlate(g, state.activePlateId) : getGroupValues(g);
      n = vals.length;
      if (n > 0) m = mean(vals);
    }
    gridHtml += '<div class="dr-grid-row">' +
      '<div><span class="legend-dot" style="background:' + g.color + '"></span>' + escapeHtml(g.name) + '</div>' +
      '<div>' + n + '</div>' +
      '<div style="font-family:ui-monospace,monospace;">' + (n > 0 ? m.toFixed(4) : '-') + '</div>' +
      '<div><input type="text" class="dr-conc-input" data-id="' + g.id + '" placeholder="1uM, 5nM, 1e-6..." /></div>' +
      '</div>';
  });
  gridHtml += '</div>';
  gridHtml += '<div class="actions" style="margin:0.5rem 0;">' +
    '<button id="fit-4pl-btn">📈 拟合 4PL 曲线</button>' +
    '<button id="extract-conc-btn" class="secondary tiny">从分组名自动提取浓度</button>' +
    '<button id="clear-conc-btn" class="secondary tiny">清空浓度</button>' +
    '</div>';
  gridHtml += '<div id="dr-result-container"></div>';
  container.innerHTML = gridHtml;
  let autoExtracted = 0;
  state.groups.forEach(function(g) {
    const inp = container.querySelector('.dr-conc-input[data-id="' + g.id + '"]');
    if (!inp) return;
    if (g._concRaw != null) { inp.value = g._concRaw; return; }
    const m = g.name.match(/([\d.]+(?:e[-+]?\d+)?)\s*(pM|nM|uM|um|μM|mM|M)\b/i);
    if (m) {
      const txt = m[1] + (m[2].toLowerCase() === 'um' ? 'uM' : m[2]);
      const parsed = parseConc(txt);
      if (parsed != null && parsed > 0) { inp.value = txt; g._concRaw = txt; autoExtracted++; }
    }
  });
  if (autoExtracted > 0 && !window._drAutoToastShown) {
    window._drAutoToastShown = true;
    setTimeout(function() { showToast('已从分组名自动提取 ' + autoExtracted + ' 个浓度（可直接编辑）', 'info'); }, 300);
  }
  document.getElementById('fit-4pl-btn').addEventListener('click', runFitAndRender);
  document.getElementById('extract-conc-btn').addEventListener('click', function() {
    let extracted = 0;
    state.groups.forEach(function(g) {
      const inp = container.querySelector('.dr-conc-input[data-id="' + g.id + '"]');
      if (!inp || inp.value) return;
      const m = g.name.match(/([\d.]+(?:e[-+]?\d+)?)\s*(pM|nM|uM|um|μM|mM|M)\b/i);
      if (m) { inp.value = m[1] + (m[2].toLowerCase() === 'um' ? 'uM' : m[2]); extracted++; }
    });
    showToast(extracted > 0 ? ('已从分组名提取 ' + extracted + ' 个浓度') : '未在分组名中找到浓度（请手动输入）', extracted > 0 ? 'success' : 'info');
  });
  document.getElementById('clear-conc-btn').addEventListener('click', function() {
    state.groups.forEach(function(g) { g._concRaw = null; });
    document.querySelectorAll('.dr-conc-input').forEach(function(i) { i.value = ''; });
    document.getElementById('dr-result-container').innerHTML = '';
  });
  const lastResult = window._lastDrResult;
  if (lastResult && lastResult.points && lastResult.points.length >= 3 && lastResult.fit) {
    document.getElementById('dr-result-container').innerHTML = renderDRChart(lastResult) + renderDRParams(lastResult) + renderDRWarnings(lastResult);
  }
}
function runFitAndRender() {
  state.groups.forEach(function(g) {
    const inp = document.querySelector('.dr-conc-input[data-id="' + g.id + '"]');
    if (inp) g._concRaw = inp.value;
  });
  autoSaveDebounced();
  const result = fitDoseResponse();
  window._lastDrResult = result;
  const rc = document.getElementById('dr-result-container');
  if (!rc) return;
  if (result.error) { rc.innerHTML = '<div class="dr-warning">⚠ ' + result.error + '</div>'; return; }
  rc.innerHTML = renderDRChart(result) + renderDRParams(result) + renderDRWarnings(result);
  showToast('拟合完成 · IC50 = ' + formatConc(result.fit.EC50) + ' · R² = ' + result.fit.R2.toFixed(4), 'success');
}
function renderDRWarnings(result) {
  const f = result.fit;
  let html = '';
  if (f.nPoints < 5) html += '<div class="dr-warning">⚠ 只有 ' + f.nPoints + ' 个数据点 · 建议至少 5-8 个浓度以稳定 IC50 估计</div>';
  if (f.HillSlope < 0.5 || f.HillSlope > 4) html += '<div class="dr-warning">⚠ HillSlope = ' + f.HillSlope.toFixed(3) + (f.HillSlope < 0.5 ? ' 偏小' : ' 偏大') + ' · 典型值 ~1.0</div>';
  if (f.Bottom < -10) html += '<div class="dr-warning">⚠ Bottom = ' + f.Bottom.toFixed(1) + '% · 低于 0% 可能表示基线扣除不完全</div>';
  if (f.Top > 120 || f.Top < 80) html += '<div class="dr-warning">⚠ Top = ' + f.Top.toFixed(1) + '% · 偏离 100% 较多</div>';
  const ciFold = f.ci95High / f.ci95Low;
  if (isFinite(ciFold) && ciFold > 10) html += '<div class="dr-warning">⚠ IC50 95% CI 很宽（' + ciFold.toFixed(1) + ' 倍）· 数据需要更多浓度点</div>';
  if (!f.converged) html += '<div class="dr-warning">⚠ 优化未完全收敛，IC50 估计可能不准确</div>';
  if (f.R2 >= 0.95) html = '<div class="dr-success">✓ 拟合优良：R² = ' + f.R2.toFixed(4) + '</div>' + html;
  return html;
}
function renderDRChart(result) {
  const points = result.points;
  const fit = result.fit;
  const nPts = points.length;
  const minLX = points[0].logConc;
  const maxLX = points[nPts - 1].logConc;
  const range = maxLX - minLX;
  const padding = Math.max(0.6, range * 0.2);
  const lxMin = Math.floor(minLX - padding);
  const lxMax = Math.ceil(maxLX + padding);
  const yDataMax = Math.max.apply(null, points.map(function(p) { return p.viability + (isFinite(p.semViab) ? p.semViab : 0); }));
  const yTop = fit ? Math.max(yDataMax, fit.Top) * 1.12 : yDataMax * 1.15;
  const yBot = fit ? Math.min(0, fit.Bottom - 5) : 0;
  const yRange = yTop - yBot;
  const width = 720, height = 330;
  const margin = { top: 30, right: 30, bottom: 56, left: 60 };
  const innerW = width - margin.left - margin.right;
  const innerH = height - margin.top - margin.bottom;
  function xScale(lx) { return margin.left + (lx - lxMin) / (lxMax - lxMin) * innerW; }
  function yScale(y) { return margin.top + innerH - (y - yBot) / yRange * innerH; }
  let svg = '<svg class="dr-chart" viewBox="0 0 ' + width + ' ' + height + '" preserveAspectRatio="xMidYMid meet">';
  svg += '<text class="title-label" x="' + margin.left + '" y="14">Dose-Response · 4PL (Hill) fit · y = B + (T-B)/(1 + 10^((lx-logEC50)*n))</text>';
  svg += '<text class="axis-label" transform="rotate(-90, 14, ' + (margin.top + innerH / 2) + ')" x="14" y="' + (margin.top + innerH / 2) + '" text-anchor="middle">Viability (%)</text>';
  const yTicks = 5;
  for (let i = 0; i <= yTicks; i++) {
    const v = yBot + yRange * i / yTicks;
    const yp = yScale(v);
    svg += '<line class="' + (i === 0 ? 'axis-line' : 'grid-line') + '" x1="' + margin.left + '" y1="' + yp + '" x2="' + (margin.left + innerW) + '" y2="' + yp + '" />';
    svg += '<text class="axis-label" x="' + (margin.left - 6) + '" y="' + (yp + 3) + '" text-anchor="end">' + v.toFixed(0) + '</text>';
  }
  for (let lx = lxMin; lx <= lxMax; lx++) {
    if (lx < lxMin || lx > lxMax) continue;
    const xp = xScale(lx);
    svg += '<line class="grid-line" x1="' + xp + '" y1="' + margin.top + '" x2="' + xp + '" y2="' + (margin.top + innerH) + '" />';
    svg += '<text class="axis-label" x="' + xp + '" y="' + (margin.top + innerH + 14) + '" text-anchor="middle">' + formatConc(Math.pow(10, lx)) + '</text>';
  }
  svg += '<text class="axis-label" x="' + (margin.left + innerW) + '" y="' + (margin.top + innerH + 32) + '" text-anchor="end">log10[concentration]</text>';
  if (fit) {
    if (isFinite(fit.ci95Low) && isFinite(fit.ci95High) && isFinite(fit.seLogEC50)) {
      const loLX = Math.log10(fit.ci95Low), hiLX = Math.log10(fit.ci95High);
      if (loLX >= lxMin && hiLX <= lxMax) {
        const yy1 = yScale(yBot), yy2 = yScale(yTop);
        svg += '<rect class="ec50-band" x="' + xScale(loLX) + '" y="' + yy1 + '" width="' + (xScale(hiLX) - xScale(loLX)) + '" height="' + (yy2 - yy1) + '" />';
      }
    }
    const xp = xScale(fit.logEC50);
    svg += '<line class="ec50-line" x1="' + xp + '" y1="' + margin.top + '" x2="' + xp + '" y2="' + (margin.top + innerH) + '" />';
    svg += '<text class="ec50-label" x="' + (xp + 6) + '" y="' + (margin.top + 12) + '" text-anchor="start">IC50 = ' + formatConc(fit.EC50) + '</text>';
  }
  if (fit) {
    const NP = 100;
    let pathD = '';
    for (let i = 0; i <= NP; i++) {
      const lx = lxMin + (lxMax - lxMin) * i / NP;
      const yh = fit.Bottom + (fit.Top - fit.Bottom) / (1 + Math.pow(10, (lx - fit.logEC50) * fit.HillSlope));
      if (yh < yBot - 50 || yh > yTop + 50) continue;
      const px = xScale(lx), py = yScale(yh);
      pathD += (pathD === '' ? 'M' : 'L') + px.toFixed(2) + ',' + py.toFixed(2) + ' ';
    }
    svg += '<path class="fit-curve" stroke="#0891b2" d="' + pathD + '" />';
  }
  points.forEach(function(p) {
    const px = xScale(p.logConc);
    const py = yScale(p.viability);
    if (isFinite(p.semViab)) {
      const yT = yScale(p.viability + p.semViab);
      const yB = yScale(p.viability - p.semViab);
      const cap = 4;
      svg += '<line class="err-bar" x1="' + px + '" y1="' + yT + '" x2="' + px + '" y2="' + yB + '" />';
      svg += '<line class="err-bar" x1="' + (px - cap) + '" y1="' + yT + '" x2="' + (px + cap) + '" y2="' + yT + '" />';
      svg += '<line class="err-bar" x1="' + (px - cap) + '" y1="' + yB + '" x2="' + (px + cap) + '" y2="' + yB + '" />';
    }
    svg += '<circle class="data-point" cx="' + px + '" cy="' + py + '" r="5" fill="' + p.color + '" />';
    svg += '<text class="data-label" x="' + px + '" y="' + (py + 4) + '" fill="white">' + p.n + '</text>';
  });
  svg += '</svg>';
  return svg;
}
function renderDRParams(result) {
  if (!result.fit) return '';
  const f = result.fit;
  const ciStr = isFinite(f.ci95Low) ? ('[' + formatConc(f.ci95Low) + ', ' + formatConc(f.ci95High) + ']') : '-';
  const items = [
    { name: 'Top', val: f.Top.toFixed(2), unit: '%' },
    { name: 'Bottom', val: f.Bottom.toFixed(2), unit: '%' },
    { name: 'HillSlope', val: f.HillSlope.toFixed(3), unit: '' },
    { name: 'EC50', val: formatConc(f.EC50, 3), unit: '', highlight: true },
    { name: 'logEC50', val: f.logEC50.toFixed(3), unit: 'log10(M)' },
    { name: 'R²', val: f.R2.toFixed(4), unit: '' },
    { name: '95% CI (EC50)', val: ciStr, unit: '' },
    { name: 'n 数据点', val: String(f.nPoints), unit: '' }
  ];
  let html = '<div class="dr-params">';
  items.forEach(function(it) {
    html += '<div class="param' + (it.highlight ? ' highlight' : '') + '">' +
      '<div class="param-name">' + it.name + '</div>' +
      '<div class="param-val">' + it.val + '<span class="param-unit">' + it.unit + '</span></div>' +
      '</div>';
  });
  html += '</div>';
  html += '<p class="hint" style="margin-top:0.3rem;">模型：4PL Hill 方程 · Levenberg-Marquardt · 95% CI：Wald · 拟合输入：' + (state.statsMode === 'currentPlate' ? '当前板数据' : '全部板聚合数据') + '</p>';
  return html;
}
document.getElementById('parse-btn').addEventListener('click', function() {
  const txt = document.getElementById('data-input').value.trim();
  if (!txt) { showToast('请先粘贴数据', 'error'); return; }
  const numbers = txt.split(/[\s,;]+/).filter(function(s) { return s.length > 0; }).map(function(s) { return parseFloat(s); }).filter(function(n) { return !isNaN(n); });
  if (numbers.length === 0) { showToast('未识别出任何数字', 'error'); return; }
  const ap = getActivePlate();
  const layout = plateLayout(ap);
  const max = layout.max, cols = layout.cols, rows = layout.rows;
  if (numbers.length < max) showToast('数据 ' + numbers.length + ' 个，少于 ' + max + ' 孔', 'info');
  if (numbers.length > max) showToast('数据 ' + numbers.length + ' 个，超出 ' + max + '，只取前 ' + max, 'info');
  ap.wellData = {};
  const used = Math.min(numbers.length, max);
  for (let i = 0; i < used; i++) {
    const r = Math.floor(i / cols) + 1;
    const c = (i % cols) + 1;
    ap.wellData[wellIdFromRowCol(r, c)] = numbers[i];
  }
  document.getElementById('parse-status').textContent = '+ 已解析 ' + used + ' 个数据点到「' + ap.name + '」（' + rows + 'x' + cols + '）';
  renderPlate(); renderPlatesBar(); autoSaveDebounced();
});
document.getElementById('load-demo-btn').addEventListener('click', function() {
  const platesData = [
    { name: '板 1 - Day 1', format: 96, assignments: [
      { group: 'Control', wells: ['A1','A2','A3','A4','A5','A6'], mean: 1.00, sd: 0.05 },
      { group: 'DMSO', wells: ['B1','B2','B3','B4','B5','B6'], mean: 1.02, sd: 0.06 },
      { group: 'Drug 1uM', wells: ['C1','C2','C3','C4','C5','C6'], mean: 0.78, sd: 0.07 },
      { group: 'Drug 5uM', wells: ['D1','D2','D3','D4','D5','D6'], mean: 0.52, sd: 0.06 },
      { group: 'Drug 10uM', wells: ['E1','E2','E3','E4','E5','E6'], mean: 0.31, sd: 0.04 }
    ]},
    { name: '板 2 - Day 2', format: 96, assignments: [
      { group: 'Control', wells: ['A1','A2','A3','A4','A5','A6'], mean: 0.98, sd: 0.06 },
      { group: 'DMSO', wells: ['B1','B2','B3','B4','B5','B6'], mean: 1.00, sd: 0.05 },
      { group: 'Drug 1uM', wells: ['C1','C2','C3','C4','C5','C6'], mean: 0.76, sd: 0.08 },
      { group: 'Drug 5uM', wells: ['D1','D2','D3','D4','D5','D6'], mean: 0.50, sd: 0.07 },
      { group: 'Drug 10uM', wells: ['E1','E2','E3','E4','E5','E6'], mean: 0.29, sd: 0.05 }
    ]}
  ];
  state.groups = []; state.activeGroupId = null; state.controlGroupId = null;
  state.plates = [];
  platesData.forEach(function(pd, idx) {
    const pid = genId('p');
    const plateObj = { id: pid, name: pd.name, format: pd.format, wellData: {} };
    for (let r = 1; r <= 8; r++) for (let c = 1; c <= 12; c++) plateObj.wellData[wellIdFromRowCol(r, c)] = +(0.3 + Math.random() * 0.6).toFixed(4);
    pd.assignments.forEach(function(asgn) {
      asgn.wells.forEach(function(w) { plateObj.wellData[w] = +(asgn.mean + (Math.random() - 0.5) * 2 * asgn.sd).toFixed(4); });
    });
    state.plates.push(plateObj);
    pd.assignments.forEach(function(asgn) {
      let g = state.groups.find(function(g2) { return g2.name === asgn.group; });
      if (!g) { g = { id: genId('g'), name: asgn.group, color: COLORS[state.groups.length % COLORS.length], wells: new Set() }; state.groups.push(g); }
      asgn.wells.forEach(function(w) { g.wells.add(makeWellKey(pid, w)); });
    });
  });
  state.controlGroupId = state.groups.find(function(g) { return g.name === 'Control'; }).id;
  state.activePlateId = state.plates[0].id;
  state.selectedWells.clear();
  document.getElementById('parse-status').textContent = '+ 已载入 2 个示例板 + 5 个跨板分组';
  renderPlatesBar(); renderPlate(); renderGroups(); renderGroupSelect(); renderResults();
  showToast('已载入 2 板示例，Control 已设为对照', 'success');
  autoSaveDebounced();
});
document.getElementById('plate-name-input').addEventListener('change', function(e) {
  const ap = getActivePlate();
  const v = e.target.value.trim();
  if (v) { ap.name = v.slice(0, 32); renderPlatesBar(); renderResults(); autoSaveDebounced(); }
});
document.getElementById('plate-format-select').addEventListener('change', function(e) {
  const ap = getActivePlate();
  const newFmt = parseInt(e.target.value, 10);
  if (newFmt === ap.format) return;
  if (Object.keys(ap.wellData).length > 0 && !confirm('切换格式将清空当前板「' + ap.name + '」的数据，确定吗？')) {
    e.target.value = ap.format;
    return;
  }
  ap.format = newFmt;
  ap.wellData = {};
  const ap_id = ap.id;
  state.groups.forEach(function(g) {
    const toDel = [];
    g.wells.forEach(function(k) { if (parsePlateFromKey(k) === ap_id) toDel.push(k); });
    toDel.forEach(function(k) { g.wells.delete(k); });
  });
  document.getElementById('parse-status').textContent = '';
  renderPlatesBar(); renderPlate(); renderGroups(); renderGroupSelect(); renderResults(); autoSaveDebounced();
});
document.querySelectorAll('input[name="chartMode"]').forEach(function(r) {
  r.addEventListener('change', function(e) { state.chartMode = e.target.value; renderResults(); autoSaveDebounced(); });
});
document.querySelectorAll('input[name="statsMode"]').forEach(function(r) {
  r.addEventListener('change', function(e) { state.statsMode = e.target.value; renderResults(); autoSaveDebounced(); });
});
document.getElementById('new-group-btn').addEventListener('click', function() {
  const g = { id: genId('g'), name: '组 ' + (state.groups.length + 1), color: COLORS[state.groups.length % COLORS.length], wells: new Set() };
  state.groups.push(g);
  state.activeGroupId = g.id;
  renderGroups(); renderGroupSelect(); renderResults(); renderPlatesBar(); autoSaveDebounced();
});
document.getElementById('assign-btn').addEventListener('click', function() {
  if (!state.activeGroupId) { showToast('请先选择或创建一个分组', 'error'); return; }
  if (state.selectedWells.size === 0) { showToast('请先选择孔位（点击当前板）', 'error'); return; }
  const group = state.groups.find(function(g) { return g.id === state.activeGroupId; });
  let replaced = 0;
  state.selectedWells.forEach(function(k) {
    state.groups.forEach(function(g) {
      if (g.id !== group.id && g.wells.has(k)) { g.wells.delete(k); replaced++; }
    });
    group.wells.add(k);
  });
  renderPlate(); renderGroups(); renderResults(); autoSaveDebounced();
  showToast('已分配 ' + state.selectedWells.size + ' 个孔到「' + group.name + '」' + (replaced ? '（' + replaced + ' 个从其他组移出）' : ''));
});
document.getElementById('unassign-btn').addEventListener('click', function() {
  if (state.selectedWells.size === 0) { showToast('请先选择孔位', 'error'); return; }
  let removed = 0;
  state.selectedWells.forEach(function(k) {
    state.groups.forEach(function(g) {
      if (g.wells.has(k)) { g.wells.delete(k); removed++; }
    });
  });
  renderPlate(); renderGroups(); renderResults(); autoSaveDebounced();
  showToast(removed > 0 ? ('已从分组移除 ' + removed + ' 个孔') : '选中孔未分配任何分组', removed > 0 ? 'info' : 'success');
});
document.getElementById('deselect-btn').addEventListener('click', function() {
  state.selectedWells.clear();
  state.lastClickedWell = null;
  renderPlate(); renderGroups();
});
document.getElementById('clear-wells-btn').addEventListener('click', function() {
  if (state.groups.length === 0) return;
  if (!confirm('清空所有分组中的孔位？数据本身保留。')) return;
  state.groups.forEach(function(g) { g.wells.clear(); });
  state.controlGroupId = null;
  renderPlate(); renderGroups(); renderGroupSelect(); renderResults(); renderPlatesBar(); autoSaveDebounced();
});
document.getElementById('control-select').addEventListener('change', function(e) {
  state.controlGroupId = e.target.value || null;
  renderGroups(); renderResults(); autoSaveDebounced();
});
function triggerDownload(blob, filename) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url; a.download = filename;
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
  setTimeout(function() { URL.revokeObjectURL(url); }, 100);
  showToast('已导出：' + filename, 'success');
}
document.getElementById('export-csv-btn').addEventListener('click', function() {
  if (state.groups.length === 0) { showToast('没有数据可导出', 'error'); return; }
  const rows = [['Plate', 'Group', 'Well', 'Absorbance']];
  state.groups.forEach(function(g) {
    Array.from(g.wells).sort().forEach(function(k) {
      const pid = parsePlateFromKey(k), wid = parseWellFromKey(k);
      const plate = state.plates.find(function(p) { return p.id === pid; });
      const plateName = plate ? plate.name : pid;
      const v = plate ? plate.wellData[wid] : undefined;
      rows.push([plateName, g.name, wid, v !== undefined ? v : '']);
    });
  });
  const stats = computeStats();
  if (stats.hasControl) {
    rows.push([]);
    rows.push(['=== Summary (' + (state.statsMode === 'currentPlate' ? 'Current Plate' : 'Aggregated All Plates') + ') ===']);
    rows.push(['Group', 'n', 'Mean', 'SD', 'SEM', 'CV(%)', 'Viability(%)', 'p_adjusted', 'Significance']);
    stats.groups.forEach(function(s) {
      const isCtrl = s.group.id === state.controlGroupId;
      rows.push([s.group.name, s.n, s.mean.toFixed(4), s.sd.toFixed(4), s.sem.toFixed(4), s.cv.toFixed(2), s.viability.toFixed(2), isCtrl ? '' : s.pAdjusted.toFixed(4), s.significance]);
    });
  }
  const csv = rows.map(function(r) {
    return r.map(function(c) {
      const s = String(c);
      return /[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
    }).join(',');
  }).join('\r\n');
  const blob = new Blob(['\ufeff' + csv], { type: 'text/csv;charset=utf-8' });
  triggerDownload(blob, 'mtt-multi-' + new Date().toISOString().slice(0, 10) + '.csv');
});
document.getElementById('export-png-btn').addEventListener('click', function() {
  const svg = document.querySelector('.bar-chart');
  if (!svg) { showToast('请先生成图表', 'error'); return; }
  const svgData = new XMLSerializer().serializeToString(svg);
  const canvas = document.createElement('canvas');
  const W = 1400, H = 800;
  canvas.width = W; canvas.height = H;
  const ctx = canvas.getContext('2d');
  ctx.fillStyle = 'white'; ctx.fillRect(0, 0, W, H);
  const img = new Image();
  img.onload = function() {
    ctx.drawImage(img, 0, 0, W, H);
    canvas.toBlob(function(blob) { triggerDownload(blob, 'mtt-chart-' + new Date().toISOString().slice(0, 10) + '.png'); }, 'image/png');
  };
  img.src = 'data:image/svg+xml;base64,' + btoa(unescape(encodeURIComponent(svgData)));
});

function serializeState() {
  return {
    plates: state.plates.map(function(p) {
      return { id: p.id, name: p.name, format: p.format, wellData: Object.assign({}, p.wellData) };
    }),
    activePlateId: state.activePlateId,
    groups: state.groups.map(function(g) {
      return { id: g.id, name: g.name, color: g.color, _concRaw: g._concRaw || null, wells: Array.from(g.wells) };
    }),
    activeGroupId: state.activeGroupId,
    controlGroupId: state.controlGroupId,
    chartMode: state.chartMode,
    statsMode: state.statsMode
  };
}
function deserializeState(data) {
  if (!data) return;
  state.plates = (data.plates || []).map(function(p) {
    return { id: p.id, name: p.name, format: p.format, wellData: Object.assign({}, p.wellData || {}) };
  });
  if (state.plates.length === 0) state.plates = [{ id: genId('p'), name: '板 1', format: 96, wellData: {} }];
  state.activePlateId = data.activePlateId || state.plates[0].id;
  state.groups = (data.groups || []).map(function(g) {
    return { id: g.id, name: g.name, color: g.color, _concRaw: g._concRaw || null, wells: new Set(g.wells || []) };
  });
  state.activeGroupId = data.activeGroupId || null;
  state.controlGroupId = data.controlGroupId || null;
  state.chartMode = data.chartMode || 'abs';
  state.statsMode = data.statsMode || 'aggregate';
  state.selectedWells.clear();
}
function saveStateAuto() {
  try {
    const data = serializeState();
    localStorage.setItem('mtt-auto', JSON.stringify({ ts: Date.now(), data: data }));
  } catch (e) { /* ignore quota */ }
  updateSessionInfo();
}
function loadStateAuto() {
  try {
    const raw = localStorage.getItem('mtt-auto');
    if (!raw) return null;
    const parsed = JSON.parse(raw);
    return (parsed && parsed.data && parsed.data.plates) ? parsed : null;
  } catch (e) { return null; }
}
function saveSnapshot(name) {
  try {
    const data = serializeState();
    const snaps = loadSnapshots();
    snaps.unshift({ id: 's' + Date.now().toString(36), name: name, ts: Date.now(), data: data });
    if (snaps.length > 6) snaps.length = 6;
    localStorage.setItem('mtt-snapshots', JSON.stringify(snaps));
    return true;
  } catch (e) { return false; }
}
function loadSnapshots() {
  try { return JSON.parse(localStorage.getItem('mtt-snapshots') || '[]'); }
  catch (e) { return []; }
}
function deleteSnapshot(id) {
  try {
    const snaps = loadSnapshots().filter(function(s) { return s.id !== id; });
    localStorage.setItem('mtt-snapshots', JSON.stringify(snaps));
  } catch (e) {}
}
let saveTimeout = null;
function autoSaveDebounced() {
  if (saveTimeout) clearTimeout(saveTimeout);
  saveTimeout = setTimeout(saveStateAuto, 700);
}
function updateSessionInfo() {
  const el = document.getElementById('session-info-text');
  if (!el) return;
  const auto = loadStateAuto();
  if (!auto) { el.textContent = '新会话（未自动保存）'; return; }
  const d = new Date(auto.ts);
  const ts = d.getHours().toString().padStart(2, '0') + ':' + d.getMinutes().toString().padStart(2, '0') + ':' + d.getSeconds().toString().padStart(2, '0');
  el.textContent = '自动保存于 ' + ts + ' · ' + auto.data.plates.length + '板/' + auto.data.groups.length + '组';
}
function showSnapshotsModal() {
  const back = document.createElement('div');
  back.className = 'snapshot-modal-back';
  const modal = document.createElement('div');
  modal.className = 'snapshot-modal';
  const auto = loadStateAuto();
  let html = '<h3>会话历史</h3>';
  if (auto) {
    const d = new Date(auto.ts);
    html += '<div class="snapshot-section-title">自动保存</div>';
    html += '<div class="snapshot-item">' +
      '<span class="name">自动保存</span>' +
      '<span class="ts">' + d.toLocaleString() + '</span>' +
      '<button class="secondary tiny" data-act="restore-auto">恢复</button>' +
      '</div>';
  }
  const snaps = loadSnapshots();
  if (snaps.length === 0) {
    html += '<div class="snapshot-section-title">命名快照</div>';
    html += '<div class="snapshot-empty">暂无命名快照。点击"保存快照"创建。</div>';
  } else {
    html += '<div class="snapshot-section-title">命名快照</div>';
    snaps.forEach(function(s) {
      const d = new Date(s.ts);
      html += '<div class="snapshot-item" data-id="' + s.id + '">' +
        '<span class="name">' + escapeHtml(s.name) + '</span>' +
        '<span class="ts">' + d.toLocaleString() + '</span>' +
        '<button class="secondary tiny" data-act="restore" data-id="' + s.id + '">恢复</button>' +
        '<button class="danger tiny" data-act="delete" data-id="' + s.id + '">删除</button>' +
        '</div>';
    });
  }
  html += '<div style="margin-top:0.85rem; display:flex; gap:0.4rem; justify-content:flex-end;">' +
    '<button class="outline tiny" data-act="close">关闭</button>' +
    '</div>';
  modal.innerHTML = html;
  document.body.appendChild(back);
  document.body.appendChild(modal);
  const close = function() { back.remove(); modal.remove(); };
  back.addEventListener('click', close);
  modal.addEventListener('click', function(e) {
    const a = e.target.dataset.act;
    if (!a) return;
    if (a === 'close') { close(); return; }
    if (a === 'restore-auto') {
      const a2 = loadStateAuto();
      if (a2) { deserializeState(a2.data); renderAll(); updateSessionInfo(); showToast('已恢复自动保存', 'success'); }
      close();
    } else {
      const id = e.target.dataset.id;
      const s = loadSnapshots().find(function(s) { return s.id === id; });
      if (!s) return;
      if (a === 'restore') { deserializeState(s.data); renderAll(); updateSessionInfo(); showToast('已恢复: ' + s.name, 'success'); close(); }
      else if (a === 'delete') {
        if (!confirm('删除快照「' + s.name + '」？')) return;
        deleteSnapshot(id);
        close();
        showSnapshotsModal();
      }
    }
  });
}
function renderAll() {
  renderPlatesBar(); renderPlate(); renderGroups(); renderGroupSelect(); renderResults(); renderDoseResponse();
}
function startNewSession() {
  if (state.plates.length > 0 && state.plates[0].wellData && Object.keys(state.plates[0].wellData).length > 0 || state.groups.length > 0) {
    saveSnapshot('会话_' + new Date().toLocaleDateString() + '_' + Date.now().toString(36).slice(-4));
  }
  if (!confirm('开始新会话将清空所有板、分组、数据。\n（当前状态已自动保存为快照）\n\n确定继续吗？')) return;
  state.plates = [{ id: genId('p'), name: '板 1', format: 96, wellData: {} }];
  state.activePlateId = state.plates[0].id;
  state.groups = []; state.activeGroupId = null; state.controlGroupId = null;
  state.selectedWells.clear();
  renderAll();
  document.getElementById('parse-status').textContent = '';
  document.getElementById('data-input').value = '';
  autoSaveDebounced();
  showToast('新会话已开始', 'success');
}
document.getElementById('save-snapshot-btn').addEventListener('click', function() {
  const def = '快照_' + new Date().toLocaleDateString() + '_' + Date.now().toString(36).slice(-4);
  const name = prompt('为快照起个名字：', def);
  if (!name) return;
  if (saveSnapshot(name)) { showToast('已保存: ' + name, 'success'); updateSessionInfo(); }
  else showToast('保存失败', 'error');
});
document.getElementById('snapshots-btn').addEventListener('click', showSnapshotsModal);
document.getElementById('restore-auto-btn').addEventListener('click', function() {
  const a = loadStateAuto();
  if (!a) { showToast('没有自动保存的会话', 'info'); return; }
  deserializeState(a.data);
  renderAll();
  updateSessionInfo();
  showToast('已恢复上次自动保存', 'success');
});
document.getElementById('new-session-btn').addEventListener('click', startNewSession);



function renderHeatmap(stats) {
  const ap = getActivePlate();
  const layout = plateLayout(ap);
  const rows = layout.rows, cols = layout.cols;
  const cellSize = Math.max(28, Math.min(48, Math.floor(720 / (cols + 2))));
  const headerSize = 22;
  const margin = { top: 30, right: 24, bottom: 50, left: 30 };
  const width = margin.left + headerSize + cols * cellSize + margin.right;
  const height = margin.top + headerSize + rows * cellSize + margin.bottom + 30;

  let ctrlMean = 0;
  if (stats.control && stats.control.mean > 0) {
    ctrlMean = stats.control.mean;
  } else {
    const allMeans = [];
    for (let i = 0; i < state.groups.length; i++) {
      const g = state.groups[i];
      const v = state.statsMode === "currentPlate" ? getGroupValuesForPlate(g, state.activePlateId) : getGroupValues(g);
      if (v.length > 0) allMeans.push(mean(v));
    }
    ctrlMean = allMeans.length > 0 ? Math.max.apply(null, allMeans) : 1;
  }
  if (ctrlMean <= 0) ctrlMean = 1;

  function colorFor(v) {
    if (v == null || isNaN(v)) return { fill: "#f1f5f9", text: "#94a3b8" };
    const cv = Math.max(0, Math.min(130, v));
    let r, g, b;
    if (cv < 50) {
      const t = cv / 50;
      r = 220; g = 38 + (202 - 38) * t; b = 38;
    } else if (cv < 100) {
      const t = (cv - 50) / 50;
      r = 220 - (220 - 22) * t; g = 202; b = 38 + (94 - 38) * t;
    } else {
      const t = Math.min(1, (cv - 100) / 30);
      r = Math.max(0, 22 - 22 * t); g = 202 - (202 - 150) * t; b = 94;
    }
    const dark = cv < 50 || cv > 110;
    return { fill: "rgb(" + Math.round(r) + "," + Math.round(g) + "," + Math.round(b) + ")", text: dark ? "#f8fafc" : "#1e293b" };
  }

  let svg = '<svg class="heatmap" viewBox="0 0 ' + width + ' ' + height + '" preserveAspectRatio="xMidYMid meet">';
  svg += '<text class="title-label" x="' + margin.left + '" y="14">孔板热图（viability % 色阶 · 0%=红 · 100%=绿）</text>';

  for (let c = 1; c <= cols; c++) {
    const x = margin.left + headerSize + (c - 1) * cellSize + cellSize / 2;
    svg += '<text class="axis-label" x="' + x + '" y="' + (margin.top - 6) + '" text-anchor="middle" style="font-weight:600;">' + c + '</text>';
  }
  for (let r = 1; r <= rows; r++) {
    const y = margin.top + headerSize + (r - 1) * cellSize + cellSize / 2 + 4;
    svg += '<text class="axis-label" x="' + (margin.left - 6) + '" y="' + y + '" text-anchor="end" style="font-weight:600;">' + String.fromCharCode(64 + r) + '</text>';
  }

  let qcWells = [];
  for (let r = 1; r <= rows; r++) {
    for (let c = 1; c <= cols; c++) {
      const wid = String.fromCharCode(64 + r) + c;
      const key = makeWellKey(ap.id, wid);
      const data = ap.wellData[wid];
      const group = state.groups.find(function(g) { return g.wells.has(key); });

      const x = margin.left + headerSize + (c - 1) * cellSize;
      const y = margin.top + headerSize + (r - 1) * cellSize;

      let viability = null;
      let color = { fill: "#f1f5f9", text: "#94a3b8" };
      let textVal = "";
      if (data !== undefined && data !== null && data !== "") {
        viability = (data / ctrlMean) * 100;
        color = colorFor(viability);
        textVal = viability.toFixed(0);
        if (group && group.id === state.controlGroupId && (viability < 70 || viability > 110)) qcWells.push(wid);
      }

      svg += '<rect x="' + x + '" y="' + y + '" width="' + (cellSize - 2) + '" height="' + (cellSize - 2) + '" fill="' + color.fill + '" stroke="' + (group ? group.color : "#cbd5e1") + '" stroke-width="' + (group ? "2.5" : "0.5") + '" rx="3" />';
      if (textVal) {
        svg += '<text class="bar-label" x="' + (x + cellSize / 2 - 1) + '" y="' + (y + cellSize / 2 + 4) + '" text-anchor="middle" fill="' + color.text + '">' + textVal + '</text>';
      }
      if (group) {
        const gName = group.name.length > 5 ? group.name.substring(0, 4) + ".." : group.name;
        svg += '<text class="axis-label" x="' + (x + cellSize / 2 - 1) + '" y="' + (y + cellSize - 4) + '" text-anchor="middle" fill="' + color.text + '" style="font-size:7px;opacity:0.7;">' + escapeHtml(gName) + '</text>';
      }
    }
  }

  if (qcWells.length > 0) {
    svg += '<rect class="qc-band" x="' + margin.left + '" y="' + (margin.top + headerSize + rows * cellSize + 8) + '" width="' + (width - 2 * margin.left) + '" height="22" rx="3" />';
    svg += '<text class="qc-text" x="' + (margin.left + 8) + '" y="' + (margin.top + headerSize + rows * cellSize + 22) + '" style="font-size:11px;">⚠ QC 报警：对照孔 (' + qcWells.join(", ") + ") viability 偏离 [70, 110] 范围</text>";
  }

  const legY = height - 14;
  const legX = margin.left;
  const legW = 22;
  const items = [{ v: 0, label: "0%" }, { v: 25, label: "25%" }, { v: 50, label: "50%" }, { v: 75, label: "75%" }, { v: 100, label: "100%" }, { v: 130, label: "130%" }];
  for (let i = 0; i < items.length; i++) {
    const c = colorFor(items[i].v);
    svg += '<rect x="' + (legX + i * (legW + 4)) + '" y="' + legY + '" width="' + legW + '" height="12" fill="' + c.fill + '" stroke="#cbd5e1" />';
    svg += '<text class="axis-label" x="' + (legX + i * (legW + 4) + legW / 2) + '" y="' + (legY + 24) + '" text-anchor="middle">' + items[i].label + '</text>';
  }
  svg += '</svg>';
  return svg;
}

const _restored = loadStateAuto();
if (_restored) deserializeState(_restored.data);
renderPlatesBar();
renderPlate();
renderGroups();
renderGroupSelect();
renderResults();
renderDoseResponse();
updateSessionInfo();

