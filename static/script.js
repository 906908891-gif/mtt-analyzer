

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
  // Uses the Numerical Recipes continued-fraction algorithm.
  // Reference: Press et al. (1992), Numerical Recipes in C, sec. 6.4.
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  // Use symmetry: I_x(a, b) = 1 - I_(1-x)(b, a)
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
  // Continued fraction for incomplete beta (NR Press et al.)
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
  // Step-down correction with running MAX (correct direction).
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
