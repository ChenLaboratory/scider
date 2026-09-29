#define R_NO_REMAP
#include "utils.h"
#include "rand.h"
#include <vector>
#include <algorithm>

// Random labelling permutation test for corDensity().

namespace {

const int BLOCK = 16;

struct Problem {
  int nct, ngrid, nroi, ngroup, npair, nstat, ncell;
  const int *goff;    // group offsets into the cells (ngroup + 1)
  const int *labs;    // cell type of the coi cells, in cell order
  const int *loff;    // group offsets into labs (ngroup + 1)
  const int *roi;     // ROI of each tested grid (0 based)
  std::vector<int> pool0;          // lattice slot of every cell
  std::vector<int> roi_n, pc, pd;
  // separable kernel
  int cs, nslot, nh, nv;
  std::vector<double> gx, gy, scal;
  std::vector<int> h_slot, h_tap, h_ntap, v_idx;
};

uint64_t splitmix64(uint64_t z) {
  z += 0x9e3779b97f4a7c15ULL;
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

// Unbiased integer in [0, range) (Lemire's multiply and reject method).
uint32_t bounded(pcg32_random_t *rng, uint32_t range) {
  uint64_t m = (uint64_t)pcg32_random_r(rng) * range;
  uint32_t l = (uint32_t)m;
  if (l < range) {
    uint32_t t = (0u - range) % range;
    while (l < t) {
      m = (uint64_t)pcg32_random_r(rng) * range;
      l = (uint32_t)m;
    }
  }
  return (uint32_t)(m >> 32);
}

// Weighted sums over n rows of a row major array with ncol columns
template <int NC>
inline void wsum(const double *w, int wstep, const double *src,
                 const int *idx, int n, int ncol, int c0, double *dst,
                 size_t dstep) {
  double a[NC] = {0.0};
  for (int t = 0; t < n; t++) {
    const double wt = w[(size_t)t * wstep];
    const double *row = src + (size_t)(idx ? idx[t] : t) * ncol + c0;
    for (int c = 0; c < NC; c++) a[c] += wt * row[c];
  }
  for (int c = 0; c < NC; c++) dst[(size_t)(c0 + c) * dstep] = a[c];
}

inline void wsumAll(const double *w, int wstep, const double *src,
                    const int *idx, int n, int ncol, double *dst,
                    size_t dstep) {
  int c0 = 0;
  for (; c0 + 4 <= ncol; c0 += 4) {
    wsum<4>(w, wstep, src, idx, n, ncol, c0, dst, dstep);
  }
  switch (ncol - c0) {
  case 3: wsum<3>(w, wstep, src, idx, n, ncol, c0, dst, dstep); break;
  case 2: wsum<2>(w, wstep, src, idx, n, ncol, c0, dst, dstep); break;
  case 1: wsum<1>(w, wstep, src, idx, n, ncol, c0, dst, dstep); break;
  default: break;
  }
}

// Kernel smoothing of per slot values cnt (nslot x ncol, row major) to the
// tested grids, D (ngrid x ncol, column major), using the workspace H.
void smooth(const Problem &P, int ncol, const double *cnt, double *H,
            double *D) {
  for (int h = 0; h < P.nh; h++) {
    wsumAll(P.gx.data() + P.h_tap[h], P.cs, cnt + (size_t)P.h_slot[h] * ncol,
            NULL, P.h_ntap[h], ncol, H + (size_t)h * ncol, 1);
  }
  for (int i = 0; i < P.ngrid; i++) {
    wsumAll(P.gy.data(), 1, H, P.v_idx.data() + (size_t)i * P.nv, P.nv, ncol,
            D + i, P.ngrid);
  }
}

// Sort key of a grid
struct KeyIdx {
  uint64_t k;
  int i;
};

inline uint64_t sortKey(double v) {
  uint64_t u;
  std::memcpy(&u, &v, sizeof u);
  return (u >> 63) ? ~u : (u | 0x8000000000000000ULL);
}

// D, W and U are column major (ngrid x nct)
struct Workspace {
  std::vector<int> pool, jrec, ord, tmp, seen;
  std::vector<KeyIdx> ki, ki2;
  std::vector<std::vector<uint32_t> > hist;
  std::vector<double> cnt, H, D, W, U, rsum, accW, ssW, accU, ssU, wl, ul,
    stat;

  explicit Workspace(const Problem &P) :
    pool(P.pool0), jrec(P.ncell), ord(P.ngrid), tmp(P.nroi), seen(P.nroi),
    ki(P.ngrid), ki2(P.ngrid),
    cnt((size_t)P.nslot * P.nct, 0.0), H((size_t)P.nh * P.nct),
    D((size_t)P.ngrid * P.nct), W((size_t)P.ngrid * P.nct),
    U((size_t)P.ngrid * P.nct), rsum(P.nroi),
    accW((size_t)P.nroi * P.npair), ssW((size_t)P.nroi * P.nct),
    accU(P.npair), ssU(P.nct), wl(P.nct), ul(P.nct), stat(P.nstat) {}
};

// Sorts a[0..n) by key, using buf[0..n) as scratch
void sortKeys(Workspace &ws, KeyIdx *a, KeyIdx *buf, int n, int depth) {
  if (n <= 32) {
    for (int t = 1; t < n; t++) {
      KeyIdx x = a[t];
      int u = t;
      while (u > 0 && a[u - 1].k > x.k) {
        a[u] = a[u - 1];
        u--;
      }
      a[u] = x;
    }
    return;
  }
  uint64_t kmin = a[0].k, kmax = a[0].k;
  for (int t = 1; t < n; t++) {
    kmin = std::min(kmin, a[t].k);
    kmax = std::max(kmax, a[t].k);
  }
  const uint64_t range = kmax - kmin;
  if (range == 0) return;
  int bits = 0, B = 1;
  while (bits < 64 && (range >> bits) != 0) bits++;
  while ((1 << B) < n && B < 16) B++;
  const int shift = bits > B ? bits - B : 0;
  const int nb = (int)(range >> shift) + 1;
  if ((int)ws.hist.size() <= depth) ws.hist.resize(depth + 1);
  std::vector<uint32_t> &hv = ws.hist[depth];
  if ((int)hv.size() < nb + 1) hv.resize(nb + 1);
  uint32_t *h = hv.data();
  std::fill(h, h + nb + 1, 0u);
  for (int t = 0; t < n; t++) h[((a[t].k - kmin) >> shift) + 1]++;
  for (int b = 0; b < nb; b++) h[b + 1] += h[b];
  for (int t = 0; t < n; t++) buf[h[(a[t].k - kmin) >> shift]++] = a[t];
  std::copy(buf, buf + n, a);
  // bucket b now spans [h[b - 1], h[b]) (h[-1] = 0)
  uint32_t lo = 0;
  for (int b = 0; b < nb; b++) {
    const uint32_t hi = h[b];
    if (hi - lo > 1) sortKeys(ws, a + lo, buf + lo, (int)(hi - lo), depth + 1);
    lo = hi;
  }
}

// Densities of one relabelling (permute = false: the observed labelling).
void density(const Problem &P, Workspace &ws, pcg32_random_t *rng,
             bool permute) {
  for (int g = 0; g < P.ngroup; g++) {
    int *pl = ws.pool.data() + P.goff[g];
    const int n = P.goff[g + 1] - P.goff[g];
    const int *lb = P.labs + P.loff[g];
    const int m = P.loff[g + 1] - P.loff[g];
    int *jr = ws.jrec.data() + P.goff[g];
    if (permute) {
      for (int k = 0; k < m; k++) {
        int j = k + (int)bounded(rng, (uint32_t)(n - k));
        jr[k] = j;
        std::swap(pl[k], pl[j]);
      }
    }
    for (int k = 0; k < m; k++) ws.cnt[(size_t)pl[k] * P.nct + lb[k]] += 1.0;
  }
  smooth(P, P.nct, ws.cnt.data(), ws.H.data(), ws.D.data());
  for (int c = 0; c < P.nct; c++) {
    double *Dc = ws.D.data() + (size_t)c * P.ngrid;
    for (int i = 0; i < P.ngrid; i++) Dc[i] *= P.scal[i];
  }
  for (int g = 0; g < P.ngroup; g++) {
    int *pl = ws.pool.data() + P.goff[g];
    const int *lb = P.labs + P.loff[g];
    const int m = P.loff[g + 1] - P.loff[g];
    const int *jr = ws.jrec.data() + P.goff[g];
    for (int k = 0; k < m; k++) ws.cnt[(size_t)pl[k] * P.nct + lb[k]] = 0.0;
    if (permute) {
      for (int k = m - 1; k >= 0; k--) std::swap(pl[k], pl[jr[k]]);
    }
  }
}

// Spearman statistics of all pairs
void statistics(const Problem &P, Workspace &ws) {
  const int nct = P.nct, ngrid = P.ngrid;
  for (int c = 0; c < nct; c++) {
    const double *col = ws.D.data() + (size_t)c * ngrid;
    double *Wc = ws.W.data() + (size_t)c * ngrid;
    double *Uc = ws.U.data() + (size_t)c * ngrid;

    // order the grids by density: zeros need no sorting
    int nzero = 0, nnz = 0, nneg = 0;
    for (int i = 0; i < ngrid; i++) {
      const double v = col[i];
      if (v == 0.0) {
        ws.ord[nzero++] = i;
      } else {
        ws.ki[nnz].k = sortKey(v);
        ws.ki[nnz++].i = i;
        nneg += (v < 0.0);
      }
    }
    sortKeys(ws, ws.ki.data(), ws.ki2.data(), nnz, 0);
    if (nneg > 0) {
      // negative densities cannot arise from the kernel, kept for safety
      std::copy_backward(ws.ord.begin(), ws.ord.begin() + nzero,
                         ws.ord.begin() + nneg + nzero);
      for (int t = 0; t < nneg; t++) ws.ord[t] = ws.ki[t].i;
    }
    for (int t = nneg; t < nnz; t++) ws.ord[nzero + t] = ws.ki[t].i;

    // walk the tie groups
    std::fill(ws.seen.begin(), ws.seen.end(), 0);
    std::fill(ws.rsum.begin(), ws.rsum.end(), 0.0);
    int s = 0;
    while (s < ngrid) {
      const double v = col[ws.ord[s]];
      int e = s + 1;
      while (e < ngrid && col[ws.ord[e]] == v) e++;
      const double gr = 0.5 * (s + 1 + e);
      if (e == s + 1) {
        const int i = ws.ord[s], r = P.roi[i];
        Uc[i] = gr;
        ws.rsum[r] += gr;
        Wc[i] = ws.seen[r] + 1.0 - 0.5 * (P.roi_n[r] + 1);
        ws.seen[r]++;
      } else {
        for (int t = s; t < e; t++) {
          const int i = ws.ord[t];
          Uc[i] = gr;
          ws.rsum[P.roi[i]] += gr;
          ws.tmp[P.roi[i]]++;
        }
        for (int t = s; t < e; t++) {
          const int i = ws.ord[t], r = P.roi[i];
          Wc[i] = ws.seen[r] + 0.5 * (ws.tmp[r] + 1) - 0.5 * (P.roi_n[r] + 1);
        }
        for (int t = s; t < e; t++) {
          const int r = P.roi[ws.ord[t]];
          ws.seen[r] += ws.tmp[r];
          ws.tmp[r] = 0;
        }
      }
      s = e;
    }
    for (int i = 0; i < ngrid; i++) {
      const int r = P.roi[i];
      Uc[i] -= ws.rsum[r] / P.roi_n[r];
    }
  }

  // cross products
  std::fill(ws.accW.begin(), ws.accW.end(), 0.0);
  std::fill(ws.ssW.begin(), ws.ssW.end(), 0.0);
  std::fill(ws.accU.begin(), ws.accU.end(), 0.0);
  std::fill(ws.ssU.begin(), ws.ssU.end(), 0.0);
  double *w = ws.wl.data(), *u = ws.ul.data();
  for (int i = 0; i < ngrid; i++) {
    const int r = P.roi[i];
    for (int c = 0; c < nct; c++) {
      w[c] = ws.W[(size_t)c * ngrid + i];
      u[c] = ws.U[(size_t)c * ngrid + i];
    }
    double *sw = ws.ssW.data() + (size_t)r * nct;
    double *aw = ws.accW.data() + (size_t)r * P.npair;
    for (int c = 0; c < nct; c++) {
      sw[c] += w[c] * w[c];
      ws.ssU[c] += u[c] * u[c];
    }
    for (int p = 0; p < P.npair; p++) {
      aw[p] += w[P.pc[p]] * w[P.pd[p]];
      ws.accU[p] += u[P.pc[p]] * u[P.pd[p]];
    }
  }

  auto corr = [](double a, double s1, double s2) {
    if (!(s1 > 0.0) || !(s2 > 0.0)) return NA_REAL;
    double r = a / std::sqrt(s1 * s2);
    return r > 1.0 ? 1.0 : (r < -1.0 ? -1.0 : r);
  };
  for (int p = 0; p < P.npair; p++) {
    ws.stat[p] = corr(ws.accU[p], ws.ssU[P.pc[p]], ws.ssU[P.pd[p]]);
  }
  for (int r = 0; r < P.nroi; r++) {
    const double *sw = ws.ssW.data() + (size_t)r * nct;
    for (int p = 0; p < P.npair; p++) {
      ws.stat[P.npair + (size_t)r * P.npair + p] =
        corr(ws.accW[(size_t)r * P.npair + p], sw[P.pc[p]], sw[P.pd[p]]);
    }
  }
}

// Lattice set up. Node coordinates a (columns) and r (rows) start at 0.
void setupLattice(Problem &P, const int *cell_a, const int *cell_r,
                  const int *grid_a, const int *grid_r, int amax, int rmax,
                  int ra, int rr, int parity) {
  const int wa = amax + 2 * ra + 1, hr = rmax + 2 * rr + 1;
  const int wc = (wa + P.cs - 1) / P.cs;
  P.nslot = wc * hr;
  P.nv = 2 * rr + 1;
  P.pool0.resize(P.ncell);
  for (int k = 0; k < P.ncell; k++) {
    P.pool0[k] = (cell_r[k] + rr) * wc + (cell_a[k] + ra) / P.cs;
  }
  // positions of the horizontal pass: column of a tested grid, rows within rr
  std::vector<int> mark((size_t)wa * hr, -1);
  P.v_idx.resize((size_t)P.ngrid * P.nv);
  P.nh = 0;
  for (int i = 0; i < P.ngrid; i++) {
    const int aa = grid_a[i] + ra;
    for (int k = 0; k < P.nv; k++) {
      const int row = grid_r[i] + k;          // padded row index
      int &id = mark[(size_t)row * wa + aa];
      if (id < 0) {
        id = P.nh++;
        // nodes of this row within ra of aa, in padded column units
        int first = aa - ra;
        if (P.cs == 2) {
          const int r0 = row - rr;             // lattice row
          const int node_par = ((parity - r0 + ra) % 2 + 2) % 2;
          if (((first % 2) + 2) % 2 != node_par) first++;
        }
        P.h_slot.push_back(row * wc + first / P.cs);
        P.h_tap.push_back(first - (aa - ra));
        P.h_ntap.push_back((aa + ra - first) / P.cs + 1);
      }
      P.v_idx[(size_t)i * P.nv + k] = id;
    }
  }
}

} // namespace

extern "C" {
  SEXP C_corDensityPerm(SEXP prob, SEXP obs_, SEXP acc_, SEXP b_start_,
                        SEXP nb_, SEXP seed_, SEXP cpu_threads) {
    Problem P;
    P.goff = INTEGER(VECTOR_ELT(prob, 0));
    P.labs = INTEGER(VECTOR_ELT(prob, 1));
    P.loff = INTEGER(VECTOR_ELT(prob, 2));
    const int *cell_a = INTEGER(VECTOR_ELT(prob, 3));
    const int *cell_r = INTEGER(VECTOR_ELT(prob, 4));
    const int *grid_a = INTEGER(VECTOR_ELT(prob, 5));
    const int *grid_r = INTEGER(VECTOR_ELT(prob, 6));
    const int *lat = INTEGER(VECTOR_ELT(prob, 7));
    SEXP gx = VECTOR_ELT(prob, 8), gy = VECTOR_ELT(prob, 9);
    SEXP overall = VECTOR_ELT(prob, 10);
    P.roi = INTEGER(VECTOR_ELT(prob, 11));
    P.nroi = Rf_asInteger(VECTOR_ELT(prob, 12));
    P.nct = Rf_asInteger(VECTOR_ELT(prob, 13));
    P.ngroup = Rf_length(VECTOR_ELT(prob, 0)) - 1;
    P.ncell = Rf_length(VECTOR_ELT(prob, 3));
    P.ngrid = Rf_length(VECTOR_ELT(prob, 5));
    P.cs = lat[4];
    P.gx.assign(REAL(gx), REAL(gx) + Rf_length(gx));
    P.gy.assign(REAL(gy), REAL(gy) + Rf_length(gy));
    setupLattice(P, cell_a, cell_r, grid_a, grid_r, lat[0], lat[1], lat[2],
                 lat[3], lat[5]);

    P.roi_n.assign(P.nroi, 0);
    for (int i = 0; i < P.ngrid; i++) P.roi_n[P.roi[i]]++;
    for (int c = 0; c < P.nct - 1; c++) {
      for (int d = c + 1; d < P.nct; d++) {
        P.pc.push_back(c);
        P.pd.push_back(d);
      }
    }
    P.npair = (int)P.pc.size();
    P.nstat = P.npair * (1 + P.nroi);

    // per grid scaling of the smoothed counts to the overall density
    P.scal.assign(P.ngrid, 1.0);
    if (!Rf_isNull(overall)) {
      std::vector<double> tot(P.nslot, 0.0), H(P.nh), base(P.ngrid);
      for (int k = 0; k < P.ncell; k++) tot[P.pool0[k]] += 1.0;
      smooth(P, 1, tot.data(), H.data(), base.data());
      const double *ov = REAL(overall);
      for (int i = 0; i < P.ngrid; i++) {
        P.scal[i] = base[i] > 0.0 ? ov[i] / base[i] : 0.0;
      }
    }

    const int nb = Rf_asInteger(nb_);
    const int b_start = Rf_asInteger(b_start_);
    const uint64_t seed = (uint64_t)Rf_asReal(seed_);
    const int n_cpu = std::max(1, Rf_asInteger(cpu_threads));
    const int nstat = P.nstat;

    if (nb == 0) {
      Workspace ws(P);
      density(P, ws, NULL, false);
      statistics(P, ws);
      SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
      SEXP st = SET_VECTOR_ELT(out, 0, Rf_allocVector(REALSXP, nstat));
      std::copy(ws.stat.begin(), ws.stat.end(), REAL(st));
      SEXP dens = SET_VECTOR_ELT(out, 1, Rf_allocMatrix(REALSXP, P.ngrid, P.nct));
      std::copy(ws.D.begin(), ws.D.end(), REAL(dens));
      UNPROTECT(1);
      return out;
    }

    const double *obs = REAL(obs_);
    const double tol = 1e-12;
    const int nblk = (nb + BLOCK - 1) / BLOCK;
    std::vector<double> blk((size_t)nblk * 4 * nstat, 0.0);

    parallel(std::min(n_cpu, nblk), nblk, [&](int first, int last) {
      Workspace ws(P);
      pcg32_random_t rng;
      for (int k = first; k <= last; k++) {
        double *ge = blk.data() + (size_t)k * 4 * nstat;
        double *le = ge + nstat, *n = le + nstat, *sum = n + nstat;
        const int b_end = std::min(nb, (k + 1) * BLOCK);
        for (int b = k * BLOCK; b < b_end; b++) {
          uint64_t id = (uint64_t)b_start + (uint64_t)b;
          pcg32_srandom_r(&rng, splitmix64(seed * 0x9e3779b97f4a7c15ULL + id));
          density(P, ws, &rng, true);
          statistics(P, ws);
          for (int s = 0; s < nstat; s++) {
            const double v = ws.stat[s];
            if (ISNAN(v)) continue;
            n[s] += 1.0;
            sum[s] += v;
            if (!ISNAN(obs[s])) {
              ge[s] += (v >= obs[s] - tol);
              le[s] += (v <= obs[s] + tol);
            }
          }
        }
      }
    });

    SEXP acc = PROTECT(Rf_duplicate(acc_));
    double *pa = REAL(acc);
    for (int k = 0; k < nblk; k++) {
      const double *bk = blk.data() + (size_t)k * 4 * nstat;
      for (int j = 0; j < 4 * nstat; j++) pa[j] += bk[j];
    }
    UNPROTECT(1);
    return acc;
  }
}
