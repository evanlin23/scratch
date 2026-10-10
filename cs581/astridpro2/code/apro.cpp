// apro: internode-distance matrices for species-tree estimation from multi-copy gene-family trees.
//
//   apro -i GENES.nwk -o OUT.phy [-M multi|pro|pros] [-u] [-R 0|1] [-s FIRSTPASS.nwk] [-T]
//
//   -M multi : ASTRID-multi. All copy pairs, all internal nodes of the unrooted tree; per-gene mean
//              over copy pairs, then mean over genes (ASTRID's IndSpeciesMapping::average).
//   -M pro   : ASTRID-Pro. Orthologous pairs only (LCA tagged speciation), counting only speciation
//              nodes on the path (the LCA included). Per-gene mean, then mean over genes.
//   -M pros  : ASTRID-Pro-S. As pro, but a speciation node u entered from child c counts
//              1/s_hat(o), where o is the child of M(u) on the other side of the path in the rooted
//              first-pass species tree (-s) and s_hat is the survival estimated by reconciliation.
//   -u       : species = gene label up to the first '_' (SimPhy); default: the whole label.
//   -R 1     : count the gene-tree root (default 1 for pro/pros, as the theorem requires; ignored for
//              multi, which never counts the degree-2 root).
//   -T       : keep the given root and read tags from internal labels 'D' (simulator output).
//   -s FILE  : rooted first-pass species tree for pros.
// Gene trees are rooted by minimising the number of duplications (ties: fewest losses, at most
// 40 candidates), and a node is a duplication iff its children's species sets overlap.
// Pairs never observed are filled with the largest observed value; the count is reported on stderr.
#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;

struct Tree {
  vector<int> par;
  vector<vector<int>> ch;
  vector<string> lab;
  int root = 0;
  int add(int p) {
    par.push_back(p); ch.emplace_back(); lab.emplace_back();
    int i = (int)par.size() - 1;
    if (p >= 0) ch[p].push_back(i);
    return i;
  }
};

static Tree parse(const string& s) {
  Tree t;
  int cur = -1;
  size_t i = 0, n = s.size();
  if (i < n && s[i] == '[') { while (i < n && s[i] != ']') i++; i++; }
  while (i < n) {
    char c = s[i];
    if (c == '(') { cur = t.add(cur); i++; }
    else if (c == ',') { i++; }
    else if (c == ')') {
      i++;
      string l;
      while (i < n && !strchr(",():;[", s[i])) l += s[i++];
      t.lab[cur] = l;
      if (i < n && s[i] == '[') { while (i < n && s[i] != ']') i++; i++; }
      if (i < n && s[i] == ':') { i++; while (i < n && !strchr(",();[", s[i])) i++; }
      if (i < n && s[i] == '[') { while (i < n && s[i] != ']') i++; i++; }
      if (t.par[cur] >= 0) cur = t.par[cur];
    } else if (c == ';') break;
    else if (isspace((unsigned char)c)) i++;
    else {
      string l;
      if (c == '\'') { i++; while (i < n && s[i] != '\'') l += s[i++]; i++; }
      else while (i < n && !strchr(",():;[", s[i])) l += s[i++];
      int v = t.add(cur);
      t.lab[v] = l;
      if (i < n && s[i] == '[') { while (i < n && s[i] != ']') i++; i++; }
      if (i < n && s[i] == ':') { i++; while (i < n && !strchr(",();[", s[i])) i++; }
      if (i < n && s[i] == '[') { while (i < n && s[i] != ']') i++; i++; }
    }
  }
  t.root = 0;
  return t;
}

// ---------------------------------------------------------------- species bitsets
static int W = 1;  // words per bitset
struct BS {
  vector<uint64_t> w;
  void init() { w.assign(W, 0); }
};
static inline bool inter(const uint64_t* a, const uint64_t* b) {
  for (int i = 0; i < W; i++) if (a[i] & b[i]) return true;
  return false;
}
static inline int popc(const uint64_t* a) {
  int c = 0;
  for (int i = 0; i < W; i++) c += __builtin_popcountll(a[i]);
  return c;
}

static bool UNDERSCORE = false;
static unordered_map<string, int> SPIDX;
static vector<string> SPNAMES;
static int spof(const string& l) {
  string s = l;
  if (UNDERSCORE) { size_t p = s.find('_'); if (p != string::npos) s = s.substr(0, p); }
  auto it = SPIDX.find(s);
  if (it != SPIDX.end()) return it->second;
  int k = (int)SPNAMES.size();
  SPIDX[s] = k; SPNAMES.push_back(s);
  return k;
}

// rooted gene tree with tags
struct RTree {
  vector<int> par;
  vector<vector<int>> ch;
  vector<int> sp;    // leaf species, -1 for internal
  vector<char> dup;  // tag
  int root;
  vector<int> post;  // postorder
};

static void postorder(RTree& r) {
  r.post.clear();
  vector<pair<int, int>> st{{r.root, 0}};
  while (!st.empty()) {
    auto& [v, k] = st.back();
    if (k < (int)r.ch[v].size()) { int c = r.ch[v][k++]; st.push_back({c, 0}); }
    else { r.post.push_back(v); st.pop_back(); }
  }
}

// species sets + tags for a rooted tree; returns loss score (sum over dup nodes of |S(v)|-|S(c)|)
static long tag(RTree& r, vector<uint64_t>& S) {
  int n = (int)r.par.size();
  S.assign((size_t)n * W, 0);
  r.dup.assign(n, 0);
  long loss = 0;
  for (int v : r.post) {
    uint64_t* sv = &S[(size_t)v * W];
    if (r.ch[v].empty()) { sv[r.sp[v] >> 6] |= 1ULL << (r.sp[v] & 63); continue; }
    bool d = false;
    for (int c : r.ch[v]) {
      uint64_t* sc = &S[(size_t)c * W];
      if (!d && inter(sv, sc)) d = true;
      for (int i = 0; i < W; i++) sv[i] |= sc[i];
    }
    r.dup[v] = d;
  }
  for (int v : r.post)
    if (r.dup[v]) {
      int pv = popc(&S[(size_t)v * W]);
      for (int c : r.ch[v]) loss += pv - popc(&S[(size_t)c * W]);
    }
  return loss;
}

// build rooted tree on edge (a,b) of the unrooted adjacency nb
static RTree build(const vector<vector<int>>& nb, const vector<int>& leafsp, int a, int b) {
  RTree r;
  r.par.push_back(-1); r.ch.emplace_back(); r.sp.push_back(-1);
  r.root = 0;
  vector<tuple<int, int, int>> st{{a, b, 0}, {b, a, 0}};
  while (!st.empty()) {
    auto [v, frm, p] = st.back(); st.pop_back();
    int i = (int)r.par.size();
    r.par.push_back(p); r.ch.emplace_back(); r.sp.push_back(nb[v].size() == 1 ? leafsp[v] : -1);
    r.ch[p].push_back(i);
    for (int q : nb[v]) if (q != frm) st.push_back({q, v, i});
  }
  postorder(r);
  return r;
}

static RTree root_and_tag(const Tree& t, vector<uint64_t>& S, bool keep) {
  int n = (int)t.par.size();
  vector<int> leafsp(n, -1);
  for (int v = 0; v < n; v++) if (t.ch[v].empty()) leafsp[v] = spof(t.lab[v]);
  if (keep) {
    RTree r;
    r.par = t.par; r.ch = t.ch; r.sp = leafsp; r.root = t.root;
    postorder(r);
    tag(r, S);
    for (int v = 0; v < n; v++) if (!t.ch[v].empty()) r.dup[v] = (t.lab[v] == "D");
    return r;
  }
  // unrooted adjacency, suppress degree-2 root
  vector<vector<int>> nb(n);
  for (int v = 0; v < n; v++) if (t.par[v] >= 0) { nb[v].push_back(t.par[v]); nb[t.par[v]].push_back(v); }
  int rt = t.root;
  if (nb[rt].size() == 2) {
    int a = nb[rt][0], b = nb[rt][1];
    for (int* x : {&a, &b}) { auto& L = nb[*x]; L.erase(find(L.begin(), L.end(), rt)); }
    nb[a].push_back(b); nb[b].push_back(a); nb[rt].clear();
  }
  int r0 = -1;
  for (int v = 0; v < n; v++) if (nb[v].size() >= 2) { r0 = v; break; }
  if (r0 < 0) {  // 2 leaves
    int a = -1, b = -1;
    for (int v = 0; v < n; v++) if (nb[v].size() == 1) { if (a < 0) a = v; else b = v; }
    RTree r = build(nb, leafsp, a, b);
    tag(r, S);
    return r;
  }
  vector<int> par(n, -2), order;
  par[r0] = -1; order.push_back(r0);
  for (size_t k = 0; k < order.size(); k++) {
    int v = order[k];
    for (int q : nb[v]) if (par[q] == -2) { par[q] = v; order.push_back(q); }
  }
  vector<uint64_t> D((size_t)n * W, 0), U((size_t)n * W, 0);
  for (int k = (int)order.size() - 1; k >= 0; k--) {
    int v = order[k];
    uint64_t* dv = &D[(size_t)v * W];
    if (nb[v].size() == 1) { dv[leafsp[v] >> 6] |= 1ULL << (leafsp[v] & 63); continue; }
    for (int q : nb[v]) if (q != par[v]) { uint64_t* dq = &D[(size_t)q * W]; for (int i = 0; i < W; i++) dv[i] |= dq[i]; }
  }
  vector<uint64_t> tmp(W);
  for (int v : order) {
    vector<int> cs;
    for (int q : nb[v]) if (q != par[v]) cs.push_back(q);
    // prefix/suffix not needed for small degree: direct
    for (int c : cs) {
      uint64_t* uc = &U[(size_t)c * W];
      const uint64_t* uv = &U[(size_t)v * W];
      for (int i = 0; i < W; i++) uc[i] = uv[i];
      for (int c2 : cs) if (c2 != c) { uint64_t* d2 = &D[(size_t)c2 * W]; for (int i = 0; i < W; i++) uc[i] |= d2[i]; }
    }
  }
  auto side = [&](int v, int q) -> const uint64_t* {  // species on q's side of edge v-q
    return par[q] == v ? &D[(size_t)q * W] : &U[(size_t)v * W];
  };
  auto dupgiven = [&](int v, int p) -> int {  // is v a duplication when its parent is p
    if (nb[v].size() == 1) return 0;
    for (int i = 0; i < W; i++) tmp[i] = 0;
    for (int q : nb[v]) if (q != p) {
      const uint64_t* s = side(v, q);
      if (inter(tmp.data(), s)) return 1;
      for (int i = 0; i < W; i++) tmp[i] |= s[i];
    }
    return 0;
  };
  long base = 0;
  for (int w : order) if (w != r0 && nb[w].size() > 1) base += dupgiven(w, par[w]);
  vector<long> g(n, 0);
  for (int k = 1; k < (int)order.size(); k++) {
    int v = order[k], p = par[v];
    if (p == r0) g[v] = base + dupgiven(r0, v);
    else g[v] = g[p] - dupgiven(p, par[p]) + dupgiven(p, v);
  }
  // rooting on edge (v, par v): score = dups of all nodes with that orientation + new root
  long best = LONG_MAX;
  vector<int> cands;
  for (int k = 1; k < (int)order.size(); k++) {
    int v = order[k];
    // g[v]: duplications of all nodes when rooted on edge (v, par v), plus the new root node
    long sc = g[v];
    sc += inter(&D[(size_t)v * W], &U[(size_t)v * W]) ? 1 : 0;
    if (sc < best) { best = sc; cands.clear(); }
    if (sc == best) cands.push_back(v);
  }
  RTree bestr;
  long bl = LONG_MAX;
  int lim = min((int)cands.size(), 40);
  for (int j = 0; j < lim; j++) {
    int v = cands[j];
    RTree r = build(nb, leafsp, v, par[v]);
    vector<uint64_t> S2;
    long l = tag(r, S2);
    if (lim == 1) { S.swap(S2); return r; }
    if (l < bl) { bl = l; bestr = std::move(r); S.swap(S2); }
  }
  return bestr;
}

// Root the first-pass species tree on the edge that minimises the total number of duplications
// (LCA reconciliation) over at most 200 gene trees (evenly spaced). Input may be rooted or not.
static Tree reroot_on(const vector<vector<int>>& nb, const vector<string>& lab, int a, int b) {
  Tree t;
  int r = t.add(-1);
  vector<tuple<int, int, int>> stck{{a, b, r}, {b, a, r}};
  while (!stck.empty()) {
    auto [v, frm, p] = stck.back(); stck.pop_back();
    int i = t.add(p);
    if (nb[v].size() == 1) t.lab[i] = lab[v];
    for (int q : nb[v]) if (q != frm) stck.push_back({q, v, i});
  }
  return t;
}
static Tree root_species(const Tree& st, const vector<RTree>& rts, vector<int>& sleaf_unused) {
  int n = (int)st.par.size();
  vector<vector<int>> nb(n);
  for (int v = 0; v < n; v++) if (st.par[v] >= 0) { nb[v].push_back(st.par[v]); nb[st.par[v]].push_back(v); }
  int rt = st.root;
  if (nb[rt].size() == 2) {
    int a = nb[rt][0], b = nb[rt][1];
    for (int* x : {&a, &b}) { auto& L = nb[*x]; L.erase(find(L.begin(), L.end(), rt)); }
    nb[a].push_back(b); nb[b].push_back(a); nb[rt].clear();
  }
  vector<size_t> use;
  size_t G = rts.size(), m = min<size_t>(G, 200);
  for (size_t j = 0; j < m; j++) use.push_back(j * G / m);
  long best = LONG_MAX;
  Tree bt;
  for (int a = 0; a < n; a++)
    for (int b : nb[a]) {
      if (b < a) continue;
      Tree t = reroot_on(nb, st.lab, a, b);
      int tn = (int)t.par.size();
      vector<int> dep(tn, 0), leafof(SPNAMES.size(), -1);
      for (int v = 1; v < tn; v++) dep[v] = dep[t.par[v]] + 1;  // parents precede children
      for (int v = 0; v < tn; v++) if (t.ch[v].empty()) leafof[spof(t.lab[v])] = v;
      auto lca = [&](int x, int y) {
        if (x < 0) return y;
        while (dep[x] > dep[y]) x = t.par[x];
        while (dep[y] > dep[x]) y = t.par[y];
        while (x != y) { x = t.par[x]; y = t.par[y]; }
        return x;
      };
      long dups = 0;
      for (size_t g : use) {
        const RTree& r = rts[g];
        vector<int> M(r.par.size(), -1);
        for (int v : r.post) {
          if (r.ch[v].empty()) { M[v] = leafof[r.sp[v]]; continue; }
          int mm = -1;
          for (int c : r.ch[v]) mm = lca(mm, M[c]);
          M[v] = mm;
          for (int c : r.ch[v]) if (M[c] == mm) { dups++; break; }
        }
      }
      if (dups < best) { best = dups; bt = t; }
    }
  fprintf(stderr, "first-pass root: %ld duplications over %zu genes\n", best, use.size());
  return bt;
}

// ---------------------------------------------------------------- distances
struct Ent { int sp; double n; double s; };  // species, #leaves, sum of path weights to them

int main(int argc, char** argv) {
  string in, out, mode = "pro", first;
  int countroot = 1;
  bool keep = false;
  for (int i = 1; i < argc; i++) {
    string a = argv[i];
    if (a == "-i") in = argv[++i];
    else if (a == "-o") out = argv[++i];
    else if (a == "-M") mode = argv[++i];
    else if (a == "-u") UNDERSCORE = true;
    else if (a == "-R") countroot = atoi(argv[++i]);
    else if (a == "-T") keep = true;
    else if (a == "-s") first = argv[++i];
  }
  auto t0 = chrono::steady_clock::now();
  vector<string> lines;
  { ifstream f(in); string l; while (getline(f, l)) if (l.find(';') != string::npos) lines.push_back(l); }
  vector<Tree> genes;
  genes.reserve(lines.size());
  for (auto& l : lines) genes.push_back(parse(l));
  for (auto& t : genes) for (size_t v = 0; v < t.par.size(); v++) if (t.ch[v].empty()) spof(t.lab[v]);
  // first-pass species tree for pros: register its species first is not needed; map by name
  Tree st;
  if (mode == "pros") { ifstream f(first); string l; getline(f, l); st = parse(l); for (size_t v = 0; v < st.par.size(); v++) if (st.ch[v].empty()) spof(st.lab[v]); }
  int k = (int)SPNAMES.size();
  W = (k + 63) / 64;
  // species tree structures for pros
  vector<int> sdepth, sleaf(k, -1);
  if (mode == "pros") {
    sdepth.assign(st.par.size(), 0);
    vector<int> ord{st.root};
    for (size_t j = 0; j < ord.size(); j++) for (int c : st.ch[ord[j]]) { sdepth[c] = sdepth[ord[j]] + 1; ord.push_back(c); }
    for (size_t v = 0; v < st.par.size(); v++) if (st.ch[v].empty()) sleaf[spof(st.lab[v])] = (int)v;
  }
  auto slca = [&](int a, int b) {
    if (a < 0) return b;
    if (b < 0) return a;
    while (sdepth[a] > sdepth[b]) a = st.par[a];
    while (sdepth[b] > sdepth[a]) b = st.par[b];
    while (a != b) { a = st.par[a]; b = st.par[b]; }
    return a;
  };
  auto t1 = chrono::steady_clock::now();
  vector<RTree> rts(genes.size());
  {
    vector<uint64_t> S;
    for (size_t g = 0; g < genes.size(); g++) rts[g] = root_and_tag(genes[g], S, keep);
  }
  auto t2 = chrono::steady_clock::now();
  bool pro = mode != "multi";
  if (mode == "pros") st = root_species(st, rts, sleaf);
  if (mode == "pros") {  // recompute depths and leaf ids for the rerooted tree
    sdepth.assign(st.par.size(), 0);
    vector<int> ord{st.root};
    for (size_t j = 0; j < ord.size(); j++) for (int c : st.ch[ord[j]]) { sdepth[c] = sdepth[ord[j]] + 1; ord.push_back(c); }
    for (size_t v = 0; v < st.par.size(); v++) if (st.ch[v].empty()) sleaf[spof(st.lab[v])] = (int)v;
  }
  // pros: reconciliation and survival estimates
  vector<double> shat;
  vector<vector<int>> maps(genes.size());
  if (mode == "pros") {
    int ns = (int)st.par.size();
    vector<double> yes(ns, 0), tot(ns, 0);
    for (size_t g = 0; g < rts.size(); g++) {
      RTree& r = rts[g];
      vector<int>& M = maps[g];
      M.assign(r.par.size(), -1);
      for (int v : r.post) {
        if (r.ch[v].empty()) { M[v] = sleaf[r.sp[v]]; continue; }
        int m = -1;
        for (int c : r.ch[v]) m = slca(m, M[c]);
        M[v] = m;
      }
      // x maximal inside clade w=M-side: the species child w of M(p) containing M(x)
      for (int x : r.post) {
        int p = r.par[x];
        if (p < 0) continue;
        if (M[x] == M[p]) continue;  // x not maximal in a strict sub-clade of M(p)
        // w = child of species node on the path from M(x) up to (not including) M(p): the
        // lineage passes through every species node strictly between; at each such node u
        // (from M(x) upward), the sibling of the child on the path was lost (speciation not seen),
        // except at M(p) if p is a speciation.
        int w = M[x];
        while (st.par[w] != M[p]) {
          int u = st.par[w];
          for (int o : st.ch[u]) if (o != w) tot[o] += 1;  // sibling o: lost
          w = u;
        }
        for (int o : st.ch[M[p]]) if (o != w) { tot[o] += 1; if (!r.dup[p]) yes[o] += 1; }
      }
    }
    shat.assign(ns, 1.0);
    for (int v = 0; v < ns; v++) if (tot[v] > 0) shat[v] = max(yes[v] / tot[v], 0.02);
    // write estimates to stderr
    for (int v = 0; v < ns; v++) if (st.ch[v].empty()) cerr << "shat " << st.lab[v] << " " << shat[v] << " n=" << tot[v] << "\n";
  }
  vector<double> SUM((size_t)k * k, 0), NG((size_t)k * k, 0), tot((size_t)k * k, 0), cnt((size_t)k * k, 0);
  vector<int> touched;
  vector<vector<Ent>> E;
  vector<double> cum, wt;
  for (size_t g = 0; g < rts.size(); g++) {
    RTree& r = rts[g];
    int n = (int)r.par.size();
    // counted weight of entering node u from child c: stored on c as wt[c]
    wt.assign(n, 0); cum.assign(n, 0);
    vector<char> counted(n, 0);
    for (int v = 0; v < n; v++) counted[v] = !r.ch[v].empty() && (!pro || !r.dup[v]);
    if (!pro || !countroot) if (r.ch[r.root].size() < 3) counted[r.root] = 0;
    for (int c = 0; c < n; c++) {
      int u = r.par[c];
      if (u < 0 || !counted[u]) continue;
      if (mode == "pros") {
        // other side: species child of M(u) not containing M(c)
        int w = maps[g][c];
        if (w == maps[g][u]) { wt[c] = 1.0; continue; }  // inferred S node whose child spans M(u)
        while (st.par[w] != maps[g][u]) w = st.par[w];
        int o = -1;
        for (int q : st.ch[maps[g][u]]) if (q != w) o = q;
        wt[c] = 1.0 / shat[o];
      } else wt[c] = 1.0;
    }
    // cum[x] = sum of wt over edges from root down to x (weights of nodes entered from below)
    for (int j = n - 1; j >= 0; j--) { int v = r.post[j]; int p = r.par[v]; cum[v] = p < 0 ? 0 : cum[p] + wt[v]; }
    // leaf a, LCA v, child c of v on a's side: path weight from a up to v excluding v = cum[a]-cum[c]
    E.assign(n, {});
    for (int v : r.post) {
      if (r.ch[v].empty()) { E[v].push_back({r.sp[v], 1.0, cum[v]}); continue; }
      bool use = !pro || !r.dup[v];
      double self = counted[v] ? 1.0 : 0.0;  // the LCA counts 1
      auto& cs = r.ch[v];
      if (use)
        for (size_t i = 0; i < cs.size(); i++)
          for (size_t j = i + 1; j < cs.size(); j++)
            for (auto& A : E[cs[i]])
              for (auto& B : E[cs[j]]) {
                if (A.sp == B.sp) continue;
                // sum over pairs of (cum[a]-cum[ci]) + (cum[b]-cum[cj]) + self
                double s = A.s * B.n + B.s * A.n - A.n * B.n * (cum[cs[i]] + cum[cs[j]] - self);
                double c2 = A.n * B.n;
                int x = min(A.sp, B.sp), y = max(A.sp, B.sp);
                size_t id = (size_t)x * k + y;
                if (cnt[id] == 0) touched.push_back((int)id);
                tot[id] += s; cnt[id] += c2;
              }
      // merge children's entries (sorted by species)
      vector<Ent> m;
      for (int c : cs) {
        vector<Ent> nm;
        nm.reserve(m.size() + E[c].size());
        size_t a = 0, b = 0;
        auto& X = E[c];
        while (a < m.size() || b < X.size()) {
          if (b == X.size() || (a < m.size() && m[a].sp < X[b].sp)) nm.push_back(m[a++]);
          else if (a == m.size() || X[b].sp < m[a].sp) nm.push_back(X[b++]);
          else { nm.push_back({m[a].sp, m[a].n + X[b].n, m[a].s + X[b].s}); a++; b++; }
        }
        m.swap(nm);
        vector<Ent>().swap(E[c]);
      }
      E[v].swap(m);
    }
    for (int id : touched) { SUM[id] += tot[id] / cnt[id]; NG[id] += 1; tot[id] = 0; cnt[id] = 0; }
    touched.clear();
  }
  auto t3 = chrono::steady_clock::now();
  double mx = 0;
  int miss = 0;
  vector<double> Dm((size_t)k * k, 0);
  for (int x = 0; x < k; x++) for (int y = x + 1; y < k; y++) {
    size_t id = (size_t)x * k + y;
    if (NG[id] > 0) { Dm[id] = SUM[id] / NG[id]; mx = max(mx, Dm[id]); } else miss++;
  }
  FILE* f = fopen(out.c_str(), "w");
  fprintf(f, "%d\n", k);
  for (int x = 0; x < k; x++) {
    fprintf(f, "%s", SPNAMES[x].c_str());
    for (int y = 0; y < k; y++) {
      double d = 0;
      if (x != y) { size_t id = (size_t)min(x, y) * k + max(x, y); d = NG[id] > 0 ? Dm[id] : mx; }
      fprintf(f, " %.6f", d);
    }
    fprintf(f, "\n");
  }
  fclose(f);
  auto sec = [](auto a, auto b) { return chrono::duration<double>(b - a).count(); };
  fprintf(stderr, "genes %zu species %d missing_pairs %d parse %.3f root_tag %.3f dist %.3f\n", genes.size(), k, miss,
          sec(t0, t1), sec(t1, t2), sec(t2, t3));
  return 0;
}
