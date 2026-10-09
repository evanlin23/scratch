// qtool: explicit quartet table from gene trees, quartet scoring and SPR hill-climbing.
//
//   qtool score  -g genes.tre -t trees.tre [-w gene_weights] [-m mode]
//   qtool search -g genes.tre -t start.tre [-w gene_weights] [-m mode] [-s seed] [-T seconds]
//   qtool genesupport -g genes.tre -t species.tre    (per-gene fraction of quartets agreeing)
//
// Gene trees are unrooted; missing taxa are allowed. Quartet q = {i<j<k<l}; topology 0 = ij|kl,
// 1 = ik|jl, 2 = il|jk. Q[q][t] = sum of gene weights of gene trees that induce topology t.
// The score of a species tree T is sum_q Q[q][top_T(q)] (the ASTRAL / MQSST objective).
// Modes transform Q before scoring/search:
//   raw       as counted
//   capminor  larger of the two minority topologies capped at the smaller one (MSC predicts
//             equal minorities; directional excess is attributed to HGT/error)
//   vote      1 for the dominant topology, 0 otherwise (unweighted dominant-quartet consensus)
#include <bits/stdc++.h>
using namespace std;

static map<string, int> taxon_id;
static vector<string> taxon_name;
static int tid(const string& s) {
  auto it = taxon_id.find(s);
  if (it != taxon_id.end()) return it->second;
  int k = taxon_name.size();
  taxon_id[s] = k;
  taxon_name.push_back(s);
  return k;
}

// Unrooted tree: nodes with adjacency lists; leaves carry a taxon id.
struct Tree {
  vector<vector<int>> adj;
  vector<int> taxon;  // -1 for internal
  int add(int t) { adj.emplace_back(); taxon.push_back(t); return adj.size() - 1; }
  void link(int a, int b) { adj[a].push_back(b); adj[b].push_back(a); }
  void unlink(int a, int b) {
    adj[a].erase(find(adj[a].begin(), adj[a].end(), b));
    adj[b].erase(find(adj[b].begin(), adj[b].end(), a));
  }
};

// Minimal Newick parser: ignores branch lengths, internal labels and [comments].
static Tree parse_newick(const string& s) {
  Tree T;
  size_t p = 0;
  function<int()> rec = [&]() -> int {
    while (isspace(s[p])) p++;
    int node;
    if (s[p] == '(') {
      p++;
      node = T.add(-1);
      while (true) {
        int c = rec();
        T.link(node, c);
        while (isspace(s[p])) p++;
        if (s[p] == ',') { p++; continue; }
        if (s[p] == ')') { p++; break; }
        throw runtime_error("bad newick near " + to_string(p));
      }
      // internal label
      while (p < s.size() && s[p] != ':' && s[p] != ',' && s[p] != ')' && s[p] != ';') p++;
    } else {
      string name;
      if (s[p] == '\'') {
        p++;
        while (s[p] != '\'') name += s[p++];
        p++;
      } else
        while (p < s.size() && s[p] != ':' && s[p] != ',' && s[p] != ')' && s[p] != ';' &&
               !isspace(s[p]))
          name += s[p++];
      node = T.add(tid(name));
    }
    if (p < s.size() && s[p] == '[') { while (s[p] != ']') p++; p++; }
    if (p < s.size() && s[p] == ':') {
      p++;
      while (p < s.size() && s[p] != ',' && s[p] != ')' && s[p] != ';') p++;
    }
    return node;
  };
  int root = rec();
  // suppress degree-2 root (rooted input)
  if (T.adj[root].size() == 2) {
    int a = T.adj[root][0], b = T.adj[root][1];
    T.unlink(root, a); T.unlink(root, b); T.link(a, b);
  }
  return T;
}

static vector<string> read_trees(const string& fn) {
  ifstream in(fn);
  if (!in) throw runtime_error("cannot open " + fn);
  vector<string> out;
  string cur, line;
  while (getline(in, line)) {
    for (char c : line) {
      cur += c;
      if (c == ';') { out.push_back(cur); cur.clear(); }
    }
  }
  return out;
}

// Leaf-to-leaf edge distances (n x n, 255 = taxon absent).
static int N;
static vector<uint8_t> leaf_dist(const Tree& T) {
  vector<uint8_t> D(N * N, 255);
  int m = T.adj.size();
  vector<int> dist(m), st;
  for (int s = 0; s < m; s++) {
    if (T.taxon[s] < 0 || T.adj[s].empty()) continue;
    fill(dist.begin(), dist.end(), -1);
    dist[s] = 0;
    st.assign(1, s);
    while (!st.empty()) {
      int u = st.back(); st.pop_back();
      for (int v : T.adj[u])
        if (dist[v] < 0) { dist[v] = dist[u] + 1; st.push_back(v); }
    }
    int a = T.taxon[s];
    for (int t = 0; t < m; t++)
      if (T.taxon[t] >= 0 && dist[t] >= 0) D[a * N + T.taxon[t]] = min(dist[t], 254);
  }
  return D;
}

static vector<array<uint8_t, 4>> quartets;
static void make_quartets() {
  quartets.clear();
  for (int i = 0; i < N; i++)
    for (int j = i + 1; j < N; j++)
      for (int k = j + 1; k < N; k++)
        for (int l = k + 1; l < N; l++) quartets.push_back({(uint8_t)i, (uint8_t)j, (uint8_t)k, (uint8_t)l});
}

// topology of quartet in tree with distance matrix D; -1 if unresolved or taxon missing
static inline int topo(const uint8_t* D, const array<uint8_t, 4>& q) {
  int i = q[0], j = q[1], k = q[2], l = q[3];
  if (D[i * N + i] == 255 || D[j * N + j] == 255 || D[k * N + k] == 255 || D[l * N + l] == 255) return -1;
  int s0 = D[i * N + j] + D[k * N + l], s1 = D[i * N + k] + D[j * N + l], s2 = D[i * N + l] + D[j * N + k];
  if (s0 < s1 && s0 < s2) return 0;
  if (s1 < s0 && s1 < s2) return 1;
  if (s2 < s0 && s2 < s1) return 2;
  return -1;
}


static double score(const Tree& T, const vector<array<float, 3>>& Q) {
  auto D = leaf_dist(T);
  const uint8_t* d = D.data();
  double s = 0;
  for (size_t q = 0; q < quartets.size(); q++) {
    int t = topo(d, quartets[q]);
    if (t >= 0) s += Q[q][t];
  }
  return s;
}

static string to_newick(const Tree& T) {
  int start = -1;
  for (int i = 0; i < (int)T.adj.size(); i++)
    if (T.taxon[i] >= 0 && !T.adj[i].empty()) { start = i; break; }
  function<string(int, int)> rec = [&](int u, int par) -> string {
    if (T.taxon[u] >= 0 && u != start) return taxon_name[T.taxon[u]];
    string s = "(";
    bool first = true;
    for (int v : T.adj[u])
      if (v != par) { if (!first) s += ","; s += rec(v, u); first = false; }
    return s + ")";
  };
  int nb = T.adj[start][0];
  return "(" + taxon_name[T.taxon[start]] + "," + rec(nb, start).substr(1);
}

int main(int argc, char** argv) {
  if (argc < 2) { cerr << "usage: qtool score|search|genesupport ...\n"; return 1; }
  string cmd = argv[1], gfile, tfile, wfile, mode = "raw";
  unsigned seed = 1;
  double tlimit = 1e9;
  for (int a = 2; a < argc; a++) {
    string o = argv[a];
    if (o == "-g") gfile = argv[++a];
    else if (o == "-t") tfile = argv[++a];
    else if (o == "-w") wfile = argv[++a];
    else if (o == "-m") mode = argv[++a];
    else if (o == "-s") seed = atoi(argv[++a]);
    else if (o == "-T") tlimit = atof(argv[++a]);
  }
  auto gstr = read_trees(gfile);
  auto tstr = read_trees(tfile);
  vector<Tree> genes;
  for (auto& s : gstr) genes.push_back(parse_newick(s));
  vector<Tree> sp;
  for (auto& s : tstr) sp.push_back(parse_newick(s));
  N = taxon_name.size();
  make_quartets();
  vector<double> w(genes.size(), 1.0);
  if (!wfile.empty()) {
    ifstream in(wfile);
    for (auto& x : w) in >> x;
  }
  auto t0 = chrono::steady_clock::now();
  vector<array<float, 3>> Q(quartets.size(), {0, 0, 0});
  vector<vector<uint8_t>> GD;
  for (size_t g = 0; g < genes.size(); g++) {
    auto D = leaf_dist(genes[g]);
    for (int i = 0; i < N; i++) D[i * N + i] = (D[i * N + i] == 255 ? 255 : 0);
    const uint8_t* d = D.data();
    for (size_t q = 0; q < quartets.size(); q++) {
      int t = topo(d, quartets[q]);
      if (t >= 0) Q[q][t] += w[g];
    }
    if (cmd == "genesupport") GD.push_back(move(D));
  }
  double total = 0;
  for (auto& x : Q) total += x[0] + x[1] + x[2];
  if (mode == "capminor" || mode == "vote") {
    for (auto& x : Q) {
      int mx = max_element(x.begin(), x.end()) - x.begin();
      int a = (mx + 1) % 3, b = (mx + 2) % 3;
      if (mode == "capminor") { float m = min(x[a], x[b]); x[a] = x[b] = m; }
      else { bool tie = x[mx] == x[a] || x[mx] == x[b] || x[mx] == 0; x = {0, 0, 0}; if (!tie) x[mx] = 1; }
    }
  }
  double tq = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
  cerr << "taxa " << N << " genes " << genes.size() << " quartets " << quartets.size()
       << " table_sec " << tq << "\n";

  if (cmd == "score") {
    // raw score and normalized (by total gene-tree quartets, as ASTRAL reports)
    for (auto& T : sp) {
      double s = score(T, Q);
      printf("%.6f\t%.8f\n", s, s / total);
    }
    return 0;
  }
  if (cmd == "genesupport") {
    auto D = leaf_dist(sp[0]);
    for (int i = 0; i < N; i++) D[i * N + i] = 0;
    for (auto& G : GD) {
      long agree = 0, tot = 0;
      for (auto& q : quartets) {
        int t = topo(G.data(), q);
        if (t < 0) continue;
        tot++;
        if (topo(D.data(), q) == t) agree++;
      }
      printf("%.6f\n", tot ? (double)agree / tot : 1.0);
    }
    return 0;
  }
  if (cmd == "search") {
    mt19937 rng(seed);
    Tree T = sp[0];
    double best = score(T, Q), start = best;
    long evals = 1, moves = 0;
    auto tstart = chrono::steady_clock::now();
    bool improved = true;
    while (improved) {
      improved = false;
      // candidate list: (u, v) prune v-side subtree hanging from internal u; regraft on edge (x, y)
      vector<pair<int, int>> dir;
      for (int u = 0; u < (int)T.adj.size(); u++)
        if (T.adj[u].size() == 3)
          for (int v : T.adj[u]) dir.push_back({u, v});
      shuffle(dir.begin(), dir.end(), rng);
      for (auto [u, v] : dir) {
        if (chrono::duration<double>(chrono::steady_clock::now() - tstart).count() > tlimit) break;
        // detach
        vector<int> others;
        for (int x : T.adj[u]) if (x != v) others.push_back(x);
        int a = others[0], b = others[1];
        // collect subtree nodes on v side
        vector<char> insub(T.adj.size(), 0);
        vector<int> st = {v};
        insub[v] = 1;
        while (!st.empty()) {
          int x = st.back(); st.pop_back();
          for (int y : T.adj[x]) if (y != u && !insub[y]) { insub[y] = 1; st.push_back(y); }
        }
        T.unlink(u, a); T.unlink(u, b); T.link(a, b);
        vector<pair<int, int>> edges;
        for (int x = 0; x < (int)T.adj.size(); x++)
          if (x != u && !insub[x])
            for (int y : T.adj[x])
              if (y > x && !insub[y] && y != u && !(min(x, y) == min(a, b) && max(x, y) == max(a, b)))
                edges.push_back({x, y});
        shuffle(edges.begin(), edges.end(), rng);
        bool done = false;
        for (auto [x, y] : edges) {
          T.unlink(x, y); T.link(x, u); T.link(u, y);
          double s = score(T, Q);
          evals++;
          if (s > best + 1e-6) { best = s; moves++; done = true; break; }
          T.unlink(x, u); T.unlink(u, y); T.link(x, y);
        }
        if (done) { improved = true; break; }
        T.unlink(a, b); T.link(u, a); T.link(u, b);
      }
    }
    double sec = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
    printf("%s;\n", to_newick(T).c_str());
    fprintf(stderr, "start %.3f final %.3f moves %ld evals %ld sec %.2f\n", start, best, moves, evals, sec);
    return 0;
  }
  return 1;
}
