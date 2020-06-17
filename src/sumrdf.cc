#include <cassert>
#include "../include/sumrdf.h"
#include "../include/util.h"

#include <random>
#include <set>

SumRDF::Type::Type() {
  classes_.clear();
  outgoing_.clear();
  incoming_.clear();
}

// Convert type to string
string SumRDF::Type::ToString() {
  string s = "";
  // C == classes_
  for (int x : classes_) {
    s = std::to_string(x) + s;
    s = "|" + s;
  }
  s = "|" + s;
  // O == outgoing_
  for (int x : outgoing_) {
    s = std::to_string(x) + s;
    s = "|" + s;
  }
  s = "|" + s;
  // I == incoming_
  for (int x : incoming_) {
    s = std::to_string(x) + s;
    s = "|" + s;
  }
  s = "|" + s;
  return s;
}

// Convert type to string using only classes (C(d) == classes_)
string SumRDF::Type::ToStringByClass() {
  string s = "";
  // C == classes_
  for (int x : classes_) {
    s = std::to_string(x) + s;
    s = "|" + s;
  }
  s = "|" + s;
  return s;
}

SumRDF::Resource::Resource() {
  type_ = Type();
}

SumRDF::Bucket::Bucket(Type& t) {
  type_.classes_ = t.classes_;
  type_.outgoing_ = t.outgoing_;
  type_.incoming_ = t.incoming_;
  resources_.clear();
}

SumRDF::Partition::Partition() { }
SumRDF::Partition::Partition(vector<int>& p, int block_cnt, vector<vector<int>>& block_list) {
  partition_ = p;
  blocks_.clear();
  for (int i = 0; i < block_cnt; i++) {
    blocks_.push_back(block_list[i]);
  }
}

void SumRDF::PrepareSummaryStructure(DataGraph& g, double ratio) {
  // 5.1 Build a typed summary
  // Compute types for each data vertex
  g_resources_.resize(0); g_resources_.resize(g.GetNumVertices());
  for (size_t i = 0; i < g.GetNumVertices(); i++) {
    for (auto rv = g.GetVLabels(i); rv.begin != rv.end; rv.begin++) {
      int vl = *rv.begin;
      // C == type_.classes_
      // For fair comparisons, we assign vertex labels to C for URIs and literals
      g_resources_[i].type_.classes_.insert(vl);
    }
    for (auto re = g.GetELabels(i, true); re.begin != re.end; re.begin++) {
      int el = *re.begin;
      // O == type_.outgoing_
      g_resources_[i].type_.outgoing_.insert(el);
    }
    for (auto re = g.GetELabels(i, false); re.begin != re.end; re.begin++) {
      int el = *re.begin;
      // I == type_.incoming_
      g_resources_[i].type_.incoming_.insert(el);
    }
  }
  double target_size = ratio * DataGraphSize(g);
  target_ = (target_size - sizeof (int) * (static_cast<size_t>(g.GetNumVertices()))) / (sizeof (int) * 3);
  target_ = std::max(target_, 10000.0);
  unordered_map<string, int> bucket_idx; // bucket_idx maps type (string) to a summary vertex
  vector<int> mu; // mu: maps a data vertex v to a summary vertex mu[v]
  w1_.resize(0);
  // Compute mu(v)
  // w1_ denotes w(b) in the paper
  for (size_t i = 0; i < g_resources_.size(); i++) {
    Resource& r = g_resources_[i];
    string s = r.type_.ToString();
    if (!bucket_idx.count(s)) { // if new type exists, create new b_{t_v}
      sm_.CreateBucket(r, i);
      bucket_idx[s] = sm_.buckets_.size() - 1;
      w1_.push_back(0);
    }
    mu.push_back(bucket_idx[s]);
    sm_.buckets_[bucket_idx[s]].resources_.push_back(i);
    w1_[bucket_idx[s]]++; // increase weight of the corresponding summary vertex
  }
  // Add edges to summary graph
  unordered_map<int, unordered_map<int, unordered_map<int, int>>> edge_idx; // edge_idx maps an edge (srcid, el, dstid) to an edge number
  sm_.adj_list_.clear();
  sm_.adj_list_.resize(sm_.buckets_.size());
  sm_.rev_adj_list_.clear();
  sm_.rev_adj_list_.resize(sm_.buckets_.size());
  for (int srcid = 0; srcid < g.GetNumVertices(); srcid++) {
    for (auto re = g.GetELabels(srcid, true); re.begin != re.end; re.begin++) {
      int el = *re.begin;
      for (auto rd = g.GetAdj(srcid, el, true); rd.begin != rd.end; rd.begin++) {
        int dstid = *rd.begin;
        AddEdgeToSummary(mu[srcid], el, mu[dstid], edge_idx);
      }
    }
  }
  // Sort edges and adjust the order of weights
  vector<int> cp(w2_.size());
  cp.swap(w2_);
  std::sort(sm_.edges_.begin(), sm_.edges_.end());
  for (size_t i = 0; i < sm_.edges_.size(); i++) {
    int src = sm_.edges_[i].src;
    int dst = sm_.edges_[i].dst;
    int el = sm_.edges_[i].el;
    w2_[i] = cp[edge_idx[src][el][dst]];
  }
  // 5.2 Summary Refinement by MinHashing
  // if a typed summary is large, we compress it by merging similar buckets and types
  if (sm_.multiplicity_ > target_) {
    // Set trivial summary
    mu_.resize(sm_.buckets_.size());
    for (size_t i = 0; i < sm_.buckets_.size(); i++) {
      mu_[i] = i;
    }
    // Summary refinement by MinHashing
    edge_idx.clear();
    SummarizeUsingMinHash(g, ratio);
    // After building a summary using MinHashing, set summary edges and weights
    sm_.multiplicity_ = 0;
    // sm_.edges_[i]: the i-th edge of the summary graph S
    sm_.edges_.clear();
    sm_.adj_list_.clear();
    sm_.adj_list_.resize(sm_.buckets_.size());
    sm_.rev_adj_list_.clear();
    sm_.rev_adj_list_.resize(sm_.buckets_.size());
    w2_.resize(0);
    w1_.clear(); w1_.resize(sm_.buckets_.size(), 0);
    unordered_map<int, int> renum;
    int cnt = 0;
    sm_.buckets_.clear();
    for (int srcid = 0; srcid < g.GetNumVertices(); srcid++) {
      for (auto re = g.GetELabels(srcid, true); re.begin != re.end; re.begin++) {
        int el = *re.begin;
        for (auto rd = g.GetAdj(srcid, el, true); rd.begin != rd.end; rd.begin++) {
          int dstid = *rd.begin;
          int b1 = Find(mu[srcid]);
          int b2 = Find(mu[dstid]);
          if (renum.count(b1) == 0) renum[b1] = cnt++;
          if (renum.count(b2) == 0) renum[b2] = cnt++;
          b1 = renum[b1]; b2 = renum[b2];
          if (sm_.buckets_.size() <= b1) sm_.buckets_.resize(b1 + 1, Bucket());
          if (sm_.buckets_.size() <= b2) sm_.buckets_.resize(b2 + 1, Bucket());
          AddEdgeToSummary(b1, el, b2, edge_idx);
          sm_.buckets_[b1].type_.outgoing_.insert(el);
          sm_.buckets_[b2].type_.incoming_.insert(el);
        }
      }
      int b = Find(mu[srcid]);
      if (renum.count(b) == 0) renum[b] = cnt++;
      b = renum[b];
      if (sm_.buckets_.size() <= b) sm_.buckets_.resize(b + 1, Bucket());
      w1_[b]++;
      sm_.buckets_[b].resources_.push_back(srcid);
      for (auto r = g.GetVLabels(srcid); r.begin != r.end; r.begin++) {
        int vl = *r.begin;
        sm_.buckets_[b].type_.classes_.insert(vl);
      }
    }
    w1_.erase(w1_.begin() + cnt, w1_.end());
    sm_.adj_list_.erase(sm_.adj_list_.begin() + cnt, sm_.adj_list_.end());
    sm_.rev_adj_list_.erase(sm_.rev_adj_list_.begin() + cnt, sm_.rev_adj_list_.end());
    // Sort edges and adjust the order of weights
    vector<int> cp(w2_.size());
    std::copy(w2_.begin(), w2_.end(), cp.begin());
    std::sort(sm_.edges_.begin(), sm_.edges_.end());
    for (size_t i = 0; i < sm_.edges_.size(); i++) {
      int src = sm_.edges_[i].src;
      int dst = sm_.edges_[i].dst;
      int el = sm_.edges_[i].el;
      w2_[i] = cp[edge_idx[src][el][dst]];
    }
  }
}

// Add edge <s, p, o> to summary graph
void SumRDF::AddEdgeToSummary(int s, int p, int o,
    unordered_map<int, unordered_map<int, unordered_map<int, int>>>& idx) {
  if (!idx[s][p].count(o)) { // if a new edge exists in summary
    sm_.edges_.push_back(Edge(s, o, p));
    sm_.multiplicity_++;
    idx[s][p][o] = sm_.edges_.size() - 1;
    w2_.push_back(0);
    sm_.adj_list_[s].push_back(make_pair(o, p));
    sm_.rev_adj_list_[o].push_back(make_pair(s, p));
  }
  w2_[idx[s][p][o]]++;
}

void SumRDF::AddEdgeToSummary(int s, int p, int o,
  std::unordered_set<Edge,EdgeHasher>& idx) {
  auto result = idx.emplace(s,o,p);
  if (result.second) {
    sm_.adj_list_[s].push_back(make_pair(o, p));
    sm_.rev_adj_list_[o].push_back(make_pair(s, p));
  }
}
// Set parameters for MinHashing
// reference reference code: summarisation/factory/minhash/MinHashConfiguration.java
void SumRDF::SummarizeUsingMinHash(DataGraph& g, double ratio) {
  scheme_rows_ = n_ = 40;
  scheme_cols_ = 2;
  m_ = scheme_cols_ + 1;
  threshold_ = std::stod(string(getenv("GCARE_SUMRDF_THRESHOLD")));
  a_.clear();
  b_.clear();
  bucket_lst_.clear();
  signatures_.clear();

  MinHashScheme();
  MakeBucketListByType();
  CreateSummary(g, ratio);
}

// Generate a MinHash scheme of size m x n, m == scheme_rows_, n == scheme_cols_
// reference code: summarisation/factory/minhash/MinHash.java
void SumRDF::MinHashScheme() {
  std::mt19937 generator(0);
  std::uniform_int_distribution<int> dis(-2147483648, 2147483647);
  for (int i = 0; i < scheme_rows_; i++) {
    a_.push_back(vector<int>(scheme_cols_, 0));
    b_.push_back(vector<int>(scheme_cols_, 0));
    for (int j = 0; j < scheme_cols_; j++) {
      a_[i][j] = dis(generator);
      b_[i][j] = dis(generator);
    }
  }
}

// refer to Algorithm 1 of the original paper
// reference code: summarisation/factory/minhash/MinHash.java #81
void SumRDF::CreateSummary(DataGraph& g, double ratio) {

  // Merge similar types (Line 20 of Algorithm 1)
  // It merges types with the same classes

  int before = sm_.multiplicity_;

  CalcSignatureSize();
  vector<bool> erased_bucket(sm_.buckets_.size(), false);
  bool merged = true;
  while (merged) {
    for (size_t t = 0; t < bucket_lst_.size(); t++) {
      auto &bucket_lst = bucket_lst_[t];
      CreateSignature(t);
      auto &signatures = signatures_;
      for (int row = 0; row < scheme_rows_; row++) {
        vector<vector<int>> bins(32768, vector<int>());
        // B_t == bucket_lst_[t]
        for (size_t i = 0; i < bucket_lst.size(); i++) {
          int bid = bucket_lst[i];
          if (erased_bucket[bid]) continue; // Do not consider already removed types
          int type_idx = sm_.buckets_[bid].type_idx_;
          // Compute M^b[i], M^b[i] == val
          int val = signatures[type_idx * n_ * m_ + row * m_ + scheme_cols_];
          // Add b to Bins[LSH(M^b[i])]
          bins[BinHash(val)].push_back(bid);
        }
        for (vector<int>& lst : bins) {
          if (lst.size() == 0) continue;
          int b1 = lst[0];
          for (size_t i = 1; i < lst.size(); ++i) {
            int b2 = lst[i];
            if (Similarity(sm_.buckets_[b1].type_idx_, sm_.buckets_[b2].type_idx_) >= threshold_) 
            {
              erased_bucket[Find(b1) < Find(b2) ? b2 : b1] = true; 
              ShallowMerge(b1, b2);
              merged = true;
              break;
            }
          }       
        }
      }
    }
    // Adjust the summary graph structure after summarization
    sm_.edges_.clear();
    int multiplicity = 0;
    vector<vector<pair<int,int>>> new_adj_list_(sm_.buckets_.size());
    for (size_t t = 0; t < bucket_lst_.size(); t++) {
      std::unordered_set<Edge,EdgeHasher> edge_idx;
      int local_multiplicity = 0;
      auto &bucket_lst = bucket_lst_[t];
      for (size_t i = 0; i < bucket_lst.size(); i++) {
        int s = bucket_lst[i];
        auto &adj_list = sm_.adj_list_[s];
        if (adj_list.size() == 0) continue;
        int ns = Find(s);
        for (size_t j = 0; j < adj_list.size(); j++) {
          auto &pair = adj_list[j];
          int no = Find(pair.first);
          int p = pair.second;
          auto result = edge_idx.emplace(ns, no, p);
          if (result.second) {
            new_adj_list_[ns].push_back(make_pair(no, p));
            ++local_multiplicity;
          }
        }
      }
      multiplicity += local_multiplicity;
    }
    new_adj_list_.swap(sm_.adj_list_);
    vector<vector<pair<int,int>>> new_rev_adj_list_(sm_.buckets_.size());
    for (size_t s = 0; s < sm_.adj_list_.size(); ++s) {
      auto &adj_list = sm_.adj_list_[s];
      if (adj_list.size() == 0) continue;
      int ns = Find(s);
      for (size_t i = 0; i < adj_list.size(); ++i) {
        auto &pair = adj_list[i];
        int no = Find(pair.first);
        int p = pair.second;
        new_rev_adj_list_[no].push_back(make_pair(ns, p));
      }
    }
    new_rev_adj_list_.swap(sm_.rev_adj_list_);
    sm_.multiplicity_ = multiplicity;
    if (multiplicity <= target_) return;
    if (before - multiplicity < 1000) break;
    before = multiplicity;
  }
  int n = static_cast<double>(SummaryGraphSize(g)) / DataGraphSize(g) / ratio;
  n = std::max(n, 1);
  for (size_t t = 0; t < bucket_lst_.size(); t += n) {
    int p = -1;
    for (size_t i = t; i < t + n && i < bucket_lst_.size(); i++) {
      for (size_t j = 0; j < bucket_lst_[i].size(); j++) {
        int bid = bucket_lst_[i][j];
        if (erased_bucket[bid]) continue;
        if (p == -1) p = bid;
        ShallowMerge(p, bid);
      }
    }
  }
}

// Merge similar types (Line 20 of Algorithm 1)
void SumRDF::MakeBucketListByType() {
  // s_vlabel[s] maps type (string)  to the new type number
  // Here, it merges types with the same classes
  unordered_map<string, int> s_vlabel;
  s_vlabel.clear();
  int cnt = 0;
  for (size_t i = 0; i < sm_.buckets_.size(); i++) {
    Bucket& bucket = sm_.buckets_[i];
    Type t = bucket.type_;
    string s = t.ToStringByClass(); // Merge types using only classes
    if (s_vlabel.count(s)) {
      bucket_lst_[s_vlabel[s]].push_back(i);
      bucket.type_idx_ = bucket_lst_[s_vlabel[s]].size() - 1;
      continue;
    }
    s_vlabel[s] = cnt++;
    bucket_lst_.push_back(vector<int>());
    bucket_lst_[cnt - 1].push_back(i);
    bucket.type_idx_ = 0;
  }
}

size_t SumRDF::DataGraphSize(DataGraph &g) {
  size_t vl_cnt = 0;
  for (size_t i = 0; i < g_resources_.size(); i++)
    vl_cnt += g_resources_[i].type_.classes_.size();
  return sizeof (int) * (static_cast<size_t>(g.GetNumVertices()) + static_cast<size_t>(g.GetNumEdges()) * 2 + vl_cnt);
}

size_t SumRDF::SummaryGraphSize(DataGraph& g) {
  size_t v_cnt = 0, vl_cnt = 0;
  for (size_t i = 0; i < sm_.buckets_.size(); i++) {
    if (Find(i) != i) {
      sm_.buckets_[Find(i)].type_.classes_.insert(sm_.buckets_[i].type_.classes_.begin(), sm_.buckets_[i].type_.classes_.end());
    }
  }
  for (size_t i = 0; i < sm_.buckets_.size(); i++) {
    if (Find(i) == i) {
      v_cnt++;
      vl_cnt += sm_.buckets_[i].type_.classes_.size();
    }
  }
  return sizeof (int) * (v_cnt * 2 + static_cast<size_t>(sm_.multiplicity_) * 3 + static_cast<size_t>(g.GetNumVertices()) + vl_cnt);
}

void SumRDF::MakeBucketList() {
  int n = std::stoi(string(getenv("N")));
  vector<vector<int>> cp = bucket_lst_;
  bucket_lst_.clear();
  for (size_t i = 0, i2 = 0; i < cp.size(); i += n, i2++) {
    bucket_lst_.push_back(vector<int>());
    for (size_t j = i; j < i + n && j < cp.size(); j++) {
      for (size_t k = 0; k < cp[j].size(); k++) {
        bucket_lst_[i2].push_back(cp[j][k]);
        sm_.buckets_[cp[j][k]].type_idx_ = bucket_lst_[i].size() - 1;
      }
    }
  }
}

// Compute the size of signature
void SumRDF::CalcSignatureSize() {
  int mx = 0;
  for (size_t t = 0; t < bucket_lst_.size(); t++)
    mx = std::max(mx, static_cast<int>(bucket_lst_[t].size()));
  signatures_.resize(mx * n_ * m_);
}

// reference code: summarisation/factory/minhash/MinHash.java #204 similarity()
void SumRDF::CreateSignature(int t) {
  auto &signatures = signatures_;

  for (size_t i = 0; i < signatures.size(); i++) {
    signatures[i] = std::numeric_limits<int>::max();
  }

  for (size_t i = 0; i < bucket_lst_[t].size(); i++) {
    int bid = bucket_lst_[t][i];
    Type& type = sm_.buckets_[bid].type_;
    for (pair<int, int> e : sm_.adj_list_[bid]) {
      int dstid = e.first;
      int elabel = e.second;
      UpdateSignature(bid, elabel, dstid);
    }
    for (pair<int, int> e : sm_.rev_adj_list_[bid]) {
      int srcid = e.first;
      int elabel = e.second;
      UpdateSignature(bid, srcid, elabel);
    }
    for (int r = 0; r < scheme_rows_; r++) {
      int index = sm_.buckets_[bid].type_idx_ * n_ * m_ + r * m_;
      int hash = 1;
      int c = 0;
      for (c = 0; c < scheme_cols_; c++)
        hash = 31 * hash + signatures[index + c];
      signatures[index + c] = hash;
    }
  }
}

// reference code: summarisation/factory/minhash/MinHash.java #239 minHash()
void SumRDF::UpdateSignature(int bid, int bid1, int bid2) {
  auto &signatures = signatures_;
  for (int r = 0; r < scheme_rows_; r++) {
    int index = sm_.buckets_[bid].type_idx_ * n_ * m_ + r * m_;
    for (int c = 0; c < scheme_cols_; c++) {
      int hash = SchemeHash(bid1, bid2, r, c);
      if (signatures[index + c] > hash) signatures[index + c] = hash;
    }
  }
}

// reference code: summarisation/factory/minhash/MinHash.java #244 in minHash()
int SumRDF::SchemeHash(int x, int y, int row, int col) {
  int val = (x * 127 + y) * 31;
  return (val * a_[row][col]) + b_[row][col]; // global variable
}

// LSH function == BinHash function
int SumRDF::BinHash(int x) {
  int res = x ^ (static_cast<uint32_t>(x) >> 16);
  return (res & 0x7fffffff) & 32767;
}

// reference code: summarisation/factory/minhash/MinHash.java #188 similarity()
double SumRDF::Similarity(int b1, int b2) {
  auto &signatures = signatures_;
  double res = 0.0;
  for (int row = 0; row < scheme_rows_; row++) {
    int idx1 = b1 * n_ * m_ + row * m_;
    int idx2 = b2 * n_ * m_ + row * m_;
    for (int col = 0; col < scheme_cols_; col++) {
      if (signatures[idx1 + col] == signatures[idx2 + col]) res += 1.0;
    }
  }
  return res / (scheme_rows_ * scheme_cols_);
}

// Merge two types
void SumRDF::ShallowMerge(int b1, int b2) {
  if (b1 > b2) {
    ShallowMerge(b2, b1);
    return;
  }
  Union(b1, b2);
}

// Initialize a summary graph structure
SumRDF::Summary::Summary() {
  buckets_.clear();
  edges_.clear();
  multiplicity_ = ori_multiplicity_ = 0;
  build_time_ = 0.0;
}

// Add a new bucket to a bucket list
void SumRDF::Summary::CreateBucket(Resource& r, int rid) {
  buckets_.push_back(Bucket(r.type_));
}

void SumRDF::WriteSummaryFile(const char* fn) {
  string fn2 = string(fn) + ".summary";
  FILE* fp = fopen(fn2.c_str(), "w");
  for (size_t i = 0; i < sm_.buckets_.size(); i++) {
    Bucket& bucket = sm_.buckets_[i];
    fprintf(fp, "b %d", i); 
    for (int rsc : bucket.resources_) {
      fprintf(fp, " %d", rsc);
    }
    fprintf(fp, "\n");
  }
  for (size_t i = 0; i < sm_.edges_.size(); ++i) {
    fprintf(fp, "e %d %d %d %d\n", sm_.edges_[i].src, sm_.edges_[i].dst, sm_.edges_[i].el, w2_[i]);
  }
  fclose(fp);
}

void SumRDF::WriteSummary(const char* fn) {
  WriteSummaryFile(fn);
  FILE* fp = fopen(fn, "w");
  fprintf(fp, "t # 0\n");
  for (size_t i = 0; i < sm_.buckets_.size(); i++) {
    Bucket& bucket = sm_.buckets_[i];
    Type t = bucket.type_;
    string s = t.ToStringByClass();
    fprintf(fp, "v %d", static_cast<int>(i));
    if (bucket.type_.classes_.size() == 0) fprintf(fp, " -1");
    for (int vl : bucket.type_.classes_) fprintf(fp, " %d", vl);
    fprintf(fp, "\n");
  }
  for (size_t i = 0; i < sm_.edges_.size(); i++) {
    fprintf(fp, "e %d %d %d\n", sm_.edges_[i].src, sm_.edges_[i].dst, sm_.edges_[i].el);
  }
  fclose(fp);
  string fn2 = string(fn) + ".weight";
  FILE* fp2 = fopen(fn2.c_str(), "w");
  fprintf(fp2, "%d\n", static_cast<int>(w1_.size()));
  fwrite(w1_.data(), sizeof (int), w1_.size(), fp2);
  fprintf(fp2, "\n%d\n", static_cast<int>(w2_.size()));
  fwrite(w2_.data(), sizeof (int), w2_.size(), fp2);
  fprintf(fp2, "\n%d\n", static_cast<int>(sm_.buckets_.size()));
  for (size_t i = 0; i < sm_.buckets_.size(); i++) {
    Bucket& bucket = sm_.buckets_[i];
    fprintf(fp2, "%d\n", static_cast<int>(bucket.resources_.size()));
    for (int r : bucket.resources_) fprintf(fp2, "%d ", r);
    fprintf(fp2, "\n");
  }
  fclose(fp2);
  fn2 = string(fn) + ".bin";
  s_.ReadText(fn);
  s_.MakeBinary();
  s_.WriteBinary(fn2.c_str());
  s_.ClearRawData();
}

void SumRDF::ReadSummary(const char* fn) {
  string fn2 = string(fn) + ".bin";
  s_.ReadBinary(fn2.c_str());
  fn2 = string(fn) + ".weight";
  FILE* fp = fopen(fn2.c_str(), "r");
  int size, size2; char c;
  fscanf(fp, "%d%c", &size, &c);
  w1_.resize(0); w1_.resize(size);
  if (fread(w1_.data(), sizeof (int), size, fp)) { }
  fscanf(fp, "%c%d%c", &c, &size, &c);
  w2_.resize(0); w2_.resize(size);
  if (fread(w2_.data(), sizeof (int), size, fp)) { }
  fscanf(fp, "%d", &size);
  s_buckets_.clear(); s_buckets_.resize(size, Bucket());
  for (size_t i = 0; i < size; i++) {
    int vl, el, rid;
    fscanf(fp, "%d", &size2);
    assert(size2 >= 0);
    while (size2--) {
      fscanf(fp, "%d", &rid);
      s_buckets_[i].resources_.push_back(rid);
    }
  }
  fclose(fp);
}

// reference code: queryanswering/SPARQLEvaluator.java #38 SPARQLEvaluator()
void SumRDF::Init() {
  // Set types for each query vertex
  for (size_t i = 0; i < s_buckets_.size(); i++) {
    for (int srcid : s_buckets_[i].resources_) {
      for (auto r = g->GetVLabels(srcid); r.begin != r.end; r.begin++) {
        int vl = *r.begin;
        s_buckets_[i].type_.classes_.insert(vl);
      }
      for (auto re = g->GetELabels(srcid, true); re.begin != re.end; re.begin++) {
        int el = *re.begin;
        s_buckets_[i].type_.outgoing_.insert(el);
      }
      for (auto re = g->GetELabels(srcid, false); re.begin != re.end; re.begin++) {
        int el = *re.begin;
        s_buckets_[i].type_.incoming_.insert(el);
      }
    }
  }
  q_resources_.clear(); q_resources_.resize(q->GetNumVertices(), Resource());
  for (int i = 0; i < q->GetNumVertices(); i++) {
    if (q->GetVLabel(i) != -1) q_resources_[i].type_.classes_.insert(q->GetVLabel(i));
    for (auto& r : q->GetAdj(i, true)) {
      int dst = r.first;
      int el = r.second;
      // O == type_.outgoing_
      q_resources_[i].type_.outgoing_.insert(el);   
      // I == type_.incoming_
      q_resources_[dst].type_.incoming_.insert(el);
    }
    q_resources_[i].bound_ = q->GetBound(i);
  }
  // Set summary graph (S == smo_, H == smo_.data_edges)
  smo_.prev.clear(); smo_.prev.resize(q->GetNumEdges(), -1);
  smo_.data_edges.clear();
  for (int srcid = 0; srcid < s_.GetNumVertices(); srcid++) {
    for (auto re = s_.GetELabels(srcid, true); re.begin != re.end; re.begin++) {
      assert(re.begin!=re.end);
      int el = *re.begin;
      for (auto rd = s_.GetAdj(srcid, el, true); rd.begin != rd.end; rd.begin++) {
        int dstid = *rd.begin;
        smo_.data_edges.push_back(Edge(srcid, dstid, el));
      }
    }
  }
  result_ = 0.0;
  emb_.resize(0);
  tau_.resize(0); tau_.resize(q->GetNumVertices(), -1);
  smo_.candidates.clear();
  partitions_.resize(0);
  CreateBasePartition();
  edges_.resize(0); edges_.resize(q->GetNumEdges());
  indexes_.resize(0);
  is_init_ = false;
  pos_ = -1;
  // Set candidate edges for each query edge for subgraph matching
  // candidates[i]: a set of summary edges which can be matched with i-th query edge
  for (size_t i = 0; i < q->GetNumEdges(); i++) {
    Edge qedge = q->GetEdge(i);
    smo_.candidates.push_back(vector<int>());
    for (size_t j = 0; j < smo_.data_edges.size(); j++) {
      Edge sedge = smo_.data_edges[j];
      bool flag = true;
      if (q->GetBound(qedge.src) != -1) {
        flag = false;
        assert(sedge.src < s_buckets_.size());
        for (int rid : s_buckets_[sedge.src].resources_) {
          if (rid == q->GetBound(qedge.src)) {
            flag = true;
            break;
          }
        }
      }
      if (q->GetBound(qedge.dst) != -1) {
        flag = false;
        assert(sedge.dst < s_buckets_.size());
        for (int rid : s_buckets_[sedge.dst].resources_) {
          if (rid == q->GetBound(qedge.dst)) {
            flag = true;
            break;
          }
        }
      }
      if (!flag) continue;
      if (Includes(s_buckets_[sedge.src].type_, q_resources_[qedge.src].type_) &&
          Includes(s_buckets_[sedge.dst].type_, q_resources_[qedge.dst].type_) &&
          sedge.el == qedge.el) { // Check conditions for subgraph matching
        smo_.candidates[i].push_back(j);
      }
    }
  }
}


int SumRDF::DecomposeQuery() {
  return 1;
}

// 4.2 Formalisation
// A partitoin groups the atoms (query edges) of q into disjoint sets, where all atoms from each such group are mapped to the same triple
// reference code: queryanswering/SPARQLEvaluator.java #186 PartitionBase()
void SumRDF::CreateBasePartition() {
  blocks_.clear(); blocks_.resize(q->GetNumEdges(), vector<int>());
  ComputePartitions();
  // a set of partitions == partitions_ and P == partition
  for (Partition& partition : partitions_) {
    for (size_t i = 0; i < partitions_.size(); i++) {
      Partition finer = partitions_[i];
      if (partition.IsFiner(finer)) partition.finer_partitions_.push_back(i);
    }
  }
}

// reference code: queryanswering/SPARQLEvaluator.java #213 computePartitions()
void SumRDF::ComputePartitions() {
  vector<int> p(q->GetNumEdges(), 0);
  vector<int> m(q->GetNumEdges(), 0);
  TrivialPartition(p);
  while (true) {
    int index = FindDecrementable(p);
    if (index < 1) break;
    p[index]--;
    m[index] = std::max(p[index - 1], m[index - 1]);
    for (size_t i = index + 1; i < p.size(); i++) {
      m[i] = std::max(p[i - 1], m[i - 1]);
      p[i] = m[i] + 1;
    }
    int count_partitions = std::max(p[p.size() - 1], m[p.size() - 1]) + 1;
    Partition partition = Partition(p, count_partitions, blocks_);
    if (IsUnifiable(partition, p, count_partitions)) { // Check whether we can make a new partition from the query
      partition.blocks_ = blocks_;
      partition.ComputeQuery(*q);
      int sum = 0; for (vector<int> b: partition.blocks_) assert(b.size()>0);
      partitions_.push_back(partition);
    }
  }
}

// Compute the query from a partition made by query
void SumRDF::Partition::ComputeQuery(QueryGraph& q) {
  q_mgu_.clear();
  for (int i = 0; i < q.GetNumEdges(); i++) {
    Edge edge = q.GetEdge(i);
    Edge new_edge = Edge(gamma_[edge.src], gamma_[edge.dst], edge.el);
    bool flag = false;
    for (int j = 0; j < q_mgu_.size(); j++) {
      if (new_edge == q_mgu_[j]) {
        flag = true;
        break;
      }
    }
    if (flag) continue;
    q_mgu_.push_back(new_edge);
  }
}

// reference code: queryanswering/SPARQLEvaluator.java #253
void SumRDF::TrivialPartition(vector<int>& p) {
  for (size_t i = 0; i < p.size(); i++) {
    p[i] = i;
    blocks_[i].push_back(i);
  }
  Partition partition = Partition(p, p.size(), blocks_);
  partition.gamma_.clear(); partition.gamma_.resize(q->GetNumVertices(), 0);
  for (int i = 0; i < q->GetNumVertices(); i++) {
    partition.gamma_[i] = i;
  }
  partition.ComputeQuery(*q);
  partitions_.push_back(partition);
}

// reference code: queryanswering/SPARQLEvaluator.java #244 findDecrementable()
int SumRDF::FindDecrementable(vector<int>& p) {
  for (int i = static_cast<int>(p.size()) - 1; i >= 0; i--)
    if (p[i] > 0) return i;
  return -1;
}

// Check whether one partition contains the other partition
// reference code: queryanswering/SPARQLEvaluator.java #418 isFiner()
bool SumRDF::Partition::IsFiner(Partition& finer) {
  if (finer.blocks_.size() < blocks_.size()) {
    for (vector<int>& block : blocks_) {
      bool flag = false;
      for (vector<int>& finer_block : finer.blocks_) {
        if (IsContained(block, finer_block)) flag = true;
      }
      if (!flag) return false;
    }
    return true;
  }
  return false;
}

// reference code: queryanswering/reasoner/SPARQLEvaluator.java isContained()
bool SumRDF::Partition::IsContained(vector<int>& block, vector<int>& finer_block) {
  for (int e : block) {
    bool flag = false;
    for (int e2 : finer_block) {
      if (e == e2) {
        flag = true;
        break;
      }
    }
    if (!flag) return false;
  }
  return true;
}

// Check whether type t1 includes type t2 (matching condition)
bool SumRDF::Includes(Type& t1, Type& t2) {
  return std::includes(t1.classes_.begin(), t1.classes_.end(),
      t2.classes_.begin(), t2.classes_.end()) &&
    std::includes(t1.outgoing_.begin(), t1.outgoing_.end(),
        t2.outgoing_.begin(), t2.outgoing_.end()) &&
    std::includes(t1.incoming_.begin(), t1.incoming_.end(),
        t2.incoming_.begin(), t2.incoming_.end());
}


// GetSubstructure return an embedding of Q in S
// reference code: queryanswering/SPARQLEvaluator.java #91 evaluate()
bool SumRDF::GetSubstructure(int subquery_index) {
  if (!is_init_) {
    is_init_ = true;
  }
  smo_.Init(*g, *q);
  if (!smo_.FindMatchedSubgraph(0)) { // if we cannot find any subgraph, then stop it
    return false;
  }
  tau_ = smo_.embedding;
  edges_ = smo_.edge_idx;
  return true;
}

// 4.2 Formalisation
// reference code: queryanswering/SPARQLEvaluator.java #109 Process()
double SumRDF::EstCard(int subquery_index) {
  result_ = 0.0;
  query_image_.clear();
  // query_image_[i]: i-th query edge
  // Here, if query edges are mapped to the same summary edge, it considers only one query edge
  for (size_t i = 0; i < q->GetNumEdges(); i++) {
    bool flag = false;
    for (size_t j = 0; j < i; j++) {
      if (edges_[i] == edges_[j]) {
        flag = true;
        break;
      }
    }
    if (flag) continue;
    query_image_.push_back(edges_[i]);
  }

  int q_size = q->GetNumVertices();
  int q_remaining = q_size % 8;
  q_size = q_remaining == 0 ? q_size : 8 + q_size - q_remaining;
  solution_chk = (bool*) malloc(q_size);
  int pid = 0;
  for (Partition partition : partitions_) {
    if (!partition.IsTauUnifiable(tau_)) continue;
    // P(m) == p
    double cnt = MinimalSolutions(partition);
    double p = QueryFactor(partition);
    // THEOREM 4.3. E_{q,s} == results_
    result_ += p * cnt;
  }
  free(solution_chk);
  return result_;
}

// reference code: queryanswering/SPARQLEvaluator.java #144 queryFactor()
double SumRDF::QueryFactor(Partition& partition) {
  double res = 1.0;
  for (size_t i = 0; i < query_image_.size(); i++) {
    res *= AtomFactor(partition, query_image_[i]);
  }
  return res;
}

// reference code: queryanswering/SPARQLEvaluator.java #126 minimalSolutions()
double SumRDF::MinimalSolutions(Partition& partition) {
  double sum = 0.0;
  for (int idx : partition.finer_partitions_) {
    Partition& finer = partitions_[idx];
    double val = MinimalSolutions(finer);
    if (finer.IsTauUnifiable(tau_)) {
      sum += val;
    }
  }
  double x = Solutions(partition);

  return x - sum;
}

// reference code: queryanswering/SPARQLEvaluator.java #136 solutions()
inline double SumRDF::Solutions(Partition& partition) {
  double res = 1;
  bool* chk = solution_chk;
  uint64_t* vec = (uint64_t*) chk;
  for (int i = 0; i < q->GetNumVertices(); i+= 8) {
      vec[i>>3] = 0;
  }
  for (Edge &edge : partition.q_mgu_) {
    int srcid = edge.src;
    int dstid = edge.dst;
    if (!chk[srcid] && q->GetBound(srcid) == -1) {
      res *= w1_[tau_[srcid]];  // s == w1_
      chk[srcid] = true;
    }
    if (!chk[dstid] && q->GetBound(dstid) == -1) {
      res *= w1_[tau_[dstid]];
      chk[dstid] = true;
    }
  }
  return (double) res;
}

// reference code: queryanswering/SPARQLEvaluator.java #152 atomFactor()
double SumRDF::AtomFactor(Partition& partition, int idx) {
  int cnt = GetWeight(partition, idx);
  if (cnt > w2_[idx]) { // w == w2_
    return 0.0;
  }
  double size = static_cast<double>(w1_[smo_.data_edges[idx].src]) *
    static_cast<double>(w1_[smo_.data_edges[idx].dst]);
  double res = 1.0;
  for (int i = 0; i < cnt; i++) {
    assert(size - i > 0);
    res *= static_cast<double>(w2_[idx] - i) / static_cast<double>(size - i);
  }
  return res;
}

// reference code: queryanswering/SPARQLEvaluator.java #164 preimage()
int SumRDF::GetWeight(Partition& partition, int idx) {
  int size = 0;
  for (Edge atom : partition.q_mgu_) {
    assert(atom.src < tau_.size());
    assert(atom.dst < tau_.size());
    assert(idx < smo_.data_edges.size());
    if (tau_[atom.src] == smo_.data_edges[idx].src && tau_[atom.dst] == smo_.data_edges[idx].dst && atom.el == smo_.data_edges[idx].el) size++;
  }
  return size;
}

// reference code: queryanswering/SPARQLEvaluator.java #400
bool SumRDF::Partition::IsTauUnifiable(vector<int>& tau) {
  for (size_t i = 0; i < gamma_.size(); i++) {
    if (tau[i] != tau[gamma_[i]]) return false;
  }
  return true;
}

// Check if we make a new partition from the query
// reference code: queryanswering/SPARQLEvaluator.java #273
bool SumRDF::IsUnifiable(Partition& partition, vector<int>& p, int count_partitions_) {
  vector<int> subject_for_partition(count_partitions_, -1);
  vector<int> predicate_for_partition(count_partitions_, -1);
  vector<int> object_for_partition(count_partitions_, -1);
  blocks_.clear();
  blocks_.resize(count_partitions_, vector<int>());
  partition.gamma_.resize(q_resources_.size(), 0);
  for (size_t i = 0; i < partition.gamma_.size(); i++)
    partition.gamma_[i] = i;
  assert(p.size() == q->GetNumEdges());
  for (size_t i = 0; i < p.size(); i++) {
    int partition_index = p[i];
    blocks_[partition_index].push_back(i);
    Edge atom = q->GetEdge(i);
    int subject_index = atom.src;
    if (subject_for_partition[partition_index] < 0)
      subject_for_partition[partition_index] = subject_index;
    int predicate_val = atom.el;
    if (predicate_for_partition[partition_index] < 0)
      predicate_for_partition[partition_index] = predicate_val;
    int object_index = atom.dst;
    if (object_for_partition[partition_index] < 0)
      object_for_partition[partition_index] = object_index;
    int subject_partition_index = partition.gamma_[subject_for_partition[partition_index]];
    Resource subject_partition_resource = q_resources_[subject_partition_index];
    int subject_atom_index = partition.gamma_[subject_index];
    Resource subject_atom_resource = q_resources_[subject_atom_index];
    if (subject_atom_resource.bound_ == -1 && subject_partition_resource.bound_ == -1) {
      ReplaceGamma(partition, subject_atom_index, subject_partition_index);
    } else if (subject_atom_resource.bound_ == -1) {
      ReplaceGamma(partition, subject_atom_index, subject_partition_index);
    } else if (subject_partition_resource.bound_ == -1) {
      ReplaceGamma(partition, subject_partition_index, subject_atom_index);
    } else {
      return false;
    }
    if (predicate_val != predicate_for_partition[partition_index]) {
      return false;
    }
    int object_partition_index = partition.gamma_[object_for_partition[partition_index]];
    Resource object_partition_resource = q_resources_[object_partition_index];
    int object_atom_index = partition.gamma_[object_index];
    Resource object_atom_resource = q_resources_[object_atom_index];
    if (object_atom_resource.bound_ == -1 && object_partition_resource.bound_ == -1) {
      ReplaceGamma(partition, object_atom_index, object_partition_index);
    } else if (object_atom_resource.bound_ == -1) {
      ReplaceGamma(partition, object_atom_index, object_partition_index);
    } else if (object_partition_resource.bound_ == -1) {
      ReplaceGamma(partition, object_partition_index, object_atom_index);
    } else {
      return false;
    }
  }
  return true;
}

// reference code: queryanswering/SPARQLEvaluator.java #352 replaceInMgu
void SumRDF::ReplaceGamma(Partition& partition, int i, int j) {
  for (size_t k = 0; k < partition.gamma_.size(); k++)
    if (partition.gamma_[k] == i) partition.gamma_[k] = j;
}

double SumRDF::AggCard() {
  double res = 0.0;
  for (double card : card_vec_)
    res += card;
  return res;
}

double SumRDF::GetSelectivity() {
  return 1.0;
}

