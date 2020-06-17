#ifndef GRAPH_OP_
#define GRAPH_OP_

#include "data_graph.h"
#include "query_graph.h"

struct SubgraphMatching {
  DataGraph* g;
  QueryGraph* q;
  vector<Edge> data_edges; // s_edges_
  vector<vector<int>> candidates; // -> iterators_
  vector<int> embedding, cur_embedding; // -> tau_
  vector<int> edge_idx, cur_edge_idx, prev, cur_idx; // -> edges_
  bool flag;
  SubgraphMatching() {
    g = NULL;
    q = NULL;
    data_edges.clear();
    candidates.clear();
    embedding.clear();
    cur_embedding.clear();
    edge_idx.clear();
    cur_edge_idx.clear();
    prev.clear();
    cur_idx.clear();
    flag = false;
  }
  void Init(DataGraph& g_, QueryGraph& q_) {
    g = &g_;
    q = &q_;
    flag = false;
    embedding.clear(); embedding.resize(q->GetNumVertices(), -1);
    cur_embedding.clear(); cur_embedding.resize(q->GetNumVertices(), -1);
    edge_idx.clear(); edge_idx.resize(candidates.size(), -1);
    cur_edge_idx.clear(); cur_edge_idx.resize(candidates.size(), -1);
    cur_idx.clear(); cur_idx.resize(candidates.size(), -1);
  }
  bool FindMatchedSubgraph(int i) {
    if (i == candidates.size()) {
      if (!flag) {
        return false;
      }
      embedding = cur_embedding;
      edge_idx = cur_edge_idx;
      prev = cur_idx;
      return true;
    }
    int src = q->GetEdge(i).src;
    int dst = q->GetEdge(i).dst;
    int el = q->GetEdge(i).el;
    for (int j = flag ? 0 : std::max(prev[i], 0); j < static_cast<int>(candidates[i].size()); j++) {
      int d_src = data_edges[candidates[i][j]].src;
      int d_dst = data_edges[candidates[i][j]].dst;
      int d_el = data_edges[candidates[i][j]].el;
      if (cur_embedding[src] != -1 && cur_embedding[src] != d_src) continue;
      if (cur_embedding[dst] != -1 && cur_embedding[dst] != d_dst) continue;
      if (el != d_el) continue;
      if (!flag && prev[i] < j) flag = true;
      int prev1 = cur_embedding[src];
      int prev2 = cur_embedding[dst];
      cur_embedding[src] = d_src;
      cur_embedding[dst] = d_dst;
      cur_edge_idx[i] = candidates[i][j];
      cur_idx[i] = j;
      if (FindMatchedSubgraph(i + 1)) return true;
      cur_embedding[src] = prev1;
      cur_embedding[dst] = prev2;
    }
    return false;
  }
};
#endif
