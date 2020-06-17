#include "../include/impr.h"
#include "../include/data_graph.h"
#include "../include/query_graph.h"
#include "../include/estimator.h"

#include <vector>
#include <unordered_map>

using std::vector;
using std::unordered_map;
using std::pair;

void Impr::PrepareSummaryStructure(DataGraph& g, double p) {
}

void Impr::WriteSummary(const char* fn) {
}

// Initialize variables
void Impr::Init() {
  // Impr can process 3,4,5 node queries
  if (q->GetNumVertices() > 5) return;
  pos_embs_.clear();
	for (int start = 0; start < q->GetNumVertices(); start++) {
		chk_.clear();
		chk_.resize(q->GetNumVertices(), false);
		qv_.clear();
		DFS(start);
	}
  std::sort(pos_embs_.begin(), pos_embs_.end());
  pos_embs_.erase(std::unique(pos_embs_.begin(), pos_embs_.end()), pos_embs_.end());
  for (int u = 0; u < q->GetNumVertices(); u++) {
    for (auto p : q->GetAdj(u, true)) {
      query_labels_.push_back(make_pair(p.second, true));
      query_labels_.push_back(make_pair(p.second, false));
    }
  }
  std::sort(query_labels_.begin(), query_labels_.end());
  query_labels_.erase(std::unique(query_labels_.begin(), query_labels_.end()), query_labels_.end());
  xk_.clear();
  el_.clear();
  case_num_.clear();
  int sum = 0;
  for (size_t u = 0; u < q->GetNumVertices(); u++) {
    for (auto& p : q->GetAdj(u, true)) {
      int elabel = p.second;
      sum += g->GetNumEdges(elabel);
    }
  }
  // Compute the target number of steps as the sampling budget
  s_num_ = sample_ratio * sum;
  chk_.clear(); chk_.resize(q->GetNumVertices(), false);
  path_.clear();
  // Compute the beta value
  beta_ = GetBeta(0);
  // step_ is the current number of steps during random walks
  step_ = 0;
}

void Impr::ReadSummary(const char* fn) {
}

// Compute beta value by traversing the query graph
int Impr::GetBeta(int cnt) {
	if (cnt == q->GetNumVertices() - 1) {
		for (int i = 1; i < cnt; i++) {
      bool flag = false;
      for (auto p : q->GetAdj(path_[i - 1], true)) {
        if (p.first == path_[i]) {
          flag = true;
          break;
        }
      }
      for (auto p : q->GetAdj(path_[i], true)) {
        if (p.first == path_[i - 1]) {
          flag = true;
          break;
        }
      }
      if (!flag) return 0;
		}
		return 1;
	}
	int res = 0;
	for (int i = 0; i < q->GetNumVertices(); i++) {
		if (chk_[i]) continue;
		path_.push_back(i);
		chk_[i] = true;
		res += GetBeta(cnt + 1);
		path_.pop_back();
		chk_[i] = false;
	}
	return res;
}

int  Impr::DecomposeQuery() {
  return 1;
}

bool Impr::GetSubstructure(int subgraph_idx) {
  // Impr can processes 3,4,5 node queries
  if (q->GetNumVertices() > 5 || beta_ == 0) {
    return false;
  }
  if (el_.size() > 0) el_.erase(el_.begin());
  if (xk_.size() > 0) xk_.erase(xk_.begin());
  if (case_num_.size() > 0) case_num_.erase(case_num_.begin());
  // Compute s_{i+1} from s_i, s_i == xk_
  if (!GetNextSample(xk_, case_num_, el_, step_, s_num_)) {
      return false;
  }
  // if it consumes all sampling budget, stop random walks
  if (step_ >= s_num_) {
      return false;
  }
  return true;
}

void Impr::DFS(int srcid) {
	if (chk_[srcid]) return;
	chk_[srcid] = true;
	qv_.push_back(srcid);
	if (static_cast<int>(qv_.size()) == q->GetNumVertices() - 1) {
		int sum = 0;
		for (int x : qv_)
			sum += x;
		pos_embs_.push_back(qv_);
		pos_embs_[pos_embs_.size() - 1].push_back((q->GetNumVertices() - 1) * q->GetNumVertices() / 2 - sum);
		return;
	}
	for (pair<int, int>& e : q->GetAdj(srcid, true)) {
		int dstid = e.first;
		int elabel = e.second;
		if (!chk_[dstid]) {
			DFS(dstid);
			chk_[dstid] = false;
			qv_.pop_back();
		}
	}
	for (pair<int, int>& e : q->GetAdj(srcid, false)) {
		int dstid = e.first;
		int elabel = e.second;
		if (!chk_[dstid]) {
			DFS(dstid);
			chk_[dstid] = false;
			qv_.pop_back();
		}
	}
}

// Choose the next edge label to walk further
pair<int, bool> Impr::ChooseELabel(int v, vector<pair<int, bool>>& cand, int& res, int& sum) {
	sum = 0;
	for (auto& p : cand) {
    if (!g->HasELabel(v, p.first, p.second)) continue;
    sum += g->GetAdjSize(v, p.first, p.second);
	}
	if (sum == 0) return make_pair(-1, -1);
  int n = rand() % sum;
	for (auto& p : cand) {
    if (!g->HasELabel(v, p.first, p.second)) continue;
    n -= g->GetAdjSize(v, p.first, p.second);
    if (n < 0) {
      res = *(g->GetAdj(v, p.first, p.second).begin + n + g->GetAdjSize(v, p.first, p.second));
      return make_pair(p.first, p.second);
    }
	}
	return make_pair(-1, -1);
}

double Impr::EstCard(int subgraph_index) {
  assert(static_cast<int>(xk_.size()) == q->GetNumVertices() - 1);
  is_duplicate_.clear();
  int f = 0;
  for (size_t i = 0; i < pos_embs_.size(); i++) {
    bool dir = false;
    int sum = 0, selected_el = -1;
    int selected_v = Select(xk_, pos_embs_[i], dir, selected_el);
    if (selected_v == -1) continue;
    assert(selected_el != -1);
    if (!g->HasELabel(selected_v, selected_el, dir)) continue;
    for (auto r = g->GetAdj(selected_v, selected_el, dir); r.begin != r.end; r.begin++) {
      int nbr = *r.begin;
      xk_.push_back(nbr);
      if (Check(xk_, pos_embs_[i], el_)) f++;
      step_++;
      xk_.pop_back();
    }
  }
  double inv_w = GetWeight(xk_) / (2 * g->GetNumEdges()) * beta_;
  assert(beta_ > 0);
  assert(inv_w > 0);
  return static_cast<double>(f) / inv_w;
}

// Compute s_{i+1} from s_i, s_i == xk
bool Impr::GetNextSample(vector<int>& xk, vector<int>& case_num, vector<pair<int, bool>>& el, int& step, int s_num) {
	if (el.size() == 0) el.clear();
	int sum = 0, cnt = 0;
  // if there is no vertex in s_i, insert a vertex chosen randomly
	if (xk.size() == 0) {
    int x = rand() % g->GetNumVertices();
		xk.push_back(x);
		case_num.push_back(-1);
		step++;
	}
	while (static_cast<int>(xk.size()) < q->GetNumVertices() - 1 && step < s_num) {
		xk.push_back(-1);
    // Choose edge label to walk further
		el.push_back(ChooseELabel(xk[xk.size() - 2], query_labels_, xk[xk.size() - 1], sum));
		case_num.push_back(sum);
		step++;
		if (el[el.size() - 1].first == -1) {
			assert(xk.size() > 1);
			xk.pop_back();
			case_num.pop_back();
			el.pop_back();
			if (xk.size() > 0) xk.erase(xk.begin());
			if (case_num.size() > 0) case_num.erase(case_num.begin());
			if (el.size() > 0) el.erase(el.begin());
			if (xk.size() == 0) {
				cnt++;
				if (cnt > 0.1 * g->GetNumVertices()) {
					return false;
				}
				assert(g->GetNumVertices() > 0);
				xk.push_back(rand() % g->GetNumVertices());
				step++;
				case_num.push_back(-1);
			}
		}
	}
	return true;
}

// Check the matching conditions
bool Impr::Check(vector<int>& v, vector<int>& u, vector<pair<int, bool>>& el) {
	vector<int> mapping(u.size());
  int cnt = 0;
  vector<bool> chk(v.size() - 2, false);
	for (size_t i = 0; i < u.size(); i++) {
		mapping[u[i]] = v[i];
	}
  if (is_duplicate_[mapping]) return false;
	for (size_t i = 0; i < q->GetNumVertices(); i++) {
    // Check binded data vertex
		if (q->GetBound(i) != -1 && mapping[i] != q->GetBound(i)) return false;
    // Check vertex label
		if (q->GetVLabel(i) != -1
      && !binary_search(g->GetVLabels(mapping[i]).begin, g->GetVLabels(mapping[i]).end, q->GetVLabel(i))) return false;
		for (auto& p : q->GetAdj(i, true)) {
			int dstid = p.first;
			int elabel = p.second;
			int from = mapping[i];
			int to = mapping[dstid];
      if (!binary_search(g->GetAdj(from, elabel, true).begin, g->GetAdj(from, elabel, true).end, to)) return false;
      for (int j = 0; j < static_cast<int>(v.size()) - 2; j++) {
        if (chk[j]) continue;
        int elabel = el[j].first;
        bool dir = el[j].second;
        if ((dir && v[j] == from && v[j + 1] == to
          && binary_search(g->GetAdj(v[j], elabel, true).begin,
          g->GetAdj(v[j], elabel, true).end, v[j + 1]))
          || (!dir && v[j] == to && v[j + 1] == from
          && binary_search(g->GetAdj(v[j + 1], elabel, true).begin,
          g->GetAdj(v[j + 1], elabel, true).end, v[j]))) {
          chk[j] = true;
          cnt++;
          break;
        }
      }
		}
	}
  if (cnt != chk.size()) return false;
  return (is_duplicate_[mapping] = true);
}

// Select the vertex to refer its adjacency list
int Impr::Select(vector<int>& xk, vector<int>& qv, bool& dir, int& selected_el) {
	int res = -1;
	int to = qv[qv.size() - 1];
	size_t min = 999999999;
	for (size_t i = 0; i < xk.size(); i++) {
		int from = qv[i];
		if (q->GetELabel(from, to) != -1) { // dir == true
			if (min > g->GetAdjSize(xk[i], q->GetELabel(from, to), true) &&
					g->GetAdjSize(xk[i], q->GetELabel(from, to), true) > 0) {
				min = g->GetAdjSize(xk[i], q->GetELabel(from, to), true);
				selected_el = q->GetELabel(from, to);
				dir = true;
				res = xk[i];
			}
		} else if (q->GetELabel(to, from) != -1) {
			if (min > g->GetAdjSize(xk[i], q->GetELabel(to, from), false) &&
					g->GetAdjSize(xk[i], q->GetELabel(to, from), false) > 0) {
				min = g->GetAdjSize(xk[i], q->GetELabel(to, from), false);
				selected_el = q->GetELabel(to, from);
				dir = false;
				res = xk[i];
			}
		}
	}
	return res;
}

// Compute the weight W(s_i)
double Impr::GetWeight(vector<int>& xk) {
  chk_.clear(); chk_.resize(xk.size(), false);
  sum_ = 0.0;
  size_ = 0;
  for (size_t start = 0; start < xk.size(); start++) {
    std::fill(chk_.begin(), chk_.end(), false);
    FindPaths(start, xk, 0, 1.0);
  }
  assert(size_ > 0);
  return static_cast<double>(sum_) / size_;
}

void Impr::FindPaths(int i, vector<int>& xk, int cnt, double pr) {
  if (cnt == q->GetNumVertices() - 2) {
    size_++;
    sum_ += pr;
    return;
  }
  chk_[i] = true;
  int num_nbrs = 0;
  for (auto p : query_labels_) {
    int el = p.first;
    bool dir = p.second;
    num_nbrs += g->GetAdjSize(xk[i], el, dir);
  }
  if (cnt == 0) num_nbrs = 1;
  for (auto p : query_labels_) {
    int el = p.first;
    bool dir = p.second;
    for (size_t j = 0; j < chk_.size(); j++) {
      if (chk_[j] || !binary_search(g->GetAdj(xk[i], el, dir).begin,
        g->GetAdj(xk[i], el, dir).end,
        xk[j])) continue;
      assert(num_nbrs > 0);
      FindPaths(j, xk, cnt + 1, pr * 1.0 / num_nbrs);
      chk_[j] = false;
    }
  }
}

// Impr uses average for aggregation
double Impr::AggCard() {
  if (q->GetNumVertices() > 5 || beta_ == 0) return -1.0;
  if (card_vec_.size() == 0) return 0.0;
  double res = 0.0;
  for (double card : card_vec_)
    res += card;
  return res / card_vec_.size();
}

double Impr::GetSelectivity() {
  return 1.0;
}

