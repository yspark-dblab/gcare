#include <cassert>
#include "../include/jsub.h"

void JSUB::PrepareSummaryStructure(DataGraph& g, double p) {
}

void JSUB::WriteSummary(const char* fn) {
}

void JSUB::ReadSummary(const char* fn) {
}

void JSUB::Init() {
    node_num_ = sample_size_ = sample_cnt_ = 0;
    dp_time_ = num_est_card_ = num_memoi_ = r1_tuple_idx_ = 0;
    offset_ = g->GetNumELabels();

    visited_.clear();
    adj_.clear();
    node_to_v_.clear();
    counterparts_.clear();
    plans_.clear();
    counterpart_.clear();
    plan_.clear();
    w_.clear();

    //set sample size
    int sum = 0;
    offset_ = 0;
    int e_cnt = q->GetNumEdges(); 
    for (int v = 0; v < q->GetNumVertices(); v++) {
        int vl = q->GetVLabel(v); 
        if (vl != -1) {
            sum += g->GetNumVertices(vl); 
            node_to_v_[offset_] = v;
            offset_++;
        }
        for (auto& e : q->GetAdj(v, true)) {
            sum += g->GetNumEdges(e.second);
        }
    }
    sample_size_ = sum / (offset_+ e_cnt); 
    sample_size_ *= sample_ratio;
    node_num_ = offset_ + e_cnt;
    adj_.resize(node_num_);
    w_.resize(node_num_);
    
    //similar to WanderJoin
    for (int i = 0; i < node_num_; i++) {
        for (int j = i + 1; j < node_num_; j++) {
            //edges i and j can have join conditions on both query vertices, 
            //i.e., i and j makes a loop
            int c1[2] = {-1, -1}; 
            int c2[2] = {-1, -1};
            int n = 0;
            if (i < offset_) {
                int v = node_to_v_[i]; 
                int vl = q->GetVLabel(v);
                assert(vl != -1);
                if (j >= offset_) {
                    auto e2 = q->GetEdge(j - offset_);
                    if (v == e2.src) {
                        c1[n] = c2[n] = 0;
                        n++;
                    }
                    if (v == e2.dst) {
                        c1[n] = 0;
                        c2[n] = 1;
                        n++;
                    }
                } else {
                    int u = node_to_v_[j]; 
                    assert(u != v);
                }
                //do not consider j >= e_cnt since no q vertex has more than 1 label
                //so c1 and c2 remains -1
            } else {
                auto e1 = q->GetEdge(i - offset_);
                if (j >= offset_) {
                    auto e2 = q->GetEdge(j - offset_);
                    if (e1.src == e2.src) {
                        c1[n] = c2[n] = 0;
                        n++;
                    }
                    if (e1.dst == e2.src) {
                        c1[n] = 1;
                        c2[n] = 0;
                        n++;
                    }
                    if (e1.src == e2.dst) {
                        c1[n] = 0;
                        c2[n] = 1;
                        n++;
                    }
                    if (e1.dst == e2.dst) {
                        c1[n] = c2[n] = 1;
                        n++;
                    }
                } else {
                    int v = node_to_v_[j]; 
                    int vl = q->GetVLabel(v);
                    assert(vl != -1);
                    if (e1.src == v) {
                        c1[n] = c2[n] = 0;
                        n++;
                    }
                    if (e1.dst == v) {
                        c1[n] = 1;
                        c2[n] = 0;
                        n++;
                    }
                }
            }
            for (int m = 0; m < n; m++) { 
				if (c1[m] != -1 && c2[m] != -1) {
					adj_[i].push_back(make_tuple(c1[m], j, c2[m]));
					adj_[j].push_back(make_tuple(c2[m], i, c1[m]));
				}
            }
        }
    }
}

void JSUB::generateWalkPlans(int cnt) {
    if (sample_cnt_ == 0) {
        return;
    }
    if (cnt == node_num_ - 1) {
        sample_cnt_--;
        plan_.clear();
        plan_.push_back(make_pair(seq_[0], -1));
        counterpart_.clear();
        counterpart_.push_back(make_pair(-1, -1));

        for (int pos1 = 1; pos1 < node_num_; pos1++) {
            int i = seq_[pos1];
            bool acyclic = true; 
            int connected_col = -1;
            pair<int, int> prev(-1, -1);
            for (int pos2 = 0; pos2 < pos1; pos2++) {
                int j = seq_[pos2];
                for (auto k : adj_[i]) {
                    if (get<1>(k) == j) { //node i and j is connected
                        if (connected_col == -1) {
                            connected_col = get<0>(k);
                            prev = make_pair(j, get<2>(k));
                        } else if (connected_col != get<0>(k))
                            acyclic = false;
                    }
                }
            }
            if (acyclic) {
				plan_.push_back(make_pair(i, connected_col));
                counterpart_.push_back(prev);
            }
        }
        plans_.push_back(plan_);
        counterparts_.push_back(counterpart_);

		return;
    }
    for (int i = 0; i < node_num_; i++) {
        if (visited_[i])
            continue;
        for (auto n : adj_[i]) {
            int j = get<1>(n); 
            if (visited_[j]) {
                visited_[i] = true;
                seq_[cnt + 1] = i;
                generateWalkPlans(cnt + 1);
                visited_[i] = false;
                seq_[cnt + 1] = -1;
                break;
            }
        }
    }
}

void JSUB::generateWalkPlans() {
    for (int s = 0; s < node_num_; s++) {
        visited_.clear();
        seq_.clear();
        visited_.resize(node_num_, false);
        seq_.resize(node_num_, -1);
        visited_[s] = true;
        seq_[0] = s;
        generateWalkPlans(0);
    }
}

bool JSUB::checkBoundedVertices(int node, vector<int> t) {
	if (node >= offset_) {
		auto e = q->GetEdge(node - offset_);
		if (q->GetBound(e.src) >= 0 && q->GetBound(e.src) != t[0])
			return false;
		if (q->GetBound(e.dst) >= 0 && q->GetBound(e.dst) != t[1])
			return false;
	} else {
		int u = node_to_v_[node];
		if (q->GetBound(u) >= 0 && q->GetBound(u) != t[0])
			return false;
	}
	return true;
}

int JSUB::nodeToOffset(int node) {
    if (node >= offset_) { 
		auto e = q->GetEdge(node - offset_);
        return e.el;
    } else {
        int u = node_to_v_[node];
        return offset_ + q->GetVLabel(u);
    }
}

bool JSUB::checkLabelStatistics(int node) { 
	if (node >= offset_) { 
		auto e = q->GetEdge(node - offset_);
		if (g->GetNumEdges(e.el) == 0) {
			return false;
        }
	} else {
		int u = node_to_v_[node];
		int vl = q->GetVLabel(u);
		if (g->GetNumVertices(vl) == 0) {
			return false;
        }
	}
	return true;
}

//sample from the first node in the walk plan
//returns inverse probability
int JSUB::sampleTuple(int node) { 
    assert(sampled_tuples_.size() == 0);
	int ret;
	if (node >= offset_) {
		auto e = q->GetEdge(node - offset_);
		vector<int> t = g->GetRandomEdge(e.el);
		sampled_tuples_.push_back(t);
		ret = g->GetNumEdges(e.el);
	} else {
		int u = node_to_v_[node];
		int vl = q->GetVLabel(u);
        assert(vl != -1);
		vector<int> t = g->GetRandomVertex(vl); 
        assert(t.size() == 1);
        t.push_back(-1);
        assert(t.size() == 2);
		sampled_tuples_.push_back(t);
		ret = g->GetNumVertices(vl);
        assert(sampled_tuples_.size() == 1);
        assert(sampled_tuples_[0].size() == 2);
	}
	return ret;
}

int JSUB::getR1TupleNum(int node) {
    int ret;
	if (node >= offset_) {
		auto e = q->GetEdge(node - offset_);
		return g->GetNumEdges(e.el);
	} else {
		int u = node_to_v_[node];
		int vl = q->GetVLabel(u);
        assert(vl != -1);
        return g->GetNumVertices(vl);
	}
}

//get next tuple from the first node in the walk plan
//returns the size of the first node
int JSUB::getNextTuple(int node) { 
    assert(sampled_tuples_.size() == 0);
	int ret;
	if (node >= offset_) {
		auto e = q->GetEdge(node - offset_);
		vector<int> t = g->GetEdge(e.el, r1_tuple_idx_);
		sampled_tuples_.push_back(t);
		ret = g->GetNumEdges(e.el);
	} else {
		int u = node_to_v_[node];
		int vl = q->GetVLabel(u);
        assert(vl != -1);
		vector<int> t = g->GetVertex(vl, r1_tuple_idx_); 
        assert(t.size() == 1);
        t.push_back(-1);
        assert(t.size() == 2);
		sampled_tuples_.push_back(t);
		ret = g->GetNumVertices(vl);
        assert(sampled_tuples_.size() == 1);
        assert(sampled_tuples_[0].size() == 2);
	}
    r1_tuple_idx_++;
	return ret;
}

//sample from subsequent nodes, using graph adj. list
//returns inverse probability
int JSUB::sampleTuple(int node, int v, int c) { 
	int ret;
	if (node >= offset_) {
		auto e = q->GetEdge(node - offset_);
		auto t = g->GetRandomEdge(v, e.el, c == 0);
		if (t.size() > 0) {
			sampled_tuples_.push_back(t);
			ret = g->GetAdjSize(v, e.el, c == 0);
		} else
			ret = 0;
	} else {
		int u = node_to_v_[node];
		int vl = q->GetVLabel(u);
		if (g->HasVLabel(v, vl)) {
			sampled_tuples_.push_back(vector<int>(1, v));
			ret = 1;
		} else 
			ret = 0;
	}
	return ret;
}

//for each possible q1 and o, use WanderJoin to estimate |q1|
//also calculate M(q1) and choose q1, o with min. |q1| * M(q1)
int JSUB::DecomposeQuery() {
    sample_cnt_ = sample_size_;
    generateWalkPlans();
    sample_cnt_ = sample_size_;
    pos_ = -1;

    if (plans_.size() == 0) {
        return 1;
    }

    double min_bound = std::numeric_limits<double>::max();

    //i iterates from 0 to n-1, try WanderJoin only once per each plan
    for (int i = 0, count = sample_size_; i < plans_.size() && count > 0; i++, count--) {
        sampled_tuples_.clear();
        assert(sampled_tuples_.size() == 0);
		inv_prob_ = 1.0;

        int start_node = plans_[i][0].first;
        assert(start_node >= 0 && start_node < node_num_);
        if (!checkLabelStatistics(start_node))
			return 1;
		inv_prob_ *= sampleTuple(start_node);

        if (!checkBoundedVertices(start_node, sampled_tuples_[0])) {
            i = i + 1;
            continue;
		}
        
        bool valid = true;
		for (int cur_order = 1; cur_order < plans_[i].size(); cur_order++) {
			int cur_node = plans_[i][cur_order].first;
			int cur_col  = plans_[i][cur_order].second;
			int prev_order = 0;
			for (int k = 1; k < cur_order; k++) 
				if (plans_[i][k].first == counterparts_[i][cur_order].first) 
					prev_order = k;
			int c = counterparts_[i][cur_order].second;
			int v = sampled_tuples_[prev_order][c];
			inv_prob_ *= sampleTuple(cur_node, v, plans_[i][cur_order].second); 
            if (inv_prob_ == 0) {
				valid = false;
				break;
			}
            if (!checkBoundedVertices(cur_node, sampled_tuples_.back())) {
				valid = false;
				break;
			}
        }

		if (!valid) {
            i = i + 1;
            continue;
        }
        
        //if sampling succeeds with the current plan and we obtain a nonzero estimate of |q1|,
        //see if the current plan gives the minimum estimate of |q1| * M(q1)
        double cur_bound = inv_prob_ * M(i);
        if (cur_bound < min_bound) {
            min_bound = cur_bound;
            pos_ = i;
        }
    }

    if (pos_ != -1) {
        sample_cnt_ *= node_num_;
#ifdef FULL_DP
        r1_tuple_num_ = getR1TupleNum(plans_[pos_][0].first); 
#endif
    }

    return 1;
}

//get the a tuple t1 from the first relation R1 in the walk plan
bool JSUB::GetSubstructure(int subquery_index) {
    //no plan is chosen
    if (pos_ == -1) {
        return false;
    }
#ifdef FULL_DP
    if (r1_tuple_idx_ >= r1_tuple_num_)
        return false;
#else
	//out of sample size
    if (sample_cnt_ <= 0) {
        return false;
    }
#endif
    int start_node = plans_[pos_][0].first;
    assert(start_node >= 0 && start_node < node_num_);
    sampled_tuples_.clear();
#ifdef FULL_DP
    inv_prob_ = getNextTuple(start_node);
#else
    inv_prob_ = sampleTuple(start_node);
#endif
    sample_cnt_--;
    
    if (!checkBoundedVertices(start_node, sampled_tuples_[0])) {
        inv_prob_ = 0;
        return true; 
    }
        
    r1_tuple_ = make_pair(sampled_tuples_[0][0], sampled_tuples_[0][1]);
    return true;
}

//use dynamic programming to get w(t1) for the sampled t1
//multiply it with the inverse probability of sampling t1
//finally, multiply with M which is always 1 in our context
double JSUB::EstCard(int subquery_index) {
    if (inv_prob_ == 0) {
        return 0;
    }
    out_of_cnt_ = false;
    num_est_card_++;
    double join = memoi(0, r1_tuple_); 
    if (out_of_cnt_) {
        return 0;
    }
    return inv_prob_ * join * M(pos_);
}

//perform dynamic programming using the remaining sample_cnt_ 
double JSUB::memoi(int order, pair<int, int> tuple) {
#ifndef FULL_DP
    if (sample_cnt_ <= 0) {
        out_of_cnt_ = true;
        return 0;
    }
#endif
    num_memoi_++;
    if (w_[order].find(tuple) != w_[order].end()) {
        return w_[order][tuple];
    }
    double res = 1;
    for (int next = order + 1; next < plans_[pos_].size(); next++) {
        if (counterparts_[pos_][next].first != plans_[pos_][order].first) {
            continue;
        }
        double sum = 0;
        int next_node = plans_[pos_][next].first; 
        int next_col  = plans_[pos_][next].second; 
        int prev_col  = counterparts_[pos_][next].second;
        int v = prev_col == 0 ? tuple.first : tuple.second;
        int c = plans_[pos_][next].second;
        if (next_node >= offset_) {
            auto e = q->GetEdge(next_node - offset_);
            auto r = g->GetAdj(v, e.el, c == 0); 
            for (; r.begin != r.end; r.begin++) {
                if (c == 0)
                    sum += memoi(next, make_pair(v, *r.begin));
                else
                    sum += memoi(next, make_pair(*r.begin, v));
                sample_cnt_--;
            }             
        } else {
            int u = node_to_v_[next_node];
            int vl = q->GetVLabel(u);
            if (g->HasVLabel(v, vl)) {
                sum += memoi(next, make_pair(v, -1));
            }         
        }
        res *= sum;
    }
    w_[order][tuple] = res;
    return res;
}
	
//avg
double JSUB::AggCard() {
	double res = 0.0;
    if (card_vec_.size() == 0)
        return res;
    for (double card : card_vec_) res += card;
    return res / card_vec_.size();
}

double JSUB::GetSelectivity() {
    return 1;
}

//maximum number of the join results between a tuple that matches q1 and the residual relations 
int JSUB::M(int p) {
	return 1;
}
