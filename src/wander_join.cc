#include <cassert>
#include <random>
#include "../include/wander_join.h"

void WanderJoin::PrepareSummaryStructure(DataGraph& g, double ratio) {

}

void WanderJoin::WriteSummary(const char* fn) {

}

void WanderJoin::ReadSummary(const char* fn) {

}

void WanderJoin::Init() {
    plans_generated_ = plan_chosen_ = false;
	walk_size_ = sample_size_ = sample_cnt_ = pos_ = inv_prob_ = 0;

    visited_.clear();
    join_from_.clear();
    join_to_.clear();
    node_to_v_.clear();
    counterparts_.clear();
    plans_.clear();
    counterpart_.clear();
    plan_.clear();
    success_cnt_.clear();
    est_.clear();
    num_idx_lookup_.clear();
    sampled_tuples_.clear();
    
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
    sample_size_ = sum / (offset_ + e_cnt); 
    sample_size_ *= sample_ratio;
    walk_size_ = offset_ + e_cnt;
    
    //set join from and to 
    for (int i = 0; i < walk_size_; i++) {
        for (int j = i + 1; j < walk_size_; j++) {
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
                    join_from_.push_back(make_pair(i, c1[m]));
                    join_to_.push_back(make_pair(j, c2[m]));
                    join_from_.push_back(make_pair(j, c2[m]));
                    join_to_.push_back(make_pair(i, c1[m]));
                }
            }
        }
    }

}

int WanderJoin::DecomposeQuery() {
    return 1;
}

//top function
void WanderJoin::generateWalkPlans() {
    for (int s = 0; s < walk_size_; s++) {
        visited_.clear();
        visited_.resize(walk_size_, false);
        visited_[s] = true;
        generateWalkPlans(0);
    }
}

//generateWalkPlanssive function
void WanderJoin::generateWalkPlans(int cnt) {
    if (sample_cnt_ == 0) 
        return;
    if (cnt == walk_size_ - 1) {
        sample_cnt_--;
        assert(plan_.size() == walk_size_ - 1);
        assert(counterpart_.size() == walk_size_ - 1);
        plans_.push_back(plan_);
        counterparts_.push_back(counterpart_);
        return;
    }
    for (int i = 0; i < join_from_.size(); i++) {
        if (visited_[join_from_[i].first] && !visited_[join_to_[i].first]) {
            visited_[join_to_[i].first] = true;
            counterpart_.push_back(join_from_[i]);
            plan_.push_back(join_to_[i]);
            generateWalkPlans(cnt + 1);
            visited_[join_to_[i].first] = false;
            counterpart_.pop_back();
            plan_.pop_back();
        }
    }
}

bool WanderJoin::checkBoundedVertices(int node, vector<int> t) {
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

bool WanderJoin::checkLabelStatistics(int node) { 
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

bool WanderJoin::checkNonTreeEdges(vector<pair<int, int>>& plan) {
	for (int i = 0; i < join_from_.size(); i++) {
		int n1 = join_from_[i].first;
		int c1 = join_from_[i].second;
		int n2 = join_to_[i].first;
		int c2 = join_to_[i].second;
		int pos1 = 0, pos2 = 0;
		for (int j = 0; j < plan.size(); j++) {
			if (plan[j].first == n1)
				pos1 = j + 1;
			if (plan[j].first == n2)
				pos2 = j + 1;
		}
		if (sampled_tuples_[pos1][c1] != sampled_tuples_[pos2][c2])
			return false;
	}
	return true;
}

//sample from the first node in the walk plan
//returns inverse probability
int WanderJoin::sampleTuple(int node) { 
	int ret;
	if (node >= offset_) {
		auto e = q->GetEdge(node - offset_);
		auto t = g->GetRandomEdge(e.el);
        sampled_tuples_.push_back(t);
        ret = g->GetNumEdges(e.el);
    } else {
        int u = node_to_v_[node];
        int vl = q->GetVLabel(u);
        auto t = g->GetRandomVertex(vl); 
        sampled_tuples_.push_back(t);
        ret = g->GetNumVertices(vl);
	}
	return ret;
}

//sample from subsequent nodes, using graph adj. list
//returns inverse probability
int WanderJoin::sampleTuple(int node, int v, int c) { 
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

//generates all walk orders
//perform random walks
//returns si and P(si)
bool WanderJoin::GetSubstructure(int subquery_index) {
	if (!plans_generated_) {
		sample_cnt_ = sample_size_;
		//this consumes sample_cnt_, preventing generating many useless walk orders
		//(an optimization)
		generateWalkPlans();
        //restore sample_cnt_
		sample_cnt_ = sample_size_;
		success_cnt_.resize(plans_.size(), 0);
		est_.resize(plans_.size());
		num_idx_lookup_.resize(plans_.size());
		plans_generated_ = true;
		if (plans_.size() == 0)
			return false;
	}
	if (sample_cnt_ <= 0)
		return false;
	while (sample_cnt_) {
		sample_cnt_--;

		sampled_tuples_.clear();
		valid_ = true;
		inv_prob_ = 1.0;
        int lookup = 0;

		int start_node = 0;
		if (counterparts_[pos_].size() > 0)
			start_node  = counterparts_[pos_][0].first;
		assert(start_node < walk_size_);
		if (!checkLabelStatistics(start_node))
			return false;
		//randomly sample an edge/vertex with edge/vertex label of start_node
		inv_prob_ *= sampleTuple(start_node);
        lookup++;
		if (!checkBoundedVertices(start_node, sampled_tuples_[0])) {
			if (!plan_chosen_) {
                est_[pos_].push_back(0);
                num_idx_lookup_[pos_].push_back(lookup);
				pos_ = (pos_ + 1) % plans_.size();
            }
            valid_ = false;
            return true;
		}

		for (int cur_order = 0; cur_order < plans_[pos_].size(); cur_order++) {
			int cur_node = plans_[pos_][cur_order].first;
			int prev_order = 0;
			for (int k = 0; k < cur_order; k++) 
				if (plans_[pos_][k].first == counterparts_[pos_][cur_order].first) 
					prev_order = k + 1;
			int c = counterparts_[pos_][cur_order].second;
			int v = sampled_tuples_[prev_order][c];
			inv_prob_ *= sampleTuple(cur_node, v, plans_[pos_][cur_order].second); 
            lookup++;
			if (inv_prob_ == 0) {
				valid_ = false;
				break;
			}
			if (!checkBoundedVertices(cur_node, sampled_tuples_.back())) {
				valid_ = false;
				break;
			}
		}
		if (!valid_ || !checkNonTreeEdges(plans_[pos_])) {
			if (!plan_chosen_) {
                est_[pos_].push_back(0);
                num_idx_lookup_[pos_].push_back(lookup);
				pos_ = (pos_ + 1) % plans_.size();
            }
            valid_ = false;
            return true;
        }
        success_cnt_[pos_]++;
        if (!plan_chosen_) {
            est_[pos_].push_back(inv_prob_);
            num_idx_lookup_[pos_].push_back(lookup);
            if (success_cnt_[pos_] >= 100) {
                int min_pos = -1;
                double min_val = std::numeric_limits<double>::max();
                for (int i = 0; i < plans_.size(); i++) {
                    if (success_cnt_[i] >= 50) {
                        double sum = 0;
                        for (double e : est_[i])
                            sum += e;
                        double mean = sum / est_[i].size();
                        double var = 0;
                        for (double e : est_[i])
                            var += (e - mean) * (e - mean);
                        //variance of estimates 
                        var = var / (est_[i].size() - 1); 

                        sum = 0;
                        for (double l : num_idx_lookup_[i])
                            sum += l;
                        //mean of index lookups
                        mean = sum / num_idx_lookup_[i].size();
                        if (var * mean < min_val) {
                            min_pos = i;
                            min_val = var * mean;
                        }
                    }
                }
                assert(min_pos != -1);
                pos_ = min_pos;
                plan_chosen_ = true;
            }
            else 
                pos_ = (pos_ + 1) % plans_.size();
        }
		return true;
	}
	//used all sample_cnt_ before choosing an order
	if (sample_cnt_ == 0)
		return false;
}

//HT estimator,
//returns 1/P(si) if valid, 0 otherwise
double WanderJoin::EstCard(int subquery_index) {
    if (valid_) {
        return inv_prob_;
    }
    else
        return 0.0;
}

double WanderJoin::AggCard() {
    double res = 0.0;
    if (card_vec_.size() == 0)
        return res;
    for (double card : card_vec_) res += card;
    return res / card_vec_.size();
}

double WanderJoin::GetSelectivity() {
    return 1;
}
