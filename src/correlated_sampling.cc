#include "../include/correlated_sampling.h"
#include <random>

const uint64_t M61 = 2305843009213693951ull;

double CorrelatedSampling::Hash(std::pair<uint64_t, uint64_t> &seed, int x) {
    uint64_t res = seed.first * x + seed.second;
    res = res % M61;
    return static_cast<double>(res) / M61;
}

double CorrelatedSampling::EstCard(int subquery_index) {
    // 1. Calculate join size
    uint64_t join_size = EstCard_(*g, *q);
    // 2. Calculate P(s_0)
    double p_s_0 = 1.0;
    for (double pmin : pmins_) {
        p_s_0 *= pmin; // pmin:  min_{A_{R_i} \ni a} p^{\frac{1}{|A_{R_i}|}}
    }
    double ans = static_cast<double>(join_size) / p_s_0;
    return ans;
}

uint64_t CorrelatedSampling::EstCard_(DataGraph& data_graph, QueryGraph& query_graph) {        
    if (samples_.size() == 0) return 0.0;
    Relation<int> &sample_l = samples_[0]; 
    if (sample_l.size() == 0) return 0.0;

    vector<bool> chk_sample(query_graph.relations_.size(), false); // if chk_sample[i] is true, then the sample i is joined to sample_l
    vector<bool> chk_attr(query_graph.num_attrs(), false); // if chk_attr[i] is true, then the sample_l contains the attribute i
    //==============================================
    // 1. Start from table 0
    // ---------------------------------------------
    chk_sample[0] = true;
    auto &relation_l = query_graph.relations_[0];
    // 1-1 make colum names for sample 0
    std::vector<int> cols_l(relation_l.size());
    for (size_t i = 0; i < cols_l.size(); ++i) {
        int cname = relation_l[i].id;
        cols_l[i] = cname;
        chk_attr[cname] = true;
    }
    // start join the remaining samples on the table 0
    std::vector<int> cols_r;
    for (int i = 0; i < query_graph.relations_.size() - 1; ++i) {
        // choose a sample using a heuristic method
        int t = 0, mx = 0;
        for (int j = 1; j < query_graph.relations_.size(); ++j) {
            if (chk_sample[j]) continue;
            if (samples_[j].size() == 0) return 0.0;
            auto &attrs = query_graph.relations_[j].attrs;
            int cnt = 0;
            for (int k = 0; k < attrs.size(); ++k) {
                auto &attr = attrs[k];
                if (chk_attr[attr.id]) cnt++;
            }
            if (mx < cnt) {
                mx = cnt;
                t = j;
            }
        }
        chk_sample[t] = true;
        auto &attrs = query_graph.relations_[t].attrs;
        cols_r.resize(attrs.size());
        Relation<int> &sample_r = samples_[t]; 
        for (size_t j = 0; j < attrs.size(); ++j) {
            int cname = attrs[j].id;
            cols_r[j] = cname; //pos;
            chk_attr[cname] = true;
        }
        sample_l.join(cols_l, cols_r, sample_r);
        if (sample_l.size() == 0) return 0.0;
    }
    return sample_l.size();
}

void CorrelatedSampling::Init() {
    samples_.clear();
    seeds_.clear();
    pmins_.clear();
    m3_ = -1;
}

bool CorrelatedSampling::GetSubstructure(int subquery_index) {
    if (!status_) return false;
    status_ = false;
    DataGraph& data = *g;
    QueryGraph& query = *q;

    int rand_val = rand();
    std::mt19937 generator_csj(rand_val);
    std::uniform_int_distribution<uint64_t> dis_csj(0, M61 - 1);

    //=============================================
    // 1. preprocess
    //---------------------------------------------
    size_t num_attrs = query.num_attrs();
    seeds_.resize(num_attrs, std::make_pair<uint64_t, uint64_t>(0, 0));
    std::vector<bool> is_join_attribute(num_attrs, false);
    pmins_.resize(num_attrs, 1.0);
    vector<double> sample_probs(query.relations_.size(), 1.0);
    for (size_t i = 0; i < query.relations_.size(); ++i) {
        auto & rel = query.relations_[i];
        // calculate p^{\frac{1}{|A_{R_i}|}}
        int num_join_attrs = 0;
        for (size_t j = 0; j < rel.attrs.size(); ++j) {
            auto &attr = rel.attrs[j];
            if (attr.is_bound) {
            } else if (attr.ref_cnt > 1) { // this is a join attribute.
                is_join_attribute[attr.id] = true;
                num_join_attrs += 1;
            }
        }
        double sample_prob = pow(sample_ratio, 1.0 / static_cast<double>(num_join_attrs));
        assert(num_join_attrs <= 2);
        assert(num_join_attrs >= 1);
        sample_probs[i] = sample_prob;
        // num_join_attrs == 2 ? sample_prob : std::sqrt(sample_prob);
        // calculate min_{A_{R_i} \ni a} p^{\frac{1}{|A_{R_i}|}}
        for (size_t j = 0; j < rel.attrs.size(); ++j) {
            auto &attr = rel.attrs[j];
            if (attr.is_bound) { // the bound attribute is not a join attribute
            } else if (attr.ref_cnt > 1) { // this is a join attribute. calculate minimum probability of it among related relations
                pmins_[attr.id] = std::min(sample_prob, pmins_[attr.id]);
            }
        }
    }
    // Prepare seed of h_a for each join attribute
    for (size_t i = 0; i < num_attrs; ++i) {
        if (!is_join_attribute[i]) continue;
        seeds_[i] = std::make_pair<uint64_t, uint64_t>(dis_csj(generator_csj), dis_csj(generator_csj));
    }
    //=============================================
    // 2. Create the sample s_0 as a list of relations <S1, ..., Sn>
    //---------------------------------------------
    samples_.resize(query.relations_.size()); // change code to resize samples_ once
    for (size_t i = 0; i < query.relations_.size(); ++i) {
        auto &rel = query.relations_[i];
        int t = data.get_table_id(rel.id); // get table id 

        auto &sample = samples_[i];
        sample.SetNumCols(rel.attrs.size());

        double sample_prob = sample_probs[i];
        vector<int> tuple(rel.attrs.size());

        // Sample the tuples for sample i
        for (size_t j = 0; j < data.table_[t].size(); j++) {
            bool pass = true;
            for (size_t k = 0; k < rel.attrs.size(); ++k) {
                auto &attr = rel.attrs[k];
                int pos = attr.pos;
                auto val = data.table_[t][j][pos];
                if (!attr.is_bound) {
                    if (attr.ref_cnt < 2) continue;
                    double h_a = Hash(seeds_[attr.id], val);
                    if (h_a >= pmins_[attr.id]) { // drop this tuple
                        pass = false;
                        break;
                    }
                    tuple[k] = val;
                } else if (attr.bound != val) { // drop this tuple
                    pass = false;
                    break;
                }
            }
            if (pass) {
                samples_[i].append(tuple);
            }
        }
    }
    //=============================================
    return true;
}
