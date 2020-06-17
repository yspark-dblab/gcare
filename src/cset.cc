#include "../include/cset.h"
#include <boost/functional/hash.hpp>

void CharacteristicSets::PrepareSummaryStructure(DataGraph& g, double ratio) {
    const int MAX = 1e8;
    map<size_t, int> idx; //hash value of a characterisic set (CSet) to an index to csets_
    
    //build characterisic sets with forward stars
    csets_.clear();
    csets_.resize(g.GetNumVertices(), CSet()); 
    int ccnt = 0;
    for (int srcid = 0; srcid < g.GetNumVertices(); srcid++) {
        vector<int> freq(g.GetNumVLabels(srcid) + g.GetNumELabels(srcid, true), 0);
        size_t hv = 0;
        auto r = g.GetVLabels(srcid);
        int i = 0;
        for (; r.begin != r.end; r.begin++) {
            int vl = *r.begin;
            boost::hash_combine(hv, vl);
            freq[i]++;
            i++;
        }
        r = g.GetELabels(srcid, true);
        i = 0;
        for (; r.begin != r.end; r.begin++) {
            int el = *r.begin;
            boost::hash_combine(hv, el + g.GetNumVLabels());
            freq[g.GetNumVLabels(srcid) + i] += g.GetAdjSize(srcid, el, true);
            i++;
        }

        if (idx.find(hv) == idx.end()) { 
            assert(csets_[ccnt].count_ == 0);
            idx[hv] = ccnt++;
        }

        i = idx[hv];
        csets_[i].count_++;
        csets_[i].vid_ = srcid;
        if (csets_[i].count_ == 1) {
            csets_[i].freq_ = freq;
        } else {
            for (size_t j = 0; j < freq.size(); j++)
                csets_[i].freq_[j] += freq[j];
        }
    }
    csets_.erase(csets_.begin() + ccnt, csets_.end());

    //build characterisic sets with backward stars
    idx.clear();
    rev_csets_.clear();
    rev_csets_.resize(g.GetNumVertices(), CSet());
    ccnt = 0;
    for (int dstid = 0; dstid < g.GetNumVertices(); dstid++) {
        vector<int> freq(g.GetNumELabels(dstid, false), 0);
        size_t hv = 0;
        auto r = g.GetELabels(dstid, false);
        int i = 0;
        for (; r.begin != r.end; r.begin++) {
            int el = *r.begin;
            boost::hash_combine(hv, el);
            freq[i] += g.GetAdjSize(dstid, el, false);
            i++;
        }

        if (idx.find(hv) == idx.end()) { 
            assert(rev_csets_[ccnt].count_ == 0);
            idx[hv] = ccnt++;
        }

        i = idx[hv];
        rev_csets_[i].count_++;
        rev_csets_[i].vid_ = dstid;
        if (rev_csets_[i].count_ == 1) {
            rev_csets_[i].freq_ = freq;
        } else {
            for (size_t j = 0; j < freq.size(); j++)
                rev_csets_[i].freq_[j] += freq[j];
        }
    }
    rev_csets_.erase(rev_csets_.begin() + ccnt, rev_csets_.end());

    //build histograms for basic join selectivity estimation
    num_buckets_ = std::min(g.GetNumVertices() + 1, (int)(csets_.size() + rev_csets_.size()));
    if (num_buckets_ > MAX / (g.GetNumVLabels() + g.GetNumELabels()))
        num_buckets_ = MAX / (g.GetNumVLabels() + g.GetNumELabels());
    bucket_size_ = (g.GetNumVertices() + g.GetNumVLabels()) / num_buckets_;
    bucket_size_ = std::max(1, bucket_size_);

    hist_.clear();
    hist_.resize(g.GetNumVLabels() + g.GetNumELabels(), vector<vector<int>>(2, vector<int>(num_buckets_, 0)));
    for (int srcid = 0; srcid < g.GetNumVertices(); srcid++) {
        auto r = g.GetVLabels(srcid);
        for (; r.begin != r.end; r.begin++) {
            int vl = *r.begin;
            hist_[vl][0][(srcid + g.GetNumVLabels()) % num_buckets_]++;
            hist_[vl][1][vl % num_buckets_]++;
        }
        r = g.GetELabels(srcid, true);
        for (; r.begin != r.end; r.begin++) {
            int el = *r.begin;
            auto adj = g.GetAdj(srcid, el, true);
            for (; adj.begin != adj.end; adj.begin++) {
                int dstid = *adj.begin; 
                hist_[g.GetNumVLabels() + el][0][(srcid + g.GetNumVLabels()) % num_buckets_]++;
                hist_[g.GetNumVLabels() + el][1][(dstid + g.GetNumVLabels()) % num_buckets_]++;
            }
        }
    }
}

void CharacteristicSets::WriteSummary(const char* fn) {
    FILE* fp = fopen(fn, "w");
    fprintf(fp, "%zu\n", csets_.size());
    for (size_t i = 0; i < csets_.size(); i++) {
        fprintf(fp, "%d %d %zu\n", csets_[i].count_, csets_[i].vid_, csets_[i].freq_.size());
        for (int cnt : csets_[i].freq_)
            fprintf(fp, "%d ", cnt);
        fprintf(fp, "\n");
    }
    fprintf(fp, "%zu\n", rev_csets_.size());
    for (size_t i = 0; i < rev_csets_.size(); i++) {
        fprintf(fp, "%d %d %zu\n", rev_csets_[i].count_, rev_csets_[i].vid_, rev_csets_[i].freq_.size());
        for (int cnt : rev_csets_[i].freq_)
            fprintf(fp, "%d ", cnt);
        fprintf(fp, "\n");
    }
    fprintf(fp, "%d %d %zu\n", num_buckets_, bucket_size_, hist_.size()); 

    string hist_fn = string(fn) + ".hist";
    FILE* hp = fopen(hist_fn.c_str(), "w");
    for (int label = 0; label < hist_.size(); label++) {
        for (int i = 0; i < 2; i++) {
            fwrite(hist_[label][i].data(), sizeof(int), num_buckets_, hp);
        }
    }
    fclose(hp);
}

void CharacteristicSets::ReadSummary(const char* fn) {
    FILE* fp = fopen(fn, "r");
    int csize;
    fscanf(fp, "%d", &csize);
    while (csize--) {
        csets_.push_back(CSet());
        size_t idx = csets_.size() - 1;
        int size;
        fscanf(fp, "%d %d %d", &csets_[idx].count_, &csets_[idx].vid_, &size);
        int tmp;
        while (size--) {
            fscanf(fp, "%d", &tmp);
            csets_[idx].freq_.push_back(tmp);
        }
    }
    fscanf(fp, "%d", &csize);
    while (csize--) {
        rev_csets_.push_back(CSet());
        size_t idx = rev_csets_.size() - 1;
        int size;
        fscanf(fp, "%d %d %d", &rev_csets_[idx].count_, &rev_csets_[idx].vid_, &size);
        int tmp;
        while (size--) {
            fscanf(fp, "%d", &tmp);
            rev_csets_[idx].freq_.push_back(tmp);
        }
    }
    int num_hist;
    fscanf(fp, "%d %d %d", &num_buckets_, &bucket_size_, &num_hist);
    fclose(fp);

    hist_.clear();
    hist_.resize(num_hist, vector<vector<int>>(2, vector<int>(num_buckets_, 0)));
    string hist_fn = string(fn) + ".hist";
    FILE* hp = fopen(hist_fn.c_str(), "r");
    for (int label = 0; label < num_hist; label++) {
        for (int i = 0; i < 2; i++) {
            fread(hist_[label][i].data(), sizeof(int), num_buckets_, hp);
        }
    }
    fclose(hp);
}

void CharacteristicSets::Init() {
    offset_ = g->GetNumVLabels();
    covered_.clear();
    dq_.clear();
    pos_ = -1;
    nodes_.clear();
    rdf_q_adj_lists_.clear();
    rdf_q_rev_adj_lists_.clear();
    join_q_adj_mat_.clear();
}

int CharacteristicSets::DecomposeQuery() {
    //build adj. lists of RDF triples for each query vertex or vertex label
    rdf_q_adj_lists_.resize(q->GetNumVertices() + offset_);
    rdf_q_rev_adj_lists_.resize(q->GetNumVertices() + offset_);

    for (int u = 0; u < q->GetNumVertices(); u++) {
        int srcid = u + offset_; 
        int vl = q->GetVLabel(u);
        if (vl != -1) {
            nodes_.emplace_back(srcid, vl, vl);
            rdf_q_adj_lists_[srcid].emplace_back(vl, vl, nodes_.size() - 1);
            rdf_q_rev_adj_lists_[vl].emplace_back(srcid, vl, nodes_.size() - 1);
        }
        auto& r = q->GetAdj(u, true);
        for (auto& p : r) {
            int dstid = p.first + offset_;
            int el = p.second + offset_;
            nodes_.emplace_back(srcid, dstid, el);
            rdf_q_adj_lists_[srcid].emplace_back(dstid, el, nodes_.size() - 1);
            rdf_q_rev_adj_lists_[dstid].emplace_back(srcid, el, nodes_.size() - 1);
        }
    }

    join_q_adj_mat_.resize(nodes_.size(), vector<bool>(nodes_.size(), false));
    for (int srcid = 0; srcid < rdf_q_adj_lists_.size(); srcid++) {
        vector<int> lst;
        for (int i = 0; i < rdf_q_adj_lists_[srcid].size(); i++) {
            int nodeid = rdf_q_adj_lists_[srcid][i].third;
            lst.push_back(nodeid);
        }
        for (int i = 0; i < rdf_q_rev_adj_lists_[srcid].size(); i++) {
            int nodeid = rdf_q_rev_adj_lists_[srcid][i].third;
            lst.push_back(nodeid);
        }
        for (int i = 0; i < lst.size(); i++) {
            for (int j = i + 1; j < lst.size(); j++) {
                join_q_adj_mat_[lst[i]][lst[j]] = join_q_adj_mat_[lst[j]][lst[i]] = true;
            }
        }
    }
    //decompose query, extracting forward stars first
    covered_.resize(nodes_.size());

    while (true) {
        int mx = -1, max_u = -1;
        for (int srcid = offset_; srcid < rdf_q_adj_lists_.size(); srcid++) {
            int cnt = 0;
            for (int i = 0; i < rdf_q_adj_lists_[srcid].size(); i++) {
                int dstid = rdf_q_adj_lists_[srcid][i].first;
                int nodeid = rdf_q_adj_lists_[srcid][i].third;
                if (covered_[nodeid]) 
                    continue;
                cnt++;
            }
            if (mx < cnt) {
                mx = cnt;
                max_u = srcid;
            }
        }
        if (mx < 2) 
            break;
        vector<int> lst;
        for (int i = 0; i < rdf_q_adj_lists_[max_u].size(); i++) {
            int dstid = rdf_q_adj_lists_[max_u][i].first;
            int nodeid = rdf_q_adj_lists_[max_u][i].third;
            covered_[nodeid] = true;
            lst.push_back(nodeid);
        }
        for (int i = 0; i < lst.size(); i++) {
            for (int j = i + 1; j < lst.size(); j++) {
                join_q_adj_mat_[lst[i]][lst[j]] = join_q_adj_mat_[lst[j]][lst[i]] = false;
            }
        }
        dq_.push_back(make_pair(max_u, true));
    }
    //then, extract backward stars
    while (true) {
        int mx = -1, max_u = -1;
        for (int srcid = offset_; srcid < rdf_q_rev_adj_lists_.size(); srcid++) {
            int cnt = 0;
            for (int i = 0; i < rdf_q_rev_adj_lists_[srcid].size(); i++) {
                int dstid = rdf_q_rev_adj_lists_[srcid][i].first;
                int nodeid = rdf_q_rev_adj_lists_[srcid][i].third;
                if (covered_[nodeid]) 
                    continue;
                cnt++;
            }
            if (mx < cnt) {
                mx = cnt;
                max_u = srcid;
            }
        }
        if (mx < 2) 
            break;
        vector<int> lst;
        for (int i = 0; i < rdf_q_rev_adj_lists_[max_u].size(); i++) {
            int dstid = rdf_q_rev_adj_lists_[max_u][i].first;
            int nodeid = rdf_q_rev_adj_lists_[max_u][i].third;
            covered_[nodeid] = true;
            lst.push_back(nodeid);
        }
        for (int i = 0; i < lst.size(); i++) {
            for (int j = i + 1; j < lst.size(); j++) {
                join_q_adj_mat_[lst[i]][lst[j]] = join_q_adj_mat_[lst[j]][lst[i]] = false;
            }
        }
        dq_.push_back(make_pair(max_u, false));
    }
    for (int i = 0; i < covered_.size(); i++) {
        if (!covered_[i]) {
            dq_.push_back(make_pair(-i, true));
        }
    }
    return dq_.size();
}

bool CharacteristicSets::GetSubstructure(int subquery_index) {
    int v = dq_[subquery_index].first;
    //an edge between unlabeled vertices
    if (v < 0) {
        if (pos_ == -2) {
            pos_ = -1;
            return false;
        }
        pos_ = -2; 
        return true;
    }
    //forward star
    if (dq_[subquery_index].second) {
        for (int i = pos_ + 1; i < csets_.size(); i++) {
            bool flag = true;
            for (int j = 0; j < rdf_q_adj_lists_[v].size(); j++) {
                int pi = rdf_q_adj_lists_[v][j].second - offset_;
                int srcid = csets_[i].vid_;
                if (pi >= 0) {
                    if (!g->HasELabel(srcid, pi, true)) {
                        flag = false;
                        break;
                    }
                } else {
                    if (!g->HasVLabel(srcid, pi + offset_)) {
                        flag = false;
                        break;
                    }
                }
            }
            if (!flag)
                continue;
            pos_ = i;
            return true; 
        }
    } else {
    //backward star
        for (int i = pos_ + 1; i < rev_csets_.size(); i++) {
            bool flag = true;
            for (int j = 0; j < rdf_q_rev_adj_lists_[v].size(); j++) {
                int pi = rdf_q_rev_adj_lists_[v][j].second - offset_; 
                assert(pi >= 0);
                int dstid = rev_csets_[i].vid_;
                int begin, end;
                if (!g->HasELabel(dstid, pi, false)) {
                    flag = false;
                    break;
                }
            }
            if (!flag) 
                continue;
            pos_ = i;
            return true; 
        }
    }
    pos_ = -1;
    return false;
}

//please refer to the original paper for details,
//especially about computing "o" values
double CharacteristicSets::EstCard(int subquery_index) {
    int v = dq_[subquery_index].first;
    //an edge between unlabeled vertices
    if (v < 0) {
        return getNodeSelectivity(-v);
    }
    double res = 0.0;
    //forward star
    if (dq_[subquery_index].second) {
        double m = 1.0, o = 1.0;
        for (int j = 0; j < rdf_q_adj_lists_[v].size(); j++) {
            int pi = rdf_q_adj_lists_[v][j].second - offset_; 
            if (pi < 0) {
                continue;
            }
            int oi = rdf_q_adj_lists_[v][j].first - offset_;
            assert(csets_[pos_].count_ > 0);
            int srcid = csets_[pos_].vid_;
            assert(g->GetAdjSize(srcid, pi, true) > 0);
            int i = g->GetELabelIndex(srcid, pi, true);
            assert(i >= 0);
            if (q->GetBound(oi) != -1) {
                int bound = q->GetBound(oi);
                o = std::min(o, (double)1.0 / csets_[pos_].freq_[g->GetNumVLabels(srcid) + i]); 
            } else {
                m *= (double)(csets_[pos_].freq_[g->GetNumVLabels(srcid) + i]) / csets_[pos_].count_;
            }
        }
        res = (double)(csets_[pos_].count_) * m * o;
        if (v - offset_ >= 0 && q->GetBound(v - offset_) != -1) {
            assert(csets_[pos_].count_ > 0);
            res /= csets_[pos_].count_;
        }
    } else {
    //backward star
        double m = 1.0, o = 1.0;
        for (int j = 0; j < rdf_q_rev_adj_lists_[v].size(); j++) {
            int pi = rdf_q_rev_adj_lists_[v][j].second - offset_; 
            assert(pi >= 0);
            int oi = rdf_q_rev_adj_lists_[v][j].first - offset_;
            assert(rev_csets_[pos_].count_ > 0);
            int dstid = rev_csets_[pos_].vid_;
            assert(g->GetAdjSize(dstid, pi, false) > 0);
            int i = g->GetELabelIndex(dstid, pi, false);
            assert(i >= 0);
            if (q->GetBound(oi) != -1) {
                int bound = q->GetBound(oi); 
                o = std::min(o, 1.0 / rev_csets_[pos_].freq_[i]); 
            } else {
                m *= (double)(rev_csets_[pos_].freq_[i]) / rev_csets_[pos_].count_;
            }
        }
        res = (double)(rev_csets_[pos_].count_) * m * o;
        if (v - offset_ >= 0 && q->GetBound(v - offset_) != -1) {
            assert(rev_csets_[pos_].count_ > 0);
            res /= rev_csets_[pos_].count_;
        }
    }
    return res;
}

//sum
double CharacteristicSets::AggCard() {
    double res = 0.0;
    for (double card : card_vec_) res += card;
    return res;
}

double CharacteristicSets::GetSelectivity() {
    double res = 1.0;
    for (int i = 0; i < nodes_.size(); i++) {
        for (int j = i + 1; j < nodes_.size(); j++) {
            if (join_q_adj_mat_[i][j]) {
                res *= getPairwiseSelectivity(i, j);
            }
        }
    }
    return res;
}

double CharacteristicSets::getNodeSelectivity(int nodeid) {
    int el = nodes_[nodeid].third - offset_;
    if (el < 0) return hist_[el + offset_][1][(el + offset_) % num_buckets_];
    int srcid = q->GetBound(nodes_[nodeid].first - offset_);
    int dstid = q->GetBound(nodes_[nodeid].second - offset_); 
    if (el >= g->GetNumELabels())
        return 0;
    if (srcid != -1 && dstid != -1)
        return g->HasEdge(srcid, dstid, el, true) ? 1 : 0;
    if (srcid != -1 && dstid == -1)
        return g->GetAdjSize(srcid, el, true);
    if (srcid == -1 && dstid != -1)
        return g->GetAdjSize(dstid, el, false);
    return g->GetNumEdges(el); 
}

double CharacteristicSets::getPairwiseSelectivity(int i, int j) {
    int el1 = nodes_[i].third;
    int el2 = nodes_[j].third;
    double cnt1 = getNodeSelectivity(i), cnt2 = getNodeSelectivity(j);
    if (cnt1 < 1e-9 || cnt2 < 1e-9) return 0.0;
    double sum = 0.0;
    int c1 = 1, c2 = 1;
    if (nodes_[i].first == nodes_[j].first) {
        c1 = c2 = 0;
    } else if (nodes_[i].first == nodes_[j].second) {
        c1 = 0;
        c2 = 1;
    } else if (nodes_[i].second  == nodes_[j].first) {
        c1 = 1;
        c2 = 0;
    } else if (nodes_[i].second == nodes_[j].second) {
        c1 = c2 = 1;
    }
    //use basic join selectivity estimation
    for (int b = 0; b < num_buckets_; b++) {
        assert(bucket_size_ > 0);
        sum += (double)(hist_[el1][c1][b]) * hist_[el2][c2][b] / bucket_size_;
    }
    return sum / cnt1 / cnt2;
}
