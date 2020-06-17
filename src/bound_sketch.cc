#include "../include/bound_sketch.h"
#include "../include/util.h"
#include <math.h>
#include <omp.h>
#include <experimental/filesystem>
#include <chrono>
#include <mutex>

BoundSketch::BoundSketch() {
    sketch_map_.clear();
    offline_skethces_.clear();
}

void BoundSketch::PrepareSummaryStructure(DataGraph& g, double ratio) {
#ifndef ONLINE
    buckets_ = ratio;
    assert(buckets_ >= 1);

    offline_skethces_.clear();
#ifdef PARALLEL_BUILD
    int num = g.table_.size();
    omp_lock_t lock;
    omp_init_lock(&lock);
    omp_set_num_threads(16);
#pragma omp parallel for num_threads(16) 
    for (int t = 0; t < num; t++) 
    {
        OfflineSketch *s = new OfflineSketch(t, buckets_, &g);
        omp_set_lock(&lock);
        offline_skethces_.push_back(s);
        omp_unset_lock(&lock);
    }
#else
    for (int t = 0; t < g.table_.size(); t++) {
        OfflineSketch *s = new OfflineSketch(t, buckets_, &g);
        offline_skethces_.push_back(s);
    }
#endif
#endif
}

void BoundSketch::WriteSummary(const char* fn) {
#ifndef ONLINE
    namespace fs = std::experimental::filesystem;
	fs::create_directory(fn);
    int num = offline_skethces_.size();
#ifdef PARALLEL_BUILD
#pragma omp parallel for num_threads(16) 
#endif
    for (int i = 0; i < num; i++) 
    {
        offline_skethces_[i]->serialize(fn);
    }
#endif
}

void BoundSketch::ReadSummary(const char* fn) {
    sketch_map_.clear();
#ifndef ONLINE
	OfflineSketch::deserialize(fn, sketch_map_, g);
    for (OfflineSketch* s : offline_skethces_)
        delete s;
    offline_skethces_.clear();
#endif
}

void BoundSketch::Init() {
    sketch_build_time_ = 0;
	buckets_ = sample_ratio;
    assert(buckets_ >= 1);
	bf_index_ = -1;
}

int BoundSketch::DecomposeQuery() {
	return 1;
}

//generates all bounding formulae as an intialization
//and returns a bounding formula (with index bf_index_) for each call
bool BoundSketch::GetSubstructure(int subquery_index) {
	if (bf_index_ == -1) {
		getJoinAttributeCovers();
        getBoundFormulae();
	}
    bf_index_++;
    if (bf_index_ == bound_formulae_.size())
        return false;
    return true;
}

void BoundSketch::getJoinAttributeCovers() {
    covers_.clear();
    covers_.resize(q->num_attrs());
    has_join_attribute_ = false;
    join_attribute_cnt_.clear();
    int mx = -1;
    for (int i = 0; i < q->num_relations(); i++) {
        auto& r = q->relations_[i];
        int cnt = 0;
        for (auto& a : r.attrs) {
            if (a.ref_cnt > 1) {
                covers_[a.id].push_back(i);
                cnt++;
            }
        }
        join_attribute_cnt_.push_back(cnt);
        has_join_attribute_ = true;
    }

    assignments_.clear();
    assignments_.resize(covers_.size());
    join_attribute_covers_.clear();
    getJoinAttributeCovers(0);

    rel_to_covered_attributes_.clear();
    for (vector<int>& join_cover : join_attribute_covers_) {
        bool safe = true;
        unordered_map<int, vector<int>> join_var_cover_map;
        for (int r = 0; r < q->num_relations(); r++) {
            vector<int> responsibilities;
            for (int a = 0; a < covers_.size(); a++) {
                if (join_cover[a] == r)
                    responsibilities.push_back(a);
            }
            join_var_cover_map[r] = responsibilities;
        }
        if (safe)
            rel_to_covered_attributes_.push_back(join_var_cover_map);
    }
}

void BoundSketch::getJoinAttributeCovers(int pos) {
    if (pos == covers_.size()) {
        bool safe = true;
        for (int i = 0; i < q->num_relations(); i++) {
            int cnt = 0;
            for (int j : assignments_) {
                if (i == j) 
                    cnt++;
            }
            if ((cnt > 0) && (cnt != join_attribute_cnt_[i]) && (cnt != join_attribute_cnt_[i] - 1)) {
                safe = false;
                break;
            }
        }
        if (safe) 
            join_attribute_covers_.push_back(assignments_);
        return;
    }
    //pos is not a join attribute
    if (covers_[pos].size() == 0) {
        assignments_[pos] = -1;
        getJoinAttributeCovers(pos + 1);
    } else {
        for (int r : covers_[pos]) {
            assignments_[pos] = r;
            getJoinAttributeCovers(pos + 1);
        }
    }
}

//generate all bounding formulae
void BoundSketch::getBoundFormulae() {
    bound_formulae_.clear();

    int curr = 0;

    for (auto& map : rel_to_covered_attributes_) {
        vector<Sketch*> uncL;
        vector<Sketch*> conL;
        vector<int> activeL;

        /* generate hash sizes for each attribute */
        unordered_map<int, bool> partitioned;
        vector<int> unconditionals;
        for (int r = 0; r < q->num_relations(); r++) {
            if (map[r].size() == join_attribute_cnt_[r]) 
                unconditionals.push_back(r);
        }

        bool covered_by_unc, covered_by_con;
        int num_partitioned = 0;
        for (int i = 0; i < covers_.size(); i++) {
            //i is not a join attribute
            if (covers_[i].size() == 0) {
                continue;
            }
            covered_by_unc = covered_by_con = false;
            for (int r : covers_[i]) {
                if (find(unconditionals.begin(), unconditionals.end(), r) != unconditionals.end()) {
                    covered_by_unc = true;
                } else {
                    covered_by_con = true;
                }
            }
            if (covered_by_unc && covered_by_con) {
                //must be a join attribute
                partitioned[i] = true;
                num_partitioned++;
            } else {
                partitioned[i] = false;
            }
        }

        unordered_map<int, int> hash_sizes_map;
        for (int i = 0; i < covers_.size(); i++) {
            hash_sizes_map[i] = 1;
            //i is not a join attribute
            if (covers_[i].size() == 0) {
                continue;
            }
            if (partitioned[i]) {
                hash_sizes_map[i] = round(pow(buckets_, 1.0 / num_partitioned));
                hash_sizes_map[i] = std::max(hash_sizes_map[i], 1);
            }
        }

        /* generate sketches for each relation in the join */ 
        for (int r = 0; r < q->num_relations(); r++) {
            auto& rel = q->relations_[r];
            Sketch* s = NULL;

            num_partitioned = 0;

            /* generate array of columns and array of actual attributes */
            vector<int> join_attrs_specific;
            vector<int> join_cols;
			vector<int> hash_sizes;
            for (auto& a : rel.attrs) {
                if (a.ref_cnt > 1) {
                    if (hash_sizes_map[a.id] > 1) {
                        join_attrs_specific.push_back(a.id);
                        join_cols.push_back(a.pos);
                        hash_sizes.push_back(hash_sizes_map[a.id]);
                        num_partitioned++;
                    }
                }
            }

            /* get the active column (if none (unconditional) this is null) */
            int active_attribute = -1;
            int active_col = -1;
            if (join_attribute_cnt_[r] != map[r].size()) {
                for (auto& a : rel.attrs) {
                    if (a.ref_cnt > 1) {
                        if (find(map[r].begin(), map[r].end(), a.id) == map[r].end()) {
                            active_attribute = a.id;
                            active_col = a.pos; //0 for src, 1 for dst
                            break;
                        }
                    }
                }
            }

            //for example, a sketch with hash sizes 1024,4 with join column src,dst 
            //is equivalent to a sketch with hash sizes 4,1024 with join columns dst,src
            if (join_attrs_specific.size() == 2) {
                if (join_cols[0] != 0) {
                    join_cols[0] = 0;
                    join_cols[1] = 1;
                    int temp = join_attrs_specific[0];
                    join_attrs_specific[0] = join_attrs_specific[1]; 
                    join_attrs_specific[1] = temp;
                    temp = hash_sizes[0];
                    hash_sizes[0] = hash_sizes[1];
                    hash_sizes[1] = temp; 
                }
            }
            
            if (hash_sizes.size() == 0)
                hash_sizes.push_back(1);
            vector<int> bounds;
            vector<int> bound_cols;

            /* build a probe to see if we already saw this particular sketch from another query */
            string probe;
            int alias = g->get_table_id(rel.id); 
            probe.append(to_string(alias));
            probe.append("[");
			//this wastes computation for TwoDimensionalSketchCon!
            probe.append(to_string(active_col));
            probe.append("][");
            for (int h : hash_sizes) {
                probe.append(to_string(h));
                probe.append(", ");
            }
            probe.append("][");
            for (int c : join_cols) {
                probe.append(to_string(c));
                probe.append(", ");
            }
            probe.append("][");
            for (auto& a : rel.attrs) {
#ifdef ONLINE
                if (a.is_bound) {
                    probe.append(to_string(a.pos)); 
                    probe.append(":"); 
                    probe.append(to_string(a.bound)); 
                    probe.append(", ");
                    bounds.push_back(a.bound);
                    bound_cols.push_back(a.pos);
                }
#endif
            }
            probe.append("]");
            if (sketch_map_.find(probe) != sketch_map_.end()) {
                s = sketch_map_[probe];
            } else {
#ifndef ONLINE
                assert(false); //should have been preprocessed offline! 
#endif
                auto ckpt = chrono::high_resolution_clock::now();
                //# of partitioned attributes
                if (join_attrs_specific.size() == 0) {
                    if (active_col == -1)
                        s = new ZeroDimensionalSketchUnc(alias, active_col, join_cols, hash_sizes, bounds, bound_cols, g, ""); 
                    else
                        s = new ZeroDimensionalSketchCon(alias, active_col, join_cols, hash_sizes, bounds, bound_cols, g, ""); 
                } else if (join_attrs_specific.size() == 1) {
                    if (active_col == -1)
                        s = new OneDimensionalSketchUnc(alias, active_col, join_cols, hash_sizes, bounds, bound_cols, g, ""); 
                    else 
                        s = new OneDimensionalSketchCon(alias, active_col, join_cols, hash_sizes, bounds, bound_cols, g, ""); 

                } else if (join_attrs_specific.size() == 2) {
                    if (active_attribute == -1) 
                        s = new TwoDimensionalSketchUnc(alias, active_col, join_cols, hash_sizes, bounds, bound_cols, g, "");
                    else 
                        s = new TwoDimensionalSketchCon(alias, active_col, join_cols, hash_sizes, bounds, bound_cols, g, "");
                }
                else {
                    cerr << "you're asking for too many attributes..." << endl;
                    exit(-1);
                }
                auto elapsed_nano = chrono::high_resolution_clock::now() - ckpt;
                sketch_build_time_ += (double) elapsed_nano.count() / 1000000; //in milliseconds
                sketch_map_[probe] = s;
            }

            for (int i = 0; i < join_attrs_specific.size(); i++) {
                s->l2gIndex[curr][i] = join_attrs_specific[i];
                s->g2lIndex[curr][join_attrs_specific[i]] = i;
            }

            /* make sure the active attribute appears in the g2l and l2g maps */
            if (join_attrs_specific.size() == 0 && active_attribute != -1) {
                s->l2gIndex[curr][0] = active_attribute;
                s->g2lIndex[curr][active_attribute] = 0;
            }

            if (map[r].size() == join_attribute_cnt_[r]) {
                uncL.push_back(s);
            } else if (map[r].size() == join_attribute_cnt_[r] - 1) {
                for (auto& a : rel.attrs){
                    if (a.ref_cnt > 1) {
                        if (find(map[r].begin(), map[r].end(), a.id) == map[r].end()) {
                            activeL.push_back(a.id);
                            break;
                        }
                    }
                }
                conL.push_back(s);
            } else {
                assert(map[r].size() == 0);
            }
        }

        vector<int> hash_sizes_global(covers_.size(), 1);
        for (size_t i = 0; i < covers_.size(); i++) {
            hash_sizes_global[i] = hash_sizes_map[i];
            assert(hash_sizes_global[i] > 0);
        }

        assert(conL.size() == activeL.size());
        
        if (uncL.size() > 0) {
            BoundFormula bf(curr, uncL, conL, activeL, hash_sizes_global);
            curr++;

            bound_formulae_.push_back(bf);
        }
    }
}

//returns the summation of the selected bounding formula (with index bf_index_) 
//instantiated with counts and maximum degrees of partitions 
double BoundSketch::EstCard(int subquery_index) {
    BoundFormula& bf = bound_formulae_[bf_index_];
    if (!has_join_attribute_) {
        ZeroDimensionalSketchUnc* u = (ZeroDimensionalSketchUnc*) bf.uncList[0];
        bf_index_ = bound_formulae_.size(); //no next GetSubstructure
        return u->unc[0];
    }
    else {
        long res = 0;
        CrossProductIterator cp(bf.hash_sizes);
        if (cp.totalBuckets > 1) {
            while (cp.hasNext()) {
                res += bf.execute(cp.next());
            }
        } else {
            vector<int> index(bf.hash_sizes.size(), 0);
            res = bf.execute(index);
        }
        if (res < 0)
            res = numeric_limits<long>::max();
        return (double) res;
    }
}

//min
double BoundSketch::AggCard() {
    if (card_vec_.size() == 0)
        return 0.0;
    double res = std::numeric_limits<double>::max();
    for (double card : card_vec_)
        if (card < res)
            res = card;
    return res;
}

double BoundSketch::GetSelectivity() {
#ifdef ONLINE
    cout << "online sketch build time: " << sketch_build_time_ << endl;
#endif
    return 1;
}

BoundSketch::~BoundSketch() {
    for (auto& p : sketch_map_)
        delete p.second;
}
