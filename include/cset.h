#ifndef CSET_H_
#define CSET_H_

#include <algorithm>
#include <vector>
#include "estimator.h" 

class CharacteristicSets : public Estimator {
public:
	//build mode
	void PrepareSummaryStructure(DataGraph&, double); 
	void WriteSummary(const char*); 
	
	//query mode
	void Init();
	void ReadSummary(const char*); 
	int  DecomposeQuery(); 
	bool GetSubstructure(int); 
	double EstCard(int);
	double AggCard();
	double GetSelectivity();
	
	struct CSet {
		int vid_;
		int count_;
		vector<int> freq_;

		CSet() : count_(0) { freq_.clear(); }
	};

private:
	double getNodeSelectivity(int);
	double getPairwiseSelectivity(int, int);

	struct Triple {
		int first;
		int second;
		int third;
		Triple(int f, int s, int t) : first(f), second(s), third(t) {}
	};
	
	//vector of (src, dst, predicate), i.e., each node is a triple in RDF query
	//here, (src, dst, predicate) = (u, v, el) for query vertices u, v and edge label el,
	//or (u, vl, vl) indicating that u has a vertex label vl;
	//we regard it has a triple from u to a virtual vertex vl with predicate label vl
	vector<Triple> nodes_;
	int offset_; //# vertex labels in data graph
	vector<bool> covered_; //whether a node is covered in DecomposeQuery
	vector<pair<int, bool>> dq_; //decomposed subqueries
	//rdf_q_adj_lists_[i] represents  i's forward adj. list
	//which is a vector of (dst, predicate, index to nodes_)
	vector<vector<Triple>> rdf_q_adj_lists_, rdf_q_rev_adj_lists_;  
	vector<vector<bool>> join_q_adj_mat_;
	
	//csets_[i]: i'th forward characteristic set, rev_csets_[i]: i'th backward
	vector<CSet> csets_, rev_csets_; 
	int pos_; //index to csets_ or rev_csets_

	int num_buckets_;
	int bucket_size_;
	vector<vector<vector<int>>> hist_; //vertex/edge label, src/dst, bucket 
};

#endif
