#ifndef JSUB_H_ 
#define JSUB_H_

#include "estimator.h"

class JSUB : public Estimator {
public:
	//build mode
	void PrepareSummaryStructure(DataGraph&, double); 
	void WriteSummary(const char*); 
	
	//query mode
	void ReadSummary(const char*); 
	void Init();
	int  DecomposeQuery(); 
	bool GetSubstructure(int); 
	double EstCard(int); 
	double AggCard();
	double GetSelectivity();
  
private:
	void generateWalkPlans();
	void generateWalkPlans(int);
	bool checkBoundedVertices(int, vector<int>);
	bool checkLabelStatistics(int); 
	int  sampleTuple(int);
	int  sampleTuple(int, int, int); 
	double memoi(int, pair<int, int>);
	int  nodeToOffset(int);
	int  M(int);
	
	int getR1TupleNum(int);
	int getNextTuple(int);

	//node is similarly as in WanderJoin, i.e., a table in relational model
	map<int, int> node_to_v_;
	//use adj_ instead of join_from_ and join_to_
	vector<vector<tuple<int, int, int>>> adj_; //node id -> vector of (column, joinable node id, joinable node column)

	vector<bool> visited_;
	vector<int> seq_;
	vector<pair<int, int>> plan_, counterpart_;
	vector<vector<pair<int, int>>> plans_, counterparts_;
	vector<vector<int>> sampled_tuples_;

	int offset_;
	int node_num_;
	int sample_cnt_, sample_size_;
	int pos_;
	double inv_prob_;
	bool out_of_cnt_;
	vector<map<pair<int, int>, double>> w_; //order -> tuple -> DP result

	pair<int, int> r1_tuple_;
	int r1_tuple_idx_, r1_tuple_num_;
};

#endif
