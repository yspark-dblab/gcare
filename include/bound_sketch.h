#ifndef BOUND_SKETCH_H_ 
#define BOUND_SKETCH_H_ 

#include "../include/estimator.h"
#include "../include/bound_formula.h"

class BoundSketch : public Estimator {
public:
	BoundSketch();
	~BoundSketch();
	
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
	void getJoinAttributeCovers();
	void getJoinAttributeCovers(int);
	void getBoundFormulae();

	unordered_map<string, Sketch*> sketch_map_;
	vector<OfflineSketch*> offline_skethces_; 
	int buckets_;
	int bf_index_;
	vector<vector<int>> covers_; //query's join attribute -> vector of covering relations
	vector<int> join_attribute_cnt_; //relation -> # covering join attributes in query
	bool has_join_attribute_;
	vector<int> assignments_; //join attribute -> a covering relation in a bounding formula 
	vector<vector<int>> join_attribute_covers_; //vector of assignments
	vector<unordered_map<int, vector<int>>> rel_to_covered_attributes_; //vector of relation -> covering join attributes mappings
	vector<BoundFormula> bound_formulae_;
	double sketch_build_time_;
};

#endif

