#ifndef MT_H_
#define MT_H_

#include <algorithm>
#include <vector>
#include <string>
#include <set>
#include "estimator.h" 

class MarkovTable : public Estimator {
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

    double EstCardAllMax(int subquery_index);
    double EstCardGreedyMax(int subquery_index);

private:
    bool getSubstructureFlag;

    // new encoding
    vector<long> mt1_;
    vector<vector<vector<vector<long>>>> mt2_;

    void insertEntryToMT(const vector<string> &entry);
    void getExtensions(vector<tuple<int, int, Edge, Edge>> &extensions, const vector<Edge> &current, const int &currentEnc);
    double calcExtRate(const tuple<int, int, Edge, Edge> &ext);
    void makeEstAndAddToQueue(vector<SubQEdgeNode> &subQEdgeNodePool, const int &currentIdx, const double &currentEst,
                              const tuple<int, int, Edge, Edge> &ext, deque<int> &queue,
                              vector<bool> &processed, vector<double> &estimates);
};

#endif
