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

private:
    bool getSubstructureFlag;
    map<string, map<string, long>> mt_;
    map<string, set<string>> ceg; // small subquery -> its next level nodes
    set<string> largestMTEntries; // longest MT entries

    void decompose(const string &vListString, int mtLen);
    void getDecom(const vector<string> &vListEdges, const int &mtLen, int depth, const string &current, const string &parent);

    set<pair<string, string>> getExtensions(const string &currentVList, const string &nextVList);
    double getMaxExt(const set<pair<string, string>> &extensions, const string &queryVList, const string &queryLabelSeq);
};

#endif
