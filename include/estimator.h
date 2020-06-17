#ifndef ESTIMATOR_H_ 
#define ESTIMATOR_H_ 

#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>


#ifdef RELATION
  #include "data_relations.h"
  #include "query_relations.h"
#else
    #include "data_graph.h"
    #include "query_graph.h"
#endif

using namespace std;

class Estimator {
public:

  enum ErrCode { NORMAL, TIMEOUT };
	unsigned long long int dp_time_;
	int num_est_card_;
	int num_memoi_;

  void Summarize(DataGraph& g, const char* fn, double p) {
		PrepareSummaryStructure(g, p);
		WriteSummary(fn);
	}

	double Run(DataGraph& gin, QueryGraph& qin, double p) {
		g = &gin;
		q = &qin;
		sample_ratio = p;
		subquery_card_.clear();
		Init();
		
		num_subqueries_ = DecomposeQuery();
		for (int j = 0; j < num_subqueries_; j++) {
			card_vec_.clear();
			while (GetSubstructure(j)) {
				double est_card = EstCard(j);
				card_vec_.push_back(est_card);
			}
			double agg_card = AggCard();
			subquery_card_.push_back(agg_card);
		}
		
		double ret = 1.0;
		for (int j = 0; j < num_subqueries_; j++) {
			ret *= subquery_card_[j];
		}
		selectivity_ = GetSelectivity();
		ret *= selectivity_;
		return ret;
	}

	//method-specific functions to be implemented
	
	//build mode
	virtual void PrepareSummaryStructure(DataGraph&, double) = 0; 
	virtual void WriteSummary(const char*) = 0; 
	
	//query mode
	virtual void ReadSummary(const char*) = 0; 
	virtual void Init() = 0;
	virtual int DecomposeQuery() = 0; 
	virtual bool GetSubstructure(int) = 0; 
	virtual double EstCard(int) = 0; 
	virtual double AggCard() = 0;
	virtual double GetSelectivity() = 0;

    void SetDataGraph(DataGraph* _g) {
        g = _g;
    }

    DataGraph* GetDataGraph() {
        return g;
    }

protected:
	DataGraph* g;
	QueryGraph *q;
	double sample_ratio;

	int num_subqueries_;
	double selectivity_;
	vector<double> subquery_card_; //for each subquery
	vector<double> card_vec_;      //for each subquery and substructure
};

#endif
