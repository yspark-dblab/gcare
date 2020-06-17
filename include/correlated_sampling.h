#ifndef CORRELATED_SAMPLING_H_
#define CORRELATED_SAMPLING_H_

#include "relation.h"
#include "estimator.h"
#include <boost/functional/hash.hpp>
#include <algorithm>
#include <vector>
#include <random>


class CorrelatedSampling: public Estimator {
public:
    void PrepareSummaryStructure(DataGraph&, double) {}
    void WriteSummary(const char*) {}

    void ReadSummary(const char*) {}

    void Init();
    int DecomposeQuery() { num_subqueries_ = 1; status_ = true; return 1; }
    bool GetSubstructure(int);
    double EstCard(int);
    double AggCard() { 
	return card_vec_[0];
    }
    double GetSelectivity() { return 1; }
private:
    std::vector<Relation<int>> samples_;
    std::vector<std::pair<uint64_t,uint64_t>> seeds_;
    bool status_;
    std::vector<int> pos_;
    std::vector<double> pmins_;
    int m3_;
    uint64_t EstCard_(DataGraph&, QueryGraph&);
    double Hash(std::pair<uint64_t, uint64_t> &, int);
};



#endif
