#ifndef IMPR_H_
#define IMPR_H_

#include "../include/estimator.h"
#include "../include/data_graph.h"
#include "../include/query_graph.h"

#include <algorithm>
#include <unordered_map>
#include <vector>
#include <boost/functional/hash.hpp>

using std::unordered_map;
using std::pair;
using std::vector;

class Impr : public Estimator {
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

  double Count(int, double);
  int getBeta(vector<int>, vector<bool>, int);
  off_t StrLen(char*);
  int StrToInt(char*);
  void Preprocess(void);
  void DFS(int);
  pair<int, bool> ChooseELabel(int, vector<pair<int, bool>>&, int&, int&);
  int Select(vector<int>&, vector<int>&, bool&, int&);
  bool Check(vector<int>&, vector<int>&, vector<pair<int, bool>>&);
  bool GetNextSample(vector<int>&, vector<int>&, vector<pair<int, bool>>&, int&, int);
  double GetWeight(vector<int>&);
  void FindPaths(int, vector<int>&, int, double);
  int GetBeta(int);

 private:
  struct ContainerHash {
    size_t operator () (vector<pair<int, bool>> const& c) const {
      return boost::hash_range(c.begin(), c.end());
    }
  };
  struct Vector1DHash {
    size_t operator () (vector<int> const& c) const {
      return boost::hash_range(c.begin(), c.end());
    }
  };
  struct Comp {
    bool operator()(const pair<int, int>& lhs, const int& rhs) { return lhs.first < rhs; }
    bool operator()(const int& lhs, const pair<int, int>& rhs) { return lhs < rhs.first; }
    bool operator()(const int& lhs, const int& rhs) { return lhs < rhs; }
  };
  unordered_map<vector<int>, bool, Vector1DHash> is_duplicate_;
  vector<int> qv_;
  vector<bool> chk_;
  vector<int> indeg_, outdeg_;
  int cnt1_, cnt2_;
  vector<double> est_;
  double sum_;
  int size_;
  vector<int> xk_, case_num_;
  vector<pair<int, bool>> el_;
  int beta_, step_, s_num_;
  vector<pair<int, bool>> query_labels_;
  vector<vector<int>> pos_embs_;
  vector<int> path_;
};

#endif
