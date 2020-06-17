#ifndef SUMRDF_H_
#define SUMRDF_H_

#include <boost/functional/hash.hpp>
#include <algorithm>
#include <vector>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include "estimator.h" 
#include "data_graph.h"
#include "util.h"
#include "graph_op.h"

class SumRDF : public Estimator {
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
  class Type {
  public:
    set<int> classes_;
    set<int> outgoing_;
    set<int> incoming_;
    Type();
    string ToString();
    string ToStringByClass();
  };
  class Resource {
   public:
    Type type_;
    int bound_;
    Resource();
  };
  class Bucket {
  public:
    Type type_;
    int type_idx_;
    vector<int> resources_;
    Bucket(Type&);
    Bucket() { };
  };
  class Summary {
   public:
    vector<Bucket> buckets_;
    vector<Edge> edges_;
    vector<vector<pair<int, int>>> adj_list_, rev_adj_list_;
    int multiplicity_, ori_multiplicity_;
    double build_time_;
    Summary();
    void CreateBucket(Resource&, int);
  };
  class Partition {
   public:
    int pid_;
    vector<int> partition_;
    vector<int> finer_partitions_;
    vector<vector<int>> blocks_;
    vector<int> gamma_;
    vector<Edge> q_mgu_;
    Partition();
    Partition(vector<int>&, int, vector<vector<int>>&);
    bool IsTauUnifiable(vector<int>&);
    bool IsFiner(Partition&);
    void ComputeQuery(QueryGraph&);
    bool IsContained(vector<int>&, vector<int>&);
  };
  Summary sm_;
  DataGraph s_;
  vector<Edge> s_edges_;
  vector<Resource> g_resources_, q_resources_;
  vector<Bucket> s_buckets_;
  vector<int> mu_;
  vector<vector<int>> emb_, indexes_;
  double target_;
  int scheme_rows_, scheme_cols_;
  double threshold_;
  int n_, m_;
  vector<vector<int>> a_, b_;
  vector<vector<int>> bucket_lst_;
  vector<int> signatures_;

  vector<int> tau_;
  vector<vector<int>> iterators_;
  vector<Partition> partitions_;
  vector<int> query_image_;
  vector<int> edges_, w1_, w2_;
  int pos_;
  double result_;
  bool is_init_;
  SubgraphMatching smo_;
  vector<vector<int>> blocks_;
  bool *solution_chk;

  void AddEdgeToSummary(int, int, int, unordered_map<int, unordered_map<int, unordered_map<int, int>>>&);
  void AddEdgeToSummary(int, int, int, unordered_set<Edge,EdgeHasher> &);
  void SummarizeUsingMinHash(DataGraph&, double);
  void MinHashScheme();
  void CreateSummary(DataGraph&, double);
  void MakeBucketListByType();
  void MakeBucketList();
  void CalcSignatureSize();
  void CreateSignature(int);
  void UpdateSignature(int, int, int);
  double Similarity(int, int);
  int SchemeHash(int, int, int, int);
  int BinHash(int);
  void ShallowMerge(int, int);
  void CreateBasePartition();
  void ComputePartitions();
  void TrivialPartition(vector<int>&);
  int FindIncrementable(vector<int>&, vector<int>&);
  int FindDecrementable(vector<int>&);
  bool Includes(Type&, Type&);
  void Evaluate(int);
  double QueryFactor(Partition&);
  double MinimalSolutions(Partition&i);
  double Solutions(Partition&);
  double AtomFactor(Partition&, int);
  int GetWeight(Partition&, int);
  bool IsUnifiable(Partition&, vector<int>&, int);
  void ReplaceGamma(Partition&, int, int);
  void WriteSummaryFile(const char*); 
  size_t DataGraphSize(DataGraph&);
  size_t SummaryGraphSize(DataGraph&);
  int Find(int v) {
    return mu_[v] == v ? v : Find(mu_[v]);
  }
  void Union(int v1, int v2) {
    v1 = Find(v1);
    v2 = Find(v2);
    if (v1 == v2) return;
    if (mu_[v1] > mu_[v2]) std::swap(v1, v2);
    mu_[v2] = v1;
  }
};

#endif

