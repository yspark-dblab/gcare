#ifndef QUERY_GRAPH_H_
#define QUERY_GRAPH_H_

#include <vector>
#include <string>
#include <unordered_map>
#include "util.h"

class QueryGraph {
private:
	int vnum_, enum_, vl_num_, el_num_; 
	vector<Edge> edge_;
	vector<vector<pair<int, int>>> adj_, in_adj_;
	vector<int> vl_;
	vector<int> bound_;
	unordered_map<Edge, int> edge_enc_;

public:
	QueryGraph() { 
		vnum_ = enum_ = vl_num_ = el_num_ = 0;
		edge_.clear();
		adj_.clear();
		in_adj_.clear();
		vl_.clear();
		bound_.clear();
        edge_enc_.clear();
	}
	void ReadText(const char*);
    void ReadText(std::vector<std::string>&);
	inline int GetNumVertices() { return vnum_; }
	inline int GetNumEdges() { return enum_; }
	Edge GetEdge(int i) { return edge_[i]; }; 
	vector<pair<int, int>>& GetAdj(int, bool);
	int GetELabel(int, int);
	inline int GetVLabel(int v) { return vl_[v]; }
	inline int GetBound(int v) { return bound_[v]; }

    void getAll2Paths(vector<tuple<int, int, Edge, Edge>> &result);
	int encodeSubQ(const vector<Edge> &edges);

  string fn_; // XXX

  pair<string, string> toVListAndLabelSeq() {
      string vList = "";
      string labelSeq = "";
      for (int i = 0; i < edge_.size(); ++i) {
          pair<string, string> e = edge_[i].toVListAndLabelSeq();
          vList += ";" + e.first;
          labelSeq += "->" + e.second;
      }
      return make_pair(vList.substr(1), labelSeq.substr(2));
  }
};

#endif
