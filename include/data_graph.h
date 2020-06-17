#ifndef GRAPH_H_
#define GRAPH_H_

#include <algorithm>
#include <vector>
#include <map>
#include <iostream>
#include <unordered_map>
#include "util.h"

using namespace std;

//use for making binary
struct RawDataGraph {
	int max_vl_, max_el_;
	vector<int> vl_cnt_, el_cnt_;

	vector<vector<int>> vlabels_;
	vector<Edge> out_edges_; 
	vector<Edge> in_edges_; 

	vector<int> offset_;
	vector<int> label_;
	vector<int> adj_offset_;
	vector<int> adj_;

	vector<int> in_offset_;
	vector<int> in_label_;
	vector<int> in_adj_offset_;
	vector<int> in_adj_;

	vector<int> vl_offset_;
	vector<int> vl_;

	vector<vector<int>> vl_rel_;
	vector<vector<pair<int, int>>> el_rel_;
};

class DataGraph {
private:
	int vnum_, enum_, vl_num_, el_num_; 

	vector<int> vl_cnt_, el_cnt_;
	
	size_t encode_size_;

	const int* offset_; //data vertex id -> offset 
	//const pair<int, int>* label_;
	const int* label_;
	const int* adj_offset_;
	const int* adj_;

	const int* in_offset_; //data vertex id -> offset 
	//const pair<int, int>* in_label_;
	const int* in_label_;
	const int* in_adj_offset_;
	const int* in_adj_;

	//const int* rel_offset_; 
	//const int* rel_;

	//const int* in_rel_offset_; 
	//const int* in_rel_;
	
	const int* vl_offset_;
	const int* vl_;

	//const int* vl_offset_; 
	//const int* vl_rel_;
	
	const int* el_rel_offset_;
	const pair<int, int>* el_rel_;

	const int* vl_rel_offset_;
	const int* vl_rel_;
	
	RawDataGraph raw_;
		
public:
	//build mode
	bool HasBinary(const char*);
	void ReadText(const char*);
	void MakeBinary();
	size_t BinarySize();
	void WriteBinary(const char*);
	void ClearRawData();

	void ReadBinary(const char*);
	int GetNumVertices();
	int GetNumVertices(int);
	int GetNumEdges();
	int GetNumEdges(int);
	int GetNumVLabels(int = -1); 
	int GetNumELabels(int = -1, bool = true); 
	range GetVLabels(int);
	range GetELabels(int, bool);
	int GetELabelIndex(int, int, bool);
	range GetAdj(int, int, bool);
	//range GetRel(int, bool);
	//range GetUni(int);
	int   GetAdjSize(int, int, bool);
	//int   GetRelSize(int, bool);
	//int   GetUniSize(int);
	bool  HasEdge(int, int, int, bool);
	bool  HasVLabel(int, int);
	bool  HasELabel(int, int, bool);
	vector<int> GetRandomVertex(int);
	vector<int> GetVertex(int, int);
	vector<int> GetRandomEdge(int);
	vector<int> GetEdge(int, int);
	vector<int> GetRandomEdge(int, int, bool);
};

#endif
