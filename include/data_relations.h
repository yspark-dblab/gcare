#ifndef DATA_RELATIONS_H_
#define DATA_RELATIONS_H_

#include "ndvector.h"

#include <algorithm>
#include <vector>
#include <unordered_map>
#include <experimental/filesystem>

using std::vector;
using std::unordered_map;
class DataGraph {

 public:
//converter--
  struct CvtDataGraph {
    Vector3D table;
    Vector4D index;

    vector<vector<vector<int>>> edge_table; // ith table, jth row, k column(=2)
    vector<vector<vector<int>>> vlabel; // ith table, jth row, k column(=1)
    vector<vector<unordered_map<int, vector<int>>>> vidx, eidx; // ith table, jth column, k value, lth element
    vector<vector<int>> vmap, emap; // i value, (t, c, v) list

    vector<vector<int>> vidx_max_key, eidx_max_key;

    int vnum;
    int max_vid, max_vlabel, max_elabel;
    int base;
  };
  
  int* container_;
  vector<CvtDataGraph> g_;
  vector<int> vec1d_table_, vec1d_index_;
  uint64_t table_cnt_, idx_cnt_, map_cnt_;

  void Make1DTable(const char*);
  void Make1DIndex(char*);
//--converter
  int vnum_;
  Vector2D map_;
  Vector3D table_;
  Vector4D index_;
  int base_;
  int max_vid_, max_vlabel_, max_elabel_;
  DataGraph(void);
  int Mapping(int, int, int);

  int get_table_id(int _id) { return _id < 0 ? base_ - _id - 1 : _id; }


  void WriteBinary(char*);

  void ReadText(const char* fn);
  bool HasBinary(const char* fn) { 
    std::string metadata = std::string(fn) + ".relation.meta";
    return std::experimental::filesystem::exists(metadata.c_str());
  }
  void MakeBinary() { }
  void WriteBinary(const char* filename);
  void ReadBinary(const char* fname);
  void ClearRawData() {}
/*
  vector<int> GetRandomTuple(int tid) {
    vector<int> ret;
    auto &table = table_[tid];
    if (table.size() == 0) return ret;
    int r = (rand() % table_[tid].size());
    auto &record = table[r];
    for (size_t i = 0; i < record.size(); ++i) {
      ret.push_back(record[i]);
    }
    return ret;
  }

  size_t GetNumTuples(int tid) {
    return table_[tid].size();
  }

  size_t GetNumTuples(int tid, int pos, int v) {
    // if it has index, use index
    
    // else, do scan
    auto &table = table_[tid];
    size_t cnt = 0;
    for (size_t i = 0; i < table.size(); ++i) {
      cnt += table[pos] == v ? 1 : 0;
    }
    return cnt;
  }

  bool HasValue(int tid, int pos, int v) {
    // if it has index, use index
    
    // else, do scan
    auto &table = table_[tid];
    size_t cnt = 0;
    for (size_t i = 0; i < table.size(); ++i) {
      if (table[pos] == v) return true;
    }
    return false;
  }
*/
};

#endif
