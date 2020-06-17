#include <fstream>
#include <cassert>
#include <iostream>
#include <string>
#include "../include/query_graph.h"

void QueryGraph::ReadText(const char* fn) {
  fn_ = string(fn);
    vnum_ = enum_ = 0;
    edge_.clear();
    adj_.clear();
    in_adj_.clear();
    vl_.clear();
    bound_.clear();
	ifstream input(fn);
	string line;
	int max_vl = -1, max_el = -1;
    while (getline(input, line)) {
        auto tok = parse(line, " ");
        if (tok[0][0] == 'v') {
			int vl = stoi(tok[2]);
			max_vl = std::max(max_vl, vl);
			vl_.push_back(vl);
			int bound = stoi(tok[3]);
			bound_.push_back(bound);
			vnum_++;
        }
        if (tok[0][0] == 'e') {
            if (adj_.size() == 0)
                adj_.resize(vnum_);
            if (in_adj_.size() == 0)
                in_adj_.resize(vnum_);
            int src = stoi(tok[1]);
            int dst = stoi(tok[2]);
            int el  = stoi(tok[3]);
            assert(src < vnum_);
            assert(dst < vnum_);
            max_el = std::max(max_el, el);
            edge_.emplace_back(src, dst, el);
            adj_[src].push_back(make_pair(dst, el));
            in_adj_[dst].push_back(make_pair(src, el));
            //cout << "Q add edge (" << src << ", " << dst << ", " << el << ")" << endl;
            enum_++;
        }
    }
    if (adj_.size() == 0)
        adj_.resize(vnum_);
    if (in_adj_.size() == 0)
        in_adj_.resize(vnum_);
	vl_num_ = max_vl + 1;
	el_num_ = max_el + 1;
}

void QueryGraph::ReadText(std::vector<std::string> &text) {
    //fn_ = string(fn);
    vnum_ = enum_ = 0;
    edge_.clear();
    adj_.clear();
    in_adj_.clear();
    vl_.clear();
    bound_.clear();
	//ifstream input(fn);
	string line;
	int max_vl = -1, max_el = -1;
    for (auto &line: text) {
    //while (getline(input, line)) {
        auto tok = parse(line, " ");
        if (tok[0][0] == 'v') {
			int vl = stoi(tok[2]);
			max_vl = std::max(max_vl, vl);
			vl_.push_back(vl);
			int bound = stoi(tok[3]);
			bound_.push_back(bound);
			vnum_++;
        }
        if (tok[0][0] == 'e') {
            if (adj_.size() == 0)
                adj_.resize(vnum_);
            if (in_adj_.size() == 0)
                in_adj_.resize(vnum_);
            int src = stoi(tok[1]);
            int dst = stoi(tok[2]);
            int el  = stoi(tok[3]);
            assert(src < vnum_);
            assert(dst < vnum_);
            max_el = std::max(max_el, el);
            edge_.emplace_back(src, dst, el);
            adj_[src].push_back(make_pair(dst, el));
            in_adj_[dst].push_back(make_pair(src, el));
            //cout << "Q add edge (" << src << ", " << dst << ", " << el << ")" << endl;
            enum_++;
        }
    }
    if (adj_.size() == 0)
        adj_.resize(vnum_);
    if (in_adj_.size() == 0)
        in_adj_.resize(vnum_);
	vl_num_ = max_vl + 1;
	el_num_ = max_el + 1;
}

vector<pair<int, int>>& QueryGraph::GetAdj(int v, bool dir) { 
    return dir ? adj_[v] : in_adj_[v]; 
}

//assume one edge label between any two query vertices
int QueryGraph::GetELabel(int u, int v) {
	for (auto& e : adj_[u]) {
		if (e.first == v)
			return e.second;
	}
	return -1;
}
