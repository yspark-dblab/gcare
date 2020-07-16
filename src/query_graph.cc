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

	for (int i = 0; i < edge_.size(); ++i) {
	    edge_enc_[edge_[i]] = i;
	}
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

void QueryGraph::getAll2Paths(vector<tuple<int, int, Edge, Edge>> &result) {
    Edge e;
    int minVId;
    for (int i = 0; i < GetNumEdges(); ++i) {
        e = GetEdge(i);
        minVId = min(e.src, e.dst);
        const vector<pair<int, int>> &srcAdj = GetAdj(e.src, true);
        for (const pair<int, int> &nbr : srcAdj) {
            if (e.dst == nbr.first && e.el == nbr.second) continue;
            if (nbr.first < minVId || nbr.first < e.dst) continue;
            if (e.el < nbr.second) {
                result.emplace_back(make_tuple(Edge::FORWARD, Edge::FORWARD, e, Edge(e.src, nbr.first, nbr.second)));
            } else {
                result.emplace_back(make_tuple(Edge::FORWARD, Edge::FORWARD, Edge(e.src, nbr.first, nbr.second), e));
            }
        }

        const vector<pair<int, int>> &srcInAdj = GetAdj(e.src, false);
        for (const pair<int, int> &nbr : srcInAdj) {
            if (nbr.first < minVId || nbr.first < e.dst) continue;
            result.emplace_back(make_tuple(Edge::FORWARD, Edge::BACKWARD, e, Edge(nbr.first, e.src, nbr.second)));
        }

        const vector<pair<int, int>> &dstAdj = GetAdj(e.dst, true);
        for (const pair<int, int> &nbr : dstAdj) {
            if (nbr.first < minVId || nbr.first < e.src) continue;
            result.emplace_back(make_tuple(Edge::FORWARD, Edge::BACKWARD, Edge(e.dst, nbr.first, nbr.second), e));
        }

        const vector<pair<int, int>> &dstInAdj = GetAdj(e.dst, false);
        for (const pair<int, int> &nbr : dstInAdj) {
            if (e.src == nbr.first && e.el == nbr.second) continue;
            if (nbr.first < minVId || nbr.first < e.src) continue;
            if (e.el < nbr.second) {
                result.emplace_back(make_tuple(Edge::BACKWARD, Edge::BACKWARD, e, Edge(nbr.first, e.dst, nbr.second)));
            } else {
                result.emplace_back(make_tuple(Edge::BACKWARD, Edge::BACKWARD, Edge(nbr.first, e.dst, nbr.second), e));
            }
        }
    }
}

int QueryGraph::encodeSubQ(const vector<Edge> &edges) {
    int enc = 0;
    for (const Edge &e : edges) {
        enc |= (1 << edge_enc_[e]);
    }
    return enc;
}

int QueryGraph::encodeSubQ(const Edge &edge) {
    return (1 << edge_enc_[edge]);
}