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
    adjE_.clear();
    in_adjE_.clear();
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
            if (adj_.size() == 0) {
                adj_.resize(vnum_);
                adjE_.resize(vnum_);
            }
            if (in_adj_.size() == 0) {
                in_adj_.resize(vnum_);
                in_adjE_.resize(vnum_);
            }
            int src = stoi(tok[1]);
            int dst = stoi(tok[2]);
            int el  = stoi(tok[3]);
            assert(src < vnum_);
            assert(dst < vnum_);
            max_el = std::max(max_el, el);
            Edge e = Edge(src, dst, el, edge_.size());
            edge_.push_back(e);
            adj_[src].push_back(make_pair(dst, el));
            in_adj_[dst].push_back(make_pair(src, el));
            adjE_[src].push_back(e);
            in_adjE_[dst].push_back(e);
            //cout << "Q add edge (" << src << ", " << dst << ", " << el << ")" << endl;
            enum_++;
        }
    }
    if (adj_.size() == 0) {
        adj_.resize(vnum_);
        adjE_.resize(vnum_);
    }
    if (in_adj_.size() == 0) {
        in_adj_.resize(vnum_);
        in_adjE_.resize(vnum_);
    }
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

vector<Edge>& QueryGraph::GetAdjE(int v, bool dir) {
    return dir ? adjE_[v] : in_adjE_[v];
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
        const vector<Edge> &srcAdj = GetAdjE(e.src, true);
        for (const Edge &extE : srcAdj) {
            if (e.dst == extE.dst && e.el == extE.el) continue;
            if (extE.dst < minVId || extE.dst < e.dst) continue;
            if (e.el < extE.el) {
                result.emplace_back(make_tuple(Edge::FORWARD, Edge::FORWARD, e, extE));
            } else {
                result.emplace_back(make_tuple(Edge::FORWARD, Edge::FORWARD, extE, e));
            }
        }

        const vector<Edge> &srcInAdj = GetAdjE(e.src, false);
        for (const Edge &extE : srcInAdj) {
            if (extE.src < minVId || extE.src < e.dst) continue;
            result.emplace_back(make_tuple(Edge::FORWARD, Edge::BACKWARD, e, extE));
        }

        const vector<Edge> &dstAdj = GetAdjE(e.dst, true);
        for (const Edge &extE : dstAdj) {
            if (extE.dst < minVId || extE.dst < e.src) continue;
            result.emplace_back(make_tuple(Edge::FORWARD, Edge::BACKWARD, extE, e));
        }

        const vector<Edge> &dstInAdj = GetAdjE(e.dst, false);
        for (const Edge &extE : dstInAdj) {
            if (e.src == extE.src && e.el == extE.el) continue;
            if (extE.src < minVId || extE.src < e.src) continue;
            if (e.el < extE.el) {
                result.emplace_back(make_tuple(Edge::BACKWARD, Edge::BACKWARD, e, extE));
            } else {
                result.emplace_back(make_tuple(Edge::BACKWARD, Edge::BACKWARD, extE, e));
            }
        }
    }
}

int QueryGraph::encodeSubQ(const vector<Edge> &edges) {
    int enc = 0;
    for (const Edge &e : edges) {
        enc |= (1 << e.id);
    }
    return enc;
}

int QueryGraph::encodeSubQ(const Edge &edge) {
    return (1 << edge.id);
}