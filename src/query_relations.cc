#include "../include/query_relations.h"
#include "../include/ndvector.h"

#include <algorithm>
#include <vector>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <string>
#include <cassert>
#include <iostream>

using std::vector;
using std::string;
using std::pair;
using std::make_pair;

/*
QueryGraph::QueryGraph(void) {
	table_.clear();
	condition_.clear();
	par_.clear();
	cnt_.clear();
}

int QueryGraph::Find(int x) {
	return (par_[x] == x) ? x : (par_[x] = Find(par_[x]));
}

void QueryGraph::Unite(int a, int b) {
	if ((a = QueryGraph::Find(a)) == (b = QueryGraph::Find(b))) return;
    (a > b) ? (par_[a] = b) : (par_[b] = a);
}

void QueryGraph::ReadText(const char* filename) {
    //std::cout << "QueryGraph::ReadText " << filename << "\n";
	FILE* fp = fopen(filename, "r");
	char buf[1111];
	int cnt = 0;
	vector<int> vlabel_vec, vid_vec;
    max_vid_ = max_vlabel_ = max_elabel_ = enum_ = 0;
	while (fgets(buf, 1111, fp) != NULL) {
		if (buf[0] == '\n') continue;
        //std::cout << "buf: " << buf;
		int buflen = 0;
		while (buf[buflen++] != '\n') {
			if (buf[buflen] == ' ') {
				buf[buflen] = 0;
            }
        }
		if (buf[0] == 'v') {
			int vid = strtol(buf + 2, 0, 0);
			max_vid_ = std::max(max_vid_, vid);
			int offset = 2 + strlen(buf + 2) + 1;
			int vlabel = strtol(buf + offset, 0, 0);
			max_vlabel_ = std::max(max_vlabel_, vlabel);
			offset += strlen(buf + offset) + 1;
			int vidattr = -1;
			if (offset < buflen)
				vidattr = strtol(buf + offset, 0, 0);
			if (vid >= static_cast<int>(vlabel_vec.size()))
				vlabel_vec.resize(vid + 1);
			if (vid >= static_cast<int>(vid_vec.size()))
				vid_vec.resize(vid + 1);
			vlabel_vec[vid] = vlabel + 1;
			vid_vec[vid] = vidattr;
			if (vlabel != -1) {
				table_.push_back(-vlabel_vec[vid]);
				condition_.push_back(vector<QueryGraph::Column>());
				int t = table_.size() - 1;
				condition_[t].push_back(QueryGraph::Column(vid, cnt, -1, -1));
				par_.push_back(cnt++);
			}
		} else if (buf[0] == 'e') {
			int srcid = strtol(buf + 2, 0, 0);
			int offset = 2 + strlen(buf + 2) + 1;
			int dstid = strtol(buf + offset, 0, 0);
			offset += strlen(buf + offset) + 1;
			int elabel = strtol(buf + offset, 0, 0);
			enum_++;
			max_elabel_ = std::max(max_elabel_, elabel);
			table_.push_back(elabel);
			condition_.push_back(vector<QueryGraph::Column>());
			int t = table_.size() - 1;
			condition_[t].push_back(
					QueryGraph::Column(srcid, cnt, vlabel_vec[srcid], vid_vec[srcid]));
			condition_[t].push_back(
					QueryGraph::Column(dstid, cnt + 1, vlabel_vec[dstid], vid_vec[dstid]));
			par_.push_back(cnt);
		    par_.push_back(cnt + 1);
			cnt += 2;
		}
	}
	for (size_t j = 0; j < condition_.size() - 1; j++) {
		for (size_t k = j + 1; k < condition_.size(); k++) {
			for (auto& x: condition_[j]) {
				for (auto& y: condition_[k]) {
					if (x.vid == y.vid) Unite(x.cls, y.cls);
                }
            }
        }
    }
	cnt = 0;
	cnt_.resize(par_.size());
	for (size_t j = 0; j < condition_.size(); j++) {
		for (size_t k = 0; k < condition_[j].size(); k++) {
			condition_[j][k].cls = par_[cnt++];
			cnt_[condition_[j][k].cls]++;
		}
	}
	fclose(fp);
    //std::cout << "~QueryGraph::ReadText " << filename << "\n";
}
*/

