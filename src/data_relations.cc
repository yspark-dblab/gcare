#include "../include/data_relations.h"

#include <algorithm>
#include <cassert>
#include <vector>
#include <string>
#include <cstring>
#include <ctime>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <map>
#include <iostream>
#include <boost/functional/hash.hpp>

using std::vector;
using std::string;
using std::pair;
using std::make_pair;

off_t StrLen(char* str) {
	off_t len = 0;
	while (*(str + len) != '\n' && *(str + len) != ' ')
		len++;
	return len;
}

int StrToInt(char* str) {
	size_t i = 0;
	int res = 0;
  bool flag = false;
  if (*str == '-') {
    flag = true;
    i++;
  }
	while (*(str + i) != '\n' && *(str + i) != ' ') {
		res *= 10;
		res += static_cast<int>(*(str + i) - '0');
		i++;
	}
  if (flag) res *= -1;
	return res;
}

DataGraph::DataGraph(void) {
	g_.clear();
	table_cnt_ = idx_cnt_ = 0;
    vnum_ = base_ = 0;
}

int DataGraph::Mapping(int v, int t, int c) {
    assert(map_[v].size() % 3 == 0);
    for (size_t i = 0; i < map_[v].size(); i += 3)
        if (map_[v][i] == t && map_[v][i + 1] == c)
            return map_[v][i + 2];
    return -1;
}

template <typename T>
struct ContainerHash {
	size_t operator()(T const& c) const {
		return boost::hash_range(c.begin(), c.end());
	}
};

void DataGraph::ReadText(const char* filename) {
    std::cout << "DataGraph::ReadText " << filename << "\n";
	int fd = open(filename, O_RDONLY, (mode_t)0600);
	struct stat fileinfo;
	if (fstat(fd, &fileinfo) == -1) {
		perror("fstat error");
		close(fd);
		exit(0);
	}
	char* map = static_cast<char*>(mmap(0, fileinfo.st_size, PROT_READ, MAP_SHARED, fd, 0));
	if (map == MAP_FAILED) {
		perror("mmap error");
		close(fd);
		exit(0);
	}
	off_t base = 0;
	size_t gnum = 0;
	unordered_map<vector<int>, int, ContainerHash<vector<int>>> h1, h2; // (t, c, v) -> #
	unordered_map<vector<int>, int, ContainerHash<vector<int>>> cnt1, cnt2;
	while (base < fileinfo.st_size) {
		if (map[base] == '\n') {
			base++;
			continue;
		}
		off_t linelen;
		for (linelen = 0; map[base + linelen] != '\n' && base + linelen < fileinfo.st_size; linelen++) { }
		if (map[base] == 't') {
			g_.push_back(CvtDataGraph());
			gnum = g_.size();
			g_[gnum - 1].max_vid = g_[gnum - 1].max_vlabel = g_[gnum - 1].max_elabel = g_[gnum - 1].base = g_[gnum - 1].vnum = 0;
			h1.clear();
			h2.clear();
			cnt1.clear();
			cnt2.clear();
		} else if (map[base] == 'v') {
			assert(gnum > 0);
			g_[gnum - 1].vnum++;
			int vid = StrToInt(map + base + 2);
			g_[gnum - 1].max_vid = std::max(g_[gnum - 1].max_vid, vid);
			for (off_t offset = 3 + StrLen(map + base + 2); offset < linelen && base + offset < fileinfo.st_size;
					offset += StrLen(map + base + offset) + 1) {
				int vlabel = StrToInt(map + base + offset);
				g_[gnum - 1].max_vlabel = std::max(g_[gnum - 1].max_vlabel, vlabel);

				if (vlabel == -1)
					continue;

				if (static_cast<int>(g_[gnum - 1].vlabel.size()) <= vlabel)
					g_[gnum - 1].vlabel.resize(vlabel + 1);
				g_[gnum - 1].vlabel[vlabel].push_back(vector<int>{ vid });
				table_cnt_++;

				if (static_cast<int>(g_[gnum - 1].vidx.size()) <= vlabel)
					g_[gnum - 1].vidx.resize(vlabel + 1);
				if (static_cast<int>(g_[gnum - 1].vidx[vlabel].size()) < 1)
					g_[gnum - 1].vidx[vlabel].resize(1);
				g_[gnum - 1].vidx[vlabel][0][vid].push_back(g_[gnum - 1].vlabel[vlabel].size() - 1);
				idx_cnt_++;

				if (static_cast<int>(g_[gnum - 1].vmap.size()) <= vid)
					g_[gnum - 1].vmap.resize(vid + 1);
			}
		} else {
			assert(gnum > 0);
			int srcid = StrToInt(map + base + 2);
			int offset = 2;
			offset += StrLen(map + base + offset) + 1;
			int dstid = StrToInt(map + base + offset);
			offset += StrLen(map + base + offset) + 1;
			int elabel = StrToInt(map + base + offset);
			g_[gnum - 1].max_elabel = std::max(g_[gnum - 1].max_elabel, elabel);

			if (static_cast<int>(g_[gnum - 1].edge_table.size()) <= elabel)
				g_[gnum - 1].edge_table.resize(elabel + 1);
			g_[gnum - 1].edge_table[elabel].push_back(vector<int>{ srcid, dstid });
			table_cnt_ += 2;

			if (static_cast<int>(g_[gnum - 1].eidx.size()) <= elabel)
				g_[gnum - 1].eidx.resize(elabel + 1);
			if (static_cast<int>(g_[gnum - 1].eidx[elabel].size()) < 2)
				g_[gnum - 1].eidx[elabel].resize(2);
			g_[gnum - 1].eidx[elabel][0][srcid].push_back(g_[gnum - 1].edge_table[elabel].size() - 1);
			g_[gnum - 1].eidx[elabel][1][dstid].push_back(g_[gnum - 1].edge_table[elabel].size() - 1);
			idx_cnt_ += 2;

			if (static_cast<int>(g_[gnum - 1].emap.size()) <= std::max(srcid, dstid))
				g_[gnum - 1].emap.resize(std::max(srcid, dstid) + 1);
		}
		base += linelen + 1;
	}

	map_cnt_ += g_.size();
	for (size_t i = 0; i < g_.size(); i++) {
		g_[i].base = g_[i].edge_table.size();
		map_cnt_ += g_[i].emap.size() + g_[i].vmap.size();
		h1.clear();
		h2.clear();
		for (size_t j = 0; j < g_[i].eidx.size(); j++) {
			g_[i].eidx_max_key.push_back(vector<int>(g_[i].eidx[j].size(), 0));
			int t = j;
			for (size_t k = 0; k < g_[i].eidx[j].size(); k++) {
				g_[i].eidx_max_key[j][k] = 0;
				int c = k;
				int cnt = 0;
				for (auto& it : g_[i].eidx[j][k]) {
					int v = it.first;
					vector<int> key{t, c, v};
					if (h1[key] == 0) {
						h1[key] = ++cnt;
						g_[i].emap[v].push_back(t);
						g_[i].emap[v].push_back(c);
						g_[i].emap[v].push_back(h1[key] - 1);
						map_cnt_ += 3;
					}
					g_[i].eidx_max_key[j][k] = std::max(g_[i].eidx_max_key[j][k], it.first);
				}
			}
		}
		for (size_t j = 0; j < g_[i].vidx.size(); j++) {
			g_[i].vidx_max_key.push_back(vector<int>(g_[i].vidx[j].size(), 0));
			int t = j + g_[i].base;
			for (size_t k = 0; k < g_[i].vidx[j].size(); k++) {
				g_[i].vidx_max_key[j][k] = 0;
				int c = k;
				int cnt = 0;
				for (auto& it : g_[i].vidx[j][k]) {
					int v = it.first;
					vector<int> key{t, c, v};
					if (h2[key] == 0) {
						h2[key] = ++cnt;
						g_[i].vmap[v].push_back(t);
						g_[i].vmap[v].push_back(c);
						g_[i].vmap[v].push_back(h2[key] - 1);
						map_cnt_ += 3;
					}
					g_[i].vidx_max_key[j][k] = std::max(g_[i].vidx_max_key[j][k], it.first);
				}
			}
		}
	}
	close(fd);
	if (munmap(map, fileinfo.st_size) == -1) {
		perror("error un-mapping the file");
		exit(0);
	}
}

void DataGraph::Make1DTable(const char* dataname) {
    std::cout << "DataGraph::Make1DTable to " << dataname << "\n";
	FILE* fp = fopen(dataname, "w");
	uint64_t rem = table_cnt_ * 2 + 10;
	uint64_t unit = 1024;
	while (rem > 0) {
		vector<int> tmp(std::min(rem, unit), 0);
		rem -= fwrite(tmp.data(), sizeof(int), tmp.size(), fp);
	}
	fclose(fp);

	int fd_table = open(dataname, O_RDWR | O_CREAT, (mode_t)0600);
	struct stat fileinfo_table;
	if (fstat(fd_table, &fileinfo_table) == -1) {
		perror("fstat error");
		close(fd_table);
		exit(0);
	}
	int* ptr_table = static_cast<int*>(mmap(0, fileinfo_table.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd_table, 0));
	if (ptr_table == MAP_FAILED) {
		close(fd_table);
		perror("Error mapping the file");
		exit(0);
	}

	off_t end = g_.size() + 1;
	off_t size = end - 1;
	for (size_t i = 0; i < g_.size(); i++) { // t
		ptr_table[i] = end - i;
		end += g_[i].edge_table.size() + g_[i].vlabel.size();
	}
	ptr_table[size] = end - size;
	end++;
	off_t pos = size + 1;
	size = end - 1;
	for (size_t i = 0; i < g_.size(); i++) { // r
		for (size_t j = 0; j < g_[i].edge_table.size(); j++) {
			assert(pos < fileinfo_table.st_size);
			ptr_table[pos] = end - pos;
			pos++;
			end += g_[i].edge_table[j].size();
		}
		for (size_t j = 0; j < g_[i].vlabel.size(); j++) {
			assert(pos < fileinfo_table.st_size);
			ptr_table[pos] = end - pos;
			pos++;
			end += g_[i].vlabel[j].size();
		}
	}
	ptr_table[size] = end - size;
	end++;
	pos = size + 1;
	size = end - 1;
	for (size_t i = 0; i < g_.size(); i++) {
		for (size_t j = 0; j < g_[i].edge_table.size(); j++) {
			for (size_t k = 0; k < g_[i].edge_table[j].size(); k++) {
				ptr_table[pos] = end - pos;
				pos++;
				for (size_t l = 0; l < g_[i].edge_table[j][k].size(); l++)
					ptr_table[end++] = g_[i].edge_table[j][k][l];
			}
		}
		for (size_t j = 0; j < g_[i].vlabel.size(); j++) {
			for (size_t k = 0; k < g_[i].vlabel[j].size(); k++) {
				ptr_table[pos] = end - pos;
				pos++;
				for (size_t l = 0; l < g_[i].vlabel[j][k].size(); l++)
					ptr_table[end++] = g_[i].vlabel[j][k][l];
			}
		}
	}
	ptr_table[size] = end - size;
	if (msync(ptr_table, fileinfo_table.st_size, MS_SYNC)) {
		perror("fail to sync the file to disk");
		exit(0);
	}
	close(fd_table);
	if (munmap(ptr_table, fileinfo_table.st_size) == -1) {
		perror("error un-mapping the file");
		exit(0);
	}
    std::cout << "~DataGraph::Make1DTable to " << dataname << "\n";
}

void DataGraph::Make1DIndex(char* dataname) {
    std::cout << "DataGraph::Make1DIndex to " << dataname << "\n";
	string fn = static_cast<string>(dataname) + ".index";

	FILE* fp = fopen(fn.c_str(), "w");
	uint64_t rem = idx_cnt_ * 2 + 1;
	uint64_t unit = 1024;
	while (rem > 0) {
		vector<int> tmp(std::min(rem, unit), 0);
		rem -= fwrite(tmp.data(), sizeof(int), tmp.size(), fp);
	}
	fclose(fp);

	int fd_idx = open(fn.c_str(), O_RDWR | O_CREAT, (mode_t)0600);
	struct stat fileinfo_idx;
	if (fstat(fd_idx, &fileinfo_idx) == -1) {
		perror("fstat error");
		close(fd_idx);
		exit(0);
	}
	int* ptr_idx = static_cast<int*>(mmap(0, fileinfo_idx.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd_idx, 0));
	if (ptr_idx == MAP_FAILED) {
		close(fd_idx);
		perror("Error mapping the file");
		exit(0);
	}

	off_t end = g_.size();
	for (size_t i = 0; i < g_.size(); i++) { // t
		ptr_idx[i] = end - i;
		end += g_[i].eidx.size() + g_[i].vidx.size();
	}
	off_t size = end;
	off_t pos = g_.size();
	for (size_t i = 0; i < g_.size(); i++) { // c
		for (size_t j = 0; j < g_[i].eidx.size(); j++) {
			assert(pos < fileinfo_idx.st_size);
			ptr_idx[pos] = end - pos;
			pos++;
			end += g_[i].eidx[j].size();
		}
		for (size_t j = 0; j < g_[i].vidx.size(); j++) {
			assert(pos < fileinfo_idx.st_size);
			ptr_idx[pos] = end - pos;
			pos++;
			end += g_[i].vidx[j].size();
		}
	}
	pos = size;
	size = end;
	for (size_t i = 0; i < g_.size(); i++) { // v
		for (size_t j = 0; j < g_[i].eidx.size(); j++) {
			for (size_t k = 0; k < g_[i].eidx[j].size(); k++) {
				ptr_idx[pos] = end - pos;
				pos++;
				end += g_[i].eidx[j][k].size();
			}
		}
		for (size_t j = 0; j < g_[i].vidx.size(); j++) {
			for (size_t k = 0; k < g_[i].vidx[j].size(); k++) {
				ptr_idx[pos] = end - pos;
				pos++;
				end += g_[i].vidx[j][k].size();
			}
		}
	}
	end++;
	pos = size;
	size = end - 1;
	for (size_t i = 0; i < g_.size(); i++) { // i
		for (size_t j = 0; j < g_[i].eidx.size(); j++) {
			for (size_t k = 0; k < g_[i].eidx[j].size(); k++) {
				for (auto& it : g_[i].eidx[j][k]) {
					ptr_idx[pos] = end - pos;
					pos++;
					for (int e : it.second)
						ptr_idx[end++] = e;
				}
			}
		}
		for (size_t j = 0; j < g_[i].vidx.size(); j++) {
			for (size_t k = 0; k < g_[i].vidx[j].size(); k++) {
				for (auto& it : g_[i].vidx[j][k]) {
					ptr_idx[pos] = end - pos;
					pos++;
					for (int e : it.second)
						ptr_idx[end++] = e;
				}
			}
		}
	}
	ptr_idx[size] = end - size;
	if (msync(ptr_idx, fileinfo_idx.st_size, MS_SYNC)) {
		perror("fail to sync the file to disk");
		exit(0);
	}
	close(fd_idx);
	if (munmap(ptr_idx, fileinfo_idx.st_size) == -1) {
		perror("error un-mapping the file");
		exit(0);
	}

	fn = static_cast<string>(dataname) + ".map";

	fp = fopen(fn.c_str(), "w");
	rem = map_cnt_ + 1;
	unit = 1024;
	std::cout << "rem: " << rem << std::endl;
	while (rem > 0) {
		vector<int> tmp(std::min(rem, unit), 0);
		rem -= fwrite(tmp.data(), sizeof(int), tmp.size(), fp);
	}
	fclose(fp);

	int fd_map = open(fn.c_str(), O_RDWR | O_CREAT, (mode_t)0600);
	struct stat fileinfo_map;
	if (fstat(fd_map, &fileinfo_map) == -1) {
		perror("fstat error");
		close(fd_map);
		exit(0);
	}
	int* ptr_map = static_cast<int*>(mmap(0, fileinfo_map.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd_map, 0));
	if (ptr_map == MAP_FAILED) {
		close(fd_map);
		perror("Error mapping the file");
		exit(0);
	}

	end = g_.size();
	for (size_t i = 0; i < g_.size(); i++) { // v
		ptr_map[i] = end - i;
		end += std::max(g_[i].emap.size(), g_[i].vmap.size());
	}
	end++;
	size = end - 1;
	pos = g_.size();
	for (size_t i = 0; i < g_.size(); i++) {
		for (size_t v = 0; v < std::max(g_[i].emap.size(), g_[i].vmap.size()); v++) {
			assert(pos < fileinfo_map.st_size);
			ptr_map[pos] = end - pos;
			pos++;
			if (v < g_[i].emap.size())
				for (size_t j = 0; j < g_[i].emap[v].size(); j++)
					ptr_map[end++] = g_[i].emap[v][j];
			if (v < g_[i].vmap.size())
				for (size_t j = 0; j < g_[i].vmap[v].size(); j++)
					ptr_map[end++] = g_[i].vmap[v][j];
		}
	}
	ptr_map[size] = end - size;
	if (msync(ptr_map, fileinfo_map.st_size, MS_SYNC)) {
		perror("fail to sync the file to disk");
		exit(0);
	}
	close(fd_map);
	if (munmap(ptr_map, fileinfo_map.st_size) == -1) {
		perror("error un-mapping the file");
		exit(0);
	}
}

void DataGraph::WriteBinary(const char* dataname) {
  string fname = string(dataname) + ".relation";
  std::cout << "DataGraph::WriteBinary to " << fname << "\n";
  Make1DTable(fname.c_str());
	string fn = fname + ".meta";
	FILE* fp = fopen(fn.c_str(), "w");
	fprintf(fp, "%zu\n", g_.size());
	for (size_t i = 0; i < g_.size(); i++)
		fprintf(fp, "%zu %d %d %d\n", g_[i].edge_table.size(), g_[i].max_vid, g_[i].max_vlabel, g_[i].max_elabel);
	fclose(fp);
    std::cout << "~DataGraph::WriteBinary to " << fname << "\n";
}

void DataGraph::ReadBinary(const char* dataname) {
  string fname = string(dataname) + ".relation";
  std::cout << "DataGraph::ReadBinary from " << fname << "\n";
	string metadata = fname + ".meta";
	FILE* fp = fopen(metadata.c_str(), "r");
	int gnum;
	fscanf(fp, "%d", &gnum);
	fscanf(fp, "%d%d%d%d", &base_, &max_vid_, &max_vlabel_, &max_elabel_);
	fclose(fp);

	int fd = open(fname.c_str(), O_RDONLY, (mode_t)0600);
	struct stat fileinfo;
	fstat(fd, &fileinfo);
	int* ptr = static_cast<int*>(mmap(0, fileinfo.st_size, PROT_READ, MAP_SHARED, fd, 0));
    container_ = static_cast<int*>(malloc(fileinfo.st_size));
    memcpy(container_, ptr, fileinfo.st_size);

    table_.array_ = container_ + container_[0];
    table_.size_ = container_[1] - container_[0] + 1;
	if (munmap(ptr, fileinfo.st_size) == -1) {
		perror("error un-mapping the file");
		exit(0);
	}
	close(fd);
    std::cout << "~DataGraph::ReadBinary from " << fname << "\n";
}
