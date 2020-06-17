#include "../include/data_graph.h"

#include <algorithm>
#include <fstream>
#include <vector>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <string>
#include <cstring>
#include <cassert>
#include <experimental/filesystem>
#include <iostream>
void DataGraph::ReadText(const char* fn) {
    std::cout << "DataGraph::ReadText " << fn << "\n";
	raw_.max_vl_ = raw_.max_el_ = -1;
	ifstream input(fn);
	string line;

    uint64_t num = 0;
	while (getline(input, line)) {
		auto tok = parse(line, " ");
		if (tok[0][0] == 'v') {
			vector<int> labels;
			for (int i = 2; i < tok.size(); i++) {
				int vl = stoi(tok[i]);
        if (vl < 0) continue;
				labels.push_back(vl);
				raw_.max_vl_ = std::max(raw_.max_vl_, vl);
			}
			raw_.vlabels_.push_back(labels);
		}
		if (tok[0][0] == 'e') {
			int src = stoi(tok[1]);
			int dst = stoi(tok[2]);
			for (int i = 3; i < tok.size(); i++) {
				int el = stoi(tok[i]);
				raw_.out_edges_.emplace_back(src, dst, el);
				raw_.in_edges_.emplace_back(dst, src, el);
				raw_.max_el_ = std::max(raw_.max_el_, el);
			}
		}
        ++num;
        if (num % 10000000 == 0) std::cout << "Parse " << num << " lines...\n";
	}
	raw_.el_cnt_.resize(raw_.max_el_ + 1);
	raw_.vl_cnt_.resize(raw_.max_vl_ + 1);
	raw_.el_rel_.resize(raw_.max_el_ + 1);
	raw_.vl_rel_.resize(raw_.max_vl_ + 1);
    std::cout << "~DataGraph::ReadText " << fn << "\n";
}

void DataGraph::MakeBinary() {
    std::cout << "DataGraph::MakeBinary\n";
	int vnum = raw_.vlabels_.size();
	auto& out = raw_.out_edges_;
	auto& in = raw_.in_edges_;

	sort(out.begin(), out.end());
	sort(in.begin(), in.end());
	out.erase(unique(out.begin(), out.end()), out.end());
	in.erase(unique(in.begin(), in.end()), in.end());

	{
		auto& offsets = raw_.offset_;
		auto& elabels = raw_.label_;
		auto& adj_offsets = raw_.adj_offset_;
		auto& adj = raw_.adj_;

		int prev_u = -1, prev_el = -1;
		for (auto& e : out) {
			if (e.src != prev_u) {
				for (int u = prev_u + 1; u <= e.src; u++) {
					offsets.push_back(elabels.size());
				}
				elabels.push_back(e.el);
				adj_offsets.push_back(adj.size()); 
				adj.push_back(e.dst);

				prev_u = e.src;
				prev_el = e.el;
			} else if (e.el != prev_el) {
				elabels.push_back(e.el);
				adj_offsets.push_back(adj.size()); 
				adj.push_back(e.dst);

				prev_el = e.el;
			} else {
				adj.push_back(e.dst);
			}
			raw_.el_cnt_[e.el]++;
		}

		for (int u = prev_u + 1; u <= vnum; u++)
			offsets.push_back(elabels.size());
		elabels.push_back(-1);
		adj_offsets.push_back(adj.size());
	}

	{
		auto& offsets = raw_.in_offset_;
		auto& elabels = raw_.in_label_;
		auto& adj_offsets = raw_.in_adj_offset_;
		auto& adj = raw_.in_adj_;

		int prev_u = -1, prev_el = -1;
		for (auto& e : in) {
			if (e.src != prev_u) {
				for (int u = prev_u + 1; u <= e.src; u++) {
					offsets.push_back(elabels.size());
				}
				elabels.push_back(e.el);
				adj_offsets.push_back(adj.size()); 
				adj.push_back(e.dst);

				prev_u = e.src;
				prev_el = e.el;
			} else if (e.el != prev_el) {
				elabels.push_back(e.el);
				adj_offsets.push_back(adj.size()); 
				adj.push_back(e.dst);

				prev_el = e.el;
			} else {
				adj.push_back(e.dst);
			}
			//raw_.el_cnt_[e.el]++;
		}

		for (int u = prev_u + 1; u <= vnum; u++)
			offsets.push_back(elabels.size());
		elabels.push_back(-1);
		adj_offsets.push_back(adj.size());
	}

	int num_vl = 0;
	for (int u = 0; u < vnum; u++) {
		raw_.vl_offset_.push_back(num_vl);
		sort(raw_.vlabels_[u].begin(), raw_.vlabels_[u].end());
		num_vl += raw_.vlabels_[u].size();
		raw_.vl_.insert(raw_.vl_.end(), raw_.vlabels_[u].begin(), raw_.vlabels_[u].end()); 
		for (auto& vl : raw_.vlabels_[u]) {
			raw_.vl_cnt_[vl]++;
			raw_.vl_rel_[vl].push_back(u);
		}
	}
	raw_.vl_offset_.push_back(num_vl);

	{
		//stable_sort(out.begin(), out.end(), [](const Edge& a, const Edge& b) -> bool { return a.el < b.el; });
		for (auto& e : out)
			raw_.el_rel_[e.el].push_back(make_pair(e.src, e.dst));
	}
    std::cout << "~DataGraph::MakeBinary\n";
}

size_t DataGraph::BinarySize() {
	size_t ret = 0;

	size_t vn = raw_.vlabels_.size();

	ret += sizeof(int) * (vn + 1); 
	ret += sizeof(int);
	ret += sizeof(int) * raw_.label_.size() * 2;
	ret += sizeof(int);
	ret += sizeof(int) * raw_.adj_.size();

	ret += sizeof(int) * (vn + 1); 
	ret += sizeof(int);
	ret += sizeof(int) * raw_.in_label_.size() * 2;
	ret += sizeof(int);
	ret += sizeof(int) * raw_.in_adj_.size();

	ret += sizeof(int) * (vn + 1); 
	ret += sizeof(int);
	ret += sizeof(int) * raw_.vl_.size();

	ret += sizeof(int) * (raw_.max_el_ + 2);
	ret += sizeof(pair<int, int>) * raw_.out_edges_.size(); 

	ret += sizeof(int) * (raw_.max_vl_ + 2); 
	ret += sizeof(int) * raw_.vl_.size();

	return ret;
}

bool DataGraph::HasBinary(const char* filename) {
  string metadata = string(filename) + ".graph.meta";
	return std::experimental::filesystem::exists(metadata.c_str());
}

void DataGraph::WriteBinary(const char* filename) {
  string fname = string(filename) + ".graph";
  std::cout << "DataGraph::WriteBinary" << fname << "\n";
	string metadata = fname + ".meta";
	FILE* fp = fopen(metadata.c_str(), "w");
	size_t encode_size = BinarySize();

	int vn = raw_.vlabels_.size();
	int en = raw_.out_edges_.size(); 

	fprintf(fp, "%d %d %d %d %zu\n", vn, en, raw_.max_vl_+1, raw_.max_el_+1, encode_size);
	for (int vl = 0; vl <= raw_.max_vl_; vl++)
		fprintf(fp, "%d ", raw_.vl_cnt_[vl]); 
	fprintf(fp, "\n");
	for (int el = 0; el <= raw_.max_el_; el++)
		fprintf(fp, "%d ", raw_.el_cnt_[el]); 
	fprintf(fp, "\n");
	fclose(fp);

	char* buffer = new char[encode_size];
	char* orig = buffer;
	int size[1] = {0};

	{
		assert(raw_.offset_.size() == vn + 1);
		assert(raw_.offset_.back() + 1 == raw_.label_.size());
		size_t offset_size = sizeof(int) * (vn + 1); 
		memcpy(buffer, raw_.offset_.data(), offset_size); 
		buffer += offset_size;

		size[0] = raw_.label_.size();
		memcpy(buffer, size, sizeof(int));
		buffer += sizeof(int);

		size_t label_size = sizeof(int) * raw_.label_.size(); 
		assert(raw_.label_.size() == raw_.adj_offset_.size());
		memcpy(buffer, raw_.label_.data(), label_size);
		buffer += label_size;
		memcpy(buffer, raw_.adj_offset_.data(), label_size);
		buffer += label_size;

		size[0] = raw_.adj_.size();
		memcpy(buffer, size, sizeof(int));
		buffer += sizeof(int);

		size_t adj_size = sizeof(int) * raw_.adj_.size();
		assert(raw_.adj_offset_.back() == raw_.adj_.size());
		memcpy(buffer, raw_.adj_.data(), adj_size);
		buffer += adj_size;
	}

	{
		assert(raw_.in_offset_.size() == vn + 1);
		assert(raw_.in_offset_.back() + 1 == raw_.in_label_.size());
		size_t offset_size = sizeof(int) * (vn + 1); 
		memcpy(buffer, raw_.in_offset_.data(), offset_size); 
		buffer += offset_size;

		size[0] = raw_.in_label_.size();
		memcpy(buffer, size, sizeof(int));
		buffer += sizeof(int);

		size_t label_size = sizeof(int) * raw_.in_label_.size(); 
		assert(raw_.in_label_.size() == raw_.in_adj_offset_.size());
		memcpy(buffer, raw_.in_label_.data(), label_size);
		buffer += label_size;
		memcpy(buffer, raw_.in_adj_offset_.data(), label_size);
		buffer += label_size;

		size[0] = raw_.in_adj_.size();
		memcpy(buffer, size, sizeof(int));
		buffer += sizeof(int);

		size_t adj_size = sizeof(int) * raw_.in_adj_.size();
		assert(raw_.in_adj_offset_.back() == raw_.in_adj_.size());
		memcpy(buffer, raw_.in_adj_.data(), adj_size);
		buffer += adj_size;
	}

	{
		assert(raw_.vl_offset_.size() == vn + 1);
		size_t offset_size = sizeof(int) * raw_.vl_offset_.size();
		memcpy(buffer, raw_.vl_offset_.data(), offset_size); 
		buffer += offset_size;

		size[0] = raw_.vl_.size();
		memcpy(buffer, size, sizeof(int));
		buffer += sizeof(int);

		size_t label_size = sizeof(int) * raw_.vl_.size();
		memcpy(buffer, raw_.vl_.data(), label_size);
		buffer += label_size;
	}

	{
		assert(raw_.el_rel_.size() == raw_.max_el_ + 1);
		size[0] = 0;  
		for (int el = 0; el <= raw_.max_el_; el++) {
			memcpy(buffer, size, sizeof(int)); 
			size[0] += raw_.el_rel_[el].size();
			buffer += sizeof(int);
		}
		memcpy(buffer, size, sizeof(int));
		buffer += sizeof(int);
		for (int el = 0; el <= raw_.max_el_; el++) {
			memcpy(buffer, raw_.el_rel_[el].data(), sizeof(pair<int, int>) * raw_.el_rel_[el].size()); 
			buffer += sizeof(pair<int, int>) * raw_.el_rel_[el].size();
		}
	}

	{
		assert(raw_.vl_rel_.size() == raw_.max_vl_ + 1);
		size[0] = 0;  
		for (int vl = 0; vl <= raw_.max_vl_; vl++) {
			memcpy(buffer, size, sizeof(int)); 
			size[0] += raw_.vl_rel_[vl].size();
			buffer += sizeof(int);
		}
		memcpy(buffer, size, sizeof(int));
		buffer += sizeof(int);
		for (int vl = 0; vl <= raw_.max_vl_; vl++) {
			memcpy(buffer, raw_.vl_rel_[vl].data(), sizeof(int) * raw_.vl_rel_[vl].size()); 
			buffer += sizeof(int) * raw_.vl_rel_[vl].size();
		}
	}

	assert((buffer - orig) == encode_size);
	FILE* f = fopen(fname.c_str(), "w");
	fwrite(orig, 1, encode_size, f);
	cout << "wrote " << (buffer - orig) << " bytes to file " << fname << endl;
	fclose(f);
    std::cout << "~DataGraph::WriteBinary" << fname << "\n";
}

void DataGraph::ReadBinary(const char* filename) {
  string fname = string(filename) + ".graph";
    std::cout << "DataGraph::ReadBinary" << fname << "\n";
	string metadata = fname + ".meta";
	FILE* fp = fopen(metadata.c_str(), "r");
	size_t encode_size;
	fscanf(fp, "%d %d %d %d %zu", &vnum_, &enum_, &vl_num_, &el_num_, &encode_size);

	vl_cnt_.resize(vl_num_);
	el_cnt_.resize(el_num_);
	for (int vl = 0; vl < vl_num_; vl++)
		fscanf(fp, "%d ", &vl_cnt_[vl]); 
	for (int el = 0; el < el_num_; el++)
		fscanf(fp, "%d ", &el_cnt_[el]); 
	fclose(fp);

	FILE* f = fopen(fname.c_str(), "r");
	char* buffer = new char[encode_size];
	fread(buffer, 1, encode_size, f);
	fclose(f);
	char* orig = buffer;

	{
		offset_ = (const int*) buffer;
		assert(offset_[0] == 0);
		buffer += sizeof(int) * (vnum_ + 1);

		int* size = (int*) buffer;
		assert(offset_[vnum_] + 1 == size[0]);
		buffer += sizeof(int);

		label_ = (const int*) buffer;
		buffer += sizeof(int) * size[0];
		adj_offset_ = (const int*) buffer;
		buffer += sizeof(int) * size[0];

		size = (int*) buffer;
		buffer += sizeof(int);

		adj_ = (const int*) buffer;
		buffer += sizeof(int) * size[0];
	}

	{
		in_offset_ = (const int*) buffer;
		assert(in_offset_[0] == 0);
		buffer += sizeof(int) * (vnum_ + 1);

		int* size = (int*) buffer;
		assert(in_offset_[vnum_] + 1 == size[0]);
		buffer += sizeof(int);

		in_label_ = (const int*) buffer;
		buffer += sizeof(int) * size[0];
		in_adj_offset_ = (const int*) buffer;
		buffer += sizeof(int) * size[0];

		size = (int*) buffer;
		buffer += sizeof(int);

		in_adj_ = (const int*) buffer;
		buffer += sizeof(int) * size[0];
	}

	{
		vl_offset_ = (const int*) buffer;
		assert(vl_offset_[0] == 0);
		buffer += sizeof(int) * (vnum_ + 1);

		int* size = (int*) buffer;
		buffer += sizeof(int);

		vl_ = (const int*) buffer;
		buffer += sizeof(int) * size[0];
	}

	{
		el_rel_offset_ = (const int*) buffer; 
		assert(el_rel_offset_[0] == 0);
		buffer += sizeof(int) * (el_num_ + 1);

		el_rel_ = (const pair<int, int>*) buffer;
		buffer += sizeof(pair<int, int>) * el_rel_offset_[el_num_];
	}

	{
		vl_rel_offset_ = (const int*) buffer; 
		assert(vl_rel_offset_[0] == 0);
		buffer += sizeof(int) * (vl_num_ + 1);

		vl_rel_ = (const int*) buffer;
		buffer += sizeof(int) * vl_rel_offset_[vl_num_];
	}
	assert((buffer - orig) == encode_size);
    std::cout << "~DataGraph::ReadBinary" << fname << "\n";
}

void DataGraph::ClearRawData() {
	raw_.vl_cnt_.clear();
	raw_.el_cnt_.clear();
	raw_.vlabels_.clear();
	raw_.out_edges_.clear();
	raw_.in_edges_.clear();
	raw_.offset_.clear();
	raw_.label_.clear();
	raw_.adj_offset_.clear();
	raw_.adj_.clear();
	raw_.in_offset_.clear();
	raw_.in_label_.clear();
	raw_.in_adj_offset_.clear();
	raw_.in_adj_.clear();
	raw_.vl_offset_.clear();
	raw_.vl_.clear();
	raw_.vl_rel_.clear();
	raw_.el_rel_.clear();
}

int DataGraph::GetNumVertices() {
	return vnum_;
}

int DataGraph::GetNumVertices(int vl) {
	return vl_cnt_[vl];
}

int DataGraph::GetNumEdges() {
	return enum_;
}

int DataGraph::GetNumEdges(int el) {
	return el_cnt_[el];
}

int DataGraph::GetNumVLabels(int v) {
	if (v == -1)
		return vl_num_; 
	else {
		range r = GetVLabels(v);
		return r.end - r.begin;
	}
}

int DataGraph::GetNumELabels(int v, bool dir) {
	if (v == -1)
		return el_num_; 
	else {
		range r = GetELabels(v, dir);
		return r.end - r.begin;
	}
}

range DataGraph::GetVLabels(int v) {
	range r;
	r.begin = vl_ + vl_offset_[v];
	r.end   = vl_ + vl_offset_[v+1];
	return r;
}

bool DataGraph::HasVLabel(int v, int vl) {
	int begin = vl_offset_[v];
	int end   = vl_offset_[v+1];

	int res = search(vl_, begin, end, vl);
	return res != -1;
}

range DataGraph::GetELabels(int v, bool dir = true) {
	const int* offset = dir ? offset_ : in_offset_;
	const int* label  = dir ? label_  : in_label_;

	range r;
	r.begin = label + offset[v]; 
	r.end   = label + offset[v+1]; 
	return r;
}

bool DataGraph::HasELabel(int v, int el, bool dir = true) {
	range r = GetAdj(v, el, dir);
	return r.end != r.begin;
}

int DataGraph::GetELabelIndex(int v, int el, bool dir = true) {
	const int* offset = dir ? offset_ : in_offset_;
	const int* label  = dir ? label_ : in_label_; 
	const int* adj_o  = dir ? adj_offset_ : in_adj_offset_; 
	const int* adj    = dir ? adj_   : in_adj_;

	int begin = offset[v];
	int end   = offset[v+1];
	int res = search(label, begin, end, el);
	return res == -1 ? -1 : res - begin;
}

range DataGraph::GetAdj(int v, int el, bool dir = true) {
	const int* offset = dir ? offset_ : in_offset_;
	//const pair<int, int>* label  = dir ? label_ : in_label_; 
	const int* label  = dir ? label_ : in_label_; 
	const int* adj_o  = dir ? adj_offset_ : in_adj_offset_; 
	const int* adj    = dir ? adj_   : in_adj_;

	range r;
	r.begin = r.end = adj;

	int begin = offset[v];
	int end   = offset[v+1];
	int res = search(label, begin, end, el);
	if (res == -1) return r;
	r.begin = adj + adj_o[res];
	r.end   = adj + adj_o[res+1];
	return r;

	/*int beginL = offset[v];
	  int endL   = offset[v+1];
	  int beginB = 0;
	  int endB   = 0;
	  if (endL - beginL > BINARY_THRESHOLD) {
	  int mid;
	  while (beginL < endL) {
	  mid = (beginL + endL) / 2;
	  if (label[mid] < el)
	  beginL = mid + 1;
	  else if (label[mid] > el)
	  endL = mid;
	  else {
	  beginB = adj_o[mid]; 
	  endB   = adj_o[mid+1]; 
	  break;
	  }
	  }
	  } else {
	  for (int i = beginL; i < endL; i++) {
	  if (label[i] == el) {
	  beginB = adj_o[i]; 
	  endB   = adj_o[i+1]; 
	  break;
	  }
	  }
	  }

	  range r;
	  r.begin = adj + beginB;
	  r.end   = adj + endB;
	  return r;*/
}

int DataGraph::GetAdjSize(int v, int el, bool dir = true) {
	range r = GetAdj(v, el, dir);
	return r.end - r.begin;
}

bool DataGraph::HasEdge(int u, int v, int el, bool dir = true) {
	int u_adj_size = GetAdjSize(u, el, dir);
	int v_adj_size = GetAdjSize(v, el, !dir);

	range r = u_adj_size < v_adj_size ? GetAdj(u, el, dir) : GetAdj(v, el, !dir); 
	int target = u_adj_size < v_adj_size ? v : u;

	int s = 0; 
	int e = r.end - r.begin;

	return search(r.begin, s, e, target) != -1;

	/*if (e - s > BINARY_THRESHOLD) {
	  int mid;
	  while (s < e) {
	  mid = (s + e) / 2;
	  if (r.begin[mid] < target)
	  s = mid + 1;
	  else if (r.begin[mid] > target)
	  e = mid;
	  else {
	  return true;
	  }
	  }
	  } else {
	  for (int i = s; i < e; i++) {
	  if (r.begin[i] == target) {
	  return true;
	  }
	  }
	  }
	  return false;*/
}

vector<int> DataGraph::GetRandomEdge(int el) {
	vector<int> ret;
  //std::cout << "GetRandomEdge(" << el << ")\n";
	int begin = el_rel_offset_[el]; 
	int end   = el_rel_offset_[el+1]; 

	if (begin == end)
		return ret;
	int r = rand();
  //std::cout << "Random num: " << r << "\n";	
  	r %= (end - begin);
  //std::cout << "Random E among " << (end - begin) << " -> " << r << "\n";
	r += begin;
	ret.resize(2);
	ret[0] = el_rel_[r].first;
	ret[1] = el_rel_[r].second;
  //int r = (rand() % (end - begin)) + begin; 
	//ret.push_back(el_rel_[r].first);
	//ret.push_back(el_rel_[r].second);
	return ret;
}

vector<int> DataGraph::GetEdge(int el, int i) { // XXX
    vector<int> ret;
    int begin = el_rel_offset_[el]; 
    int end   = el_rel_offset_[el+1]; 

    if (begin == end)
        return ret;
    assert(i < end - begin);
    int r = i + begin;
    ret.resize(2);
    ret[0] = el_rel_[r].first;
    ret[1] = el_rel_[r].second;
    //ret.push_back(el_rel_[r].first);
    //ret.push_back(el_rel_[r].second);
    return ret;
}

vector<int> DataGraph::GetRandomEdge(int v, int el, bool dir) {
    vector<int> ret;
  //std::cout << "GetRandomEdge(" << el << "," << dir << "," << v << ")\n";
	range r = GetAdj(v, el, dir);
	if (r.begin == r.end)
		return ret;
	int rv = rand();
  //std::cout << "Random num: " << rv << "\n";
 	rv %= (r.end - r.begin);
  //std::cout << "Random E among " << (r.end - r.begin) << " -> " << rv << "\n";
	int other = r.begin[rv];
	ret.resize(2);
	ret[0] = dir ? v : other;
	ret[1] = dir ? other : v;
	//int other = r.begin[rand() % (r.end - r.begin)];
	//ret.push_back(dir ? v : other);
	//ret.push_back(dir ? other : v); 
	return ret;
}

vector<int> DataGraph::GetRandomVertex(int vl) {
	vector<int> ret;

	int begin = vl_rel_offset_[vl]; 
	int end   = vl_rel_offset_[vl+1]; 

	if (begin == end)
		return ret;
	int r = rand();
  //std::cout << "Random num: " << r << "\n";
	r %= (end - begin);
  //std::cout << "Random V among " << (end - begin) << " -> " << r << "\n";
	r += begin; 
	ret.push_back(vl_rel_[r]);
	return ret;
}

vector<int> DataGraph::GetVertex(int vl, int i) {
    vector<int> ret;

    int begin = vl_rel_offset_[vl]; 
    int end   = vl_rel_offset_[vl+1]; 

    if (begin == end)
        return ret;
    assert(i < end - begin);
    int r = i + begin; 
    ret.push_back(vl_rel_[r]);
    return ret;
}

/*range DataGraph::GetRel(int el, bool dir = true) {
  const int* offset = dir ? rel_offset_ : in_rel_offset_;
  const int* re     = dir ? rel_        : in_rel_;

  range r;
  r.begin = re + offset[el];
  r.end   = re + offset[el+1];
  return r;
  }

  int DataGraph::GetRelSize(int el, bool dir = true) {
  const int* offset = dir ? rel_offset_ : in_rel_offset_;
  const int* re     = dir ? rel_        : in_rel_;

  range r;
  r.begin = re + offset[el];
  r.end   = re + offset[el+1];
  return r.end - r.begin;
  }

  range DataGraph::GetUni(int vl) {
  range r;
  r.begin  = vl_rel_ + vl_offset_[vl];
  r.end    = vl_rel_ + vl_offset_[vl+1];
  return r;
  }

  int DataGraph::GetUniSize(int vl) {
  range r = GetUni(vl);
  return r.end - r.begin;
  }*/
