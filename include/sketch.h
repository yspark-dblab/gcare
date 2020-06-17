#ifndef SKETCH_H_
#define SKETCH_H_

#include <map>
#include <experimental/filesystem>
#include <set>
#include <math.h>
#include "util.h"
#include "data_relations.h"
#include "query_relations.h"

class Sketch {
public:
	int t; //data table id
	int active_col; //0 for src, 1 for dst
	vector<int> join_cols; //0 for src, 1 for dst, vector of partitioned (hashed) attributes
	vector<int> hash_sizes; //same length with join_cols, vector of hash sizes
	vector<int> bound_cols; //bounded attributes
	vector<int> bounds; //bounded data vertices
	DataGraph* g;
	map<int, map<int, int>> l2gIndex;
	map<int, map<int, int>> g2lIndex;

	Sketch(int _t, int _active_col, vector<int>& _join_cols, vector<int>& _hash_sizes, vector<int>& _bounds, vector<int>& _bound_cols, DataGraph* _g)
		: t(_t), active_col(_active_col), join_cols(_join_cols), hash_sizes(_hash_sizes), bounds(_bounds), bound_cols(_bound_cols), g(_g) {}
	
	virtual int access(int boundID, int gVarIndex, vector<int>& arr) = 0;
};


class ZeroDimensionalSketchUnc : public Sketch {
public:
	int unc[1];

	ZeroDimensionalSketchUnc(int _t, int _active_col, vector<int>& _join_cols, vector<int>& _hash_sizes, vector<int>& _bounds, vector<int>& _bound_cols, DataGraph* _g, string file) : 
		Sketch(_t, _active_col, _join_cols, _hash_sizes, _bounds, _bound_cols, _g)
	{
		if (file.length() == 0) {
			unc[0] = 0;
			populate();
		}
		else {
			deserialize(file);
		}
	}

	void deserialize(string file) {
		ifstream in(file.c_str());
		string line;
		while (getline(in, line))
			unc[0] = stoi(line);
		in.close();
	}

	void populate() {
		if (g->table_[t].size() == 0)
			return;

		if (bounds.size() == 0) {
			unc[0] = g->table_[t].size();
			return;
		}
		
		for (int i = 0; i < g->table_[t].size(); i++) {
			bool match = true;
			for (int j = 0; j < bounds.size(); j++) {
				if (bounds[j] != g->table_[t][i][bound_cols[j]])
					match = false;
			}
			if (match)
				unc[0]++;
		}
	}

	int access(int boundID, int gVarIndex, vector<int>& arr) {
		return unc[0];
	}
};

class ZeroDimensionalSketchCon : public Sketch {
public:
	map<int, vector<int>> con;

	ZeroDimensionalSketchCon(int _t, int _active_col, vector<int>& _join_cols, vector<int>& _hash_sizes, vector<int>& _bounds, vector<int>& _bound_cols, DataGraph* _g, string file) : 
		Sketch(_t, _active_col, _join_cols, _hash_sizes, _bounds, _bound_cols, _g) 
	{
		assert(hash_sizes.size() == 1);
		assert(hash_sizes[0] == 1);
		
		if (file.length() == 0) {
			con[0].resize(hash_sizes[0]);
			populate();
		}
		else {
			deserialize(file);
		}
	}

	void deserialize(string file) {
		ifstream in(file.c_str());
		string line;
		while (getline(in, line))
			con[0].push_back(stoi(line));
		assert(con[0].size() == hash_sizes[0]);
		in.close();
	}

	void populate() {
		if (g->table_[t].size() == 0)
			return;

		map<int, int> cnt;
		for (int i = 0; i < g->table_[t].size(); i++) {
			bool match = true;
			for (int j = 0; j < bounds.size(); j++) {
				if (bounds[j] != g->table_[t][i][bound_cols[j]])
					match = false;
			}
			if (match) {
				int v = g->table_[t][i][active_col];
				cnt[v]++;
			}
		}
		for (auto p : cnt) {
			if (p.second > con[0][0])
				con[0][0] = p.second;
		}
	}

	int access(int boundID, int gVarIndex, vector<int>& arr) {
		assert(arr[gVarIndex] == 0);

		return con[g2lIndex[boundID][gVarIndex]][arr[l2gIndex[boundID][0]]];
	}
};

class OneDimensionalSketchUnc : public Sketch {
public:
	vector<int> unc;
	
	OneDimensionalSketchUnc(int _t, int _active_col, vector<int>& _join_cols, vector<int>& _hash_sizes, vector<int>& _bounds, vector<int>& _bound_cols, DataGraph* _g, string file) : 
		Sketch(_t, _active_col, _join_cols, _hash_sizes, _bounds, _bound_cols, _g) 
	{
		if (file.length() == 0) {
			unc.resize(hash_sizes[0], 0);
			populate();
		}
		else {
			deserialize(file);
		}
	}

	void deserialize(string file) {
		unc.clear();
		ifstream in(file.c_str());
		string line;
		while (getline(in, line))
			unc.push_back(stoi(line));
		assert(unc.size() == hash_sizes[0]);
		in.close();
	}

	void populate() {
		assert(join_cols.size() == 1);
		assert(hash_sizes.size() == 1);

		if (g->table_[t].size() == 0)
			return;

		for (int i = 0; i < g->table_[t].size(); i++) {
			bool match = true;
			for (int j = 0; j < bounds.size(); j++) {
				if (bounds[j] != g->table_[t][i][bound_cols[j]])
					match = false;
			}
			if (match) {
				int h = g->table_[t][i][join_cols[0]] % hash_sizes[0];
				unc[h]++;
			}
		}
	}

	int access(int boundID, int gVarIndex, vector<int>& arr) {
		if (l2gIndex[boundID].empty())
			return unc[0];
		else {
			assert(gVarIndex == -1);
			assert(l2gIndex[boundID][0] < arr.size());
			assert(arr[l2gIndex[boundID][0]] < hash_sizes[0]);

			return unc[arr[l2gIndex[boundID][0]]];
		}
	}
};

class OneDimensionalSketchCon : public Sketch {
public:
	map<int, vector<int>> con;
	
	OneDimensionalSketchCon(int _t, int _active_col, vector<int>& _join_cols, vector<int>& _hash_sizes, vector<int>& _bounds, vector<int>& _bound_cols, DataGraph* _g, string file) : 
		Sketch(_t, _active_col, _join_cols, _hash_sizes, _bounds, _bound_cols, _g) 
	{
		if (file.length() == 0) {
			con[0].resize(hash_sizes[0]);
			populate();
		}
		else {
			deserialize(file);
		}
	}

	void deserialize(string file) {
		ifstream in(file.c_str());
		string line;
		while (getline(in, line))
			con[0].push_back(stoi(line));
		assert(con[0].size() == hash_sizes[0]);
		in.close();
	}

	void populate() {
		assert(join_cols.size() == 1);
		assert(hash_sizes.size() == 1);

		if (g->table_[t].size() == 0)
			return;

		vector<map<int, int>> cnt(hash_sizes[0]);
		for (int i = 0; i < g->table_[t].size(); i++) {
			bool match = true;
			for (int j = 0; j < bounds.size(); j++) {
				if (bounds[j] != g->table_[t][i][bound_cols[j]])
					match = false;
			}
			if (match) {
				int v = g->table_[t][i][active_col];
				int h = g->table_[t][i][join_cols[0]] % hash_sizes[0];
				cnt[h][v]++;
			}
		}
		for (int h = 0; h < hash_sizes[0]; h++) {
			for (auto p : cnt[h]) {
				if (p.second > con[0][h])
					con[0][h] = p.second;
			}
		}
	}

	int access(int boundID, int gVarIndex, vector<int>& arr) {
		int res = con[g2lIndex[boundID][gVarIndex]][arr[l2gIndex[boundID][0]]];
		return res; 
	}
};

class TwoDimensionalSketchUnc : public Sketch {
public:
	vector<vector<int>> unc;
	
	TwoDimensionalSketchUnc(int _t, int _active_col, vector<int>& _join_cols, vector<int>& _hash_sizes, vector<int>& _bounds, vector<int>& _bound_cols, DataGraph* _g, string file) : 
		Sketch(_t, _active_col, _join_cols, _hash_sizes, _bounds, _bound_cols, _g) 
	{
		unc.resize(hash_sizes[0], vector<int>(hash_sizes[1], 0));
		if (file.length() == 0) {
			populate();
		}
		else {
			deserialize(file);
		}
	}

	void deserialize(string file) {
		ifstream in(file.c_str());
		string line;
		int i = 0;
		while (getline(in, line)) {
			unc[i / hash_sizes[1]][i % hash_sizes[1]] = stoi(line);
			i++;
		}
		assert(i == hash_sizes[0] * hash_sizes[1]);
		in.close();
	}

	void populate() {
		assert(join_cols.size() == 2);
		assert(hash_sizes.size() == 2);

		if (g->table_[t].size() == 0)
			return;

		for (int i = 0; i < g->table_[t].size(); i++) {
			bool match = true;
			for (int j = 0; j < bounds.size(); j++) {
				if (bounds[j] != g->table_[t][i][bound_cols[j]])
					match = false;
			}
			if (match) {
				int h0 = g->table_[t][i][join_cols[0]] % hash_sizes[0];
				int h1 = g->table_[t][i][join_cols[1]] % hash_sizes[1];
				unc[h0][h1]++;
			}
		}
	}

	int access(int boundID, int gVarIndex, vector<int>& arr) {
		assert(gVarIndex == -1);
		assert(l2gIndex[boundID][0] < arr.size());
		assert(l2gIndex[boundID][1] < arr.size());
		assert(arr[l2gIndex[boundID][0]] < hash_sizes[0]);
		assert(arr[l2gIndex[boundID][1]] < hash_sizes[1]);

		return unc
			[arr[l2gIndex[boundID][0]]]
			[arr[l2gIndex[boundID][1]]];
	}
};

class TwoDimensionalSketchCon : public Sketch {
public:
	map<int, vector<vector<int>>> con;
	
	TwoDimensionalSketchCon(int _t, int _active_col, vector<int>& _join_cols, vector<int>& _hash_sizes, vector<int>& _bounds, vector<int>& _bound_cols, DataGraph* _g, string file) : 
		Sketch(_t, _active_col, _join_cols, _hash_sizes, _bounds, _bound_cols, _g) 
	{
		con[0].resize(hash_sizes[0], vector<int>(hash_sizes[1], 0));
		con[1].resize(hash_sizes[0], vector<int>(hash_sizes[1], 0));
		if (file.length() == 0) {
			populate();
		}
		else {
			deserialize(file);
		}
	}

	void deserialize(string file) {
		ifstream in(file.c_str());
		string line;
		int i = 0;
		while (getline(in, line)) {
			con[0][i / hash_sizes[1]][i % hash_sizes[1]] = stoi(line);
			i++;
		}
		assert(i == hash_sizes[0] * hash_sizes[1]);
		in.close();
	}

	void populate() {
		assert(join_cols.size() == 2);
		assert(hash_sizes.size() == 2);

		if (g->table_[t].size() == 0)
			return;

		vector<vector<map<int, int>>> cnt(hash_sizes[0], vector<map<int, int>>(hash_sizes[1]));
		for (int i = 0; i < g->table_[t].size(); i++) {
			bool match = true;
			for (int j = 0; j < bounds.size(); j++) {
				if (bounds[j] != g->table_[t][i][bound_cols[j]])
					match = false;
			}
			if (match) {
				int v = g->table_[t][i][active_col];
				int h0 = g->table_[t][i][join_cols[0]] % hash_sizes[0];
				int h1 = g->table_[t][i][join_cols[1]] % hash_sizes[1];
				cnt[h0][h1][v]++;
			}
		}
		for (int h0 = 0; h0 < hash_sizes[0]; h0++) {
			for (int h1 = 0; h1 < hash_sizes[1]; h1++) {
				for (auto p : cnt[h0][h1]) {
					if (p.second > con[0][h0][h1])
						con[0][h0][h1] = p.second;
				}
			}
		}
	}

	int access(int boundID, int gVarIndex, vector<int>& arr) {
		assert(g2lIndex[boundID][gVarIndex] == 0 || g2lIndex[boundID][gVarIndex] == 1);
		assert(l2gIndex[boundID][0] < arr.size());
		assert(l2gIndex[boundID][1] < arr.size());
		assert(arr[l2gIndex[boundID][0]] < hash_sizes[0]);
		assert(arr[l2gIndex[boundID][1]] < hash_sizes[1]);

		return con[g2lIndex[boundID][gVarIndex]]
			[arr[l2gIndex[boundID][0]]]
			[arr[l2gIndex[boundID][1]]];
	}
};

class OfflineSketch {
public:
	int t; //data table id
	int buckets;
	DataGraph* g;

	set<int> hash_sizes; //candidate hash values (larger than 1)
	vector<int> h_sizes;
	int num_hash;
	vector<vector<int>> unc_2D;
	vector<vector<vector<int>>> con_2D; //active col

	OfflineSketch(int _t, int _buckets, DataGraph* _g) : t(_t), buckets(_buckets), g(_g) {
		int h = buckets;
		for (int p = 1; ; p++) {
			h = round(pow(buckets, 1.0 / p));
			if (h < 1)
				h = 1;
			hash_sizes.insert(h);
			if (h == 1)
				break;
		}
		h_sizes.clear();
		for (int h : hash_sizes)
			h_sizes.push_back(h);
		num_hash = h_sizes.size();

		populate();
	}

	void populate() {
		
		//edge label table
		if (t < g->base_) {

			unc_2D.resize(num_hash * num_hash); //partition on both columns
			con_2D.resize(2, vector<vector<int>>(num_hash * num_hash));

			for (int a = 0; a < num_hash; a++) {
				for (int b = 0; b < num_hash; b++) {
					int hs0 = h_sizes[a];
					int hs1 = h_sizes[b];

					if ((hs0-1) * (hs1-1) > buckets)
						continue;

					vector<int>& unc  = unc_2D[a * num_hash + b];
					vector<int>& con0 = con_2D[0][a * num_hash + b];
					vector<int>& con1 = con_2D[1][a * num_hash + b];
					unc.resize(hs0 * hs1, 0);
					con0.resize(hs0 * hs1, 0);
					con1.resize(hs0 * hs1, 0);

					vector<map<int, int>> cnt0;
					vector<map<int, int>> cnt1;
					cnt0.resize(hs0 * hs1);
					cnt1.resize(hs0 * hs1);

					for (int i = 0; i < g->table_[t].size(); i++) {
						int v0 = g->table_[t][i][0];
						int v1 = g->table_[t][i][1];
						int h0 = v0 % hs0; 
						int h1 = v1 % hs1; 
						cnt0[h0 * hs1 + h1][v0]++;
						cnt1[h0 * hs1 + h1][v1]++;
						 unc[h0 * hs1 + h1]++;
					}
					for (int h0 = 0; h0 < hs0; h0++) {
						for (int h1 = 0; h1 < hs1; h1++) {
							for (auto p : cnt0[h0 * hs1 + h1]) {
								if (p.second > con0[h0 * hs1 + h1])
									con0[h0 * hs1 + h1] = p.second;
							}
							for (auto p : cnt1[h0 * hs1 + h1]) {
								if (p.second > con1[h0 * hs1 + h1])
									con1[h0 * hs1 + h1] = p.second;
							}
						}
					}
				}
			}
		}
		//vertex label table
		else {

			unc_2D.resize(num_hash); //partition on both columns
			con_2D.resize(1, vector<vector<int>>(num_hash * num_hash));

			for (int a = 0; a < num_hash; a++) {
				int hs = h_sizes[a];

				vector<int>& unc  = unc_2D[a];
				vector<int>& con  = con_2D[0][a];
				unc.resize(hs, 0);
				con.resize(hs, 0);

				vector<map<int, int>> cnt(hs);

				for (int i = 0; i < g->table_[t].size(); i++) {
					int v = g->table_[t][i][0];
					int h = v % hs; 
					cnt[h][v]++;
					unc[h]++;
				}
				for (int h = 0; h < hs; h++) {
					for (auto p : cnt[h]) {
						if (p.second > con[h])
							con[h] = p.second;
					}
				}
			}
		}
	}
	
	//write sketches into different files 
	void serialize(const char* dir) {
		string fn(dir); 
		fn += "/";

		if (t < g->base_) {
			for (int a = 0; a < num_hash; a++) {
				for (int b = 0; b < num_hash; b++) {
					int hs0 = h_sizes[a];
					int hs1 = h_sizes[b];
					
					if ((hs0-1) * (hs1-1) > buckets)
						continue;

					//unc
					string name(fn);
					if (hs0 > 1 && hs1 > 1)
						name += "2d_" + to_string(t) + "_-1_" + to_string(hs0) + "_" + to_string(hs1) + "_0_1_.txt";
					else if (hs0 > 1 && hs1 == 1)
						name += "1d_" + to_string(t) + "_-1_" + to_string(hs0) + "_0_.txt";
					else if (hs0 == 1 && hs1 > 1)
						name += "1d_" + to_string(t) + "_-1_" + to_string(hs1) + "_1_.txt";
					else
						name += "0d_" + to_string(t) + "_-1_.txt";
					//cout << "serializing a sketch to " << name << endl;
					ofstream out(name.c_str());

					vector<int>& unc  = unc_2D[a * num_hash + b];
					for (int h0 = 0; h0 < hs0; h0++) { 
						for (int h1 = 0; h1 < hs1; h1++) {
							out << unc[h0 * hs1 + h1] << endl;
						}
					}
					out.close();

					//con
					for (int active_col = 0; active_col <= 1; active_col++) {
						string name(fn);
						if (hs0 > 1 && hs1 > 1)
							name += "2d_" + to_string(t) + "_" + to_string(active_col) + "_" + to_string(hs0) + "_" + to_string(hs1) + "_0_1_.txt";
						else if (hs0 > 1 && hs1 == 1)
							name += "1d_" + to_string(t) + "_" + to_string(active_col) + "_" + to_string(hs0) + "_0_.txt";
						else if (hs0 == 1 && hs1 > 1)
							name += "1d_" + to_string(t) + "_" + to_string(active_col) + "_" + to_string(hs1) + "_1_.txt";
						else
							name += "0d_" + to_string(t) + "_" + to_string(active_col) + "_.txt";
						//cout << "serializing a sketch to " << name << endl;
						ofstream out(name.c_str());

						vector<int>& con  = con_2D[active_col][a * num_hash + b];
						for (int h0 = 0; h0 < hs0; h0++) {
							for (int h1 = 0; h1 < hs1; h1++) {
								out << con[h0 * hs1 + h1] << endl;
							}
						}
						out.close();
					}
				}
			}
		}
		else {
			for (int a = 0; a < num_hash; a++) {
				int hs = h_sizes[a];

				//unc
				string name(fn);
				if (hs > 1)
					name += "1d_" + to_string(t) + "_-1_" + to_string(hs) + "_0_.txt";
				else 
					name += "0d_" + to_string(t) + "_-1_.txt";
				//cout << "serializing a sketch to " << name << endl;
				ofstream out(name.c_str());

				vector<int>& unc  = unc_2D[a];
				for (int h = 0; h < hs; h++) { 
					out << unc[h] << endl;
				}
				out.close();

				//con
				for (int active_col = 0; active_col <= 0; active_col++) {
					string name(fn);
					if (hs > 1)
						name += "1d_" + to_string(t) + "_" + to_string(active_col) + "_" + to_string(hs) + "_0_.txt";
					else
						name += "0d_" + to_string(t) + "_" + to_string(active_col) + "_.txt";
					//cout << "serializing a sketch to " << name << endl;
					ofstream out(name.c_str());

					vector<int>& con  = con_2D[active_col][a];
					for (int h = 0; h < hs; h++) {
						out << con[h] << endl;
					}
					out.close();
				}
			}
		}
	}

	static void deserialize(const char* dir, unordered_map<string, Sketch*>& sketch_map, DataGraph* g) {
		namespace fs = std::experimental::filesystem;
		for (auto& dir_entry : fs::recursive_directory_iterator(dir)) {
			Sketch* s;
			string last = dir_entry.path().filename();
			string name = dir_entry.path().string(); //full path
			//cout << "deserializing from " << name << endl;
			auto tokens = tokenize(last, "_");
			int t = stoi(tokens[1]);
			int active_col = stoi(tokens[2]);
			assert(active_col >= -1 && active_col <= 1);
			vector<int> join_cols, hash_sizes, bounds, bound_cols;

			if (tokens[0][0] == '0') {
				assert(tokens.size() == 4);
				hash_sizes.push_back(1);
				if (active_col == -1)
					s = new ZeroDimensionalSketchUnc(t, active_col, join_cols, hash_sizes, bounds, bound_cols, g, name); 
				else
					s = new ZeroDimensionalSketchCon(t, active_col, join_cols, hash_sizes, bounds, bound_cols, g, name); 
			} else if (tokens[0][0] == '1') {
				assert(tokens.size() == 6);
				hash_sizes.push_back(stoi(tokens[3]));
				join_cols.push_back(stoi(tokens[4]));
				if (active_col == -1)
					s = new OneDimensionalSketchUnc(t, active_col, join_cols, hash_sizes, bounds, bound_cols, g, name); 
				else
					s = new OneDimensionalSketchCon(t, active_col, join_cols, hash_sizes, bounds, bound_cols, g, name); 
			} else if (tokens[0][0] == '2') {
				assert(tokens.size() == 8);
				hash_sizes.push_back(stoi(tokens[3]));
				hash_sizes.push_back(stoi(tokens[4]));
				join_cols.push_back(stoi(tokens[5]));
				join_cols.push_back(stoi(tokens[6]));
				if (active_col == -1)
					s = new TwoDimensionalSketchUnc(t, active_col, join_cols, hash_sizes, bounds, bound_cols, g, name); 
				else
					s = new TwoDimensionalSketchCon(t, active_col, join_cols, hash_sizes, bounds, bound_cols, g, name); 
			} else {
				cout << "is " << name << " a sketch?" << endl;
				continue;
			}

			string probe;
			probe.append(to_string(t));
			probe.append("[");
			probe.append(to_string(active_col));
			probe.append("][");
			for (int h : hash_sizes) {
				probe.append(to_string(h));
				probe.append(", ");
			}
			probe.append("][");
			for (int c : join_cols) {
				probe.append(to_string(c));
				probe.append(", ");
			}
			probe.append("][");
			probe.append("]");

			//cout << "deserializing to " << probe << endl;

			sketch_map[probe] = s;
		}
	}
};


#endif
