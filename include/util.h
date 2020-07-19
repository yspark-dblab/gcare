#ifndef UTIL_H_
#define UTIL_H_

#include <boost/functional/hash.hpp>
#include <string>
#include <vector>
#include <utility>

#define BINARY_THRESHOLD 4

using namespace std;

vector<string> parse(string line, string del);

inline int search(const int* arr, int begin, int end, int target) {
	if (end - begin > BINARY_THRESHOLD) {
		int mid;
		while (begin < end) {
			mid = (begin + end) / 2;
			int val = arr[mid];
			if (val < target) begin = mid + 1;
			else if (val > target) end = mid;
			else return mid;
		}
	} else {
		for (int i = begin; i < end; i++) {
			int val = arr[i];
			if (val < target) continue;
			return val == target ? i : -1;
		}
	}
	return -1;
}

vector<string> tokenize(const string& line, const char* delim);

struct range {
	const int* begin;
	const int* end;

	bool operator< (const range& other) const {
		int size = this->end - this->begin;
		int osize = other.end - other.begin;
		if (size == osize)
			return this->begin < other.begin;
		else
			return size < osize;
	}

	bool operator== (const range& other) const {
		return this->begin == other.begin && this->end == other.end; 
	}

	bool good() {
		for (int i = 0; i < (end-begin)-1; i++)
			if (begin[i] >= begin[i+1])
				return false;
		return true;
	}

	int size() {
		return end - begin;
	}
};

struct Edge {
	int src, dst, el, id;

	Edge(int s, int d, int e) : src(s), dst(d), el(e) {}
    Edge(int s, int d, int e, int id) : src(s), dst(d), el(e), id(id) {}
	Edge() {}

	static const int FORWARD = 0;
    static const int BACKWARD = 1;

	bool operator<(const Edge& other) const {
		if (src < other.src)
			return true;
		else if (src == other.src) {
			if (el < other.el)
				return true;
			else if (el == other.el) {
				return dst < other.dst;
			}
		}
		return false;
	}
	bool operator==(const Edge& other) const {
		return src == other.src && dst == other.dst && el == other.el;
	}

  pair<string, string> toVListAndLabelSeq() {
    return make_pair(to_string(src) + "-" + to_string(dst), to_string(el));
  }
};

namespace std {
    template<>
    struct hash<Edge> {
        size_t operator()(const Edge &e) const {
            return e.src + e.dst + e.el;
        }
    };
}

struct EdgeHasher {
	std::size_t operator () (const Edge &key) const 
	{
		// The following line is a stright forward implementation. But it can be
		// hard to design a good hash function if KeyData is complex.
		//return (key.id << 32 | key.age); // suppose size_t is 64-bit and int is 32-bit
		// A commonly used way is to use boost
		std::size_t seed = 0;
		boost::hash_combine(seed, boost::hash_value(key.src));
		boost::hash_combine(seed, boost::hash_value(key.dst));
		boost::hash_combine(seed, boost::hash_value(key.el));
		return seed;
	}
};

string sortVList(const string &vList);
bool isAcyclicConnected(const string &vList);
string extractLabelSeq(const string &subVList, const string &queryVList, const string &queryLabelSeq);

#endif
