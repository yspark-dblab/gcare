#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include "../include/util.h"
#include <boost/algorithm/string.hpp>

vector<string> parse(string line, string del)
{
	vector<string> ret;

	size_t pos = 0;
	string token;
	while((pos = line.find(del)) != string::npos)
	{
		token = line.substr(0, pos);
		ret.push_back(token);
		line.erase(0, pos + del.length());
	}
	ret.push_back(line);
	return ret;
}
/*
int search(const int* arr, int begin, int end, int target) {
	if (end - begin > BINARY_THRESHOLD) {
		int mid;
		while (begin < end) {
			mid = (begin + end) / 2;
			int val = arr[mid];
			if (val < target)
				begin = mid + 1;
			else if (val > target)
				end = mid;
			else {
				return mid;
			}
		}
	} else {
		for (int i = begin; i < end; i++) {
			int val = arr[i];
			if (val < target) continue;
			return val == target ? i : -1;

//			if (arr[i] == target) {
//				return i; 
//			}
		}
	}
	return -1;
}
*/

vector<string> tokenize(const string& line, const char* delim)
{
	vector<string> tokens;
	char *dup = strdup(line.c_str());
	char *tok = strtok(dup, delim);
	while(tok != NULL)
	{
		tokens.push_back(string(tok));
		tok = strtok(NULL, delim);
	}
	free(dup);
	
	return tokens;
}

string sortVList(const string &vList) {
    vector<string> sorted;
    boost::split(sorted, vList, boost::is_any_of(";"));
    sort(sorted.begin(), sorted.end());
    string result;
    for (const string &e : sorted) {
        result.append(";").append(e);
    }
    return result.substr(1);
}

bool isAcyclicConnected(const string &vListStr) {
    set<string> component;
    vector<string> vList;
    boost::split(vList, vListStr, boost::is_any_of(";-"));
    for (int i = 0; i < vList.size(); i += 2) {
        if (!component.empty()) {
            bool hasSrc = component.count(vList[i]);
            bool hasDest = component.count(vList[i + 1]);
            if ((hasSrc && hasDest) || (!hasSrc && !hasDest)) {
                return false;
            }
        }

        component.insert(vList[i]);
        component.insert(vList[i + 1]);
    }
    return true;
}

string extractLabelSeq(const string &subVListStr, const string &queryVListStr, const string &queryLabelSeqStr) {
    vector<string> subVList, queryVList, queryLabelSeq;
    boost::split(subVList, subVListStr, boost::is_any_of(";"));
    boost::split(queryVList, queryVListStr, boost::is_any_of(";"));
    boost::split(queryLabelSeq, queryLabelSeqStr, boost::is_any_of("->"));

    string result;
    int i = 0;
    for (const string &edge : subVList) {
        for (; i < queryVList.size(); ++i) {
            if (queryVList[i] == edge) {
                result.append("->").append(queryLabelSeq[i]);
                break;
            }
        }
    }
    return result.substr(2);
}
