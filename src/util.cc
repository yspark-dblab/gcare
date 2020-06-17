#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include "../include/util.h"

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

