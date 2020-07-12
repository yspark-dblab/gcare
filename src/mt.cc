#include "../include/mt.h"
#include <boost/algorithm/string.hpp>
#include "../include/util.h"
#include <queue>
#include <limits>

void MarkovTable::Init() {
    getSubstructureFlag = true;
    mt_.clear();
    ceg.clear();
    largestMTEntries.clear();
}

void MarkovTable::PrepareSummaryStructure(DataGraph& g, double ratio) {
}

void MarkovTable::WriteSummary(const char* fn) {
}

void MarkovTable::ReadSummary(const char* fn) {
    string line;
    ifstream catalogueFile(fn);
    if (catalogueFile.is_open()) {
        while (getline(catalogueFile, line)) {
            vector<string> entry;
            boost::split(entry, line, boost::is_any_of(","));
            if (mt_.count(entry[1]) == 0) {
                mt_.insert(pair<string, map<string, long>>(entry[1], map<string, long>()));
            }
            mt_[entry[1]].insert(pair<string, long>(entry[2], stol(entry[3])));
        }
        catalogueFile.close();
    }
}

int MarkovTable::DecomposeQuery() {
    return 1;
}

bool MarkovTable::GetSubstructure(int subquery_index) {
    if (getSubstructureFlag) {
        getSubstructureFlag = false;
        return true;
    }
    return false;
}

double MarkovTable::EstCard(int subquery_index) {
    pair<string, string> vListAndLabelSeq = q->toVListAndLabelSeq();
    const string &queryVList = vListAndLabelSeq.first;
    const string &queryLabelSeq = vListAndLabelSeq.second;
    decompose(vListAndLabelSeq.first, 2);

    map<string, double> estimates;
    queue<string> queue;
    for (const string &startVList : largestMTEntries) {
        queue.push(startVList);
        string labelSeq = extractLabelSeq(startVList, queryVList, queryLabelSeq);
        estimates.insert(pair<string, double>(startVList, mt_[startVList][labelSeq]));
    }

    // perform BFS
    set<string> inQueue;
    while (!queue.empty()) {
        string currentVList = queue.front();
        queue.pop();

        for (const string &nextVList : ceg[currentVList]) {
            set<pair<string, string>> extensions = getExtensions(currentVList, nextVList);
            estimates.insert(pair<string, double>(nextVList, getMaxExt(extensions, queryVList, queryLabelSeq)));
            if (inQueue.count(nextVList)) continue;
            queue.push(nextVList);
            inQueue.insert(nextVList);
        }
    }

    return estimates[queryVList];
}

double MarkovTable::AggCard() {
    return card_vec_[0];
}

double MarkovTable::GetSelectivity() {
    return 1;
}

void MarkovTable::getDecom(const vector<string> &vListEdges, const int &mtLen, int depth, const string &current, const string &parent) {
    if (!parent.empty() || depth == mtLen) {
        string sortedCurrent = sortVList(current);
        if (depth == mtLen) {
            largestMTEntries.insert(sortedCurrent);
        }
        if (!parent.empty()) {
            ceg[sortVList(parent)].insert(sortedCurrent);
            if (depth == vListEdges.size()) {
                ceg.insert(pair<string, set<string>>(sortedCurrent, set<string>()));
                return;
            }
        }
    }

    for (int i = 0; i < vListEdges.size(); i++) {
        string next = vListEdges[i];
        string updated;
        if (current.empty()) {
            updated = next;
        } else {
            updated.append(current);
            updated.append(";");
            updated.append(next);
        }
        if (depth > 0 && !isAcyclicConnected(updated)) continue;
        getDecom(vListEdges, mtLen, depth + 1, updated, current);
    }
}

void MarkovTable::decompose(const string &vListString, int mtLen) {
    vector<string> vListEdges;
    boost::split(vListEdges, vListString, boost::is_any_of(";"));
    getDecom(vListEdges, mtLen, 0, "", "");
}

set<pair<string, string>> MarkovTable::getExtensions(const string &currentVListStr, const string &nextVListStr) {
    vector<string> currentVList, nextVList;
    boost::split(currentVList, currentVListStr, boost::is_any_of(";"));
    boost::split(nextVList, nextVListStr, boost::is_any_of(";"));

    set<string> currentEdges;
    string ext;

    int j = 0;
    for (int i = 0; i < nextVList.size(); ++i) {
        if (nextVList[i] == currentVList[j]) {
            currentEdges.insert(currentVList[j]);
            ++j;
        } else {
            ext = nextVList[i];
            break;
        }
    }

    set<pair<string, string>> numerAndDenom;
    string numerator;
    for (int i = 0; i < currentVList.size(); ++i) {
        if (i < j) {
            numerator = currentVList[i] + ";" + ext;
        } else {
            numerator = ext + ";" + currentVList[i];
        }

        if (isAcyclicConnected(numerator)) {
            numerAndDenom.insert(pair<string, string>(numerator, currentVList[i]));
        }
    }

    return numerAndDenom;
}

double MarkovTable::getMaxExt(const set<pair<string, string>> &extensions, const string &queryVList, const string &queryLabelSeq) {
    double maxExt = numeric_limits<double>::min();
    string numerLabelSeq, denomLabelSeq;
    for (const pair<string, string> &ext : extensions) {
        numerLabelSeq = extractLabelSeq(ext.first, queryVList, queryLabelSeq);
        denomLabelSeq = extractLabelSeq(ext.second, queryVList, queryLabelSeq);
        maxExt = max(maxExt, ((double) mt_[ext.first][numerLabelSeq]) / mt_[ext.second][denomLabelSeq]);
    }
    return maxExt;
}
