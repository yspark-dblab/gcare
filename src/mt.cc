#include "../include/mt.h"
#include <boost/algorithm/string.hpp>
#include "../include/util.h"
#include <queue>
#include <limits>
#include <cmath>
#include <deque>

void MarkovTable::Init() {
    getSubstructureFlag = true;
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
            insertEntryToMT(entry);
        }
        catalogueFile.close();
    }
}

void MarkovTable::insertEntryToMT(const vector<string> &entry) {
    const long count = stol(entry[3]);
    vector<string> edges;
    boost::split(edges, entry[1], boost::is_any_of(";"));
    if (edges.size() == 1) {
        mt1_[stoi(entry[2])] = count;
    } else if (edges.size() == 2) {
        int vList[4];
        for (int i = 0; i < edges.size(); ++i) {
            vector<string> srcDest;
            boost::split(srcDest, edges[i], boost::is_any_of("-"));
            vList[i * 2] = stoi(srcDest[0]);
            vList[i * 2 + 1] = stoi(srcDest[1]);
        }

        vector<string> labelStrs;
        boost::split(labelStrs, entry[2], boost::is_any_of("->"));
        int labels[2] = {stoi(labelStrs[0]), stoi(labelStrs[1])};

        if (vList[0] == vList[2] || vList[0] == vList[3]) {
            if (vList[0] == vList[2]) {
                if (labels[0] < labels[1]) {
                    mt2_[Edge::FORWARD][Edge::FORWARD][labels[0]][labels[1]] = count;
                } else {
                    mt2_[Edge::FORWARD][Edge::FORWARD][labels[1]][labels[0]] = count;
                }
            } else {
                mt2_[Edge::FORWARD][Edge::BACKWARD][labels[0]][labels[1]] = count;
            }
        } else {
            if (vList[1] == vList[2]) {
                mt2_[Edge::FORWARD][Edge::BACKWARD][labels[1]][labels[0]] = count;
            } else {
                if (labels[0] < labels[1]) {
                    mt2_[Edge::BACKWARD][Edge::BACKWARD][labels[0]][labels[1]] = count;
                } else {
                    mt2_[Edge::BACKWARD][Edge::BACKWARD][labels[1]][labels[0]] = count;
                }
            }
        }
    }
}

void MarkovTable::OldReadSummary(const char *fn) {
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
    vector<tuple<int, int, Edge, Edge>> twoPaths;
    q->getAll2Paths(twoPaths);

    int subQEnc;
    double estimates[1 << q->GetNumEdges()];
    deque<vector<Edge>> queue;
    for (const tuple<int, int, Edge, Edge> &twoPath : twoPaths) {
        vector<Edge> starter;
        starter.emplace_back(get<2>(twoPath));
        starter.emplace_back(get<3>(twoPath));
        queue.emplace_back(starter);
        subQEnc = q->encodeSubQ(starter);
        estimates[subQEnc] = mt2_[get<0>(twoPath)][get<1>(twoPath)][get<2>(twoPath).el][get<3>(twoPath).el];
    }

    double extRate, currentEst;
    bool processed[1 << q->GetNumEdges()];
    while (!queue.empty()) {
        const vector<Edge> &current = queue.front();
        queue.pop_front();
        currentEst = estimates[q->encodeSubQ(current)];

        vector<tuple<int, int, Edge, Edge>> extensions;
        extensions.reserve(pow(q->GetNumEdges(), 2));
        getExtensions(extensions, current);
        for (const tuple<int, int, Edge, Edge> &ext : extensions) {
            extRate = calcExtRate(ext);
            vector<Edge> nextSubQ;
            nextSubQ.reserve(current.size() + 1);
            for (const Edge &e : current) {
                nextSubQ.push_back(e);
            }
            nextSubQ.push_back(get<3>(ext));
            subQEnc = q->encodeSubQ(nextSubQ);
            if (processed[subQEnc]) {
                estimates[subQEnc] = max(estimates[subQEnc], currentEst * extRate);
            } else {
                estimates[subQEnc] = currentEst * extRate;
                queue.emplace_back(nextSubQ);
            }
        }
    }

    return estimates[(1 << q->GetNumEdges()) - 1];
}

double MarkovTable::EstCardGreedyMax(int subquery_index) {
    pair<string, string> vListAndLabelSeq = q->toVListAndLabelSeq();
    const string &queryVList = vListAndLabelSeq.first;
    const string &queryLabelSeq = vListAndLabelSeq.second;
    decompose(vListAndLabelSeq.first, 2);

    map<string, double> estimates;
    queue<string> queue;

    // get first level node (size-2)
    pair<string, double> greedyMax = make_pair<>("", numeric_limits<double>::min());
    for (const string &startVList : largestMTEntries) {
        string labelSeq = extractLabelSeq(startVList, queryVList, queryLabelSeq);
        if (mt_[startVList][labelSeq] > greedyMax.second) {
            greedyMax.first = startVList;
            greedyMax.second = mt_[startVList][labelSeq];
        }
    }
    queue.push(greedyMax.first);
    estimates.insert(greedyMax);

    // extend
    while (!queue.empty()) {
        string currentVList = queue.front();
        queue.pop();

        greedyMax.first = "";
        greedyMax.second = numeric_limits<double>::min();
        for (const string &nextVList : ceg[currentVList]) {
            set<pair<string, string>> extensions = getExtensions(currentVList, nextVList);
            double nextEst = estimates[currentVList] * getMaxExt(extensions, queryVList, queryLabelSeq);
            if (nextEst > greedyMax.second) {
                greedyMax.first = nextVList;
                greedyMax.second = nextEst;
            }
        }

        if (!greedyMax.first.empty()) {
            queue.push(greedyMax.first);
            estimates.insert(greedyMax);
        }
    }

    return estimates[queryVList];
}

double MarkovTable::EstCardAllMax(int subquery_index) {
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
            double nextEst = estimates[currentVList] * getMaxExt(extensions, queryVList, queryLabelSeq);
            if (estimates.count(nextVList)) {
                estimates[nextVList] = max(estimates[nextVList], nextEst);
            } else {
                estimates[nextVList] = nextEst;
            }
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

void MarkovTable::getExtensions(vector<tuple<int, int, Edge, Edge>> &extensions, const vector<Edge> &current) {
    int minVId;
    for (const Edge &e : current) {
        minVId = min(e.src, e.dst);
        const vector<pair<int, int>> &srcAdj = q->GetAdj(e.src, true);
        for (const pair<int, int> &nbr : srcAdj) {
            if (nbr.first < minVId) continue;
            extensions.emplace_back(make_tuple(Edge::FORWARD, Edge::FORWARD, e, Edge(e.src, nbr.first, nbr.second)));
        }

        const vector<pair<int, int>> &srcInAdj = q->GetAdj(e.src, false);
        for (const pair<int, int> &nbr : srcInAdj) {
            if (nbr.first < minVId) continue;
            extensions.emplace_back(make_tuple(Edge::FORWARD, Edge::BACKWARD, e, Edge(nbr.first, e.src, nbr.second)));
        }

        const vector<pair<int, int>> &dstAdj = q->GetAdj(e.dst, true);
        for (const pair<int, int> &nbr : dstAdj) {
            if (nbr.first < minVId) continue;
            extensions.emplace_back(make_tuple(Edge::BACKWARD, Edge::FORWARD, e, Edge(e.dst, nbr.first, nbr.second)));
        }

        const vector<pair<int, int>> &dstInAdj = q->GetAdj(e.dst, false);
        for (const pair<int, int> &nbr : dstInAdj) {
            if (nbr.first < minVId) continue;
            extensions.emplace_back(make_tuple(Edge::BACKWARD, Edge::BACKWARD, e, Edge(nbr.first, e.dst, nbr.second)));
        }
    }
}

double MarkovTable::calcExtRate(const tuple<int, int, Edge, Edge> &ext) {
    const int &baseDir = get<0>(ext);
    const int &extDir = get<1>(ext);
    const int &baseEl = get<2>(ext).el;
    const int &extEl = get<3>(ext).el;
    const double &denom = mt1_[baseEl];
    if (baseDir < extDir) {
        return mt2_[baseDir][extDir][baseEl][extEl] / denom;
    } else if (baseDir > extDir) {
        return mt2_[extDir][baseDir][extEl][baseEl] / denom;
    } else if (baseEl < extEl) {
        return mt2_[baseDir][extDir][baseEl][extEl] / denom;
    } else {
        return mt2_[extDir][baseDir][extEl][baseEl] / denom;
    }
}