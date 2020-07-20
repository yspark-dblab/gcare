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
        mt1_ = vector<long>(g->GetNumELabels());
        for (int i = 0; i < 2; ++i) {
            mt2_.emplace_back(vector<vector<vector<long>>>(2));
            for (int j = 0; j < 2; ++j) {
                mt2_[i].emplace_back(vector<vector<long>>(g->GetNumELabels()));
                for (int k = 0; k < g->GetNumELabels(); ++k) {
                    mt2_[i][j].emplace_back(vector<long>(g->GetNumELabels(), 0));
                }
            }
        }

        while (getline(catalogueFile, line)) {
            vector<string> entry;
            boost::split(entry, line, boost::is_any_of(","));
            insertEntryToMT(entry);
        }
        catalogueFile.close();
    }
}

void MarkovTable::insertEntryToMT(const vector<string> &entry) {
    const long count = stol(entry[2]);
    vector<string> directions;
    boost::split(directions, entry[0], boost::is_any_of(";"));
    vector<string> labels;
    boost::split(labels, entry[1], boost::is_any_of(";"));
    if (labels.size() == 1) {
        mt1_[stoi(labels[0])] = count;
    } else if (labels.size() == 2) {
        mt2_[stoi(directions[0])][stoi(directions[1])][stoi(labels[0])][stoi(labels[1])] = count;
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
    return EstCardAllMax(subquery_index);
}

double MarkovTable::EstCardAllMax(int subquery_index) {
    vector<tuple<int, int, Edge, Edge>> twoPaths;
    q->getAll2Paths(twoPaths);

    const int NUM_SUBQ = 1 << q->GetNumEdges();
    int subQEnc, currentEnc;
    vector<double> estimates(NUM_SUBQ);
    deque<int> queue;
    vector<SubQEdgeNode> subQEdgeNodePool;
    subQEdgeNodePool.reserve(NUM_SUBQ / 2);
    for (const tuple<int, int, Edge, Edge> &twoPath : twoPaths) {
        subQEnc = (1 << get<2>(twoPath).id) | (1 << get<3>(twoPath).id);
        subQEdgeNodePool.emplace_back(0, get<2>(twoPath), 1 << get<2>(twoPath).id, -1);
        subQEdgeNodePool.emplace_back(1, get<3>(twoPath), subQEnc, subQEdgeNodePool.size() - 1);
        queue.push_back(subQEdgeNodePool.size() - 1);
        estimates[subQEnc] = mt2_[get<0>(twoPath)][get<1>(twoPath)][get<2>(twoPath).el][get<3>(twoPath).el];
    }

    double currentEst;
    vector<bool> processed(NUM_SUBQ, false);
    while (!queue.empty()) {
        int currentIdx = queue.front();
        queue.pop_front();
        SubQEdgeNode *current = &subQEdgeNodePool[currentIdx];
        currentEnc = current->currentEnc;
        currentEst = estimates[currentEnc];

        int i = current->count;
        while (i >= 0) {
            vector<Edge> srcAdj = q->GetAdjE(current->edge.src, true);
            for (const Edge &extE : srcAdj) {
                if (current->edge.dst == extE.dst && current->edge.el == extE.el) continue;
                if ((q->encodeSubQ(extE) & currentEnc) != 0) continue;
                makeEstAndAddToQueue(subQEdgeNodePool, currentIdx, currentEst,
                        make_tuple(Edge::FORWARD, Edge::FORWARD, current->edge, extE), queue, processed, estimates);
            }

            vector<Edge> srcInAdj = q->GetAdjE(current->edge.src, false);
            for (const Edge &extE : srcInAdj) {
                if ((q->encodeSubQ(extE) & currentEnc) != 0) continue;
                makeEstAndAddToQueue(subQEdgeNodePool, currentIdx, currentEst,
                        make_tuple(Edge::FORWARD, Edge::BACKWARD, current->edge, extE), queue, processed, estimates);
            }

            vector<Edge> dstAdj = q->GetAdjE(current->edge.dst, true);
            for (const Edge &extE : dstAdj) {
                if ((q->encodeSubQ(extE) & currentEnc) != 0) continue;
                makeEstAndAddToQueue(subQEdgeNodePool, currentIdx, currentEst,
                        make_tuple(Edge::BACKWARD, Edge::FORWARD, current->edge, extE), queue, processed, estimates);
            }

            vector<Edge> dstInAdj = q->GetAdjE(current->edge.dst, false);
            for (const Edge &extE : dstInAdj) {
                if (current->edge.src == extE.src && current->edge.el == extE.el) continue;
                if ((q->encodeSubQ(extE) & currentEnc) != 0) continue;
                makeEstAndAddToQueue(subQEdgeNodePool, currentIdx, currentEst,
                        make_tuple(Edge::BACKWARD, Edge::BACKWARD, current->edge, extE), queue, processed, estimates);
            }

            --i;
            current = &subQEdgeNodePool[current->prevIdx];
        }
    }

    return estimates[NUM_SUBQ - 1];
}

void MarkovTable::makeEstAndAddToQueue(vector<SubQEdgeNode> &subQEdgeNodePool, const int &currentIdx, const double &currentEst,
        const tuple<int, int, Edge, Edge> &ext, deque<int> &queue,
        vector<bool> &processed, vector<double> &estimates) {
    double extRate = calcExtRate(ext);
    const SubQEdgeNode &current = subQEdgeNodePool[currentIdx];
    int subQEnc = current.currentEnc | (1 << get<3>(ext).id);

    if (processed[subQEnc]) {
        estimates[subQEnc] = max(estimates[subQEnc], currentEst * extRate);
    } else {
        estimates[subQEnc] = currentEst * extRate;
        subQEdgeNodePool.emplace_back(current.count + 1, get<3>(ext), subQEnc, currentIdx);
        queue.emplace_back(subQEdgeNodePool.size() - 1);
        processed[subQEnc] = true;
    }
}

double MarkovTable::EstCardGreedyMax(int subquery_index) {
    vector<tuple<int, int, Edge, Edge>> twoPaths;
    q->getAll2Paths(twoPaths);

    double est;
    pair<vector<Edge>, double> current = pair<vector<Edge>, double>(NULL, numeric_limits<double>::min());
    for (const tuple<int, int, Edge, Edge> &twoPath : twoPaths) {
        est = mt2_[get<0>(twoPath)][get<1>(twoPath)][get<2>(twoPath).el][get<3>(twoPath).el];
        if (est > current.second) {
            vector<Edge> starter(2);
            starter[0] = get<2>(twoPath);
            starter[1] = get<3>(twoPath);
            current.first = starter;
            current.second = est;
        }
    }

    int subQEnc;
    double extRate;
    const int NUM_Q_EDGES = q->GetNumEdges();
    pair<vector<Edge>, double> maxNext;
    while (current.first.size() < NUM_Q_EDGES) {
        subQEnc = q->encodeSubQ(current.first);
        vector<tuple<int, int, Edge, Edge>> extensions;
        extensions.reserve(pow(q->GetNumEdges(), 2));
        getExtensions(extensions, current.first, subQEnc);

        maxNext = pair<vector<Edge>, double>(NULL, numeric_limits<double>::min());
        for (const tuple<int, int, Edge, Edge> &ext : extensions) {
            extRate = calcExtRate(ext);
            if (extRate > maxNext.second) {
                vector<Edge> nextSubQ;
                nextSubQ.reserve(current.first.size() + 1);
                for (const Edge &e : current.first) {
                    nextSubQ.push_back(e);
                }
                nextSubQ.push_back(get<3>(ext));
                maxNext.first = nextSubQ;
                maxNext.second = extRate;
            }
        }

        maxNext.second *= current.second;
        current = maxNext;
    }

    return current.second;
}

double MarkovTable::AggCard() {
    return card_vec_[0];
}

double MarkovTable::GetSelectivity() {
    return 1;
}

void MarkovTable::getExtensions(vector<tuple<int, int, Edge, Edge>> &extensions, const vector<Edge> &current, const int &currentEnc) {
    Edge extE;
    for (const Edge &e : current) {
        vector<pair<int, int>> srcAdj = q->GetAdj(e.src, true);
        for (const pair<int, int> &nbr : srcAdj) {
            if (e.dst == nbr.first && e.el == nbr.second) continue;
            extE = Edge(e.src, nbr.first, nbr.second);
            if ((q->encodeSubQ(extE) & currentEnc) != 0) continue;
            extensions.emplace_back(make_tuple(Edge::FORWARD, Edge::FORWARD, e, extE));
        }

        vector<pair<int, int>> srcInAdj = q->GetAdj(e.src, false);
        for (const pair<int, int> &nbr : srcInAdj) {
            extE = Edge(nbr.first, e.src, nbr.second);
            if ((q->encodeSubQ(extE) & currentEnc) != 0) continue;
            extensions.emplace_back(make_tuple(Edge::FORWARD, Edge::BACKWARD, e, extE));
        }

        vector<pair<int, int>> dstAdj = q->GetAdj(e.dst, true);
        for (const pair<int, int> &nbr : dstAdj) {
            extE = Edge(e.dst, nbr.first, nbr.second);
            if ((q->encodeSubQ(extE) & currentEnc) != 0) continue;
            extensions.emplace_back(make_tuple(Edge::BACKWARD, Edge::FORWARD, e, extE));
        }

        vector<pair<int, int>> dstInAdj = q->GetAdj(e.dst, false);
        for (const pair<int, int> &nbr : dstInAdj) {
            if (e.src == nbr.first && e.el == nbr.second) continue;
            extE = Edge(nbr.first, e.dst, nbr.second);
            if ((q->encodeSubQ(extE) & currentEnc) != 0) continue;
            extensions.emplace_back(make_tuple(Edge::BACKWARD, Edge::BACKWARD, e, extE));
        }
    }
}

double MarkovTable::calcExtRate(const tuple<int, int, Edge, Edge> &ext) {
    int baseDir = get<0>(ext);
    int extDir = get<1>(ext);
    int baseEl = get<2>(ext).el;
    int extEl = get<3>(ext).el;
    double denom = mt1_[baseEl];
    if (baseDir < extDir) {
        return mt2_[baseDir][extDir][baseEl][extEl] / denom;
    } else if (baseDir > extDir) {
        return mt2_[extDir][baseDir][extEl][baseEl] / denom;
    } else {
        if (baseEl < extEl) {
            return mt2_[baseDir][extDir][baseEl][extEl] / denom;
        } else {
            return mt2_[extDir][baseDir][extEl][baseEl] / denom;
        }
    }
}