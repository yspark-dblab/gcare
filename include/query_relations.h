#ifndef QUERY_RELATIONS_H_
#define QUERY_RELATIONS_H_

#include <vector>
#include <string>
#include <cstring>
#include <unordered_map>
#include <cassert>
#include <iostream>

using std::vector;
class QueryGraph {
public:

    struct Attribute {
        int id; // id of the attribute. 
        bool is_bound = false;
        int bound; //
        int pos; // position at the data table
        int ref_cnt; // reference count
        Attribute(int _id, bool _is_bound, int _bound, int _pos, int _ref_cnt)
        : id(_id), is_bound(_is_bound), bound(_bound), pos(_pos), ref_cnt(_ref_cnt) {}
        Attribute() {}
        void print() { std::cout << "a " << id << " " << is_bound << " " << bound << " " << pos << " " << ref_cnt << "\n"; }
    };

    struct Relation {
        int id;
        std::vector<Attribute> attrs;

        Relation(int _id, std::vector<Attribute> _attrs): id(_id), attrs(_attrs) {}
        Relation() {}
        Attribute &operator[](size_t i) {
            return attrs[i];
        }
        void print() {
            std::cout << "r " << id << "\n";
            for (int i = 0; i < attrs.size(); ++i) attrs[i].print();
        }
        size_t size() { return attrs.size(); }
    };

    std::vector<Relation> relations_;
    size_t num_attrs_;

    size_t num_attrs() { return num_attrs_; }
    size_t num_relations() { return relations_.size(); }

    QueryGraph() {
        relations_.clear();
        num_attrs_ = 0;
    }

    struct Edge { 
        int srcid, dstid, lbl;
        Edge(int _srcid, int _dstid, int _lbl): srcid(_srcid), dstid(_dstid), lbl(_lbl) {}
        Edge() {}
        void print() { std::cout << "e " << srcid << " " << dstid << " " << lbl << "\n"; }
    };

    struct Vertex { 
        int id, lbl, bound; 
        Vertex(int _id, int _lbl, int _bound): id(_id), lbl(_lbl), bound(_bound) {} 
        Vertex() {}
        void print() { std::cout << "v " << id << " " << lbl << " " << bound <<  "\n"; }
    };

    void ParseFile(const char* filename, std::vector<Vertex> &vertices, std::vector<Edge> &edges) {
        char buf[1111];
        FILE* fp = fopen(filename, "r");
        while (fgets(buf, 1111, fp) != NULL) {
            if (buf[0] == '\n') continue;
            int buflen = 0;
            while (buf[buflen++] != '\n') {
                if (buf[buflen] == ' ') {
                    buf[buflen] = 0;
                }
            }
            if (buf[0] == 'v') { // parse vertex
                int id = strtol(buf + 2, 0, 0);
                int offset = 2 + strlen(buf + 2) + 1;
                int lbl = strtol(buf + offset, 0, 0);
                offset += strlen(buf + offset) + 1;
                int bound = -1;
                if (offset < buflen)
                    bound = strtol(buf + offset, 0, 0);
                vertices.emplace_back(id, lbl, bound);
            } else if (buf[0] == 'e') { // parse edge
                int srcid = strtol(buf + 2, 0, 0);
                int offset = 2 + strlen(buf + 2) + 1;
                int dstid = strtol(buf + offset, 0, 0);
                offset += strlen(buf + offset) + 1;
                int lbl = strtol(buf + offset, 0, 0);
                edges.emplace_back(srcid, dstid, lbl); 
            }
        }
        fclose(fp);
    }

    void ParseText(std::vector<std::string> &text, std::vector<Vertex> &vertices, std::vector<Edge> &edges) {
        for (auto &line : text) {
            char* buf = (char*) line.c_str();
            if (buf[0] == '\n') continue;
            int buflen = 0;
            while (buf[buflen++] != '\n') {
                if (buf[buflen] == ' ') {
                    buf[buflen] = 0;
                }
            }
            if (buf[0] == 'v') { // parse vertex
                int id = strtol(buf + 2, 0, 0);
                int offset = 2 + strlen(buf + 2) + 1;
                int lbl = strtol(buf + offset, 0, 0);
                offset += strlen(buf + offset) + 1;
                int bound = -1;
                if (offset < buflen)
                    bound = strtol(buf + offset, 0, 0);
                vertices.emplace_back(id, lbl, bound);
            } else if (buf[0] == 'e') { // parse edge
                int srcid = strtol(buf + 2, 0, 0);
                int offset = 2 + strlen(buf + 2) + 1;
                int dstid = strtol(buf + offset, 0, 0);
                offset += strlen(buf + offset) + 1;
                int lbl = strtol(buf + offset, 0, 0);
                edges.emplace_back(srcid, dstid, lbl); 
            }
        }
    }

    void ConvertToRelations(std::vector<Vertex> &vertices, std::vector<Edge> &edges) {
        //=============================================
        //1. Find valid attributes (i.e. join attributes and bound attributes)
	    //---------------------------------------------
        std::unordered_map<int, int> bounds;
        std::unordered_map<int, int> ref_cnts;
        for (size_t i = 0; i < vertices.size(); ++i) {
            Vertex &v = vertices[i];
            bounds[v.id] = v.bound;
            if (v.bound != -1) {
            } else if (v.lbl != -1) {
                if (ref_cnts.count(v.id) == 0) {
                    ref_cnts[v.id] = 1;
                } else {
                    ref_cnts[v.id] += 1;
                }
            }
        }
        for (size_t i = 0; i < edges.size(); ++i) {
            Edge &e = edges[i];
            ref_cnts[e.srcid] += 1;
            ref_cnts[e.dstid] += 1;
        }
        //=============================================
        // 2. Create releations using the valid attributes
        //---------------------------------------------
        relations_.clear();
        // 2-1 Create relations for vertices
        for (size_t i = 0; i < vertices.size(); ++i) {
            Vertex &v = vertices[i];
            if (v.lbl != -1) { // add relation
                std::vector<Attribute> attrs;
                attrs.emplace_back(v.id, v.bound != -1, v.bound, 0, ref_cnts[v.id]);
                relations_.emplace_back(-(v.lbl + 1), attrs);
            }
        }
        // 2-2 Create relations for edges
        for (size_t i = 0; i < edges.size(); ++i) {
            Edge &e = edges[i];
            std::vector<Attribute> attrs;
            if (bounds[e.srcid] != -1 || ref_cnts[e.srcid] > 1) {
                attrs.emplace_back(e.srcid, bounds[e.srcid] != -1, bounds[e.srcid], 0, ref_cnts[e.srcid]);
            } 
            if (bounds[e.dstid] != -1 || ref_cnts[e.dstid] > 1) {
                attrs.emplace_back(e.dstid, bounds[e.dstid] != -1, bounds[e.dstid], 1, ref_cnts[e.dstid]);
            }
            assert(attrs.size() > 0);
            relations_.emplace_back(e.lbl, attrs);
        }
        //=============================================
        //=============================================
        // 4. Change attribute ids to contiguous integers from 0
        //---------------------------------------------
        std::unordered_map<int, int> vmap;
        num_attrs_ = 0;
        for (size_t i = 0; i < relations_.size(); ++i) {
            Relation &r = relations_[i];
            for (size_t j = 0; j < r.attrs.size(); ++j) {
                Attribute &a = r.attrs[j];
                if (vmap.count(a.id) == 0) {
                    vmap[a.id] = num_attrs_;
                    a.id = num_attrs_++;
                } else {
                    a.id = vmap[a.id];
                }
            }
        }
    }


    void ReadText(const char* filename) {
        std::vector<Vertex> vertices;
        std::vector<Edge> edges;
        //=============================================
        //1. Parse the query graph
	    //---------------------------------------------
        ParseFile(filename, vertices, edges);
	    //---------------------------------------------
        //2. Convert the query graph to a query relations
        ConvertToRelations(vertices, edges);
        //=============================================
    }

    void ReadText(std::vector<std::string> &text) {
        std::vector<Vertex> vertices;
        std::vector<Edge> edges;
        //=============================================
        //1. Parse the query graph
	    //---------------------------------------------
        ParseText(text, vertices, edges);
	    //---------------------------------------------
        //2. Convert the query graph to a query relations
        ConvertToRelations(vertices, edges);
        //=============================================
    }

/*
    vector<int> table_;

    struct Column {
        Column(int x, int y, int z, int w): vid(x), cls(y), vlabel(z), vidattr(w) { }
        int vid, cls, vlabel, vidattr;
    };

    vector<vector<Column>> condition_;
    vector<int> par_, cnt_;
    int max_vid_, max_vlabel_, max_elabel_, enum_;

    int Find(int);
    void Unite(int, int);
    QueryGraph(void);

    void ReadText(const char*);
*/
};
#endif
