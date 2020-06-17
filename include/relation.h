#ifndef RELATION_H_
#define RELATION_H_

#include <iostream>
#include <cstdint>
#include <vector>
#include <stdlib.h> 
#include <cassert>
#include <cstring>
#include <algorithm>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include <sstream>

template<typename CellType>
struct Relation {

    enum ErrCode { MEMORY };

    uint64_t num_cols_;
    uint64_t num_rows_;
    uint64_t max_rows_;

    CellType* data_;

    struct Row {
        CellType* data_;
        uint64_t num_cols_;
        CellType &operator[](uint64_t c) {
            assert(c >= 0);
            assert(c < num_cols_);
            assert(data_);
            return data_[c];
        }
        Row(CellType* data=NULL, uint64_t num_cols=0): data_(data), num_cols_(num_cols) {}
        uint64_t size() { return num_cols_; }
    };

    Relation()
    : num_rows_(0), data_(NULL), max_rows_(0), num_cols_(0)
    {}
    
    Relation(uint64_t num_cols)
    : num_rows_(0), data_(NULL), max_rows_(0), num_cols_(num_cols)
    {}

    void SetNumCols(uint64_t num_cols) { num_cols_ = num_cols; }

    Relation(uint64_t num_cols, uint64_t num_rows)
    : num_rows_(num_rows), data_(NULL), max_rows_(num_rows), num_cols_(num_cols)
    {
        data_ = (CellType*) malloc(sizeof(CellType) * num_cols_ * num_rows);
    }

    ~Relation() {
        if (data_) free(data_);
    }
    void clear() {
        if (data_) free(data_);
        data_ = NULL;
        max_rows_ = 0;
        num_rows_ = 0;
    }

    Row operator[](uint64_t r) {
        assert(r >= 0);
        assert(r < num_rows_);
        return Row(data_ + (num_cols_ * r), num_cols_);
    }

    void handleOverflow() {
        uint64_t new_max_rows = max_rows_ == 0 ? 1 : max_rows_ * 2;
        CellType* new_data = (CellType*) malloc(sizeof(CellType) * num_cols_ * new_max_rows);

        if (!new_data) throw ErrCode::MEMORY;

        memcpy((void*) new_data, (void*) data_, sizeof(CellType) * num_cols_ * max_rows_);
        free((void*) data_);
        max_rows_ = new_max_rows;
        data_ = new_data;
    }
    
    void append(std::vector<CellType> &vals) {
        if (num_rows_ == max_rows_) handleOverflow();
        CellType* cells = data_ + (num_cols_ * num_rows_);
        for (uint64_t i = 0; i < num_cols_; ++i) {
            cells[i] = vals[i];
        }
        num_rows_++;
    }

    void append(CellType* vals) {
        if (num_rows_ == max_rows_) handleOverflow();
        CellType* cells = data_ + (num_cols_ * num_rows_);
        for (uint64_t i = 0; i < num_cols_; ++i) {
            cells[i] = vals[i];
        }
        num_rows_++;
    }

    void swap(Relation &other) {
        std::swap(data_, other.data_);        
        std::swap(num_cols_, other.num_cols_);
        std::swap(num_rows_, other.num_rows_);
        std::swap(max_rows_, other.max_rows_);
    }

    uint64_t size() { return num_rows_; }

    template <typename T>
    struct ContainerHash {
        size_t operator()(T const&c) const {
            return boost::hash_range(c.begin(), c.end());
        }
    };

    void join(std::vector<int> &lcols, std::vector<int> &rcols, Relation<CellType> &rr) {

        Relation<CellType> &lr = *this;
        if (lr.size() == 0 || rr.size() == 0) {
            lr.clear();
            rr.clear();
            return;
        }
        
        if (lr.size() > rr.size()) {
            lr.swap(rr);
            lcols.swap(rcols);
        }
       
        std::vector<int> ncols = lcols;
        std::vector<uint64_t> l_pos, r_pos, nr_pos;

        // 1. Detect join columns and fine the columns of the new relation
        for (size_t i = 0; i < rcols.size(); ++i) {
            bool pass = true;
            for (size_t j = 0; j < lcols.size(); ++j) {
                if (lcols[j] == rcols[i]) {
                    r_pos.push_back(i);
                    l_pos.push_back(j);
                    pass = false;
                    break;
                }
            }
            if (pass) {
                ncols.push_back(rcols[i]);
                nr_pos.push_back(i);
            }
        }

        // 2. Build hash table of the left relation
        std::unordered_multimap<std::vector<CellType>, uint64_t, ContainerHash<std::vector<CellType>>> hash_table;
        
        std::vector<CellType> key(l_pos.size());
        for (uint64_t i = 0; i < lr.size(); ++i) {
            for (size_t j = 0; j < key.size(); ++j) {
                key[j] = lr[i][l_pos[j]];
            }
            hash_table.emplace(key, i);
        }

        std::vector<CellType> tuple(ncols.size());
        Relation<CellType> nr(ncols.size());
        // 3. Do join
        for (uint64_t i = 0; i < rr.size(); ++i) {
            for (size_t j = 0; j < key.size(); ++j) {
                key[j] = rr[i][r_pos[j]];
            }
            auto range = hash_table.equal_range(key);
            for (auto it = range.first; it != range.second; ++it) {
                uint64_t idx = it->second;
                for (size_t j = 0; j < lcols.size(); ++j) {
                    tuple[j] = lr[idx][j];
                }
                for (size_t j = 0; j < nr_pos.size(); ++j) {
                    tuple[j + lcols.size()] = rr[i][nr_pos[j]];
                }
                nr.append(tuple);
            }
        }
        rr.clear();
        lcols.swap(ncols);
        lr.swap(nr);
    }
};

#endif
