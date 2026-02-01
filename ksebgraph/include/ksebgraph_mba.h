#ifndef KSEBGRAPH_MBA_H_
#define KSEBGRAPH_MBA_H_
#include <iostream>
#include <list>
#include <unordered_map>
#include<fstream>
#include <map>
#include <cstdint>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>
#include <optional>
#include "../../infrastructure/tgraph.h"
#include "ksebgraph.h"

using namespace std;

struct pair_hash {
    size_t operator()(const pair<int, int>& p) const {
        return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
    }
};
typedef unordered_set<int64_t> EdgeHash;

typedef struct final {
    uint32_t vid;
    uint32_t eid;
  } ArrayEntry;
// triangle type
typedef std::tuple<uint32_t, uint32_t, uint32_t> Triangle;
typedef std::unordered_set<int64_t> EdgeHash;
struct SuperEdge {
    int t, sp, x, y;
};
typedef std::vector<uint32_t> VI;

class MbaG final {
public:
    MbaG(const uint32_t n, const uint32_t m);
    MbaG(const uint32_t n, const uint32_t m, const std::string& fn);
    MbaG(const MbaG&) = delete;
    MbaG& operator=(const MbaG&) = delete;
    ~MbaG() {}
    void loadGraph(const std::string& dataset_path);
    void trussDecomp();
    void KdeltaTrussDecomp();
    void onePass();
    void twoPhase();
    void constructKSEBGraphForK(int k);
    vector<vector<int>> searchCommunityByKSEBraphs(int query, int k, uint32_t delta);
    void writeTGraph(const std::string& filename);
    void writeKSEBGraphsIndexToFile(const std::string& filename);
    void readTGraphFromFile(const std::string& filename);
    void readKSEBGraphsIndexFromFile(const std::string& filename);

    inline bool tryInsertEdge(std::unordered_set<int64_t>& seen, int u, int v) {
        int64_t key = (int64_t(u) << 32) | int64_t(v);
        return seen.insert(key).second;
    }

        
    void initKSEBGraphs(int k, int m) {
       ksebgraph.emplace(k, m);
    }

    void getKSpanDict(int k, vector<vector<int>>& kspan_dict){
        for (int i = 0; i < tg_.m_; i++) {
            uint32_t sp = kspan_[k][i];
            if (sp != UINT32_MAX) {
                kspan_dict[sp].push_back(i);
            }   
        }
    }

    inline void getminmax(uint32_t a, uint32_t b, uint32_t c){
        if (a > b) a ^= b ^= a ^= b;
        if (a > c) a ^= c ^= a ^= c;
        if (b > c){maxe = b; mine = a;}
        else {maxe = c; mine = a;}
    }

    inline Triangle sortTri(uint32_t a, uint32_t b, uint32_t c){
        if (a > b) a ^= b ^= a ^= b;
        if (a > c) a ^= c ^= a ^= c;
        if (b > c) b ^= c ^= b ^= c;
        return {a, b, c};
    }

    void print_info(uint32_t t_, uint32_t n_, uint32_t m_) {
        std::cout << "-----------info--------------" << endl;
        std::cout << "Number of timestamps: " << t_ << endl;
        std::cout << "Number of vertices: " << n_ << endl;
        std::cout << "Number of temporal edges: " << m_ << endl;
        std::cout << "-----------------------------" << endl;
    }
    
    int kmax() { return kmax_; };
private:
    const uint32_t l_;
    const uint32_t n_;
    optional<KSEBGraphs> ksebgraph;
    TGraph tg_;
    vector<vector<uint32_t>> kspan_;
    // vector<int> unoccupiedTreeId;
    vector<uint32_t> trn_;
    uint32_t m_;
    // data members
    uint32_t kmax_ = 0;	//the # maximum k value
    uint32_t tmax_ = 0; //the #maximum timestamps
    // the k-support
    VI ks_;
    // the edge peeling order
    VI ord_;
    //the position of edges in order 
    VI pos;
    uint32_t mine, maxe;
};
#endif
