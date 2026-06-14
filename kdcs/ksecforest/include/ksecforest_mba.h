#ifndef KSECFOREST_MBAF_H
#define KSECFOREST_MBAF_H
#include <iostream>
#include <unordered_map>
#include <cstdint>
#include <string>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <cstring>
#include <list>
#include <optional>
#include "../../infrastructure/include/tgraph.h"
#include "../../infrastructure/include/anchor_union_find.h"
#include "../../infrastructure/include/KDeltaTrussDecomp.h"

namespace KDeltaTrussCommunity {

typedef std::unordered_set<int64_t> EdgeHash;
typedef std::pair<uint32_t, uint32_t> TriId;
typedef vector<uint32_t> KSEB;

struct TN {
    uint32_t id_; uint32_t delta_;
    uint32_t father_;              // father TN id, -1 if none
    vector<uint32_t> children_;   // children TN ids
    // elements_:
    // - build phase: stores SGN ids
    // - finalize phase: stores edge ids
    vector<uint32_t> elements_;
    TN() : id_(UINT32_MAX), father_(UINT32_MAX), delta_(UINT32_MAX) {}
    TN(uint32_t delta, uint32_t sgnid) : id_(sgnid), father_(UINT32_MAX), delta_(delta){};
};


struct Forest {
    uint32_t maxBlockId;
    vector<uint32_t> etoTreeNID;  // edge ID → Tree node ID
    unordered_map<uint32_t, TN> idTN; // Super node ID → Tree node object
    Forest(uint32_t m) : etoTreeNID(m, UINT32_MAX), maxBlockId(0) {};
};


/***************************************************
 * Abstract Base Class: KDEquiForests
 ***************************************************/
struct KDEquiForests{
public:
    uint32_t num_; // number of forest 
    uint32_t m_; // number of edges
    // Original vertex → set of super nodes
    std::vector<Forest> forests;
    KDEquiForests() = default; 
    KDEquiForests(uint32_t num, uint32_t m) : num_(num), m_(m), forests(num_, Forest(m)){};
    KDEquiForests(const KDEquiForests&) = delete;
    KDEquiForests& operator=(const KDEquiForests&) = delete;
    ~KDEquiForests() {}

};


class MbaF final {
public:
    MbaF(uint32_t n, uint32_t m);
    MbaF(uint32_t n, uint32_t m, const std::string& fn);
    MbaF(const MbaF&) = delete;
    MbaF& operator=(const MbaF&) = delete;
    ~MbaF() {}
    // load graph from file, the file format is: t n m, followed by m lines of t u v
    void loadGraph(const std::string& file_name);
    // load graph index from file
    void loadGraphIndex(const std::string& fn);
    // write graph index to file
    void writeGraphToFile(const std::string& filename);
    // write KSECForest index to file
    void writeKSECForestIndexToFile(const std::string& filename);
    // load KSECForest index from file
    void readKSECForestIndexFromFile(const std::string& filename);
    // two-phase construction of KSECForest index
    void KSECForestConstruction();
    // from bottom to up merging of KSECForest index
    void constructForestForK(uint32_t k, vector<uint32_t>& edgeProcessed);
    // serial merging of tree nodes in KSECForest index
    void Serialmerge(uint32_t seed, uint32_t k, unordered_set<uint32_t>& treeNodeToMerge);
    // merging of tree nodes with the same k-span in KSECForest index
    uint32_t sameKSpanMerging(const std::vector<uint32_t> &toMerge, uint32_t k);
    // insert edge into the graph and update the KSECForest index
    void addEdge(uint32_t u, uint32_t v, uint32_t timestamps);
    // search k-delta-truss community by KSpanIndex
    vector<vector<int>> searchCommunityByKSpanIndex(int query, int k, uint32_t delta);
    // search k-delta-truss community by KSECForest index
    vector<vector<uint32_t>> searchCommunityByKSECForest(uint32_t query, uint32_t k, uint32_t delta);
    // construct forest by graph
    void constructForestByGraph(Forest& forest, vector<vector<pair<uint32_t, uint32_t>>>& waitList, vector<vector<pair<uint32_t, uint32_t>>>& superE);

    KDeltaTruss& getKDeltaTruss() { return kdt_; }

    inline void dealwithWaitedSE(uint32_t rx, uint32_t ry, vector<uint32_t>& BidtoCid, Forest& forest, uint32_t delta);
    inline bool tryInsertEdge(std::unordered_set<int64_t>& seen, uint32_t u, uint32_t v) {
        int64_t key = (int64_t(u) << 32) | int64_t(v);
        return seen.insert(key).second;
    }
    inline pair<uint32_t, uint32_t> findAncestorAndFatherLessD(unordered_map<uint32_t, TN>& idtn, uint32_t nodeId, uint32_t delta) {
        uint32_t fatherId = idtn[nodeId].father_;
        while (fatherId != UINT32_MAX && idtn[fatherId].delta_ < delta) {
            nodeId = fatherId;
            fatherId = idtn[fatherId].father_;
        }  
        return {nodeId, fatherId};
    }

    inline uint32_t findAncestorEqualOrLargeD(unordered_map<uint32_t, TN>& idtn, uint32_t nodeId, uint32_t delta) {
        if (idtn[nodeId].delta_ >= delta) return nodeId;
        uint32_t fatherId = idtn[nodeId].father_;
        while (fatherId != UINT32_MAX && idtn[fatherId].delta_ < delta) {
            fatherId = idtn[fatherId].father_;
        }  
        return fatherId;
    }
     
    inline void mergeNode(uint32_t xid, uint32_t yid, unordered_map<uint32_t, TN>& idtn, vector<uint32_t>& BidToCid) {
        TN &nodex = idtn[xid], &nodey = idtn[yid];
        if (nodex.children_.size() < nodey.children_.size()) {
            for (uint32_t id : nodex.elements_) {
                nodey.elements_.push_back(id); BidToCid[id] = yid;    
            }
            for (uint32_t cid : nodex.children_) {
                idtn[cid].father_ = yid;
                nodey.children_.push_back(cid);
            }
            idtn.erase(xid);
        } else {
            for (uint32_t id : nodey.elements_) {
                nodex.elements_.push_back(id); BidToCid[id] = xid;    
            }
            for (uint32_t cid : nodey.children_) {
                idtn[cid].father_ = xid;
                nodex.children_.push_back(cid);
            }
            idtn.erase(yid);
        }
        
    }

    inline bool isSuperEdge(TN& father, TN& child) {
        return child.father_ == father.id_;
    }

    inline void addSuperEdge(TN& father, TN& child) {
        father.children_.push_back(child.id_);
        child.father_ = father.id_;
    }

    inline void delSuperEdge(TN& father, TN& child) {
        father.children_.erase(std::remove(father.children_.begin(), father.children_.end(), child.id_), father.children_.end());
        child.father_ = UINT32_MAX;
        return;
    }
    
    void initKSECForest(uint32_t k) {
        ksecforest.emplace(k, l_);
    }

    void updateKSECForest() {
        ksecforest->forests.push_back(Forest(l_));
    }

    void initAUF(uint32_t m) {
        auf.emplace(m);
    }

    void getKSpanDict(uint32_t k, vector<vector<uint32_t>>& kspan_dict){
        auto& kspan = kdt_.span_[k];
        uint32_t m = tg_.m_;
        for (uint32_t i = 0; i < m; i++) {
            uint32_t sp = kspan[i];
            if (sp != UINT32_MAX) {
                kspan_dict[sp].push_back(i);
            }   
        }
    }

    void print_info(uint32_t t_, uint32_t n_, uint32_t m_) {
        cout << "-----------info--------------" << endl;
        cout << "Number of timestamps: " << t_ << endl;
        cout << "Number of vertices: " << n_ << endl;
        cout << "Number of temporal edges: " << m_ << endl;
        cout << "-----------------------------" << endl;
    }


    
    uint32_t kmax() { return kmax_; };
private:
    uint32_t l_;
    uint32_t n_;
    uint32_t m_;
    uint32_t kmax_ = 0;	//the # maximum k value
    uint32_t tmax_ = 0; //the #maximum timestamps
    vector<uint32_t> k_;
    unordered_map<TriId, uint32_t, pairHash> triMTSpan;
    optional<AnchorUnionFind> auf;
    optional<KDEquiForests> ksecforest;
    KDeltaTruss kdt_;
    TGraph tg_;
    vector<uint32_t> unoccupiedTreeId;
};
} // namespace KDeltaTrussCommunity
#endif
