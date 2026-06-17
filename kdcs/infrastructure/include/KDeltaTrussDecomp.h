#ifndef KDELTA_TRUSS_KDeltaTruss_H_
#define KDELTA_TRUSS_KDeltaTruss_H_
#include "tgraph.h"
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <map>

namespace KDeltaTrussCommunity{

typedef std::tuple<uint32_t, uint32_t, uint32_t> Triangle;
struct pairHash {
    std::size_t operator()(const std::pair<uint32_t, uint32_t>& p) const {
        return std::hash<uint32_t>{}(p.first) ^ (std::hash<uint32_t>{}(p.second) << 1);
    }
};
class MbaF; // forward declaration
class KDeltaTruss final {
public:
    friend class KDeltaTrussCommunity::MbaF;
    using TriIdMap = std::unordered_map<std::pair<uint32_t, uint32_t>, bool, pairHash>;
    KDeltaTruss(uint32_t n, uint32_t l, TGraph& tg, std::vector<uint32_t>& k, std::unordered_map<EdgT, uint32_t, pairHash>& triMTSpan);
    KDeltaTruss(uint32_t n, uint32_t l, TGraph& tg, std::vector<uint32_t>& k, std::unordered_map<EdgT, uint32_t, pairHash>& triMTSpan, const std::string& fn);
    KDeltaTruss(const KDeltaTruss&) = delete;
    KDeltaTruss& operator=(const KDeltaTruss&) = delete;
    ~KDeltaTruss() {}
    void trussDecomp();
    void kDeltaTrussDecomp();
    void loadKSpanIndex(const std::string& fn);
    void loadTriToMtsIndex(const std::string& fn);
    void writeKspanIndexToFile(const std::string& filename);
    void writeTriToMtsToFile(const std::string& filename);
    void insertT(const uint32_t e, const uint32_t timestamp);
    void insertE(const uint32_t e, const uint32_t timestamp);
    void assignKspanForEdges(std::vector<std::vector<uint32_t>>& inc, uint32_t k);
    void expand(TriIdMap& validTri, std::vector<uint32_t>& seeds, std::vector<uint32_t>& visited, uint32_t level, uint32_t k, uint32_t visittoken);
    std::vector<std::vector<uint32_t>> XH_truss_maintenance(const uint32_t e);
    void maintenance(uint32_t k, uint32_t seed);
    std::vector<std::unordered_set<uint32_t>> addEdgeWithVariety(uint32_t e, uint32_t timestamp, bool variety);
    std::vector<std::unordered_set<uint32_t>> addEdge(uint32_t u, uint32_t v, uint32_t timestamp);
    inline uint32_t getTriMts(uint32_t e1, uint32_t e2, uint32_t e3)  {
        if (e1 > e2) std::swap(e1, e2);
        if (e1 > e3) std::swap(e1, e3);
        if (e2 > e3) { return triMTSpan_[{e1, e2}];}
            else { return triMTSpan_[{e1, e3}];}
    }

    inline void setTriMts(uint32_t e1, uint32_t e2, uint32_t e3) {
        if (e1 > e2) std::swap(e1, e2);
        if (e1 > e3) std::swap(e1, e3);
        if (e2 > e3) { triMTSpan_[{e1, e2}] = tg_.GetMts(tg_.tau_[e1], tg_.tau_[e2], tg_.tau_[e3]);}
            else { triMTSpan_[{e1, e3}] = tg_.GetMts(tg_.tau_[e1], tg_.tau_[e2], tg_.tau_[e3]);}
    }
    
    inline Triangle sortTri(uint32_t e1, uint32_t e2, uint32_t e3){
        if (e1 > e2) std::swap(e1, e2);
        if (e1 > e3) std::swap(e1, e3);
        if (e2 > e3) std::swap(e2, e3);
        return {e1, e2, e3};
    }

    inline void getTriId(uint32_t a, uint32_t b, uint32_t c){
        if (a > b) a ^= b ^= a ^= b;
        if (a > c) a ^= c ^= a ^= c;
        if (b > c){maxeid = b; mineid = a;}
        else {maxeid = c; mineid = a;}
    }
    
    inline void addCandidate(uint32_t level, uint32_t eid) {
        if (candidateEdgeBuf_[level].empty()) {
            touchedCandidateLevels_.push_back(level);
        }
        candidateEdgeBuf_[level].push_back(eid);
    }

    inline void addDeltaTriangle(uint32_t level, Triangle& tri) {
        if (deltaTriListBuf_[level].empty()) {
            touchedDeltaLevels_.push_back(level);
        }
        deltaTriListBuf_[level].push_back(tri);
    }

    inline void clearMaintenanceBuffers(uint32_t level) {
        for (uint32_t mts : touchedDeltaLevels_) {
            if (mts < level) {
                for (auto& [e1, e2, e3] : deltaTriListBuf_[mts]) {
                    sup_[e1] = 0; sup_[e2] = 0; sup_[e3] = 0;
                } 
            }
            deltaTriListBuf_[mts].clear();
        }
        for (uint32_t ksp : touchedCandidateLevels_) {
            candidateEdgeBuf_[ksp].clear();
        } 
        touchedCandidateLevels_.clear();
        touchedDeltaLevels_.clear();
    }



    private:
    // members
    uint32_t l_;
    uint32_t n_;
    uint32_t m_;
    // the min and max edge ID in a triangle
    uint32_t mineid, maxeid;
    // graph and trussnesses
    std::vector<uint32_t>& k_;
    TGraph& tg_;
    std::vector<std::vector<uint32_t>> span_;
    std::vector<int32_t> sup_;
    // data members
    uint32_t kmax_ = 0;	//the # maximum k value
    uint32_t tmax_ = 0; //the #maximum minimum time span
    // the k-support
    std::vector<uint32_t> ks_;
    // the edge peeling order
    std::vector<uint32_t> ord_;
    //the position of edges in order 
    std::vector<uint32_t> pos;
    std::unordered_map<EdgT, uint32_t, pairHash>& triMTSpan_;
    std::vector<std::vector<Triangle>> time2triangle;
    std::vector<std::vector<uint32_t>> candidateEdgeBuf_;
    std::vector<std::vector<Triangle>> deltaTriListBuf_;
    std::vector<uint32_t> touchedCandidateLevels_;
    std::vector<uint32_t> touchedDeltaLevels_;
    std::vector<std::unordered_set<uint32_t>> affected;
    std::vector<uint32_t> visited, removed, ins;
	uint32_t token = 1, visittoken = 0;
    
};
}
#endif
