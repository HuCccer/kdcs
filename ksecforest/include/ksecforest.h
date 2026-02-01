#ifndef KSECFOREST_H_
#define KSECFOREST_H_
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <list>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdint>

using namespace std;

typedef vector<int> KSEB;

struct TN {
    int id_; uint32_t delta_;
    int father_;              // father TN id, -1 if none
    vector<int> children_;   // children TN ids
    // elements_:
    // - build phase: stores SGN ids
    // - finalize phase: stores edge ids
    vector<int> elements_;
    TN() : id_(-1), father_(-1), delta_(UINT32_MAX) {}
    TN(uint32_t delta, int sgnid) : id_(sgnid), delta_(delta), father_(-1) {};
};


struct Forest {
    int maxBlockId;
    vector<int> etoTreeNID;  // edge ID → Tree node ID
    unordered_map<int, TN> idTN; // Super node ID → Tree node object
    Forest(int m) : etoTreeNID(m, -1), maxBlockId(0) {};
};


/***************************************************
 * Abstract Base Class: KDEquiForests
 ***************************************************/
class KDEquiForests{
public:
    int num_; // number of forest 
    int m_; // number of edges
    // Original vertex → set of super nodes
    vector<Forest> forests;
    KDEquiForests() = default; 
    KDEquiForests(int num, int m) : num_(num), m_(m), forests(num_, Forest(m)){};
    KDEquiForests(const KDEquiForests&) = delete;
    KDEquiForests& operator=(const KDEquiForests&) = delete;
    ~KDEquiForests() {}

};
#endif