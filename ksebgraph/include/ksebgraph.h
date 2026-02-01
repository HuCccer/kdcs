#ifndef KSEBGRAPH_H_
#define KSEBGRAPH_H_
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

struct GN {
    int id_; uint32_t delta_;
    vector<pair<int, int>> adj;   // children TN ids
    // elements_:
    // - build phase: stores SGN ids
    // - finalize phase: stores edge ids
    vector<int> elements_;
    GN() : id_(-1), delta_(UINT32_MAX) {}
    GN(uint32_t delta, int sgnid) : id_(sgnid), delta_(delta){};
};


struct Graph {
    int maxBlockId;
    vector<int> etoGraphNID;  // edge ID → Tree node ID
    unordered_map<int, GN> idGN; // Super node ID → Tree node object
    Graph(int m) : etoGraphNID(m, -1), maxBlockId(0) {};
};


/***************************************************
 * Abstract Base Class: KDEquiForests
 ***************************************************/
class KSEBGraphs{
public:
    int num_; // number of forest 
    int m_; // number of edges
    // Original vertex → set of super nodes
    vector<Graph> graphs;
    KSEBGraphs() = default; 
    KSEBGraphs(int num, int m) : num_(num), m_(m), graphs(num_, Graph(m)){};
    KSEBGraphs(const KSEBGraphs&) = delete;
    KSEBGraphs& operator=(const KSEBGraphs&) = delete;
    ~KSEBGraphs() {}

};

#endif
