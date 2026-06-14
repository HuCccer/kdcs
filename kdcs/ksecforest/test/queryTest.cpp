#include <chrono>
#include <fstream>
#include <string>
#include <filesystem>
#include "../include/ksecforest_mba.h"

// a sample program
int main(int argc, char** argv) {
    // load graph index
    const std::string graph_index_file = argv[1];
    // load kspan index
    const std::string kspan_index_file = argv[2];
    // load ksecforest index
    const std::string ksecforest_index_file = argv[3];
    // load queries
    const std::string queries_file = argv[4];
    // parameters for query
    int k = stoi(argv[5]), delta = stoi(argv[6]);
    std::ifstream infile(graph_index_file, std::ios::binary);
    uint32_t t, n, m;
    infile.read(reinterpret_cast<char*>(&t), sizeof(t))
          .read(reinterpret_cast<char*>(&n), sizeof(n))
          .read(reinterpret_cast<char*>(&m), sizeof(m));
    KDeltaTrussCommunity::MbaF mbaF(n, m * 1.5);
    infile.close();
    mbaF.loadGraphIndex(graph_index_file);
    auto& kdt = mbaF.getKDeltaTruss();
    kdt.loadKSpanIndex(kspan_index_file);
    mbaF.readKSECForestIndexFromFile(ksecforest_index_file);
    // read the queries
    std::vector<int> query_vertices;
    std::ifstream qf(queries_file);
    int v;
    while (qf >> v) query_vertices.push_back(v);
    auto beg = std::chrono::steady_clock::now();
    for (uint32_t vertex : query_vertices) {
        //search (k, delta)-Truss Community by KSpan index
        // mbaF.searchCommunityByKSpanIndex(vertex, k, delta);
        //search (k, delta)-Truss Community by KSECForest index
        mbaF.searchCommunityByKSECForest(vertex, k, delta);
    }
    auto end = std::chrono::steady_clock::now();
    const auto dif = end - beg;
    printf("Applying the query costs %f ms.\n",
          std::chrono::duration<double, std::milli>(dif).count());
    return 0;
}
