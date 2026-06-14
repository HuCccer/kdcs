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
    // load edges to update
    const std::string update_file = argv[4]; 
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
    // read the edges to update
    uint32_t num_updates = 0;
    std::ifstream update_in(update_file);
    update_in >> num_updates;
    std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> updates;
    while (num_updates > 0) { 
        uint32_t t, u, v;
        update_in >> t >> u >> v;
        updates.emplace_back(t, u, v);
        num_updates--;
    }
    const auto beg = std::chrono::steady_clock::now();
    for (const auto& [t, u, v] : updates) {  
        // edge insertion maintainance
        mbaF.addEdge(u, v, t);
    }
    const auto end = std::chrono::steady_clock::now();
    const auto dif = end - beg;
    printf("Applying the edge insertion maintenance costs %f ms.\n",
          std::chrono::duration<double, std::milli>(dif).count());
    return 0;
}
