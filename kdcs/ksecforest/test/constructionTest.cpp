#include <chrono>
#include <fstream>
#include <string>
#include <filesystem>
#include "../include/ksecforest_mba.h"

// a sample program
int main(int argc, char** argv) {
    std::string graph_file = argv[1];
    std::ifstream infile(graph_file);
    uint32_t t, n, m;
    infile >> t >> n >> m;
    infile.close();
    KDeltaTrussCommunity::MbaF mbaF(n, 1.2 * m, graph_file);
    std::string base = std::filesystem::path(graph_file).filename().string();
    auto& kdt = mbaF.getKDeltaTruss();
    auto beg = std::chrono::steady_clock::now();  
    kdt.kDeltaTrussDecomp();
    mbaF.KSECForestConstruction();
    auto end = std::chrono::steady_clock::now();
    printf("Index construction costs %.5f ms\n",
            std::chrono::duration<double, std::milli>(end - beg).count());
    mbaF.writeGraphToFile("../index/graph_" + base);
    kdt.writeKspanIndexToFile("../index/kspan_" + base);
    mbaF.writeKSECForestIndexToFile("../index/forest_" + base); 
    return 0;
}
