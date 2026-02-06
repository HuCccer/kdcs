#include <unistd.h>
#include <chrono>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem>

#include "../include/ksecforest_mba.h"

// a sample program
int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr,
            "Usage:\n"
            "  Build: ./main build graph.txt\n"
            "  Query: ./main query graph.idx kspan.idx kde.idx queries.txt k delta\n");
        return 1;
    }
    std::string mode = argv[1];

    /* ================= BUILD ================= */
    if (mode == "build") {
        if (argc < 3) {
            fprintf(stderr, "Build mode requires graph file\n");
            return 1;
        }

        std::string graph_file = argv[2];
        std::string op = argv[3];
        std::ifstream infile(graph_file);
        if (!infile) {
            perror("open graph file");
            return 1;
        }

        uint32_t t, n, m;
        infile >> t >> n >> m;
        infile.close();

        MbaF mbaF(n, m, graph_file);
        if (stoi(op) == 0) {
            auto beg = std::chrono::steady_clock::now();
            mbaF.twoPhase();
            auto end = std::chrono::steady_clock::now();
            printf("Index construction costs %.3f ms\n",
            std::chrono::duration<double, std::milli>(end - beg).count());
            std::string base = std::filesystem::path(graph_file).filename().string();
            mbaF.writeTGraph("./" + base);
            mbaF.writeKSECForestIndexToFile("./ksecforest_" + base);
        } else {
            auto beg = std::chrono::steady_clock::now();
            mbaF.onePass();
            auto end = std::chrono::steady_clock::now();
            printf("Index construction costs %.3f ms\n",
            std::chrono::duration<double, std::milli>(end - beg).count());
        }
    }

    /* ================= QUERY ================= */
    else if (mode == "query") {
        if (argc < 7) {
            fprintf(stderr,
                "Query mode requires 5 arguments\n");
            return 1;
        }

        std::string graph_file = argv[2];
        std::string ksecforest_file = argv[3];
        std::string queries_file = argv[4];
        int k = stoi(argv[5]);
        int delta = stoi(argv[6]);

        uint32_t n, m;
        std::ifstream infile(graph_file, std::ios::binary);
        infile.read((char*)&n, sizeof(n));
        infile.read((char*)&m, sizeof(m));
        infile.close();
        MbaF mbaF(n, m * 2);
        std::vector<int> query_vertices;
        std::ifstream qf(queries_file);
        int v;
        while (qf >> v) query_vertices.push_back(v);
        auto beg = std::chrono::steady_clock::now();
        for (int i = 0; i < query_vertices.size(); i++) {
            mbaF.searchCommunityByKSECForest(query_vertices[i], k, delta);
        }
        auto end = std::chrono::steady_clock::now();
        printf("KSpan query costs %.3f ms\n",
            std::chrono::duration<double, std::milli>(end - beg).count());
    }

    else {
        fprintf(stderr, "Unknown mode: %s\n", mode.c_str());
        return 1;
    }

    return 0;
}
