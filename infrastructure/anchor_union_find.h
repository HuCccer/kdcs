#include <vector>
#include <algorithm>
using namespace std;

class AnchorUnionFind {
public:
    vector<int> parent_;  
    vector<int> rank_;   
    vector<int> anchor_; 

public:
    AnchorUnionFind(int m) {
        parent_.resize(m);
        rank_.resize(m, 0);
        anchor_.resize(m);
        for (int i = 0; i < m; i++) {
            parent_[i] = i; anchor_[i] = i;
        }
    }

    int find(int x) {
        int root = x;
        while (root != parent_[root]) {
            root = parent_[root];
        }
        while (x != root) {
            int px = parent_[x];
            parent_[x] = root;
            x = px;
        }
        return root;
    }

    void update(int x, int y, int anchor) {
        int px = find(x), py = find(y);
        if (px == py) return;
        int new_root = -1;
        if (rank_[px] < rank_[py]) {
            parent_[px] = py;
            new_root = py;
        } else {
            if (rank_[px] == rank_[py]) rank_[px]++;
            parent_[py] = px;
            new_root = px;
        }
        anchor_[new_root] = anchor;
    }

    void reset() {
        int m = parent_.size();
        for (int i = 0; i < m; i++) {
            parent_[i] = i; anchor_[i] = i; 
        }
        fill(rank_.begin(), rank_.end(), 0);
    }

};
