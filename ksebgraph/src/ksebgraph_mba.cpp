#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <utility>
#include <queue>
#include <unordered_set>
#include <chrono>
#include <random>
#include "../include/ksebgraph_mba.h"

#define ASSERT(truth) \
    if (!(truth)) { \
      std::cerr << "\x1b[1;31mASSERT\x1b[0m: " \
                << "LINE " << __LINE__ \
                << ", " << __FILE__ \
                << std::endl; \
      std::exit(EXIT_FAILURE); \
    } else

#define ASSERT_MSG(truth, msg) \
    if (!(truth)) { \
      std::cerr << "\x1b[1;31mASSERT\x1b[0m: " \
                << "LINE " << __LINE__ \
                << ", " << __FILE__ << '\n' \
                << "\x1b[1;32mINFO\x1b[0m: " << msg \
                << std::endl; \
      std::exit(EXIT_FAILURE); \
    } else

vector<unordered_map<uint32_t, uint32_t>> Mc;
vector<vector<Triangle>> time2triangle;

MbaG::MbaG(const uint32_t n, const uint32_t l)
    : n_(n), l_(l), tg_(n, l, trn_) {
    ASSERT_MSG(1 <= l_ && l_ < (static_cast<uint32_t>(1) << 29),
                "it is required 64 <= l <= 2^29 for the ease of implementation");
    trn_   = std::vector<uint32_t>(l_, 0);
}

MbaG::MbaG(const uint32_t n, const uint32_t l, const std::string& fn)
    : n_(n), l_(l), tg_(n, l, trn_) {
    ASSERT_MSG(1 <= l_ && l_ < (static_cast<uint32_t>(1) << 29),
                "it is required 64 <= l <= 2^29 for the ease of implementation");
    trn_   = std::vector<uint32_t>(l_, 0);
    loadGraph(fn);
}

void MbaG::loadGraph(const std::string& file_name) {
    std::ifstream infile(file_name, std::ios::in);
    ASSERT_MSG(infile.is_open(), "cannot open the file");
    // read the size of the graph
    uint32_t num_of_timestamp = 0, num_of_vertex = 0, temporal_edges = 0;
    infile >> num_of_timestamp >> num_of_vertex >> temporal_edges;
	tg_.ttshreld_ = num_of_timestamp; time2triangle.resize(tg_.ttshreld_ + 1);
    ASSERT(temporal_edges <= l_ && n_ == num_of_vertex);
    print_info(num_of_timestamp, num_of_vertex, temporal_edges);
    cout << endl << "Load Graph.." << endl;
    ASSERT_MSG(!infile.eof(), "invalid graph file");
    // read the edges
    map<EdgT,vector<uint32_t>> ET;
    for (int te = 0; te < temporal_edges; te++){
        uint32_t t, v1, v2;
        infile >> t >> v1 >> v2;
		if (t < 5) continue;
        if (v1 > v2) std::swap(v1, v2);
        ET[{v1,v2}].emplace_back(t);
    }
    infile.close();
    uint32_t eid = 0;
    for(auto& [edge, times] : ET){
		tg_.tau_[eid].assign(times.begin(), times.end());
        ASSERT(tg_.LazyInsert(edge.first, edge.second) == eid++);
    }
	m_ = tg_.m();
    Mc.resize(m_);
     // rectify the graph
    tg_.Rectify();
    // clear ET
    decltype(ET)().swap(ET);
    cout << "Load graph successfully!" << endl;
}

void MbaG::trussDecomp() {
 // truss decomposition
    // 1. compute the support of each edge by triangle listing
    // 1.1. define a total order over the vertices
    const auto pred = [this](const uint32_t v1, const uint32_t v2) {
    const size_t deg1 = tg_.adj_[v1].size();
    const size_t deg2 = tg_.adj_[v2].size();
        if (deg1 != deg2) return deg1 > deg2;
        else return v1 > v2;
    };
    VI verts(n_);
  	iota(verts.begin(), verts.end(), 0);
  	sort(verts.begin(), verts.end(), pred);
  	// 1.2. call the "forward" algorithm to list triangles
  	trn_.resize(m_, 0);
  	vector<vector<ArrayEntry>> A(n_);
    //delta-triangle list
	for (const uint32_t v : verts) {
    	for (const auto& ae : tg_.adj_[v]) {
			uint32_t u = ae.vid;
			uint32_t e = ae.eid;
			if (!pred(v, u)) continue;
			size_t pv = 0, pu = 0;
			while (pv < A[v].size() && pu < A[u].size()) {
				if (A[v][pv].vid == A[u][pu].vid) {
					++trn_[A[v][pv].eid]; ++trn_[A[u][pu].eid];
					++trn_[e];
					uint32_t interval = tg_.GetMst(tg_.tau_[e],tg_.tau_[A[v][pv].eid],tg_.tau_[A[u][pu].eid]);
					if(interval > tmax_) tmax_ = interval;
					time2triangle[interval].emplace_back(sortTri(e,A[v][pv].eid,A[u][pu].eid));
					++pv; ++pu;
				} else if (pred(A[v][pv].vid, A[u][pu].vid)) {
				++pv;
				} else {
				++pu;
				}
			}
      		A[u].emplace_back(ArrayEntry{v, e});
		}
	}
	// 2. decomposition
	// 2.1. sort the edges according to their supports
	const uint32_t max_sup = *max_element(trn_.cbegin(), trn_.cend());
	VI bin(max_sup + 1, 0);
  	for (uint32_t eid = 0; eid < m_; ++eid) ++bin[trn_[eid]];
	for (uint32_t i = 0, start = 0; i <= max_sup; ++i) {
		start += bin[i];
		bin[i] = start - bin[i];
	}

  	ord_.resize(m_, 0);
	pos.resize(m_, 0);
  	for (uint32_t eid = 0; eid < m_; ++eid) {
    	pos[eid] = bin[trn_[eid]];
    	ord_[pos[eid]] = eid;
    	++bin[trn_[eid]];
  	}
	rotate(bin.rbegin(), bin.rbegin() + 1, bin.rend());
  	bin[0] = 0;
  	// 2.2. peeling
  	ks_.resize(m_, 0);
	vector<bool> removed(m_, false);
  	uint32_t k = 0;
  	for (uint32_t i = 0; i < m_; ++i) {
    	k = max(k, trn_[ord_[i]]);
    	// ASSERT(bin[k] == i);
    	const uint32_t eid = ord_[i];
    	++bin[trn_[eid]];
    	removed[eid] = true;
    	// find triangles containing the edge with ID eid
    	const auto tris = tg_.GetTriangles(eid);
		// update ks_[eid]
		for (const auto& tri : tris) {
			const uint32_t e1 = tri.first;
			const uint32_t e2 = tri.second;
			if (trn_[e1] >= k && trn_[e2] >= k) ++ks_[eid];
			if (removed[e1] || removed[e2]) continue;
			for (const uint32_t e : {e1, e2}) {
				if (trn_[e] > k) {
				const uint32_t pe3 = bin[trn_[e]];
				const uint32_t pe = pos[e];
				if (pe3 != pe) {
					const uint32_t e3 = ord_[pe3];
					ord_[pe] = e3;
					pos[e3] = pe;
					ord_[pe3] = e;
					pos[e] = pe3;
				}
				++bin[trn_[e]];
				--trn_[e];
				}
			}
		}
  	}
	kmax_ = k;
}

void MbaG::KdeltaTrussDecomp() {
    trussDecomp();
	vector<uint32_t> trussness(m_);
	memcpy(trussness.data(),trn_.data(),trussness.size() * sizeof(uint32_t));
    queue<uint32_t> q;
	vector<bool> Ins(m_,false);
    kspan_.resize(kmax_ + 1, vector<uint32_t>(m_, UINT32_MAX));
	for(uint32_t t = tmax_ ; t > 0 ; --t){
		if(time2triangle[t].size() == 0) continue;
		for(auto [e1, e2, e3] : time2triangle[t]){
			Mc[e1][e3] = t;
            uint32_t k1 = trussness[e1], k2 = trussness[e2], k3 = trussness[e3];
			uint32_t mink = min(k1, min(k2, k3));
			//update k-support
			if(k1 == mink) { if(ks_[e1]-- == k1) { q.push(e1); Ins[e1] = true;} } 
			if(k2 == mink) { if(ks_[e2]-- == k2) { q.push(e2); Ins[e2] = true;} } 
			if(k3 == mink) { if(ks_[e3]-- == k3) { q.push(e3); Ins[e3] = true;} }
			//update the trussness of edge
			while(!q.empty()){
				uint32_t e = q.front(); q.pop();
                Ins[e] = false;
				uint32_t k = trussness[e];
				trussness[e]--; kspan_[k][e] = t;
				if(trussness[e] == 0) continue;
				uint32_t kse = 0;
				const auto tris = tg_.GetTriangles(e);
				for (const auto [e1, e2] : tris) {
					if (trussness[e1] < k - 1 || trussness [e2] < k - 1) continue;
					getminmax(e, e1, e2);
					auto pos = Mc[mine].find(maxe);
					if(pos == Mc[mine].end()) {
						kse++;
						if (trussness[e1] >= k && trussness[e2] >= k) {
							if(trussness[e1] == k && Ins[e1] == false && ks_[e1]-- == trussness[e1]) {
								q.push(e1); Ins[e1] = true;
							}
							if(trussness[e2] == k && Ins[e2] == false && ks_[e2]-- == trussness[e2]) {
								q.push(e2); Ins[e2] = true;
							}	
						}
					}
				}
				ks_[e] = kse;
			}
		}//for time2triangle[t]
	}
    for(uint32_t e = 0; e < m_; ++e){
		while(trussness[e] > 0){
			kspan_[trussness[e]][e] = 0;
			trussness[e]--; 
		}  
    }
	return;
}

void MbaG::twoPhase() {
	KdeltaTrussDecomp();
	initKSEBGraphs(kmax_ + 1, m_);
	for (int k = 1; k <= kmax_; k++) {
		constructKSEBGraphForK(k);
	}
	return;
}

void MbaG::onePass() {
    trussDecomp();
	initKSEBGraphs(kmax_ + 1, m_);
	vector<vector<unordered_map<int, int>>> edgeMap(kmax_ + 1, vector<unordered_map<int, int>>(m_)); 
    vector<vector<tuple<int, int, int>>> missedgeMap (kmax_ + 1, vector<tuple<int, int, int>>(m_, {-1, -1, -1})); 
    vector<vector<int>> RidtoBid(kmax_ + 1, vector<int>(m_)); 
	vector<bool> edgeProcessed(m_, false);
	vector<vector<vector<pair<int, int>>>> superE(kmax_ + 1, vector<vector<pair<int, int>>>(tmax_ + 1));
	vector<vector<vector<int>>> superN(kmax_ + 1);
	vector<vector<pair<int, vector<int>>>> R(kmax_ + 1);
	queue<int> q;
    vector<int> CID(kmax_ + 1, -1); 
	vector<int> trussness(m_);
	memcpy(trussness.data(),trn_.data(),trussness.size() * sizeof(int));

	auto tryInsert = [](auto& mp, const auto& key, const auto& value) {
		auto [it, inserted] = mp.try_emplace(key, value);
		if (!inserted) {
			if (value < it->second) it->second = value;
		}
	};
	for(int delta = tmax_ ; delta > 0 ; --delta){
		if(time2triangle[delta].size() == 0) continue;
		for(auto [ea, eb, ec] : time2triangle[delta]){
			Mc[ea][ec] = delta;
            int ka = trussness[ea], kb = trussness[eb], kc = trussness[ec];
			int mink = min(ka, min(kb, kc));
			//update k-support
			if(ka == mink) { if(ks_[ea]-- == ka) { q.push(ea); edgeProcessed[ea] = true;} } 
			if(kb == mink) { if(ks_[eb]-- == kb) { q.push(eb); edgeProcessed[eb] = true;} } 
			if(kc == mink) { if(ks_[ec]-- == kc) { q.push(ec); edgeProcessed[ec] = true;} }
			if (q.empty()) continue;
			// generate a node for mink-span equal to t
			int curid = ++CID[mink];
			superN[mink].emplace_back(vector<int>{curid});
			RidtoBid[mink][curid] = curid;
			auto& RE = R[mink].emplace_back(delta, vector<int>{}).second;
			//update the trussness of edge
			while(!q.empty()){
				int e = q.front(); q.pop();
				int k = trussness[e], kse = 0;
				RE.emplace_back(e);
                edgeProcessed[e] = false;
                int missede = get<0>(missedgeMap[k][e]);
                if (missede != -1 && trussness[missede] >= k) {
					int missedt = get<1>(missedgeMap[k][e]), missedsp = get<2>(missedgeMap[k][e]);
					if (missedsp == delta) tryInsert(edgeMap[k][missede], curid, missedt);
                }
				// check the super edges among super nodes
				if (!edgeMap[k][e].empty()) {
					//preid: the id of super nodes with kspan no less than t, delta: the weight of two super nodes
                    for (auto [exploredRid, t] : edgeMap[k][e]) {
                        if (exploredRid == curid) continue; 	   
                        int RKspan = R[k][exploredRid].first;
						if (RKspan == delta && t <= RKspan) {
							int preBid = RidtoBid[k][exploredRid], curBid = RidtoBid[k][curid]; 
							if (preBid == curBid) continue;
							if (superN[k][preBid].size() <= superN[k][curBid].size()) {
								for (int rid : superN[k][preBid]) {
									superN[k][curBid].push_back(rid);
									RidtoBid[k][rid] = curBid;
								}
								superN[k][preBid].clear();
							} else {
								for (int rid : superN[k][curBid]) {
									superN[k][preBid].push_back(rid);
									RidtoBid[k][rid] = preBid;
								}
								superN[k][curBid].clear();
							}
						} else {
							superE[k][t].push_back({exploredRid, curid});
						}
                    }
                }
				//check the k-triangles containing e
				const auto tris = tg_.GetTriangles(e);
				for (const auto [e1, e2] : tris) {
                    if (trussness[e1] < k - 1 || trussness[e2] < k - 1) continue;
					getminmax(e, e1, e2);
					auto pos = Mc[mine].find(maxe);
					if(pos == Mc[mine].end()) {
                        int mst = tg_.GetMst(tg_.tau_[e], tg_.tau_[e1], tg_.tau_[e2]);
						if (trussness[e1] >= k && trussness[e2] >= k) {
							if (!edgeProcessed[e1]) {
								if (trussness[e1] == k && ks_[e1]-- == k) {
									q.push(e1); edgeProcessed[e1] = true;
								} else {
									tryInsert(edgeMap[k][e1], curid, mst);
								}
							} 
							if (!edgeProcessed[e2]) {
								if (trussness[e2] == k && ks_[e2]-- == k) {
									q.push(e2); edgeProcessed[e2] = true;
								} else {
									tryInsert(edgeMap[k][e2], curid, mst);
								}
							} 
						}
						kse++;
					} else {
                        // if (pos->second == delta) continue;
						if (trussness[e1] >= k && trussness[e2] >= k) {
							if (!edgeProcessed[e1]) { tryInsert(edgeMap[k][e1], curid, pos->second);}
							if (!edgeProcessed[e2]) { tryInsert(edgeMap[k][e2], curid, pos->second);}
                            if (pos->second == delta) continue;
                            missedgeMap[k][e1] = {e2, pos->second, delta};
                            missedgeMap[k][e2] = {e1, pos->second, delta};
						}
					}
				}
				ks_[e] = kse;
				trussness[e]--;
			}
		}//for time2triangle[t]	
	}
	// deal with kspan equal to 0 specially
	vector<vector<int>> kspan_dict(kmax_ + 1);  
    for(int e = 0; e < m_; ++e){
		if (trussness[e] > 0) kspan_dict[trussness[e]].push_back(e);
    }
	vector<int> tmp;
	int tmpSize = 0;
	for (int k = kmax_; k >= 1; k--) {
		auto& graph  = ksebgraph->graphs[k];
		tmpSize += kspan_dict[k].size();
		tmp.reserve(tmpSize); tmp.insert(tmp.end(),  std::make_move_iterator(kspan_dict[k].begin()), std::make_move_iterator(kspan_dict[k].end()));
		int i = 0;
		while (i < tmpSize){
			int e0 = tmp[i]; i++;
			if (edgeProcessed[e0]) continue;
 			CID[k]++; 
			int curid = CID[k];
			superN[k].emplace_back(vector<int>{curid});
			RidtoBid[k][curid] = curid;
			auto& eR = R[k].emplace_back(0, vector<int>{}).second;
			q.push(e0);
			edgeProcessed[e0] = true;
			while (!q.empty()) {
                int e = q.front(); q.pop();
                eR.emplace_back(e);
                if (!edgeMap[k][e].empty()) {
                    for (const auto [exploredRId, t] : edgeMap[k][e]) {
                        if (exploredRId == curid) continue; 	   
						superE[k][t].push_back({exploredRId, curid});
                    }
                } 
				const auto  tris = tg_.GetTriangles(e);
				for (const auto& [e1, e2] : tris) {
					if (trussness[e1] < k || trussness[e2] < k) continue;
					getminmax(e, e1, e2);
					auto pos = Mc[mine].find(maxe);
					if(pos == Mc[mine].end()) {
						if (!edgeProcessed[e1]) { q.push(e1); edgeProcessed[e1] = true; }
						if (!edgeProcessed[e2]) { q.push(e2); edgeProcessed[e2] = true; }
					} else {
						if (!edgeProcessed[e1]) { tryInsert(edgeMap[k][e1], curid, pos->second);} 
						if (!edgeProcessed[e2]) { tryInsert(edgeMap[k][e2], curid, pos->second);}
					}
				}
			}
		}
		fill(edgeProcessed.begin(), edgeProcessed.end(), false);
	}
	// generategraphs for each k
	for (int k = 1; k <= kmax_; k++) {
		auto& graph  = ksebgraph->graphs[k];
		auto& snk = superN[k];
		auto& Rk = R[k];
		auto& Mk = RidtoBid[k];
		int snkSize = snk.size();
		// the position of blocks in super node[k] 
		vector<int> bidPos;
		int& Blockid = graph.maxBlockId;
		for (int idx = 0; idx < snkSize; idx++) {
			if (snk[idx].size() == 0) continue;
			for (int Rid : snk[idx]) {
				Mk[Rid] = Blockid;
			}
			bidPos.push_back(idx);
			int Bkspan= R[k][snk[idx][0]].first;
            graph.idGN.emplace(Blockid, GN(Bkspan, Blockid));
			Blockid++;
		}
        EdgeHash seenSuperE;
		vector<vector<pair<int, int>>> blockAdj(Blockid);
        auto& SE = superE[k];
		for (int t = 0; t <= tmax_; t++) {
            for (auto [x, y] : SE[t]) {
                int bidx = Mk[x], bidy = Mk[y];
                if (bidx == bidy) continue;
				if (tryInsertEdge(seenSuperE, x, y)) {
					blockAdj[bidx].push_back({bidy,t});
					blockAdj[bidy].push_back({bidx,t});
				}
            }		
		} 
		
		for(auto& [bid, gn] : graph.idGN) {
			size_t edge_cnt = 0;
			int pos = bidPos[bid];
            for (int rid : snk[pos]) edge_cnt += Rk[rid].second.size();
			gn.elements_.reserve(edge_cnt);
            for (int rid : snk[pos]) {
                gn.elements_.insert(gn.elements_.end(),  
                std::make_move_iterator(Rk[rid].second.begin()), 
                std::make_move_iterator(Rk[rid].second.end()));
            }
			gn.adj = std::move(blockAdj[bid]);
			for (int e : gn.elements_) { graph.etoGraphNID[e] = bid; }
		}
	}
	return;
}

void MbaG::constructKSEBGraphForK(int k) {
	vector<vector<int>> kspan_dict(tmax_ + 1);
	getKSpanDict(k, kspan_dict);
	vector<unordered_map<int, int>> edgeMap(m_); 
	vector<vector<pair<int, int>>> superE(tmax_ + 1);
    vector<bool> edgeProcessed(m_);
    queue<int> q;
	auto& graph = ksebgraph->graphs[k];
	auto& ksp = kspan_[k];
	auto& etoBlockId = graph.etoGraphNID;
	int& graphNodeId = graph.maxBlockId;
	for (int delta = 0; delta <= tmax_; delta++) {
		if (kspan_dict[delta].size() == 0) continue;
		for (int e0 : kspan_dict[delta]) {
			if (edgeProcessed[e0]) continue;
			edgeProcessed[e0] = true;
			auto [it, inserted] = graph.idGN.emplace(graphNodeId, GN(delta, graphNodeId));
			KSEB& curKseb = it->second.elements_;
			q.push(e0);
			while (!q.empty()){
				int e = q.front(); q.pop();
				curKseb.emplace_back(e);
                etoBlockId[e] = graphNodeId;
				if (!edgeMap[e].empty()) {
					for (auto [exploredGraphId, t] : edgeMap[e]) {
						if (exploredGraphId == graphNodeId) continue;
						superE[t].push_back({exploredGraphId, graphNodeId});
					}
				}
				// find triangles containing the edge with ID eid
				const auto tris = tg_.GetTriangles(e, k);
				for (const auto [e1, e2] : tris) {
                    if (edgeProcessed[e1] && edgeProcessed[e2]) continue;
					int sp1 = ksp[e1], sp2 = ksp[e2], triMst = 0;
					getminmax(e, e1, e2); 
					if (auto pos = Mc[mine].find(maxe); pos != Mc[mine].end()) {
						triMst = pos->second;
					}
					int maxsp = max(sp1, sp2);
					// function to process_edge
					// process edge e1
                    if (!edgeProcessed[e1] && sp1 == maxsp) {
                        if (sp1 == delta && triMst <= delta) {
                            q.push(e1); edgeProcessed[e1] = true;
                        } else {
                            auto [it, inserted] = edgeMap[e1].try_emplace(graphNodeId, triMst);
                            if (!inserted) {
                                if (triMst < it->second) it->second = triMst;
                            } 
                        }
                    }
                      if (!edgeProcessed[e2] && sp2 == maxsp) {
                        if (sp2 == delta && triMst <= delta) {
                            q.push(e2); edgeProcessed[e2] = true;
                        } else {
                            auto [it, inserted] = edgeMap[e2].try_emplace(graphNodeId, triMst);
                            if (!inserted) {
                                if (triMst < it->second) it->second = triMst;
                            } 
                        }
                    }
				}
			}
			graphNodeId++;
		}
	}

	EdgeHash seenSuperE;
	vector<vector<pair<int, int>>> blockAdj(graphNodeId);
    for (int t = 0; t <= tmax_; t++) {
        for (auto [x, y] : superE[t]) {
            if (tryInsertEdge(seenSuperE, x, y)) {
				blockAdj[x].push_back({y,t});
				blockAdj[y].push_back({x,t});
            }
        }		
    } 
	for (auto& [id, gn] : graph.idGN) {
		gn.adj = std::move(blockAdj[id]);
	}
	return;
}

vector<vector<int>> MbaG::searchCommunityByKSEBraphs(int query, int k, uint32_t delta) {
	Graph& graph = ksebgraph->graphs[k];
    vector<bool> Ins(graph.maxBlockId, false);
    queue<GN*> q;
    vector<vector<int>> result;
    for (const auto& ae : tg_.adj_[query]) {
        int e = ae.eid;
        int gnid = graph.etoGraphNID[e];
        if (gnid == -1 || Ins[gnid]) continue;
        GN* start = &graph.idGN[gnid];
        if (start->delta_ > delta) continue;
        vector<int> Ai;
        Ins[start->id_] = true;
        q.push(start);
        while (!q.empty()) {
            GN* cur = q.front(); q.pop();
            Ai.insert(Ai.end(), cur->elements_.begin(), cur->elements_.end());
            for (auto [adjId, edgeWeight] : cur->adj) {
                if (Ins[adjId] || edgeWeight > delta) continue;
                GN* adjNode = &graph.idGN[adjId];
                if (adjNode->delta_ > delta) continue;
                Ins[adjNode->id_] = true;
                q.push(adjNode);
            }
        }
        result.emplace_back(std::move(Ai));
    }
    // int size = 0;
    // for (auto& community : result) {
	// 	size += community.size();
	// 	for (int& e : community) {
	// 		cout<<e<<" "<<kspan_[k][e]<<endl;
	// 	}
	// }
	// cout<<size<<endl;
    return result;
}

void MbaG::writeTGraph(const std::string& filename) {
	std::ofstream outfile(filename, std::ios::binary);
    if (!outfile) {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
        return;
    }
	double fsz = 0.0;
	/* ---------- write temporal graph ---------- */
    outfile.write((char*)&n_, sizeof(n_));
    outfile.write((char*)&m_, sizeof(m_));
	fsz += sizeof(n_);
    fsz += sizeof(m_);

    for (uint32_t i = 0; i < m_; i++) {
        uint32_t u = tg_.edge_info_[i].first;
        uint32_t v = tg_.edge_info_[i].second;
        uint32_t tsize = tg_.tau_[i].size();

        outfile.write((char*)&u, sizeof(u));
        outfile.write((char*)&v, sizeof(v));
        outfile.write((char*)&tsize, sizeof(tsize));
        outfile.write((char*)tg_.tau_[i].data(), tsize * sizeof(uint32_t));
		fsz += sizeof(uint32_t) * 3;               // u, v, tsize
        fsz += tsize * sizeof(uint32_t);
    }
	outfile.close();
	std::cout << "Index successfully written to: " << filename << std::endl;
	double mb = fsz / (1024.0 * 1024.0);
	printf("Index size: %.2f MB\n", mb);
	return;
}


void MbaG::writeKSEBGraphsIndexToFile(const std::string& filename) {
    std::ofstream outfile(filename, std::ios::binary | std::ios::trunc);
    if (!outfile) {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
        return;
    }

    double fsz = 0.0;

    // ---- kmax ----
    outfile.write(reinterpret_cast<const char*>(&kmax_), sizeof(kmax_));
    fsz += sizeof(kmax_);

    // ---- graphs ----
    for (uint32_t k = 1; k <= kmax_; ++k) {
        auto& graph = ksebgraph->graphs[k];

        outfile.write(reinterpret_cast<const char*>(&graph.maxBlockId),
                      sizeof(graph.maxBlockId));
        fsz += sizeof(graph.maxBlockId);

        int gn_cnt = static_cast<int>(graph.idGN.size());
        outfile.write(reinterpret_cast<const char*>(&gn_cnt), sizeof(gn_cnt));
        fsz += sizeof(gn_cnt);

        for (auto& [gid, gn] : graph.idGN) {
            outfile.write(reinterpret_cast<const char*>(&gid), sizeof(gid));
            fsz += sizeof(gid);

            outfile.write(reinterpret_cast<const char*>(&gn.delta_),
                          sizeof(gn.delta_));
            fsz += sizeof(gn.delta_);

            int adj_cnt = static_cast<int>(gn.adj.size());
            outfile.write(reinterpret_cast<const char*>(&adj_cnt),
                          sizeof(adj_cnt));
            fsz += sizeof(adj_cnt);

            // ---- SAFE adj write ----
            for (auto& [to, w] : gn.adj) {
                outfile.write(reinterpret_cast<const char*>(&to), sizeof(int));
                outfile.write(reinterpret_cast<const char*>(&w), sizeof(int));
                fsz += sizeof(int) * 2;
            }

            int edge_cnt = static_cast<int>(gn.elements_.size());
            outfile.write(reinterpret_cast<const char*>(&edge_cnt),
                          sizeof(edge_cnt));
            fsz += sizeof(edge_cnt);

            outfile.write(reinterpret_cast<const char*>(gn.elements_.data()),
                          edge_cnt * sizeof(int));
            fsz += edge_cnt * sizeof(int);
        }
    }

    // ---- etoGraphNID ----
    for (int e = 0; e < m_; ++e) {
        int maxk_to_e = trn_[e];
        outfile.write(reinterpret_cast<const char*>(&maxk_to_e),
                      sizeof(maxk_to_e));
        fsz += sizeof(maxk_to_e);

        for (int k = 1; k <= maxk_to_e; ++k) {
            int gid = ksebgraph->graphs[k].etoGraphNID[e];
            outfile.write(reinterpret_cast<const char*>(&gid), sizeof(gid));
            fsz += sizeof(gid);
        }
    }

    outfile.close();

    std::cout << "Index successfully written to: " << filename << std::endl;
    printf("Index size: %.2f MB\n", fsz / (1024.0 * 1024.0));
}

void MbaG::readTGraphFromFile(const std::string& filename) {
	std::ifstream infile(filename, std::ios::binary);
    if (!infile) {
        std::cerr << "Failed to open file for reading: " << filename << std::endl;
        return;
    }
	/* -------- read temporal graph -------- */
	uint32_t n = 0, m = 0;
    infile.read((char*)&n, sizeof(n));
    infile.read((char*)&m, sizeof(m));
	// ASSERT(m <= l_ && n_ == n);
	for (uint32_t e = 0; e < m; ++e) {
      	uint32_t u, v, tsize;
		infile.read((char*)&u, sizeof(u));
		infile.read((char*)&v, sizeof(v));
		tg_.LazyInsert(u, v);
		infile.read((char*)&tsize, sizeof(tsize));
        tg_.tau_[e].resize(tsize);
		infile.read((char*)tg_.tau_[e].data(),
                    tsize * sizeof(uint32_t));
    } 
	m_ = tg_.m_;
	infile.close();
	//  std::cout << "Index successfully read from: " << filename << std::endl;
	return;
}



void MbaG::readKSEBGraphsIndexFromFile(const std::string& filename) {
    std::ifstream infile(filename, std::ios::binary);
    if (!infile) {
        std::cerr << "Failed to open file for reading: " << filename << std::endl;
        return;
    }

    // ---- kmax ----
    infile.read(reinterpret_cast<char*>(&kmax_), sizeof(kmax_));

    initKSEBGraphs(kmax_ + 1, m_);  // 你已有的初始化函数

    // ---- graphs ----
    for (uint32_t k = 1; k <= kmax_; ++k) {
        auto& graph = ksebgraph->graphs[k];

        // maxBlockId
        infile.read(reinterpret_cast<char*>(&graph.maxBlockId),
                    sizeof(graph.maxBlockId));

        // GN count
        int gn_cnt;
        infile.read(reinterpret_cast<char*>(&gn_cnt), sizeof(gn_cnt));

        graph.idGN.reserve(gn_cnt);

        for (int i = 0; i < gn_cnt; ++i) {
            int gid;
            infile.read(reinterpret_cast<char*>(&gid), sizeof(gid));

            GN gn;
            gn.id_ = gid;

            infile.read(reinterpret_cast<char*>(&gn.delta_),
                        sizeof(gn.delta_));

            // ---- adj ----
            int adj_cnt;
            infile.read(reinterpret_cast<char*>(&adj_cnt), sizeof(adj_cnt));
            gn.adj.resize(adj_cnt);

            for (int j = 0; j < adj_cnt; ++j) {
                int to, w;
                infile.read(reinterpret_cast<char*>(&to), sizeof(int));
                infile.read(reinterpret_cast<char*>(&w), sizeof(int));
                gn.adj[j] = {to, w};
            }

            // ---- elements ----
            int edge_cnt;
            infile.read(reinterpret_cast<char*>(&edge_cnt), sizeof(edge_cnt));
            gn.elements_.resize(edge_cnt);

            infile.read(reinterpret_cast<char*>(gn.elements_.data()),
                        edge_cnt * sizeof(int));

            graph.idGN.emplace(gid, std::move(gn));
        }
    }

    // ---- etoGraphNID ----
    for (int e = 0; e < m_; ++e) {
        int maxk_to_e;
        infile.read(reinterpret_cast<char*>(&maxk_to_e),
                    sizeof(maxk_to_e));

        for (int k = 1; k <= maxk_to_e; ++k) {
            int gid;
            infile.read(reinterpret_cast<char*>(&gid), sizeof(gid));
            ksebgraph->graphs[k].etoGraphNID[e] = gid;
        }
    }

    infile.close();
    // std::cout << "Index successfully read from: " << filename << std::endl;
}
