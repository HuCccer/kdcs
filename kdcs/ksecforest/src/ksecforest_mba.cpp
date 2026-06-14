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
#include "../include/ksecforest_mba.h"

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

using namespace std;
using namespace KDeltaTrussCommunity;

MbaF::MbaF(uint32_t n, uint32_t l)
    : n_(n), l_(l), k_(l), tg_(n, l, k_), kdt_(n, l, tg_, k_, triMTSpan){
    ASSERT_MSG(1 <= l_ && l_ < (static_cast<uint32_t>(1) << 29),
                "it is required 64 <= l <= 2^29 for the ease of implementation");
}

MbaF::MbaF(uint32_t n, uint32_t l, const string& fn)
    : n_(n), l_(l), k_(l), tg_(n, l, k_), kdt_(n, l, tg_, k_, triMTSpan){
    ASSERT_MSG(1 <= l_ && l_ < (static_cast<uint32_t>(1) << 29),
                "it is required 64 <= l <= 2^29 for the ease of implementation");
	loadGraph(fn);
}

void MbaF::loadGraph(const std::string& file_name) {
    std::ifstream infile(file_name, std::ios::in);
    ASSERT_MSG(infile.is_open(), "cannot open the file");
    // read the size of the graph
    uint32_t num_of_timestamp = 0, num_of_vertex = 0, temporal_edges = 0;
    infile >> num_of_timestamp >> num_of_vertex >> temporal_edges;
    ASSERT(temporal_edges <= l_ && n_ == num_of_vertex);
	tg_.maxtimestamps_ = num_of_timestamp;
    print_info(num_of_timestamp, num_of_vertex, temporal_edges);
    cout << endl << "Load Graph.." << endl;
    ASSERT_MSG(!infile.eof(), "invalid graph file");
    // read the edges
    map<EdgT,vector<uint32_t>> ET;
    for (uint32_t te = 0; te < temporal_edges; te++){
        uint32_t t, v1, v2;
        infile >> t >> v1 >> v2;
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
     // rectify the graph
    tg_.Rectify();
    // clear ET
    decltype(ET)().swap(ET);
    cout << "Load graph successfully!" << endl;
}


void MbaF::Serialmerge(uint32_t seed, uint32_t k, unordered_set<uint32_t>& treeNodeToInsert) {
	auto& tn = ksecforest->forests[k].idTN;
	auto comparer = [&](uint32_t i, uint32_t j) {
        return tn.at(i).delta_ > tn.at(j).delta_;
    };
	priority_queue<uint32_t, vector<uint32_t>, decltype(comparer)> pq(comparer);
	unordered_set<uint32_t> inpq;
	unordered_set<uint32_t> doneNode;
	pq.push(seed); inpq.insert(seed);
	for (uint32_t nodeId : treeNodeToInsert) {
		pq.push(nodeId); inpq.insert(nodeId);
	}
	uint32_t lastNodeId = UINT32_MAX;
	while (!pq.empty()) {
		uint32_t currMinDelta = tn.at(pq.top()).delta_;
		vector<uint32_t> toMerge;
		while (!pq.empty() && tn.at(pq.top()).delta_ == currMinDelta) {
			uint32_t nodeId = pq.top(); pq.pop();
			toMerge.push_back(nodeId);
			uint32_t fatherId = tn.at(nodeId).father_;
			if (fatherId != UINT32_MAX && !inpq.count(fatherId)) {
				inpq.insert(fatherId);
				pq.push(fatherId);
			}
		}
		uint32_t id = sameKSpanMerging(toMerge, k);
		auto& children = tn.at(id).children_;
		uint32_t i = 0;
		while (i < children.size()) {
			uint32_t childId = children[i];
			if (doneNode.count(childId) || !tn.count(childId)) {
				swap(children[i], children.back());
				children.pop_back();
			} else {
				i++;
			}
		}
		if (lastNodeId != UINT32_MAX) {
			if (tn.at(id).elements_.size() == 0 && children.size() == 0) {
				tn.erase(id); doneNode.insert(id);
				tn.at(lastNodeId).father_ = UINT32_MAX;
				continue;
			}
			tn.at(lastNodeId).father_ = id;
			children.push_back(lastNodeId);
		} 
		doneNode.insert(id);
		lastNodeId = id;
	}
    return;
}


uint32_t MbaF::sameKSpanMerging(const std::vector<uint32_t> &toMerge, uint32_t k) {
    if (toMerge.size() < 2) {
        return toMerge.empty() ? UINT32_MAX : toMerge[0];
    }
	auto& tn = ksecforest->forests[k].idTN;
	auto& etoTN = ksecforest->forests[k].etoTreeNID;
	uint32_t finalNodeId = toMerge[0];
    auto& finalNode = tn.at(finalNodeId);
    for (uint32_t i = 1; i < toMerge.size(); i++) {
		auto& nodei = tn.at(toMerge[i]);
		for (uint32_t tid : nodei.children_) {
			if (tn.count(tid) == 0) continue;
			finalNode.children_.push_back(tid);
			tn.at(tid).father_ = finalNodeId;
		}
		auto& edges = nodei.elements_; // edge IDs contained in the node to be merged
		copy(edges.begin(), edges.end(), back_inserter(finalNode.elements_));
		for (uint32_t eid : edges) {
			etoTN[eid] = finalNodeId;
		}
		tn.erase(toMerge[i]);
    }
    return finalNodeId;
}


 void MbaF::addEdge(uint32_t u, uint32_t v, uint32_t timestamp) {
	bool variety = false;
	uint32_t e0 = tg_.InsertEdge(u, v, variety);
	vector<uint32_t> oldtaue0 = tg_.tau_[e0];
	auto aff = kdt_.addEdgeWithVariety(e0, timestamp, variety); 
	uint32_t maxk = aff.size() - 1;
	if (maxk == 0) return;
	if (maxk > kmax_) {
		updateKSECForest(); kmax_ = maxk;
	} 
	for (uint32_t k = 1; k <= maxk; k++) {
		// cout<<"k"<<k<<endl;
		auto& forest = ksecforest->forests[k];
		auto& kspan = kdt_.span_[k];
		auto& etoTID = forest.etoTreeNID;
		auto& idtn = forest.idTN;
		map<uint32_t, vector<uint32_t>> splitaff;
		unordered_map<uint32_t, unordered_set<uint32_t>> origNodes;
		// collect the original nodes contained in the super nodes affected by the insertion of the new edge
		for (uint32_t e : aff[k]) {		
			uint32_t ksp = kspan[e]; uint32_t tid = etoTID[e];
			splitaff[ksp].push_back(e);
			if (tid != UINT32_MAX) {
				origNodes[tid].insert(e);
			}
		}
		// delete edges with kspan increased to ksp
		for (auto& [tid, edges] : origNodes) {
			auto& nodeEdges = idtn[tid].elements_;
			// uint32_t delta =  idtn[tid].delta_, size = nodeEdges.size();
			for (uint32_t i = 0; i < nodeEdges.size();) {
				uint32_t e = nodeEdges[i];
				if (edges.count(e)) {
					swap(nodeEdges[i], nodeEdges.back()); nodeEdges.pop_back();
					etoTID[e] = UINT32_MAX;
				} else {
					i++;
				}
			}
			if (nodeEdges.empty()) {
				uint32_t fatherId = idtn[tid].father_;
				if (fatherId == UINT32_MAX) {
					for(uint32_t childId : idtn[tid].children_) idtn[childId].father_ = UINT32_MAX;
				} else {
					auto& children = idtn[fatherId].children_;
					children.erase(remove(children.begin(), children.end(), tid), idtn[fatherId].children_.end());
					for (uint32_t childId : idtn[tid].children_) {
						idtn[childId].father_ = fatherId; idtn[fatherId].children_.push_back(childId);
					}
				}
				idtn.erase(tid);
			}
		}
		unordered_map<uint32_t, uint32_t> treeNode;
		unordered_set<uint32_t> treeNodeToInsert;
		for (auto& [ksp, edges] : splitaff) {
			treeNode.clear(); treeNodeToInsert.clear();
			uint32_t newNodeId = forest.maxBlockId++;
			idtn.emplace(newNodeId, TN(ksp, newNodeId));
			idtn[newNodeId].elements_ = edges;
			// update the mapping from edge ID to tree node ID
			for (uint32_t e : edges) etoTID[e] = newNodeId;
			vector<uint32_t> q;
			unordered_set<uint32_t> seen;
			q.insert(q.end(), edges.begin(), edges.end());
			while (!q.empty()){
				uint32_t e = q.back(); q.pop_back();
				seen.insert(e);
				tg_.ForEachTriangle(e, k, [&](uint32_t e1, uint32_t e2) {
					uint32_t tid1 = etoTID[e1], tid2 = etoTID[e2];
					if (tid1 == UINT32_MAX || tid2 == UINT32_MAX || (tid1 == newNodeId && tid2 == newNodeId)) return;
					if (seen.count(e1) || seen.count(e2)) return;
					uint32_t sp1 = kspan[e1], sp2 = kspan[e2], spe = kspan[e];
					uint32_t triMst = tg_.GetMts(tg_.tau_[e], tg_.tau_[e1], tg_.tau_[e2]);
					uint32_t maxsp = max(sp1, max(sp2, spe));
					uint32_t weight = max(maxsp, triMst);
					if (tid1 != newNodeId) {
						auto [it, inserted] = treeNode.try_emplace(tid1, weight);
						if (!inserted) if (weight < it->second) it->second = weight;
					}
					if (tid2 != newNodeId) {
						auto [it, inserted] = treeNode.try_emplace(tid2, weight);
						if (!inserted) if (weight < it->second) it->second = weight;
					}
				});
			};
			uint32_t newnodekspan = idtn[newNodeId].delta_;
			for (auto& [tid, weight] : treeNode) {
				uint32_t delta = idtn[tid].delta_;
				if (delta < weight) {
					auto [ancestor, fatherofancestor] = findAncestorAndFatherLessD(idtn, tid, weight);
					if (fatherofancestor == newNodeId || ancestor == newNodeId) continue;
					if (weight == newnodekspan) {
						if (fatherofancestor != UINT32_MAX) {
							delSuperEdge(idtn[fatherofancestor], idtn[ancestor]);
						} 
						addSuperEdge(idtn[newNodeId], idtn[ancestor]);
						if (fatherofancestor != UINT32_MAX) treeNodeToInsert.insert(fatherofancestor);
					} else{
						if (treeNodeToInsert.count(ancestor)) continue;
						if (fatherofancestor == UINT32_MAX) {
							uint32_t newemptyNodeId = forest.maxBlockId++;
							auto [tn, inserted]	= idtn.emplace(newemptyNodeId, TN(weight, newemptyNodeId));
							// add super-edge between ancestor and new emptyNodeId
							tn->second.children_.push_back(ancestor);
							idtn[ancestor].father_ = newemptyNodeId;
							treeNodeToInsert.insert(newemptyNodeId);
						} else {
							if (idtn[fatherofancestor].delta_ != weight) {
								uint32_t newemptyNodeId = forest.maxBlockId++;
								auto [tn, inserted]	= idtn.emplace(newemptyNodeId, TN(weight, newemptyNodeId));
								delSuperEdge(idtn[fatherofancestor], idtn[ancestor]);
								tn->second.children_.push_back(ancestor);
								idtn[ancestor].father_ = newemptyNodeId;
								treeNodeToInsert.insert(newemptyNodeId);
								treeNodeToInsert.insert(fatherofancestor);
							} else {
								treeNodeToInsert.insert(fatherofancestor);
							}
						}
					}
				} else{
					treeNodeToInsert.insert(tid);
				} 
			}
			Serialmerge(newNodeId, k, treeNodeToInsert);
		}
		if (aff[k].count(e0)) continue;
		treeNode.clear(); treeNodeToInsert.clear();
		uint32_t e0TID = etoTID[e0];
		tg_.ForEachTriangle(e0, k, [&](uint32_t e1, uint32_t e2) {
			if (aff[k].count(e1) || aff[k].count(e2)) return;
			uint32_t oldmts = tg_.GetMts(oldtaue0, tg_.tau_[e1], tg_.tau_[e2]);
			uint32_t newMst = tg_.GetMts(tg_.tau_[e0], tg_.tau_[e1], tg_.tau_[e2]);
			if (oldmts == newMst) return;
			uint32_t sp1 = kspan[e1], sp2 = kspan[e2], spe = kspan[e0];
			uint32_t maxsp = max(sp1, max(sp2, spe));
			if (oldmts < maxsp) return;
			uint32_t weight = max(maxsp, newMst);
			uint32_t tid1 = etoTID[e1], tid2 = etoTID[e2];
			if (tid1 != e0TID) {
				auto [it, inserted] = treeNode.try_emplace(tid1, weight);
				if (!inserted) if (weight < it->second) it->second = weight;
			}
			if (tid2 != e0TID) {
				auto [it, inserted] = treeNode.try_emplace(tid2, weight);
				if (!inserted) if (weight < it->second) it->second = weight;
			}
		});
		if (treeNode.empty()) continue;
		uint32_t e0kspan = idtn[e0TID].delta_;
		for (auto& [tid, weight] : treeNode) {
			uint32_t ancestorofe0 = findAncestorEqualOrLargeD(idtn, e0TID, weight);
			uint32_t delta = idtn[tid].delta_;
			if (delta < weight) {
				auto [ancestor, fatherofancestor] = findAncestorAndFatherLessD(idtn, tid, weight);
				if (fatherofancestor == e0TID || ancestor == e0TID || fatherofancestor == ancestorofe0) continue;
				if (weight == e0kspan) {
					if (fatherofancestor != UINT32_MAX) {
						delSuperEdge(idtn[fatherofancestor], idtn[ancestor]);
					} 
					addSuperEdge(idtn[e0TID], idtn[ancestor]);
					if (fatherofancestor != UINT32_MAX) treeNodeToInsert.insert(fatherofancestor);
				} else{
					if (treeNodeToInsert.count(ancestor)) continue;
					if (fatherofancestor == UINT32_MAX) {
						uint32_t newemptyNodeId = forest.maxBlockId++;
						auto [tn, inserted]	= idtn.emplace(newemptyNodeId, TN(weight, newemptyNodeId));
						tn->second.children_.push_back(ancestor);
						idtn[ancestor].father_ = newemptyNodeId;
						treeNodeToInsert.insert(newemptyNodeId);
					} else {
						if (idtn[fatherofancestor].delta_ != weight) {
							uint32_t newemptyNodeId = forest.maxBlockId++;
							auto [tn, inserted]	= idtn.emplace(newemptyNodeId, TN(weight, newemptyNodeId));
							delSuperEdge(idtn[fatherofancestor], idtn[ancestor]);
							tn->second.children_.push_back(ancestor);
							idtn[ancestor].father_ = newemptyNodeId;
							treeNodeToInsert.insert(newemptyNodeId);
							treeNodeToInsert.insert(fatherofancestor);
						} else {
							treeNodeToInsert.insert(fatherofancestor);
						}
					}
				}
			} else {
				if (ancestorofe0 != tid) treeNodeToInsert.insert(tid);
			}
		}
		Serialmerge(e0TID, k, treeNodeToInsert);
	}
}


void MbaF::KSECForestConstruction() {
	kmax_ = kdt_.kmax_ ; 
	ASSERT(m_ == tg_.m_); ASSERT(l_ == tg_.l_);
	initKSECForest(kmax_ + 1);
	initAUF(m_);
	vector<uint32_t> edgeProcessed(m_);
	for (uint32_t k = 1; k <= kmax_; k++) {
		constructForestForK(k, edgeProcessed);
		auf->reset();
	}
	return;
}

void MbaF::constructForestForK(uint32_t k, vector<uint32_t>& edgeProcessed) {
	auto& forest = ksecforest->forests[k];
	auto& idtn = forest.idTN;
	auto& kspan = kdt_.span_[k];
	tmax_ = tg_.maxtimestamps_;
	vector<vector<uint32_t>> kspan_dict(tmax_ + 1);
	getKSpanDict(k, kspan_dict);
	vector<unordered_map<uint32_t, uint32_t>> edgeMap(m_); 
	vector<vector<tuple<uint32_t, uint32_t, uint32_t>>> superE(tmax_ + 1);
	vector<KSEB> ksebSet;
    vector<uint32_t> q;
	uint32_t& treeNodeId = forest.maxBlockId;
	for (uint32_t delta = 0; delta <= tmax_; delta++) {
		if (kspan_dict[delta].size() == 0) continue;
		for (uint32_t e0 : kspan_dict[delta]) {
			if (edgeProcessed[e0] == k) continue;
			edgeProcessed[e0] = k;
			ksebSet.emplace_back(KSEB());
			auto [tn, inserted] = idtn.emplace(treeNodeId, TN(delta, treeNodeId));
			tn->second.elements_.push_back(treeNodeId);
			KSEB& curKseb = ksebSet.back();
			q.push_back(e0);
			while (!q.empty()){
				uint32_t e = q.back(); q.pop_back();
				curKseb.emplace_back(e);
				if (!edgeMap[e].empty()) {
					for (auto& [exploredTreeId, t] : edgeMap[e]) {
						if (exploredTreeId == treeNodeId) continue;
						superE[t].push_back({delta, exploredTreeId, treeNodeId});
					}
				}
				tg_.ForEachTriangle(e, k, [&](uint32_t e1, uint32_t e2) {
					if (edgeProcessed[e1] == k && edgeProcessed[e2] == k) return;
					uint32_t sp1 = kspan[e1], sp2 = kspan[e2];
					uint32_t triMst = tg_.GetMts(tg_.tau_[e], tg_.tau_[e1], tg_.tau_[e2]);
					uint32_t maxsp = max(sp1, max(sp2, delta));
					if (edgeProcessed[e1] != k && sp1 == maxsp) {
						if (sp1 == delta && triMst <= delta) {
							q.push_back(e1); edgeProcessed[e1] = k;
						} else {
							auto [it, inserted] = edgeMap[e1].try_emplace(treeNodeId, triMst);
							if (!inserted) {
								if (triMst < it->second) it->second = triMst;
							} 
						}
					}
					// process edge e2
					if (edgeProcessed[e2] != k && sp2 == maxsp) {
						if (sp2 == delta && triMst <= delta) {
							q.push_back(e2); edgeProcessed[e2] = k;
						} else {
							auto [it, inserted] = edgeMap[e2].try_emplace(treeNodeId, triMst);
							if (!inserted) {
								if (triMst < it->second) it->second = triMst;
							}  
						}
					}
				});
			}
			treeNodeId++;
		}
	}
	if (treeNodeId == 0) return;
	uint32_t deltaKMax = idtn[treeNodeId - 1].delta_;
	vector<vector<pair<uint32_t, uint32_t>>> sortedSuperEdges(deltaKMax + 1);
	vector<vector<pair<uint32_t,uint32_t>>> waitlist(tmax_ + 1);
	vector<EdgeHash> seenWait(tmax_ + 1); 
	EdgeHash seenSuperE;
	for (uint32_t t = 0; t <= tmax_; t++) {
		for (auto [sp, bidx, bidy] : superE[t]) {
			if (t > sp) {
				if (tryInsertEdge(seenSuperE, bidx, bidy)) {
					waitlist[t].push_back({bidy, bidx});
				}
			}
			else {
				if (tryInsertEdge(seenSuperE, bidx, bidy)) {
					sortedSuperEdges[sp].push_back({bidy, bidx});
				}
			}
			
		}		
	} 
	constructForestByGraph(forest, waitlist, sortedSuperEdges);
	uint32_t size = ksebSet.size();
	for(auto& [tid, tn] : idtn) {
		vector<uint32_t> BIds = std::move(tn.elements_);
		tn.elements_.clear(); 
		size_t edge_cnt = 0;
		for (uint32_t bid : BIds) {
			if (bid >= size) continue;
			q.push_back(bid); edge_cnt += ksebSet[bid].size();
		}
		tn.elements_.reserve(edge_cnt);
		while (!q.empty()) {
			uint32_t sid = q.back(); q.pop_back();
			tn.elements_.insert(tn.elements_.end(), std::make_move_iterator(ksebSet[sid].begin()), std::make_move_iterator(ksebSet[sid].end()));
		}
		for (uint32_t e : tn.elements_) { 
			forest.etoTreeNID[e] = tid; 
		}
	}
	return;
}

void MbaF::constructForestByGraph(Forest& forest, vector<vector<pair<uint32_t, uint32_t>>>& waitlist, vector<vector<pair<uint32_t, uint32_t>>>& superEdges) {
	auto& idtn = forest.idTN;
	uint32_t numOfBlocks = idtn.size();
	vector<uint32_t> BidtoCid(2 * numOfBlocks);
	iota(BidtoCid.begin(), BidtoCid.end(), 0);
	vector<pair<uint32_t, uint32_t>> Buffer;
	uint32_t maxKspan = superEdges.size() - 1, maxmts = waitlist.size() - 1;
	for (uint32_t delta = 0; delta <= maxKspan; delta++) {
		auto& superE = superEdges[delta];
		for (auto [x, y] : superE) {
			uint32_t ry = auf->find(y); uint32_t anchorid = auf->anchor_[ry];
			uint32_t lowLevelId = BidtoCid[anchorid]; uint32_t highLevelId = BidtoCid[x];
			TN& lowLevelTN = idtn[lowLevelId], &highLevelTN = idtn[highLevelId];
			if (lowLevelTN.father_ == highLevelId) continue;
			if (lowLevelTN.father_ == UINT32_MAX) {
				lowLevelTN.father_ = highLevelId;
				highLevelTN.children_.push_back(lowLevelId);
			} else {
				mergeNode(lowLevelTN.father_, highLevelId, idtn, BidtoCid);
			} 
			Buffer.push_back({ry, x}); 
		}

		for (auto [x, y] : Buffer) {
			auf->update(x, y, y);
		}
		Buffer.clear();
		for (auto [x, y] : waitlist[delta]) {
			uint32_t rx = auf->find(x), ry = auf->find(y);
			if (rx == ry) continue;
			dealwithWaitedSE(rx, ry, BidtoCid, forest, delta);
		}
	}	
	
	for (uint32_t delta = maxKspan + 1; delta <= maxmts; delta++) {
		if (waitlist[delta].size() == 0) continue;
		for (auto [x, y] : waitlist[delta]) {
			uint32_t rx = auf->find(x), ry = auf->find(y);
            if (rx == ry) continue;
            dealwithWaitedSE(rx, ry, BidtoCid, forest, delta);
		}
    }
	return;
}

inline void MbaF::dealwithWaitedSE(uint32_t rx, uint32_t ry, vector<uint32_t>& BidtoCid, Forest& forest, uint32_t delta) {
	auto& idtn = forest.idTN;
	uint32_t ax = auf->anchor_[rx], ay = auf->anchor_[ry];
	uint32_t tidx = BidtoCid[ax], tidy = BidtoCid[ay];
	uint32_t deltax = idtn[tidx].delta_, deltay = idtn[tidy].delta_;
	uint32_t maxdelta = max(deltax, deltay);
	if (maxdelta < delta) {
		uint32_t& treeNodeId = forest.maxBlockId;
		auto [tn, inserted] = idtn.emplace(treeNodeId, TN(delta, treeNodeId));
		tn->second.elements_.push_back(treeNodeId);
		idtn[treeNodeId].children_.push_back(tidx);
		idtn[tidx].father_ = treeNodeId;
		idtn[treeNodeId].children_.push_back(tidy);
		idtn[tidy].father_ = treeNodeId;
		auf->update(rx, treeNodeId, treeNodeId);
		auf->update(ry, treeNodeId, treeNodeId);
		treeNodeId++;
	} else {
		if (maxdelta == deltax && maxdelta == deltay) {
			mergeNode(tidx, tidy, idtn, BidtoCid);
			auf->update(rx, ry, ax);
		} else if (maxdelta == deltax) {
			idtn[tidy].father_ = tidx;
			idtn[tidx].children_.push_back(tidy);
			auf->update(rx, ry, ax);
		} else if (maxdelta == deltay) {
			idtn[tidx].father_ = tidy;
			idtn[tidy].children_.push_back(tidx);
			auf->update(rx, ry, ay);
		}
	}
	return ;
}

vector<vector<int>> MbaF::searchCommunityByKSpanIndex(int query, int k, uint32_t delta) {
	auto& kspan = kdt_.span_[k];
    vector<bool> Ins(m_);
    queue<int> q;
    vector<vector<int>> result;
    for (const auto& ae : tg_.adj_[query]) {
        int e = ae.eid;
		if (kspan[e] <= delta && !Ins[e]) {
			q.push(e); Ins[e] = true; 
			vector<int> Ai;
			while (!q.empty()) {
            	int eid = q.front(); q.pop(); Ai.push_back(eid);
				const auto tris =  tg_.GetTriangles(eid);
				for (const auto& tri : tris) {
					uint32_t e1 = tri.first, e2 = tri.second;
					if (kspan[e1] == UINT32_MAX || kspan[e2] == UINT32_MAX) continue;
					if (kspan[e1] > delta || kspan[e2] > delta) continue;
					uint32_t deltatri = tg_.GetMts(tg_.tau_[eid], tg_.tau_[e1], tg_.tau_[e2]);
					if (deltatri > delta) continue;
					if (!Ins[e1]) {
						q.push(e1); Ins[e1] = true;
					}
					if (!Ins[e2]) {
						q.push(e2); Ins[e2] = true;
					}
				}
        	}
			result.emplace_back(std::move(Ai));
		}
    }
	int size = 0;
    for (auto& community : result) {
		size += community.size();
		// for (int& e : community) {
		// 	cout<<e<<" ";
		// }
	}
	cout<<size<<endl;
    return result;
 }




vector<vector<uint32_t>> MbaF::searchCommunityByKSECForest(uint32_t query, uint32_t k, uint32_t delta) {
	Forest& forest = ksecforest->forests[k];
	auto& idtn = forest.idTN;
    unordered_set<uint32_t> Ins;
    queue<uint32_t> q;
    vector<vector<uint32_t>> result;
    for (const auto& ae : tg_.adj_[query]) {
        uint32_t e = ae.eid;
        uint32_t treeid = forest.etoTreeNID[e];
        if (treeid == UINT32_MAX || Ins.count(treeid)) continue;
        const TN* node = &forest.idTN[treeid];
        if (node->delta_ > delta) continue;
        while (node->father_ != UINT32_MAX && idtn[node->father_].delta_ <= delta) {
            node = &idtn[node->father_];
        }
        vector<uint32_t> Ai;
        q.push(node->id_);
        while (!q.empty()) {
            uint32_t treeid = q.front(); q.pop();
            Ins.insert(treeid);
            const TN* node = &idtn[treeid];
            Ai.insert(Ai.end(), node->elements_.begin(), node->elements_.end());
            for (uint32_t cid : node->children_) {
                q.push(cid);
            }
        }
        result.emplace_back(std::move(Ai));
    }
	int size = 0;
    for (auto& community : result) {
		size += community.size();
		// for (int& e : community) {
		// 	cout<<e<<" ";
		// }
	}
	cout<<size<<endl;
    return result;
}


void MbaF::writeKSECForestIndexToFile(const std::string& filename) {
    std::ofstream outfile(filename, std::ios::binary);
    if (!outfile) {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
        return;
    }
	double fsz = 0.0;
    /* ---------- write KDE forests ---------- */
    outfile.write((char*)&kmax_, sizeof(kmax_));
	outfile.write((char*)&m_, sizeof(m_));
	fsz += sizeof(kmax_); fsz += sizeof(m_);
    for (uint32_t k = 1; k <= kmax_; ++k) {
        auto& forest = ksecforest->forests[k];
		/* ---- maxBlockid ---- */
        outfile.write((char*)&forest.maxBlockId, sizeof(forest.maxBlockId));
		fsz += sizeof(uint32_t); // maxBlockId
        uint32_t tn_cnt = forest.idTN.size();
        outfile.write((char*)&tn_cnt, sizeof(tn_cnt));
		fsz += sizeof(uint32_t); // tn_cnt
        for (auto& [tid, tn] : forest.idTN) {
            outfile.write((char*)&tid, sizeof(tid));
            uint32_t delta = tn.delta_;
            outfile.write((char*)&delta, sizeof(delta));
            uint32_t father = tn.father_;
            outfile.write((char*)&father, sizeof(father));
            uint32_t ch_cnt = tn.children_.size();
            outfile.write((char*)&ch_cnt, sizeof(ch_cnt));
            outfile.write((char*)tn.children_.data(), ch_cnt * sizeof(uint32_t));
            uint32_t edge_cnt = tn.elements_.size();
            outfile.write((char*)&edge_cnt, sizeof(edge_cnt));
            outfile.write((char*)tn.elements_.data(), edge_cnt * sizeof(uint32_t));
			fsz += (5 * sizeof(uint32_t));                 // tid, delta, father, ch_cnt, edge_cnt
            fsz += tn.children_.size() * sizeof(uint32_t);
            fsz += tn.elements_.size() * sizeof(uint32_t);
        }
    }
	/* ---- etoTreeNID ---- */
	
	for (uint32_t e = 0; e < m_; e++) {
		uint32_t maxk_to_e = k_[e];
		outfile.write((char*)&maxk_to_e, sizeof(maxk_to_e));
		fsz += sizeof(uint32_t);                     // maxk
		for (uint32_t k = 1; k <= maxk_to_e; k++) {
			uint32_t tid = ksecforest->forests[k].etoTreeNID[e];
			outfile.write((char*)&tid, sizeof(tid));
			fsz += sizeof(tid); // tid
		}
	}
		
    outfile.close();
    std::cout << "Index successfully written to: " << filename << std::endl;
	double mb = fsz / (1024.0 * 1024.0);
	printf("Index size: %.2f MB\n", mb);

}

void MbaF::readKSECForestIndexFromFile(const std::string& filename) {
    std::ifstream infile(filename, std::ios::binary);
    if (!infile) {
        std::cerr << "Failed to open file for reading: " << filename << std::endl;
        return;
    }

    /* -------- read forests -------- */
    infile.read((char*)&kmax_, sizeof(kmax_));
	infile.read((char*)&m_, sizeof(m_));
    initKSECForest(kmax_ + 1);

    for (uint32_t k = 1; k <= kmax_; ++k) {
        auto& forest = ksecforest->forests[k];
		/* ---- fid ---- */
        infile.read((char*)&forest.maxBlockId, sizeof(forest.maxBlockId));

        /* -------- idTN -------- */
        uint32_t tn_cnt;
        infile.read((char*)&tn_cnt, sizeof(tn_cnt));
		forest.idTN.reserve(tn_cnt);
        for (uint32_t i = 0; i < tn_cnt; i++) {
            uint32_t tid;
            infile.read((char*)&tid, sizeof(tid));

            TN tn;  
			tn.id_ = tid;
            infile.read((char*)&tn.delta_, sizeof(uint32_t));
            infile.read((char*)&tn.father_, sizeof(uint32_t));

            uint32_t ch_cnt;
            infile.read((char*)&ch_cnt, sizeof(ch_cnt));
            tn.children_.resize(ch_cnt);
            infile.read((char*)tn.children_.data(),
                        ch_cnt * sizeof(uint32_t));

            uint32_t edge_cnt;
            infile.read((char*)&edge_cnt, sizeof(edge_cnt));
            tn.elements_.resize(edge_cnt);
            infile.read((char*)tn.elements_.data(),
                        edge_cnt * sizeof(uint32_t));

			forest.idTN.emplace(tid, std::move(tn));
        }
    }

	/* -------- etoTreeNID -------- */
	for (uint32_t e = 0; e < tg_.m_; e++) {
        infile.read((char*)(&k_[e]), sizeof(k_[e]));
        for (uint32_t k = 1; k <= k_[e]; k++) {
            uint32_t tid;
            infile.read((char*)(&tid), sizeof(tid));
            ksecforest->forests[k].etoTreeNID[e] = tid;
        }
    }

    infile.close();
    // std::cout << "Index successfully read from: " << filename << std::endl;
}

void MbaF::writeGraphToFile(const std::string& filename) {
    std::ofstream outfile(filename, std::ios::binary);
    if (!outfile) {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
        return;
    }
	double fsz = 0.0;
	outfile.write(reinterpret_cast<const char*>(&tg_.maxtimestamps_), sizeof(tg_.maxtimestamps_));
	outfile.write(reinterpret_cast<const char*>(&tg_.n_), sizeof(tg_.n_));
    outfile.write(reinterpret_cast<const char*>(&tg_.m_), sizeof(tg_.m_));
	fsz += sizeof(tg_.maxtimestamps_); fsz += sizeof(tg_.n_); fsz += sizeof(tg_.m_);
    for (uint32_t i = 0; i < tg_.m_; i++) {
        uint32_t tsize = tg_.tau_[i].size();
        uint32_t header[3] = {
            tg_.edge_info_[i].first,
            tg_.edge_info_[i].second,
            tsize
        };
        outfile.write(reinterpret_cast<const char*>(header), sizeof(header));
        if (tsize > 0) {
            outfile.write(
                reinterpret_cast<const char*>(tg_.tau_[i].data()),
                tsize * sizeof(uint32_t)
            );
        }
		fsz += sizeof(header);
		fsz += tsize * sizeof(uint32_t);
    }
    outfile.close();
    std::cout << "Index successfully written to: " << filename << std::endl;
    double mb = fsz / (1024.0 * 1024.0);
    printf("Index size: %.2f MB\n", mb);
}

void MbaF::loadGraphIndex(const std::string& fn) {
    std::ifstream infile(fn, std::ios::binary);
    if (!infile) {
        std::cerr << "Failed to open index file: " << fn << std::endl;
        return;
    }
	uint32_t maxtimestamps, n, m;
    infile.read(reinterpret_cast<char*>(&maxtimestamps), sizeof(maxtimestamps))
          .read(reinterpret_cast<char*>(&n), sizeof(n))
          .read(reinterpret_cast<char*>(&m), sizeof(m));
	tg_.maxtimestamps_ = maxtimestamps;
    ASSERT(m <= l_ && n == tg_.n_);
    for (uint32_t e = 0, buf[3]; e < m; ++e) {
        infile.read(reinterpret_cast<char*>(buf), sizeof(buf));
        ASSERT(tg_.LazyInsert(buf[0], buf[1]) == e);
        uint32_t taue_size = buf[2];
        tg_.tau_[e].resize(taue_size);
		infile.read(reinterpret_cast<char*>(tg_.tau_[e].data()),
                taue_size * sizeof(uint32_t));
    }
    tg_.Rectify();
    ASSERT(tg_.m() == m);
    infile.close();
}


