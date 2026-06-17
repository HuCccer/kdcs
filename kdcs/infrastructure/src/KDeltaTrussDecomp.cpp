#include "../include/KDeltaTrussDecomp.h"
#include <fstream>
#include <cstring>
#include <unordered_set>

using namespace std;
using namespace KDeltaTrussCommunity;


KDeltaTruss::KDeltaTruss(const uint32_t n, const uint32_t l, TGraph& tg, std::vector<uint32_t>& k, std::unordered_map<EdgT, uint32_t, pairHash>& triMTSpan)
    : l_(l), n_(n), tg_(tg), k_(k), sup_(l_, 0), triMTSpan_(triMTSpan) {
    ASSERT_MSG(1 <= l_ && l_ < (static_cast<uint32_t>(1) << 29),
                "it is required 64 <= l <= 2^29 for the ease of implementation");
    visited    = std::vector<uint32_t>(l_, 0);
    removed   = std::vector<uint32_t>(l_, 0);
    ins   = std::vector<uint32_t>(l_, 0);
}

KDeltaTruss::KDeltaTruss(const uint32_t n, const uint32_t l, TGraph& tg, std::vector<uint32_t>& k, std::unordered_map<EdgT, uint32_t, pairHash>& triMTSpan, const std::string& fn)
    : KDeltaTruss(n, l, tg, k, triMTSpan) {
    loadKSpanIndex(fn);
}

void KDeltaTruss::loadKSpanIndex(const std::string& fn) {
    std::ifstream infile(fn, std::ios::binary);
    if (!infile) {
        std::cerr << "Failed to open index file: " << fn << std::endl;
        return;
    }
    infile.read(reinterpret_cast<char*>(&kmax_), sizeof(kmax_));
    span_.resize(kmax_ + 1);
    for (uint32_t k = 1; k <= kmax_; k++) {
        span_[k].assign(tg_.l_, UINT32_MAX);
    }
	infile.read(reinterpret_cast<char*>(&m_), sizeof(m_));
    for (uint32_t e = 0; e < tg_.m_; e++) {
        infile.read(reinterpret_cast<char*>(&k_[e]), sizeof(k_[e]));
        for (int k = 1; k <= k_[e]; k++) {
            uint32_t ksp = 0;
            infile.read(reinterpret_cast<char*>(&ksp), sizeof(ksp));
            span_[k][e] = ksp;
        }
    }
    infile.close();
}


void KDeltaTruss::trussDecomp() {
 // truss decomposition
    // 1. compute the support of each edge by triangle listing
    // 1.1. define a total order over the vertices
	m_ = tg_.m_;
    const auto pred = [this](const uint32_t v1, const uint32_t v2) {
    const size_t deg1 = tg_.adj_[v1].size();
    const size_t deg2 = tg_.adj_[v2].size();
        if (deg1 != deg2) return deg1 > deg2;
        else return v1 > v2;
    };
    vector<uint32_t> verts(n_);
  	iota(verts.begin(), verts.end(), 0);
  	sort(verts.begin(), verts.end(), pred);
  	// 1.2. call the "forward" algorithm to list triangles
    vector<uint32_t> trn_(m_, 0);
  	vector<vector<ArrayEntry>> A(n_);
    //delta-triangle list
	time2triangle.resize(tg_.maxtimestamps_ + 1);
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
					uint32_t interval = tg_.GetMts(tg_.tau_[e],tg_.tau_[A[v][pv].eid],tg_.tau_[A[u][pu].eid]);
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
	vector<uint32_t> bin(max_sup + 1, 0);
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
        tg_.ForEachTriangle(eid, [&](uint32_t e1, uint32_t e2) {
            if (trn_[e1] >= k && trn_[e2] >= k) ++ks_[eid];
            if (removed[e1] || removed[e2]) return;
            for (const uint32_t e : {e1, e2}) {
                if (trn_[e] > k) {  
                    const uint32_t pe3 = bin[trn_[e]], pe = pos[e];
                    if (pe3 != pe) {
                        const uint32_t e3 = ord_[pe3];
                        ord_[pe] = e3; pos[e3] = pe;
                        ord_[pe3] = e; pos[e] = pe3;
                    }
                    ++bin[trn_[e]];
                    --trn_[e];
                }
            }
        });
  	}
	kmax_ = k;
    for (uint32_t i = 0; i < m_; ++i) k_[i] = trn_[i];
}

void KDeltaTruss::kDeltaTrussDecomp() {
    trussDecomp();
	vector<uint32_t> trussness(m_);
	for (uint32_t i = 0; i < m_; ++i) trussness[i] = k_[i];
    vector<uint32_t> q;
	vector<bool> Ins(m_,false);
    span_.resize(kmax_ + 1, vector<uint32_t>(l_, UINT32_MAX));
    unordered_set<EdgT, pairHash> Mc;
	for(uint32_t t = tmax_ ; t > 0 ; --t){
		if(time2triangle[t].size() == 0) continue;
		for(const auto& [e1, e2, e3] : time2triangle[t]){
			Mc.insert({e1, e3}); 
            uint32_t k1 = trussness[e1], k2 = trussness[e2], k3 = trussness[e3];
			uint32_t mink = min(k1, min(k2, k3));
			//update k-support
			if(k1 == mink) { if(ks_[e1]-- == k1) { q.push_back(e1); Ins[e1] = true;} } 
			if(k2 == mink) { if(ks_[e2]-- == k2) { q.push_back(e2); Ins[e2] = true;} } 
			if(k3 == mink) { if(ks_[e3]-- == k3) { q.push_back(e3); Ins[e3] = true;} }
			//update the trussness of edge
			while(!q.empty()){
				uint32_t e = q.back(); q.pop_back();
                Ins[e] = false;
				uint32_t k = trussness[e];
				trussness[e]--; span_[k][e] = t;
				if(trussness[e] == 0) continue;
				uint32_t kse = 0;
				tg_.ForEachTriangle(e, [&](uint32_t e1, uint32_t e2) {
                    if (trussness[e1] < k - 1 || trussness [e2] < k - 1) return;
                    getTriId(e, e1, e2);
                    if (!Mc.count({mineid, maxeid})) {
                        kse++;
                        if (trussness[e1] >= k && trussness[e2] >= k) {
							if(trussness[e1] == k && Ins[e1] == false && ks_[e1]-- == trussness[e1]) {
								q.push_back(e1); Ins[e1] = true;
							}
							if(trussness[e2] == k && Ins[e2] == false && ks_[e2]-- == trussness[e2]) {
								q.push_back(e2); Ins[e2] = true;
							}	
						}
                    }

                });
				ks_[e] = kse;
			}
		}//for time2triangle[t]
	}
    for(uint32_t e = 0; e < m_; ++e){
		while(trussness[e] > 0){
			span_[trussness[e]][e] = 0; trussness[e]--; 
		}  
    }

	return;
}


vector<std::unordered_set<uint32_t>> KDeltaTruss::addEdgeWithVariety(uint32_t e, uint32_t timestamp, bool variety) {
	if (variety == true) insertE(e, timestamp);
	else insertT(e, timestamp);
	return affected;
}

std::vector<std::unordered_set<uint32_t>> KDeltaTruss::addEdge(uint32_t u, uint32_t v, uint32_t timestamp) {
	bool variety = false;
	uint32_t e = tg_.InsertEdge(u, v, variety);
	if (variety == true) insertE(e, timestamp);
	else insertT(e, timestamp);
	return affected;
}

vector<vector<uint32_t>>  KDeltaTruss::XH_truss_maintenance(const uint32_t eid) {
	vector<uint32_t> level(kmax_ + 1);
	uint32_t max_e_tri_level = 0;
	tg_.ForEachTriangle(eid, [&](uint32_t e1, uint32_t e2) {
		uint32_t e_tri_level = std::min(k_[e1],k_[e2]);
		max_e_tri_level = std::max(max_e_tri_level, e_tri_level);
		level[e_tri_level]++;
	});
	uint32_t k1 = 0, k2 = 0, count_k = 0;
    for (int k = max_e_tri_level; k >= 0; --k) {
		count_k += level[k];
        if (k1 == 0 && count_k >= k) k1 = k; 	
        if (k2 == 0 && count_k >= k + 1) k2 = k + 1;
        if(k1 && k2) break;
    }
	// the maximum trussness of e is 2;
	if (k2 == 0) return{};
	k_[eid] = k1;
	uint32_t m = tg_.m();
	// the maximum trussness of the edge that trussness can increase is k2 - 1
	vector<vector<uint32_t>> L(k2);
	vector<uint32_t> status(m, UINT32_MAX);
	vector<vector<uint32_t>> inc(k2 + 1);
	if (k1 == k2 - 1) {
		L[k1].push_back(eid); ins[eid] = token;
	} 
	for (uint32_t k = k1; k > 0; k--) {
		inc[k].push_back(eid);
	}
	tg_.ForEachTriangle(eid, [&](uint32_t e1, uint32_t e2) {
		uint32_t min_k = std::min(k_[e1], k_[e2]);
		if (min_k <= k2 - 1) {
			if (k_[e1] == min_k) {
				L[min_k].push_back(e1); ins[e1] = token;
			} 
			if (k_[e2] == min_k) {
				L[min_k].push_back(e2); ins[e2] = token;
			} 
		}
	});
	vector<uint32_t> subq;
	//truss maintenance
	vector<uint32_t> tmpEdge;
	for (int k = k2 - 1; k >= 0; k--) {
		uint32_t Lidx = 0;
		while (Lidx < L[k].size()) {
			uint32_t leid = L[k][Lidx++];
			tmpEdge.clear();
			sup_[leid] = 0;
			tg_.ForEachTriangle(leid, [&](uint32_t e1, uint32_t e2) {
				if (k_[e1] < k || k_[e2] < k) return;
				if (k_[e1] == k && status[e1] == k) return;
				if (k_[e2] == k && status[e2] == k) return;
				sup_[leid]++;
				if (k_[e1] == k && ins[e1] != token) tmpEdge.push_back(e1);
				if (k_[e2] == k && ins[e2] != token) tmpEdge.push_back(e2);
			});
			if (sup_[leid] > k) {
				status[leid] = k + 1;
				for (uint32_t teid : tmpEdge) {
					L[k].push_back(teid); ins[teid] = token;
				}
			} else {
				subq.clear();
				uint32_t subqidx = 0;
				subq.push_back(leid);
				while (subqidx < subq.size()){
					uint32_t subqeid = subq[subqidx++]; status[subqeid] = k;
					tg_.ForEachTriangle(subqeid, [&](uint32_t e1, uint32_t e2) {
						if (k_[e1] < k || k_[e2] < k) return;
						if (k_[e1] == k && status[e1] == k) return;
						if (k_[e2] == k && status[e2] == k) return;
						if (k_[e1] == k && status[e1] == k + 1) {
							if (--sup_[e1] == k) subq.push_back(e1);
						} 
						if (k_[e2] == k && status[e2] == k + 1) {
							if (--sup_[e2] == k) subq.push_back(e2);
						}
					});
				}
			}
		}
		for (uint32_t leid : L[k]) {
			if (status[leid] == k + 1) {
				k_[leid]++; inc[k + 1].push_back(leid);
			} 
			sup_[leid] = 0;
		}
	}
	token++;
	return inc;
}

void KDeltaTruss::assignKspanForEdges(vector<vector<uint32_t>>& inc, uint32_t k) {
	auto& kspan = span_[k];
	uint32_t level = 0;
	for (uint32_t eid : inc[k]) {
		tg_.ForEachTriangle(eid, [&](uint32_t e1, uint32_t e2) {
			if (k_[e1] < k || k_[e2] < k) return;
			uint32_t mts = tg_.GetMts(tg_.tau_[eid], tg_.tau_[e1], tg_.tau_[e2]);
			level = max(level, mts);
			if (kspan[e1] < UINT32_MAX) level = max(level, kspan[e1]);
			if (kspan[e2] < UINT32_MAX) level = max(level, kspan[e2]);
		});
	}
	for (uint32_t eid : inc[k]) {
		kspan[eid] = level;
		affected[k].insert(eid);
	}
	return;
}



void KDeltaTruss::insertE(const uint32_t e, const uint32_t timestamp) {
	tg_.InsertTimestamp(e, timestamp);
	const auto& tris = tg_.GetTriangles(e);
	vector<vector<uint32_t>> inc = XH_truss_maintenance(e);
	uint32_t ke = k_[e];
	affected.resize(ke + 1);
	for (uint32_t k = 1; k <= ke; k++) affected[k].clear();
	if (ke == 0) return;
	if (ke > kmax_) {
		kmax_ = ke; span_.resize(kmax_ + 1);
		span_[kmax_].resize(l_, UINT32_MAX);
	}
	uint32_t m = tg_.m_;
	for (uint32_t k = ke; k > 1; k--) {
	 	assignKspanForEdges(inc, k);
		maintenance(k, e);
		// clear maintenance buffers
        for (uint32_t mts : touchedDeltaLevels_) deltaTriListBuf_[mts].clear();
		for (uint32_t ksp : touchedCandidateLevels_) candidateEdgeBuf_[ksp].clear();
        touchedCandidateLevels_.clear();
        touchedDeltaLevels_.clear();
		// update label token
		token++;
		visittoken += 2;
	}
	// especially maintenance for k = 1
	auto& kspan = span_[1];
	uint32_t min_mts = UINT32_MAX;
	for (const auto& [e1, e2] : tris) {
		uint32_t kspane1 = kspan[e1], kspane2 = kspan[e2];
		uint32_t new_mts = tg_.GetMts(tg_.tau_[e], tg_.tau_[e1], tg_.tau_[e2]);
		if (new_mts < kspane1) {
			kspan[e1] = new_mts; affected[1].insert(e1);
		} 
		if (new_mts < kspane2) {
			kspan[e2] = new_mts; affected[1].insert(e2);
		} 
		min_mts = min(min_mts, new_mts);
	}
	if (min_mts < kspan[e]) {
		kspan[e] = min_mts; affected[1].insert(e);
	} 
}
	


void KDeltaTruss::insertT(const uint32_t e, const uint32_t timestamp) {
	vector<uint32_t> oldtaue = tg_.tau_[e];
	const auto& trise = tg_.GetTriangles(e);
	tg_.InsertTimestamp(e, timestamp);
	uint32_t ke = k_[e], m = tg_.m_;
	affected.resize(ke + 1);
	for (uint32_t k = 1; k <= ke; k++) affected[k].clear();
	if (ke == 0) return;
	for (uint32_t k = ke; k > 1; k--) {
		auto& kspan = span_[k];
		for (const auto& [e1, e2] : trise) {
			if (kspan[e1] == UINT32_MAX || kspan[e2] == UINT32_MAX) continue;
			uint32_t ori_tri_mts = tg_.GetMts(oldtaue, tg_.tau_[e1], tg_.tau_[e2]);
			uint32_t cur_tri_mts = tg_.GetMts(tg_.tau_[e], tg_.tau_[e1], tg_.tau_[e2]);
			uint32_t maxKSpan = max(kspan[e], max(kspan[e1], kspan[e2]));
			if (cur_tri_mts >= maxKSpan || ori_tri_mts < maxKSpan || ori_tri_mts == cur_tri_mts) continue;
			uint32_t seed;
			if (kspan[e1] == maxKSpan) seed = e1;
			else if (kspan[e2] == maxKSpan) seed = e2;
			else seed = e;
			maintenance(k, seed);
			if (kspan[seed] < kspan[e]) kspan[e] = kspan[seed];
			clearMaintenanceBuffers(kspan[seed]);
			token++; visittoken += 2;
		}
	}
	// especially maintenance for k = 1
	auto& kspan = span_[1];
	uint32_t min_mts = UINT32_MAX;
	for (const auto& [e1, e2] : trise) {
		uint32_t kspane1 = kspan[e1], kspane2 = kspan[e2];
		uint32_t new_mts = tg_.GetMts(tg_.tau_[e], tg_.tau_[e1], tg_.tau_[e2]);
		if (new_mts < kspane1) {
			kspan[e1] = new_mts; affected[1].insert(e1);
		} 
		if (new_mts < kspane2) {
			kspan[e2] = new_mts; affected[1].insert(e2);
		}
		min_mts = min(min_mts, new_mts);
	}
	if (min_mts < kspan[e]) {
		kspan[e] = min_mts; affected[1].insert(e);
	} 
}

void KDeltaTruss::maintenance(uint32_t k, uint32_t seed) {
	auto& kspan = span_[k];
	uint32_t level = kspan[seed];
	TriIdMap validTri;
	if (candidateEdgeBuf_.size() <= level) {
    	candidateEdgeBuf_.resize(level + 1);
	}
	if (deltaTriListBuf_.size() <= level) {
		deltaTriListBuf_.resize(level + 1);
	}
	vector<uint32_t> seeds, q;
	seeds.push_back(seed); visited[seed] = visittoken + 1;
	while (true) {
		if (candidateEdgeBuf_[level].size() || seeds.size()) {
			expand(validTri, seeds, visited, level, k, visittoken);
		}
		for (auto& [e1, e2, e3] : deltaTriListBuf_[level]) {
			if (removed[e1] == token || removed[e2] == token || removed[e3] == token) continue;
			validTri[{e1, e3}] = false;
			sup_[e1]--; sup_[e2]--; sup_[e3]--;
			if(sup_[e1] < k && kspan[e1] >= level) { q.push_back(e1); ins[e1] = token;}
			if(sup_[e2] < k && kspan[e2] >= level) { q.push_back(e2); ins[e2] = token;}
			if(sup_[e3] < k && kspan[e3] >= level) { q.push_back(e3); ins[e3] = token;}
			while(!q.empty()){
				uint32_t eid = q.back(); q.pop_back();
				removed[eid] = token;
				if (kspan[eid] != level) {
					affected[k].insert(eid); kspan[eid] = level;
				} 
				tg_.ForEachTriangle(eid, [&](uint32_t e1, uint32_t e2) {
					if(kspan[e1] == UINT32_MAX || kspan[e2] == UINT32_MAX) return;
					if (removed[e1] == token || removed[e2] == token) return;
					auto tri = sortTri(eid, e1, e2);
					uint32_t mintid = get<0>(tri), maxtid = get<2>(tri);
					uint32_t mts = tg_.GetMts(tg_.tau_[eid], tg_.tau_[e1], tg_.tau_[e2]);
					if (mts > level) return;
					auto it = validTri.find({mintid, maxtid});
					if(it == validTri.end() || it->second == false) return;
					it->second = false;
					sup_[eid]--; sup_[e1]--; sup_[e2]--;
					if(sup_[e1] < k && kspan[e1] >= level && ins[e1] != token){
						q.push_back(e1); ins[e1] = token;
					}
					if(sup_[e2] < k && kspan[e2] >= level && ins[e2] != token){
						q.push_back(e2); ins[e2] = token;
					} 	
				});
			}
		}	
		if (removed[seed] == token) break;
		if (level == 0) break;
		level--;
	}
}

void KDeltaTruss::expand(TriIdMap& validTri, vector<uint32_t>& seeds, vector<uint32_t>& visited, uint32_t level, uint32_t k, uint32_t visittoken) {
	if (level == 0) return;
	const uint32_t seenMark = visittoken + 1;
	const uint32_t doneMark = visittoken + 2;
	for (uint32_t eid : candidateEdgeBuf_[level]) {
		if (sup_[eid] == 0) visited[eid] = visittoken;
		else seeds.push_back(eid);
	}
	auto& kspan = span_[k];
	while (!seeds.empty()) {
		uint32_t eid = seeds.back(); seeds.pop_back();
		visited[eid] = doneMark;
		tg_.ForEachTriangle(eid, [&](uint32_t e1, uint32_t e2) {
			if (kspan[e1] > level || kspan[e2] > level) return;
			auto tri = sortTri(eid, e1, e2);
			uint32_t mintid = get<0>(tri), maxtid = get<2>(tri);
			uint32_t mts = tg_.GetMts(tg_.tau_[eid], tg_.tau_[e1], tg_.tau_[e2]);
			if (mts > level) return;
			if (visited[e1] != doneMark && visited[e2] != doneMark) {
				addDeltaTriangle(mts, tri);
				validTri[{mintid, maxtid}] = true;
				sup_[eid]++; sup_[e1]++; sup_[e2]++;
			}
			if (visited[e1] != seenMark && visited[e1] != doneMark) {
				if (kspan[e1] == level) {
					seeds.push_back(e1);  visited[e1] = seenMark;
				} else {
					if (mts < kspan[e1]) {
						addCandidate(kspan[e1], e1); visited[e1] = seenMark;
					}
				} 
			}
			if (visited[e2] != seenMark && visited[e2] != doneMark) {
				if (kspan[e2] == level) {
					seeds.push_back(e2);  visited[e2] = seenMark;
				} else {
					if (mts < kspan[e2]) {
						addCandidate(kspan[e2], e2); visited[e2] = seenMark;
					}
				}
			}
		});
	}
}



void KDeltaTruss::writeKspanIndexToFile(const std::string& filename) {
    std::ofstream outfile(filename, std::ios::binary);
    if (!outfile) {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
        return;
    }
    double fsz = 0.0;
    outfile.write(reinterpret_cast<const char*>(&kmax_), sizeof(kmax_));
	outfile.write(reinterpret_cast<const char*>(&tg_.m_), sizeof(tg_.m_));
    fsz += sizeof(kmax_);
    fsz += sizeof(m_);
    for (uint32_t e = 0; e < m_; e++) {
        int maxk_to_e = k_[e];
        outfile.write(reinterpret_cast<const char*>(&maxk_to_e), sizeof(maxk_to_e));
        fsz += sizeof(maxk_to_e);
        for (int k = 1; k <= maxk_to_e; k++) {
            uint32_t ksp = span_[k][e];
            outfile.write(reinterpret_cast<const char*>(&ksp), sizeof(ksp));
            fsz += sizeof(ksp);
        }
    }
    outfile.close();
    std::cout << "Index successfully written to: " << filename << std::endl;
    double mb = fsz / (1024.0 * 1024.0);
    printf("Index size: %.2f MB\n", mb);
}







        
