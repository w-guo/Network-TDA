/*
 * utility.h
 *
 *  Created on: May 9, 2019
 *      Author: robustor
 */

#ifndef UTILITY_H_
#define UTILITY_H_

#include "Trie.h"

typedef pair<int, int> edge;
typedef vector<edge> EdgeSet;
typedef map<int, list<int>> mapNBVer;
typedef map<int, Trie *> mapTrie;
typedef list<MCVec> MCSet; // each MC is sorted

MCPair * getNeighborMCs(MCVec * MC, int m, mapNBVer * mLV, mapTrie * mTrie) {
	// create a set to store neighboring MCs of the given MC
	MCPair * nbSet = new MCPair();

	for (auto vid : *MC) {
		if (mTrie->find(vid) != mTrie->end()) {
			Trie * T = (*mTrie)[vid];
			MCPair * MCP = T->outputMC(T->root);
			nbSet->insert(MCP->begin(), MCP->end());
			MCP->clear();
		}
	}
	nbSet->erase(m);

	for (auto vid : *MC) {
		if (mLV->find(vid) != mLV->end()) {
			list<int> Lv = (*mLV)[vid];
			for (auto vit = MC->begin(); *vit<vid; vit++) {
				Lv.remove(*vit);
			}
			for (auto uid : Lv) {
				if (mTrie->find(uid) != mTrie->end()) {
					Trie * T = (*mTrie)[uid];
					if (T->NodeList.find(vid) != T->NodeList.end()) {
						for (auto lit : (T->NodeList)[vid]) {
							MCPair * MCP = T->outputMC(lit);
							nbSet->insert(MCP->begin(), MCP->end());
							MCP->clear();
						}
					}
				}

			}
		}
	}

	return nbSet;
}

MCSet * listSubsumedMCsEdgeDeletion(EdgeSet * delEdges, mapNBVer * mLV, mapNBVer * mHV, MCSet * delMCs) {
	MCSet * newMcs = new MCSet();
	MCSet * tempMCs = new MCSet();;
	for (auto  MCvec : *delMCs) {
		tempMCs->push_back(MCvec);
	}
}

#endif /* UTILITY_H_ */
