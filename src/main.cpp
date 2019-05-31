//============================================================================
// Name        : TDA_networks.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description :
//============================================================================

#include <iostream>
#include "CliqueGraph.h"
using namespace std;



//PD * buildCT(mapNBVer * mLV, mapNBVer * mHV);
int main() {

	MCVec mc1 = {1, 2, 3, 6};
	MCVec mc2 = {2, 3, 4, 5, 6};
	MCVec mc3 = {4, 5, 7, 12};
	MCVec mc4 = {5, 7, 8, 9};
	MCVec mc5 = {6, 10, 11};
	MCVec mc6 = {1, 2, 6};
	MCVec mc7 = {4, 5, 7};
	MCVec mc8 = {6, 10};
	MCVec mc9 = {11};
	MCVec mc10 = {12};
	MCVec mc11 = {5, 6, 9};

//	Trie* T1 = new Trie(1);
//	T1->addMC(&mc1, 10);
//	Trie* T2 = new Trie(2);
//	T2->addMC(&mc2, 11);
//	Trie* T3 = new Trie(4);
//	T3->addMC(&mc3, 12);
//	Trie* T4 = new Trie(5);
//	T4->addMC(&mc4, 13);
//	Trie* T5 = new Trie(6);
//	T5->addMC(&mc5, 14);

	mapTrie * mTrie = new mapTrie();
//	mTrie->insert(make_pair(T1->root->key, T1));
//	mTrie->insert(make_pair(T2->root->key, T2));
//	mTrie->insert(make_pair(T3->root->key, T3));
//	mTrie->insert(make_pair(T4->root->key, T4));
//	mTrie->insert(make_pair(T5->root->key, T5));

	mapNBVer * mLV = new mapNBVer();
	list<int> L1({1});
	list<int> L2({1, 2});
	list<int> L3({2, 3});
	list<int> L4({2, 3, 4});
	list<int> L5({1, 2, 3, 4, 5});
	list<int> L6({4, 5});
	list<int> L7({5, 7});
	list<int> L8({5, 7, 8});
	list<int> L9({6});
	list<int> L10({6, 10});
	list<int> L11({4, 5, 7});

	mLV->insert(make_pair(2, L1));
	mLV->insert(make_pair(3, L2));
	mLV->insert(make_pair(4, L3));
	mLV->insert(make_pair(5, L4));
	mLV->insert(make_pair(6, L5));
	mLV->insert(make_pair(7, L6));
	mLV->insert(make_pair(8, L7));
	mLV->insert(make_pair(9, L8));
	mLV->insert(make_pair(10, L9));
	mLV->insert(make_pair(11, L10));
	mLV->insert(make_pair(12, L11));

	MCSet newMCs;
	mapMC mMC;
	int nlevels = 5;

	CliqueGraph * CG = new CliqueGraph(nlevels);

	newMCs.push_back(mc1);
	newMCs.push_back(mc2);
	newMCs.push_back(mc3);
	newMCs.push_back(mc4);
	newMCs.push_back(mc5);

//	for (auto it : newMCs) {
//		for (auto iter : it) {
//			cout << iter << " ";
//		}
//		cout << '\n';
//	}

	CG->addMCs(&newMCs, mLV, mTrie);

	for (int k = nlevels; k > 1; k--) {
		CG->checkSpanTreeConnection(k);
	}

	PD * P = CG->generatePD();

	for (auto it : *P) {
        cout << "(b,d): "<< "(" << it.first << "," << it.second << ")" << endl;
	}

	MCSet delMCs;
	delMCs.push_back(mc1);
	delMCs.push_back(mc3);
	delMCs.push_back(mc5);

	CG->removeMCs(&delMCs, mTrie);

	newMCs.clear();
	newMCs.push_back(mc6);
	newMCs.push_back(mc7);
	newMCs.push_back(mc8);
	//newMCs.push_back(mc9);
	//newMCs.push_back(mc10);

	EdgeSet * delEdges = new EdgeSet();

	delEdges->push_back(make_pair(1, 3));
	delEdges->push_back(make_pair(4, 12));
	delEdges->push_back(make_pair(5, 12));
	delEdges->push_back(make_pair(7, 12));
	delEdges->push_back(make_pair(6, 11));
	delEdges->push_back(make_pair(10, 11));

	for (auto e : *delEdges) {
		list<int> & Lv = (*mLV)[e.second];
		Lv.remove(e.first);
	}

	for (auto it : *mLV) {
		cout << "vertex: "<<it.first << endl;
		for (auto iter : it.second) {
			cout << iter<< " ";
		}
		cout << endl;
	}

	cout << "Number of Tries: " << mTrie->size() << endl;
	for (auto it : *mTrie) {
		cout << "Trie " << it.first << endl;
		if (!it.second->root->chi.empty()) {
			MCPair * MCP = it.second->outputMC(it.second->root);
			for (auto iter : *MCP) {
				cout << "MC id: " << iter.first << endl;
				for (auto vit : iter.second) {
					cout << vit << " ";
				}
				cout << endl;
			}
		}
	}

	CG->addMCs(&newMCs, mLV, mTrie);

	newMCs.clear();
	newMCs.push_back(mc11);

	EdgeSet * newEdges = new EdgeSet();
	newEdges->push_back(make_pair(6, 9));

	for (auto e : *newEdges) {
			list<int> & Lv = (*mLV)[e.second];
			Lv.push_back(e.first);
	}

	for (auto it : *mLV) {
			cout << "vertex: "<<it.first << endl;
			for (auto iter : it.second) {
				cout << iter<< " ";
			}
			cout << endl;
		}
	CG->addMCs(&newMCs, mLV, mTrie);

	for (int k = nlevels; k > 1; k--) {
		CG->checkSpanTreeConnection(k);
	}

	P = CG->generatePD();

	for (auto it : *P) {
	        cout << "(b,d): "<< "(" << it.first << "," << it.second << ")" << endl;
	}
	return 0;

}

//PD * buildCT(mapNBVer * mLV, mapNBVer * mHV) {
//	int m0 = 0;
//	mapTrie * mTrie;
//	int nlevels = 5;
//	CliqueGraph * CG = new CliqueGraph(nlevels);
//	CG->addMCs(&newMCs, mLV, mTrie);
//	for (int k = nlevels; k > 1; k--) {
//		CG->checkSpanTreeConnection(k);
//	}
//	PD * P = CG->generatePD();
//  return P;
//}

