/*
 * CliqueGraph.h
 */

#ifndef CLIQUEGRAPH_H_
#define CLIQUEGRAPH_H_

#include "EulerTourTree.h"
#include "utility.h"
#include <vector>
#include <map>
#include <list>

typedef map<int, ETNode *> ETForest;
typedef vector<ETForest> ETForestVec;
typedef map<int, CGNodeArr *> CGNodes;

typedef pair<int, int> bdPair;
typedef vector<bdPair> PD;

struct MC {

	MC(int m, int MCsize)
		: MCid(m)
		, birthT(MCsize)
		, deathT(-1)
		, visited(false)
		, child()
		, pCG(NULL)
		{}

	~MC() {
		for (auto cg: *pCG) {
			cg->adjNontreeCGEdges.clear();
			cg->adjTreeCGEdges.clear();
			delete cg;
		}
	}

public:

	int MCid, birthT, deathT;
	bool visited;
	vector<MC *> child;
	CGNodeArr * pCG;
};

typedef map<int, MC *> mapMC;

// a node in the Clique graph is a representation of a MC and
// corresponds an ET Node (and there is an ET forest for each level)
class CliqueGraph {

public:
	CliqueGraph(int nlevels)
		: m0(0)
	{
		forests.resize(nlevels);
	}

	~CliqueGraph() {
		clear();
	};

public:
	void init(MC * M) {
		int m = M->MCid;
		int nlevels = M->birthT;
		CGNodeArr * cga = (CGNodeArr *)malloc((nlevels-1)*sizeof(CGNode));
		for (int i = 0; i < nlevels-1; i++) {
			CGNode * cg = new CGNode(m, nlevels);
			EulerTourTree::createTour(cg);
			cga->emplace_back(cg);
			forests[i+1].insert(make_pair(m, cg->loopNode));
		}
		M->pCG = cga;
		//cout << "cga address: " << cga << endl;
	};


	void clear() {
		for (auto f : forests) {
			for (auto fit : f) {
				EulerTourTree::destroyTour(fit.second->findRoot());
			}
			f.clear();
		}
		forests.clear();

//		for (auto it : CGnodes) {
//			for (CGNode * cg : *it.second) {
//				cg->adjNontreeCGEdges.clear();
//				cg->adjTreeCGEdges.clear();
//			}
//		}
//		CGnodes.clear();
		for (auto Mit : mMC) {
			delete Mit.second;
		}
		mMC.clear();
	}


	void addMCs(MCSet * newMCs, mapNBVer * mLV, mapTrie * mTrie) {
		int m = m0;
		// initialize the clique graph at each level
		for (auto MCvec : *newMCs) {
			m++;
			// add each MC to Trie map
			int min_vid = *min_element(MCvec.begin(), MCvec.end());
			if (mTrie->find(min_vid) != mTrie->end()) {
				Trie * T = (*mTrie)[min_vid];
				T->addMC(&MCvec, m);
			}
			else {
				Trie * T = new Trie(min_vid);
				T->addMC(&MCvec, m);
				mTrie->insert(make_pair(min_vid, T));
			}
			if (MCvec.size() > 1) {
				// add each MC to MC map
				MC * M = new MC(m, MCvec.size());
				mMC.insert(make_pair(m, M));
				init(M);
//				for (int k = MCvec.size()-1; k > 0; k--) {
//					//	cout << "cga address: " << M->pCG << endl;
//					//	cout << "cg address: " << (*M->pCG)[k-1] << endl;
//					//	cout << "cg address2: " << (*M->pCG)[k-1]->loopNode->cg1 << endl;
//					forests[k].insert(make_pair(m, (*M->pCG)[k-1]->loopNode));
//				}
			}
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

		for (auto MCvec : *newMCs) {
			int min_vid = *min_element(MCvec.begin(), MCvec.end());
			Trie * T = (*mTrie)[min_vid];
			m = T->getID(&MCvec);
			MC * M = mMC[m];
			cout << "add MC id: " << M->MCid << endl;
			MCPair * nbSet = getNeighborMCs(&MCvec, m, mLV, mTrie);
			for (auto iter : *nbSet) {
				cout << "neighbor MC id: " << iter.first << endl;
				for (auto vit : iter.second) {
					cout << vit << " ";
				}
				cout << endl;
			}
			if (nbSet->size() > 0) {
				for (auto nbMCPair : *nbSet) {
					MC * nbM = mMC[nbMCPair.first];
					if (nbM->visited and nbMCPair.first > m0) {
						continue;
					}
					else {
						int k = min(M->birthT, nbM->birthT);
						if (k == 2) {
							insertCGEdges(M->pCG, nbM->pCG, 1, 2);

						}
						else {
							MCVec nbMC = nbMCPair.second, interMC;
							set_intersection(MCvec.begin(), MCvec.end(), nbMC.begin(), nbMC.end(), back_inserter(interMC));
							insertCGEdges(M->pCG, nbM->pCG, interMC.size(), k);
						}
					}
				}
			}
			M->visited = true;
		}
		m0 = m;
	}


	void removeMCs(MCSet * delMCs, mapTrie * mTrie) {
		for (auto MCvec : *delMCs) {
			int min_vid = *min_element(MCvec.begin(), MCvec.end());
			Trie * T = (*mTrie)[min_vid];
			if (MCvec.size() > 1) {
				int m = T->getID(&MCvec);
				MC * M = mMC[m];
				cout << "remove MC id: " << M->MCid << endl;
				deleteCGEdges(M->pCG);
				for (int k = MCvec.size()-1; k > 0; k--) {
					forests[k].erase(m);
				}
				if (M->deathT == 1) {
					forests[0].clear();
				}
				// remove this MC from the map and delete it
				mMC.erase(m);
				for (auto cg : *M->pCG) {
					delete cg->loopNode;
				}
				delete M;
				M = NULL;
			}
			T->removeMC(&MCvec);
			if (T->root->chi.empty()) {
				mTrie->erase(min_vid);
			}
		}
	}


	void insertCGEdges(CGNodeArr * cga_u, CGNodeArr * cga_v, int w, int k) {

//		CGEdge * cg_e = NULL;
//		if (!edgepool.empty()) {
//			cg_e = edgepool.back();
//			cg_e->set(cga_u, cga_v, w);
//			edgepool.pop_back();
//		}
//		else {
//			cg_e = new CGEdge(cga_u, cga_v, w);
//			edgepool.push_back(cg_e);
//		}
		CGEdge * cg_e = new CGEdge(cga_u, cga_v, w);
		//edgepool.push_back(cg_e);

		CGNode * cg_u = (*cga_u)[w-1];
		CGNode * cg_v = (*cga_v)[w-1];

		// if they're already connected or should not be connected by a tree edge for now
		if (EulerTourTree::connected(cg_u, cg_v) || w != k-1) {
			// add edge to non-tree adjacency lists and update weights
			for (int l = w; l > 0; l--) {
				addNontreeEdge((*cga_u)[l-1], (*cga_v)[l-1], cg_e);
			}
		}
		else {
			ETNode * et_ru = cg_u->loopNode->splay();
			ETNode * et_rv = cg_v->loopNode->splay();
			updateTreeEdgeInsertion(et_ru, et_rv, w);
			for (int l = w; l > 0; l--) {
				addTreeEdge((*cga_u)[l-1], (*cga_v)[l-1], cg_e, l-1);
//				EulerTourTree::print(newroot);
//				cout << endl;
			}
		}

	}


	void deleteCGEdges(CGNodeArr * cga) {
		cout << "Number of non-tree edges at level 0: " <<(*cga)[0]->adjNontreeCGEdges.size() << endl;
		if (!(*cga)[0]->adjNontreeCGEdges.empty()) {
			for (auto e : (*cga)[0]->adjNontreeCGEdges) {
				for (int l = e->weight; l > 0; l--) {
					deleteNontreeEdge((*e->cga1)[l-1], (*e->cga2)[l-1], e);
				}
				delete e;
			}
		}

		cout << "Number of tree edges at level 0: " <<(*cga)[0]->adjTreeCGEdges.size() << endl;
		if (!(*cga)[0]->adjTreeCGEdges.empty()) {
			for (auto e : (*cga)[0]->adjTreeCGEdges) {
				int w = e->weight;
				for (int l = w; l > 0; l--) {
					deleteTreeEdge((*e->cga1)[l-1], (*e->cga2)[l-1], e, l-1);
				}
				ETNode * et_ru = (*e->cga1)[w-1]->loopNode->splay();
				ETNode * et_rv = (*e->cga2)[w-1]->loopNode->splay();
				//cout << "root ETNode of the first ET tree: "; et_ru->print(); cout << endl;
				//cout << "root ETNode of the second ET tree: "; et_rv->print(); cout << endl;
				updateTreeEdgeDeletion(et_ru, et_rv, w);
				delete e;
			}
		}
	};


	void addNontreeEdge(CGNode * cg_u, CGNode * cg_v, CGEdge * cg_e) {
		// put edge into the adjacency lists
		cg_u->adjNontreeCGEdges.push_back(cg_e);
		cg_v->adjNontreeCGEdges.push_back(cg_e);

		// update the weight of the relevant node (weight == number of non-tree edges)
		cg_u->loopNode->addWeight(1);
		cg_v->loopNode->addWeight(1);

		CGNode * cg_1 = (*cg_e->cga1)[0];
		CGNode * cg_2 = (*cg_e->cga2)[0];
		cout << "added non-tree edge (" << cg_1->cgid << "," << cg_2->cgid << ")" << endl;
	}


	ETNode * addTreeEdge(CGNode * cg_u, CGNode * cg_v, CGEdge * cg_e, int l) {
		// put edge into the adjacency lists (for quick iteration)
		cg_u->adjTreeCGEdges.push_back(cg_e);
		cg_v->adjTreeCGEdges.push_back(cg_e);

		ETNode * neuv = NULL, * nevu = NULL;
		ETNode * newroot = EulerTourTree::link(cg_u, cg_v, &neuv, &nevu);
		CGNode * cg_1 = (*cg_e->cga1)[0];
		CGNode * cg_2 = (*cg_e->cga2)[0];
		cout << "new ET tree after adding edge (" << cg_1->cgid << "," << cg_2->cgid << "): ";
		EulerTourTree::print(newroot); cout << endl;

		// store the pointers in the edge
		if ( l >= cg_e->arcs.size() ) {
			cg_e->arcs.resize(l+1);
		}
		cg_e->arcs[l].first = neuv;
		cg_e->arcs[l].second = nevu;

		return newroot;
	}


	void deleteNontreeEdge(CGNode * cg_u, CGNode * cg_v, CGEdge * cg_e) {

		// remove edge from adjacency lists
		removeFromAdjList(cg_u->adjNontreeCGEdges, cg_u->cgid, cg_v->cgid);
		removeFromAdjList(cg_v->adjNontreeCGEdges, cg_u->cgid, cg_v->cgid);

		// update the weight of the relevant node (weight == number of non-tree edges)
		cg_u->loopNode->addWeight(-1);
		cg_v->loopNode->addWeight(-1);

		CGNode * cg_1 = (*cg_e->cga1)[0];
		CGNode * cg_2 = (*cg_e->cga2)[0];
		cout << "removed non-tree edge (" << cg_1->cgid << "," << cg_2->cgid << ")" << endl;
	}


	void deleteTreeEdge(CGNode * cg_u, CGNode * cg_v, CGEdge * cg_e, int l) {
		// remove edge from adjacency lists
		removeFromAdjList(cg_u->adjTreeCGEdges, cg_u->cgid, cg_v->cgid);
		removeFromAdjList(cg_v->adjTreeCGEdges, cg_u->cgid, cg_v->cgid);

		ETNode * et_ru = NULL, * et_rv = NULL;
		EulerTourTree::cut(cg_e->arcs[l].first, cg_e->arcs[l].second, &et_ru, &et_rv);

		CGNode * cg_1 = (*cg_e->cga1)[0];
		CGNode * cg_2 = (*cg_e->cga2)[0];
		cout << "split ET trees after removing edge (" << cg_1->cgid << "," << cg_2->cgid << ")" << endl;
		EulerTourTree::print(et_ru); cout << endl;
		EulerTourTree::print(et_rv); cout << endl;

	}


	void removeFromAdjList(AdjacencyList & L, int m1, int m2) {
		for (int i = 0; i < L.size(); i++) {
			CGNode * cg_1 = (*L[i]->cga1)[0];
			CGNode * cg_2 = (*L[i]->cga2)[0];
			if (cg_1->cgid == m1 and cg_2->cgid == m2) {
				L.erase(L.begin()+i);
				break;
			}
		}
	}


	void updateTreeEdgeInsertion(ETNode * et_ru, ETNode * et_rv, int w) {
		// representative ETNode IDs
		int uid = et_ru->dw[2];
		int vid = et_rv->dw[2];
		cout << "first rep ETNode ID: " << uid << endl;
		cout << "second rep ETNode ID: " << vid << endl;

		ETNode * et_u = forests[w][uid];
		ETNode * et_v = forests[w][vid];

		if (et_u->cg1->size >= et_v->cg1->size) {
			for (int l = w; l > 0; l--) {
				et_rv = forests[l][vid]->splay();
				et_rv->dw[2] = uid;
				forests[l].erase(vid);
			}
		}
		else {
			for (int l = w; l > 0; l--) {
				et_ru = forests[l][uid]->splay();
				et_ru->dw[2] = vid;
				forests[l].erase(uid);
			}
		}

		for (int l = w; l > 0; l--) {
			cout << "At level = " << l << endl;
			for (auto it : forests[l]) {
				cout << "ETNode: "; it.second->print(); cout << endl;
			}
		}

	}

	void updateTreeEdgeDeletion(ETNode * et_ru, ETNode * et_rv, int w) {
			ETNode * et_v = NULL;

			if (forests[w].find(et_ru->dw[2]) != forests[w].end()) {
				et_v = EulerTourTree::findRep(et_rv);
			}
			else {
				et_v = EulerTourTree::findRep(et_ru);
			}

			int vid = et_v->cg1->cgid;
			CGNodeArr * cga_v = mMC[vid]->pCG;

			for (int l = w; l > 0; l--) {
				et_v = (*cga_v)[l-1]->loopNode;
				et_v->dw[2] = vid;
				forests[l].insert(make_pair(vid, et_v));
			}
		}

//	void updateTreeEdgeDeletion(ETNode * et_ru, ETNode * et_rv, int w) {
//		ETNode * et_u = forests[w][et_ru->dw[2]];
//		ETNode * et_v = NULL;
//
//		if (et_u->findRoot() == et_ru) {
//			et_v = EulerTourTree::findRep(et_rv);
//		}
//		else {
//			et_v = EulerTourTree::findRep(et_ru);
//		}
//		et_v->dw[2] = et_v->cg1->cgid;
//
//		for (int l = w; l > 0; l--) {
//			forests[l].insert(make_pair(et_v->cg1->cgid, et_v));
//		}
//	}


	void checkSpanTreeConnection(int k) {
		int i = 0;
		ETNodeStack ETnodes;

		while (i < forests[k-1].size()-1) {
			cout << "forest size is " << forests[k-1].size() << " at level " << k-1 << endl;
			ETnodes.clear();
			for (auto it : forests[k-1]) {
				ETnodes.push_back(it.second);
				cout << "ETNode: "; it.second->print(); cout << endl;
			}
			for (int j = i+1; j < forests[k-1].size(); j++) {
				cout << "i = " << i << ", " << "j = " << j << endl;
				ETNode * et_ru = ETnodes[i]->splay();
				ETNode * et_rv = ETnodes[j]->splay();
				cout << "ET tree rooted at "; EulerTourTree::print(et_ru); cout << " has " << et_ru->stweight << " non-tree edges"<< endl;
				cout << "ET tree rooted at "; EulerTourTree::print(et_rv); cout << " has " << et_rv->stweight << " non-tree edges"<< endl;
				if (et_ru->stweight > 0 and et_rv->stweight > 0) {
					if (et_ru->stweight > et_rv->stweight) {
						swap(et_ru, et_rv);
					}
					if (findSpanTreeConnection(et_ru, et_rv, k)) {
						i--;
						break;
					}
				}
			}
			i++;
		}
	}


	bool findSpanTreeConnection(ETNode * et_ru, ETNode * et_rv, int k) {
		ETNodeStack ETnodes;
		ETNode * et_m = et_ru;
		while (et_m) {
			if (et_m->isLoop() and et_m->stweight > 0) {
				cout << "ETNode "; et_m->print(); cout << " with edge number = " << et_m->stweight <<endl;
				for (auto e : et_m->cg1->adjNontreeCGEdges) {
					cout << "weight = " << e->weight << endl;
					if (e->weight == k-1) {
						if ((*e->cga2)[k-2]->cgid == et_m->cg1->cgid) {
							swap(e->cga1, e->cga2);
						}
						if ((*e->cga2)[k-2]->loopNode->findRoot() == et_rv) {
							updateTreeEdgeInsertion(et_ru, et_rv, k-1);
							for (int l = k-1; l > 0; l--) {
								CGNode * cg_u = (*e->cga1)[l-1];
								CGNode * cg_v = (*e->cga2)[l-1];
								deleteNontreeEdge(cg_u, cg_v, e);
								addTreeEdge(cg_u, cg_v, e, l-1);
							}
							return true;
						}
					}
				}
			}
			et_m = et_m->nextPreOrder(ETnodes, true);
		}
		return false;
	}

	//generate persistence diagram
	PD * generatePD() {
		ETNode * et_m = forests[1].begin()->second;
		for (auto it : forests[1]) {
			if (it.second->cg1->size > et_m->cg1->size) {
				et_m = it.second;
			}
		}
		cout << "ETNode at level 0: "; et_m->print(); cout << endl;
		forests[0].insert(make_pair(et_m->cg1->cgid, et_m));

		for (int k = forests.size()-1; k > 0; k--) {
			recordDeathTime(k);
		}

		MC * M = mMC[et_m->cg1->cgid];
		M->deathT = 1;
		PD * P = new PD();
		preOrder(P, M);

		return P;
	}


	void recordDeathTime(int k) {
		ETForest D = forests[k];
		for (auto fit : forests[k-1]) {
			D.erase(fit.first);
		}
		cout << "k = " << k << endl;
		if (D.size() > 0) {
			for (auto dit : D) {
				cout << "child MC id: " << dit.first << endl;
				MC * M = mMC[dit.first];
				M->deathT = k;
				for (auto fit : forests[k-1]) {
					if (EulerTourTree::connected((*M->pCG)[k-2]->loopNode, fit.second)) {
						cout << "parent MC id: " << fit.first << endl;
						MC * M_par = mMC[fit.first];
						M_par->child.push_back(M);
						break;
					}
				}
			}
		}
	}


	void preOrder(PD * P, MC * M) {
		if (!M) return;
		else {
			P->push_back(make_pair(M->birthT, M->deathT));
			while (!M->child.empty()) {
				MC * M_chi = M->child.back();
				M->child.pop_back();
				preOrder(P, M_chi);
			}
		}
	}
public:
	// current maximum id of MCs
	int m0;

	// the set of MCs that each one is represented by an array of CG nodes (one for each level)
	mapMC mMC;

	// the set of forests (one for each level)
	ETForestVec forests;
};

#endif // CLIQUEGRAPH_H_
