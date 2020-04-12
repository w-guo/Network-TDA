/*
 * CliqueGraph.h
 */

#ifndef CLIQUEGRAPH_H_
#define CLIQUEGRAPH_H_

#include "EulerTourTree.h"
#include "utility.h"

using namespace std;

typedef THash<TInt, ETNode *> ETForest;
typedef TVec<ETForest> ETForestVec;

struct MC {

  MC(int m, int MCsize)
      : MCid(m), birthT(MCsize), deathT(-1), visited(false), child() {
    pCG.Reserve(MCsize - 1);
  }

  ~MC() {
    for (TVec<CGNode *>::TIter it = pCG.BegI(); it < pCG.EndI(); it++) {
      CGNode *cg = *it;
      delete cg->loopNode;
      delete cg;
    }
  }

public:
  int MCid, birthT, deathT;
  bool visited;
  TVec<MC *> child;
  CGNodeArr pCG;
};

typedef THash<TInt, MC *> TIntMCH;

// a node in the Clique graph is a representation of a MC and
// corresponds an ET Node (and there is an ET forest for each level)
class CliqueGraph {

public:
  CliqueGraph(int nlevels) : m0(0) { forests.Gen(nlevels); }

  ~CliqueGraph() { clear(); };

public:
  void init(MC *M) {
    int m = M->MCid;
    int nlevels = M->birthT;

    for (int i = 0; i < nlevels - 1; i++) {
      CGNode *cg = new CGNode(m, nlevels);
      EulerTourTree::createTour(cg);
      (M->pCG)[i] = cg;
      forests[i + 1].AddDat(m, cg->loopNode);
    }
  };

  void clear() {
    for (TVec<ETForest>::TIter it = forests.BegI(); it < forests.EndI(); it++) {
      ETForest f = *it;
      for (THashKeyDatI<TInt, ETNode *> fIt = f.BegI(); !fIt.IsEnd(); fIt++) {
        EulerTourTree::destroyTour(fIt.GetDat()->findRoot());
      }
      f.Clr();
    }
    forests.Clr();

    for (THashKeyDatI<TInt, MC *> MIt = MCH.BegI(); !MIt.IsEnd(); MIt++) {
      delete MIt.GetDat();
    }
    MCH.Clr();
  }

  void addMCs(PUNGraph &G, TVec<TIntV> &newMCs, TIntTrieH &TrieH,
              TVec<TIntH> &prevCG) {
    int m = m0;
    // initialize the clique graph at each level
    for (TVec<TIntV>::TIter it = newMCs.BegI(); it != newMCs.EndI(); it++) {
      TIntV MCV = *it;
      m++;
      // add each MC to the hash table of Trie
      int min_vId = MCV[0]; // each MC is sorted
      if (TrieH.IsKey(min_vId)) {
        Trie *T = TrieH.GetDat(min_vId);
        T->addMC(MCV, m);
      } else {
        Trie *T = new Trie(min_vId);
        T->addMC(MCV, m);
        TrieH.AddDat(min_vId, T);
      }
      // add each MC to MC map
      MC *M = new MC(m, MCV.Len());
      MCH.AddDat(m, M);
      init(M);
    }

    for (TVec<TIntV>::TIter it = newMCs.BegI(); it != newMCs.EndI(); it++) {
      TIntV MCV = *it;
      int min_vId = MCV[0]; // each MC is sorted
      Trie *T = TrieH.GetDat(min_vId);
      m = T->getID(MCV);
      MC *M = MCH.GetDat(m);

      TIntIntVH nbMCs;
      getNeighborMCs(G, MCV, m, TrieH, nbMCs);

      if (nbMCs.Len() > 0) {
        for (THashKeyDatI<TInt, TIntV> mIt = nbMCs.BegI(); !mIt.IsEnd();
             mIt++) {
          MC *nbM = MCH.GetDat(mIt.GetKey());
          if (nbM->visited and mIt.GetKey() > m0) {
            continue;
          } else {
            if (TMath::Mn(M->birthT, nbM->birthT) == 2) {
              insertCGEdges(&M->pCG, &nbM->pCG, 1, prevCG);
            } else {
              TIntV nbMCV = mIt.GetDat();
              TIntV interMCV;
              MCV.Intrs(nbMCV, interMCV);
              insertCGEdges(&M->pCG, &nbM->pCG, interMCV.Len(), prevCG);
            }
          }
        }
      }
      M->visited = true;
    }
    m0 = m;
  }

  void removeMC(TIntV &MCV, TIntTrieH &TrieH) {
    int min_vId = MCV[0]; // each MC is sorted
    Trie *T = TrieH.GetDat(min_vId);
    int m = T->getID(MCV);
    MC *M = MCH.GetDat(m);
    deleteCGEdges(M);

    for (int k = MCV.Len() - 1; k > 0; k--) {
      forests[k].DelIfKey(m);
    }
    // remove this MC from the hash table and delete it
    MCH.DelKey(m);
    delete M;
    M = NULL;

    T->removeMC(MCV);
    if (T->NodeList.Empty()) {
      TrieH.DelKey(min_vId);
    }
  }

  void removeMCs(TVec<TIntV> &delMCs, TIntTrieH &TrieH) {
    for (TVec<TIntV>::TIter it = delMCs.BegI(); it != delMCs.EndI(); it++) {
      TIntV MCV = *it;
      removeMC(MCV, TrieH);
    }
  }

  void insertCGEdges(CGNodeArr *cga_u, CGNodeArr *cga_v, int w,
                     TVec<TIntH> &prevCG) {
    CGEdge *cg_e = new CGEdge(cga_u, cga_v, w);

    for (int l = w; l > 0; l--) {
      CGNode *cg_u = (*cga_u)[l - 1];
      CGNode *cg_v = (*cga_v)[l - 1];
      // if they are already connected
      if (EulerTourTree::connected(cg_u, cg_v)) {
        return;
      } else {
        ETNode *et_ru = cg_u->loopNode->splay();
        ETNode *et_rv = cg_v->loopNode->splay();
        // ETNode *et_u = cg_u->loopNode;
        // ETNode *et_v = cg_v->loopNode;
        updateTreeEdgeInsertion(et_ru, et_rv, l, prevCG);
        addTreeEdge(cg_u, cg_v, cg_e, l - 1);
      }
    }
  }

  void deleteCGEdges(MC *M) {
    CGNodeArr *cga = &M->pCG;
    for (int l = M->birthT - 1; l > 0; l--) {
      while (!(*cga)[l - 1]->adjTreeCGEdges.Empty()) {
        CGEdge *cg_e = (*cga)[l - 1]->adjTreeCGEdges.Last();
        deleteTreeEdge((*cg_e->cga1)[l - 1], (*cg_e->cga2)[l - 1], cg_e, l - 1);
        if (l == 1) {
          delete cg_e;
        }
      }
    }
  }

  ETNode *addTreeEdge(CGNode *cg_u, CGNode *cg_v, CGEdge *cg_e, int l) {
    // put edge into the adjacency lists (for quick iteration)
    cg_u->adjTreeCGEdges.Add(cg_e);
    cg_v->adjTreeCGEdges.Add(cg_e);

    ETNode *neuv = NULL, *nevu = NULL;
    ETNode *newroot = EulerTourTree::link(cg_u, cg_v, &neuv, &nevu);
    // CGNode *cg_1 = (*cg_e->cga1)[0];
    // CGNode *cg_2 = (*cg_e->cga2)[0];
    // printf("new ET tree after adding edge (%d,%d): ", cg_1->cgId,
    // cg_2->cgId); EulerTourTree::print(newroot); cout << endl;

    // store the pointers in the edge
    cg_e->arcs[l].Val1 = neuv;
    cg_e->arcs[l].Val2 = nevu;

    return newroot;
  }

  void deleteTreeEdge(CGNode *cg_u, CGNode *cg_v, CGEdge *cg_e, int l) {
    // remove edge from adjacency lists
    removeFromAdjList(cg_u->adjTreeCGEdges, cg_u->cgId, cg_v->cgId);
    removeFromAdjList(cg_v->adjTreeCGEdges, cg_u->cgId, cg_v->cgId);

    ETNode *et_ru = NULL, *et_rv = NULL;
    EulerTourTree::cut(cg_e->arcs[l].Val1, cg_e->arcs[l].Val2, &et_ru, &et_rv);

    // CGNode *cg_1 = (*cg_e->cga1)[0];
    // CGNode *cg_2 = (*cg_e->cga2)[0];
    // printf("split ET trees after removing edge (%d,%d)\n", cg_1->cgId,
    // cg_2->cgId); EulerTourTree::print(et_ru); cout << endl;
    // EulerTourTree::print(et_rv);
    // cout << endl;
  }

  void removeFromAdjList(AdjacencyList &L, int m1, int m2) {
    for (int i = 0; i < L.Len(); i++) {
      CGNode *cg_1 = (*L[i]->cga1)[0];
      CGNode *cg_2 = (*L[i]->cga2)[0];
      if ((cg_1->cgId == m1 && cg_2->cgId == m2) ||
          (cg_1->cgId == m2 && cg_2->cgId == m1)) {
        L.Del(i);
        break;
      }
    }
  }

  void updateTreeEdgeInsertion(ETNode *et_ru, ETNode *et_rv, int w,
                               TVec<TIntH> &prevCG) {
    // representative ETNode IDs
    int uId = et_ru->dw;
    int vId = et_rv->dw;
    // cout << "first representative node ID: " << uId << endl;
    // cout << "second representative node ID: " << vId << endl;
    int uSize, vSize;
    if (uId != vId) {
      if (forests[w].IsKey(uId)) {
        ETNode *et_u = forests[w].GetDat(uId);
        uSize = et_u->cg1->size;
      } else {
        uSize = prevCG[w].GetDat(uId);
      }
      if (forests[w].IsKey(vId)) {
        ETNode *et_v = forests[w].GetDat(vId);
        vSize = et_v->cg1->size;
      } else {
        vSize = prevCG[w].GetDat(vId);
      }
      if ((uSize > vSize) || (uSize == vSize && uId < vId)) {
        et_rv->dw = uId;
        // cout << "removed representative node ID: " << vId << endl;
        forests[w].DelIfKey(vId);
      } else {
        et_ru->dw = vId;
        // cout << "removed representative node ID: " << uId << endl;
        forests[w].DelIfKey(uId);
      }
    }
  }

  // generate persistence diagram
  TIntPrV *generatePD() {
    THashKeyDatI<TInt, ETNode *> it = forests[1].BegI();
    ETNode *et_m = it.GetDat();

    for (THashKeyDatI<TInt, ETNode *> it = forests[1].BegI();
         it < forests[1].EndI(); it++) {
      if (it.GetDat()->cg1->size > et_m->cg1->size) {
        et_m = it.GetDat();
      }
    }

    // cout << "ETNode at level 0: ";
    // et_m->print();
    // cout << endl;
    forests[0].Clr();
    forests[0].AddDat(et_m->cg1->cgId, et_m);

    for (int k = forests.Len() - 1; k > 0; k--) {
      recordDeathTime(k);
    }

    MC *M = MCH.GetDat(et_m->cg1->cgId);
    M->deathT = 1;

    TIntPrV *P = new TIntPrV();
    preOrder(P, M);

    return P;
  }

  void recordDeathTime(int k) {
    ETForest D = forests[k];
    for (THashKeyDatI<TInt, ETNode *> fIt = forests[k - 1].BegI(); !fIt.IsEnd();
         fIt++) {
      D.DelIfKey(fIt.GetKey());
    }

    if (D.Len() > 0) {
      for (THashKeyDatI<TInt, ETNode *> dIt = D.BegI(); !dIt.IsEnd(); dIt++) {
        MC *M = MCH.GetDat(dIt.GetKey());
        // printf("child MC id: %d merged at level %d\n", M->MCid, k);
        M->deathT = k;
        if (k > 1) {
          for (THashKeyDatI<TInt, ETNode *> fIt = forests[k - 1].BegI();
               !fIt.IsEnd(); fIt++) {
            if (EulerTourTree::connected((M->pCG)[k - 2]->loopNode,
                                         fIt.GetDat())) {
              MC *M_par = MCH.GetDat(fIt.GetKey());
              // printf("its parent id: %d\n", M_par->MCid);
              M_par->child.Add(M);
              break;
            }
          }
        } else {
          THashKeyDatI<TInt, ETNode *> fIt = forests[0].BegI();
          MC *M_par = MCH.GetDat(fIt.GetKey());
          M_par->child.Add(M);
        }
      }
    }
  }

  void preOrder(TIntPrV *P, MC *M) {
    if (!M)
      return;
    else {
      P->Add(TIntPr(M->birthT, M->deathT));
      while (!M->child.Empty()) {
        MC *M_chi = M->child.Last();
        M->child.DelLast();
        preOrder(P, M_chi);
      }
    }
  }

  void enumVerticesInCommunity(ETNode *n, TIntIntVH &allMCH, TIntSet &CmtySet) {
    n = (ETNode *)n->findRoot()->subtree_minimum();
    while (n != NULL) {
      if (n->isLoop()) {
        TIntV MCV = allMCH.GetDat(n->cg1->cgId);
        CmtySet.AddKeyV(MCV);
      }
      n = (ETNode *)n->next();
    }
  }

  void convertCommunities(TIntTrieH &TrieH, THash<TInt, TVec<TIntV>> &cmtyH) {
    TIntIntVH allMCH;
    convertTrieH(TrieH, allMCH);

    for (int l = 1; l < forests.Len(); l++) {
      TVec<TIntV> NIdCmtyVV;
      TIntSet CmtySet;
      for (THashKeyDatI<TInt, ETNode *> fIt = forests[l].BegI();
           fIt < forests[l].EndI(); fIt++) {
        CmtySet.Clr(false);
        ETNode *n = fIt.GetDat();
        enumVerticesInCommunity(n, allMCH, CmtySet);
        NIdCmtyVV.Add();
        CmtySet.GetKeyV(NIdCmtyVV.Last());
        NIdCmtyVV.Last().Sort();
      }
      cmtyH.AddDat(l + 1, NIdCmtyVV);
    }
  }

private:
  // current maximum id of MCs
  int m0;

public:
  // the set of MCs that each one is represented by an array of CG nodes (one
  // for each level)
  TIntMCH MCH;

  // the set of forests (one for each level)
  ETForestVec forests;
};

#endif // CLIQUEGRAPH_H_