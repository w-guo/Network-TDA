/*
 * utility.h
 *
 */

#ifndef UTILITY_H_
#define UTILITY_H_

#include "Trie.h"

typedef THash<TInt, Trie *> TIntTrieH;

bool isKeyVIn(TIntV &KeyV, TIntSet &A) {
  for (int i = 0; i < KeyV.Len(); i++) {
    if (!A.IsKey(KeyV[i])) {
      return false;
    }
  }
  return true;
}

void getRelativeComplement(const TIntSet &A, const TIntSet &B,
                           TIntSet &Complement) {
  for (TIntSet::TIter it = A.BegI(); it < A.EndI(); it++) {
    const int nId = it.GetKey();
    if (!B.IsKey(nId))
      Complement.AddKey(nId);
  }
}

void getIntersection(const TIntSet &A, const TIntSet &B, TIntSet &C) {
  if (A.Len() < B.Len()) {
    for (THashSetKeyI<TInt> it = A.BegI(); it < A.EndI(); it++)
      if (B.IsKey(it.GetKey()))
        C.AddKey(it.GetKey());
  } else {
    for (THashSetKeyI<TInt> it = B.BegI(); it < B.EndI(); it++)
      if (A.IsKey(it.GetKey()))
        C.AddKey(it.GetKey());
  }
}

void getUnion(const TIntSet &A, const TIntSet &B, TIntSet &C) {
  if (A.Len() < B.Len()) {
    C = B;
    for (THashSetKeyI<TInt> it = A.BegI(); it < A.EndI(); it++)
      if (!B.IsKey(it.GetKey()))
        C.AddKey(it.GetKey());
  } else {
    C = A;
    for (THashSetKeyI<TInt> it = B.BegI(); it < B.EndI(); it++)
      if (!A.IsKey(it.GetKey()))
        C.AddKey(it.GetKey());
  }
}

void getPairCombination(TIntV &K, TIntPrV &PrV) {
  K.Sort();
  for (int i = 0; i < K.Len() - 1; i++) {
    for (int j = i + 1; j < K.Len(); j++) {
      PrV.Add(TIntPr(K[i], K[j]));
    }
  }
}

void getNbrs(PUNGraph &G, int NId, TIntSet &Nbrs) {
  TUNGraph::TNodeI node = G->GetNI(NId);
  int deg = node.GetDeg();
  for (int i = 0; i < deg; i++)
    Nbrs.AddKey(node.GetNbrNId(i));
}

void getLNbrs(PUNGraph &G, int NId, TIntSet &LNbrs) {
  TUNGraph::TNodeI node = G->GetNI(NId);
  node.SortNIdV();
  for (int i = 0; i < node.GetDeg(); i++) {
    int nId = node.GetNbrNId(i);
    if (nId > NId) {
      return;
    }
    LNbrs.AddKey(nId);
  }
}

void getHNbrs(PUNGraph &G, int NId, TIntSet &HNbrs) {
  TUNGraph::TNodeI node = G->GetNI(NId);
  node.SortNIdV();
  for (int i = node.GetDeg() - 1; i >= 0; i--) {
    int nId = node.GetNbrNId(i);
    if (nId < NId) {
      return;
    }
    HNbrs.AddKey(nId);
  }
}

void enumPotentialSubsumedMCs(TIntV &C, TIntPrV &EdgeV, TVec<TIntV> &S) {
  S.Add(C); // C is sorted
  TIntPrV CliqueEdgeV, IntrsEdgeV;
  getPairCombination(C, CliqueEdgeV);
  CliqueEdgeV.Intrs(EdgeV, IntrsEdgeV);

  for (int i = 0; i < IntrsEdgeV.Len(); i++) {
    int u = IntrsEdgeV[i].Val1;
    int v = IntrsEdgeV[i].Val2;
    TVec<TIntV> Sp;
    for (int j = 0; j < S.Len(); j++) {
      TIntV Cp = S[j];
      if (Cp.IsIn(u) && Cp.IsIn(v)) {
        TIntV C1 = Cp;
        TIntV C2 = Cp;
        C1.DelIfIn(u);
        C2.DelIfIn(v);
        Sp.Add(C1);
        Sp.Add(C2);
      } else {
        Sp.Add(Cp);
      }
    }
    S = Sp;
  }
}

// store neighboring MCs of the given MC in nbMCs
// each MC is sorted
void getNeighborMCs(PUNGraph &G, TIntV &MC, int m, TIntTrieH &TrieH,
                    TIntIntVH &nbMCs) {
  for (TIntV::TIter it = MC.BegI(); it != MC.EndI(); it++) {
    int vId = *it;
    if (TrieH.IsKey(vId)) {
      Trie *T = TrieH.GetDat(vId);
      TIntIntVH *MCVH = T->outputMC(T->root);
      for (THashKeyDatI<TInt, TIntV> mIt = MCVH->BegI(); !mIt.IsEnd(); mIt++) {
        nbMCs.AddDat(mIt.GetKey(), mIt.GetDat());
      }
      MCVH->Clr();
    }
  }
  nbMCs.DelKey(m);

  for (TIntV::TIter it = MC.BegI(); it != MC.EndI(); it++) {
    int vId = *it;
    TIntSet vLNbrs;
    getLNbrs(G, vId, vLNbrs);
    if (vLNbrs.Len() > 0) {
      for (TIntV::TIter vIt = MC.BegI(); *vIt < vId; vIt++) {
        vLNbrs.DelKey(*vIt);
      }
      for (TIntSet::TIter uIt = vLNbrs.BegI(); uIt < vLNbrs.EndI(); uIt++) {
        int uId = *uIt;
        if (TrieH.IsKey(uId)) {
          Trie *T = TrieH.GetDat(uId);
          if (T->NodeList.IsKey(vId)) {
            TrieNode *node = T->NodeList.GetDat(vId);
            while (node != NULL) {
              TIntIntVH *MCVH = T->outputMC(node);
              for (THashKeyDatI<TInt, TIntV> mIt = MCVH->BegI(); !mIt.IsEnd();
                   mIt++) {
                nbMCs.AddDat(mIt.GetKey(), mIt.GetDat());
              }
              MCVH->Clr();
              node = node->next;
            }
          }
        }
      }
    }
  }
}

// convert MCs saved in TrieH to Hashtable with MC id as the key
void convertTrieH(TIntTrieH &TrieH, TIntIntVH &allMCH) {
  allMCH.Clr();
  for (THashKeyDatI<TInt, Trie *> TIt = TrieH.BegI(); !TIt.IsEnd(); TIt++) {
    Trie *T = TIt.GetDat();
    TIntIntVH *MCVH = T->outputMC(T->root);
    for (THashKeyDatI<TInt, TIntV> mIt = MCVH->BegI(); !mIt.IsEnd(); mIt++) {
      allMCH.AddDat(mIt.GetKey(), mIt.GetDat());
    }
    MCVH->Clr();
  }
}

// convert a vector of edges to a undirected graph
PUNGraph convertEdgeVToGraph(TIntPrV &EdgeV) {
  PUNGraph G = TUNGraph::New();
  for (int e = 0; e < EdgeV.Len(); e++) {
    if (!G->IsNode(EdgeV[e].Val1)) {
      G->AddNode(EdgeV[e].Val1);
    }
    if (!G->IsNode(EdgeV[e].Val2)) {
      G->AddNode(EdgeV[e].Val2);
    }
    G->AddEdge(EdgeV[e].Val1, EdgeV[e].Val2);
  }
  return G;
}

/*
 * Implementation of the 2-opt algorithm for a minimum weighted vertex cover by
 * R. Bar-Yehuda and S. Even. A linear time approximation algorithm for the
 * weighted vertex cover problem. J. of Algorithms 2:198-203, 1981. The solution
 * is guaranteed to be within 2 times the optimum solution.
 */
int BarYehudaEvenTwoApproxVC(PUNGraph &G) {
  int cnt = 0;

  while (G->GetEdges() > 0) {
    TUNGraph::TEdgeI EI = G->BegEI();
    int u = EI.GetSrcNId();
    // printf("include vertex id %d in the cover set\n", u);
    cnt++;
    // remove all the edges incident to u
    G->DelNode(u);
  }
  // return the size of approx minimum vertex cover set
  return cnt;
}

// compute an upper bound on the TSN
int TSNUpperBound(TIntPrV &NewEdgeV) {
  // PUNGraph G1 = convertEdgeVToGraph(DelEdgeV);
  PUNGraph G2 = convertEdgeVToGraph(NewEdgeV);

  // upper bound for RSN
  // int RSNUB = BarYehudaEvenTwoApproxVC(G1);
  // upper bound for ASN
  int ASNUB = BarYehudaEvenTwoApproxVC(G2);
  // int TSNUB = RSNUB + ASNUB;
  // G1->Defrag();
  G2->Defrag();

  return ASNUB;
}

void saveTxt(const TVec<TIntV> &cmtyVV, const TStr &FNm) {
  FILE *F = fopen(FNm.CStr(), "wt");

  for (int i = 0; i < cmtyVV.Len(); i++) {
    for (int j = 0; j < cmtyVV[i].Len(); j++) {
      fprintf(F, "%d\t", cmtyVV[i][j]);
    }
    fprintf(F, "\n");
  }
  fclose(F);
}

#endif /* UTILITY_H_ */