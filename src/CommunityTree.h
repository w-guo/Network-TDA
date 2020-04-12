/*
 * CommunityTree.h
 */

#ifndef COMMUNITYTREE_H_
#define COMMUNITYTREE_H_

#include "CliqueGraph.h"

class CommunityTree {
public:
  CommunityTree() : copyG(), MCV(), MCVs(NULL), maxMCSize(1), time(0) {}
  ~CommunityTree() {}

private:             // recursion variables only
  PUNGraph copyG;    // a copy of the graph
  TIntV MCV;         // vector of a maximum clique
  TVec<TIntV> *MCVs; // vector of maximum cliques
  int maxMCSize;     // max size of maximum cliques
public:
  int time;

private:
  int maxNbrsInCandNodeId(const TIntSet &SUBG, const TIntSet &CAND) const;
  int getNodeIdWithMaxDeg(const TIntSet &Set) const;
  void expandClique(const TIntSet &SUBG, TIntSet &CAND);
  bool excludeEdgesInNonMaxClique(TIntV &K, TIntPrV &ExcludeEdgeV);
  void TTTExcludeEdges(PUNGraph &G, TIntV &K, TIntSet &Cand, TIntSet &Fini,
                       TIntPrV &ExcludeEdgeV);
  bool isMC(TIntV &C);

public:
  void getMaximalCliques(const PUNGraph &G, TVec<TIntV> &initMCs);
  static void getMCs(const PUNGraph &G, TVec<TIntV> &initMCs);
  void enumNewMCs(const PUNGraph &newG, TIntPrV &NewEdgeV, TVec<TIntV> &newMCs);
  void enumSubsumedMCsEdgeInsertion(TIntPrV &NewEdgeV, TIntTrieH &TrieH,
                                    TVec<TIntV> &newMCs, TVec<TIntV> &delMCs);
  void enumSubsumedMCsEdgeDeletion(PUNGraph &G, TIntPrV &DelEdgeV,
                                   TVec<TIntV> &delMCs, TVec<TIntV> &newMCs);
  CliqueGraph *buildCT(PUNGraph &G, TIntTrieH &TrieH);
  void updateCTEI(PUNGraph &G, TIntPrV &newEdgeV, TIntTrieH &TrieH,
                  CliqueGraph *CG);
};

int CommunityTree::maxNbrsInCandNodeId(const TIntSet &SUBG,
                                       const TIntSet &CAND) const {
  int id = -1;
  int maxIntersection = -1;
  //
  for (THashSetKeyI<TInt> it = SUBG.BegI(); it < SUBG.EndI(); it++) {
    int nId = it.GetKey();
    TUNGraph::TNodeI nIt = copyG->GetNI(nId);
    int deg = nIt.GetDeg();
    //
    int curIntersection = 0;
    for (int i = 0; i < deg; i++) {
      int nbrId = nIt.GetNbrNId(i);
      if (CAND.IsKey(nbrId))
        curIntersection++;
    }
    //
    if (maxIntersection < curIntersection) {
      maxIntersection = curIntersection;
      id = nId;
    }
  }
  return id;
}

// ExcludeEdgeV is sorted
bool CommunityTree::excludeEdgesInNonMaxClique(TIntV &K,
                                               TIntPrV &ExcludeEdgeV) {
  TIntPrV CliqueEdgeV;
  getPairCombination(K, CliqueEdgeV);
  if (CliqueEdgeV.IntrsLen(ExcludeEdgeV) > 0) {
    return true;
  } else {
    return false;
  }
}

void CommunityTree::TTTExcludeEdges(PUNGraph &G, TIntV &K, TIntSet &Cand,
                                    TIntSet &Fini, TIntPrV &ExcludeEdgeV) {
  if (Cand.Len() == 0 && Fini.Len() == 0) {
    MCVs->Add(K);
    if (K.Len() > maxMCSize) {
      maxMCSize = K.Len();
    }
    return;
  }
  // union of Cand and Fini
  TIntSet CFUnion(TMath::Mx(Cand.Len(), Fini.Len()));
  getUnion(Cand, Fini, CFUnion);
  // Get u in CFUnion that maximaze Cand intersection with neighbors of u
  int u = maxNbrsInCandNodeId(CFUnion, Cand);
  // Get neighbors of u
  TIntSet uNbrs;
  getNbrs(G, u, uNbrs);
  // Get relative complement of uNbrs in Cand
  TIntSet Ext;
  getRelativeComplement(Cand, uNbrs, Ext);

  for (THashSetKeyI<TInt> it = Ext.BegI(); it < Ext.EndI(); it++) {
    int q = it.GetKey();
    TIntV Kq = K;
    Kq.Add(q);
    // only enumerate the cliques within the graph that do not contain any edge
    // from ExcludeEdgeV
    if (excludeEdgesInNonMaxClique(Kq, ExcludeEdgeV)) {
      Cand.DelKey(q);
      Fini.AddKey(q);
      continue;
    }
    TIntSet qNbrs, Candq, Finiq;
    getNbrs(G, q, qNbrs);
    getIntersection(Cand, qNbrs, Candq);
    getIntersection(Fini, qNbrs, Finiq);
    TTTExcludeEdges(G, Kq, Candq, Finiq, ExcludeEdgeV);
    Cand.DelKey(q);
    Fini.AddKey(q);
  }
}

int CommunityTree::getNodeIdWithMaxDeg(const TIntSet &Set) const {
  int id = -1;
  int maxDeg = -1;
  //
  for (THashSetKeyI<TInt> it = Set.BegI(); it < Set.EndI(); it++) {
    int nId = it.GetKey();
    int deg = copyG->GetNI(nId).GetDeg();
    if (maxDeg < deg) {
      maxDeg = deg;
      id = nId;
    }
  }
  return id;
}

void CommunityTree::expandClique(const TIntSet &SUBG, TIntSet &CAND) {
  if (SUBG.Len() == 0) {
    TIntV copyMCV = MCV;
    copyMCV.Sort();
    MCVs->Add(copyMCV);
    if (MCV.Len() > maxMCSize) {
      maxMCSize = MCV.Len();
    }
    return;
  }
  if (CAND.Len() == 0)
    return;
  // Get u that maximaze Cand intersection with neighbours of vertex u
  int u = maxNbrsInCandNodeId(SUBG, CAND);
  // Get neighbours of node u
  TIntSet nbrsU;
  getNbrs(copyG, u, nbrsU);
  // Get relative complement of nbrsU in CAND
  TIntSet EXT;
  getRelativeComplement(CAND, nbrsU, EXT);
  while (EXT.Len() != 0) {
    int q = getNodeIdWithMaxDeg(EXT);
    //
    MCV.Add(q);
    //
    THashSet<TInt> nbrsQ;
    getNbrs(copyG, q, nbrsQ);
    //
    THashSet<TInt> SUBGq;
    getIntersection(SUBG, nbrsQ, SUBGq);
    //
    THashSet<TInt> CANDq;
    getIntersection(CAND, nbrsQ, CANDq);
    //
    expandClique(SUBGq, CANDq);
    //
    CAND.DelKey(q);
    MCV.DelLast();
    //
    EXT.Clr();
    getRelativeComplement(CAND, nbrsU, EXT);
  }
}

void CommunityTree::getMaximalCliques(const PUNGraph &G, TVec<TIntV> &initMCs) {
  if (G->GetNodes() == 0)
    return;
  copyG = G;
  MCVs = &initMCs;
  MCV.Clr();

  TIntSet SUBG, CAND;
  for (TUNGraph::TNodeI NI = copyG->BegNI(); NI < copyG->EndNI(); NI++) {
    TInt nId = NI.GetId();
    SUBG.AddKey(nId);
    CAND.AddKey(nId);
  }
  expandClique(SUBG, CAND);
  initMCs = *MCVs;
}

void CommunityTree::getMCs(const PUNGraph &G, TVec<TIntV> &initMCs) {
  CommunityTree CT;
  initMCs.Clr(false);
  CT.getMaximalCliques(G, initMCs);
};

void CommunityTree::enumNewMCs(const PUNGraph &G, TIntPrV &NewEdgeV,
                               TVec<TIntV> &newMCs) {
  copyG = G;
  TIntPrV ExcludeEdgeV;
  for (int e = 0; e < NewEdgeV.Len(); e++) {
    if (!copyG->IsNode(NewEdgeV[e].Val1)) {
      copyG->AddNode(NewEdgeV[e].Val1);
    }
    if (!copyG->IsNode(NewEdgeV[e].Val2)) {
      copyG->AddNode(NewEdgeV[e].Val2);
    }
    copyG->AddEdge(NewEdgeV[e].Val1, NewEdgeV[e].Val2);
  }

  TIntV K, NbrV, SubGraphV;
  for (int e = 0; e < NewEdgeV.Len(); e++) {
    const int u = NewEdgeV[e].Val1;
    const int v = NewEdgeV[e].Val2;
    // find common neighbors between u and v and store the node ids in NbrV
    TSnap::GetCmnNbrs(copyG, u, v, NbrV);
    K = TIntV::GetV(u, v);
    NbrV.Union(K, SubGraphV);
    PUNGraph SubG = TSnap::ConvertSubGraph<PUNGraph>(copyG, SubGraphV);
    TIntSet Cand, Fini;
    Cand.AddKeyV(NbrV);
    TVec<TIntV> S;
    MCVs = &S;
    TTTExcludeEdges(SubG, K, Cand, Fini, ExcludeEdgeV);
    newMCs.Union(S);
    ExcludeEdgeV.Add(TIntPr(u, v));
  }
}

void CommunityTree::enumSubsumedMCsEdgeInsertion(TIntPrV &NewEdgeV,
                                                 TIntTrieH &TrieH,
                                                 TVec<TIntV> &newMCs,
                                                 TVec<TIntV> &delMCs) {

  for (int i = 0; i < newMCs.Len(); i++) {
    TVec<TIntV> S;
    enumPotentialSubsumedMCs(newMCs[i], NewEdgeV, S);
    for (int j = 0; j < S.Len(); j++) {
      TIntV C = S[j]; // C is sorted
      if (TrieH.IsKey(C[0])) {
        Trie *T = TrieH.GetDat(C[0]);
        if (T->getID(C) != -1) {
          delMCs.AddUnique(C);
        }
      }
    }
  }
}

bool CommunityTree::isMC(TIntV &C) {
  C.Sort();
  int v = C[0];
  TIntSet vLNbrs;
  getLNbrs(copyG, v, vLNbrs);
  vLNbrs.AddKey(v);
  for (THashSetKeyI<TInt> wIt = vLNbrs.BegI(); wIt < vLNbrs.EndI(); wIt++) {
    int w = wIt.GetKey();
    TIntSet wHNbrs;
    getHNbrs(copyG, w, wHNbrs);
    wHNbrs.AddKey(w);
    if (isKeyVIn(C, wHNbrs)) {
      TIntSet C_HS, comp;
      C_HS.AddKeyV(C);
      getRelativeComplement(wHNbrs, C_HS, comp);
      for (THashSetKeyI<TInt> uIt = comp.BegI(); uIt < comp.EndI(); uIt++) {
        TIntSet uNbrs;
        getNbrs(copyG, uIt.GetKey(), uNbrs);
        if (isKeyVIn(C, uNbrs)) {
          return false;
        }
      }
    }
  }
  return true;
}

void CommunityTree::enumSubsumedMCsEdgeDeletion(PUNGraph &G, TIntPrV &DelEdgeV,
                                                TVec<TIntV> &delMCs,
                                                TVec<TIntV> &newMCs) {
  copyG = G;
  for (int e = 0; e < DelEdgeV.Len(); e++) {
    copyG->DelEdge(DelEdgeV[e].Val1, DelEdgeV[e].Val2);
  }

  for (int i = 0; i < delMCs.Len(); i++) {
    TVec<TIntV> S;
    enumPotentialSubsumedMCs(delMCs[i], DelEdgeV, S);
    for (int j = 0; j < S.Len(); j++) {
      TIntV C = S[j];
      if (isMC(C)) {
        newMCs.AddUnique(C);
      }
    }
  }
}

CliqueGraph *CommunityTree::buildCT(PUNGraph &G, TIntTrieH &TrieH) {
  TVec<TIntV> initMCs;
  getMaximalCliques(G, initMCs);
  printf("number of init MCs: %d\n", initMCs.Len());
  // sort the size of the MCs in a descending order before added to CG
  initMCs.Sort(false);
  // generate a clique graph with levels determined by the max size of MCs
  CliqueGraph *CG = new CliqueGraph(maxMCSize);
  TVec<TIntH> prevCG;
  CG->addMCs(G, initMCs, TrieH, prevCG);

  return CG;
}

void CommunityTree::updateCTEI(PUNGraph &G, TIntPrV &newEdgeV, TIntTrieH &TrieH,
                               CliqueGraph *CG) {
  TVec<TIntV> newMCs, delMCs;
  TVec<TIntH> prevCG;
  prevCG.Gen(CG->forests.Len());
  for (int i = 0; i < CG->forests.Len(); i++) {
    for (THashKeyDatI<TInt, ETNode *> fIt = CG->forests[i].BegI();
         fIt < CG->forests[i].EndI(); fIt++) {
      prevCG[i].AddDat(fIt.GetKey(), fIt.GetDat()->cg1->size);
    }
  }

  TTmProfiler Profiler;
  int TimerId = Profiler.AddTimer("Profiler");

  Profiler.ResetTimer(TimerId);
  Profiler.StartTimer(TimerId);
  enumNewMCs(G, newEdgeV, newMCs);
  printf("(V, E) = (%d, %d)\n", G->GetNodes(), G->GetEdges());
  Profiler.StopTimer(TimerId);
  printf("number of new MCs: %d\n", newMCs.Len());
  printf("run time for enumerating new MCs: %f\n",
         Profiler.GetTimerSec(TimerId));

  if (maxMCSize > CG->forests.Len()) {
    ETForestVec extraLevels;
    extraLevels.Gen(maxMCSize - CG->forests.Len());
    CG->forests.AddV(extraLevels);
  }

  Profiler.ResetTimer(TimerId);
  Profiler.StartTimer(TimerId);
  enumSubsumedMCsEdgeInsertion(newEdgeV, TrieH, newMCs, delMCs);
  Profiler.StopTimer(TimerId);
  printf("number of deleted MCs: %d\n", delMCs.Len());
  printf("run time for enumerating deleted MCs: %f\n",
         Profiler.GetTimerSec(TimerId));

  newMCs.Sort(false);
  Profiler.ResetTimer(TimerId);
  Profiler.StartTimer(TimerId);
  CG->removeMCs(delMCs, TrieH);
  Profiler.StopTimer(TimerId);
  printf("run time for removing deleted MCs: %f\n",
         Profiler.GetTimerSec(TimerId));

  Profiler.ResetTimer(TimerId);
  Profiler.StartTimer(TimerId);
  CG->addMCs(G, newMCs, TrieH, prevCG);
  Profiler.StopTimer(TimerId);
  printf("run time for adding MCs: %f\n", Profiler.GetTimerSec(TimerId));
}

#endif /* COMMUNITYTREE_H_ */