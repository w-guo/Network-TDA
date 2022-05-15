#ifndef CPM_H_
#define CPM_H_

#include "CommunityTree.h"

PUNGraph calculateOverlapMtxLevel(const TVec<TIntV> &MCs, int l) {
  const int n = MCs.Len();
  PUNGraph overlapMtxLevel = TUNGraph::New();
  for (int i = 0; i < n; i++) {
    overlapMtxLevel->AddNode(i);
  }
  // Calculate clique-clique overlap matrix
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      if (MCs[i].IntrsLen(MCs[j]) >= l - 1) {
        overlapMtxLevel->AddEdge(i, j);
      }
    }
  }
  return overlapMtxLevel;
}

void getCPMCommunities(const PUNGraph &G, int l, TVec<TIntV> &NIdCmtyVV) {
  // printf("\nClique Percolation Method\n");
  TVec<TIntV> allMCs;
  CommunityTree::getMCs(G, allMCs);
  allMCs.Sort(false);

  TVec<TIntV> MCs;
  for (int i = 0; i < allMCs.Len(); i++) {
    if (allMCs[i].Len() >= l) {
      MCs.Add(allMCs[i]);
    } else {
      break;
    }
  }

  // get the overlap matrix at level i
  PUNGraph overlapMtxLevel = calculateOverlapMtxLevel(MCs, l);
  // connected components are communities
  TCnComV CnComV;
  TSnap::GetWccs(overlapMtxLevel, CnComV);
  overlapMtxLevel->Defrag();

  NIdCmtyVV.Clr(false);
  TIntSet CmtySet;
  for (int c = 0; c < CnComV.Len(); c++) {
    CmtySet.Clr(false);
    for (int i = 0; i < CnComV[c].Len(); i++) {
      const TIntV &CliqueNIdV = MCs[CnComV[c][i]];
      CmtySet.AddKeyV(CliqueNIdV);
    }
    NIdCmtyVV.Add();
    CmtySet.GetKeyV(NIdCmtyVV.Last());
    NIdCmtyVV.Last().Sort();
  }
}

void getCPMAllLvlCommunities(const PUNGraph &G,
                             THash<TInt, TVec<TIntV>> &cmtyH) {
  printf("\nClique Percolation Method\n");
  TTmProfiler Profiler;
  int TimerId = Profiler.AddTimer("Profiler");
  int TimerAll = Profiler.AddTimer("ProfilerAll");

  Profiler.ResetTimer(TimerAll);
  Profiler.StartTimer(TimerAll);

  Profiler.ResetTimer(TimerId);
  Profiler.StartTimer(TimerId);

  TVec<TIntV> allMCs;
  CommunityTree::getMCs(G, allMCs);
  Profiler.StopTimer(TimerId);
  printf("run time for calculating all MCs: %f\n",
         Profiler.GetTimerSec(TimerId));
  allMCs.Sort(false);
  int maxMCSize = allMCs[0].Len();
  int minMCSize = 2;

  for (int l = maxMCSize; l >= minMCSize; l--) {
    printf("at level %d:\n", l);

    TVec<TIntV> MCs;
    if (l == minMCSize) {
      MCs = allMCs;
    } else {
      for (int i = 0; i < allMCs.Len(); i++) {
        if (allMCs[i].Len() >= l) {
          MCs.Add(allMCs[i]);
        } else {
          break;
        }
      }
    }

    Profiler.ResetTimer(TimerId);
    Profiler.StartTimer(TimerId);
    // get the overlap matrix at level i
    PUNGraph overlapMtxLevel = calculateOverlapMtxLevel(MCs, l);
    Profiler.StopTimer(TimerId);
    printf("run time for calculating overlapping matrix: %f\n",
           Profiler.GetTimerSec(TimerId));

    Profiler.ResetTimer(TimerId);
    Profiler.StartTimer(TimerId);
    // connected components are communities
    TCnComV CnComV;
    TSnap::GetWccs(overlapMtxLevel, CnComV);
    overlapMtxLevel->Defrag();

    TVec<TIntV> NIdCmtyVV;
    TIntSet CmtySet;
    for (int c = 0; c < CnComV.Len(); c++) {
      CmtySet.Clr(false);
      for (int i = 0; i < CnComV[c].Len(); i++) {
        const TIntV &CliqueNIdV = MCs[CnComV[c][i]];
        CmtySet.AddKeyV(CliqueNIdV);
      }
      NIdCmtyVV.Add();
      CmtySet.GetKeyV(NIdCmtyVV.Last());
      NIdCmtyVV.Last().Sort();
    }
    cmtyH.AddDat(l, NIdCmtyVV);
    Profiler.StopTimer(TimerId);
    printf("run time for converting CCs to communities: %f\n",
           Profiler.GetTimerSec(TimerId));
  }
  Profiler.StopTimer(TimerAll);
  printf("done [%f].\n", Profiler.GetTimerSec(TimerAll));
}

#endif /* CPM_H_ */