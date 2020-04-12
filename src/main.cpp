//============================================================================
// Author      :
// Description :
//============================================================================

#include <iostream>
#include "Snap.h"
#include "CPM.h"

using namespace std;

int main(int argc, char *argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(
      TStr::Fmt("Network community detection. build: %s, %s. Time: % s ",
                __TIME__, __DATE__, TExeTm::GetCurTm()));
  TTmProfiler Profiler;

  int TimerId = Profiler.AddTimer("Profiler");
  int TimerAll = Profiler.AddTimer("ProfilerAll");

  const TStr InFNm = Env.GetIfArgPrefixStr("-i:", "../data/sms01.txt",
                                           "Input graph (undirected graph) ");

  // First load the temporal evolving network
  PTimeNENet FullNENet = TTimeNENet::LoadEdgeTm(InFNm, 1, 2, 3, ssfWhiteSep);
  TTimeNet::TTmBucketV TmBucketV;
  TTmUnit TmUnit = tmuWeek;
  FullNENet->GetEdgeTmBuckets(TmUnit, TmBucketV);
  int startTm = 7;
  int endTm = 9;
  // int endTm = TmBucketV.Len();

  TIntPrV PrevEdgeV, CurEdgeV, NewEdgeV;
  TIntV EdgeLstV;
  for (int t = 0; t <= startTm; t++) {
    printf("edgelist length: %d\n", EdgeLstV.Len());
    EdgeLstV.AddV(TmBucketV[t].NIdV);
  }
  // get the initial graph G_0
  PUNGraph PrevGraph = TSnap::ConvertESubGraph<PUNGraph>(FullNENet, EdgeLstV);
  // delete self-edges
  TSnap::DelSelfEdges(PrevGraph);
  TSnap::DelZeroDegNodes(PrevGraph);
  if (PrevGraph->GetEdges() == 0) {
    printf(
        "No edge exists in the network. Choose another starting time period.");
    return 0;
  }
  printf("***%d/%d: by %s (V,E) = (%d,%d)\n", startTm + 1, TmBucketV.Len(),
         TmBucketV[startTm].BegTm.GetStr(TmUnit).CStr(), PrevGraph->GetNodes(),
         PrevGraph->GetEdges());

  // reserve memory for edge values in the vector
  PrevEdgeV.Gen(PrevGraph->GetEdges(), 0);
  for (TUNGraph::TEdgeI EI = PrevGraph->BegEI(); EI < PrevGraph->EndEI();
       EI++) {
    PrevEdgeV.Add(TIntPr(EI.GetSrcNId(), EI.GetDstNId()));
  }

  // sort the vector
  PrevEdgeV.Sort();

  Profiler.ResetTimer(TimerAll);
  Profiler.StartTimer(TimerAll);

  // build the initial community tree
  TIntTrieH TrieH;
  CommunityTree *CT = new CommunityTree();
  CT->time = startTm;

  CliqueGraph *CG = CT->buildCT(PrevGraph, TrieH);
  TIntPrV *P = CG->generatePD();
  THash<TInt, TVec<TIntV>> cmtyH;
  CG->convertCommunities(TrieH, cmtyH);
  Profiler.StopTimer(TimerAll);
  printf("run time: %f\n", Profiler.GetTimerSec(TimerAll));
  cmtyH.Defrag();

  int level = 5;
  TVec<TIntV> cmtyVV = cmtyH.GetDat(level);
  saveTxt(cmtyVV, TStr::Fmt("t%du.txt", startTm + 1));
  for (int i = 0; i < CG->forests.Len(); i++) {
    printf("%d communities at level %d\n", CG->forests[i].Len(), i + 1);
    for (THashKeyDatI<TInt, ETNode *> fIt = CG->forests[i].BegI();
         fIt < CG->forests[i].EndI(); fIt++) {
      cout << fIt.GetKey() << " ";
    }
    cout << endl;
  }
  // for (TIntPrV::TIter it = P->BegI(); it < P->EndI(); it++) {
  //   printf("(b,d): (%d,%d)\n", it->Val1, it->Val2);
  // }
  delete P;
  PrevGraph->Defrag();

  int windowLen = 3;

  if (endTm - startTm - 1 > windowLen) {
    TSnapQueue<int> TSNUpBdQ;
    TSNUpBdQ.Gen(windowLen, 0);
    int TSNUpBdQSum = 0;

    for (int t = startTm + 1; t < endTm; t++) {
      printf("edgelist length: %d\n", EdgeLstV.Len());
      EdgeLstV.AddV(TmBucketV[t].NIdV);

      PUNGraph CurGraph =
          TSnap::ConvertESubGraph<PUNGraph>(FullNENet, EdgeLstV);
      TSnap::DelSelfEdges(CurGraph);
      TSnap::DelZeroDegNodes(CurGraph);
      printf("***%d/%d: by %s (V,E) = (%d,%d) ", t + 1, TmBucketV.Len(),
             TmBucketV[t].BegTm.GetStr(TmUnit).CStr(), CurGraph->GetNodes(),
             CurGraph->GetEdges());

      CurEdgeV.Gen(CurGraph->GetEdges(), 0);
      for (TUNGraph::TEdgeI EI = CurGraph->BegEI(); EI < CurGraph->EndEI();
           EI++) {
        CurEdgeV.Add(TIntPr(EI.GetSrcNId(), EI.GetDstNId()));
      }
      CurEdgeV.Sort();
      // new edges are stored in NewEdgeV
      CurEdgeV.Diff(PrevEdgeV, NewEdgeV);

      int TSNUpBd = TSNUpperBound(NewEdgeV);
      cout << "estimated upper bound for TSN: " << TSNUpBd << endl;

      if (TSNUpBdQ.Len() < windowLen) {
        TSNUpBdQ.Push(TSNUpBd);
        TSNUpBdQSum += TSNUpBd;

        cmtyVV.Clr();
        getCPMCommunities(CurGraph, level, cmtyVV);
        saveTxt(cmtyVV, TStr::Fmt("t%d.txt", t + 1));
      } else {
        // if the current estimated upper bound for TSN is no less the moving
        // mean of the previous window
        if (TSNUpBd >= TSNUpBdQSum / windowLen) {
          int preUpdateT = CT->time;
          cout << "previous update time: " << preUpdateT + 1 << endl;

          TIntV lastEdgeLstV;
          for (int i = 0; i <= preUpdateT; i++) {
            lastEdgeLstV.AddV(TmBucketV[i].NIdV);
          }

          PUNGraph LastGraph =
              TSnap::ConvertESubGraph<PUNGraph>(FullNENet, lastEdgeLstV);
          TSnap::DelSelfEdges(LastGraph);
          TSnap::DelZeroDegNodes(LastGraph);
          printf("last graph (V,E) = (%d,%d) ", LastGraph->GetNodes(),
                 LastGraph->GetEdges());

          Profiler.ResetTimer(TimerAll);
          Profiler.StartTimer(TimerAll);
          TIntPrV LastEdgeV, NewEdgeVFrLastG;
          LastEdgeV.Gen(LastGraph->GetEdges(), 0);
          for (TUNGraph::TEdgeI EI = LastGraph->BegEI();
               EI < LastGraph->EndEI(); EI++) {
            LastEdgeV.Add(TIntPr(EI.GetSrcNId(), EI.GetDstNId()));
          }
          LastEdgeV.Sort();
          CurEdgeV.Diff(LastEdgeV, NewEdgeVFrLastG);

          // update the community tree
          CT->updateCTEI(LastGraph, NewEdgeVFrLastG, TrieH, CG);
          CT->time = t;

          TIntPrV *P = CG->generatePD();

          Profiler.ResetTimer(TimerId);
          Profiler.StartTimer(TimerId);
          THash<TInt, TVec<TIntV>> cmtyH;
          CG->convertCommunities(TrieH, cmtyH);
          Profiler.StopTimer(TimerId);
          printf("run time for converting to communities: %f\n",
                 Profiler.GetTimerSec(TimerId));
          Profiler.StopTimer(TimerAll);
          printf("run time: %f\n", Profiler.GetTimerSec(TimerAll));
          cmtyH.Defrag();

          cmtyVV.Clr();
          cmtyVV = cmtyH.GetDat(level);
          saveTxt(cmtyVV, TStr::Fmt("t%du.txt", t + 1));
          for (int i = 0; i < CG->forests.Len(); i++) {
            printf("%d communities at level %d\n", CG->forests[i].Len(), i + 1);
            for (THashKeyDatI<TInt, ETNode *> fIt = CG->forests[i].BegI();
                 fIt < CG->forests[i].EndI(); fIt++) {
              cout << fIt.GetKey() << " ";
            }
            cout << endl;
          }

          delete P;
          LastGraph->Defrag();
        } else {
          cmtyVV.Clr();
          getCPMCommunities(CurGraph, level, cmtyVV);
          saveTxt(cmtyVV, TStr::Fmt("t%d.txt", t + 1));
        }
        // update the sum of the queue
        TSNUpBdQSum = TSNUpBdQSum - TSNUpBdQ.Top();
        TSNUpBdQ.Pop();
        TSNUpBdQ.Push(TSNUpBd);
        TSNUpBdQSum += TSNUpBd;
      }
      CurGraph->Defrag();
      PrevEdgeV = CurEdgeV;
    }
  }
}