//////////////////////////////////////////////////////////////////////
// 
// TPiEG2.h
//
//////////////////////////////////////////////////////////////////////

#ifndef __TPiEG2__
#define __TPiEG2__

#include "TClasTool.h"
#include "TPartSieve.h"
#include "TPartSieveHists.h"
#include "TInterrupt.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TTree.h"  
#include "TStopwatch.h"
//#include "haprad_constants.h"
//#include "TRadCor.h"
//#include "THapradConfig.h"

class TPiEG2: public TClasTool
{

 public:

  Bool_t gFiltered_data;

  TPartSieve *Sieve;
  TPartSieveHists *SHists;

  TInterrupt *Interrupt;
  TStopwatch *Time;
  TObjArray  *H; // Stores the histogram objects for operations on all histos.

//  TRadCor *rc;
//const THapradConfig* fC;
  // Pointers for the histograms. Also stored in the H object array.
  // This is for easy access to the histogram objects.
  
  TH1D *Elec_vert_z;
  TH2D *Elec_vert_xy;
  TH1F *SampFrac;
  TH1F *EMass;

  TH1D *PiPMom;
  TH1D *PiPCos;
  TH1D *PiMMom;
  TH1D *PiMCos;

  // Pointers for tree 
  TTree *treeZCut;  

 public:
  TPiEG2(); // Initialize code
  ~TPiEG2()
  {
    // Cleanup to avoid memory leakage.
    delete Time;
    DeleteHistos();
    delete SHists;
    delete Sieve;
  }
  Int_t Run(Int_t Runn=0, Int_t Nevt=2147483647);
  void  InitHistos(void);
  void  DeleteHistos(void);
  void  ClearHistos(void);
  void  Write(void);
  TVector3 vertex_corr(TVector3 v, TVector3 p);
  int sect(Double_t phi);

  ClassDef(TPiEG2,1); // A Class for checking the PID.
};

#endif
