// -*- C++ -*-
//
// 
/**\class TauPolInterface.h TauPolInterface.cc
 Description: 
*/
//
// Original Author:  Vladimir Cherepanov 
//         Created:  Sun Sep 10 17:49:02 CET 2017
//         Description of observables can be found in
//         https://arxiv.org/pdf/1805.10552.pdf and
//         http://inspirehep.net/record/318161
//



#ifndef TauPolInterface_h
#define TauPolInterface_h
#include <vector>
#include "TLorentzVector.h"
#include "TComplex.h"
#include "TMatrixT.h"
#include "TVectorT.h"
#include "TMatrixTSym.h"
#include <string.h>
#include "rhoHelper.h"
#include "a1Helper.h"

using namespace std;
class TauPolInterface {
 public:
  TauPolInterface();
  TauPolInterface(vector<TLorentzVector> TauAndProd, string type, int TauCharge=1);
  TauPolInterface(vector<TLorentzVector> TauAndProd1, string type1,  vector<TLorentzVector> TauAndProd2, string type2,int TauCharge1=1, int TauCharge2=1);
  ~TauPolInterface();


  void Configure(vector<TLorentzVector> TauAndProd, string type, int TauCharge = 1 );
  void ConfigurePair(vector<TLorentzVector> TauAndProd1, string type1,  vector<TLorentzVector> TauAndProd2, string type2, int TauCharge1 = 1,int TauCharge2 = 1);
  bool  isConfigured();
  bool  isPairConfigured();
  TLorentzVector Boost(TLorentzVector pB, TLorentzVector frame);
  
  double getVisibleOmega(string which="");
  double getOmega(string which="");
  double getOmegabar(string which="");
  double getCombOmega();
  double getCombOmegaBar();
  double getCombVisibleOmega();
  TLorentzVector getVisiblePairLV();
 private:

  void SetupLeg(string which="");


  bool debug;
  TMatrixT<double> convertToMatrix(TVectorT<double> V);

  TLorentzVector TauLV;
  TLorentzVector ProductLV;
  TLorentzVector TauRhoPi;
  TLorentzVector TauRhoPi0;

  TLorentzVector TauA1OSPi;
  TLorentzVector TauA1SSPi1;
  TLorentzVector TauA1SSPi2;

  TLorentzVector InvisibleLV;
  TLorentzVector DPF_TauLV;
  TLorentzVector DPF_TauRhoPi;
  TLorentzVector DPF_TauRhoPi0;
  TLorentzVector DPF_InvisibleLV;
  string type_;
  int taucharge_;

  TLorentzVector TauLV1,TauLV2;
  TLorentzVector ProductLV1,ProductLV2;

  TLorentzVector TauRhoPi1,TauRhoPi2;
  TLorentzVector TauRhoPi01,TauRhoPi02;
  TLorentzVector TauA1OSPi1,TauA1OSPi2;
  TLorentzVector TauA1SSPi11,TauA1SSPi12;
  TLorentzVector TauA1SSPi21,TauA1SSPi22;
  TLorentzVector InvisibleLV1,InvisibleLV2;
  TLorentzVector DPF_TauLV1,DPF_TauLV2;
  TLorentzVector DPF_TauRhoPi1,DPF_TauRhoPi2;
  TLorentzVector DPF_TauRhoPi01,DPF_TauRhoPi02;
  TLorentzVector DPF_InvisibleLV1,DPF_InvisibleLV2;
  string type1_,type2_;
  int taucharge1_, taucharge2_;
};
#endif
