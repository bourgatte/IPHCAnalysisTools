// -*- C++ -*-
//
// 
/**\class SCalculator.h SCalculator.cc
 Description: 
*/
//
// Original Author:  Vladimir Cherepanov 
//         Created:  Mon Sep 4 13:49:02 CET 2017
//
//

#ifndef SCalculator_h
#define SCalculator_h


#include <vector>
#include "TLorentzVector.h"
#include "TComplex.h"
#include "TMatrixT.h"
#include "TVectorT.h"
#include "TMatrixTSym.h"
#include <string.h>
#include <vector>
#include "TLorentzVector.h"
#include "PolarimetricA1.h"
using namespace std;


class SCalculator {
 
 public:
  SCalculator();
  SCalculator(string type);
  ~SCalculator();
  void Configure(vector<TLorentzVector> TauAndProd, TLorentzVector Frame, int charge=1);
  bool isConfigured();
  std::vector<TLorentzVector> getBoosted(){return TauAndProd_LF;}
  TLorentzVector Boost(TLorentzVector pB, TLorentzVector frame);
  TVector3 Rotate(TVector3 LVec, TVector3 Rot);
  TVector3 pv();
  void SortPions(std::vector<TLorentzVector >& pionsvec, std::vector<double>& charges);

  //====================


 private:
  string type_;
  int charge_;
  TVector3 polvec();
  vector<TLorentzVector> TauAndProd_LF;
  vector<TLorentzVector> TauAndProd_HRF;
  TLorentzVector TauLV;
  bool debug;
  TMatrixT<double> convertToMatrix(TVectorT<double> V);
};
#endif
