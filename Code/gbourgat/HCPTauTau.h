#ifndef HCPTauTau_h
#define HCPTauTau_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "boost/functional/hash.hpp"
//#include "SVFitStorage.h"
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
#include "TVector3.h"
//#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"
#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/MultiProngTauSolver.h"
#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"
#include "SimpleFits/FitSoftware/interface/TauA1NuConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/DiTauConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/GlobalEventFit.h"
#include "ReferenceScaleFactors.h"
#include "ScaleFactor.h"
#include "Objects.h"
//#include "PUReweight.h"
#include "tauTrigSFreader.h"
#include "DataMCCorrections.h"
#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"
#include "RecoilCorrector.h"

#include "RooWorkspace.h"
#include "RooFunctor.h"
#include <memory>

class HCPTauTau : public Selection {

 public:
  HCPTauTau(TString Name_, TString id_);
  virtual ~HCPTauTau();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {//Trigger=0,
	     //Id_and_Kin, 
	     /* NPairsFound, */
             METFiltersAndBTagVeto=0,
	     TausIsolation, 
	     //Tau2Isolation,
	     AgainstEleMu,
	     LeptonVeto,
	     PairCharge, 
	     //PairMass,
	     //MTM,
	     NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();
  ReferenceScaleFactors *RSF;
  int TriggerOkDummy, selVertexDummy, selMuon_IsoDummy, selMuon_AntiIsoDummy, selTauDummy, ChargeSumDummy;
  double MTDummy, MvisDummy, TauFLSigmaDummy;

  int Charge;

  //PUReweight reweight;//(PUReweight::RUN2ANALYSIS);
  //cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
  
  //DataMCCorrections DataMC_Corr;
  //tauTrigSFreader tauTrgSF;

  TString InputNtuplePath=Ntp->GetInputNtuplePath();
  bool isDY1050=(InputNtuplePath.Contains("10to50"));


  TFile *WorkSpaceFF2016;
  RooWorkspace *wFF2016;
  
  ClassicSVfit svfitAlgo1;
  //ClassicSVfit svfitAlgo2;
  //  SVFitStorage svfitstorage;

  
 private:
  // Selection Variables and Histos
  
  std::vector<TH1D> Tau1PT;
  std::vector<TH1D> Tau1E;
  std::vector<TH1D> Tau1Mass;
  std::vector<TH1D> Tau1Phi;
  std::vector<TH1D> Tau1Eta;
  std::vector<TH1D> Tau1dz;

  std::vector<TH1D> Tau2PT;
  std::vector<TH1D> Tau2E;
  std::vector<TH1D> Tau2Mass;
  std::vector<TH1D> Tau2Phi;
  std::vector<TH1D> Tau2Eta;
  std::vector<TH1D> Tau2dz;
  
  /* std::vector<TH1D> Tau1isolation; */
  /* std::vector<TH1D> Tau2isolation; */
  /*
  std::vector<TH1D> againstElectronVLooseMVA6_Tau1;
  std::vector<TH1D> againstElectronLooseMVA6_Tau1;
  std::vector<TH1D> againstElectronMediumMVA6_Tau1;
  std::vector<TH1D> againstElectronTightMVA6_Tau1;
  std::vector<TH1D> againstElectronVTightMVA6_Tau1;
  std::vector<TH1D> againstMuonLoose3_Tau1;
  std::vector<TH1D> againstMuonTight3_Tau1;
  std::vector<TH1D> byCombinedIsolationDeltaBetaCorrRaw3Hits_Tau1;

  std::vector<TH1D> againstElectronVLooseMVA6_Tau2;
  std::vector<TH1D> againstElectronLooseMVA6_Tau2;
  std::vector<TH1D> againstElectronMediumMVA6_Tau2;
  std::vector<TH1D> againstElectronTightMVA6_Tau2;
  std::vector<TH1D> againstElectronVTightMVA6_Tau2;
  std::vector<TH1D> againstMuonLoose3_Tau2;
  std::vector<TH1D> againstMuonTight3_Tau2;
  std::vector<TH1D> byCombinedIsolationDeltaBetaCorrRaw3Hits_Tau2;
  */
  std::vector<TH1D> ExtraLeptonVeto;
  std::vector<TH1D> Tau2HPSDecayMode;
  std::vector<TH1D> Tau1HPSDecayMode;
  std::vector<TH1D> Tau2MVADecayMode;
  std::vector<TH1D> Tau1MVADecayMode;

  std::vector<TH1D> TauTauVisMass;
  //std::vector<TH1D> TauTauTruthMass;
  std::vector<TH1D> TauTauFullMass;

  std::vector<TH1D> TauTauVisPT;
  std::vector<TH1D> Mjj;
 
  std::vector<TH1D> dRTauTau;
  std::vector<TH1D> QCDShape;

  std::vector<TH1D> NQCD;
  std::vector<TH1D> TauTauFullMass_B;
  std::vector<TH1D> TauTauFullMass_C;
  std::vector<TH1D> TauTauFullMass_D;
  
  std::vector<TH1D> NFFData;
  std::vector<TH1D> NFFLeadMC;
  std::vector<TH1D> W_res;
  
  std::vector<TH1D> MET;
  std::vector<TH1D> METphi;
  std::vector<TH1D> PUPPImet;
  std::vector<TH1D> PUPPImetphi;
  std::vector<TH1D> PUPPImetcorr;
  std::vector<TH1D> PUPPImetcorrphi;
  std::vector<TH1D> TransverseMass;
  
  std::vector<TH1D> NPrimeVtx;
  std::vector<TH1D> NPU;
  std::vector<TH1D> RHO;
  
  std::vector<TH1D> NbJets;
 
  std::vector<TH1D> ZPtVis;

  std::vector<TH1D> IstauminusvisPhysical;
  std::vector<TH1D> IstauplusvisPhysical;
  std::vector<TH1D> IsPairPhysical;

  std::vector<TH1D> ResolPullTauTauFroma1a1MZMomentum;
  std::vector<TH1D> ResolPullTauminusFroma1a1MZMomentum;
  std::vector<TH1D> ResolPullTauplusFroma1a1MZMomentum;
  std::vector<TH1D> ResolPullTauFroma1a1MZMomentum;

  std::vector<TH1D> ResolPullXVtxIna1a1;
  std::vector<TH1D> ResolPullYVtxIna1a1;
  std::vector<TH1D> ResolPullZVtxIna1a1;


  std::vector<TH1D> tauminusa1a1MomentumVis;                        
  std::vector<TH1D> tauplusa1a1MomentumVis;
  std::vector<TH1D> InvariantMasstausa1a1Vis;

  std::vector<TH1D> tauminusa1a1MomentumPairConstraint;
  std::vector<TH1D> tauplusa1a1MomentumPairConstraint;
  std::vector<TH1D> InvariantMasstausa1a1PairConstraint;

  std::vector<TH1D> polarimetricAcopAngle;

  std::vector<TH1D> polarimetricAcopAnglePVRefitNoBS;
  std::vector<TH1D> polarimetricAcopAnglePVRefitBS;
  std::vector<TH1D> polarimetricAcopAnglePVRefitNoBSZNominal;
  std::vector<TH1D> polarimetricAcopAnglePVRefitBSZNominal;

  std::vector<TH1D> polarimetricAcopAnglePVRefitOnlyNoBS;
  std::vector<TH1D> polarimetricAcopAnglePVRefitOnlyBS;
  std::vector<TH1D> polarimetricAcopAnglePVRefitOnlyNoBSZNominal;
  std::vector<TH1D> polarimetricAcopAnglePVRefitOnlyBSZNominal;

  std::vector<TH1D> polarimetricAcopAngleMVA;

  std::vector<TH1D> polarimetricAcopAnglePVRefitNoBSNewMVA;
  std::vector<TH1D> polarimetricAcopAnglePVRefitBSNewMVA;
  std::vector<TH1D> polarimetricAcopAnglePVRefitNoBSZNominalNewMVA;
  std::vector<TH1D> polarimetricAcopAnglePVRefitBSZNominalNewMVA;

  std::vector<TH1D> polarimetricAcopAnglePVRefitOnlyNoBSNewMVA;
  std::vector<TH1D> polarimetricAcopAnglePVRefitOnlyBSNewMVA;
  std::vector<TH1D> polarimetricAcopAnglePVRefitOnlyNoBSZNominalNewMVA;
  std::vector<TH1D> polarimetricAcopAnglePVRefitOnlyBSZNominalNewMVA;
   
  std::vector<TH1D> polarimetricAcopAngleTruthA1;

  std::vector<TH1D> polarimetricAcopAngleDecayPlane;

  std::vector<TH1D> polarimetricAcopAngleSVFit; 
  std::vector<TH1D> polarimetricAcopAngleMVASVFit; 

  std::vector<TH1D> test;  
  
  std::vector<TH1D> PurityDM;
  std::vector<TH1D> PurityNewMVA;

  std::vector<TH1D> TauSVFitPxResPull;
  std::vector<TH1D> TauSVFitPyResPull;
  std::vector<TH1D> TauSVFitPzResPull;
 

  std::vector<TH1D> TauPxResPull;
  std::vector<TH1D> TauPyResPull;
  std::vector<TH1D> TauPzResPull;

  std::vector<TH1D> TauSVFitPxResPullMVA;
  std::vector<TH1D> TauSVFitPyResPullMVA;
  std::vector<TH1D> TauSVFitPzResPullMVA;
 

  std::vector<TH1D> TauPxResPullMVA;
  std::vector<TH1D> TauPyResPullMVA;
  std::vector<TH1D> TauPzResPullMVA;

  std::vector<TH1D> PVXResol;

  std::vector<TH1D> PVXNoBSResol;
  std::vector<TH1D> PVXBSResol;

  std::vector<TH1D> PVYResol;
  std::vector<TH1D> PVYNoBSResol;
  std::vector<TH1D> PVYBSResol;
  
  std::vector<TH1D> PVZResol;
  std::vector<TH1D> PVZNoBSResol;
  std::vector<TH1D> PVZBSResol;
  
  std::vector<TH1D> PVXNoBSOnlyResol;
  std::vector<TH1D> PVXBSOnlyResol;
  
  std::vector<TH1D> PVYNoBSOnlyResol;
  std::vector<TH1D> PVYBSOnlyResol;
  
  std::vector<TH1D> PVZNoBSOnlyResol;
  std::vector<TH1D> PVZBSOnlyResol;

  
  //---------------------------------------------------------------------
  //Histos MC for QCD
  
  std::vector<TH1D> PUPPImetcorrQCDMC;
  
  std::vector<TH1D> Tau1PTQCDMC;
  std::vector<TH1D> Tau1EQCDMC;
  std::vector<TH1D> Tau1MassQCDMC;
  std::vector<TH1D> Tau1PhiQCDMC;
  std::vector<TH1D> Tau1EtaQCDMC;
  std::vector<TH1D> Tau1dzQCDMC;
  std::vector<TH1D> Tau1HPSDecayModeQCDMC;
  std::vector<TH1D> Tau1MVADecayModeQCDMC;
  
  std::vector<TH1D> Tau2PTQCDMC;
  std::vector<TH1D> Tau2EQCDMC;
  std::vector<TH1D> Tau2MassQCDMC;
  std::vector<TH1D> Tau2PhiQCDMC;
  std::vector<TH1D> Tau2EtaQCDMC;
  std::vector<TH1D> Tau2dzQCDMC;
  std::vector<TH1D> Tau2HPSDecayModeQCDMC;
  std::vector<TH1D> Tau2MVADecayModeQCDMC;

  std::vector<TH1D> NbJetsQCDMC;
  std::vector<TH1D> TauTauVisMassQCDMC;

  std::vector<TH1D> TauTauVisPTQCDMC;
  std::vector<TH1D> MjjQCDMC;
};

#endif
