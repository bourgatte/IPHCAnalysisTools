#ifndef ZMuTau_h
#define ZMuTau_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "SVFitStorage.h"
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
#include "TVector3.h"
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
#include "PUReweight.h"
#include "tauTrigSFreader.h"
#include "DataMCCorrections.h"


class ZMuTau : public Selection {

 public:
  ZMuTau(TString Name_, TString id_);
  virtual ~ZMuTau();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {Trigger=0,
	     Id_and_Kin,
	     NPairsFound, 
	     MuonIsolation,
	     TauIsolation, 
	     DiMuon_Veto,
	     LeptonVeto,
	     PairCharge,
	     MTM,
	     // PairMass,
	     NCuts};
 
 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();
  ReferenceScaleFactors *RSF;

  int TriggerOkDummy, selVertexDummy, selMuon_IsoDummy, selMuon_AntiIsoDummy, selTauDummy, ChargeSumDummy;
  double MTDummy, MvisDummy, TauFLSigmaDummy;

  int Charge;

  double  cMu_pt,  cMu_eta,  cTau_pt,  cTau_eta;
  PUReweight reweight;//(PUReweight::RUN2ANALYSIS);
  //  tauTrigSFreader tauTrgSF;
  DataMCCorrections DataMC_Corr;
  DataMCCorrections DataMC_CorrLeptonIso;

 private:
  // Selection Variables and Histos

  std::vector<TH1D> TauPT;
  std::vector<TH1D> TauE;
  std::vector<TH1D> TauMass;
  std::vector<TH1D> TauPhi;
  std::vector<TH1D> TauEta;
  std::vector<TH1D> Taudz;

  std::vector<TH1D> MuonPT;
  std::vector<TH1D> MuonE;
  std::vector<TH1D> MuonMass;
  std::vector<TH1D> MuonPhi;
  std::vector<TH1D> MuonEta;
  std::vector<TH1D> Muondz;
  std::vector<TH1D> Muondxy;
  std::vector<TH1D> MuonIsol;
  /*
  std::vector<TH1D> againstElectronVLooseMVA6;
  std::vector<TH1D> againstElectronLooseMVA6;
  std::vector<TH1D> againstElectronMediumMVA6;
  std::vector<TH1D> againstElectronTightMVA6;
  std::vector<TH1D> againstElectronVTightMVA6;
  std::vector<TH1D> againstMuonLoose3;
  std::vector<TH1D> againstMuonTight3;
  std::vector<TH1D> byCombinedIsolationDeltaBetaCorrRaw3Hits;

  std::vector<TH1D> DiMuonVeto;
  std::vector<TH1D> ExtraLeptonVeto;
  */
  std::vector<TH1D> TauHPSDecayMode;
  std::vector<TH1D> TauTauVisMass;
  std::vector<TH1D> TauTauTruthMass;

  std::vector<TH1D> dRTauTau;
  std::vector<TH1D> QCDShape;

  std::vector<TH1D> NQCD;
  std::vector<TH1D> NWJets;
  std::vector<TH1D> NWJetsRelaxed;
  std::vector<TH1D> MTAHigh;
  std::vector<TH1D> MTBHigh;
  std::vector<TH1D> MTBLow;
  std::vector<TH1D> TauTauVisMassAHigh;
  std::vector<TH1D> TauTauVisMassBHigh;
  std::vector<TH1D> TauTauVisMassBLow;

  std::vector<TH1D> METAHigh;

  std::vector<TH1D> MET;
  std::vector<TH1D> METphi;


  std::vector<TH1D> PUPPImet;
  std::vector<TH1D> PUPPImetphi;
  std::vector<TH1D> TransverseMass;
  
  std::vector<TH1D> NPrimeVtx;
  std::vector<TH1D> NPU;
  std::vector<TH1D> RHO;
  
  std::vector<TH1D> SVChi2;
  std::vector<TH1D> SVQuality;
  std::vector<TH2D> SVQualityVsSignificance;
  std::vector<TH1D> PVSVSignificance;
  std::vector<TH1D> NbJets;

 
  std::vector<TH1D> Etavis;
  // std::vector<TH1D> Phivispipi;
  //std::vector<TH1D> Phivispirho;
  std::vector<TH1D> Phivislpi;
  std::vector<TH1D> Phivislrho;
  //std::vector<TH1D> Phivispia1;
  //std::vector<TH1D> Phivisrhorho;
  //std::vector<TH1D> Phivisrhoa1;
  std::vector<TH1D> Phivisla1;
  // std::vector<TH1D> Phivisa1a1;
  std::vector<TH1D> Thetavis;
  
  //std::vector<TH1D> PhiInf8GeVvis;
  //std::vector<TH1D> PhiSup8GeVvis;

  std::vector<TH1D> Etatruth;
  //std::vector<TH1D> Phitruthpipi;
  // std::vector<TH1D> Phitruthpirho;
   std::vector<TH1D> Phitruthlpi;
   std::vector<TH1D> Phitruthlrho;
   // std::vector<TH1D> Phitruthpia1;
   //std::vector<TH1D> Phitruthrhorho;
   //std::vector<TH1D> Phitruthrhoa1;
  std::vector<TH1D> Phitruthla1;
  //std::vector<TH1D> Phitrutha1a1;
  std::vector<TH1D> Thetatruth;
  /* 
  std::vector<TH1D> PhiVisRespipi;
  std::vector<TH1D> PhiVisRespirho;
  std::vector<TH1D> PhiVisReslpi;
  std::vector<TH1D> PhiVisReslrho;
  std::vector<TH1D> PhiVisRespia1;
  std::vector<TH1D> PhiVisResrhorho;
  std::vector<TH1D> PhiVisResrhoa1;
  std::vector<TH1D> PhiVisResla1;
  std::vector<TH1D> PhiVisResa1a1;  
  */

  /*
  std::vector<TH1D> TauTauFullPtRes;
  std::vector<TH1D> TauTauFullEtaRes;
  std::vector<TH1D> TauTauFullPhiRes;

  std::vector<TH1D> TauplusFullPtRes;
  std::vector<TH1D> TauplusFullEtaRes;
  std::vector<TH1D> TauplusFullPhiRes;

  std::vector<TH1D> TauminusFullPtRes;
  std::vector<TH1D> TauminusFullEtaRes;
  std::vector<TH1D> TauminusFullPhiRes;
  
  std::vector<TH1D> TauTauVisPtRes;
  std::vector<TH1D> TauTauVisEtaRes;
  std::vector<TH1D> TauTauVisPhiRes;

  std::vector<TH1D> TauplusVisPtRes;
  std::vector<TH1D> TauplusVisEtaRes;
  std::vector<TH1D> TauplusVisPhiRes;
    
  std::vector<TH1D> TauminusVisPtRes;
  std::vector<TH1D> TauminusVisEtaRes;
  std::vector<TH1D> TauminusVisPhiRes;

  std::vector<TH1D> TauTauFullPtResPull;
  std::vector<TH1D> TauTauFullEtaResPull;
  std::vector<TH1D> TauTauFullPhiResPull;

  std::vector<TH1D> TauplusFullPtResPull;
  std::vector<TH1D> TauplusFullEtaResPull;
  std::vector<TH1D> TauplusFullPhiResPull;

  std::vector<TH1D> TauminusFullPtResPull;
  std::vector<TH1D> TauminusFullEtaResPull;
  std::vector<TH1D> TauminusFullPhiResPull;

  std::vector<TH1D> TauTauVisPtResPull;
  std::vector<TH1D> TauTauVisEtaResPull;
  std::vector<TH1D> TauTauVisPhiResPull;

  std::vector<TH1D> TauplusVisPtResPull;
  std::vector<TH1D> TauplusVisEtaResPull;
  std::vector<TH1D> TauplusVisPhiResPull;
    
  std::vector<TH1D> TauminusVisPtResPull;
  std::vector<TH1D> TauminusVisEtaResPull;
  std::vector<TH1D> TauminusVisPhiResPull;
*/
  //std::vector<TH1D> DRTruth;
  //std::vector<TH1D> DRFull;
  // std::vector<TH1D> DRFullTruth;
  //std::vector<TH1D> DRVisTruth;
  
  std::vector<TH1D> Pi0EnergyRes;
  std::vector<TH1D> Pi0EnergyResPull;

  std::vector<TH1D> ZPtVis;

  std::vector<TH2D> NewPhivsDeltaPhi;
  std::vector<TH2D> NewPhivsDeltaEta;
  //std::vector<TH2D> NewPhivsPhiproton;
  std::vector<TH2D> NewPhivsPhiTauplus;
  //std::vector<TH2D> NewPhivsEtaproton;
  std::vector<TH2D> NewPhivsEtaTauplus;
  std::vector<TH2D> NewPhivsZPt;
  std::vector<TH1D> NewPhiSignal;
  std::vector<TH1D> NewPhiQCD;

  std::vector<TH2D> NewPhivsPhiTauminus;
  std::vector<TH2D> NewPhivsEtaTauminus;
  //std::vector<TH2D> NewPhivsTauminusPt;
  //std::vector<TH2D> NewPhivsTauplusPt;
  std::vector<TH2D> NewPhivsDeltaEtaplusminus;

  std::vector<TH2D> EtaplusvsEtaminus;
  std::vector<TH2D> EtaTau1vsEtaTau2;

  std::vector<TH2D> PhiTau1vsPhiTau2;
  
  std::vector<TH2D> DzTau1vsDzTau2;
  
};
#endif
