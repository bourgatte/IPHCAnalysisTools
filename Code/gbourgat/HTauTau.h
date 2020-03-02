#ifndef HTauTau_h
#define HTauTau_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "boost/functional/hash.hpp"
//#include "SVFitStorage.h"
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
#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"

class HTauTau : public Selection {

 public:
  HTauTau(TString Name_, TString id_);
  virtual ~HTauTau();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {//Trigger=0,
	      Id_and_Kin=0, 
	     /* NPairsFound, */
	     /* Tau1Isolation, */
	     /* Tau2Isolation, */
	     /* LeptonVeto, */
	      PairCharge, 
	     /* PairMass, */
	     //MTM,
	     NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();
  ReferenceScaleFactors *RSF;
  int TriggerOkDummy, selVertexDummy, selMuon_IsoDummy, selMuon_AntiIsoDummy, selTauDummy, ChargeSumDummy;
  double MTDummy, MvisDummy, TauFLSigmaDummy;

  int Charge;

  PUReweight reweight;//(PUReweight::RUN2ANALYSIS);
  DataMCCorrections DataMC_Corr;
  tauTrigSFreader tauTrgSF;

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
  /*
  std::vector<TH1D> Tau1isolation;
  std::vector<TH1D> Tau2isolation;
  
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

  std::vector<TH1D> TauTauVisMass;
  std::vector<TH1D> TauTauTruthMass;
  std::vector<TH1D> TauTauFullMass;
 
  std::vector<TH1D> dRTauTau;
  std::vector<TH1D> QCDShape;

  std::vector<TH1D> NQCD;
  std::vector<TH1D> TauTauFullMass_B;
  std::vector<TH1D> TauTauFullMass_C;
  std::vector<TH1D> TauTauFullMass_D;
  
  std::vector<TH1D> MET;
  std::vector<TH1D> METphi;
  std::vector<TH1D> PUPPImet;
  std::vector<TH1D> PUPPImetphi;
  std::vector<TH1D> TransverseMass;
  
  std::vector<TH1D> NPrimeVtx;
  std::vector<TH1D> NPU;
  std::vector<TH1D> RHO;
  
  std::vector<TH1D> NbJets;
  
  std::vector<TH1D> h_SVFitMass; 
  std::vector<TH1D> h_SVFitStatus;
  /* std::vector<TH1D> svfTau1E; */
  /* std::vector<TH1D> svfTau2E; */
  
  /* std::vector<TH1D> PhiDatasvfitpipi; */
  /* std::vector<TH1D> PhiDatasvfitpirho; */
  /* std::vector<TH1D> PhiDatasvfitpia1; */
  /* std::vector<TH1D> PhiDatasvfitrhorho; */
  /* std::vector<TH1D> PhiDatasvfitrhoa1; */
  /* std::vector<TH1D> PhiDatasvfita1a1; */
  /*
  std::vector<TH1D> PhiDatavispipi;
  std::vector<TH1D> PhiDatavispirho;
  std::vector<TH1D> PhiDatavispia1;
  std::vector<TH1D> PhiDatavisrhorho;
  std::vector<TH1D> PhiDatavisrhoa1;
  std::vector<TH1D> PhiDatavisa1a1;
  */
  
  /* std::vector<TH1D> Etasvfit; */
  /* std::vector<TH1D> Phisvfitpipi; */
  /* std::vector<TH1D> Phisvfitpirho; */
  /* //std::vector<TH1D> Phisvfitlpi; */
  /* //std::vector<TH1D> Phisvfitlrho; */
  /* std::vector<TH1D> Phisvfitpia1; */
  /* std::vector<TH1D> Phisvfitrhorho; */
  /* std::vector<TH1D> Phisvfitrhoa1; */
  /* // std::vector<TH1D> Phisvfitla1; */
  /* std::vector<TH1D> Phisvfita1a1; */
  /* std::vector<TH1D> Thetasvfit; */
  
  std::vector<TH1D> Etavis;
  std::vector<TH1D> Phivispipi;
  std::vector<TH1D> Phivispirho;
  // std::vector<TH1D> Phivislpi;
  // std::vector<TH1D> Phivislrho;
  std::vector<TH1D> Phivispia1;
  std::vector<TH1D> Phivisrhorho;
  std::vector<TH1D> Phivisrhoa1;
  //std::vector<TH1D> Phivisla1;
  std::vector<TH1D> Phivisa1a1;
  std::vector<TH1D> Thetavis;
  
  std::vector<TH1D> Etatruth;
  std::vector<TH1D> Phitruthpipi;
  std::vector<TH1D> Phitruthpirho;
  // std::vector<TH1D> Phitruthlpi;
  // std::vector<TH1D> Phitruthlrho;
  std::vector<TH1D> Phitruthpia1;
  std::vector<TH1D> Phitruthrhorho;
  std::vector<TH1D> Phitruthrhoa1;
  // std::vector<TH1D> Phitruthla1;
  std::vector<TH1D> Phitrutha1a1;
  std::vector<TH1D> Thetatruth;
  /*
  std::vector<TH1D> PhiSvFitRespipi;
  std::vector<TH1D> PhiSvFitRespirho;
  std::vector<TH1D> PhiSvFitReslpi;
  std::vector<TH1D> PhiSvFitReslrho;
  std::vector<TH1D> PhiSvFitRespia1;
   std::vector<TH1D> PhiSvFitResrhorho;
  std::vector<TH1D> PhiSvFitResrhoa1;
  std::vector<TH1D> PhiSvFitResla1;
  std::vector<TH1D> PhiSvFitResa1a1;
  
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
  std::vector<TH2D> NewPhivsPhiproton;
  std::vector<TH2D> NewPhivsPhiTauplus;
  std::vector<TH2D> NewPhivsEtaproton;
  std::vector<TH2D> NewPhivsEtaTauplus;
  std::vector<TH2D> NewPhivsZPt;
  std::vector<TH1D> NewPhiSignal;
  std::vector<TH1D> NewPhiQCD;

  std::vector<TH1D> IstauminusvisPhysical;
  std::vector<TH1D> IstauplusvisPhysical;
  std::vector<TH1D> IsPairPhysical;
  
  //std::vector<TH1D> ResolPullTauTauFroma1a1MeanEnergy;
  //std::vector<TH1D> ResolPullTauTauFroma1a1MZEnergy;
  std::vector<TH1D> ResolPullTauTauFroma1a1MeanMomentum;
  std::vector<TH1D> ResolPullTauTauFroma1a1MZMomentum;
  
  //std::vector<TH1D> ResolPullTauminusFroma1a1MeanEnergy;
  //std::vector<TH1D> ResolPullTauminusFroma1a1MZEnergy;
  std::vector<TH1D> ResolPullTauminusFroma1a1MeanMomentum;
  std::vector<TH1D> ResolPullTauminusFroma1a1MZMomentum;

  //std::vector<TH1D> ResolPullTauplusFroma1a1MeanEnergy;
  //std::vector<TH1D> ResolPullTauplusFroma1a1MZEnergy;
  std::vector<TH1D> ResolPullTauplusFroma1a1MeanMomentum;
  std::vector<TH1D> ResolPullTauplusFroma1a1MZMomentum;

  //std::vector<TH1D> ResolPullTauFroma1a1MeanEnergy;
  //std::vector<TH1D> ResolPullTauFroma1a1MZEnergy;
  std::vector<TH1D> ResolPullTauFroma1a1MeanMomentum;
  std::vector<TH1D> ResolPullTauFroma1a1MZMomentum;

  //std::vector<TH1D> ResolPullTauminusFroma1XMeanEnergy;
  std::vector<TH1D> ResolPullTauminusFroma1XMeanMomentum;
  //std::vector<TH1D> ResolPullXminusFroma1XMeanEnergy;
  std::vector<TH1D> ResolPullXminusFroma1XMeanMomentum;

  //std::vector<TH1D> ResolPullTauplusFroma1XMeanEnergy;
  std::vector<TH1D> ResolPullTauplusFroma1XMeanMomentum; 
  //std::vector<TH1D> ResolPullXplusFroma1XMeanEnergy;
  std::vector<TH1D> ResolPullXplusFroma1XMeanMomentum;


  std::vector<TH1D> ResolPullXVtxIna1a1;
  std::vector<TH1D> ResolPullYVtxIna1a1;
  std::vector<TH1D> ResolPullZVtxIna1a1;


  std::vector<TH1D> tauminusa1a1MomentumVis;                        
  std::vector<TH1D> tauplusa1a1MomentumVis;
  std::vector<TH1D> InvariantMasstausa1a1Vis;
		    
  std::vector<TH1D> tauminusa1a1MomentumMean;
  std::vector<TH1D> tauplusa1a1MomentumMean;
  std::vector<TH1D> InvariantMasstausa1a1Mean;

  std::vector<TH1D> tauminusa1a1MomentumPairConstraint;
  std::vector<TH1D> tauplusa1a1MomentumPairConstraint;
  std::vector<TH1D> InvariantMasstausa1a1PairConstraint;
	
  std::vector<TH1D> tauminusa1XMomentumVis;
  std::vector<TH1D> tauplusa1XMomentumVis;
  std::vector<TH1D> InvariantMasstausa1XVis;
		    
  std::vector<TH1D> tauminusa1XMomentumMean;
  std::vector<TH1D> tauplusa1XMomentumMean;
  std::vector<TH1D> InvariantMasstausa1XMean;


  std::vector<TH1D> polarimetricAcopAngle;

  std::vector<TH1D> polarimetricAcopAnglePVRefitNoBSOld;
  std::vector<TH1D> polarimetricAcopAnglePVRefitBSOld;
  std::vector<TH1D> polarimetricAcopAnglePVRefitNoBSZNominalOld;
  std::vector<TH1D> polarimetricAcopAnglePVRefitBSZNominalOld;
   
  std::vector<TH1D> polarimetricAcopAnglePVRefitNoBSOneTrackRemovedOld;
  std::vector<TH1D> polarimetricAcopAnglePVRefitBSOneTrackRemovedOld;
  std::vector<TH1D> polarimetricAcopAnglePVRefitNoBSOneTrackRemovedZNominalOld;
  std::vector<TH1D> polarimetricAcopAnglePVRefitBSOneTrackRemovedZNominalOld;

  std::vector<TH1D> polarimetricAcopAnglePVRefitNoBSTracksRemovedOld;
  std::vector<TH1D> polarimetricAcopAnglePVRefitBSTracksRemovedOld;
  std::vector<TH1D> polarimetricAcopAnglePVRefitNoBSTracksRemovedZNominalOld;
  std::vector<TH1D> polarimetricAcopAnglePVRefitBSTracksRemovedZNominalOld;

  std::vector<TH1D> polarimetricAcopAnglePVRefitNoBSNew;
  std::vector<TH1D> polarimetricAcopAnglePVRefitBSNew;
  std::vector<TH1D> polarimetricAcopAnglePVRefitNoBSZNominalNew;
  std::vector<TH1D> polarimetricAcopAnglePVRefitBSZNominalNew;

  std::vector<TH1D> polarimetricAcopAnglePVRefitOnlyNoBSTracksRemovedOld;
  std::vector<TH1D> polarimetricAcopAnglePVRefitOnlyBSTracksRemovedOld;
  std::vector<TH1D> polarimetricAcopAnglePVRefitOnlyNoBSTracksRemovedZNominalOld;
  std::vector<TH1D> polarimetricAcopAnglePVRefitOnlyBSTracksRemovedZNominalOld;

  std::vector<TH1D> polarimetricAcopAnglePVRefitOnlyNoBSNew;
  std::vector<TH1D> polarimetricAcopAnglePVRefitOnlyBSNew;
  std::vector<TH1D> polarimetricAcopAnglePVRefitOnlyNoBSZNominalNew;
  std::vector<TH1D> polarimetricAcopAnglePVRefitOnlyBSZNominalNew;
   
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitNoBSTIP; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitBSTIP; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitNoBSZNominalTIP; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitBSZNominalTIP; */

  std::vector<TH1D> polarimetricAcopAngleMVA;

  std::vector<TH1D> polarimetricAcopAnglePVRefitNoBSOldMVA;
  std::vector<TH1D> polarimetricAcopAnglePVRefitBSOldMVA;
  std::vector<TH1D> polarimetricAcopAnglePVRefitNoBSZNominalOldMVA;
  std::vector<TH1D> polarimetricAcopAnglePVRefitBSZNominalOldMVA;
   
  std::vector<TH1D> polarimetricAcopAnglePVRefitNoBSOneTrackRemovedOldMVA;
  std::vector<TH1D> polarimetricAcopAnglePVRefitBSOneTrackRemovedOldMVA;
  std::vector<TH1D> polarimetricAcopAnglePVRefitNoBSOneTrackRemovedZNominalOldMVA;
  std::vector<TH1D> polarimetricAcopAnglePVRefitBSOneTrackRemovedZNominalOldMVA;

  std::vector<TH1D> polarimetricAcopAnglePVRefitNoBSTracksRemovedOldMVA;
  std::vector<TH1D> polarimetricAcopAnglePVRefitBSTracksRemovedOldMVA;
  std::vector<TH1D> polarimetricAcopAnglePVRefitNoBSTracksRemovedZNominalOldMVA;
  std::vector<TH1D> polarimetricAcopAnglePVRefitBSTracksRemovedZNominalOldMVA;

  std::vector<TH1D> polarimetricAcopAnglePVRefitNoBSNewMVA;
  std::vector<TH1D> polarimetricAcopAnglePVRefitBSNewMVA;
  std::vector<TH1D> polarimetricAcopAnglePVRefitNoBSZNominalNewMVA;
  std::vector<TH1D> polarimetricAcopAnglePVRefitBSZNominalNewMVA;

  std::vector<TH1D> polarimetricAcopAnglePVRefitOnlyNoBSTracksRemovedOldMVA;
  std::vector<TH1D> polarimetricAcopAnglePVRefitOnlyBSTracksRemovedOldMVA;
  std::vector<TH1D> polarimetricAcopAnglePVRefitOnlyNoBSTracksRemovedZNominalOldMVA;
  std::vector<TH1D> polarimetricAcopAnglePVRefitOnlyBSTracksRemovedZNominalOldMVA;

  std::vector<TH1D> polarimetricAcopAnglePVRefitOnlyNoBSNewMVA;
  std::vector<TH1D> polarimetricAcopAnglePVRefitOnlyBSNewMVA;
  std::vector<TH1D> polarimetricAcopAnglePVRefitOnlyNoBSZNominalNewMVA;
  std::vector<TH1D> polarimetricAcopAnglePVRefitOnlyBSZNominalNewMVA;
   
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitNoBSTIPMVA; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitBSTIPMVA; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitNoBSZNominalTIPMVA; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitBSZNominalTIPMVA; */
   
  std::vector<TH1D> polarimetricAcopAngleTruthA1;
  
  /* std::vector<TH1D> polarimetricAcopAngleTruthRho; */
  /* std::vector<TH1D> polarimetricAcopAngleTruthPi; */

  std::vector<TH1D> polarimetricAcopAngleDecayPlane;

  /* std::vector<TH1D> polarimetricAcopAnglePVRefitNoBS; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitBS; */
  /* std::vector<TH1D> polarimetricAcopAngleMVA; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitNoBSMVA; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitBSMVA; */
  
  /* std::vector<TH1D> polarimetricAcopAngle30; */
  /* std::vector<TH1D> polarimetricAcopAngle25; */
  /* std::vector<TH1D> polarimetricAcopAngle20; */
  /* std::vector<TH1D> polarimetricAcopAngle15; */
  /* std::vector<TH1D> polarimetricAcopAngle10; */
  /* std::vector<TH1D> polarimetricAcopAngle5; */
  
  /* std::vector<TH1D> AcolAngle; */
  /* std::vector<TH1D> AcolAngleSVFit; */
  /* std::vector<TH1D> AcolAngleTruth; */
  
  /* std::vector<TH1D> polarimetricAcopAnglePtTruthA1;   */
  /* std::vector<TH1D> polarimetricAcopAngleMVAPtTruthA1; */
 

  /* std::vector<TH1D> polarimetricAcopAnglePtTruthRho; */
  /* std::vector<TH1D> polarimetricAcopAngleMVAPtTruthRho; */

  /* std::vector<TH1D> polarimetricAcopAnglePtTruthPi; */
  /* std::vector<TH1D> polarimetricAcopAngleMVAPtTruthPi; */

   std::vector<TH1D> polarimetricAcopAngleSVFit; 
   std::vector<TH1D> polarimetricAcopAngleMVASVFit; 

  /* std::vector<TH1D> polarimetricAcopAngleSVFitRho; */
  /* std::vector<TH1D> polarimetricAcopAngleMVASVFitRho; */
  
  /* std::vector<TH1D> polarimetricAcopAngleSVFitPi; */
  /* std::vector<TH1D> polarimetricAcopAngleMVASVFitPi; */
  
  
  /* std::vector<TH1D> polarimetricAcopAngleBackground; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitNoBSBackground; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitBSBackground; */
  /* std::vector<TH1D> polarimetricAcopAngleMVABackground; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitNoBSMVABackground; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitBSMVABackground; */

  /* std::vector<TH1D> polarimetricAcopAnglePVRefitNoBSWithoutWSpin; */  

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
  std::vector<TH1D> PVXNoBSOldResol;
  std::vector<TH1D> PVXNoBSTracksRemovedOldResol;
  std::vector<TH1D> PVXNoBSNewResol;
  std::vector<TH1D> PVXBSOldResol;
  std::vector<TH1D> PVXBSTracksRemovedOldResol;
  std::vector<TH1D> PVXBSNewResol;

  std::vector<TH1D> PVYResol;
  std::vector<TH1D> PVYNoBSOldResol;
  std::vector<TH1D> PVYNoBSTracksRemovedOldResol;
  std::vector<TH1D> PVYNoBSNewResol;
  std::vector<TH1D> PVYBSOldResol;
  std::vector<TH1D> PVYBSTracksRemovedOldResol;
  std::vector<TH1D> PVYBSNewResol;
  
  std::vector<TH1D> PVZResol;
  std::vector<TH1D> PVZNoBSOldResol;
  std::vector<TH1D> PVZNoBSTracksRemovedOldResol;
  std::vector<TH1D> PVZNoBSNewResol;
  std::vector<TH1D> PVZBSOldResol;
  std::vector<TH1D> PVZBSTracksRemovedOldResol;
  std::vector<TH1D> PVZBSNewResol;
  
  std::vector<TH1D> PVXNoBSTracksRemovedOldOnlyResol;
  std::vector<TH1D> PVXNoBSNewOnlyResol;
  std::vector<TH1D> PVXBSTracksRemovedOldOnlyResol;
  std::vector<TH1D> PVXBSNewOnlyResol;
  
  std::vector<TH1D> PVYNoBSTracksRemovedOldOnlyResol;
  std::vector<TH1D> PVYNoBSNewOnlyResol;
  std::vector<TH1D> PVYBSTracksRemovedOldOnlyResol;
  std::vector<TH1D> PVYBSNewOnlyResol;
  
  std::vector<TH1D> PVZNoBSTracksRemovedOldOnlyResol;
  std::vector<TH1D> PVZNoBSNewOnlyResol;
  std::vector<TH1D> PVZBSTracksRemovedOldOnlyResol;
  std::vector<TH1D> PVZBSNewOnlyResol;
};
#endif
