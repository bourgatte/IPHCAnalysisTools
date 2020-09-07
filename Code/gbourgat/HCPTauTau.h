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
#include "MEtSys.h"
#include "BDTClassification.h"

#include "RooWorkspace.h"
#include "RooFunctor.h"
#include <memory>

class HCPTauTau : public Selection {

 public:
  HCPTauTau(TString Name_, TString id_);
  virtual ~HCPTauTau();

  virtual void  Configure();
  virtual void  Finish();


  enum cuts {Trigger=0,
	     //Id_and_Kin, 
	     /* NPairsFound, */
             //ZTTMC=0,
             //METFilters,
	     Id_and_Kin,
	     //genmatch,
	     TausIsolation, 
	     //Tau2Isolation,
	     AgainstEleMu,
	     LeptonVeto,
	     PairCharge, 
	     PairMass,
	     //MTM,
	     NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();
  ReferenceScaleFactors *RSF;
  int TriggerOkDummy, selVertexDummy, selMuon_IsoDummy, selMuon_AntiIsoDummy, selTauDummy, ChargeSumDummy;
  double MTDummy, MvisDummy, TauFLSigmaDummy;

  int boolEvents,boola1a1, boolSV, boolThreePions, boolacop0, boolacop1;

  int Charge;

  //PUReweight reweight;//(PUReweight::RUN2ANALYSIS);
  //cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
  
  //DataMCCorrections DataMC_Corr;
  //tauTrigSFreader tauTrgSF;

  TString InputNtuplePath=Ntp->GetInputNtuplePath();
  bool isDY1050=(InputNtuplePath.Contains("10to50"));

  TFile *WorkSpaceFF2016;
  RooWorkspace *wFF2016;
  TFile *WorkSpaceFF2017;
  RooWorkspace *wFF2017;
  TFile *WorkSpaceFF2018;
  RooWorkspace *wFF2018;
  
  //Int_t year;
  BDTClassification *BDT;
  
  //ClassicSVfit svfitAlgo1;
  //ClassicSVfit svfitAlgo2;
  //  SVFitStorage svfitstorage;

  
 private:
  // Selection Variables and Histos
  
  /* std::vector<TH1D> Tau1PT; */
  /* std::vector<TH1D> Tau1E; */
  /* std::vector<TH1D> Tau1Mass; */
  /* std::vector<TH1D> Tau1Phi; */
  /* std::vector<TH1D> Tau1Eta; */
  /* std::vector<TH1D> Tau1dz; */

  /* std::vector<TH1D> Tau2PT; */
  /* std::vector<TH1D> Tau2E; */
  /* std::vector<TH1D> Tau2Mass; */
  /* std::vector<TH1D> Tau2Phi; */
  /* std::vector<TH1D> Tau2Eta; */
  /* std::vector<TH1D> Tau2dz; */
  
  std::vector<TH1D> Tau1isolation;
  std::vector<TH1D> Tau2isolation;
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
  /* std::vector<TH1D> ExtraLeptonVeto; */
  /* std::vector<TH1D> Tau2HPSDecayMode; */
  /* std::vector<TH1D> Tau1HPSDecayMode; */
  /* std::vector<TH1D> Tau2MVADecayMode; */
  /* std::vector<TH1D> Tau1MVADecayMode; */
  /* std::vector<TH1D> Tau2GenMatch; */
  /* std::vector<TH1D> Tau1GenMatch; */
  /* std::vector<TH1D> TauTauVisMass; */
  /* std::vector<TH1D> TauTauFullMass; */
  /* std::vector<TH1D> TauTauVisPT; */
  /* std::vector<TH1D> TauTauFullPT; */
  /* std::vector<TH1D> Mjj; */
  /* std::vector<TH1D> dijetpt; */
  /* std::vector<TH1D> dijetphi; */
  /* std::vector<TH1D> jdeta; */
  /* std::vector<TH1D> jdphi; */
  /* std::vector<TH1D> jpt_2; */
  /* std::vector<TH1D> jeta_2; */
  /* std::vector<TH1D> jphi_2; */
  /* std::vector<TH1D> jpt_1; */
  /* std::vector<TH1D> jeta_1; */
  /* std::vector<TH1D> jphi_1; */
 
  /* std::vector<TH1D> dRTauTau; */
  /* std::vector<TH1D> QCDShape; */

  /* std::vector<TH1D> NQCD; */
  /* std::vector<TH1D> TauTauFullMass_B; */
  /* std::vector<TH1D> TauTauFullMass_C; */
  /* std::vector<TH1D> TauTauFullMass_D; */
  
  /* std::vector<TH1D> NFFData; */
  /* std::vector<TH1D> NFFLeadMC; */
  std::vector<TH2D> WfakesHiggs;
  std::vector<TH2D> WfakesJetFakes;
  std::vector<TH2D> WfakesZTT;

  std::vector<TH2D> WfakesHiggs_DP;
  std::vector<TH2D> WfakesJetFakes_DP;
  std::vector<TH2D> WfakesZTT_DP;

  
  /* std::vector<TH1D> MET; */
  /* std::vector<TH1D> METphi; */
  /* std::vector<TH1D> PUPPImet; */
  /* std::vector<TH1D> PUPPImetphi; */
  /* std::vector<TH1D> PUPPImetcorr; */
  /* std::vector<TH1D> PUPPImetcorrphi; */
  /* std::vector<TH1D> TransverseMass; */
  
  /* std::vector<TH1D> NPrimeVtx; */
  /* std::vector<TH1D> NPU; */
  /* std::vector<TH1D> RHO; */
  
  /* std::vector<TH1D> NbJets; */
 
  /* std::vector<TH1D> HPtVis; */

  /* std::vector<TH1D> IstauminusvisPhysical; */
  /* std::vector<TH1D> IstauplusvisPhysical; */
  /* std::vector<TH1D> IsPairPhysical; */

  /* std::vector<TH1D> ResolPullTauTauFroma1a1MZMomentum; */
  /* std::vector<TH1D> ResolPullTauminusFroma1a1MZMomentum; */
  /* std::vector<TH1D> ResolPullTauplusFroma1a1MZMomentum; */
  /* std::vector<TH1D> ResolPullTauFroma1a1MZMomentum; */

  /* std::vector<TH1D> ResolPullXVtxIna1a1; */
  /* std::vector<TH1D> ResolPullYVtxIna1a1; */
  /* std::vector<TH1D> ResolPullZVtxIna1a1; */


  /* std::vector<TH1D> tauminusa1a1MomentumVis;                         */
  /* std::vector<TH1D> tauplusa1a1MomentumVis; */
  /* std::vector<TH1D> InvariantMasstausa1a1Vis; */

  /* std::vector<TH1D> tauminusa1a1MomentumPairConstraint; */
  /* std::vector<TH1D> tauplusa1a1MomentumPairConstraint; */
  /* std::vector<TH1D> InvariantMasstausa1a1PairConstraint; */

  /* std::vector<TH1D> polarimetricAcopAngle; */

  //std::vector<TH1D> polarimetricAcopAnglePVRefitNoBS;
  //std::vector<TH1D> polarimetricAcopAnglePVRefitBS;
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitNoBSZNominal; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitBSZNominal; */

  //std::vector<TH1D> polarimetricAcopAngleMVADM;

  /* std::vector<TH1D> polarimetricAcopAnglePVRefitNoBSMVADM; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitBSMVADM; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitNoBSZNominalMVADM; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitBSZNominalMVADM; */


  //std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSMVADM;
  // std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADM;
   
  /* std::vector<TH1D> polarimetricAcopAngleMVADMHiggs; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitBSMVADMHiggs; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitBSZNominalMVADMHiggs; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggs; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMHiggs; */

  /* std::vector<TH1D> polarimetricAcopAngleMVADMJetFakes; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitBSMVADMJetFakes; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitBSZNominalMVADMJetFakes; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakes; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMJetFakes; */

  /* std::vector<TH1D> polarimetricAcopAngleMVADMZTTEmbed; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitBSMVADMZTTEmbed; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitBSZNominalMVADMZTTEmbed; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTEmbed; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMZTTEmbed; */

  std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSMVADM;
  std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSMVADMQCDMC;

  std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled;

  std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled;

  std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled;
  std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSMVADM_DP;
 std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSMVADMQCDMC_DP;
 
 std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled_DP;
 std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled_DP;
 std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled_DP;
  /* std::vector<TH2D> polarimetricAcopAngleMVADMHiggs; */
  /* std::vector<TH2D> polarimetricAcopAnglePVRefitBSMVADMHiggs; */
  /* std::vector<TH2D> polarimetricAcopAnglePVRefitBSZNominalMVADMHiggs; */
  std::vector<TH2D> polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggs;
  std::vector<TH2D> polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggs_DP;
  //  std::vector<TH2D> polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMHiggs;

  /* std::vector<TH2D> polarimetricAcopAngleMVADMJetFakes; */
  /* std::vector<TH2D> polarimetricAcopAnglePVRefitBSMVADMJetFakes; */
  /* std::vector<TH2D> polarimetricAcopAnglePVRefitBSZNominalMVADMJetFakes; */
  std::vector<TH2D> polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakes;
  std::vector<TH2D> polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakes_DP;
  //  std::vector<TH2D> polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMJetFakes;
  
  /* std::vector<TH2D> polarimetricAcopAngleMVADMZTT; */
  /* std::vector<TH2D> polarimetricAcopAnglePVRefitBSMVADMZTT; */
  /* std::vector<TH2D> polarimetricAcopAnglePVRefitBSZNominalMVADMZTT; */
  std::vector<TH2D> polarimetricAcopAnglePVRefitWithTracksBSMVADMZTT;
  std::vector<TH2D> polarimetricAcopAnglePVRefitWithTracksBSMVADMZTT_DP;
  //std::vector<TH2D> polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMZTT;
  
  
  /* std::vector<TH1D> polarimetricAcopAngleTruthA1; */

  /* std::vector<TH1D> polarimetricAcopAngleDecayPlane; */

  /* std::vector<TH1D> polarimetricAcopAngleSVFit;  */
  /* std::vector<TH1D> polarimetricAcopAngleMVASVFit;  */

  /* std::vector<TH1D> test;   */
  
  /* std::vector<TH1D> PurityDM; */
  /* std::vector<TH1D> PurityNewMVA; */

  std::vector<TH1D> TauSVFitPxResPull;
  std::vector<TH1D> TauSVFitPyResPull;
  std::vector<TH1D> TauSVFitPzResPull;
 
  std::vector<TH1D> TauSVFitPhiResPull;
  std::vector<TH1D> TauSVFitEtaResPull;

  std::vector<TH1D> TauPxResPull;
  std::vector<TH1D> TauPyResPull;
  std::vector<TH1D> TauPzResPull;
 
  std::vector<TH1D> TauPhiResPull;
  std::vector<TH1D> TauEtaResPull;

  std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSSVFitMVADM;

  /* std::vector<TH1D> TauPxResPull; */
  /* std::vector<TH1D> TauPyResPull; */
  /* std::vector<TH1D> TauPzResPull; */

  /* std::vector<TH1D> TauSVFitPxResPullMVA; */
  /* std::vector<TH1D> TauSVFitPyResPullMVA; */
  /* std::vector<TH1D> TauSVFitPzResPullMVA; */
 

  /* std::vector<TH1D> TauPxResPullMVA; */
  /* std::vector<TH1D> TauPyResPullMVA; */
  /* std::vector<TH1D> TauPzResPullMVA; */

  /* std::vector<TH1D> PVXResol; */

  /* std::vector<TH1D> PVXNoBSResol; */
  /* std::vector<TH1D> PVXBSResol; */

  /* std::vector<TH1D> PVYResol; */
  /* std::vector<TH1D> PVYNoBSResol; */
  /* std::vector<TH1D> PVYBSResol; */
  
  /* std::vector<TH1D> PVZResol; */
  /* std::vector<TH1D> PVZNoBSResol; */
  /* std::vector<TH1D> PVZBSResol; */
  
  /* std::vector<TH1D> PVXNoBSOnlyResol; */
  /* std::vector<TH1D> PVXBSOnlyResol; */
  
  /* std::vector<TH1D> PVYNoBSOnlyResol; */
  /* std::vector<TH1D> PVYBSOnlyResol; */
  
  /* std::vector<TH1D> PVZNoBSOnlyResol; */
  /* std::vector<TH1D> PVZBSOnlyResol; */

  /* std::vector<TH1D> HiggsBDTScore; */
  /* std::vector<TH1D> JetFakesBDTScore; */
  /* std::vector<TH1D> ZTTBDTScore; */

  /* //--------------------------------------------------------------------- */
  /* //Histos MC for a1-a1 */
  
  std::vector<TH1D> Tau1PTa1a1;
  /* std::vector<TH1D> Tau1Ea1a1; */
  /* std::vector<TH1D> Tau1Massa1a1; */
  /* std::vector<TH1D> Tau1Phia1a1; */
  /* std::vector<TH1D> Tau1Etaa1a1; */
  /* std::vector<TH1D> Tau1dza1a1; */

  /* std::vector<TH1D> Tau2PTa1a1; */
  /* std::vector<TH1D> Tau2Ea1a1; */
  /* std::vector<TH1D> Tau2Massa1a1; */
  /* std::vector<TH1D> Tau2Phia1a1; */
  /* std::vector<TH1D> Tau2Etaa1a1; */
  /* std::vector<TH1D> Tau2dza1a1; */
  
  std::vector<TH1D> TauTauVisMassa1a1;
  std::vector<TH1D> TauTauFullMassa1a1;
  std::vector<TH1D> TauTauVisPTa1a1;
  std::vector<TH1D> TauTauFullPTa1a1;
  std::vector<TH1D> Mjja1a1;
  /* std::vector<TH1D> dijetpta1a1; */
  /* std::vector<TH1D> dijetphia1a1; */
  std::vector<TH1D> jdetaa1a1; 
  /* std::vector<TH1D> jdphia1a1; */
  /* std::vector<TH1D> jpt_2a1a1; */
  /* std::vector<TH1D> jeta_2a1a1; */
  /* std::vector<TH1D> jphi_2a1a1; */
  std::vector<TH1D> jpt_1a1a1; 
  /* std::vector<TH1D> jeta_1a1a1; */
  /* std::vector<TH1D> jphi_1a1a1; */
  
  std::vector<TH1D> PUPPImetcorra1a1; 
  /* std::vector<TH1D> PUPPImetcorrphia1a1; */
  
  std::vector<TH1D> NbJetsa1a1; 
  
  std::vector<TH1D> HiggsBDTScorea1a1;
  std::vector<TH1D> JetFakesBDTScorea1a1;
  std::vector<TH1D> ZTTBDTScorea1a1;
  
  std::vector<TH1D> HiggsBDTScorea1a1_DP;
  std::vector<TH1D> JetFakesBDTScorea1a1_DP;
  std::vector<TH1D> ZTTBDTScorea1a1_DP;

   std::vector<TH1D> Tau1PTa1a1QCDMC; 
  /* std::vector<TH1D> Tau1Ea1a1QCDMC; */
  /* std::vector<TH1D> Tau1Massa1a1QCDMC; */
  /* std::vector<TH1D> Tau1Phia1a1QCDMC; */
  /* std::vector<TH1D> Tau1Etaa1a1QCDMC; */
  /* std::vector<TH1D> Tau1dza1a1QCDMC; */

  /* std::vector<TH1D> Tau2PTa1a1QCDMC; */
  /* std::vector<TH1D> Tau2Ea1a1QCDMC; */
  /* std::vector<TH1D> Tau2Massa1a1QCDMC; */
  /* std::vector<TH1D> Tau2Phia1a1QCDMC; */
  /* std::vector<TH1D> Tau2Etaa1a1QCDMC; */
  /* std::vector<TH1D> Tau2dza1a1QCDMC; */
  
  /* std::vector<TH1D> Tau2GenMatchQCDMC; */
  /* std::vector<TH1D> Tau1GenMatchQCDMC; */
  
  std::vector<TH1D> TauTauVisMassa1a1QCDMC;
  std::vector<TH1D> TauTauFullMassa1a1QCDMC;
  std::vector<TH1D> TauTauVisPTa1a1QCDMC;
  std::vector<TH1D> TauTauFullPTa1a1QCDMC;
  std::vector<TH1D> Mjja1a1QCDMC; 
  /* std::vector<TH1D> dijetpta1a1QCDMC; */
  /* std::vector<TH1D> dijetphia1a1QCDMC; */
   std::vector<TH1D> jdetaa1a1QCDMC; 
  /* std::vector<TH1D> jdphia1a1QCDMC; */
  /* std::vector<TH1D> jpt_2a1a1QCDMC; */
  /* std::vector<TH1D> jeta_2a1a1QCDMC; */
  /* std::vector<TH1D> jphi_2a1a1QCDMC; */
   std::vector<TH1D> jpt_1a1a1QCDMC; 
  /* std::vector<TH1D> jeta_1a1a1QCDMC; */
  /* std::vector<TH1D> jphi_1a1a1QCDMC; */
  
   std::vector<TH1D> PUPPImetcorra1a1QCDMC; 
  /* std::vector<TH1D> PUPPImetcorrphia1a1QCDMC; */

   std::vector<TH1D> NbJetsa1a1QCDMC; 

  std::vector<TH1D> HiggsBDTScorea1a1QCDMC;
  std::vector<TH1D> JetFakesBDTScorea1a1QCDMC;
  std::vector<TH1D> ZTTBDTScorea1a1QCDMC;
  std::vector<TH1D> HiggsBDTScorea1a1QCDMC_DP;
  std::vector<TH1D> JetFakesBDTScorea1a1QCDMC_DP;
  std::vector<TH1D> ZTTBDTScorea1a1QCDMC_DP;

  /* //--------------------------------------------------------------------- */
  /* //Histos MC for QCD */
  
  /* std::vector<TH1D> PUPPImetcorrQCDMC; */
  /* std::vector<TH1D> PUPPImetcorrphiQCDMC; */

  /* std::vector<TH1D> Tau1PTQCDMC; */
  /* std::vector<TH1D> Tau1EQCDMC; */
  /* std::vector<TH1D> Tau1MassQCDMC; */
  /* std::vector<TH1D> Tau1PhiQCDMC; */
  /* std::vector<TH1D> Tau1EtaQCDMC; */
  /* std::vector<TH1D> Tau1dzQCDMC; */
  /* std::vector<TH1D> Tau1HPSDecayModeQCDMC; */
  /* std::vector<TH1D> Tau1MVADecayModeQCDMC; */
  
  /* std::vector<TH1D> Tau2PTQCDMC; */
  /* std::vector<TH1D> Tau2EQCDMC; */
  /* std::vector<TH1D> Tau2MassQCDMC; */
  /* std::vector<TH1D> Tau2PhiQCDMC; */
  /* std::vector<TH1D> Tau2EtaQCDMC; */
  /* std::vector<TH1D> Tau2dzQCDMC; */
  /* std::vector<TH1D> Tau2HPSDecayModeQCDMC; */
  /* std::vector<TH1D> Tau2MVADecayModeQCDMC; */

  /* std::vector<TH1D> NbJetsQCDMC; */
  /* std::vector<TH1D> TauTauVisMassQCDMC; */
  /* std::vector<TH1D> TauTauFullMassQCDMC; */
  /* std::vector<TH1D> TauTauFullPTQCDMC; */
  /* std::vector<TH1D> TauTauVisPTQCDMC; */

  /* std::vector<TH1D> MjjQCDMC; */
  /* std::vector<TH1D> dijetptQCDMC; */
  /* std::vector<TH1D> dijetphiQCDMC; */
  /* std::vector<TH1D> jdetaQCDMC; */
  /* std::vector<TH1D> jdphiQCDMC; */
  /* std::vector<TH1D> jpt_2QCDMC; */
  /* std::vector<TH1D> jeta_2QCDMC; */
  /* std::vector<TH1D> jphi_2QCDMC; */
  /* std::vector<TH1D> jpt_1QCDMC; */
  /* std::vector<TH1D> jeta_1QCDMC; */
  /* std::vector<TH1D> jphi_1QCDMC; */

  /* std::vector<TH1D> HiggsBDTScoreQCDMC; */
  /* std::vector<TH1D> JetFakesBDTScoreQCDMC; */
  /* std::vector<TH1D> ZTTBDTScoreQCDMC; */

  /* std::vector<TH1D> polarimetricAcopAngleMVADMQCDMC; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitBSMVADMQCDMC; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitBSZNominalMVADMQCDMC; */

  //std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMQCDMC;

  /* std::vector<TH1D> polarimetricAcopAngleMVADMHiggsQCDMC; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitBSMVADMHiggsQCDMC; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitBSZNominalMVADMHiggsQCDMC; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsQCDMC; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMHiggsQCDMC; */

  /* std::vector<TH1D> polarimetricAcopAngleMVADMJetFakesQCDMC; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitBSMVADMJetFakesQCDMC; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitBSZNominalMVADMJetFakesQCDMC; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesQCDMC; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMJetFakesQCDMC; */

  /* std::vector<TH1D> polarimetricAcopAngleMVADMZTTEmbedQCDMC; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitBSMVADMZTTEmbedQCDMC; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitBSZNominalMVADMZTTEmbedQCDMC; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTEmbedQCDMC; */
  /* std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMZTTEmbedQCDMC; */

  std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolledQCDMC;

  std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolledQCDMC;

  std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolledQCDMC;
std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolledQCDMC_DP;
  std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolledQCDMC_DP;
  std::vector<TH1D> polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolledQCDMC_DP;
  /* std::vector<TH2D> polarimetricAcopAngleMVADMHiggsQCDMC; */
  /* std::vector<TH2D> polarimetricAcopAnglePVRefitBSMVADMHiggsQCDMC; */
  /* std::vector<TH2D> polarimetricAcopAnglePVRefitBSZNominalMVADMHiggsQCDMC; */
  std::vector<TH2D> polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsQCDMC;
  std::vector<TH2D> polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsQCDMC_DP;
  //  std::vector<TH2D> polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMHiggsQCDMC;

  /* std::vector<TH2D> polarimetricAcopAngleMVADMJetFakesQCDMC; */
  /* std::vector<TH2D> polarimetricAcopAnglePVRefitBSMVADMJetFakesQCDMC; */
  /* std::vector<TH2D> polarimetricAcopAnglePVRefitBSZNominalMVADMJetFakesQCDMC; */
  std::vector<TH2D> polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesQCDMC;
  std::vector<TH2D> polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesQCDMC_DP;
  //  std::vector<TH2D> polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMJetFakesQCDMC;
  
  /* std::vector<TH2D> polarimetricAcopAngleMVADMZTTQCDMC; */
  /* std::vector<TH2D> polarimetricAcopAnglePVRefitBSMVADMZTTQCDMC; */
  /* std::vector<TH2D> polarimetricAcopAnglePVRefitBSZNominalMVADMZTTQCDMC; */
  std::vector<TH2D> polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTQCDMC;
  std::vector<TH2D> polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTQCDMC_DP;
  //  std::vector<TH2D> polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMZTTQCDMC;
  
  
  ////////////////// Shape Uncertainties

  std::vector<TH2D> ShapeSystPVRefitWithTracksBSHiggs;
  std::vector<TH2D> ShapeSystPVRefitWithTracksBSWfakesHiggs;
  std::vector<TH2D> ShapeSystPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> ShapeSystPVRefitWithTracksBSWfakesJetFakes;
  std::vector<TH2D> ShapeSystPVRefitWithTracksBSZTT;
  std::vector<TH2D> ShapeSystPVRefitWithTracksBSWfakesZTT;
  
  std::vector<TH2D> ShapeSystPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> ShapeSystPVRefitWithTracksBSWfakesHiggs_DP;
  std::vector<TH2D> ShapeSystPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> ShapeSystPVRefitWithTracksBSWfakesJetFakes_DP;
  std::vector<TH2D> ShapeSystPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> ShapeSystPVRefitWithTracksBSWfakesZTT_DP;

  std::vector<TH2D> CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSHiggs;
  std::vector<TH2D> CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSHiggs;
  std::vector<TH2D> CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSHiggs;
  std::vector<TH2D> CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSHiggs;
  std::vector<TH2D> CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSHiggs;
  std::vector<TH2D> CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSHiggs;
  /* std::vector<TH2D> CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> CMS_res_j_13TeVUpPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> CMS_res_j_13TeVDownPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSHiggs; */
  std::vector<TH2D> CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSHiggs;
  std::vector<TH2D> CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSHiggs;
  std::vector<TH2D> CMS_scale_gg_13TeVUpPVRefitWithTracksBSHiggs;
  std::vector<TH2D> CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSHiggs;
  std::vector<TH2D> CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSHiggs;
  std::vector<TH2D> CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSHiggs;
  std::vector<TH2D> CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSHiggs;
  std::vector<TH2D> CMS_scale_gg_13TeVDownPVRefitWithTracksBSHiggs;
  std::vector<TH2D> CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSHiggs;
  std::vector<TH2D> CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSHiggs;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggs;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggs;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggs;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggs;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggs;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggs;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggs;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggs;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggs;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggs;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggs;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggs;
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSHiggs; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSHiggs; */

  std::vector<TH2D> PrefiringUpPVRefitWithTracksBSHiggs;
  std::vector<TH2D> PrefiringDownPVRefitWithTracksBSHiggs;
  
  std::vector<TH2D> polarimetricAcopAnglePVRefitWithTracksBSMVADMWfakesHiggs;
  std::vector<TH2D> CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSWfakesHiggs;
  std::vector<TH2D> CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSWfakesHiggs;
  std::vector<TH2D> CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSWfakesHiggs;
  std::vector<TH2D> CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSWfakesHiggs;
  std::vector<TH2D> CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSWfakesHiggs;
  std::vector<TH2D> CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSWfakesHiggs;
  /* std::vector<TH2D> CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSWfakesHiggs; */
  /* std::vector<TH2D> CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSWfakesHiggs; */
  /* std::vector<TH2D> CMS_res_j_13TeVUpPVRefitWithTracksBSWfakesHiggs; */
  /* std::vector<TH2D> CMS_res_j_13TeVDownPVRefitWithTracksBSWfakesHiggs; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSWfakesHiggs; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSWfakesHiggs; */
  /* std::vector<TH2D> CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSWfakesHiggs; */
  /* std::vector<TH2D> CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSWfakesHiggs; */
  /* std::vector<TH2D> CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSWfakesHiggs; */
  /* std::vector<TH2D> CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSWfakesHiggs; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs; */
  /* std::vector<TH2D> CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs; */
  /* std::vector<TH2D> CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs; */
  /* std::vector<TH2D> CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSWfakesHiggs; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSWfakesHiggs; */
  /* std::vector<TH2D> CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSWfakesHiggs; */
  /* std::vector<TH2D> CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSWfakesHiggs; */
  /* std::vector<TH2D> CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSWfakesHiggs; */
  /* std::vector<TH2D> CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSWfakesHiggs; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs; */
  /* std::vector<TH2D> CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs; */
  /* std::vector<TH2D> CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs; */
  /* std::vector<TH2D> CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs; */
  /* std::vector<TH2D> CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSWfakesHiggs; */
  /* std::vector<TH2D> CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSWfakesHiggs; */
  /* std::vector<TH2D> CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSWfakesHiggs; */
  /* std::vector<TH2D> CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSWfakesHiggs; */
  /* std::vector<TH2D> CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSWfakesHiggs; */
  /* std::vector<TH2D> CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSWfakesHiggs; */
  std::vector<TH2D> CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesHiggs;
  std::vector<TH2D> CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSWfakesHiggs;
  std::vector<TH2D> CMS_scale_gg_13TeVUpPVRefitWithTracksBSWfakesHiggs;
  std::vector<TH2D> CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSWfakesHiggs;
  std::vector<TH2D> CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSWfakesHiggs;
  std::vector<TH2D> CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesHiggs;
  std::vector<TH2D> CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSWfakesHiggs;
  std::vector<TH2D> CMS_scale_gg_13TeVDownPVRefitWithTracksBSWfakesHiggs;
  std::vector<TH2D> CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSWfakesHiggs;
  std::vector<TH2D> CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSWfakesHiggs;
  std::vector<TH2D> PrefiringUpPVRefitWithTracksBSWfakesHiggs;
  std::vector<TH2D> PrefiringDownPVRefitWithTracksBSWfakesHiggs;


  std::vector<TH2D> CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSJetFakes;
  /* std::vector<TH2D> CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> CMS_res_j_13TeVUpPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> CMS_res_j_13TeVDownPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSJetFakes; */
  std::vector<TH2D> CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> CMS_scale_gg_13TeVUpPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> CMS_scale_gg_13TeVDownPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakes;
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSJetFakes; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSJetFakes; */
  
  std::vector<TH2D> PrefiringUpPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> PrefiringDownPVRefitWithTracksBSJetFakes;

  std::vector<TH2D> polarimetricAcopAnglePVRefitWithTracksBSMVADMWfakesJetFakes;
  std::vector<TH2D> CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSWfakesJetFakes;
  std::vector<TH2D> CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSWfakesJetFakes;
  std::vector<TH2D> CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSWfakesJetFakes;
  std::vector<TH2D> CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSWfakesJetFakes;
  std::vector<TH2D> CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSWfakesJetFakes;
  std::vector<TH2D> CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSWfakesJetFakes;
  /* std::vector<TH2D> CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSWfakesJetFakes; */
  /* std::vector<TH2D> CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSWfakesJetFakes; */
  /* std::vector<TH2D> CMS_res_j_13TeVUpPVRefitWithTracksBSWfakesJetFakes; */
  /* std::vector<TH2D> CMS_res_j_13TeVDownPVRefitWithTracksBSWfakesJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSWfakesJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSWfakesJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSWfakesJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSWfakesJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSWfakesJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSWfakesJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSWfakesJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSWfakesJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSWfakesJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSWfakesJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSWfakesJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSWfakesJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes; */
  /* std::vector<TH2D> CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes; */
  /* std::vector<TH2D> CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSWfakesJetFakes; */
  /* std::vector<TH2D> CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSWfakesJetFakes; */
  /* std::vector<TH2D> CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSWfakesJetFakes; */
  /* std::vector<TH2D> CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSWfakesJetFakes; */
  /* std::vector<TH2D> CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSWfakesJetFakes; */
  /* std::vector<TH2D> CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSWfakesJetFakes; */
  std::vector<TH2D> CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesJetFakes;
  std::vector<TH2D> CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSWfakesJetFakes;
  std::vector<TH2D> CMS_scale_gg_13TeVUpPVRefitWithTracksBSWfakesJetFakes;
  std::vector<TH2D> CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSWfakesJetFakes;
  std::vector<TH2D> CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSWfakesJetFakes;
  std::vector<TH2D> CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesJetFakes;
  std::vector<TH2D> CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSWfakesJetFakes;
  std::vector<TH2D> CMS_scale_gg_13TeVDownPVRefitWithTracksBSWfakesJetFakes;
  std::vector<TH2D> CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSWfakesJetFakes;
  std::vector<TH2D> CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSWfakesJetFakes;
  std::vector<TH2D> PrefiringUpPVRefitWithTracksBSWfakesJetFakes;
  std::vector<TH2D> PrefiringDownPVRefitWithTracksBSWfakesJetFakes;
  
  

  std::vector<TH2D> CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSZTT;
  std::vector<TH2D> CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSZTT;
  std::vector<TH2D> CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSZTT;
  std::vector<TH2D> CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSZTT;
  std::vector<TH2D> CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSZTT;
  std::vector<TH2D> CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSZTT;
  /* std::vector<TH2D> CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> CMS_res_j_13TeVUpPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> CMS_res_j_13TeVDownPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSZTT; */
  std::vector<TH2D> CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSZTT;
  std::vector<TH2D> CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSZTT;
  std::vector<TH2D> CMS_scale_gg_13TeVUpPVRefitWithTracksBSZTT;
  std::vector<TH2D> CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSZTT;
  std::vector<TH2D> CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSZTT;
  std::vector<TH2D> CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSZTT;
  std::vector<TH2D> CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSZTT;
  std::vector<TH2D> CMS_scale_gg_13TeVDownPVRefitWithTracksBSZTT;
  std::vector<TH2D> CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSZTT;
  std::vector<TH2D> CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSZTT;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTT;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTT;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTT;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTT;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTT;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTT;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTT;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTT;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTT;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTT;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTT;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTT;
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSZTT; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSZTT; */

  std::vector<TH2D> PrefiringUpPVRefitWithTracksBSZTT;
  std::vector<TH2D> PrefiringDownPVRefitWithTracksBSZTT;
  
  std::vector<TH2D> polarimetricAcopAnglePVRefitWithTracksBSMVADMWfakesZTT;
  std::vector<TH2D> CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSWfakesZTT;
  std::vector<TH2D> CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSWfakesZTT;
  std::vector<TH2D> CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSWfakesZTT;
  std::vector<TH2D> CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSWfakesZTT;
  std::vector<TH2D> CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSWfakesZTT;
  std::vector<TH2D> CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSWfakesZTT;
  /* std::vector<TH2D> CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSWfakesZTT; */
  /* std::vector<TH2D> CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSWfakesZTT; */
  /* std::vector<TH2D> CMS_res_j_13TeVUpPVRefitWithTracksBSWfakesZTT; */
  /* std::vector<TH2D> CMS_res_j_13TeVDownPVRefitWithTracksBSWfakesZTT; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSWfakesZTT; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSWfakesZTT; */
  /* std::vector<TH2D> CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSWfakesZTT; */
  /* std::vector<TH2D> CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSWfakesZTT; */
  /* std::vector<TH2D> CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSWfakesZTT; */
  /* std::vector<TH2D> CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSWfakesZTT; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSWfakesZTT; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSWfakesZTT; */
  /* std::vector<TH2D> CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSWfakesZTT; */
  /* std::vector<TH2D> CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSWfakesZTT; */
  /* std::vector<TH2D> CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSWfakesZTT; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSWfakesZTT; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSWfakesZTT; */
  /* std::vector<TH2D> CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSWfakesZTT; */
  /* std::vector<TH2D> CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSWfakesZTT; */
  /* std::vector<TH2D> CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSWfakesZTT; */
  /* std::vector<TH2D> CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSWfakesZTT; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSWfakesZTT; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSWfakesZTT; */
  /* std::vector<TH2D> CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSWfakesZTT; */
  /* std::vector<TH2D> CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSWfakesZTT; */
  /* std::vector<TH2D> CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSWfakesZTT; */
  /* std::vector<TH2D> CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSWfakesZTT; */
  /* std::vector<TH2D> CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSWfakesZTT; */
  /* std::vector<TH2D> CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSWfakesZTT; */
  /* std::vector<TH2D> CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSWfakesZTT; */
  /* std::vector<TH2D> CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSWfakesZTT; */
  /* std::vector<TH2D> CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSWfakesZTT; */
  std::vector<TH2D> CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesZTT;
  std::vector<TH2D> CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSWfakesZTT;
  std::vector<TH2D> CMS_scale_gg_13TeVUpPVRefitWithTracksBSWfakesZTT;
  std::vector<TH2D> CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSWfakesZTT;
  std::vector<TH2D> CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSWfakesZTT;
  std::vector<TH2D> CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesZTT;
  std::vector<TH2D> CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSWfakesZTT;
  std::vector<TH2D> CMS_scale_gg_13TeVDownPVRefitWithTracksBSWfakesZTT;
  std::vector<TH2D> CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSWfakesZTT;
  std::vector<TH2D> CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSWfakesZTT;
  std::vector<TH2D> PrefiringUpPVRefitWithTracksBSWfakesZTT;
  std::vector<TH2D> PrefiringDownPVRefitWithTracksBSWfakesZTT;


  std::vector<TH2D> ttbarcontaminationHiggs;
  std::vector<TH2D> ttbarcontaminationJetFakes;
  std::vector<TH2D> ttbarcontaminationZTT;

  std::vector<TH2D> ttbarcontaminationHiggs_DP;
  std::vector<TH2D> ttbarcontaminationJetFakes_DP;
  std::vector<TH2D> ttbarcontaminationZTT_DP;
  
  std::vector<TH2D> ttbarcontaminationWfakesHiggs;
  std::vector<TH2D> ttbarcontaminationWfakesJetFakes;
  std::vector<TH2D> ttbarcontaminationWfakesZTT;

  std::vector<TH2D> ttbarcontaminationWfakesHiggs_DP;
  std::vector<TH2D> ttbarcontaminationWfakesJetFakes_DP;
  std::vector<TH2D> ttbarcontaminationWfakesZTT_DP;

  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC;
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSHiggsQCDMC; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSHiggsQCDMC; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSHiggsQCDMC; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSHiggsQCDMC; */

  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC;
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSJetFakesQCDMC; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSJetFakesQCDMC; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSJetFakesQCDMC; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSJetFakesQCDMC; */

  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC;
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSZTTQCDMC; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSZTTQCDMC; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSZTTQCDMC; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSZTTQCDMC; */


    std::vector<TH2D> CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSHiggs_DP;
  /* std::vector<TH2D> CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> CMS_res_j_13TeVUpPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> CMS_res_j_13TeVDownPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSHiggs_DP; */
  std::vector<TH2D> CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> CMS_scale_gg_13TeVUpPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> CMS_scale_gg_13TeVDownPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggs_DP;
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSHiggs_DP; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSHiggs_DP; */
 
  std::vector<TH2D> PrefiringUpPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> PrefiringDownPVRefitWithTracksBSHiggs_DP;
  
  std::vector<TH2D> polarimetricAcopAnglePVRefitWithTracksBSMVADMWfakesHiggs_DP;
  std::vector<TH2D> CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP;
  std::vector<TH2D> CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP;
  std::vector<TH2D> CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP;
  std::vector<TH2D> CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP;
  std::vector<TH2D> CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP;
  std::vector<TH2D> CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP;
  /* std::vector<TH2D> CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP; */
  /* std::vector<TH2D> CMS_res_j_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP; */
  /* std::vector<TH2D> CMS_res_j_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP; */
  /* std::vector<TH2D> CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP; */
  /* std::vector<TH2D> CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP; */
  /* std::vector<TH2D> CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP; */
  /* std::vector<TH2D> CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP; */
  /* std::vector<TH2D> CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP; */
  std::vector<TH2D> CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP;
  std::vector<TH2D> CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP;
  std::vector<TH2D> CMS_scale_gg_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP;
  std::vector<TH2D> CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP;
  std::vector<TH2D> CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP;
  std::vector<TH2D> CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP;
  std::vector<TH2D> CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP;
  std::vector<TH2D> CMS_scale_gg_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP;
  std::vector<TH2D> CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP;
  std::vector<TH2D> CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP;
  std::vector<TH2D> PrefiringUpPVRefitWithTracksBSWfakesHiggs_DP;
  std::vector<TH2D> PrefiringDownPVRefitWithTracksBSWfakesHiggs_DP;
  
  std::vector<TH2D> CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSJetFakes_DP;
  /* std::vector<TH2D> CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> CMS_res_j_13TeVUpPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> CMS_res_j_13TeVDownPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSJetFakes_DP; */
  std::vector<TH2D> CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> CMS_scale_gg_13TeVUpPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> CMS_scale_gg_13TeVDownPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakes_DP;
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSJetFakes_DP; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSJetFakes_DP; */
  
  std::vector<TH2D> PrefiringUpPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> PrefiringDownPVRefitWithTracksBSJetFakes_DP;
 
  std::vector<TH2D> polarimetricAcopAnglePVRefitWithTracksBSMVADMWfakesJetFakes_DP;
  std::vector<TH2D> CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP;
  std::vector<TH2D> CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP;
  std::vector<TH2D> CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP;
  std::vector<TH2D> CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP;
  std::vector<TH2D> CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP;
  std::vector<TH2D> CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP;
  /* std::vector<TH2D> CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP; */
  /* std::vector<TH2D> CMS_res_j_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP; */
  /* std::vector<TH2D> CMS_res_j_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP; */
  /* std::vector<TH2D> CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP; */
  /* std::vector<TH2D> CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP; */
  /* std::vector<TH2D> CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP; */
  /* std::vector<TH2D> CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP; */
  /* std::vector<TH2D> CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP; */
  std::vector<TH2D> CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP;
  std::vector<TH2D> CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP;
  std::vector<TH2D> CMS_scale_gg_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP;
  std::vector<TH2D> CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP;
  std::vector<TH2D> CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP;
  std::vector<TH2D> CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP;
  std::vector<TH2D> CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP;
  std::vector<TH2D> CMS_scale_gg_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP;
  std::vector<TH2D> CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP;
  std::vector<TH2D> CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP;
  std::vector<TH2D> PrefiringUpPVRefitWithTracksBSWfakesJetFakes_DP;
  std::vector<TH2D> PrefiringDownPVRefitWithTracksBSWfakesJetFakes_DP;
  
                    
  std::vector<TH2D> CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSZTT_DP;
  /* std::vector<TH2D> CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> CMS_res_j_13TeVUpPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> CMS_res_j_13TeVDownPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSZTT_DP; */
  std::vector<TH2D> CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> CMS_scale_gg_13TeVUpPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> CMS_scale_gg_13TeVDownPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTT_DP;
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSZTT_DP; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSZTT_DP; */
 
  std::vector<TH2D> PrefiringUpPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> PrefiringDownPVRefitWithTracksBSZTT_DP;
  
  std::vector<TH2D> polarimetricAcopAnglePVRefitWithTracksBSMVADMWfakesZTT_DP;
  std::vector<TH2D> CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSWfakesZTT_DP;
  std::vector<TH2D> CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSWfakesZTT_DP;
  std::vector<TH2D> CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSWfakesZTT_DP;
  std::vector<TH2D> CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSWfakesZTT_DP;
  std::vector<TH2D> CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSWfakesZTT_DP;
  std::vector<TH2D> CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSWfakesZTT_DP;
  /* std::vector<TH2D> CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSWfakesZTT_DP; */
  /* std::vector<TH2D> CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSWfakesZTT_DP; */
  /* std::vector<TH2D> CMS_res_j_13TeVUpPVRefitWithTracksBSWfakesZTT_DP; */
  /* std::vector<TH2D> CMS_res_j_13TeVDownPVRefitWithTracksBSWfakesZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSWfakesZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSWfakesZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSWfakesZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSWfakesZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSWfakesZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSWfakesZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSWfakesZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSWfakesZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSWfakesZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSWfakesZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSWfakesZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSWfakesZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSWfakesZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSWfakesZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSWfakesZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSWfakesZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSWfakesZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSWfakesZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSWfakesZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSWfakesZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSWfakesZTT_DP; */
  /* std::vector<TH2D> CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSWfakesZTT_DP; */
  /* std::vector<TH2D> CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSWfakesZTT_DP; */
  /* std::vector<TH2D> CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSWfakesZTT_DP; */
  /* std::vector<TH2D> CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSWfakesZTT_DP; */
  /* std::vector<TH2D> CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSWfakesZTT_DP; */
  /* std::vector<TH2D> CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSWfakesZTT_DP; */
  /* std::vector<TH2D> CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSWfakesZTT_DP; */
  std::vector<TH2D> CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesZTT_DP;
  std::vector<TH2D> CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSWfakesZTT_DP;
  std::vector<TH2D> CMS_scale_gg_13TeVUpPVRefitWithTracksBSWfakesZTT_DP;
  std::vector<TH2D> CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSWfakesZTT_DP;
  std::vector<TH2D> CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSWfakesZTT_DP;
  std::vector<TH2D> CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesZTT_DP;
  std::vector<TH2D> CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSWfakesZTT_DP;
  std::vector<TH2D> CMS_scale_gg_13TeVDownPVRefitWithTracksBSWfakesZTT_DP;
  std::vector<TH2D> CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSWfakesZTT_DP;
  std::vector<TH2D> CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSWfakesZTT_DP;
  std::vector<TH2D> PrefiringUpPVRefitWithTracksBSWfakesZTT_DP;
  std::vector<TH2D> PrefiringDownPVRefitWithTracksBSWfakesZTT_DP;
  
 
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP;
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSHiggsQCDMC_DP; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSHiggsQCDMC_DP; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSHiggsQCDMC_DP; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSHiggsQCDMC_DP; */
 
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP;
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSJetFakesQCDMC_DP; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSJetFakesQCDMC_DP; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSJetFakesQCDMC_DP; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSJetFakesQCDMC_DP; */
 
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP;
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSZTTQCDMC_DP; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSZTTQCDMC_DP; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSZTTQCDMC_DP; */
  /* std::vector<TH2D> jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSZTTQCDMC_DP; */


  std::vector<TH2D> jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSHiggs;
  std::vector<TH2D> jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSHiggs;
  std::vector<TH2D> jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSZTT;
  std::vector<TH2D> jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSZTT;
  
  std::vector<TH2D> jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSHiggsQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSHiggsQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSJetFakesQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSJetFakesQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSZTTQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSZTTQCDMC;
 
  std::vector<TH2D> jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSZTT_DP;
  
  std::vector<TH2D> jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSHiggsQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSHiggsQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSJetFakesQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSJetFakesQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSZTTQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSZTTQCDMC_DP;

  
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSHiggs;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSHiggs;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSHiggs;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSHiggs;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSHiggs;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSHiggs;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSHiggs;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSHiggs;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSHiggs;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSHiggs;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSHiggs;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSHiggs;
  
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSJetFakes;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSJetFakes;
  
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSZTT;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSZTT;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSZTT;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSZTT;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSZTT;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSZTT;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSZTT;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSZTT;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSZTT;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSZTT;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSZTT;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSZTT;
  
  
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSHiggsQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSHiggsQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSHiggsQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSHiggsQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSHiggsQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSHiggsQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSHiggsQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSHiggsQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSHiggsQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSHiggsQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSHiggsQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSHiggsQCDMC;
  
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSJetFakesQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSJetFakesQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSJetFakesQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSJetFakesQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSJetFakesQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSJetFakesQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSJetFakesQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSJetFakesQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSJetFakesQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSJetFakesQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSJetFakesQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSJetFakesQCDMC;
  
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSZTTQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSZTTQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSZTTQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSZTTQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSZTTQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSZTTQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSZTTQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSZTTQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSZTTQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSZTTQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSZTTQCDMC;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSZTTQCDMC;
  
  
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSHiggs_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSHiggs_DP;
  
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSJetFakes_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSJetFakes_DP;
  
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSZTT_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSZTT_DP;
  
  
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSHiggsQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSHiggsQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSHiggsQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSHiggsQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSHiggsQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSHiggsQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSHiggsQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSHiggsQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSHiggsQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSHiggsQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSHiggsQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSHiggsQCDMC_DP;
  
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSJetFakesQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSJetFakesQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSJetFakesQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSJetFakesQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSJetFakesQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSJetFakesQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSJetFakesQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSJetFakesQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSJetFakesQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSJetFakesQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSJetFakesQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSJetFakesQCDMC_DP;
  
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSZTTQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSZTTQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSZTTQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSZTTQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSZTTQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSZTTQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSZTTQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSZTTQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSZTTQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSZTTQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSZTTQCDMC_DP;
  std::vector<TH2D> jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSZTTQCDMC_DP;



  std::vector<TH1D> Tau1SVXminusPVX;
  std::vector<TH1D> Tau1SVYminusPVY;
  std::vector<TH1D> Tau1SVZminusPVZ;
  std::vector<TH1D> Tau1SVXminusPVXRelative;
  std::vector<TH1D> Tau1SVYminusPVYRelative;
  std::vector<TH1D> Tau1SVZminusPVZRelative;
  
  std::vector<TH1D> Tau1SVXminusPVXQCDMC;
  std::vector<TH1D> Tau1SVYminusPVYQCDMC;
  std::vector<TH1D> Tau1SVZminusPVZQCDMC;
  std::vector<TH1D> Tau1SVXminusPVXRelativeQCDMC;
  std::vector<TH1D> Tau1SVYminusPVYRelativeQCDMC;
  std::vector<TH1D> Tau1SVZminusPVZRelativeQCDMC;
  

  std::vector<TH1D> Tau2SVXminusPVX;
  std::vector<TH1D> Tau2SVYminusPVY;
  std::vector<TH1D> Tau2SVZminusPVZ;
  std::vector<TH1D> Tau2SVXminusPVXRelative;
  std::vector<TH1D> Tau2SVYminusPVYRelative;
  std::vector<TH1D> Tau2SVZminusPVZRelative;
  
  std::vector<TH1D> Tau2SVXminusPVXQCDMC;
  std::vector<TH1D> Tau2SVYminusPVYQCDMC;
  std::vector<TH1D> Tau2SVZminusPVZQCDMC;
  std::vector<TH1D> Tau2SVXminusPVXRelativeQCDMC;
  std::vector<TH1D> Tau2SVYminusPVYRelativeQCDMC;
  std::vector<TH1D> Tau2SVZminusPVZRelativeQCDMC;
};

#endif
