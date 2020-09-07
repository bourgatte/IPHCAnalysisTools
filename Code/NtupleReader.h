//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon May 15 17:49:09 2017 by ROOT version 5.34/18
// from TTree HTauTauTree/HTauTauTree
// found on file: /home-pbs/vcherepa/cms_work/CMSSW_8_0_25/src/LLRHiggsTauTau/NtupleProducer/test/HTauTauAnalysis.root
//////////////////////////////////////////////////////////

#ifndef NtupleReader_h
#define NtupleReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>
#include <vector>
#include <vector>
#include <vector>
#include <vector>
#include <vector>
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

// Fixed size dimensions of array or collections stored in the TTree if any.

class NtupleReader {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Declaration of leaf types
  ULong64_t       EventNumber;
  Int_t           RunNumber;
  Int_t           lumi;
  Int_t           year;
  //Int_t           NBadMu;
  Bool_t          passecalBadCalibFilterUpdate;
  Float_t         prefiringweight;
  Float_t         prefiringweightup;
  Float_t         prefiringweightdown;
  Long64_t        triggerbit;
  Int_t           metfilterbit;
  //Float_t         met;
  Float_t         met_er;
  Float_t         met_er_phi;
  //Float_t         metphi;
  //Float_t         PUPPImet;
  //Float_t         PUPPImetphi;
  vector<float>   *daughters_IetaIeta;
  vector<float>   *daughters_full5x5_IetaIeta;
  vector<float>   *daughters_hOverE;
  vector<float>   *daughters_deltaEtaSuperClusterTrackAtVtx;
  vector<float>   *daughters_deltaPhiSuperClusterTrackAtVtx;
  vector<float>   *daughters_IoEmIoP;
  vector<float>   *daughters_IoEmIoP_ttH;
  //Float_t         PFMETCov00;
  //Float_t         PFMETCov01;
  //Float_t         PFMETCov10;
  //Float_t         PFMETCov11;
  Float_t         PFMETsignif;
  Int_t           npv;
  Float_t         npu;
  //Float_t         PUReweight;
  Float_t         rho;
  vector<float>   *mothers_px;
  vector<float>   *mothers_py;
  vector<float>   *mothers_pz;
  vector<float>   *mothers_e;
  vector<Long64_t> *mothers_trgSeparateMatch;
  vector<string>  *trigger_name;
  vector<int>     *trigger_accept;
  vector<float>   *daughters_px;
  vector<float>   *daughters_py;
  vector<float>   *daughters_pz;
  vector<float>   *daughters_e;
  vector<int>     *daughters_charge;
  vector<float>   *daughters_charged_px;
  vector<float>   *daughters_charged_py;
  vector<float>   *daughters_charged_pz;
  vector<float>   *daughters_charged_e;
  vector<float>   *daughters_neutral_px;
  vector<float>   *daughters_neutral_py;
  vector<float>   *daughters_neutral_pz;
  vector<float>   *daughters_neutral_e;
  vector<int>     *daughters_hasTES;
  /* vector<float>   *daughters_px_TauUp; */
  /* vector<float>   *daughters_py_TauUp; */
  /* vector<float>   *daughters_pz_TauUp; */
  /* vector<float>   *daughters_e_TauUp; */
  vector<int>     *daughters_hasEES;
  /* vector<float>   *daughters_px_TauDown; */
  /* vector<float>   *daughters_py_TauDown; */
  /* vector<float>   *daughters_pz_TauDown; */
  /* vector<float>   *daughters_e_TauDown; */
  /* vector<float>   *daughters_px_EleUp; */
  /* vector<float>   *daughters_py_EleUp; */
  /* vector<float>   *daughters_pz_EleUp; */
  /* vector<float>   *daughters_e_EleUp; */
  /* vector<float>   *daughters_px_EleDown; */
  /* vector<float>   *daughters_py_EleDown; */
  /* vector<float>   *daughters_pz_EleDown; */
  /* vector<float>   *daughters_e_EleDown; */
  vector<int>     *daughters_isTauMatched;
  Int_t           PUNumInteractions;
  vector<int>     *daughters_genindex;
  Double_t         MC_weight;
  Double_t         MC_weight_scale_muF0p5;
  Double_t         MC_weight_scale_muF2;
  Double_t         MC_weight_scale_muR0p5;
  Double_t         MC_weight_scale_muR2;
  Float_t         nominal_wt;
  vector<double> *TheoreticalPSUnc;
  Double_t          TheoreticalScaleUnc1005;
  Double_t          TheoreticalScaleUnc1009;
  Double_t          TheoreticalScaleUnc5;
  Double_t          TheoreticalScaleUnc9;
  Float_t         lheHt;
  Int_t           lheNOutPartons;
  Int_t           lheNOutB;
  Int_t           lheNOutC;
  //  Float_t         aMCatNLOweight;
  vector<float>   *genpart_px;
  vector<float>   *genpart_py;
  vector<float>   *genpart_pz;
  vector<float>   *genpart_e;
  Int_t           DataMC_Type_idx;
  bool           Event_isRealData;
  vector<float>   *genpart_pca_x;
  vector<float>   *genpart_pca_y;
  vector<float>   *genpart_pca_z;
  vector<int>     *genpart_pdg;
  vector<int>     *genpart_status;
  vector<int>     *genpart_HMothInd;
  vector<int>     *genpart_MSSMHMothInd;
  vector<int>     *genpart_TopMothInd;
  vector<int>     *genpart_TauMothInd;
  vector<int>     *genpart_ZMothInd;
  vector<int>     *genpart_WMothInd;
  vector<int>     *genpart_bMothInd;
  vector<int>     *genpart_HZDecayMode;
  vector<int>     *genpart_TopDecayMode;
  vector<int>     *genpart_WDecayMode;
  vector<int>     *genpart_TauGenDecayMode;
  vector<int>     *genpart_TauGenDetailedDecayMode;
  vector<int>     *genpart_flags;
  vector<float>   *genjet_px;
  vector<float>   *genjet_py;
  vector<float>   *genjet_pz;
  vector<float>   *genjet_e;
  vector<int>     *genjet_partonFlavour;
  vector<int>     *genjet_hadronFlavour;
  Int_t           NUP;
  /* vector<float>   *SVfit_fitMETPhiTauUp; */
  /* vector<float>   *SVfit_fitMETPhiTauDown; */
  /* vector<float>   *SVfit_fitMETRhoTauUp; */
  /* vector<float>   *SVfit_fitMETRhoTauDown; */
  /* vector<float>   *SVfit_phiUncTauUp; */
  /* vector<float>   *SVfit_phiUncTauDown; */
  /* vector<float>   *SVfit_phiTauUp; */
  /* vector<float>   *SVfit_phiTauDown; */
  /* vector<float>   *SVfit_etaUncTauUp; */
  /* vector<float>   *SVfit_etaUncTauDown; */
  /* vector<float>   *SVfit_etaTauUp; */
  /* vector<float>   *SVfit_etaTauDown; */
  /* vector<float>   *SVfit_ptUncTauUp; */
  /* vector<float>   *SVfit_ptUncTauDown; */
  /* vector<float>   *SVfit_ptTauUp; */
  /* vector<float>   *SVfit_ptTauDown; */
  /* vector<float>   *SVfitTransverseMassTauUp; */
  /* vector<float>   *SVfitTransverseMassTauDown; */
  /* vector<float>   *SVfitMassTauUp; */
  /* vector<float>   *SVfitMassTauDown; */
  vector<vector<double> > *MCSignalParticle_p4;
  vector<int>     *MCSignalParticle_pdgid;
  vector<int>     *MCSignalParticle_charge;
  vector<vector<double> > *MCSignalParticle_Poca;
  vector<vector<unsigned int> > *MCSignalParticle_Tauidx;
  vector<vector<vector<double> > > *MCTauandProd_p4;
  vector<vector<vector<double> > > *MCTauandProd_Vertex;
  vector<vector<int> > *MCTauandProd_pdgid;
  vector<vector<unsigned int> > *MCTauandProd_midx;
  vector<vector<int> > *MCTauandProd_charge;
  vector<unsigned int> *MCTau_JAK;
  vector<unsigned int> *MCTau_DecayBitMask;
  vector<vector<float> > *MC_p4;
  vector<int>     *MC_pdgid;
  vector<int>     *MC_charge;
  vector<int>     *MC_midx;
  vector<vector<int> > *MC_childpdgid;
  vector<vector<int> > *MC_childidx;
  vector<int>     *MC_status;
  /* vector<float>   *SVfitMass; */
  /* vector<float>   *SVfitTransverseMass; */
  /* vector<float>   *SVfit_pt; */
  /* vector<float>   *SVfit_ptUnc; */
  /* vector<float>   *SVfit_eta; */
  /* vector<float>   *SVfit_etaUnc; */
  /* vector<float>   *SVfit_phi; */
  /* vector<float>   *SVfit_phiUnc; */
  /* vector<float>   *SVfit_fitMETRho; */
  /* vector<float>   *SVfit_fitMETPhi; */
  //  vector<bool>    *isOSCand;
  vector<float>   *METx;
  vector<float>   *METy;
  /* vector<float> *metx_up_jes; */
  /* vector<float> *mety_up_jes; */
  /* vector<float> *metx_down_jes; */
  /* vector<float> *mety_down_jes; */
  /* vector<float> *metx_up_tes; */
  /* vector<float> *mety_up_tes; */
  /* vector<float> *metx_down_tes; */
  /* vector<float> *mety_down_tes; */
  /* vector<float> *metx_up_ees; */
  /* vector<float> *mety_up_ees; */
  /* vector<float> *metx_down_ees; */
  /* vector<float> *mety_down_ees; */
  /* //vector<float>   *uncorrMETx; */
  /* //vector<float>   *uncorrMETy; */
  /* vector<float>   *MET_cov00; */
  /* vector<float>   *MET_cov01; */
  /* vector<float>   *MET_cov10; */
  /* vector<float>   *MET_cov11; */
  /* vector<float>   *MET_significance; */
  /* vector<float>   *mT_Dau1; */
  /* vector<float>   *mT_Dau2; */
  vector<int>     *PDGIdDaughters;
  vector<int>     *indexDau1;
  vector<int>     *indexDau2;
  vector<int>     *particleType;
  //vector<float>   *discriminator;
  vector<int>     *daughters_muonID;
  vector<int>     *daughters_typeOfMuon;
  vector<float>   *dxy;
  vector<float>   *dz;
  vector<float>   *dxy_innerTrack;
  //  vector<float>   *dz_innerTrack;
  vector<float>   *daughters_rel_error_trackpt;
  vector<float>   *SIP;
  //vector<bool>    *daughters_iseleBDT;
  std::vector<bool> *daughters_iseleWPLoose;
  vector<bool>    *daughters_iseleWP80;
  vector<bool>    *daughters_iseleWP90;
  std::vector<bool> *daughters_iseleNoIsoWPLoose;
  std::vector<bool> *daughters_iseleNoIsoWP80;
  std::vector<bool> *daughters_iseleNoIsoWP90;
  vector<float>   *daughters_eleMVAnt;
  vector<float>   *daughters_eleMVA_HZZ;
  vector<bool>    *daughters_passConversionVeto;
  vector<int>     *daughters_eleMissingHits;
  vector<bool>    *daughters_iseleChargeConsistent;
  //vector<int>     *daughters_eleCUTID;
  vector<int>     *decayMode;
  vector<int>     *genmatch;
  vector<Long64_t> *tauID;
  //  vector<int> *MVADM2016;
  vector<int> *MVADM2017;
  vector<float>   *combreliso;
  vector<float>   *combreliso03;
  vector<float>   *daughters_depositR03_tracker;
  vector<float>   *daughters_depositR03_ecal;
  vector<float>   *daughters_depositR03_hcal;
  vector<int>     *daughters_decayModeFindingOldDMs;
  //vector<float>   *daughters_SCeta;
  //vector<float>   *againstElectronMVA5category;
  //vector<float>   *againstElectronMVA5raw;
  //vector<float>   *byPileupWeightedIsolationRaw3Hits;
  vector<float>   *footprintCorrection;
  vector<float>   *neutralIsoPtSumWeight;
  vector<float>   *photonPtSumOutsideSignalCone;
  vector<int>     *daughters_decayModeFindingNewDMs;
  //vector<float>   *daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits;
  //vector<float>   *daughters_byIsolationMVArun2017v2DBoldDMwLTraw2017;
  //vector<float>   *daughters_byIsolationMVArun2017v1DBoldDMwLTraw2017;
  //vector<float>   *daughters_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017;
  //vector<int>     *daughters_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017;
  //vector<float>   *daughters_byIsolationMVA3oldDMwoLTraw;
  //vector<float>   *daughters_byIsolationMVA3oldDMwLTraw;
  vector<float>   *daughters_byIsolationMVA3newDMwoLTraw;
  vector<float>   *daughters_byIsolationMVA3newDMwLTraw;
  //  vector<float>   *daughters_byIsolationMVArun2v1DBoldDMwLTraw;
  vector<float>   *daughters_byDeepTau2017v2p1VSjetraw;
  vector<float>   *daughters_byDeepTau2017v2p1VSeraw;
  vector<float>   *daughters_byDeepTau2017v2p1VSmuraw;
  vector<float>   *daughters_chargedIsoPtSum;
  vector<float>   *daughters_neutralIsoPtSum;
  vector<float>   *daughters_puCorrPtSum;
  vector<int>     *daughters_numChargedParticlesSignalCone;
  vector<int>     *daughters_numNeutralHadronsSignalCone;
  vector<int>     *daughters_numPhotonsSignalCone;
  vector<int>     *daughters_daughters_numParticlesSignalCone;
  vector<int>     *daughters_numChargedParticlesIsoCone;
  vector<int>     *daughters_numNeutralHadronsIsoCone;
  vector<int>     *daughters_numPhotonsIsoCone;
  vector<int>     *daughters_numParticlesIsoCone;
  vector<float>   *daughters_leadChargedParticlePt;
  vector<float>   *daughters_trackRefPt;
  //vector<int>     *daughters_isLastTriggerObjectforPath;
  vector<Long64_t> *daughters_trgMatched;
  //vector<int>     *daughters_isTriggerObjectforPath;
  vector<Long64_t> *daughters_FilterFired;
  vector<Long64_t> *daughters_isGoodTriggerType;
  vector<Long64_t> *daughters_L3FilterFired;
  vector<Long64_t> *daughters_L3FilterFiredLast;
  vector<float>   *daughters_HLTpt;
  //vector<bool>    *daughters_isL1IsoTau28Matched;
  vector<int>     *Muon_trackCharge;
  vector<int>     *Muon_pdgid;
  vector<double>  *Muon_B;
  vector<double>  *Muon_M;
  vector<vector<double> > *Muon_par;
  //  vector<vector<double> > *Muon_cov;
  /* std::vector<std::vector<double> > *PFTau_Track_par; */
  /* std::vector<std::vector<double> > *PFTau_Track_cov; */
  /* std::vector<int>  *PFTau_Track_charge; */
  /* std::vector<int>  *PFTau_Track_pdgid; */
  /* std::vector<double>  *PFTau_Track_B; */
  /* std::vector<double>  *PFTau_Track_M; */
  vector<vector<double> > *PFTauSVPos;
  vector<vector<double> > *PFTauSVCov;
  vector<vector<vector<double> > > *PFTauPionsP4;
  vector<vector<vector<double> > > *PFTauRefitPionsP4;
  vector<vector<double> > *PFTauPionsCharge;
  vector<vector<double> > *PFTauRefitPionsCharge;
  vector<vector<double> > *PFTauSVChi2NDofMatchingQuality;
  std::vector<Float_t> *TauFLSignificance;
  std::vector<double>  *PFTauGEOMFlightLenght;
  std::vector<double>  *PFTauGEOMFlightLenghtSignificance;
  vector<int>     *daughters_jetNDauChargedMVASel;
  vector<float>   *daughters_miniRelIsoCharged;
  vector<float>   *daughters_miniRelIsoNeutral;
  vector<float>   *daughters_jetPtRel;
  vector<float>   *daughters_jetPtRatio;
  vector<float>   *daughters_jetBTagCSV;
  std::vector<Float_t> *daughters_jetBTagDeepCSV;
  std::vector<Float_t> *daughters_jetBTagDeepFlavor;
  //vector<float>   *daughters_lepMVA_mvaId;
  vector<float>   *daughters_pca_x;
  vector<float>   *daughters_pca_y;
  vector<float>   *daughters_pca_z;
  vector<float>   *daughters_pcaRefitPV_x;
  vector<float>   *daughters_pcaRefitPV_y;
  vector<float>   *daughters_pcaRefitPV_z;
  vector<float>   *daughters_pcaGenPV_x;
  vector<float>   *daughters_pcaGenPV_y;
  vector<float>   *daughters_pcaGenPV_z;
  vector<vector<double> > *PFTau_a1_lvp;
  //  vector<vector<double> > *PFTau_a1_cov;
  vector<int>     *PFTau_a1_charge;
  vector<int>     *PFTau_a1_pdgid;
  vector<double>  *PFTau_a1_B;
  vector<double>  *PFTau_a1_M;
  std::vector<float> *PFTauTrack_deltaR;
  std::vector<std::vector<double > > *PFTauLeadTrackLV;
  //Int_t           JetsNumber;
  std::vector<Long64_t> *jets_VBFleadFilterMatch;
  std::vector<Long64_t> *jets_VBFsubleadFilterMatch;
  vector<float>   *jets_px;
  vector<float>   *jets_py;
  vector<float>   *jets_pz;
  vector<float>   *jets_e;
  vector<float>   *jetsDown_px;
  vector<float>   *jetsDown_py;
  vector<float>   *jetsDown_pz;
  vector<float>   *jetsDown_e;
  vector<float>   *jetsUp_px;
  vector<float>   *jetsUp_py;
  vector<float>   *jetsUp_pz;
  vector<float>   *jetsUp_e;
  
  vector<float>   *jetRawf;
  vector<float>   *jets_area;
  std::vector<Float_t> *jets_JER;
  //  std::vector<Float_t> *jets_QGdiscr;
  vector<float>   *jets_mT;
  vector<int>     *jets_Flavour; 
  vector<int>     *jets_HadronFlavour;
  vector<int>     *jets_genjetIndex;
  vector<float>   *jets_PUJetID;
  vector<float>   *jets_PUJetIDupdated;
  vector<Int_t>   *jets_PUJetID_WP;
  vector<Int_t>   *jets_PUJetIDupdated_WP;
  vector<float>   *jets_vtxPt;
  vector<float>   *jets_vtxMass;
  vector<float>   *jets_vtx3dL;
  vector<float>   *jets_vtxNtrk;
  vector<float>   *jets_vtx3deL;
  vector<float>   *jets_leadTrackPt;
  vector<float>   *jets_leptonPtRel;
  vector<float>   *jets_leptonPt;
  vector<float>   *jets_leptonDeltaR;
  vector<float>   *_jets_chEmEF;
  vector<float>   *_jets_chHEF;
  vector<float>   *_jets_nEmEF;
  vector<float>   *_jets_nHEF;
  vector<float>   *jets_MUF;
  vector<int>     *_jets_neMult;
  vector<int>     *_jets_chMult;
  vector<float>   *jets_jecUnc;
  vector<float>   *jets_jetUnc_Absolute_up;
  vector<float>   *jets_jetUnc_FlavorQCD_up;
  vector<float>   *jets_jetUnc_RelativeBal_up;
  vector<float>   *jets_jetUnc_HF_up;
  vector<float>   *jets_jetUnc_BBEC1_up;
  vector<float>   *jets_jetUnc_EC2_up;
  vector<float>   *jets_jetUnc_BBEC1_YEAR_up;
  vector<float>   *jets_jetUnc_EC2_YEAR_up;
  vector<float>   *jets_jetUnc_Absolute_YEAR_up;
  vector<float>   *jets_jetUnc_HF_YEAR_up;
  vector<float>   *jets_jetUnc_RelativeSample_YEAR_up;
  vector<float>   *jets_jetUnc_Absolute_dw;
  vector<float>   *jets_jetUnc_FlavorQCD_dw;
  vector<float>   *jets_jetUnc_RelativeBal_dw;
  vector<float>   *jets_jetUnc_HF_dw;
  vector<float>   *jets_jetUnc_BBEC1_dw;
  vector<float>   *jets_jetUnc_EC2_dw;
  vector<float>   *jets_jetUnc_BBEC1_YEAR_dw;
  vector<float>   *jets_jetUnc_EC2_YEAR_dw;
  vector<float>   *jets_jetUnc_Absolute_YEAR_dw;
  vector<float>   *jets_jetUnc_HF_YEAR_dw;
  vector<float>   *jets_jetUnc_RelativeSample_YEAR_dw;
  vector<float>   *bDiscriminator;
  vector<float>   *bCSVscore;
  vector<float>   *pfCombinedMVAV2BJetTags;
  std::vector<Float_t> *bDeepCSV_probb; //DeepCSV_probb
  std::vector<Float_t> *bDeepCSV_probbb; //DeepCSV_probbb
  std::vector<Float_t> *bDeepCSV_probudsg; //DeepCSV_probudsg
  std::vector<Float_t> *bDeepCSV_probc; //DeepCSV_probc
  std::vector<Float_t> *bDeepCSV_probcc; //DeepCSV_probcc
  std::vector<Float_t> *bDeepFlavor_probb;  //DeepFlavor_probb
  std::vector<Float_t> *bDeepFlavor_probbb; //DeepFlavor_probbb
  std::vector<Float_t> *bDeepFlavor_problepb; //DeepFlavor_problepb
  std::vector<Float_t> *bDeepFlavor_probc; //DeepFlavor_probc
  std::vector<Float_t> *bDeepFlavor_probuds; //DeepFlavor_probuds
  std::vector<Float_t> *bDeepFlavor_probg; //DeepFlavor_probg
  std::vector<Bool_t> *looseJetID;
  std::vector<Bool_t> *tightJetID;
  std::vector<Bool_t> *tightLepVetoJetID;
  //vector<int>     *PFjetID;
  vector<float>   *ak8jets_px;
  vector<float>   *ak8jets_py;
  vector<float>   *ak8jets_pz;
  vector<float>   *ak8jets_e;
  vector<float>   *ak8jets_SoftDropMass;
  /* vector<float>   *ak8jets_PrunedMass; */
  /* vector<float>   *ak8jets_TrimmedMass; */
  /* vector<float>   *ak8jets_FilteredMass; */
  vector<float>   *ak8jets_tau1;
  vector<float>   *ak8jets_tau2;
  vector<float>   *ak8jets_tau3;
  vector<float>   *ak8jets_tau4;
  vector<float>   *ak8jets_CSV;
  std::vector<Float_t> *ak8jets_deepCSV_probb; // CSV score
  std::vector<Float_t> *ak8jets_deepCSV_probbb; // CSV score
  /* std::vector<Float_t> *ak8jets_deepFlavor_probb; // Flavor score */
  /* std::vector<Float_t> *ak8jets_deepFlavor_probbb; // Flavor score */
  /* std::vector<Float_t> *ak8jets_deepFlavor_problepb; // Flavor score */
  vector<int>     *ak8jets_nsubjets;
  vector<float>   *subjets_px;
  vector<float>   *subjets_py;
  vector<float>   *subjets_pz;
  vector<float>   *subjets_e;
  vector<float>   *subjets_CSV;
  std::vector<Float_t> *subjets_deepCSV_probb;
  std::vector<Float_t> *subjets_deepCSV_probbb;
  /* std::vector<Float_t> *subjets_deepFlavor_probb; */
  /* std::vector<Float_t> *subjets_deepFlavor_probbb; */
  /* std::vector<Float_t> *subjets_deepFlavor_problepb; */
  vector<int>     *subjets_ak8MotherIdx;
  Float_t         pv_x;
  Float_t         pv_y;
  Float_t         pv_z;
  //  vector<double>  *pv_cov;

  vector<Float_t>         *RefitPVNoBS_x;
  vector<Float_t>         *RefitPVNoBS_y;
  vector<Float_t>         *RefitPVNoBS_z;
  vector<Float_t>         *RefitPVBS_x;
  vector<Float_t>         *RefitPVBS_y;
  vector<Float_t>         *RefitPVBS_z;

  Float_t          RefitPVWithTracksBS_x;
  Float_t          RefitPVWithTracksBS_y;
  Float_t          RefitPVWithTracksBS_z;

  Float_t          RefitPVWithTracksBS_xError;
  Float_t          RefitPVWithTracksBS_yError;
  Float_t          RefitPVWithTracksBS_zError;

  std::vector<Float_t> *RefitPVBS_xError,  *RefitPVBS_yError, *RefitPVBS_zError;
  std::vector<Float_t> *RefitPVNoBS_xError,  *RefitPVNoBS_yError, *RefitPVNoBS_zError;
   

  Float_t         pvGen_x;
  Float_t         pvGen_y;
  Float_t         pvGen_z;
  Bool_t          isRefitPV;
   
  vector<size_t>  *LeptonHash;
  vector<size_t>  *VertexHashNoBS1;
  vector<size_t>  *VertexHashNoBS2;
  vector<size_t>   *VertexHashBS1;
  vector<size_t>   *VertexHashBS2;

  Int_t SelectedPairs;
  Bool_t eleveto;
  Bool_t muonveto;
  std::vector<Bool_t> *trg_doubletau;
  
  std::vector<Float_t> *puppimt_1;
  std::vector<Float_t> *gen_match_1;
  
  std::vector<Float_t> *trigweight_1;
  std::vector<Float_t> *idisoweight_1;
  std::vector<Float_t> *antieweight_1;
  std::vector<Float_t> *antimuweight_1;

  std::vector<Float_t> *puppimt_2;
  std::vector<Float_t> *gen_match_2;

  std::vector<Float_t> *trigweight_2;
  std::vector<Float_t> *idisoweight_2;
  std::vector<Float_t> *antieweight_2;
  std::vector<Float_t> *antimuweight_2;
  
  /* std::vector<Float_t> *pt_tt; */
  /* std::vector<Float_t> *pt_vis; */
  std::vector<Float_t> *mt_tot;
  //  std::vector<Float_t> *m_vis;

  Float_t met;
  Float_t metphi;
  Float_t PUPPImet;
  Float_t puppimet_ex_UnclusteredEnUp;
  Float_t puppimet_ey_UnclusteredEnUp;
  Float_t puppimet_ex_UnclusteredEnDown;
  Float_t puppimet_ey_UnclusteredEnDown;
  Float_t PUPPImetphi;
  Float_t PFMETCov00;
  Float_t PFMETCov01;
  Float_t PFMETCov10;
  Float_t PFMETCov11;
  Float_t PUPPIMETCov00;
  Float_t PUPPIMETCov01;
  Float_t PUPPIMETCov10;
  Float_t PUPPIMETCov11;

  std::vector<Float_t> *mjj;
  std::vector<Float_t> *jdeta;
  std::vector<Float_t> *mjjDown;
  std::vector<Float_t> *jdetaDown;
  std::vector<Float_t> *mjjUp;
  std::vector<Float_t> *jdetaUp;
  
  std::vector<Int_t> *njetingap;
  std::vector<Int_t> *njetingap20;
  std::vector<Float_t> *jdphi;
  std::vector<Float_t> *dijetpt;
  std::vector<Float_t> *dijetphi;
  //  std::vector<Float_t> *ptvis;

  std::vector<Int_t> *nbtag;
  std::vector<Int_t> *njets;
  std::vector<Int_t> *njetspt20;
  std::vector<Int_t> *njetsUp;
  std::vector<Int_t> *njetspt20Up;
  std::vector<Int_t> *njetsDown;
  std::vector<Int_t> *njetspt20Down;
  std::vector<Float_t> *jpt_1;
  std::vector<Float_t> *jptDown_1;
  std::vector<Float_t> *jptUp_1;
  
  std::vector<Float_t> *jeta_1;
  std::vector<Float_t> *jphi_1;
  std::vector<Float_t> *jcsv_1;
  std::vector<Float_t> *jpt_2;
  std::vector<Float_t> *jeta_2;
  std::vector<Float_t> *jphi_2;
  std::vector<Float_t> *jcsv_2;
  std::vector<Float_t> *bpt_1;
  std::vector<Float_t> *beta_1;
  std::vector<Float_t> *bphi_1;
  std::vector<Float_t> *bcsv_1;
  std::vector<Float_t> *bpt_2;
  std::vector<Float_t> *beta_2;
  std::vector<Float_t> *bphi_2;
  std::vector<Float_t> *bcsv_2;

  Float_t puweight;

  std::vector<Float_t> *weight;

  /* std::vector<Float_t> *jpfid_1; */
  /* std::vector<Float_t> *jpuid_1; */
  /* std::vector<Float_t> *jpfid_2; */
  /* std::vector<Float_t> *jpuid_2; */
  /* std::vector<Float_t> *bpfid_1; */
  /* std::vector<Float_t> *bpuid_1; */
  /* std::vector<Float_t> *bpfid_2; */
  /* std::vector<Float_t> *bpuid_2; */
  
  //  std::vector<Float_t> *byDeepTau2017v2p1VSjetraw_1;
  std::vector<Float_t> *byVVVLooseDeepTau2017v2p1VSjet_1;
  std::vector<Float_t> *byVVLooseDeepTau2017v2p1VSjet_1; 
  std::vector<Float_t> *byVLooseDeepTau2017v2p1VSjet_1;  
  std::vector<Float_t> *byLooseDeepTau2017v2p1VSjet_1;   
  std::vector<Float_t> *byMediumDeepTau2017v2p1VSjet_1;  
  std::vector<Float_t> *byTightDeepTau2017v2p1VSjet_1;   
  std::vector<Float_t> *byVTightDeepTau2017v2p1VSjet_1;  
  std::vector<Float_t> *byVVTightDeepTau2017v2p1VSjet_1; 
  
  //  std::vector<Float_t> *byDeepTau2017v2p1VSeleraw_1;
  std::vector<Float_t> *byVVVLooseDeepTau2017v2p1VSe_1;  
  std::vector<Float_t> *byVVLooseDeepTau2017v2p1VSe_1; 
  std::vector<Float_t> *byVLooseDeepTau2017v2p1VSe_1;   
  std::vector<Float_t> *byLooseDeepTau2017v2p1VSe_1;	
  std::vector<Float_t> *byMediumDeepTau2017v2p1VSe_1;   
  std::vector<Float_t> *byTightDeepTau2017v2p1VSe_1;	
  std::vector<Float_t> *byVTightDeepTau2017v2p1VSe_1;   
  std::vector<Float_t> *byVVTightDeepTau2017v2p1VSe_1;
  
  //  std::vector<Float_t> *byDeepTau2017v2p1VSmuraw_1;
  std::vector<Float_t> *byVLooseDeepTau2017v2p1VSmu_1; 
  std::vector<Float_t> *byLooseDeepTau2017v2p1VSmu_1; 
  std::vector<Float_t> *byMediumDeepTau2017v2p1VSmu_1;   
  std::vector<Float_t> *byTightDeepTau2017v2p1VSmu_1;

  //  std::vector<Float_t> *byDeepTau2017v2p1VSjetraw_2;
  std::vector<Float_t> *byVVVLooseDeepTau2017v2p1VSjet_2;
  std::vector<Float_t> *byVVLooseDeepTau2017v2p1VSjet_2; 
  std::vector<Float_t> *byVLooseDeepTau2017v2p1VSjet_2;  
  std::vector<Float_t> *byLooseDeepTau2017v2p1VSjet_2;   
  std::vector<Float_t> *byMediumDeepTau2017v2p1VSjet_2;  
  std::vector<Float_t> *byTightDeepTau2017v2p1VSjet_2;   
  std::vector<Float_t> *byVTightDeepTau2017v2p1VSjet_2;  
  std::vector<Float_t> *byVVTightDeepTau2017v2p1VSjet_2;
  
  //  std::vector<Float_t> *byDeepTau2017v2p1VSeleraw_2;
  std::vector<Float_t> *byVVVLooseDeepTau2017v2p1VSe_2;  
  std::vector<Float_t> *byVVLooseDeepTau2017v2p1VSe_2; 
  std::vector<Float_t> *byVLooseDeepTau2017v2p1VSe_2;   
  std::vector<Float_t> *byLooseDeepTau2017v2p1VSe_2;	
  std::vector<Float_t> *byMediumDeepTau2017v2p1VSe_2;   
  std::vector<Float_t> *byTightDeepTau2017v2p1VSe_2;	
  std::vector<Float_t> *byVTightDeepTau2017v2p1VSe_2;   
  std::vector<Float_t> *byVVTightDeepTau2017v2p1VSe_2;
  
  //  std::vector<Float_t> *byDeepTau2017v2p1VSmuraw_2;
  std::vector<Float_t> *byVLooseDeepTau2017v2p1VSmu_2; 
  std::vector<Float_t> *byLooseDeepTau2017v2p1VSmu_2; 
  std::vector<Float_t> *byMediumDeepTau2017v2p1VSmu_2;   
  std::vector<Float_t> *byTightDeepTau2017v2p1VSmu_2;

  /* std::vector<Float_t> *pvx; */
  /* std::vector<Float_t> *pvy; */
  /* std::vector<Float_t> *pvz; */

  std::vector<Float_t> *dm_1;
  std::vector<Float_t> *dmMVA_1;

  std::vector<Float_t> *dm_2;
  std::vector<Float_t> *dmMVA_2;

  /* std::vector<Float_t> *svx_1; */
  /* std::vector<Float_t> *svy_1; */
  /* std::vector<Float_t> *svz_1; */

  /* std::vector<Float_t> *svx_2; */
  /* std::vector<Float_t> *svy_2; */
  /* std::vector<Float_t> *svz_2; */

  std::vector<Int_t> *tau1IndexVect;
  std::vector<Int_t> *tau2IndexVect;


  // List of branches
  TBranch        *b_EventNumber;   //!
  TBranch        *b_RunNumber;   //!
  TBranch        *b_lumi;   //!
  TBranch        *b_year;
  //TBranch        *b_NBadMu;   //!
  TBranch        *b_passecalBadCalibFilterUpdate;
  TBranch        *b_prefiringweight;
  TBranch        *b_prefiringweightup;
  TBranch        *b_prefiringweightdown;
     
  TBranch        *b_triggerbit;   //!
  TBranch        *b_metfilterbit;   //!
  //TBranch        *b_met;   //!
  TBranch        *b_met_er;
  TBranch        *b_met_er_phi;
  //TBranch        *b_metphi;   //!
  //TBranch        *b_PUPPImet;   //!
  //TBranch        *b_PUPPImetphi;   //!
  TBranch        *b_daughters_IetaIeta;   //!
  TBranch        *b_daughters_full5x5_IetaIeta;   //!
  TBranch        *b_daughters_hOverE;   //!
  TBranch        *b_daughters_deltaEtaSuperClusterTrackAtVtx;   //!
  TBranch        *b_daughters_deltaPhiSuperClusterTrackAtVtx;   //!
  TBranch        *b_daughters_IoEmIoP;   //!
  TBranch        *b_daughters_IoEmIoP_ttH;   //!
  //TBranch        *b_PFMETCov00;   //!
  //TBranch        *b_PFMETCov01;   //!
  //TBranch        *b_PFMETCov10;   //!
  //TBranch        *b_PFMETCov11;   //!
  TBranch        *b_PFMETsignif;   //!
  //TBranch        *b_npv;   //!
  //TBranch        *b_npu;   //!
  //TBranch        *b_PUReweight;   //!
  //TBranch        *b_rho;   //!
  TBranch        *b_mothers_px;   //!
  TBranch        *b_mothers_py;   //!
  TBranch        *b_mothers_pz;   //!
  TBranch        *b_mothers_e;   //!
  TBranch        *b_mothers_trgSeparateMatch;   //!
  TBranch        *b_trigger_name;   //!
  TBranch        *b_trigger_accept;   //!
  TBranch        *b_daughters_px;   //!
  TBranch        *b_daughters_py;   //!
  TBranch        *b_daughters_pz;   //!
  TBranch        *b_daughters_e;   //!
  TBranch        *b_daughters_charge;   //!
  TBranch        *b_daughters_charged_px;   //!
  TBranch        *b_daughters_charged_py;   //!
  TBranch        *b_daughters_charged_pz;   //!
  TBranch        *b_daughters_charged_e;   //!
  TBranch        *b_daughters_neutral_px;   //!
  TBranch        *b_daughters_neutral_py;   //!
  TBranch        *b_daughters_neutral_pz;   //!
  TBranch        *b_daughters_neutral_e;   //!
  TBranch        *b_daughters_hasTES;   //!
  /* TBranch        *b_daughters_px_TauUp;   //! */
  /* TBranch        *b_daughters_py_TauUp;   //! */
  /* TBranch        *b_daughters_pz_TauUp;   //! */
  /* TBranch        *b_daughters_e_TauUp;   //! */
  TBranch        *b_daughters_hasEES;   //!
  /* TBranch        *b_daughters_px_TauDown;   //! */
  /* TBranch        *b_daughters_py_TauDown;   //! */
  /* TBranch        *b_daughters_pz_TauDown;   //! */
  /* TBranch        *b_daughters_e_TauDown;   //! */
  /* TBranch        *b_daughters_px_EleUp; */
  /* TBranch        *b_daughters_py_EleUp; */
  /* TBranch        *b_daughters_pz_EleUp; */
  /* TBranch        *b_daughters_e_EleUp; */
  /* TBranch        *b_daughters_px_EleDown; */
  /* TBranch        *b_daughters_py_EleDown; */
  /* TBranch        *b_daughters_pz_EleDown; */
  /* TBranch        *b_daughters_e_EleDown; */
  TBranch        *b_daughters_isTauMatched;
  TBranch        *b_PUNumInteractions;   //!
  TBranch        *b_daughters_genindex;   //!
  TBranch        *b_MC_weight;   //!
  TBranch        *b_MC_weight_scale_muF0p5;   //!
  TBranch        *b_MC_weight_scale_muF2;   //!
  TBranch        *b_MC_weight_scale_muR0p5;   //!
  TBranch        *b_MC_weight_scale_muR2;   //!
  TBranch        *b_nominal_wt;
  TBranch        *b_TheoreticalPSUnc;
  TBranch        *b_TheoreticalScaleUnc1005;
  TBranch        *b_TheoreticalScaleUnc1009;
  TBranch        *b_TheoreticalScaleUnc5;
  TBranch        *b_TheoreticalScaleUnc9;
  
  TBranch        *b_lheHt;   //!
  TBranch        *b_lheNOutPartons;   //!
  TBranch        *b_lheNOutB;   //!
  TBranch        *b_lheNOutC;   //!
  //  TBranch        *b_aMCatNLOweight;   //!
  TBranch        *b_genpart_px;   //!
  TBranch        *b_genpart_py;   //!
  TBranch        *b_genpart_pz;   //!
  TBranch        *b_genpart_e;   //!
  TBranch        *b_DataMC_Type_idx;   //!
  TBranch        *b_Event_isRealData;   //!
  TBranch        *b_genpart_pca_x;   //!
  TBranch        *b_genpart_pca_y;   //!
  TBranch        *b_genpart_pca_z;   //!
  TBranch        *b_genpart_pdg;   //!
  TBranch        *b_genpart_status;   //!
  TBranch        *b_genpart_HMothInd;   //!
  TBranch        *b_genpart_MSSMHMothInd;   //!
  TBranch        *b_genpart_TopMothInd;   //!
  TBranch        *b_genpart_TauMothInd;   //!
  TBranch        *b_genpart_ZMothInd;   //!
  TBranch        *b_genpart_WMothInd;   //!
  TBranch        *b_genpart_bMothInd;   //!
  TBranch        *b_genpart_HZDecayMode;   //!
  TBranch        *b_genpart_TopDecayMode;   //!
  TBranch        *b_genpart_WDecayMode;   //!
  TBranch        *b_genpart_TauGenDecayMode;   //!
  TBranch        *b_genpart_TauGenDetailedDecayMode;   //!
  TBranch        *b_genpart_flags;   //!
  TBranch        *b_genjet_px;   //!
  TBranch        *b_genjet_py;   //!
  TBranch        *b_genjet_pz;   //!
  TBranch        *b_genjet_e;   //!
  TBranch        *b_genjet_partonFlavour;   //!
  TBranch        *b_genjet_hadronFlavour;   //!
  TBranch        *b_NUP;   //!
  /* TBranch        *b_SVfit_fitMETPhiTauUp;   //! */
  /* TBranch        *b_SVfit_fitMETPhiTauDown;   //! */
  /* TBranch        *b_SVfit_fitMETRhoTauUp;   //! */
  /* TBranch        *b_SVfit_fitMETRhoTauDown;   //! */
  /* TBranch        *b_SVfit_phiUncTauUp;   //! */
  /* TBranch        *b_SVfit_phiUncTauDown;   //! */
  /* TBranch        *b_SVfit_phiTauUp;   //! */
  /* TBranch        *b_SVfit_phiTauDown;   //! */
  /* TBranch        *b_SVfit_etaUncTauUp;   //! */
  /* TBranch        *b_SVfit_etaUncTauDown;   //! */
  /* TBranch        *b_SVfit_etaTauUp;   //! */
  /* TBranch        *b_SVfit_etaTauDown;   //! */
  /* TBranch        *b_SVfit_ptUncTauUp;   //! */
  /* TBranch        *b_SVfit_ptUncTauDown;   //! */
  /* TBranch        *b_SVfit_ptTauUp;   //! */
  /* TBranch        *b_SVfit_ptTauDown;   //! */
  /* TBranch        *b_SVfitTransverseMassTauUp;   //! */
  /* TBranch        *b_SVfitTransverseMassTauDown;   //! */
  /* TBranch        *b_SVfitMassTauUp;   //! */
  /* TBranch        *b_SVfitMassTauDown;   //! */
  TBranch        *b_MCSignalParticle_p4;   //!
  TBranch        *b_MCSignalParticle_pdgid;   //!
  TBranch        *b_MCSignalParticle_charge;   //!
  TBranch        *b_MCSignalParticle_Poca;   //!
  TBranch        *b_MCSignalParticle_Tauidx;   //!
  TBranch        *b_MCTauandProd_p4;   //!
  TBranch        *b_MCTauandProd_Vertex;   //!
  TBranch        *b_MCTauandProd_pdgid;   //!
  TBranch        *b_MCTauandProd_midx;   //!
  TBranch        *b_MCTauandProd_charge;   //!
  TBranch        *b_MCTau_JAK;   //!
  TBranch        *b_MCTau_DecayBitMask;   //!
  TBranch        *b_MC_p4;   //!
  TBranch        *b_MC_pdgid;   //!
  TBranch        *b_MC_charge;   //!
  TBranch        *b_MC_midx;   //!
  TBranch        *b_MC_childpdgid;   //!
  TBranch        *b_MC_childidx;   //!
  TBranch        *b_MC_status;   //!
  /* TBranch        *b_SVfitMass;   //! */
  /* TBranch        *b_SVfitTransverseMass;   //! */
  /* TBranch        *b_SVfit_pt;   //! */
  /* TBranch        *b_SVfit_ptUnc;   //! */
  /* TBranch        *b_SVfit_eta;   //! */
  /* TBranch        *b_SVfit_etaUnc;   //! */
  /* TBranch        *b_SVfit_phi;   //! */
  /* TBranch        *b_SVfit_phiUnc;   //! */
  /* TBranch        *b_SVfit_fitMETRho;   //! */
  /* TBranch        *b_SVfit_fitMETPhi;   //! */
  //  TBranch        *b_isOSCand;   //!
  TBranch        *b_METx;   //!
  TBranch        *b_METy;   //!
  /* TBranch        *b_metx_up_jes; */
  /* TBranch        *b_mety_up_jes; */
  /* TBranch        *b_metx_down_jes; */
  /* TBranch        *b_mety_down_jes; */
  /* TBranch        *b_metx_up_tes; */
  /* TBranch        *b_mety_up_tes; */
  /* TBranch        *b_metx_down_tes; */
  /* TBranch        *b_mety_down_tes; */
  /* TBranch        *b_metx_up_ees; */
  /* TBranch        *b_mety_up_ees; */
  /* TBranch        *b_metx_down_ees; */
  /* TBranch        *b_mety_down_ees; */
  /* //TBranch        *b_uncorrMETx;   //! */
  /* //TBranch        *b_uncorrMETy;   //! */
  /* TBranch        *b_MET_cov00;   //! */
  /* TBranch        *b_MET_cov01;   //! */
  /* TBranch        *b_MET_cov10;   //! */
  /* TBranch        *b_MET_cov11;   //! */
  /* TBranch        *b_MET_significance;   //! */
  /* TBranch        *b_mT_Dau1;   //! */
  /* TBranch        *b_mT_Dau2;   //! */
  TBranch        *b_PDGIdDaughters;   //!
  TBranch        *b_indexDau1;   //!
  TBranch        *b_indexDau2;   //!
  TBranch        *b_particleType;   //!
  //TBranch        *b_discriminator;   //!
  TBranch        *b_daughters_muonID;   //!
  TBranch        *b_daughters_typeOfMuon;   //!
  TBranch        *b_dxy;   //!
  TBranch        *b_dz;   //!
  TBranch        *b_dxy_innerTrack;   //!
  //  TBranch        *b_dz_innerTrack;   //!
  TBranch        *b_daughters_rel_error_trackpt;   //!
  TBranch        *b_SIP;   //!
  //TBranch        *b_daughters_iseleBDT;   //!
  TBranch        *b_daughters_iseleWPLoose;
  TBranch        *b_daughters_iseleWP80;   //!
  TBranch        *b_daughters_iseleWP90;   //!
  TBranch        *b_daughters_iseleNoIsoWPLoose;
  TBranch        *b_daughters_iseleNoIsoWP80;
  TBranch        *b_daughters_iseleNoIsoWP90;
  TBranch        *b_daughters_eleMVAnt;   //!
  TBranch        *b_daughters_eleMVA_HZZ;   //!
  TBranch        *b_daughters_passConversionVeto;   //!
  TBranch        *b_daughters_eleMissingHits;   //!
  TBranch        *b_daughters_iseleChargeConsistent;   //!
  //TBranch        *b_daughters_eleCUTID;   //!
  TBranch        *b_decayMode;   //!
  TBranch        *b_genmatch;
  TBranch        *b_tauID;   //!
  //  TBranch        *b_MVADM2016;   //!
  TBranch        *b_MVADM2017;   //!
  TBranch        *b_combreliso;   //!
  TBranch        *b_combreliso03;   //!
  TBranch        *b_daughters_depositR03_tracker;   //!
  TBranch        *b_daughters_depositR03_ecal;   //!
  TBranch        *b_daughters_depositR03_hcal;   //!
  TBranch        *b_daughters_decayModeFindingOldDMs;   //!
  //TBranch        *b_daughters_SCeta;   //!
  //TBranch        *b_againstElectronMVA5category;   //!
  //TBranch        *b_againstElectronMVA5raw;   //!
  //TBranch        *b_byPileupWeightedIsolationRaw3Hits;   //!
  TBranch        *b_footprintCorrection;   //!
  TBranch        *b_neutralIsoPtSumWeight;   //!
  TBranch        *b_photonPtSumOutsideSignalCone;   //!
  TBranch        *b_daughters_decayModeFindingNewDMs;   //!
  //TBranch        *b_daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits;   //!
  //TBranch        *b_daughters_byIsolationMVArun2017v2DBoldDMwLTraw2017;
  //  TBranch        *b_daughters_byIsolationMVArun2017v1DBoldDMwLTraw2017;
  //TBranch        *b_daughters_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017;
  //TBranch        *b_daughters_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017;
  //TBranch        *b_daughters_byIsolationMVA3oldDMwoLTraw;   //!
  //TBranch        *b_daughters_byIsolationMVA3oldDMwLTraw;   //!
  TBranch        *b_daughters_byIsolationMVA3newDMwoLTraw;   //!
  TBranch        *b_daughters_byIsolationMVA3newDMwLTraw;   //!
  //  TBranch        *b_daughters_byIsolationMVArun2v1DBoldDMwLTraw;   //!
  TBranch        *b_daughters_byDeepTau2017v2p1VSjetraw;   //!
  TBranch        *b_daughters_byDeepTau2017v2p1VSeraw;   //!
  TBranch        *b_daughters_byDeepTau2017v2p1VSmuraw;   //!
  TBranch        *b_daughters_chargedIsoPtSum;   //!
  TBranch        *b_daughters_neutralIsoPtSum;   //!
  TBranch        *b_daughters_puCorrPtSum;   //!
  TBranch        *b_daughters_numChargedParticlesSignalCone;   //!
  TBranch        *b_daughters_numNeutralHadronsSignalCone;   //!
  TBranch        *b_daughters_numPhotonsSignalCone;   //!
  TBranch        *b_daughters_daughters_numParticlesSignalCone;   //!
  TBranch        *b_daughters_numChargedParticlesIsoCone;   //!
  TBranch        *b_daughters_numNeutralHadronsIsoCone;   //!
  TBranch        *b_daughters_numPhotonsIsoCone;   //!
  TBranch        *b_daughters_numParticlesIsoCone;   //!
  TBranch        *b_daughters_leadChargedParticlePt;   //!
  TBranch        *b_daughters_trackRefPt;   //!
  //TBranch        *b_daughters_isLastTriggerObjectforPath;   //!
  TBranch        *b_daughters_trgMatched;   //!
  //TBranch        *b_daughters_isTriggerObjectforPath;   //!
  TBranch        *b_daughters_FilterFired;   //!
  TBranch        *b_daughters_isGoodTriggerType;   //!
  TBranch        *b_daughters_L3FilterFired;   //!
  TBranch        *b_daughters_L3FilterFiredLast;   //!
  TBranch        *b_daughters_HLTpt;   //!
  //TBranch        *b_daughters_isL1IsoTau28Matched;   //!
  TBranch        *b_Muon_trackCharge;   //!
  TBranch        *b_Muon_pdgid;   //!
  TBranch        *b_Muon_B;   //!
  TBranch        *b_Muon_M;   //!
  TBranch        *b_Muon_par;   //!
  //  TBranch        *b_Muon_cov;   //!
  /* TBranch        *b_PFTau_Track_par; */
  /* TBranch        *b_PFTau_Track_cov; */
  /* TBranch        *b_PFTau_Track_charge; */
  /* TBranch        *b_PFTau_Track_pdgid; */
  /* TBranch        *b_PFTau_Track_B; */
  /* TBranch        *b_PFTau_Track_M; */
  TBranch        *b_PFTauSVPos;   //!
  TBranch        *b_PFTauSVCov;   //!
  TBranch        *b_PFTauPionsP4;   //!
  TBranch        *b_PFTauRefitPionsP4;   //!
  TBranch        *b_PFTauRefitPionsCharge;   //!
  TBranch        *b_PFTauPionsCharge;   //!
  TBranch        *b_PFTauSVChi2NDofMatchingQuality;   //!
  TBranch        *b_TauFLSignificance;
  TBranch        *b_PFTauGEOMFlightLenght;
  TBranch        *b_PFTauGEOMFlightLenghtSignificance;
  TBranch        *b_daughters_jetNDauChargedMVASel;   //!
  TBranch        *b_daughters_miniRelIsoCharged;   //!
  TBranch        *b_daughters_miniRelIsoNeutral;   //!
  TBranch        *b_daughters_jetPtRel;   //!
  TBranch        *b_daughters_jetPtRatio;   //!
  TBranch        *b_daughters_jetBTagCSV;   //!
  TBranch        *b_daughters_jetBTagDeepCSV;
  TBranch        *b_daughters_jetBTagDeepFlavor;
  //TBranch        *b_daughters_lepMVA_mvaId;   //!
  TBranch        *b_daughters_pca_x;   //!
  TBranch        *b_daughters_pca_y;   //!
  TBranch        *b_daughters_pca_z;   //!
  TBranch        *b_daughters_pcaRefitPV_x;   //!
  TBranch        *b_daughters_pcaRefitPV_y;   //!
  TBranch        *b_daughters_pcaRefitPV_z;   //!
  TBranch        *b_daughters_pcaGenPV_x;   //!
  TBranch        *b_daughters_pcaGenPV_y;   //!
  TBranch        *b_daughters_pcaGenPV_z;   //!
  TBranch        *b_PFTau_a1_lvp;   //!
  //  TBranch        *b_PFTau_a1_cov;   //!
  TBranch        *b_PFTau_a1_charge;   //!
  TBranch        *b_PFTau_a1_pdgid;   //!
  TBranch        *b_PFTau_a1_B;   //!
  TBranch        *b_PFTau_a1_M;   //!
  TBranch        *b_PFTauTrack_deltaR;
  TBranch        *b_PFTauLeadTrackLV;
  //TBranch        *b_JetsNumber;   //!
  TBranch        *b_jets_VBFleadFilterMatch;
  TBranch        *b_jets_VBFsubleadFilterMatch;
  TBranch        *b_jets_px;   //!
  TBranch        *b_jets_py;   //!
  TBranch        *b_jets_pz;   //!
  TBranch        *b_jets_e;   //!
  TBranch        *b_jetsDown_px;
  TBranch        *b_jetsDown_py;
  TBranch        *b_jetsDown_pz;
  TBranch        *b_jetsDown_e;
  TBranch        *b_jetsUp_px;
  TBranch        *b_jetsUp_py;
  TBranch        *b_jetsUp_pz;
  TBranch        *b_jetsUp_e;
  TBranch        *b_jetRawf;   //!
  TBranch        *b_jets_area;   //!
  TBranch        *b_jets_JER;
  //  TBranch        *b_jets_QGdiscr;
  TBranch        *b_jets_mT;   //!
  TBranch        *b_jets_Flavour;   //!
  TBranch        *b_jets_HadronFlavour;   //!
  TBranch        *b_jets_genjetIndex;   //!
  TBranch        *b_jets_PUJetID;   //!
  TBranch        *b_jets_PUJetIDupdated;   //!
  TBranch        *b_jets_PUJetID_WP;
  TBranch        *b_jets_PUJetIDupdated_WP;
  TBranch        *b_jets_vtxPt;   //!
  TBranch        *b_jets_vtxMass;   //!
  TBranch        *b_jets_vtx3dL;   //!
  TBranch        *b_jets_vtxNtrk;   //!
  TBranch        *b_jets_vtx3deL;   //!
  TBranch        *b_jets_leadTrackPt;   //!
  TBranch        *b_jets_leptonPtRel;   //!
  TBranch        *b_jets_leptonPt;   //!
  TBranch        *b_jets_leptonDeltaR;   //!
  TBranch        *b_jets_chEmEF;   //!
  TBranch        *b_jets_chHEF;   //!
  TBranch        *b_jets_nEmEF;   //!
  TBranch        *b_jets_nHEF;   //!
  TBranch        *b_jets_MUF;   //!
  TBranch        *b_jets_neMult;   //!
  TBranch        *b_jets_chMult;   //!
  TBranch        *b_jets_jecUnc;   //!
  TBranch   *b_jets_jetUnc_Absolute_up;
  TBranch   *b_jets_jetUnc_FlavorQCD_up;
  TBranch   *b_jets_jetUnc_RelativeBal_up;
  TBranch   *b_jets_jetUnc_HF_up;
  TBranch   *b_jets_jetUnc_BBEC1_up;
  TBranch   *b_jets_jetUnc_EC2_up;
  TBranch   *b_jets_jetUnc_BBEC1_YEAR_up;
  TBranch   *b_jets_jetUnc_EC2_YEAR_up;
  TBranch   *b_jets_jetUnc_Absolute_YEAR_up;
  TBranch   *b_jets_jetUnc_HF_YEAR_up;
  TBranch   *b_jets_jetUnc_RelativeSample_YEAR_up;
  TBranch   *b_jets_jetUnc_Absolute_dw;
  TBranch   *b_jets_jetUnc_FlavorQCD_dw;
  TBranch   *b_jets_jetUnc_RelativeBal_dw;
  TBranch   *b_jets_jetUnc_HF_dw;
  TBranch   *b_jets_jetUnc_BBEC1_dw;
  TBranch   *b_jets_jetUnc_EC2_dw;
  TBranch   *b_jets_jetUnc_BBEC1_YEAR_dw;
  TBranch   *b_jets_jetUnc_EC2_YEAR_dw;
  TBranch   *b_jets_jetUnc_Absolute_YEAR_dw;
  TBranch   *b_jets_jetUnc_HF_YEAR_dw;
  TBranch   *b_jets_jetUnc_RelativeSample_YEAR_dw;
  TBranch        *b_PFJet_chargedMultiplicity;   //!
  TBranch        *b_PFJet_neutralMultiplicity;   //!
  TBranch        *b_PFJet_chargedEmEnergyFraction;   //!
  TBranch        *b_PFJet_chargedHadronEnergyFraction;   //!
  TBranch        *b_PFJet_neutralHadronEnergyFraction;   //!
  TBranch        *b_PFJet_neutralEmEnergyFraction;   //!
  TBranch        *b_bDiscriminator;   //!
  TBranch        *b_bCSVscore;   //!
  TBranch        *b_pfCombinedMVAV2BJetTags;   //!
  TBranch        *b_bDeepCSV_probb; //DeepCSV_probb
  TBranch        *b_bDeepCSV_probbb; //DeepCSV_probbb
  TBranch        *b_bDeepCSV_probudsg; //DeepCSV_probudsg
  TBranch        *b_bDeepCSV_probc; //DeepCSV_probc
  TBranch        *b_bDeepCSV_probcc; //DeepCSV_probcc
  TBranch        *b_bDeepFlavor_probb;  //DeepFlavor_probb
  TBranch        *b_bDeepFlavor_probbb; //DeepFlavor_probbb
  TBranch        *b_bDeepFlavor_problepb; //DeepFlavor_problepb
  TBranch        *b_bDeepFlavor_probc; //DeepFlavor_probc
  TBranch        *b_bDeepFlavor_probuds; //DeepFlavor_probuds
  TBranch        *b_bDeepFlavor_probg; //DeepFlavor_probg
  TBranch        *b_looseJetID;
  TBranch        *b_tightJetID;
  TBranch        *b_tightLepVetoJetID;
  //TBranch        *b_PFjetID;   //!
  TBranch        *b_ak8jets_px;   //!
  TBranch        *b_ak8jets_py;   //!
  TBranch        *b_ak8jets_pz;   //!
  TBranch        *b_ak8jets_e;   //!
  TBranch        *b_ak8jets_SoftDropMass;   //!
  /* TBranch        *b_ak8jets_PrunedMass;   //! */
  /* TBranch        *b_ak8jets_TrimmedMass;   //! */
  /* TBranch        *b_ak8jets_FilteredMass;   //! */
  TBranch        *b_ak8jets_tau1;   //!
  TBranch        *b_ak8jets_tau2;   //!
  TBranch        *b_ak8jets_tau3;   //!
  TBranch        *b_ak8jets_tau4;   //!
  TBranch        *b_ak8jets_CSV;   //!
  TBranch        *b_ak8jets_deepCSV_probb; // CSV score
  TBranch        *b_ak8jets_deepCSV_probbb; // CSV score
  /* TBranch        *b_ak8jets_deepFlavor_probb; // Flavor score */
  /* TBranch        *b_ak8jets_deepFlavor_probbb; // Flavor score */
  /* TBranch        *b_ak8jets_deepFlavor_problepb; // Flavor score */
  TBranch        *b_ak8jets_nsubjets;   //!
  TBranch        *b_subjets_px;   //!
  TBranch        *b_subjets_py;   //!
  TBranch        *b_subjets_pz;   //!
  TBranch        *b_subjets_e;   //!
  TBranch        *b_subjets_CSV;   //!
  TBranch        *b_subjets_deepCSV_probb;
  TBranch        *b_subjets_deepCSV_probbb;
  /* TBranch        *b_subjets_deepFlavor_probb; */
  /* TBranch        *b_subjets_deepFlavor_probbb; */
  /* TBranch        *b_subjets_deepFlavor_problepb; */
  TBranch        *b_subjets_ak8MotherIdx;   //!
  TBranch        *b_pv_x;   //!
  TBranch        *b_pv_y;   //!
  TBranch        *b_pv_z;   //!
  //  TBranch        *b_pv_cov;   //!
  TBranch        *b_RefitPVNoBS_x;
  TBranch        *b_RefitPVNoBS_y;
  TBranch        *b_RefitPVNoBS_z;
  TBranch        *b_RefitPVBS_x;
  TBranch        *b_RefitPVBS_y;
  TBranch        *b_RefitPVBS_z;
  TBranch        *b_RefitPVNoBS_xError;
  TBranch        *b_RefitPVNoBS_yError;
  TBranch        *b_RefitPVNoBS_zError;
  TBranch        *b_RefitPVBS_xError;
  TBranch        *b_RefitPVBS_yError;
  TBranch        *b_RefitPVBS_zError;
  
  TBranch        *b_RefitPVWithTracksBS_x;
  TBranch        *b_RefitPVWithTracksBS_y;
  TBranch        *b_RefitPVWithTracksBS_z;
  TBranch        *b_RefitPVWithTracksBS_xError;
  TBranch        *b_RefitPVWithTracksBS_yError;
  TBranch        *b_RefitPVWithTracksBS_zError;
  
  TBranch        *b_pvGen_x;   //!
  TBranch        *b_pvGen_y;   //!
  TBranch        *b_pvGen_z;   //!
  TBranch        *b_isRefitPV;   //!
   
  TBranch  *b_LeptonHash;
  TBranch *b_VertexHashNoBS1;
  TBranch *b_VertexHashNoBS2;
  TBranch  *b_VertexHashBS1;
  TBranch  *b_VertexHashBS2;
  
  TBranch *b_SelectedPairs;

  TBranch  *b_eleveto;
  TBranch  *b_muonveto;

  TBranch  *b_trg_doubletau;
  

  TBranch  *b_puppimt_1;
  TBranch  *b_gen_match_1;
  
  TBranch  *b_trigweight_1;
  TBranch  *b_idisoweight_1;
  TBranch  *b_antieweight_1;
  TBranch  *b_antimuweight_1;

  TBranch  *b_puppimt_2;
  TBranch  *b_gen_match_2;

  TBranch  *b_trigweight_2;
  TBranch  *b_idisoweight_2;
  TBranch  *b_antieweight_2;
  TBranch  *b_antimuweight_2;
  
  /* TBranch  *b_pt_tt; */
  /* TBranch  *b_pt_vis; */
  TBranch  *b_mt_tot;
  //  TBranch  *b_m_vis;

  TBranch *b_met;
  TBranch *b_metphi;
  TBranch *b_PUPPImet;
  TBranch *b_puppimet_ex_UnclusteredEnUp;
  TBranch *b_puppimet_ey_UnclusteredEnUp;
  TBranch *b_puppimet_ex_UnclusteredEnDown;
  TBranch *b_puppimet_ey_UnclusteredEnDown;
  TBranch *b_PUPPImetphi;
  TBranch *b_PFMETCov00;
  TBranch *b_PFMETCov01;
  TBranch *b_PFMETCov10;
  TBranch *b_PFMETCov11;
  TBranch *b_PUPPIMETCov00;
  TBranch *b_PUPPIMETCov01;
  TBranch *b_PUPPIMETCov10;
  TBranch *b_PUPPIMETCov11;
  
  TBranch  *b_mjj;
  TBranch  *b_jdeta;
  TBranch  *b_mjjDown;
  TBranch  *b_jdetaDown;
  TBranch  *b_mjjUp;
  TBranch  *b_jdetaUp;
  
  TBranch  *b_njetingap;
  TBranch  *b_njetingap20;
  TBranch  *b_jdphi;
  TBranch  *b_dijetpt;
  TBranch  *b_dijetphi;
  //  TBranch  *b_ptvis;

  TBranch  *b_nbtag;
  TBranch  *b_njets;
  TBranch  *b_njetspt20;
  TBranch  *b_njetsUp;
  TBranch  *b_njetspt20Up;
  TBranch  *b_njetsDown;
  TBranch  *b_njetspt20Down;
  TBranch  *b_jpt_1;
  TBranch  *b_jptDown_1;
  TBranch  *b_jptUp_1;
  
  TBranch  *b_jeta_1;
  TBranch  *b_jphi_1;
  TBranch  *b_jcsv_1;
  TBranch  *b_jpt_2;
  TBranch  *b_jeta_2;
  TBranch  *b_jphi_2;
  TBranch  *b_jcsv_2;
  TBranch  *b_bpt_1;
  TBranch  *b_beta_1;
  TBranch  *b_bphi_1;
  TBranch  *b_bcsv_1;
  TBranch  *b_bpt_2;
  TBranch  *b_beta_2;
  TBranch  *b_bphi_2;
  TBranch  *b_bcsv_2;

  TBranch *b_puweight;

  TBranch  *b_weight;

  /* TBranch  *b_jpfid_1; */
  /* TBranch  *b_jpuid_1; */
  /* TBranch  *b_jpfid_2; */
  /* TBranch  *b_jpuid_2; */
  /* TBranch  *b_bpfid_1; */
  /* TBranch  *b_bpuid_1; */
  /* TBranch  *b_bpfid_2; */
  /* TBranch  *b_bpuid_2; */

  TBranch *b_npv;
  TBranch *b_npu;
  TBranch *b_rho;
  
  //  TBranch  *b_byDeepTau2017v2p1VSjetraw_1;
  TBranch  *b_byVVVLooseDeepTau2017v2p1VSjet_1;
  TBranch  *b_byVVLooseDeepTau2017v2p1VSjet_1; 
  TBranch  *b_byVLooseDeepTau2017v2p1VSjet_1;  
  TBranch  *b_byLooseDeepTau2017v2p1VSjet_1;   
  TBranch  *b_byMediumDeepTau2017v2p1VSjet_1;  
  TBranch  *b_byTightDeepTau2017v2p1VSjet_1;   
  TBranch  *b_byVTightDeepTau2017v2p1VSjet_1;  
  TBranch  *b_byVVTightDeepTau2017v2p1VSjet_1; 
  
  //  TBranch  *b_byDeepTau2017v2p1VSeleraw_1;
  TBranch  *b_byVVVLooseDeepTau2017v2p1VSe_1;  
  TBranch  *b_byVVLooseDeepTau2017v2p1VSe_1; 
  TBranch  *b_byVLooseDeepTau2017v2p1VSe_1;   
  TBranch  *b_byLooseDeepTau2017v2p1VSe_1;	
  TBranch  *b_byMediumDeepTau2017v2p1VSe_1;   
  TBranch  *b_byTightDeepTau2017v2p1VSe_1;	
  TBranch  *b_byVTightDeepTau2017v2p1VSe_1;   
  TBranch  *b_byVVTightDeepTau2017v2p1VSe_1;
  
  //  TBranch  *b_byDeepTau2017v2p1VSmuraw_1;
  TBranch  *b_byVLooseDeepTau2017v2p1VSmu_1; 
  TBranch  *b_byLooseDeepTau2017v2p1VSmu_1; 
  TBranch  *b_byMediumDeepTau2017v2p1VSmu_1;   
  TBranch  *b_byTightDeepTau2017v2p1VSmu_1;

  //  TBranch  *b_byDeepTau2017v2p1VSjetraw_2;
  TBranch  *b_byVVVLooseDeepTau2017v2p1VSjet_2;
  TBranch  *b_byVVLooseDeepTau2017v2p1VSjet_2; 
  TBranch  *b_byVLooseDeepTau2017v2p1VSjet_2;  
  TBranch  *b_byLooseDeepTau2017v2p1VSjet_2;   
  TBranch  *b_byMediumDeepTau2017v2p1VSjet_2;  
  TBranch  *b_byTightDeepTau2017v2p1VSjet_2;   
  TBranch  *b_byVTightDeepTau2017v2p1VSjet_2;  
  TBranch  *b_byVVTightDeepTau2017v2p1VSjet_2;
  
  //  TBranch  *b_byDeepTau2017v2p1VSeleraw_2;
  TBranch  *b_byVVVLooseDeepTau2017v2p1VSe_2;  
  TBranch  *b_byVVLooseDeepTau2017v2p1VSe_2; 
  TBranch  *b_byVLooseDeepTau2017v2p1VSe_2;   
  TBranch  *b_byLooseDeepTau2017v2p1VSe_2;	
  TBranch  *b_byMediumDeepTau2017v2p1VSe_2;   
  TBranch  *b_byTightDeepTau2017v2p1VSe_2;	
  TBranch  *b_byVTightDeepTau2017v2p1VSe_2;   
  TBranch  *b_byVVTightDeepTau2017v2p1VSe_2;
  
  //  TBranch  *b_byDeepTau2017v2p1VSmuraw_2;
  TBranch  *b_byVLooseDeepTau2017v2p1VSmu_2; 
  TBranch  *b_byLooseDeepTau2017v2p1VSmu_2; 
  TBranch  *b_byMediumDeepTau2017v2p1VSmu_2;   
  TBranch  *b_byTightDeepTau2017v2p1VSmu_2;

  /* TBranch  *b_pvx; */
  /* TBranch  *b_pvy; */
  /* TBranch  *b_pvz; */

  TBranch  *b_dm_1;
  TBranch  *b_dmMVA_1;

  TBranch  *b_dm_2;
  TBranch  *b_dmMVA_2;

  /* TBranch  *b_svx_1; */
  /* TBranch  *b_svy_1; */
  /* TBranch  *b_svz_1; */

  /* TBranch  *b_svx_2; */
  /* TBranch  *b_svy_2; */
  /* TBranch  *b_svz_2; */

  TBranch   *b_tau1IndexVect;
  TBranch   *b_tau2IndexVect;

  NtupleReader(TTree *tree=0);
  virtual ~NtupleReader();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef NtupleReader_cxx
NtupleReader::NtupleReader(TTree *tree) : fChain(0) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {


#ifdef SINGLE_TREE
    // The following code should be used if you want this class to access
    // a single tree instead of a chain
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
    if (!f || !f->IsOpen()) {
      f = new TFile("Memory Directory");
    }
    f->GetObject("HTauTauTree/HTauTauTree",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
    TChain * chain = new TChain("HTauTauTree/HTauTauTree","");
    chain->Add("/home-pbs/vcherepa/cms_work/CMSSW_8_0_25/src/LLRHiggsTauTau/NtupleProducer/test/HTauTauAnalysis.root/HTauTauTree/HTauTauTree");
    tree = chain;
#endif // SINGLE_TREE

  }
  Init(tree);
}



NtupleReader::~NtupleReader()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t NtupleReader::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t NtupleReader::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void NtupleReader::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set object pointer
  daughters_IetaIeta = 0;
  daughters_full5x5_IetaIeta = 0;
  daughters_hOverE = 0;
  daughters_deltaEtaSuperClusterTrackAtVtx = 0;
  daughters_deltaPhiSuperClusterTrackAtVtx = 0;
  daughters_IoEmIoP = 0;
  daughters_IoEmIoP_ttH = 0;
  mothers_px = 0;
  mothers_py = 0;
  mothers_pz = 0;
  mothers_e = 0;
  mothers_trgSeparateMatch = 0;
  trigger_name=0;
  trigger_accept = 0;
  daughters_px = 0;
  daughters_py = 0;
  daughters_pz = 0;
  daughters_e = 0;
  daughters_charge = 0;
  daughters_charged_px = 0;
  daughters_charged_py = 0;
  daughters_charged_pz = 0;
  daughters_charged_e = 0;
  daughters_neutral_px = 0;
  daughters_neutral_py = 0;
  daughters_neutral_pz = 0;
  daughters_neutral_e = 0;
  daughters_hasTES = 0;
  /* daughters_px_TauUp = 0; */
  /* daughters_py_TauUp = 0; */
  /* daughters_pz_TauUp = 0; */
  /* daughters_e_TauUp = 0; */
  daughters_hasEES = 0;
  /* daughters_px_TauDown = 0; */
  /* daughters_py_TauDown = 0; */
  /* daughters_pz_TauDown = 0; */
  /* daughters_e_TauDown = 0; */
  /* daughters_px_EleUp = 0; */
  /* daughters_py_EleUp = 0; */
  /* daughters_pz_EleUp = 0; */
  /* daughters_e_EleUp = 0; */
  /* daughters_px_EleDown = 0; */
  /* daughters_py_EleDown = 0; */
  /* daughters_pz_EleDown = 0; */
  /* daughters_e_EleDown = 0; */
  daughters_isTauMatched =0;
  daughters_genindex = 0;
  TheoreticalPSUnc=0;
  genpart_px = 0;
  genpart_py = 0;
  genpart_pz = 0;
  genpart_e = 0;
  genpart_pca_x = 0;
  genpart_pca_y = 0;
  genpart_pca_z = 0;
  genpart_pdg = 0;
  genpart_status = 0;
  genpart_HMothInd = 0;
  genpart_MSSMHMothInd = 0;
  genpart_TopMothInd = 0;
  genpart_TauMothInd = 0;
  genpart_ZMothInd = 0;
  genpart_WMothInd = 0;
  genpart_bMothInd = 0;
  genpart_HZDecayMode = 0;
  genpart_TopDecayMode = 0;
  genpart_WDecayMode = 0;
  genpart_TauGenDecayMode = 0;
  genpart_TauGenDetailedDecayMode = 0;
  genpart_flags = 0;
  genjet_px = 0;
  genjet_py = 0;
  genjet_pz = 0;
  genjet_e = 0;
  genjet_partonFlavour = 0;
  genjet_hadronFlavour = 0;
  /* SVfit_fitMETPhiTauUp = 0; */
  /* SVfit_fitMETPhiTauDown = 0; */
  /* SVfit_fitMETRhoTauUp = 0; */
  /* SVfit_fitMETRhoTauDown = 0; */
  /* SVfit_phiUncTauUp = 0; */
  /* SVfit_phiUncTauDown = 0; */
  /* SVfit_phiTauUp = 0; */
  /* SVfit_phiTauDown = 0; */
  /* SVfit_etaUncTauUp = 0; */
  /* SVfit_etaUncTauDown = 0; */
  /* SVfit_etaTauUp = 0; */
  /* SVfit_etaTauDown = 0; */
  /* SVfit_ptUncTauUp = 0; */
  /* SVfit_ptUncTauDown = 0; */
  /* SVfit_ptTauUp = 0; */
  /* SVfit_ptTauDown = 0; */
  /* SVfitTransverseMassTauUp = 0; */
  /* SVfitTransverseMassTauDown = 0; */
  /* SVfitMassTauUp = 0; */
  /* SVfitMassTauDown = 0; */
  MCSignalParticle_p4 = 0;
  MCSignalParticle_pdgid = 0;
  MCSignalParticle_charge = 0;
  MCSignalParticle_Poca = 0;
  MCSignalParticle_Tauidx = 0;
  MCTauandProd_p4 = 0;
  MCTauandProd_Vertex = 0;
  MCTauandProd_pdgid = 0;
  MCTauandProd_midx = 0;
  MCTauandProd_charge = 0;
  MCTau_JAK = 0;
  MCTau_DecayBitMask = 0;
  MC_p4 = 0;
  MC_pdgid = 0;
  MC_charge = 0;
  MC_midx = 0;
  MC_childpdgid = 0;
  MC_childidx = 0;
  MC_status = 0;
  /* SVfitMass = 0; */
  /* SVfitTransverseMass = 0; */
  /* SVfit_pt = 0; */
  /* SVfit_ptUnc = 0; */
  /* SVfit_eta = 0; */
  /* SVfit_etaUnc = 0; */
  /* SVfit_phi = 0; */
  /* SVfit_phiUnc = 0; */
  /* SVfit_fitMETRho = 0; */
  /* SVfit_fitMETPhi = 0; */
  //  isOSCand = 0;
  METx = 0;
  METy = 0;
  /* metx_up_jes=0; */
  /* mety_up_jes=0; */
  /* metx_down_jes=0; */
  /* mety_down_jes=0; */
  /* metx_up_tes=0; */
  /* mety_up_tes=0; */
  /* metx_down_tes=0; */
  /* mety_down_tes=0; */
  /* metx_up_ees=0; */
  /* mety_up_ees=0; */
  /* metx_down_ees=0; */
  /* mety_down_ees=0; */
  /* //uncorrMETx = 0; */
  /* //uncorrMETy = 0; */
  /* MET_cov00 = 0; */
  /* MET_cov01 = 0; */
  /* MET_cov10 = 0; */
  /* MET_cov11 = 0; */
  /* MET_significance = 0; */
  /* mT_Dau1 = 0; */
  /* mT_Dau2 = 0; */
  PDGIdDaughters = 0;
  indexDau1 = 0;
  indexDau2 = 0;
  particleType = 0;
  //discriminator = 0;
  daughters_muonID = 0;
  daughters_typeOfMuon = 0;
  dxy = 0;
  dz = 0;
  dxy_innerTrack = 0;
  //  dz_innerTrack = 0;
  daughters_rel_error_trackpt = 0;
  SIP = 0;
  //daughters_iseleBDT = 0;
  daughters_iseleWPLoose=0;
  daughters_iseleWP80 = 0;
  daughters_iseleWP90 = 0;
  daughters_iseleNoIsoWPLoose=0;
  daughters_iseleNoIsoWP80=0;
  daughters_iseleNoIsoWP90=0;
  daughters_eleMVAnt = 0;
  daughters_eleMVA_HZZ = 0;
  daughters_passConversionVeto = 0;
  daughters_eleMissingHits = 0;
  daughters_iseleChargeConsistent = 0;
  //daughters_eleCUTID = 0;
  decayMode = 0;
  genmatch =0;
  tauID = 0;
  //  MVADM2016=0;
  MVADM2017=0;
  combreliso = 0;
  combreliso03 = 0;
  daughters_depositR03_tracker = 0;
  daughters_depositR03_ecal = 0;
  daughters_depositR03_hcal = 0;
  daughters_decayModeFindingOldDMs = 0;
  //daughters_SCeta = 0;
  //againstElectronMVA5category = 0;
  //againstElectronMVA5raw = 0;
  //byPileupWeightedIsolationRaw3Hits = 0;
  footprintCorrection = 0;
  neutralIsoPtSumWeight = 0;
  photonPtSumOutsideSignalCone = 0;
  daughters_decayModeFindingNewDMs = 0;
  //daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits = 0;
  // daughters_byIsolationMVArun2017v2DBoldDMwLTraw2017 = 0;
  //  daughters_byIsolationMVArun2017v1DBoldDMwLTraw2017 = 0;
  //daughters_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017 = 0;
  //daughters_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017 = 0;
  //daughters_byIsolationMVA3oldDMwoLTraw = 0;
  //daughters_byIsolationMVA3oldDMwLTraw = 0;
  daughters_byIsolationMVA3newDMwoLTraw = 0;
  daughters_byIsolationMVA3newDMwLTraw = 0;
  //  daughters_byIsolationMVArun2v1DBoldDMwLTraw = 0;
  daughters_byDeepTau2017v2p1VSjetraw = 0;
  daughters_byDeepTau2017v2p1VSeraw = 0;
  daughters_byDeepTau2017v2p1VSmuraw = 0;
  daughters_chargedIsoPtSum = 0;
  daughters_neutralIsoPtSum = 0;
  daughters_puCorrPtSum = 0;
  daughters_numChargedParticlesSignalCone = 0;
  daughters_numNeutralHadronsSignalCone = 0;
  daughters_numPhotonsSignalCone = 0;
  daughters_daughters_numParticlesSignalCone = 0;
  daughters_numChargedParticlesIsoCone = 0;
  daughters_numNeutralHadronsIsoCone = 0;
  daughters_numPhotonsIsoCone = 0;
  daughters_numParticlesIsoCone = 0;
  daughters_leadChargedParticlePt = 0;
  daughters_trackRefPt = 0;
  //daughters_isLastTriggerObjectforPath = 0;
  daughters_trgMatched = 0;
  //daughters_isTriggerObjectforPath = 0;
  daughters_FilterFired = 0;
  daughters_isGoodTriggerType = 0;
  daughters_L3FilterFired = 0;
  daughters_L3FilterFiredLast = 0;
  daughters_HLTpt = 0;
  //daughters_isL1IsoTau28Matched = 0;
  Muon_trackCharge = 0;
  Muon_pdgid = 0;
  Muon_B = 0;
  Muon_M = 0;
  Muon_par = 0;
  //  Muon_cov = 0;
  /* PFTau_Track_par = 0; */
  /* PFTau_Track_cov = 0; */
  /* PFTau_Track_charge = 0; */
  /* PFTau_Track_pdgid = 0; */
  /* PFTau_Track_B = 0; */
  /* PFTau_Track_M = 0; */
  PFTauSVPos = 0;
  PFTauSVCov = 0;
  PFTauPionsP4 = 0;
  PFTauRefitPionsP4 = 0;
  //PFTau_TIP_PVPosNoBS = 0 ;
  //PFTau_TIP_PVPosBS = 0 ;
  PFTauPionsCharge = 0;
  PFTauRefitPionsCharge = 0;
  PFTauSVChi2NDofMatchingQuality = 0;
  TauFLSignificance = 0;
  PFTauGEOMFlightLenght = 0;
  PFTauGEOMFlightLenghtSignificance = 0;
  daughters_jetNDauChargedMVASel = 0;
  daughters_miniRelIsoCharged = 0;
  daughters_miniRelIsoNeutral = 0;
  daughters_jetPtRel = 0;
  daughters_jetPtRatio = 0;
  daughters_jetBTagCSV = 0;
  daughters_jetBTagDeepCSV = 0;
  daughters_jetBTagDeepFlavor = 0;
  //daughters_lepMVA_mvaId = 0;
  daughters_pca_x = 0;
  daughters_pca_y = 0;
  daughters_pca_z = 0;
  daughters_pcaRefitPV_x = 0;
  daughters_pcaRefitPV_y = 0;
  daughters_pcaRefitPV_z = 0;
  daughters_pcaGenPV_x = 0;
  daughters_pcaGenPV_y = 0;
  daughters_pcaGenPV_z = 0;
  PFTau_a1_lvp = 0;
  //  PFTau_a1_cov = 0;
  PFTau_a1_charge = 0;
  PFTau_a1_pdgid = 0;
  PFTau_a1_B = 0;
  PFTau_a1_M = 0;
  PFTauLeadTrackLV = 0;
  PFTauTrack_deltaR = 0;
  //JetsNumber = 0;
  jets_VBFleadFilterMatch = 0;
  jets_VBFsubleadFilterMatch = 0;
  jets_px = 0;
  jets_py = 0;
  jets_pz = 0;
  jets_e = 0;
  jetsDown_px = 0;
  jetsDown_py = 0;
  jetsDown_pz = 0;
  jetsDown_e = 0;
  jetsUp_px = 0;
  jetsUp_py = 0;
  jetsUp_pz = 0;
  jetsUp_e = 0;
  jetRawf = 0;
  jets_area = 0;
  jets_JER = 0;
  //  jets_QGdiscr = 0;
  jets_mT = 0;
  jets_Flavour = 0;
  jets_HadronFlavour = 0;
  jets_genjetIndex = 0;
  jets_PUJetID = 0;
  jets_PUJetIDupdated = 0;
  jets_PUJetID_WP = 0;
  jets_PUJetIDupdated_WP = 0;
  jets_vtxPt = 0;
  jets_vtxMass = 0;
  jets_vtx3dL = 0;
  jets_vtxNtrk = 0;
  jets_vtx3deL = 0;
  jets_leadTrackPt = 0;
  jets_leptonPtRel = 0;
  jets_leptonPt = 0;
  jets_leptonDeltaR = 0;
  _jets_chEmEF = 0;
  _jets_chHEF = 0;
  _jets_nEmEF = 0;
  _jets_nHEF = 0;
  jets_MUF = 0;
  _jets_neMult = 0;
  _jets_chMult = 0;
  jets_jecUnc = 0;
  jets_jetUnc_Absolute_up=0;
  jets_jetUnc_FlavorQCD_up=0;
  jets_jetUnc_RelativeBal_up=0;
  jets_jetUnc_HF_up=0;
  jets_jetUnc_BBEC1_up=0;
  jets_jetUnc_EC2_up=0;
  jets_jetUnc_BBEC1_YEAR_up=0;
  jets_jetUnc_EC2_YEAR_up=0;
  jets_jetUnc_Absolute_YEAR_up=0;
  jets_jetUnc_HF_YEAR_up=0;
  jets_jetUnc_RelativeSample_YEAR_up=0;
  jets_jetUnc_Absolute_dw=0;
  jets_jetUnc_FlavorQCD_dw=0;
  jets_jetUnc_RelativeBal_dw=0;
  jets_jetUnc_HF_dw=0;
  jets_jetUnc_BBEC1_dw=0;
  jets_jetUnc_EC2_dw=0;
  jets_jetUnc_BBEC1_YEAR_dw=0;
  jets_jetUnc_EC2_YEAR_dw=0;
  jets_jetUnc_Absolute_YEAR_dw=0;
  jets_jetUnc_HF_YEAR_dw=0;
  jets_jetUnc_RelativeSample_YEAR_dw=0;
  bDiscriminator = 0;
  bCSVscore = 0;
  pfCombinedMVAV2BJetTags = 0;
  bDeepCSV_probb=0; //DeepCSV_probb
  bDeepCSV_probbb=0; //DeepCSV_probbb
  bDeepCSV_probudsg=0; //DeepCSV_probudsg
  bDeepCSV_probc=0; //DeepCSV_probc
  bDeepCSV_probcc=0; //DeepCSV_probcc
  bDeepFlavor_probb=0;  //DeepFlavor_probb
  bDeepFlavor_probbb=0; //DeepFlavor_probbb
  bDeepFlavor_problepb=0; //DeepFlavor_problepb
  bDeepFlavor_probc=0; //DeepFlavor_probc
  bDeepFlavor_probuds=0; //DeepFlavor_probuds
  bDeepFlavor_probg=0; //DeepFlavor_probg
  looseJetID=0;
  tightJetID=0;
  tightLepVetoJetID=0;
  //PFjetID = 0;
  ak8jets_px = 0;
  ak8jets_py = 0;
  ak8jets_pz = 0;
  ak8jets_e = 0;
  ak8jets_SoftDropMass = 0;
  /* ak8jets_PrunedMass = 0; */
  /* ak8jets_TrimmedMass = 0; */
  /* ak8jets_FilteredMass = 0; */
  ak8jets_tau1 = 0;
  ak8jets_tau2 = 0;
  ak8jets_tau3 = 0;
  ak8jets_tau4 = 0;
  ak8jets_CSV = 0;
  ak8jets_deepCSV_probb = 0;
  ak8jets_deepCSV_probbb = 0;
  /* ak8jets_deepFlavor_probb = 0; */
  /* ak8jets_deepFlavor_probbb = 0; */
  /* ak8jets_deepFlavor_problepb = 0; */
  ak8jets_nsubjets = 0;
  subjets_px = 0;
  subjets_py = 0;
  subjets_pz = 0;
  subjets_e = 0;
  subjets_CSV = 0;
  subjets_deepCSV_probb = 0;
  subjets_deepCSV_probbb = 0;
  /* subjets_deepFlavor_probb = 0; */
  /* subjets_deepFlavor_probbb = 0; */
  /* subjets_deepFlavor_problepb = 0; */
  subjets_ak8MotherIdx = 0;
  //  pv_cov = 0;

  RefitPVNoBS_x = 0;
  RefitPVNoBS_y = 0;
  RefitPVNoBS_z = 0;
   
  RefitPVBS_x = 0;
  RefitPVBS_y = 0;
  RefitPVBS_z = 0;

  RefitPVNoBS_xError = 0;
  RefitPVNoBS_yError = 0;
  RefitPVNoBS_zError = 0;
   
  RefitPVBS_xError = 0;
  RefitPVBS_yError = 0;
  RefitPVBS_zError = 0;
  
  LeptonHash = 0;
  VertexHashNoBS1 = 0;
  VertexHashNoBS2 = 0;
  VertexHashBS1 = 0;
  VertexHashBS2 = 0;

  trg_doubletau=0;

  puppimt_1=0;
  gen_match_1=0;
  
  trigweight_1=0;
  idisoweight_1=0;
  antieweight_1=0;
  antimuweight_1=0;

  puppimt_2=0;
  gen_match_2=0;

  trigweight_2=0;
  idisoweight_2=0;
  antieweight_2=0;
  antimuweight_2=0;
  
  /* pt_tt=0; */
  /* pt_vis=0; */
  mt_tot=0;
  //  m_vis=0;
  
  mjj=0;
  jdeta=0;
  mjjDown=0;
  jdetaDown=0;
  mjjUp=0;
  jdetaUp=0;
  
  njetingap=0;
  njetingap20=0;
  jdphi=0;
  dijetpt=0;
  dijetphi=0;
  //  ptvis=0;

  nbtag=0;
  njets=0;
  njetspt20=0;
  njetsUp=0;
  njetspt20Up=0;
  njetsDown=0;
  njetspt20Down=0;
  jpt_1=0;
  jptDown_1=0;
  jptUp_1=0;
  
  jeta_1=0;
  jphi_1=0;
  jcsv_1=0;
  jpt_2=0;
  jeta_2=0;
  jphi_2=0;
  jcsv_2=0;
  bpt_1=0;
  beta_1=0;
  bphi_1=0;
  bcsv_1=0;
  bpt_2=0;
  beta_2=0;
  bphi_2=0;
  bcsv_2=0;

  weight=0;

  /* jpfid_1=0; */
  /* jpuid_1=0; */
  /* jpfid_2=0; */
  /* jpuid_2=0; */
  /* bpfid_1=0; */
  /* bpuid_1=0; */
  /* bpfid_2=0; */
  /* bpuid_2=0; */
  
  //  byDeepTau2017v2p1VSjetraw_1=0;
  byVVVLooseDeepTau2017v2p1VSjet_1=0;
  byVVLooseDeepTau2017v2p1VSjet_1=0; 
  byVLooseDeepTau2017v2p1VSjet_1=0;  
  byLooseDeepTau2017v2p1VSjet_1=0;   
  byMediumDeepTau2017v2p1VSjet_1=0;  
  byTightDeepTau2017v2p1VSjet_1=0;   
  byVTightDeepTau2017v2p1VSjet_1=0;  
  byVVTightDeepTau2017v2p1VSjet_1=0; 
  
  //  byDeepTau2017v2p1VSeleraw_1=0;
  byVVVLooseDeepTau2017v2p1VSe_1=0;  
  byVVLooseDeepTau2017v2p1VSe_1=0; 
  byVLooseDeepTau2017v2p1VSe_1=0;   
  byLooseDeepTau2017v2p1VSe_1=0;	
  byMediumDeepTau2017v2p1VSe_1=0;   
  byTightDeepTau2017v2p1VSe_1=0;	
  byVTightDeepTau2017v2p1VSe_1=0;   
  byVVTightDeepTau2017v2p1VSe_1=0;
  
  //  byDeepTau2017v2p1VSmuraw_1=0;
  byVLooseDeepTau2017v2p1VSmu_1=0; 
  byLooseDeepTau2017v2p1VSmu_1=0; 
  byMediumDeepTau2017v2p1VSmu_1=0;   
  byTightDeepTau2017v2p1VSmu_1=0;

  //  byDeepTau2017v2p1VSjetraw_2=0;
  byVVVLooseDeepTau2017v2p1VSjet_2=0;
  byVVLooseDeepTau2017v2p1VSjet_2=0; 
  byVLooseDeepTau2017v2p1VSjet_2=0;  
  byLooseDeepTau2017v2p1VSjet_2=0;   
  byMediumDeepTau2017v2p1VSjet_2=0;  
  byTightDeepTau2017v2p1VSjet_2=0;   
  byVTightDeepTau2017v2p1VSjet_2=0;  
  byVVTightDeepTau2017v2p1VSjet_2=0;
  
  //  byDeepTau2017v2p1VSeleraw_2=0;
  byVVVLooseDeepTau2017v2p1VSe_2=0;  
  byVVLooseDeepTau2017v2p1VSe_2=0; 
  byVLooseDeepTau2017v2p1VSe_2=0;   
  byLooseDeepTau2017v2p1VSe_2=0;	
  byMediumDeepTau2017v2p1VSe_2=0;   
  byTightDeepTau2017v2p1VSe_2=0;	
  byVTightDeepTau2017v2p1VSe_2=0;   
  byVVTightDeepTau2017v2p1VSe_2=0;
  
  //  byDeepTau2017v2p1VSmuraw_2=0;
  byVLooseDeepTau2017v2p1VSmu_2=0; 
  byLooseDeepTau2017v2p1VSmu_2=0; 
  byMediumDeepTau2017v2p1VSmu_2=0;   
  byTightDeepTau2017v2p1VSmu_2=0;

  //pvx=0;
  //pvy=0;
  //pvz=0;

  dm_1=0;
  dmMVA_1=0;

  dm_2=0;
  dmMVA_2=0;

  /* svx_1=0; */
  /* svy_1=0; */
  /* svz_1=0; */

  /* svx_2=0; */
  /* svy_2=0; */
  /* svz_2=0; */

  tau1IndexVect=0;
  tau2IndexVect=0;
   
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("evt", &EventNumber, &b_EventNumber);
  fChain->SetBranchAddress("run", &RunNumber, &b_RunNumber);
  fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
  fChain->SetBranchAddress("year", &year, &b_year);
  //fChain->SetBranchAddress("NBadMu", &NBadMu, &b_NBadMu);
  fChain->SetBranchAddress("passecalBadCalibFilterUpdate",&passecalBadCalibFilterUpdate,&b_passecalBadCalibFilterUpdate);
  fChain->SetBranchAddress("prefiringweight",&prefiringweight,&b_prefiringweight);
  fChain->SetBranchAddress("prefiringweightup",&prefiringweightup,&b_prefiringweightup);
  fChain->SetBranchAddress("prefiringweightdown",&prefiringweightdown,&b_prefiringweightdown);
  fChain->SetBranchAddress("triggerbit", &triggerbit, &b_triggerbit);
  fChain->SetBranchAddress("metfilterbit", &metfilterbit, &b_metfilterbit);
  fChain->SetBranchAddress("met", &met, &b_met);
  fChain->SetBranchAddress("met_er",&met_er,&b_met_er);
  fChain->SetBranchAddress("met_er_phi",&met_er_phi,&b_met_er_phi);
  fChain->SetBranchAddress("metphi", &metphi, &b_metphi);
  //fChain->SetBranchAddress("PUPPImet", &PUPPImet, &b_PUPPImet);
  //fChain->SetBranchAddress("PUPPImetphi", &PUPPImetphi, &b_PUPPImetphi);
  fChain->SetBranchAddress("daughters_IetaIeta", &daughters_IetaIeta, &b_daughters_IetaIeta);
  fChain->SetBranchAddress("daughters_full5x5_IetaIeta", &daughters_full5x5_IetaIeta, &b_daughters_full5x5_IetaIeta);
  fChain->SetBranchAddress("daughters_hOverE", &daughters_hOverE, &b_daughters_hOverE);
  fChain->SetBranchAddress("daughters_deltaEtaSuperClusterTrackAtVtx", &daughters_deltaEtaSuperClusterTrackAtVtx, &b_daughters_deltaEtaSuperClusterTrackAtVtx);
  fChain->SetBranchAddress("daughters_deltaPhiSuperClusterTrackAtVtx", &daughters_deltaPhiSuperClusterTrackAtVtx, &b_daughters_deltaPhiSuperClusterTrackAtVtx);
  fChain->SetBranchAddress("daughters_IoEmIoP", &daughters_IoEmIoP, &b_daughters_IoEmIoP);
  fChain->SetBranchAddress("daughters_IoEmIoP_ttH", &daughters_IoEmIoP_ttH, &b_daughters_IoEmIoP_ttH);
  /* fChain->SetBranchAddress("PFMETCov00", &PFMETCov00, &b_PFMETCov00); */
  /* fChain->SetBranchAddress("PFMETCov01", &PFMETCov01, &b_PFMETCov01); */
  /* fChain->SetBranchAddress("PFMETCov10", &PFMETCov10, &b_PFMETCov10); */
  /* fChain->SetBranchAddress("PFMETCov11", &PFMETCov11, &b_PFMETCov11); */
  fChain->SetBranchAddress("PFMETsignif", &PFMETsignif, &b_PFMETsignif);
  fChain->SetBranchAddress("npv", &npv, &b_npv);
  fChain->SetBranchAddress("npu", &npu, &b_npu);
  //fChain->SetBranchAddress("PUReweight", &PUReweight, &b_PUReweight);
  fChain->SetBranchAddress("rho", &rho, &b_rho);
  fChain->SetBranchAddress("mothers_px", &mothers_px, &b_mothers_px);
  fChain->SetBranchAddress("mothers_py", &mothers_py, &b_mothers_py);
  fChain->SetBranchAddress("mothers_pz", &mothers_pz, &b_mothers_pz);
  fChain->SetBranchAddress("mothers_e", &mothers_e, &b_mothers_e);
  fChain->SetBranchAddress("mothers_trgSeparateMatch", &mothers_trgSeparateMatch, &b_mothers_trgSeparateMatch);
  fChain->SetBranchAddress("trigger_name", &trigger_name, &b_trigger_name);
  fChain->SetBranchAddress("trigger_accept", &trigger_accept, &b_trigger_accept);
  fChain->SetBranchAddress("daughters_px", &daughters_px, &b_daughters_px);
  fChain->SetBranchAddress("daughters_py", &daughters_py, &b_daughters_py);
  fChain->SetBranchAddress("daughters_pz", &daughters_pz, &b_daughters_pz);
  fChain->SetBranchAddress("daughters_e", &daughters_e, &b_daughters_e);
  fChain->SetBranchAddress("daughters_charge", &daughters_charge, &b_daughters_charge);
  fChain->SetBranchAddress("daughters_charged_px", &daughters_charged_px, &b_daughters_charged_px);
  fChain->SetBranchAddress("daughters_charged_py", &daughters_charged_py, &b_daughters_charged_py);
  fChain->SetBranchAddress("daughters_charged_pz", &daughters_charged_pz, &b_daughters_charged_pz);
  fChain->SetBranchAddress("daughters_charged_e", &daughters_charged_e, &b_daughters_charged_e);
  fChain->SetBranchAddress("daughters_neutral_px", &daughters_neutral_px, &b_daughters_neutral_px);
  fChain->SetBranchAddress("daughters_neutral_py", &daughters_neutral_py, &b_daughters_neutral_py);
  fChain->SetBranchAddress("daughters_neutral_pz", &daughters_neutral_pz, &b_daughters_neutral_pz);
  fChain->SetBranchAddress("daughters_neutral_e", &daughters_neutral_e, &b_daughters_neutral_e);
  fChain->SetBranchAddress("daughters_hasTES", &daughters_hasTES, &b_daughters_hasTES);
  /* fChain->SetBranchAddress("daughters_px_TauUp", &daughters_px_TauUp, &b_daughters_px_TauUp); */
  /* fChain->SetBranchAddress("daughters_py_TauUp", &daughters_py_TauUp, &b_daughters_py_TauUp); */
  /* fChain->SetBranchAddress("daughters_pz_TauUp", &daughters_pz_TauUp, &b_daughters_pz_TauUp); */
  /* fChain->SetBranchAddress("daughters_e_TauUp", &daughters_e_TauUp, &b_daughters_e_TauUp); */
  fChain->SetBranchAddress("daughters_hasEES", &daughters_hasEES, &b_daughters_hasEES);
  /* fChain->SetBranchAddress("daughters_px_TauDown", &daughters_px_TauDown, &b_daughters_px_TauDown); */
  /* fChain->SetBranchAddress("daughters_py_TauDown", &daughters_py_TauDown, &b_daughters_py_TauDown); */
  /* fChain->SetBranchAddress("daughters_pz_TauDown", &daughters_pz_TauDown, &b_daughters_pz_TauDown); */
  /* fChain->SetBranchAddress("daughters_e_TauDown", &daughters_e_TauDown, &b_daughters_e_TauDown); */
  /* fChain->SetBranchAddress("daughters_px_EleUp", &daughters_px_EleUp,&b_daughters_px_EleUp); */
  /* fChain->SetBranchAddress("daughters_py_EleUp", &daughters_py_EleUp,&b_daughters_py_EleUp); */
  /* fChain->SetBranchAddress("daughters_pz_EleUp", &daughters_pz_EleUp,&b_daughters_pz_EleUp); */
  /* fChain->SetBranchAddress("daughters_e_EleUp", &daughters_e_EleUp,&b_daughters_e_EleUp); */
  /* fChain->SetBranchAddress("daughters_px_EleDown", &daughters_px_EleDown,&b_daughters_px_EleDown); */
  /* fChain->SetBranchAddress("daughters_py_EleDown", &daughters_py_EleDown,&b_daughters_py_EleDown); */
  /* fChain->SetBranchAddress("daughters_pz_EleDown", &daughters_pz_EleDown,&b_daughters_pz_EleDown); */
  /* fChain->SetBranchAddress("daughters_e_EleDown", &daughters_e_EleDown,&b_daughters_e_EleDown); */
  fChain->SetBranchAddress("daughters_isTauMatched",&daughters_isTauMatched,&b_daughters_isTauMatched);
  fChain->SetBranchAddress("PUNumInteractions", &PUNumInteractions, &b_PUNumInteractions);
  fChain->SetBranchAddress("daughters_genindex", &daughters_genindex, &b_daughters_genindex);
  fChain->SetBranchAddress("MC_weight", &MC_weight, &b_MC_weight);
  fChain->SetBranchAddress("MC_weight_scale_muF0p5", &MC_weight_scale_muF0p5, &b_MC_weight_scale_muF0p5);
  fChain->SetBranchAddress("MC_weight_scale_muF2", &MC_weight_scale_muF2, &b_MC_weight_scale_muF2);
  fChain->SetBranchAddress("MC_weight_scale_muR0p5", &MC_weight_scale_muR0p5, &b_MC_weight_scale_muR0p5);
  fChain->SetBranchAddress("MC_weight_scale_muR2", &MC_weight_scale_muR2, &b_MC_weight_scale_muR2);
  fChain->SetBranchAddress("_nominal_wt", &nominal_wt, &b_nominal_wt);
  fChain->SetBranchAddress("TheoreticalPSUnc", &TheoreticalPSUnc, &b_TheoreticalPSUnc);
  fChain->SetBranchAddress("TheoreticalScaleUnc1005", &TheoreticalScaleUnc1005, &b_TheoreticalScaleUnc1005);
  fChain->SetBranchAddress("TheoreticalScaleUnc1009", &TheoreticalScaleUnc1009, &b_TheoreticalScaleUnc1009);
  fChain->SetBranchAddress("TheoreticalScaleUnc5", &TheoreticalScaleUnc5, &b_TheoreticalScaleUnc5);
  fChain->SetBranchAddress("TheoreticalScaleUnc9", &TheoreticalScaleUnc9, &b_TheoreticalScaleUnc9);
  
  fChain->SetBranchAddress("lheHt", &lheHt, &b_lheHt);
  fChain->SetBranchAddress("lheNOutPartons", &lheNOutPartons, &b_lheNOutPartons);
  fChain->SetBranchAddress("lheNOutB", &lheNOutB, &b_lheNOutB);
  fChain->SetBranchAddress("lheNOutC", &lheNOutC, &b_lheNOutC);
  //  fChain->SetBranchAddress("aMCatNLOweight", &aMCatNLOweight, &b_aMCatNLOweight);
  fChain->SetBranchAddress("genpart_px", &genpart_px, &b_genpart_px);
  fChain->SetBranchAddress("genpart_py", &genpart_py, &b_genpart_py);
  fChain->SetBranchAddress("genpart_pz", &genpart_pz, &b_genpart_pz);
  fChain->SetBranchAddress("genpart_e", &genpart_e, &b_genpart_e);
  fChain->SetBranchAddress("DataMC_Type_idx", &DataMC_Type_idx, &b_DataMC_Type_idx);
  fChain->SetBranchAddress("Event_isRealData", &Event_isRealData, &b_Event_isRealData);
  fChain->SetBranchAddress("genpart_pca_x", &genpart_pca_x, &b_genpart_pca_x);
  fChain->SetBranchAddress("genpart_pca_y", &genpart_pca_y, &b_genpart_pca_y);
  fChain->SetBranchAddress("genpart_pca_z", &genpart_pca_z, &b_genpart_pca_z);
  fChain->SetBranchAddress("genpart_pdg", &genpart_pdg, &b_genpart_pdg);
  fChain->SetBranchAddress("genpart_status", &genpart_status, &b_genpart_status);
  fChain->SetBranchAddress("genpart_HMothInd", &genpart_HMothInd, &b_genpart_HMothInd);
  fChain->SetBranchAddress("genpart_MSSMHMothInd", &genpart_MSSMHMothInd, &b_genpart_MSSMHMothInd);
  fChain->SetBranchAddress("genpart_TopMothInd", &genpart_TopMothInd, &b_genpart_TopMothInd);
  fChain->SetBranchAddress("genpart_TauMothInd", &genpart_TauMothInd, &b_genpart_TauMothInd);
  fChain->SetBranchAddress("genpart_ZMothInd", &genpart_ZMothInd, &b_genpart_ZMothInd);
  fChain->SetBranchAddress("genpart_WMothInd", &genpart_WMothInd, &b_genpart_WMothInd);
  fChain->SetBranchAddress("genpart_bMothInd", &genpart_bMothInd, &b_genpart_bMothInd);
  fChain->SetBranchAddress("genpart_HZDecayMode", &genpart_HZDecayMode, &b_genpart_HZDecayMode);
  fChain->SetBranchAddress("genpart_TopDecayMode", &genpart_TopDecayMode, &b_genpart_TopDecayMode);
  fChain->SetBranchAddress("genpart_WDecayMode", &genpart_WDecayMode, &b_genpart_WDecayMode);
  fChain->SetBranchAddress("genpart_TauGenDecayMode", &genpart_TauGenDecayMode, &b_genpart_TauGenDecayMode);
  fChain->SetBranchAddress("genpart_TauGenDetailedDecayMode", &genpart_TauGenDetailedDecayMode, &b_genpart_TauGenDetailedDecayMode);
  fChain->SetBranchAddress("genpart_flags", &genpart_flags, &b_genpart_flags);
  fChain->SetBranchAddress("genjet_px", &genjet_px, &b_genjet_px);
  fChain->SetBranchAddress("genjet_py", &genjet_py, &b_genjet_py);
  fChain->SetBranchAddress("genjet_pz", &genjet_pz, &b_genjet_pz);
  fChain->SetBranchAddress("genjet_e", &genjet_e, &b_genjet_e);
  fChain->SetBranchAddress("genjet_partonFlavour", &genjet_partonFlavour, &b_genjet_partonFlavour);
  fChain->SetBranchAddress("genjet_hadronFlavour", &genjet_hadronFlavour, &b_genjet_hadronFlavour);
  fChain->SetBranchAddress("NUP", &NUP, &b_NUP);
  /* fChain->SetBranchAddress("SVfit_fitMETPhiTauUp", &SVfit_fitMETPhiTauUp, &b_SVfit_fitMETPhiTauUp); */
  /* fChain->SetBranchAddress("SVfit_fitMETPhiTauDown", &SVfit_fitMETPhiTauDown, &b_SVfit_fitMETPhiTauDown); */
  /* fChain->SetBranchAddress("SVfit_fitMETRhoTauUp", &SVfit_fitMETRhoTauUp, &b_SVfit_fitMETRhoTauUp); */
  /* fChain->SetBranchAddress("SVfit_fitMETRhoTauDown", &SVfit_fitMETRhoTauDown, &b_SVfit_fitMETRhoTauDown); */
  /* fChain->SetBranchAddress("SVfit_phiUncTauUp", &SVfit_phiUncTauUp, &b_SVfit_phiUncTauUp); */
  /* fChain->SetBranchAddress("SVfit_phiUncTauDown", &SVfit_phiUncTauDown, &b_SVfit_phiUncTauDown); */
  /* fChain->SetBranchAddress("SVfit_phiTauUp", &SVfit_phiTauUp, &b_SVfit_phiTauUp); */
  /* fChain->SetBranchAddress("SVfit_phiTauDown", &SVfit_phiTauDown, &b_SVfit_phiTauDown); */
  /* fChain->SetBranchAddress("SVfit_etaUncTauUp", &SVfit_etaUncTauUp, &b_SVfit_etaUncTauUp); */
  /* fChain->SetBranchAddress("SVfit_etaUncTauDown", &SVfit_etaUncTauDown, &b_SVfit_etaUncTauDown); */
  /* fChain->SetBranchAddress("SVfit_etaTauUp", &SVfit_etaTauUp, &b_SVfit_etaTauUp); */
  /* fChain->SetBranchAddress("SVfit_etaTauDown", &SVfit_etaTauDown, &b_SVfit_etaTauDown); */
  /* fChain->SetBranchAddress("SVfit_ptUncTauUp", &SVfit_ptUncTauUp, &b_SVfit_ptUncTauUp); */
  /* fChain->SetBranchAddress("SVfit_ptUncTauDown", &SVfit_ptUncTauDown, &b_SVfit_ptUncTauDown); */
  /* fChain->SetBranchAddress("SVfit_ptTauUp", &SVfit_ptTauUp, &b_SVfit_ptTauUp); */
  /* fChain->SetBranchAddress("SVfit_ptTauDown", &SVfit_ptTauDown, &b_SVfit_ptTauDown); */
  /* fChain->SetBranchAddress("SVfitTransverseMassTauUp", &SVfitTransverseMassTauUp, &b_SVfitTransverseMassTauUp); */
  /* fChain->SetBranchAddress("SVfitTransverseMassTauDown", &SVfitTransverseMassTauDown, &b_SVfitTransverseMassTauDown); */
  /* fChain->SetBranchAddress("SVfitMassTauUp", &SVfitMassTauUp, &b_SVfitMassTauUp); */
  /* fChain->SetBranchAddress("SVfitMassTauDown", &SVfitMassTauDown, &b_SVfitMassTauDown); */
  fChain->SetBranchAddress("MCSignalParticle_p4", &MCSignalParticle_p4, &b_MCSignalParticle_p4);
  fChain->SetBranchAddress("MCSignalParticle_pdgid", &MCSignalParticle_pdgid, &b_MCSignalParticle_pdgid);
  fChain->SetBranchAddress("MCSignalParticle_charge", &MCSignalParticle_charge, &b_MCSignalParticle_charge);
  fChain->SetBranchAddress("MCSignalParticle_Poca", &MCSignalParticle_Poca, &b_MCSignalParticle_Poca);
  fChain->SetBranchAddress("MCSignalParticle_Tauidx", &MCSignalParticle_Tauidx, &b_MCSignalParticle_Tauidx);
  fChain->SetBranchAddress("MCTauandProd_p4", &MCTauandProd_p4, &b_MCTauandProd_p4);
  fChain->SetBranchAddress("MCTauandProd_Vertex", &MCTauandProd_Vertex, &b_MCTauandProd_Vertex);
  fChain->SetBranchAddress("MCTauandProd_pdgid", &MCTauandProd_pdgid, &b_MCTauandProd_pdgid);
  fChain->SetBranchAddress("MCTauandProd_midx", &MCTauandProd_midx, &b_MCTauandProd_midx);
  fChain->SetBranchAddress("MCTauandProd_charge", &MCTauandProd_charge, &b_MCTauandProd_charge);
  fChain->SetBranchAddress("MCTau_JAK", &MCTau_JAK, &b_MCTau_JAK);
  fChain->SetBranchAddress("MCTau_DecayBitMask", &MCTau_DecayBitMask, &b_MCTau_DecayBitMask);
  fChain->SetBranchAddress("MC_p4", &MC_p4, &b_MC_p4);
  fChain->SetBranchAddress("MC_pdgid", &MC_pdgid, &b_MC_pdgid);
  fChain->SetBranchAddress("MC_charge", &MC_charge, &b_MC_charge);
  fChain->SetBranchAddress("MC_midx", &MC_midx, &b_MC_midx);
  fChain->SetBranchAddress("MC_childpdgid", &MC_childpdgid, &b_MC_childpdgid);
  fChain->SetBranchAddress("MC_childidx", &MC_childidx, &b_MC_childidx);
  fChain->SetBranchAddress("MC_status", &MC_status, &b_MC_status);
  /* fChain->SetBranchAddress("SVfitMass", &SVfitMass, &b_SVfitMass); */
  /* fChain->SetBranchAddress("SVfitTransverseMass", &SVfitTransverseMass, &b_SVfitTransverseMass); */
  /* fChain->SetBranchAddress("SVfit_pt", &SVfit_pt, &b_SVfit_pt); */
  /* fChain->SetBranchAddress("SVfit_ptUnc", &SVfit_ptUnc, &b_SVfit_ptUnc); */
  /* fChain->SetBranchAddress("SVfit_eta", &SVfit_eta, &b_SVfit_eta); */
  /* fChain->SetBranchAddress("SVfit_etaUnc", &SVfit_etaUnc, &b_SVfit_etaUnc); */
  /* fChain->SetBranchAddress("SVfit_phi", &SVfit_phi, &b_SVfit_phi); */
  /* fChain->SetBranchAddress("SVfit_phiUnc", &SVfit_phiUnc, &b_SVfit_phiUnc); */
  /* fChain->SetBranchAddress("SVfit_fitMETRho", &SVfit_fitMETRho, &b_SVfit_fitMETRho); */
  /* fChain->SetBranchAddress("SVfit_fitMETPhi", &SVfit_fitMETPhi, &b_SVfit_fitMETPhi); */
  //  fChain->SetBranchAddress("isOSCand", &isOSCand, &b_isOSCand);
  fChain->SetBranchAddress("METx", &METx, &b_METx);
  fChain->SetBranchAddress("METy", &METy, &b_METy);
  //fChain->SetBranchAddress("uncorrMETx", &uncorrMETx, &b_uncorrMETx);
  //fChain->SetBranchAddress("uncorrMETy", &uncorrMETy, &b_uncorrMETy);
  /* fChain->SetBranchAddress("METx_UP_JES",&metx_up_jes,&b_metx_up_jes); */
  /* fChain->SetBranchAddress("METy_UP_JES",&mety_up_jes,&b_mety_up_jes); */
  /* fChain->SetBranchAddress("METx_DOWN_JES",&metx_down_jes,&b_metx_down_jes); */
  /* fChain->SetBranchAddress("METy_DOWN_JES",&mety_down_jes,&b_mety_down_jes); */
  /* fChain->SetBranchAddress("METx_UP_TES",&metx_up_tes,&b_metx_up_tes); */
  /* fChain->SetBranchAddress("METy_UP_TES",&mety_up_tes,&b_mety_up_tes); */
  /* fChain->SetBranchAddress("METx_DOWN_TES",&metx_down_tes,&b_metx_down_tes); */
  /* fChain->SetBranchAddress("METy_DOWN_TES",&mety_down_tes,&b_mety_down_tes); */
  /* fChain->SetBranchAddress("METx_UP_EES",&metx_up_ees,&b_metx_up_ees); */
  /* fChain->SetBranchAddress("METy_UP_EES",&mety_up_ees,&b_mety_up_ees); */
  /* fChain->SetBranchAddress("METx_DOWN_EES",&metx_down_ees,&b_metx_down_ees); */
  /* fChain->SetBranchAddress("METy_DOWN_EES",&mety_down_ees,&b_mety_down_ees); */
  /* fChain->SetBranchAddress("MET_cov00", &MET_cov00, &b_MET_cov00); */
  /* fChain->SetBranchAddress("MET_cov01", &MET_cov01, &b_MET_cov01); */
  /* fChain->SetBranchAddress("MET_cov10", &MET_cov10, &b_MET_cov10); */
  /* fChain->SetBranchAddress("MET_cov11", &MET_cov11, &b_MET_cov11); */
  /* fChain->SetBranchAddress("MET_significance", &MET_significance, &b_MET_significance); */
  /* fChain->SetBranchAddress("mT_Dau1", &mT_Dau1, &b_mT_Dau1); */
  /* fChain->SetBranchAddress("mT_Dau2", &mT_Dau2, &b_mT_Dau2); */
  fChain->SetBranchAddress("PDGIdDaughters", &PDGIdDaughters, &b_PDGIdDaughters);
  fChain->SetBranchAddress("indexDau1", &indexDau1, &b_indexDau1);
  fChain->SetBranchAddress("indexDau2", &indexDau2, &b_indexDau2);
  fChain->SetBranchAddress("particleType", &particleType, &b_particleType);
  //fChain->SetBranchAddress("discriminator", &discriminator, &b_discriminator);
  fChain->SetBranchAddress("daughters_muonID", &daughters_muonID, &b_daughters_muonID);
  fChain->SetBranchAddress("daughters_typeOfMuon", &daughters_typeOfMuon, &b_daughters_typeOfMuon);
  fChain->SetBranchAddress("dxy", &dxy, &b_dxy);
  fChain->SetBranchAddress("dz", &dz, &b_dz);
  fChain->SetBranchAddress("dxy_innerTrack", &dxy_innerTrack, &b_dxy_innerTrack);
  //  fChain->SetBranchAddress("dz_innerTrack", &dz_innerTrack, &b_dz_innerTrack);
  fChain->SetBranchAddress("daughters_rel_error_trackpt", &daughters_rel_error_trackpt, &b_daughters_rel_error_trackpt);
  fChain->SetBranchAddress("SIP", &SIP, &b_SIP);
  //fChain->SetBranchAddress("daughters_iseleBDT", &daughters_iseleBDT, &b_daughters_iseleBDT);
  fChain->SetBranchAddress("daughters_iseleWPLoose",&daughters_iseleWPLoose, &b_daughters_iseleWPLoose);
  fChain->SetBranchAddress("daughters_iseleWP80", &daughters_iseleWP80, &b_daughters_iseleWP80);
  fChain->SetBranchAddress("daughters_iseleWP90", &daughters_iseleWP90, &b_daughters_iseleWP90);
  fChain->SetBranchAddress("daughters_iseleNoIsoWPLoose",&daughters_iseleNoIsoWPLoose, &b_daughters_iseleNoIsoWPLoose);
  fChain->SetBranchAddress("daughters_iseleNoIsoWP80",&daughters_iseleNoIsoWP80, &b_daughters_iseleNoIsoWP80);
  fChain->SetBranchAddress("daughters_iseleNoIsoWP90",&daughters_iseleNoIsoWP90, &b_daughters_iseleNoIsoWP90);
  fChain->SetBranchAddress("daughters_eleMVAnt", &daughters_eleMVAnt, &b_daughters_eleMVAnt);
  fChain->SetBranchAddress("daughters_eleMVA_HZZ", &daughters_eleMVA_HZZ, &b_daughters_eleMVA_HZZ);
  fChain->SetBranchAddress("daughters_passConversionVeto", &daughters_passConversionVeto, &b_daughters_passConversionVeto);
  fChain->SetBranchAddress("daughters_eleMissingHits", &daughters_eleMissingHits, &b_daughters_eleMissingHits);
  fChain->SetBranchAddress("daughters_iseleChargeConsistent", &daughters_iseleChargeConsistent, &b_daughters_iseleChargeConsistent);
  //fChain->SetBranchAddress("daughters_eleCUTID", &daughters_eleCUTID, &b_daughters_eleCUTID);
  fChain->SetBranchAddress("decayMode", &decayMode, &b_decayMode);
  fChain->SetBranchAddress("genmatch",&genmatch,&b_genmatch);
  fChain->SetBranchAddress("tauID", &tauID, &b_tauID);
  //  fChain->SetBranchAddress("MVADM2016v1", &MVADM2016, &b_MVADM2016);
  fChain->SetBranchAddress("MVADM2017v1", &MVADM2017, &b_MVADM2017);
  fChain->SetBranchAddress("combreliso", &combreliso, &b_combreliso);
  fChain->SetBranchAddress("combreliso03", &combreliso03, &b_combreliso03);
  fChain->SetBranchAddress("daughters_depositR03_tracker", &daughters_depositR03_tracker, &b_daughters_depositR03_tracker);
  fChain->SetBranchAddress("daughters_depositR03_ecal", &daughters_depositR03_ecal, &b_daughters_depositR03_ecal);
  fChain->SetBranchAddress("daughters_depositR03_hcal", &daughters_depositR03_hcal, &b_daughters_depositR03_hcal);
  fChain->SetBranchAddress("daughters_decayModeFindingOldDMs", &daughters_decayModeFindingOldDMs, &b_daughters_decayModeFindingOldDMs);
  //fChain->SetBranchAddress("daughters_SCeta", &daughters_SCeta, &b_daughters_SCeta);
  //fChain->SetBranchAddress("againstElectronMVA5category", &againstElectronMVA5category, &b_againstElectronMVA5category);
  //fChain->SetBranchAddress("againstElectronMVA5raw", &againstElectronMVA5raw, &b_againstElectronMVA5raw);
  //fChain->SetBranchAddress("byPileupWeightedIsolationRaw3Hits", &byPileupWeightedIsolationRaw3Hits, &b_byPileupWeightedIsolationRaw3Hits);
  fChain->SetBranchAddress("footprintCorrection", &footprintCorrection, &b_footprintCorrection);
  fChain->SetBranchAddress("neutralIsoPtSumWeight", &neutralIsoPtSumWeight, &b_neutralIsoPtSumWeight);
  fChain->SetBranchAddress("photonPtSumOutsideSignalCone", &photonPtSumOutsideSignalCone, &b_photonPtSumOutsideSignalCone);
  fChain->SetBranchAddress("daughters_decayModeFindingNewDMs", &daughters_decayModeFindingNewDMs, &b_daughters_decayModeFindingNewDMs);
  //fChain->SetBranchAddress("daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits", &daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits, &b_daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits);
  //fChain->SetBranchAddress("daughters_byIsolationMVA3oldDMwoLTraw", &daughters_byIsolationMVA3oldDMwoLTraw, &b_daughters_byIsolationMVA3oldDMwoLTraw);
  //fChain->SetBranchAddress("daughters_byIsolationMVArun2017v2DBoldDMwLTraw2017",&daughters_byIsolationMVArun2017v2DBoldDMwLTraw2017, &b_daughters_byIsolationMVArun2017v2DBoldDMwLTraw2017);
  //  fChain->SetBranchAddress("daughters_byIsolationMVArun2017v1DBoldDMwLTraw2017",&daughters_byIsolationMVArun2017v1DBoldDMwLTraw2017, &b_daughters_byIsolationMVArun2017v1DBoldDMwLTraw2017);
  //fChain->SetBranchAddress("daughters_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017",&daughters_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017, &b_daughters_byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017);
  //fChain->SetBranchAddress("daughters_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017", &daughters_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017, &b_daughters_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017);
  //fChain->SetBranchAddress("daughters_byIsolationMVA3oldDMwLTraw", &daughters_byIsolationMVA3oldDMwLTraw, &b_daughters_byIsolationMVA3oldDMwLTraw);
  //fChain->SetBranchAddress("daughters_byIsolationMVA3newDMwoLTraw", &daughters_byIsolationMVA3newDMwoLTraw, &b_daughters_byIsolationMVA3newDMwoLTraw);
  //fChain->SetBranchAddress("daughters_byIsolationMVA3newDMwLTraw", &daughters_byIsolationMVA3newDMwLTraw, &b_daughters_byIsolationMVA3newDMwLTraw);
  //fChain->SetBranchAddress("daughters_byIsolationMVArun2v1DBoldDMwLTraw", &daughters_byIsolationMVArun2v1DBoldDMwLTraw, &b_daughters_byIsolationMVArun2v1DBoldDMwLTraw);
  fChain->SetBranchAddress("daughters_byDeepTau2017v2p1VSjetraw", &daughters_byDeepTau2017v2p1VSjetraw, &b_daughters_byDeepTau2017v2p1VSjetraw);
  fChain->SetBranchAddress("daughters_byDeepTau2017v2p1VSeraw", &daughters_byDeepTau2017v2p1VSeraw, &b_daughters_byDeepTau2017v2p1VSeraw);
  fChain->SetBranchAddress("daughters_byDeepTau2017v2p1VSmuraw", &daughters_byDeepTau2017v2p1VSmuraw, &b_daughters_byDeepTau2017v2p1VSmuraw);
  fChain->SetBranchAddress("daughters_chargedIsoPtSum", &daughters_chargedIsoPtSum, &b_daughters_chargedIsoPtSum);
  fChain->SetBranchAddress("daughters_neutralIsoPtSum", &daughters_neutralIsoPtSum, &b_daughters_neutralIsoPtSum);
  fChain->SetBranchAddress("daughters_puCorrPtSum", &daughters_puCorrPtSum, &b_daughters_puCorrPtSum);
  fChain->SetBranchAddress("daughters_numChargedParticlesSignalCone", &daughters_numChargedParticlesSignalCone, &b_daughters_numChargedParticlesSignalCone);
  fChain->SetBranchAddress("daughters_numNeutralHadronsSignalCone", &daughters_numNeutralHadronsSignalCone, &b_daughters_numNeutralHadronsSignalCone);
  fChain->SetBranchAddress("daughters_numPhotonsSignalCone", &daughters_numPhotonsSignalCone, &b_daughters_numPhotonsSignalCone);
  fChain->SetBranchAddress("daughters_daughters_numParticlesSignalCone", &daughters_daughters_numParticlesSignalCone, &b_daughters_daughters_numParticlesSignalCone);
  fChain->SetBranchAddress("daughters_numChargedParticlesIsoCone", &daughters_numChargedParticlesIsoCone, &b_daughters_numChargedParticlesIsoCone);
  fChain->SetBranchAddress("daughters_numNeutralHadronsIsoCone", &daughters_numNeutralHadronsIsoCone, &b_daughters_numNeutralHadronsIsoCone);
  fChain->SetBranchAddress("daughters_numPhotonsIsoCone", &daughters_numPhotonsIsoCone, &b_daughters_numPhotonsIsoCone);
  fChain->SetBranchAddress("daughters_numParticlesIsoCone", &daughters_numParticlesIsoCone, &b_daughters_numParticlesIsoCone);
  fChain->SetBranchAddress("daughters_leadChargedParticlePt", &daughters_leadChargedParticlePt, &b_daughters_leadChargedParticlePt);
  fChain->SetBranchAddress("daughters_trackRefPt", &daughters_trackRefPt, &b_daughters_trackRefPt);
  //fChain->SetBranchAddress("daughters_isLastTriggerObjectforPath", &daughters_isLastTriggerObjectforPath, &b_daughters_isLastTriggerObjectforPath);
  fChain->SetBranchAddress("daughters_trgMatched", &daughters_trgMatched, &b_daughters_trgMatched);
  //fChain->SetBranchAddress("daughters_isTriggerObjectforPath", &daughters_isTriggerObjectforPath, &b_daughters_isTriggerObjectforPath);
  fChain->SetBranchAddress("daughters_FilterFired", &daughters_FilterFired, &b_daughters_FilterFired);
  fChain->SetBranchAddress("daughters_isGoodTriggerType", &daughters_isGoodTriggerType, &b_daughters_isGoodTriggerType);
  fChain->SetBranchAddress("daughters_L3FilterFired", &daughters_L3FilterFired, &b_daughters_L3FilterFired);
  fChain->SetBranchAddress("daughters_L3FilterFiredLast", &daughters_L3FilterFiredLast, &b_daughters_L3FilterFiredLast);
  fChain->SetBranchAddress("daughters_HLTpt", &daughters_HLTpt, &b_daughters_HLTpt);
  //fChain->SetBranchAddress("daughters_isL1IsoTau28Matched", &daughters_isL1IsoTau28Matched, &b_daughters_isL1IsoTau28Matched);
  fChain->SetBranchAddress("Muon_trackCharge", &Muon_trackCharge, &b_Muon_trackCharge);
  fChain->SetBranchAddress("Muon_pdgid", &Muon_pdgid, &b_Muon_pdgid);
  fChain->SetBranchAddress("Muon_B", &Muon_B, &b_Muon_B);
  fChain->SetBranchAddress("Muon_M", &Muon_M, &b_Muon_M);
  fChain->SetBranchAddress("Muon_par", &Muon_par, &b_Muon_par);
  //  fChain->SetBranchAddress("Muon_cov", &Muon_cov, &b_Muon_cov);
  /* fChain->SetBranchAddress("PFTau_Track_par", &PFTau_Track_par, &b_PFTau_Track_par); */
  /* fChain->SetBranchAddress("PFTau_Track_cov", &PFTau_Track_cov, &b_PFTau_Track_cov); */
  /* fChain->SetBranchAddress("PFTau_Track_charge", &PFTau_Track_charge, &b_PFTau_Track_charge); */
  /* fChain->SetBranchAddress("PFTau_Track_pdgid", &PFTau_Track_pdgid, &b_PFTau_Track_pdgid); */
  /* fChain->SetBranchAddress("PFTau_Track_B", &PFTau_Track_B, &b_PFTau_Track_B); */
  /* fChain->SetBranchAddress("PFTau_Track_M", &PFTau_Track_M, &b_PFTau_Track_M); */
  fChain->SetBranchAddress("PFTauSVPos", &PFTauSVPos, &b_PFTauSVPos);
  fChain->SetBranchAddress("PFTauSVCov", &PFTauSVCov, &b_PFTauSVCov);
  fChain->SetBranchAddress("PFTauPionsP4", &PFTauPionsP4, &b_PFTauPionsP4);
  fChain->SetBranchAddress("PFTauRefitPionsP4", &PFTauRefitPionsP4, &b_PFTauRefitPionsP4);
  fChain->SetBranchAddress("PFTauPionsCharge", &PFTauPionsCharge, &b_PFTauPionsCharge);
  fChain->SetBranchAddress("PFTauRefitPionsCharge", &PFTauRefitPionsCharge, &b_PFTauRefitPionsCharge);
  fChain->SetBranchAddress("PFTauSVChi2NDofMatchingQuality", &PFTauSVChi2NDofMatchingQuality, &b_PFTauSVChi2NDofMatchingQuality);
  fChain->SetBranchAddress("TauFLSignificance", &TauFLSignificance, &b_TauFLSignificance);
  fChain->SetBranchAddress("PFTauGEOMFlightLenght", &PFTauGEOMFlightLenght, &b_PFTauGEOMFlightLenght);
  fChain->SetBranchAddress("PFTauGEOMFlightLenghtSignificance", &PFTauGEOMFlightLenghtSignificance, &b_PFTauGEOMFlightLenghtSignificance);
  fChain->SetBranchAddress("daughters_jetNDauChargedMVASel", &daughters_jetNDauChargedMVASel, &b_daughters_jetNDauChargedMVASel);
  fChain->SetBranchAddress("daughters_miniRelIsoCharged", &daughters_miniRelIsoCharged, &b_daughters_miniRelIsoCharged);
  fChain->SetBranchAddress("daughters_miniRelIsoNeutral", &daughters_miniRelIsoNeutral, &b_daughters_miniRelIsoNeutral);
  fChain->SetBranchAddress("daughters_jetPtRel", &daughters_jetPtRel, &b_daughters_jetPtRel);
  fChain->SetBranchAddress("daughters_jetPtRatio", &daughters_jetPtRatio, &b_daughters_jetPtRatio);
  fChain->SetBranchAddress("daughters_jetBTagCSV", &daughters_jetBTagCSV, &b_daughters_jetBTagCSV);
  fChain->SetBranchAddress("daughters_jetBTagDeepCSV",&daughters_jetBTagDeepCSV, &b_daughters_jetBTagDeepCSV);
  fChain->SetBranchAddress("daughters_jetBTagDeepFlavor",&daughters_jetBTagDeepFlavor, &b_daughters_jetBTagDeepFlavor);
  //fChain->SetBranchAddress("daughters_lepMVA_mvaId", &daughters_lepMVA_mvaId, &b_daughters_lepMVA_mvaId);
  fChain->SetBranchAddress("daughters_pca_x", &daughters_pca_x, &b_daughters_pca_x);
  fChain->SetBranchAddress("daughters_pca_y", &daughters_pca_y, &b_daughters_pca_y);
  fChain->SetBranchAddress("daughters_pca_z", &daughters_pca_z, &b_daughters_pca_z);
  fChain->SetBranchAddress("daughters_pcaRefitPV_x", &daughters_pcaRefitPV_x, &b_daughters_pcaRefitPV_x);
  fChain->SetBranchAddress("daughters_pcaRefitPV_y", &daughters_pcaRefitPV_y, &b_daughters_pcaRefitPV_y);
  fChain->SetBranchAddress("daughters_pcaRefitPV_z", &daughters_pcaRefitPV_z, &b_daughters_pcaRefitPV_z);
  fChain->SetBranchAddress("daughters_pcaGenPV_x", &daughters_pcaGenPV_x, &b_daughters_pcaGenPV_x);
  fChain->SetBranchAddress("daughters_pcaGenPV_y", &daughters_pcaGenPV_y, &b_daughters_pcaGenPV_y);
  fChain->SetBranchAddress("daughters_pcaGenPV_z", &daughters_pcaGenPV_z, &b_daughters_pcaGenPV_z);
  //fChain->SetBranchAddress("JetsNumber", &JetsNumber, &b_JetsNumber);
  fChain->SetBranchAddress("PFTau_a1_lvp", &PFTau_a1_lvp, &b_PFTau_a1_lvp);
  //  fChain->SetBranchAddress("PFTau_a1_cov", &PFTau_a1_cov, &b_PFTau_a1_cov);
  fChain->SetBranchAddress("PFTau_a1_charge", &PFTau_a1_charge, &b_PFTau_a1_charge);
  fChain->SetBranchAddress("PFTau_a1_pdgid", &PFTau_a1_pdgid, &b_PFTau_a1_pdgid);
  fChain->SetBranchAddress("PFTau_a1_B", &PFTau_a1_B, &b_PFTau_a1_B);
  fChain->SetBranchAddress("PFTau_a1_M", &PFTau_a1_M, &b_PFTau_a1_M);
  fChain->SetBranchAddress("PFTauTrack_deltaR", &PFTauTrack_deltaR, &b_PFTauTrack_deltaR);
  fChain->SetBranchAddress("PFTauLeadTrackLV", &PFTauLeadTrackLV, &b_PFTauLeadTrackLV);
  fChain->SetBranchAddress("jets_VBFleadFilterMatch", &jets_VBFleadFilterMatch, &b_jets_VBFleadFilterMatch);
  fChain->SetBranchAddress("jets_VBFsubleadFilterMatch", &jets_VBFsubleadFilterMatch, &b_jets_VBFsubleadFilterMatch);
  fChain->SetBranchAddress("jets_px", &jets_px, &b_jets_px);
  fChain->SetBranchAddress("jets_py", &jets_py, &b_jets_py);
  fChain->SetBranchAddress("jets_pz", &jets_pz, &b_jets_pz);
  fChain->SetBranchAddress("jets_e", &jets_e, &b_jets_e);
  fChain->SetBranchAddress("jetsDown_px", &jetsDown_px, &b_jetsDown_px);
  fChain->SetBranchAddress("jetsDown_py", &jetsDown_py, &b_jetsDown_py);
  fChain->SetBranchAddress("jetsDown_pz", &jetsDown_pz, &b_jetsDown_pz);
  fChain->SetBranchAddress("jetsDown_e", &jetsDown_e, &b_jetsDown_e);
  fChain->SetBranchAddress("jetsUp_px", &jetsUp_px, &b_jetsUp_px);
  fChain->SetBranchAddress("jetsUp_py", &jetsUp_py, &b_jetsUp_py);
  fChain->SetBranchAddress("jetsUp_pz", &jetsUp_pz, &b_jetsUp_pz);
  fChain->SetBranchAddress("jetsUp_e", &jetsUp_e, &b_jetsUp_e);
  fChain->SetBranchAddress("jetRawf", &jetRawf, &b_jetRawf);
  fChain->SetBranchAddress("jets_area", &jets_area, &b_jets_area);
  fChain->SetBranchAddress("jets_JER",&jets_JER, &b_jets_JER);
  //  fChain->SetBranchAddress("jets_QGdiscr" , &jets_QGdiscr);
  fChain->SetBranchAddress("jets_mT", &jets_mT, &b_jets_mT);
  fChain->SetBranchAddress("jets_Flavour", &jets_Flavour, &b_jets_Flavour);
  fChain->SetBranchAddress("jets_HadronFlavour", &jets_HadronFlavour, &b_jets_HadronFlavour);
  fChain->SetBranchAddress("jets_genjetIndex", &jets_genjetIndex, &b_jets_genjetIndex);
  fChain->SetBranchAddress("jets_PUJetID", &jets_PUJetID, &b_jets_PUJetID);
  fChain->SetBranchAddress("jets_PUJetIDupdated", &jets_PUJetIDupdated, &b_jets_PUJetIDupdated);
  fChain->SetBranchAddress("jets_PUJetID_WP",&jets_PUJetID_WP, &b_jets_PUJetID_WP);
  fChain->SetBranchAddress("jets_PUJetIDupdated_WP",&jets_PUJetIDupdated_WP, &b_jets_PUJetIDupdated_WP);
  fChain->SetBranchAddress("jets_vtxPt", &jets_vtxPt, &b_jets_vtxPt);
  fChain->SetBranchAddress("jets_vtxMass", &jets_vtxMass, &b_jets_vtxMass);
  fChain->SetBranchAddress("jets_vtx3dL", &jets_vtx3dL, &b_jets_vtx3dL);
  fChain->SetBranchAddress("jets_vtxNtrk", &jets_vtxNtrk, &b_jets_vtxNtrk);
  fChain->SetBranchAddress("jets_vtx3deL", &jets_vtx3deL, &b_jets_vtx3deL);
  fChain->SetBranchAddress("jets_leadTrackPt", &jets_leadTrackPt, &b_jets_leadTrackPt);
  fChain->SetBranchAddress("jets_leptonPtRel", &jets_leptonPtRel, &b_jets_leptonPtRel);
  fChain->SetBranchAddress("jets_leptonPt", &jets_leptonPt, &b_jets_leptonPt);
  fChain->SetBranchAddress("jets_leptonDeltaR", &jets_leptonDeltaR, &b_jets_leptonDeltaR);
  fChain->SetBranchAddress("jets_chEmEF", &_jets_chEmEF, &b_jets_chEmEF);
  fChain->SetBranchAddress("jets_chHEF", &_jets_chHEF, &b_jets_chHEF);
  fChain->SetBranchAddress("jets_nEmEF", &_jets_nEmEF, &b_jets_nEmEF);
  fChain->SetBranchAddress("jets_nHEF", &_jets_nHEF, &b_jets_nHEF);
  fChain->SetBranchAddress("jets_MUF", &jets_MUF, &b_jets_MUF);
  fChain->SetBranchAddress("jets_neMult", &_jets_neMult, &b_jets_neMult);
  fChain->SetBranchAddress("jets_chMult", &_jets_chMult, &b_jets_chMult);
  fChain->SetBranchAddress("jets_jecUnc", &jets_jecUnc, &b_jets_jecUnc);
  // JEC Uncertainty sources
  fChain->SetBranchAddress("jets_jetUncRegrouped_Absolute_up"  , &jets_jetUnc_Absolute_up, &b_jets_jetUnc_Absolute_up); // up variations
  fChain->SetBranchAddress("jets_jetUncRegrouped_FlavorQCD_up"        , &jets_jetUnc_FlavorQCD_up,&b_jets_jetUnc_FlavorQCD_up);
  fChain->SetBranchAddress("jets_jetUncRegrouped_RelativeBal_up"      , &jets_jetUnc_RelativeBal_up, &b_jets_jetUnc_RelativeBal_up);
  fChain->SetBranchAddress("jets_jetUncRegrouped_HF_up"      , &jets_jetUnc_HF_up, &b_jets_jetUnc_HF_up);
  fChain->SetBranchAddress("jets_jetUncRegrouped_BBEC1_up"      , &jets_jetUnc_BBEC1_up, &b_jets_jetUnc_BBEC1_up);
  fChain->SetBranchAddress("jets_jetUncRegrouped_EC2_up"      , &jets_jetUnc_EC2_up, &b_jets_jetUnc_EC2_up);
  fChain->SetBranchAddress("jets_jetUncRegrouped_BBEC1_YEAR_up"      , &jets_jetUnc_BBEC1_YEAR_up, &b_jets_jetUnc_BBEC1_YEAR_up);
  fChain->SetBranchAddress("jets_jetUncRegrouped_EC2_YEAR_up"      , &jets_jetUnc_EC2_YEAR_up, &b_jets_jetUnc_EC2_YEAR_up);
  fChain->SetBranchAddress("jets_jetUncRegrouped_Absolute_YEAR_up"      , &jets_jetUnc_Absolute_YEAR_up, &b_jets_jetUnc_Absolute_YEAR_up);
  fChain->SetBranchAddress("jets_jetUncRegrouped_HF_YEAR_up"      , &jets_jetUnc_HF_YEAR_up, &b_jets_jetUnc_HF_YEAR_up);
  fChain->SetBranchAddress("jets_jetUncRegrouped_RelativeSample_YEAR_up"      , &jets_jetUnc_RelativeSample_YEAR_up, &b_jets_jetUnc_RelativeSample_YEAR_up);

  fChain->SetBranchAddress("jets_jetUncRegrouped_Absolute_dw"  , &jets_jetUnc_Absolute_dw, &b_jets_jetUnc_Absolute_dw); // up variations
  fChain->SetBranchAddress("jets_jetUncRegrouped_FlavorQCD_dw"        , &jets_jetUnc_FlavorQCD_dw,&b_jets_jetUnc_FlavorQCD_dw);
  fChain->SetBranchAddress("jets_jetUncRegrouped_RelativeBal_dw"      , &jets_jetUnc_RelativeBal_dw, &b_jets_jetUnc_RelativeBal_dw);
  fChain->SetBranchAddress("jets_jetUncRegrouped_HF_dw"      , &jets_jetUnc_HF_dw, &b_jets_jetUnc_HF_dw);
  fChain->SetBranchAddress("jets_jetUncRegrouped_BBEC1_dw"      , &jets_jetUnc_BBEC1_dw, &b_jets_jetUnc_BBEC1_dw);
  fChain->SetBranchAddress("jets_jetUncRegrouped_EC2_dw"      , &jets_jetUnc_EC2_dw, &b_jets_jetUnc_EC2_dw);
  fChain->SetBranchAddress("jets_jetUncRegrouped_BBEC1_YEAR_dw"      , &jets_jetUnc_BBEC1_YEAR_dw, &b_jets_jetUnc_BBEC1_YEAR_dw);
  fChain->SetBranchAddress("jets_jetUncRegrouped_EC2_YEAR_dw"      , &jets_jetUnc_EC2_YEAR_dw, &b_jets_jetUnc_EC2_YEAR_dw);
  fChain->SetBranchAddress("jets_jetUncRegrouped_Absolute_YEAR_dw"      , &jets_jetUnc_Absolute_YEAR_dw, &b_jets_jetUnc_Absolute_YEAR_dw);
  fChain->SetBranchAddress("jets_jetUncRegrouped_HF_YEAR_dw"      , &jets_jetUnc_HF_YEAR_dw, &b_jets_jetUnc_HF_YEAR_dw);
  fChain->SetBranchAddress("jets_jetUncRegrouped_RelativeSample_YEAR_dw"      , &jets_jetUnc_RelativeSample_YEAR_dw, &b_jets_jetUnc_RelativeSample_YEAR_dw);

  /* fChain->SetBranchAddress("jets_jetUnc_Fragmentation_up"    , &_SourceUncVal_up["Fragmentation"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_PileUpDataMC_up"     , &_SourceUncVal_up["PileUpDataMC"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_PileUpPtBB_up"       , &_SourceUncVal_up["PileUpPtBB"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_PileUpPtEC1_up"      , &_SourceUncVal_up["PileUpPtEC1"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_PileUpPtEC2_up"      , &_SourceUncVal_up["PileUpPtEC2"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_PileUpPtHF_up"       , &_SourceUncVal_up["PileUpPtHF"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_PileUpPtRef_up"      , &_SourceUncVal_up["PileUpPtRef"]); */

  /* fChain->SetBranchAddress("jets_jetUnc_RelativeFSR_up"      , &_SourceUncVal_up["RelativeFSR"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_RelativeJEREC1_up"   , &_SourceUncVal_up["RelativeJEREC1"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_RelativeJEREC2_up"   , &_SourceUncVal_up["RelativeJEREC2"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_RelativeJERHF_up"    , &_SourceUncVal_up["RelativeJERHF"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_RelativePtBB_up"     , &_SourceUncVal_up["RelativePtBB"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_RelativePtEC1_up"    , &_SourceUncVal_up["RelativePtEC1"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_RelativePtEC2_up"    , &_SourceUncVal_up["RelativePtEC2"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_RelativePtHF_up"     , &_SourceUncVal_up["RelativePtHF"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_RelativeSample_up"   , &_SourceUncVal_up["RelativeSample"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_RelativeStatEC_up"   , &_SourceUncVal_up["RelativeStatEC"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_RelativeStatFSR_up"  , &_SourceUncVal_up["RelativeStatFSR"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_RelativeStatHF_up"   , &_SourceUncVal_up["RelativeStatHF"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_SinglePionECAL_up"   , &_SourceUncVal_up["SinglePionECAL"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_SinglePionHCAL_up"   , &_SourceUncVal_up["SinglePionHCAL"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_TimePtEta_up"        , &_SourceUncVal_up["TimePtEta"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_AbsoluteFlavMap_dw"  , &_SourceUncVal_dw["AbsoluteFlavMap"]); // down variations */
  /* fChain->SetBranchAddress("jets_jetUnc_AbsoluteMPFBias_dw"  , &_SourceUncVal_dw["AbsoluteMPFBias"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_AbsoluteSample_dw"   , &_SourceUncVal_dw["AbsoluteSample"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_AbsoluteScale_dw"    , &_SourceUncVal_dw["AbsoluteScale"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_AbsoluteStat_dw"     , &_SourceUncVal_dw["AbsoluteStat"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_FlavorQCD_dw"        , &_SourceUncVal_dw["FlavorQCD"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_Fragmentation_dw"    , &_SourceUncVal_dw["Fragmentation"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_PileUpDataMC_dw"     , &_SourceUncVal_dw["PileUpDataMC"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_PileUpPtBB_dw"       , &_SourceUncVal_dw["PileUpPtBB"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_PileUpPtEC1_dw"      , &_SourceUncVal_dw["PileUpPtEC1"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_PileUpPtEC2_dw"      , &_SourceUncVal_dw["PileUpPtEC2"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_PileUpPtHF_dw"       , &_SourceUncVal_dw["PileUpPtHF"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_PileUpPtRef_dw"      , &_SourceUncVal_dw["PileUpPtRef"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_RelativeBal_dw"      , &_SourceUncVal_dw["RelativeBal"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_RelativeFSR_dw"      , &_SourceUncVal_dw["RelativeFSR"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_RelativeJEREC1_dw"   , &_SourceUncVal_dw["RelativeJEREC1"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_RelativeJEREC2_dw"   , &_SourceUncVal_dw["RelativeJEREC2"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_RelativeJERHF_dw"    , &_SourceUncVal_dw["RelativeJERHF"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_RelativePtBB_dw"     , &_SourceUncVal_dw["RelativePtBB"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_RelativePtEC1_dw"    , &_SourceUncVal_dw["RelativePtEC1"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_RelativePtEC2_dw"    , &_SourceUncVal_dw["RelativePtEC2"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_RelativePtHF_dw"     , &_SourceUncVal_dw["RelativePtHF"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_RelativeSample_dw"   , &_SourceUncVal_dw["RelativeSample"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_RelativeStatEC_dw"   , &_SourceUncVal_dw["RelativeStatEC"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_RelativeStatFSR_dw"  , &_SourceUncVal_dw["RelativeStatFSR"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_RelativeStatHF_dw"   , &_SourceUncVal_dw["RelativeStatHF"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_SinglePionECAL_dw"   , &_SourceUncVal_dw["SinglePionECAL"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_SinglePionHCAL_dw"   , &_SourceUncVal_dw["SinglePionHCAL"]); */
  /* fChain->SetBranchAddress("jets_jetUnc_TimePtEta_dw"        , &_SourceUncVal_dw["TimePtEta"]); */
   
  //fChain->SetBranchAddress("PFJet_chargedMultiplicity", &CHM, &b_PFJet_chargedMultiplicity);
  //fChain->SetBranchAddress("PFJet_neutralMultiplicity", &NumNeutralParticles, &b_PFJet_neutralMultiplicity);
  //fChain->SetBranchAddress("PFJet_chargedEmEnergyFraction", &CEMF, &b_PFJet_chargedEmEnergyFraction);
  //fChain->SetBranchAddress("PFJet_chargedHadronEnergyFraction", &CHF, &b_PFJet_chargedHadronEnergyFraction);
  //fChain->SetBranchAddress("PFJet_neutralHadronEnergyFraction", &NHF, &b_PFJet_neutralHadronEnergyFraction);
  //fChain->SetBranchAddress("PFJet_neutralEmEnergyFraction", &NEMF, &b_PFJet_neutralEmEnergyFraction);
  fChain->SetBranchAddress("bDiscriminator", &bDiscriminator, &b_bDiscriminator);
  fChain->SetBranchAddress("bCSVscore", &bCSVscore, &b_bCSVscore);
  fChain->SetBranchAddress("pfCombinedMVAV2BJetTags", &pfCombinedMVAV2BJetTags, &b_pfCombinedMVAV2BJetTags);
  fChain->SetBranchAddress("bDeepCSV_probb",&bDeepCSV_probb, &b_bDeepCSV_probb);
  fChain->SetBranchAddress("bDeepCSV_probbb",&bDeepCSV_probbb, &b_bDeepCSV_probbb);
  fChain->SetBranchAddress("bDeepCSV_probudsg",&bDeepCSV_probudsg, &b_bDeepCSV_probudsg);
  fChain->SetBranchAddress("bDeepCSV_probc",&bDeepCSV_probc, &b_bDeepCSV_probc);
  fChain->SetBranchAddress("bDeepCSV_probcc",&bDeepCSV_probcc, &b_bDeepCSV_probcc);
  fChain->SetBranchAddress("bDeepFlavor_probb",&bDeepFlavor_probb, &b_bDeepFlavor_probb);
  fChain->SetBranchAddress("bDeepFlavor_probbb",&bDeepFlavor_probbb, &b_bDeepFlavor_probbb);
  fChain->SetBranchAddress("bDeepFlavor_problepb",&bDeepFlavor_problepb, &b_bDeepFlavor_problepb);
  fChain->SetBranchAddress("bDeepFlavor_probc",&bDeepFlavor_probc, &b_bDeepFlavor_probc);
  fChain->SetBranchAddress("bDeepFlavor_probuds",&bDeepFlavor_probuds, &b_bDeepFlavor_probuds);
  fChain->SetBranchAddress("bDeepFlavor_probg",&bDeepFlavor_probg, &b_bDeepFlavor_probg);
  fChain->SetBranchAddress("looseJetID",&looseJetID, &b_looseJetID);
  fChain->SetBranchAddress("tightJetID",&tightJetID, &b_tightJetID);
  fChain->SetBranchAddress("tightLepVetoJetID",&tightLepVetoJetID, &b_tightLepVetoJetID);
  fChain->SetBranchAddress("ak8jets_px", &ak8jets_px, &b_ak8jets_px);
  fChain->SetBranchAddress("ak8jets_py", &ak8jets_py, &b_ak8jets_py);
  fChain->SetBranchAddress("ak8jets_pz", &ak8jets_pz, &b_ak8jets_pz);
  fChain->SetBranchAddress("ak8jets_e", &ak8jets_e, &b_ak8jets_e);
  fChain->SetBranchAddress("ak8jets_SoftDropMass", &ak8jets_SoftDropMass, &b_ak8jets_SoftDropMass);
  /* fChain->SetBranchAddress("ak8jets_PrunedMass", &ak8jets_PrunedMass, &b_ak8jets_PrunedMass); */
  /* fChain->SetBranchAddress("ak8jets_TrimmedMass", &ak8jets_TrimmedMass, &b_ak8jets_TrimmedMass); */
  /* fChain->SetBranchAddress("ak8jets_FilteredMass", &ak8jets_FilteredMass, &b_ak8jets_FilteredMass); */
  fChain->SetBranchAddress("ak8jets_tau1", &ak8jets_tau1, &b_ak8jets_tau1);
  fChain->SetBranchAddress("ak8jets_tau2", &ak8jets_tau2, &b_ak8jets_tau2);
  fChain->SetBranchAddress("ak8jets_tau3", &ak8jets_tau3, &b_ak8jets_tau3);
  fChain->SetBranchAddress("ak8jets_tau4", &ak8jets_tau4, &b_ak8jets_tau4);
  fChain->SetBranchAddress("ak8jets_CSV", &ak8jets_CSV, &b_ak8jets_CSV);
  fChain->SetBranchAddress("ak8jets_deepCSV_probb", &ak8jets_deepCSV_probb, &b_ak8jets_deepCSV_probb);
  fChain->SetBranchAddress("ak8jets_deepCSV_probbb", &ak8jets_deepCSV_probbb, &b_ak8jets_deepCSV_probbb);
  /* fChain->SetBranchAddress("ak8jets_deepFlavor_probb", &ak8jets_deepFlavor_probb, &b_ak8jets_deepFlavor_probb); */
  /* fChain->SetBranchAddress("ak8jets_deepFlavor_probbb", &ak8jets_deepFlavor_probbb, &b_ak8jets_deepFlavor_probbb); */
  /* fChain->SetBranchAddress("ak8jets_deepFlavor_problepb", &ak8jets_deepFlavor_problepb, &b_ak8jets_deepFlavor_problepb); */
  fChain->SetBranchAddress("ak8jets_nsubjets", &ak8jets_nsubjets, &b_ak8jets_nsubjets);
  fChain->SetBranchAddress("subjets_px", &subjets_px, &b_subjets_px);
  fChain->SetBranchAddress("subjets_py", &subjets_py, &b_subjets_py);
  fChain->SetBranchAddress("subjets_pz", &subjets_pz, &b_subjets_pz);
  fChain->SetBranchAddress("subjets_e", &subjets_e, &b_subjets_e);
  fChain->SetBranchAddress("subjets_CSV", &subjets_CSV, &b_subjets_CSV);
  fChain->SetBranchAddress("subjets_deepCSV_probb", &subjets_deepCSV_probb, &b_subjets_deepCSV_probb);
  fChain->SetBranchAddress("subjets_deepCSV_probbb", &subjets_deepCSV_probbb, &b_subjets_deepCSV_probbb);
  /* fChain->SetBranchAddress("subjets_deepFlavor_probb", &subjets_deepFlavor_probb, &b_subjets_deepFlavor_probb); */
  /* fChain->SetBranchAddress("subjets_deepFlavor_probbb", &subjets_deepFlavor_probbb, &b_subjets_deepFlavor_probbb); */
  /* fChain->SetBranchAddress("subjets_deepFlavor_problepb", &subjets_deepFlavor_problepb, &b_subjets_deepFlavor_problepb); */
  fChain->SetBranchAddress("subjets_ak8MotherIdx", &subjets_ak8MotherIdx, &b_subjets_ak8MotherIdx);
  fChain->SetBranchAddress("pv_x", &pv_x, &b_pv_x);
  fChain->SetBranchAddress("pv_y", &pv_y, &b_pv_y);
  fChain->SetBranchAddress("pv_z", &pv_z, &b_pv_z);
  //  fChain->SetBranchAddress("pv_cov", &pv_cov, &b_pv_cov);
   
  fChain->SetBranchAddress("RefitPVNoBS_x", &RefitPVNoBS_x, &b_RefitPVNoBS_x);
  fChain->SetBranchAddress("RefitPVNoBS_y", &RefitPVNoBS_y, &b_RefitPVNoBS_y);
  fChain->SetBranchAddress("RefitPVNoBS_z", &RefitPVNoBS_z, &b_RefitPVNoBS_z);
   
  fChain->SetBranchAddress("RefitPVBS_x", &RefitPVBS_x, &b_RefitPVBS_x);
  fChain->SetBranchAddress("RefitPVBS_y", &RefitPVBS_y, &b_RefitPVBS_y);
  fChain->SetBranchAddress("RefitPVBS_z", &RefitPVBS_z, &b_RefitPVBS_z);
   
  fChain->SetBranchAddress("RefitPVBS_xError", &RefitPVBS_xError, &b_RefitPVBS_xError);
  fChain->SetBranchAddress("RefitPVBS_yError", &RefitPVBS_yError, &b_RefitPVBS_yError);
  fChain->SetBranchAddress("RefitPVBS_zError", &RefitPVBS_zError, &b_RefitPVBS_zError);

  fChain->SetBranchAddress("RefitPVNoBS_xError", &RefitPVNoBS_xError, &b_RefitPVNoBS_xError);
  fChain->SetBranchAddress("RefitPVNoBS_yError", &RefitPVNoBS_yError, &b_RefitPVNoBS_yError);
  fChain->SetBranchAddress("RefitPVNoBS_zError", &RefitPVNoBS_zError, &b_RefitPVNoBS_zError);

  fChain->SetBranchAddress("pvRefit_x", &RefitPVWithTracksBS_x, &b_RefitPVWithTracksBS_x);
  fChain->SetBranchAddress("pvRefit_y", &RefitPVWithTracksBS_y, &b_RefitPVWithTracksBS_y);
  fChain->SetBranchAddress("pvRefit_z", &RefitPVWithTracksBS_z, &b_RefitPVWithTracksBS_z);
  fChain->SetBranchAddress("pvRefit_xError", &RefitPVWithTracksBS_xError, &b_RefitPVWithTracksBS_xError);
  fChain->SetBranchAddress("pvRefit_yError", &RefitPVWithTracksBS_yError, &b_RefitPVWithTracksBS_yError);
  fChain->SetBranchAddress("pvRefit_zError", &RefitPVWithTracksBS_zError, &b_RefitPVWithTracksBS_zError);
  
  fChain->SetBranchAddress("pvGen_x", &pvGen_x, &b_pvGen_x);
  fChain->SetBranchAddress("pvGen_y", &pvGen_y, &b_pvGen_y);
  fChain->SetBranchAddress("pvGen_z", &pvGen_z, &b_pvGen_z);
  fChain->SetBranchAddress("isRefitPV", &isRefitPV, &b_isRefitPV);
   
  fChain->SetBranchAddress("LeptonHash", &LeptonHash, &b_LeptonHash);
  fChain->SetBranchAddress("VertexHashNoBS1", &VertexHashNoBS1, &b_VertexHashNoBS1);
  fChain->SetBranchAddress("VertexHashNoBS2", &VertexHashNoBS2, &b_VertexHashNoBS2);
  fChain->SetBranchAddress("VertexHashBS1", &VertexHashBS1, &b_VertexHashBS1);
  fChain->SetBranchAddress("VertexHashBS2", &VertexHashBS2, &b_VertexHashBS2);

  fChain->SetBranchAddress("SelectedPairs",&SelectedPairs,&b_SelectedPairs);
  fChain->SetBranchAddress("extraelec_veto", &eleveto,&b_eleveto);
  fChain->SetBranchAddress("extramuon_veto", &muonveto,&b_muonveto);
  fChain->SetBranchAddress("trg_doubletau", &trg_doubletau,&b_trg_doubletau);
  fChain->SetBranchAddress("puppimt_1", &puppimt_1,&b_puppimt_1);
  fChain->SetBranchAddress("gen_match_1", &gen_match_1,&b_gen_match_1);
  fChain->SetBranchAddress("trigweight_1", &trigweight_1,&b_trigweight_1);
  fChain->SetBranchAddress("idisoweight_1", &idisoweight_1,&b_idisoweight_1);
  fChain->SetBranchAddress("antieweight_1", &antieweight_1,&b_antieweight_1);
  fChain->SetBranchAddress("antimuweight_1", &antimuweight_1,&b_antimuweight_1);
  fChain->SetBranchAddress("puppimt_2", &puppimt_2,&b_puppimt_2);
  fChain->SetBranchAddress("gen_match_2", &gen_match_2,&b_gen_match_2);
  fChain->SetBranchAddress("trigweight_2", &trigweight_2,&b_trigweight_2);
  fChain->SetBranchAddress("idisoweight_2", &idisoweight_2,&b_idisoweight_2);
  fChain->SetBranchAddress("antieweight_2", &antieweight_2,&b_antieweight_2);
  fChain->SetBranchAddress("antimuweight_2", &antimuweight_2,&b_antimuweight_2);
  /* fChain->SetBranchAddress("pt_tt", &pt_tt,&b_pt_tt); */
  /* fChain->SetBranchAddress("pt_vis", &pt_vis,&b_pt_vis); */
  fChain->SetBranchAddress("mt_tot", &mt_tot,&b_mt_tot);
  //  fChain->SetBranchAddress("m_vis", &m_vis,&b_m_vis);
  fChain->SetBranchAddress("PUPPImet",&PUPPImet,&b_PUPPImet);
  fChain->SetBranchAddress("puppimet_ex_UnclusteredEnUp",&puppimet_ex_UnclusteredEnUp,&b_puppimet_ex_UnclusteredEnUp);
  fChain->SetBranchAddress("puppimet_ey_UnclusteredEnUp",&puppimet_ey_UnclusteredEnUp,&b_puppimet_ey_UnclusteredEnUp);
  fChain->SetBranchAddress("puppimet_ex_UnclusteredEnDown",&puppimet_ex_UnclusteredEnDown,&b_puppimet_ex_UnclusteredEnDown);
  fChain->SetBranchAddress("puppimet_ey_UnclusteredEnDown",&puppimet_ey_UnclusteredEnDown,&b_puppimet_ey_UnclusteredEnDown);
  fChain->SetBranchAddress("PUPPImetphi",&PUPPImetphi,&b_PUPPImetphi);
  fChain->SetBranchAddress("metcov00",&PFMETCov00,&b_PFMETCov00);
  fChain->SetBranchAddress("metcov01",&PFMETCov01,&b_PFMETCov01);
  fChain->SetBranchAddress("metcov10",&PFMETCov10,&b_PFMETCov10);
  fChain->SetBranchAddress("metcov11",&PFMETCov11,&b_PFMETCov11);
  fChain->SetBranchAddress("puppimetcov00",&PUPPIMETCov00,&b_PUPPIMETCov00);
  fChain->SetBranchAddress("puppimetcov01",&PUPPIMETCov01,&b_PUPPIMETCov01);
  fChain->SetBranchAddress("puppimetcov10",&PUPPIMETCov10,&b_PUPPIMETCov10);
  fChain->SetBranchAddress("puppimetcov11",&PUPPIMETCov11,&b_PUPPIMETCov11);
  fChain->SetBranchAddress("mjj", &mjj,&b_mjj);
  fChain->SetBranchAddress("jdeta", &jdeta,&b_jdeta);
  fChain->SetBranchAddress("mjjDown", &mjjDown,&b_mjjDown);
  fChain->SetBranchAddress("jdetaDown", &jdetaDown,&b_jdetaDown);
  fChain->SetBranchAddress("mjjUp", &mjjUp,&b_mjjUp);
  fChain->SetBranchAddress("jdetaUp", &jdetaUp,&b_jdetaUp);
  
  fChain->SetBranchAddress("njetingap", &njetingap,&b_njetingap);
  fChain->SetBranchAddress("njetingap20", &njetingap20,&b_njetingap20);
  fChain->SetBranchAddress("jdphi", &jdphi,&b_jdphi);
  fChain->SetBranchAddress("dijetpt", &dijetpt,&b_dijetpt);
  fChain->SetBranchAddress("dijetphi", &dijetphi,&b_dijetphi);
  //  fChain->SetBranchAddress("ptvis", &ptvis,&b_ptvis);
  fChain->SetBranchAddress("nbtag",&nbtag,&b_nbtag);
  fChain->SetBranchAddress("njets", &njets,&b_njets);
  fChain->SetBranchAddress("njetspt20", &njetspt20,&b_njetspt20);
  fChain->SetBranchAddress("njetsUp", &njetsUp,&b_njetsUp);
  fChain->SetBranchAddress("njetspt20Up", &njetspt20Up,&b_njetspt20Up);
  fChain->SetBranchAddress("njetsDown", &njetsDown,&b_njetsDown);
  fChain->SetBranchAddress("njetspt20Down", &njetspt20Down,&b_njetspt20Down);
  fChain->SetBranchAddress("jpt_1", &jpt_1,&b_jpt_1);
  fChain->SetBranchAddress("jptDown_1", &jptDown_1,&b_jptDown_1);
  fChain->SetBranchAddress("jptUp_1", &jptUp_1,&b_jptUp_1);
  
  fChain->SetBranchAddress("jeta_1", &jeta_1,&b_jeta_1);
  fChain->SetBranchAddress("jphi_1", &jphi_1,&b_jphi_1);
  fChain->SetBranchAddress("jcsv_1", &jcsv_1,&b_jcsv_1);
  fChain->SetBranchAddress("jpt_2", &jpt_2,&b_jpt_2);
  fChain->SetBranchAddress("jeta_2", &jeta_2,&b_jeta_2);
  fChain->SetBranchAddress("jphi_2", &jphi_2,&b_jphi_2);
  fChain->SetBranchAddress("jcsv_2", &jcsv_2,&b_jcsv_2);
  fChain->SetBranchAddress("bpt_1", &bpt_1,&b_bpt_1);
  fChain->SetBranchAddress("beta_1", &beta_1,&b_beta_1);
  fChain->SetBranchAddress("bphi_1", &bphi_1,&b_bphi_1);
  fChain->SetBranchAddress("bcsv_1", &bcsv_1,&b_bcsv_1);
  fChain->SetBranchAddress("bpt_2", &bpt_2,&b_bpt_2);
  fChain->SetBranchAddress("beta_2", &beta_2,&b_beta_2);
  fChain->SetBranchAddress("bphi_2", &bphi_2,&b_bphi_2);
  fChain->SetBranchAddress("bcsv_2 ", &bcsv_2,&b_bcsv_2);
  fChain->SetBranchAddress("puweight", &puweight,&b_puweight);
  fChain->SetBranchAddress("weight", &weight,&b_weight);
  /* fChain->SetBranchAddress("jpfid_1", &jpfid_1,&b_jpfid_1); */
  /* fChain->SetBranchAddress("jpuid_1", &jpuid_1,&b_jpuid_1); */
  /* fChain->SetBranchAddress("jpfid_2", &jpfid_2,&b_jpfid_2); */
  /* fChain->SetBranchAddress("jpuid_2", &jpuid_2,&b_jpuid_2); */
  /* fChain->SetBranchAddress("bpfid_1", &bpfid_1,&b_bpfid_1); */
  /* fChain->SetBranchAddress("bpuid_1", &bpuid_1,&b_bpuid_1); */
  /* fChain->SetBranchAddress("bpfid_2", &bpfid_2,&b_bpfid_2); */
  /* fChain->SetBranchAddress("bpuid_2", &bpuid_2,&b_bpuid_2); */
  //fChain->SetBranchAddress("deepTauVsJetRaw_1", &byDeepTau2017v2p1VSjetraw_1,&b_byDeepTau2017v2p1VSjetraw_1);
  fChain->SetBranchAddress("byVVVLooseDeepTau2017v2p1VSjet_1", &byVVVLooseDeepTau2017v2p1VSjet_1,&b_byVVVLooseDeepTau2017v2p1VSjet_1);
  fChain->SetBranchAddress("byVVLooseDeepTau2017v2p1VSjet_1", &byVVLooseDeepTau2017v2p1VSjet_1,&b_byVVLooseDeepTau2017v2p1VSjet_1);
  fChain->SetBranchAddress("byVLooseDeepTau2017v2p1VSjet_1", &byVLooseDeepTau2017v2p1VSjet_1,&b_byVLooseDeepTau2017v2p1VSjet_1);
  fChain->SetBranchAddress("byLooseDeepTau2017v2p1VSjet_1", &byLooseDeepTau2017v2p1VSjet_1,&b_byLooseDeepTau2017v2p1VSjet_1);
  fChain->SetBranchAddress("byMediumDeepTau2017v2p1VSjet_1", &byMediumDeepTau2017v2p1VSjet_1,&b_byMediumDeepTau2017v2p1VSjet_1);
  fChain->SetBranchAddress("byTightDeepTau2017v2p1VSjet_1", &byTightDeepTau2017v2p1VSjet_1,&b_byTightDeepTau2017v2p1VSjet_1);
  fChain->SetBranchAddress("byVTightDeepTau2017v2p1VSjet_1", &byVTightDeepTau2017v2p1VSjet_1,&b_byVTightDeepTau2017v2p1VSjet_1);
  fChain->SetBranchAddress("byVVTightDeepTau2017v2p1VSjet_1", &byVVTightDeepTau2017v2p1VSjet_1,&b_byVVTightDeepTau2017v2p1VSjet_1);
  
  //fChain->SetBranchAddress("deepTauVsEleRaw_1", &byDeepTau2017v2p1VSeleraw_1,&b_byDeepTau2017v2p1VSeleraw_1);
  fChain->SetBranchAddress("byVVVLooseDeepTau2017v2p1VSe_1",&byVVVLooseDeepTau2017v2p1VSe_1,&b_byVVVLooseDeepTau2017v2p1VSe_1);
  fChain->SetBranchAddress("byVVLooseDeepTau2017v2p1VSe_1",&byVVLooseDeepTau2017v2p1VSe_1,&b_byVVLooseDeepTau2017v2p1VSe_1);
  fChain->SetBranchAddress("byVLooseDeepTau2017v2p1VSe_1",&byVLooseDeepTau2017v2p1VSe_1,&b_byVLooseDeepTau2017v2p1VSe_1);
  fChain->SetBranchAddress("byLooseDeepTau2017v2p1VSe_1",&byLooseDeepTau2017v2p1VSe_1,&b_byLooseDeepTau2017v2p1VSe_1);
  fChain->SetBranchAddress("byMediumDeepTau2017v2p1VSe_1",&byMediumDeepTau2017v2p1VSe_1,&b_byMediumDeepTau2017v2p1VSe_1);
  fChain->SetBranchAddress("byTightDeepTau2017v2p1VSe_1",&byTightDeepTau2017v2p1VSe_1,&b_byTightDeepTau2017v2p1VSe_1);
  fChain->SetBranchAddress("byVTightDeepTau2017v2p1VSe_1",&byVTightDeepTau2017v2p1VSe_1,&b_byVTightDeepTau2017v2p1VSe_1);
  fChain->SetBranchAddress("byVVTightDeepTau2017v2p1VSe_1",&byVVTightDeepTau2017v2p1VSe_1,&b_byVVTightDeepTau2017v2p1VSe_1);

  //fChain->SetBranchAddress("deepTauVsMuRaw_1", &byDeepTau2017v2p1VSmuraw_1,&b_byDeepTau2017v2p1VSmuraw_1);
  //fChain->SetBranchAddress("byDeepTau2017v2p1VSmuraw_1",&byDeepTau2017v2p1VSmuraw_1,&b_byDeepTau2017v2p1VSmuraw_1);
  fChain->SetBranchAddress("byVLooseDeepTau2017v2p1VSmu_1",&byVLooseDeepTau2017v2p1VSmu_1,&b_byVLooseDeepTau2017v2p1VSmu_1);
  fChain->SetBranchAddress("byLooseDeepTau2017v2p1VSmu_1",&byLooseDeepTau2017v2p1VSmu_1,&b_byLooseDeepTau2017v2p1VSmu_1);
  fChain->SetBranchAddress("byMediumDeepTau2017v2p1VSmu_1",&byMediumDeepTau2017v2p1VSmu_1,&b_byMediumDeepTau2017v2p1VSmu_1);
  fChain->SetBranchAddress("byTightDeepTau2017v2p1VSmu_1",&byTightDeepTau2017v2p1VSmu_1,&b_byTightDeepTau2017v2p1VSmu_1);
  
  //fChain->SetBranchAddress("deepTauVsJetRaw_2", &byDeepTau2017v2p1VSjetraw_2,&b_byDeepTau2017v2p1VSjetraw_2);
  fChain->SetBranchAddress("byVVVLooseDeepTau2017v2p1VSjet_2", &byVVVLooseDeepTau2017v2p1VSjet_2,&b_byVVVLooseDeepTau2017v2p1VSjet_2);
  fChain->SetBranchAddress("byVVLooseDeepTau2017v2p1VSjet_2", &byVVLooseDeepTau2017v2p1VSjet_2,&b_byVVLooseDeepTau2017v2p1VSjet_2);
  fChain->SetBranchAddress("byVLooseDeepTau2017v2p1VSjet_2", &byVLooseDeepTau2017v2p1VSjet_2,&b_byVLooseDeepTau2017v2p1VSjet_2);
  fChain->SetBranchAddress("byLooseDeepTau2017v2p1VSjet_2", &byLooseDeepTau2017v2p1VSjet_2,&b_byLooseDeepTau2017v2p1VSjet_2);
  fChain->SetBranchAddress("byMediumDeepTau2017v2p1VSjet_2", &byMediumDeepTau2017v2p1VSjet_2,&b_byMediumDeepTau2017v2p1VSjet_2);
  fChain->SetBranchAddress("byTightDeepTau2017v2p1VSjet_2", &byTightDeepTau2017v2p1VSjet_2,&b_byTightDeepTau2017v2p1VSjet_2);
  fChain->SetBranchAddress("byVTightDeepTau2017v2p1VSjet_2", &byVTightDeepTau2017v2p1VSjet_2,&b_byVTightDeepTau2017v2p1VSjet_2);
  fChain->SetBranchAddress("byVVTightDeepTau2017v2p1VSjet_2", &byVVTightDeepTau2017v2p1VSjet_2,&b_byVVTightDeepTau2017v2p1VSjet_2);
  
  //fChain->SetBranchAddress("deepTauVsEleRaw_2", &byDeepTau2017v2p1VSeleraw_2,&b_byDeepTau2017v2p1VSeleraw_2);
  fChain->SetBranchAddress("byVVVLooseDeepTau2017v2p1VSe_2",&byVVVLooseDeepTau2017v2p1VSe_2,&b_byVVVLooseDeepTau2017v2p1VSe_2);
  fChain->SetBranchAddress("byVVLooseDeepTau2017v2p1VSe_2",&byVVLooseDeepTau2017v2p1VSe_2,&b_byVVLooseDeepTau2017v2p1VSe_2);
  fChain->SetBranchAddress("byVLooseDeepTau2017v2p1VSe_2",&byVLooseDeepTau2017v2p1VSe_2,&b_byVLooseDeepTau2017v2p1VSe_2);
  fChain->SetBranchAddress("byLooseDeepTau2017v2p1VSe_2",&byLooseDeepTau2017v2p1VSe_2,&b_byLooseDeepTau2017v2p1VSe_2);
  fChain->SetBranchAddress("byMediumDeepTau2017v2p1VSe_2",&byMediumDeepTau2017v2p1VSe_2,&b_byMediumDeepTau2017v2p1VSe_2);
  fChain->SetBranchAddress("byTightDeepTau2017v2p1VSe_2",&byTightDeepTau2017v2p1VSe_2,&b_byTightDeepTau2017v2p1VSe_2);
  fChain->SetBranchAddress("byVTightDeepTau2017v2p1VSe_2",&byVTightDeepTau2017v2p1VSe_2,&b_byVTightDeepTau2017v2p1VSe_2);
  fChain->SetBranchAddress("byVVTightDeepTau2017v2p1VSe_2",&byVVTightDeepTau2017v2p1VSe_2,&b_byVVTightDeepTau2017v2p1VSe_2);

  //fChain->SetBranchAddress("deepTauVsMuRaw_2", &byDeepTau2017v2p1VSmuraw_2,&b_byDeepTau2017v2p1VSmuraw_2);
  fChain->SetBranchAddress("byVLooseDeepTau2017v2p1VSmu_2",&byVLooseDeepTau2017v2p1VSmu_2,&b_byVLooseDeepTau2017v2p1VSmu_2);
  fChain->SetBranchAddress("byLooseDeepTau2017v2p1VSmu_2",&byLooseDeepTau2017v2p1VSmu_2,&b_byLooseDeepTau2017v2p1VSmu_2);
  fChain->SetBranchAddress("byMediumDeepTau2017v2p1VSmu_2",&byMediumDeepTau2017v2p1VSmu_2,&b_byMediumDeepTau2017v2p1VSmu_2);
  fChain->SetBranchAddress("byTightDeepTau2017v2p1VSmu_2",&byTightDeepTau2017v2p1VSmu_2,&b_byTightDeepTau2017v2p1VSmu_2);

  
  /* fChain->SetBranchAddress("pvx", &pvx,&b_pvx); */
  /* fChain->SetBranchAddress("pvy ", &pvy,&b_pvy); */
  /* fChain->SetBranchAddress("pvz", &pvz,&b_pvz); */

  fChain->SetBranchAddress("dm_1", &dm_1,&b_dm_1);
  fChain->SetBranchAddress("dmMVA_1", &dmMVA_1,&b_dmMVA_1);

  fChain->SetBranchAddress("dm_2", &dm_2,&b_dm_2);
  fChain->SetBranchAddress("dmMVA_2", &dmMVA_2,&b_dmMVA_2);
  
  /* fChain->SetBranchAddress("svx_1", &svx_1,&b_svx_1); */
  /* fChain->SetBranchAddress("svy_1", &svy_1,&b_svy_1); */
  /* fChain->SetBranchAddress("svz_1", &svz_1,&b_svz_1); */

  /* fChain->SetBranchAddress("svx_2", &svx_2,&b_svx_2); */
  /* fChain->SetBranchAddress("svy_2", &svy_2,&b_svy_2); */
  /* fChain->SetBranchAddress("svz_2", &svz_2,&b_svz_2); */

  fChain->SetBranchAddress("tau1IndexVect",&tau1IndexVect,&b_tau1IndexVect);
  fChain->SetBranchAddress("tau2IndexVect",&tau2IndexVect,&b_tau2IndexVect);
  
  Notify();
}

Bool_t NtupleReader::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void NtupleReader::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t NtupleReader::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef NtupleReader_cxx
