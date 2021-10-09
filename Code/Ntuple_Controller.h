//Ntuple_Controller.h HEADER FILE

#ifndef Ntuple_Controller_h
#define Ntuple_Controller_h


// Root include files
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TMath.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TH1.h"
#include "TBits.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMatrixT.h"
#include "TMatrixTSym.h"
#include "TVectorT.h"
#include "TSystem.h"

// Include files (C & C++ libraries)
#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <utility>      // std::pair
#include <tuple>
#include <functional>

#include "NtupleReader.h"

#include "HistoConfig.h"
#ifdef USE_TauSpinner
#include "TauSpinerInterface.h"
#endif
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
#include "TauDataFormat/TauNtuple/interface/TauDecay.h"

#ifdef USE_SVfit
#include "DataFormats/SVFitObject.h"
#include "SVFitStorage.h"
#include "SVfitProvider.h"
#endif


#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/MultiProngTauSolver.h"
#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"
#include "SimpleFits/FitSoftware/interface/TauA1NuConstrainedFitter.h"

#include "PileUp.h"
//#include "TauTriggerSFs/interface/TauTriggerSFs2017.h"
//#include "TauTriggerSFs/interface/SFProvider.h"
#include "TauIDSFs/interface/TauIDSFTool.h"
//#include "TauIDSFTool.h"
#include "RecoilCorrector.h"
#include "MEtSys.h"

#include "RooWorkspace.h"
#include "RooFunctor.h"
#include <memory>
#include <unistd.h>

// Rochester muon momentum correction
//#include "CommonFiles/rochcor2012jan22.h"

// small struct needed to allow sorting indices by some value
struct sortIdxByValue {
    bool operator()(const std::pair<int,double> &left, const std::pair<int,double> &right) {
        return left.second > right.second;
    }
};

///////////////////////////////////////////////////////////////////////////////
//*****************************************************************************
//*
//*   Class: Ntuple_Controller
//*   
//*   Purpose: The purpose of this class is to provide a interface to the
//*            Ntuple
//*
//*   Designed by: Vladimir Cherepanov
//*
//*
//*****************************************************************************
///////////////////////////////////////////////////////////////////////////////


class Ntuple_Controller{
 private:
  NtupleReader *Ntp;
  TFile *newfile;
  TTree *SkimmedTree;
  int nbytes;
  int jentry;
  int nb;
  bool copyTree;

  int currentEvent;

  
  TauIDSFTool *tauSFTool2016=new TauIDSFTool("2016Legacy","DeepTau2017v2p1VSjet","Medium",true);
  TauIDSFTool *antiEleSFTool2016=new TauIDSFTool("2016Legacy","DeepTau2017v2p1VSe","VVLoose");
  TauIDSFTool *antiMuSFTool2016=new TauIDSFTool("2016Legacy","DeepTau2017v2p1VSmu","VLoose");
  TauIDSFTool *tauEmbedSFTool2016=new TauIDSFTool("2016Legacy","DeepTau2017v2p1VSjet","Medium",true,true);
  //TauIDSFTool *antiEleEmbedSFTool2016=new TauIDSFTool("2016Legacy","DeepTau2017v2p1VSe","VVLoose",false,true);
  //TauIDSFTool *antiMuEmbedSFTool2016=new TauIDSFTool("2016Legacy","DeepTau2017v2p1VSmu","VLoose",false,true);
  
  TauIDSFTool *tauSFTool2017=new TauIDSFTool("2017ReReco","DeepTau2017v2p1VSjet","Medium",true);
  TauIDSFTool *antiEleSFTool2017=new TauIDSFTool("2017ReReco","DeepTau2017v2p1VSe","VVLoose");
  TauIDSFTool *antiMuSFTool2017=new TauIDSFTool("2017ReReco","DeepTau2017v2p1VSmu","VLoose");
  TauIDSFTool *tauEmbedSFTool2017=new TauIDSFTool("2017ReReco","DeepTau2017v2p1VSjet","Medium",true,true);
  TauIDSFTool *tauSFTool2018=new TauIDSFTool("2018ReReco","DeepTau2017v2p1VSjet","Medium",true);
  TauIDSFTool *antiEleSFTool2018=new TauIDSFTool("2018ReReco","DeepTau2017v2p1VSe","VVLoose");
  TauIDSFTool *antiMuSFTool2018=new TauIDSFTool("2018ReReco","DeepTau2017v2p1VSmu","VLoose");
  TauIDSFTool *tauEmbedSFTool2018=new TauIDSFTool("2018ReReco","DeepTau2017v2p1VSjet","Medium",true,true);

  RecoilCorrector *recoilPuppiMetCorrector2016;
  RecoilCorrector *recoilPuppiMetCorrector2017;
  RecoilCorrector *recoilPuppiMetCorrector2018;
  
  MEtSys *recoilPuppiMetShifter2016;
  MEtSys *recoilPuppiMetShifter2017;
  MEtSys *recoilPuppiMetShifter2018;

  RooWorkspace *w2016;
  RooWorkspace *w2017;
  RooWorkspace *w2018;
  

  /* TString SWorkSpace2017 = (std::string)std::getenv("workdir")+"Code/LegacyCorrectionsWorkspace/output/htt_scalefactors_legacy_2017.root"; */
  /* TFile * WorkSpace2017 = new TFile(SWorkSpace2017, "read"); */
  /* TString SWorkSpace2018 = (std::string)std::getenv("workdir")+"Code/LegacyCorrectionsWorkspace/output/htt_scalefactors_legacy_2018.root"; */
  /* TFile * WorkSpace2018 = new TFile(SWorkSpace2018, "read"); */

  //std::shared_ptr<RooWorkspace> w2017_= std::shared_ptr<RooWorkspace>((RooWorkspace*)gDirectory->Get("w"));
  //std::shared_ptr<RooWorkspace> w2018_= std::shared_ptr<RooWorkspace>((RooWorkspace*)gDirectory->Get("w"));
  
  TFile *filePUdistribution2016_data;
  TFile *filePUdistribution2016_MC;
  TFile *filePUdistribution2017_data;
  TFile *filePUdistribution2017_MC;
  TFile *filePUdistribution2018_data;
  TFile *filePUdistribution2018_MC;
  
  /* TString rfilename2017_data = (std::string)std::getenv("workdir")+"Code/CommonFiles/weights/PU/data_pileup_pudistributions_data_2017.root"; */
  /* TString rfilename2017_mc = (std::string)std::getenv("workdir")+"Code/CommonFiles/weights/PU/data_pileup_pudistributions_mc_2017.root"; */
  /* TFile *filePUdistribution2017_data = new TFile(rfilename2017_data, "read"); */
  /* TFile *filePUdistribution2017_MC = new TFile(rfilename2017_mc, "read"); */
  /* TString rfilename2018_data = (std::string)std::getenv("workdir")+"Code/CommonFiles/weights/PU/pileUp_data_Autumn18.root"; */
  /* TString rfilename2018_mc = (std::string)std::getenv("workdir")+"Code/CommonFiles/weights/PU/pileUp_MC_Autumn18.root"; */
  /* TFile *filePUdistribution2018_data = new TFile(rfilename2018_data, "read"); */
  /* TFile *filePUdistribution2018_MC = new TFile(rfilename2018_mc, "read"); */

  TFile * TES2016;
  TFile * FES2016;
  TFile * TES2017;
  TFile * FES2017;
  TFile * TES2018;
  TFile * FES2018;
  
  TH1* histTES2016;
  TGraph* histFES2016;
  TH1* histTES2017;
  TGraph* histFES2017;
  TH1* histTES2018;
  TGraph* histFES2018;

  bool cannotObtainHiggsMass; // avoid repeated printing of warning when running locally

  int EmbedID;

  // Ntuple Access Functions
  virtual void Branch_Setup(TString B_Name, int type);
  virtual void Branch_Setup(){}

  // Functions to configure objects
  virtual void ConfigureObjects(); 
  void doElectrons();
  void doPhotons();
  void doJets();
  void doMuons();
  void doTaus();
  void doMET();
  unsigned int ObjEvent;

  // helper functions for internal calculations
  void printMCDecayChain(unsigned int par, unsigned int level = 0, bool printStatus = false, bool printPt = false, bool printEtaPhi = false, bool printQCD = false);

  // Object Variables
  std::vector<TLorentzVector> electrons_default;
  std::vector<TLorentzVector> photons_default;
  std::vector<TLorentzVector> jets_default;
  std::vector<TLorentzVector> muons_default;
  std::vector<TLorentzVector> taus_default;
  TLorentzVector              met_default;
  std::vector<TLorentzVector> electrons;
  std::vector<TLorentzVector> photons;
  std::vector<TLorentzVector> jets;
  std::vector<TLorentzVector> muons;
  std::vector<TLorentzVector> taus;
  TLorentzVector              met;

  // TString flags for object corrections
  TString tauCorrection;
  TString muonCorrection;
  TString elecCorrection;
  TString jetCorrection;

  // Systematic controls variables
  int theSys;
  HistoConfig HConfig;

  // Interfaces
#ifdef USE_TauSpinner  
  TauSpinerInterface TauSpinerInt;
#endif
  HistoConfig HistoC;

  // Fit Variables
  LorentzVectorParticle               theTau;
  std::vector<LorentzVectorParticle>  daughter;
  double                              LC_chi2;
  double                              ndof;
  bool                                fitStatus;
  bool                                isInit;

  // muon correction related objects
  //  rochcor2012*   rmcor;
  std::vector<TLorentzVector> Muon_corrected_p4;
  void           CorrectMuonP4();
  bool           Muon_isCorrected;

  // helpers for SVFit
#ifdef USE_SVfit
  // create SVFitObject from standard muon and standard tau_h
  void runAndSaveSVFit_MuTauh(SVFitObject* svfObj, SVFitStorage& svFitStor, const TString& metType, unsigned muIdx, unsigned tauIdx, double scaleMu, double scaleTau, bool save = true);
  // create SVFitObject from standard muon and fully reconstructed 3prong tau
  void runAndSaveSVFit_MuTau3p(SVFitObject* svfObj, SVFitStorage& svFitStor, const TString& metType, unsigned muIdx, TLorentzVector tauLV, LorentzVectorParticle neutrino, double scaleMu, double scaleTau, bool save = true);
  // create SVFitObject from standard tau_h and standard tau_h
  void runAndSaveSVFit_TauhTauh(SVFitObject* svfObj, SVFitStorage& svFitStor, const TString& metType, unsigned tauIdx1, unsigned tauIdx2, double scaleTau1, double scaleTau2, bool save = false);
#endif

 public:
  // Constructor
  Ntuple_Controller(std::vector<TString> RootFiles);

  // Destructor
  virtual ~Ntuple_Controller() ;

  // Event initializer
  void InitEvent();

  //TauSpiner function
  double TauSpinerGet(int SpinType);
  void TauSpinerSetSignal(int signalcharge){
#ifdef USE_TauSpinner
TauSpinerInt.SetTauSignalCharge(signalcharge);
#endif
}
  
   enum beamspot{BS_x0,BS_y0,BS_z0,BS_sigmaZ,BS_dxdz,BS_dydz,BS_BeamWidthX,NBS_par};
   enum TrackQuality {
     undefQuality = -1, loose = 0, tight = 1, highPurity = 2,
     confirmed = 3, goodIterative = 4, looseSetWithPV = 5, highPuritySetWithPV = 6,
     qualitySize = 7
   };
   enum TrackPar{i_qoverp = 0, i_lambda, i_phi, i_dxy,i_dsz};

   enum GenParticleFlag{Bit_isPrompt=0,
			Bit_isDecayedLeptonHadron,
			Bit_isTauDecayProduct,
			Bit_isPromptTauDecayProduct,
			Bit_isDirectTauDecayProduct,
			Bit_isDirectPromptTauDecayProduct,
			Bit_isDirectHadronDecayProduct,
			Bit_isHardProcess,
			Bit_fromHardProcess,
			Bit_isHardProcessTauDecayProduct,
			Bit_isDirectHardProcessTauDecayProduct,
			Bit_fromHardProcessBeforeFSR,
			Bit_isFirstCopy,
			Bit_isLastCopy,
			Bit_isLastCopyBeforeFSR,
			Bit_isVBFParton};




   enum TauQualityBitMask{Bit_byLooseCombinedIsolationDeltaBetaCorr3Hits=0,
			  Bit_byMediumCombinedIsolationDeltaBetaCorr3Hits,
			  Bit_byTightCombinedIsolationDeltaBetaCorr3Hits,
			  Bit_againstMuonLoose3,
			  Bit_againstMuonTight3,
			  Bit_againstElectronVLooseMVA6,
			  Bit_againstElectronLooseMVA6,
			  Bit_againstElectronMediumMVA6,
			  Bit_againstElectronTightMVA6,
			  Bit_againstElectronVTightMVA6,
			  Bit_byVLooseIsolationMVArun2v1DBoldDMwLT,
			  Bit_byLooseIsolationMVArun2v1DBoldDMwLT,
			  Bit_byMediumIsolationMVArun2v1DBoldDMwLT,
			  Bit_byTightIsolationMVArun2v1DBoldDMwLT,
			  Bit_byVTightIsolationMVArun2v1DBoldDMwLT,
			  Bit_byVLooseIsolationMVArun2v1DBnewDMwLT,
			  Bit_byLooseIsolationMVArun2v1DBnewDMwLT,
			  Bit_byMediumIsolationMVArun2v1DBnewDMwLT,
			  Bit_byTightIsolationMVArun2v1DBnewDMwLT,
			  Bit_byVTightIsolationMVArun2v1DBnewDMwLT,
			  Bit_byLooseIsolationMVArun2v1DBdR03oldDMwLT,
			  Bit_byMediumIsolationMVArun2v1DBdR03oldDMwLT,
			  Bit_byTightIsolationMVArun2v1DBdR03oldDMwLT,
			  Bit_byVTightIsolationMVArun2v1DBdR03oldDMwLT,
			  Bit_byVLooseIsolationMVArun2017v1DBoldDMwLT2017, //FRA syncApr2018
			  Bit_byLooseIsolationMVArun2017v1DBoldDMwLT2017,  //FRA syncApr2018
			  Bit_byMediumIsolationMVArun2017v1DBoldDMwLT2017, //FRA syncApr2018
			  Bit_byTightIsolationMVArun2017v1DBoldDMwLT2017,  //FRA syncApr2018
			  Bit_byVTightIsolationMVArun2017v1DBoldDMwLT2017, //FRA syncApr2018
			  Bit_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017, //FRA syncApr2018
			  Bit_byVLooseIsolationMVArun2017v2DBoldDMwLT2017, //FRA syncApr2018
			  Bit_byLooseIsolationMVArun2017v2DBoldDMwLT2017,  //FRA syncApr2018
			  Bit_byMediumIsolationMVArun2017v2DBoldDMwLT2017, //FRA syncApr2018
			  Bit_byTightIsolationMVArun2017v2DBoldDMwLT2017,  //FRA syncApr2018
			  Bit_byVTightIsolationMVArun2017v2DBoldDMwLT2017, //FRA syncApr2018
			  Bit_byVVTightIsolationMVArun2017v2DBoldDMwLT2017, //FRA syncApr2018
			  Bit_byVLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017, //FRA syncApr2018
			  Bit_byLooseIsolationMVArun2017v2DBoldDMdR0p3wLT2017,  //FRA syncApr2018
			  Bit_byMediumIsolationMVArun2017v2DBoldDMdR0p3wLT2017, //FRA syncApr2018
			  Bit_byTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017,  //FRA syncApr2018
			  Bit_byVTightIsolationMVArun2017v2DBoldDMdR0p3wLT2017, //FRA syncApr2018
			  Bit_byVVVLooseDeepTau2017v2p1VSjet,
			  Bit_byVVLooseDeepTau2017v2p1VSjet, 
			  Bit_byVLooseDeepTau2017v2p1VSjet,  
			  Bit_byLooseDeepTau2017v2p1VSjet,   
			  Bit_byMediumDeepTau2017v2p1VSjet,  
			  Bit_byTightDeepTau2017v2p1VSjet,   
			  Bit_byVTightDeepTau2017v2p1VSjet,  
			  Bit_byVVTightDeepTau2017v2p1VSjet, 
			  Bit_byVVVLooseDeepTau2017v2p1VSe,  
			  Bit_byVVLooseDeepTau2017v2p1VSe, 
			  Bit_byVLooseDeepTau2017v2p1VSe,   
			  Bit_byLooseDeepTau2017v2p1VSe,	
			  Bit_byMediumDeepTau2017v2p1VSe,   
			  Bit_byTightDeepTau2017v2p1VSe,	
			  Bit_byVTightDeepTau2017v2p1VSe,   
			  Bit_byVVTightDeepTau2017v2p1VSe,   
			  Bit_byVLooseDeepTau2017v2p1VSmu, 
			  Bit_byLooseDeepTau2017v2p1VSmu, 
			  Bit_byMediumDeepTau2017v2p1VSmu, 
			  Bit_byTightDeepTau2017v2p1VSmu};





  enum MuonQualityBitMask{Bit_MuonLoose=0,
			 Bit_MuonSoft,
			 Bit_MuonMedium,
			 Bit_MuonTight,
			 Bit_MuonHighPt,
			 Bit_MuonTight_noVtx};



  enum particleType {
  MUON = 0,
  ELECTRON = 1,
  TAU =2
  };
  float m_MVAEleIDCuts[2][2][3] ;

  enum pairType {
    MuHad  = 0,
    EHad   = 1,
    HadHad = 2,
    MuMu   = 3,
    EE     = 4,
    EMu    = 5,
    EEPrompt = 6, // prompt Z->ee/mumu decays
    MuMuPrompt = 7,
    Other  = 8 // for e.g. h->bb
  };

  enum eleMVAIDWP {
    EMVATight = 0, // 80% eff
    EMVALoose = 1  // 90% eff
  };

  enum muIDWP {
    MuLoose  = 0,
    MuSoft   = 1,
    MuMedium = 2,
    MuTight  = 3,
    MuHighPt = 4
  };

  enum aeleWP {
    aeleVLoose = 0,
    aeleLoose  = 1,
    aeleMedium = 2,
    aeleTight  = 3,
    aeleVTight = 4
  };

  enum amuWP {
    amuLoose = 0,
    amuTight = 1
  };

  typedef std::vector<float> tauPair_t; // pt1 - iso1 - idx1 - pt2 - iso2 - idx2 - idxoriginalPair
  typedef std::tuple <float, float, int, float, float, int, int> tauPair_tuple; // pt1 - iso1 - idx1 - pt2 - iso2 - idx2 - idxoriginalPair
  // access to SVFit
  #ifdef USE_SVfit
  SVFitObject* getSVFitResult_MuTauh(SVFitStorage& svFitStor, TString metType, unsigned muIdx, unsigned tauIdx, unsigned rerunEvery = 5000, TString suffix = "", double scaleMu = 1 , double scaleTau = 1);
  SVFitObject* getSVFitResult_MuTau3p(SVFitStorage& svFitStor, TString metType, unsigned muIdx, TLorentzVector tauLV, LorentzVectorParticle neutrino, TString suffix = "", double scaleMu = 1, double scaleTau = 1);
  SVFitObject* getSVFitResult_TauhTauh(SVFitStorage& svFitStor, TString metType, unsigned tauIdx1, unsigned tauIdx2, unsigned rerunEvery  = 5000 , TString suffix  ="" , double scaleTau1  =1 , double scaleTau2  =1 );


  #endif


  // Ntuple Access Functions
  virtual Int_t Get_Entries();
  virtual void Get_Event(int _jentry);
  virtual Int_t Get_EventIndex();
  virtual TString Get_File_Name();

  //Ntuple Cloning Functions
  virtual void CloneTree(TString n);
  virtual void SaveCloneTree();
  inline void AddEventToCloneTree(){if(copyTree)SkimmedTree->Fill();}

  // Systematic controls
  enum    Systematic {Default=0,NSystematics};

  int     SetupSystematics(TString sys_);
  void    SetSysID(int sysid){theSys=sysid;}


  // Data/MC switch and thin
  TH1F* hLLRCounters;
  bool isData()  {return (bool)Ntp->Event_isRealData;}
  void ThinTree();

  // Set object corrections to be applied
  void SetTauCorrections(TString tauCorr){tauCorrection = tauCorr;}
  void SetMuonCorrections(TString muonCorr){muonCorrection = muonCorr;}
  void SetElecCorrections(TString elecCorr){elecCorrection = elecCorr;}
  void SetJetCorrections(TString jetCorr){jetCorrection = jetCorr;}
  // corresponding getters
  const TString& GetTauCorrections() const {return tauCorrection;}
  const TString& GetMuonCorrections() const {return muonCorrection;}
  const TString& GetElecCorrections() const {return elecCorrection;}
  const TString& GetJetCorrections() const {return jetCorrection;}

  // Information from input Ntuple path name
  TString GetInputNtuplePath();
  TString GetInputDatasetName();
  TString GetInputPublishDataName();
  int getSampleHiggsMass();
  int readHiggsMassFromString(TString input);

  // resonance mass
  int getHiggsSampleMassFromGenInfo();
  double getResonanceMassFromGenInfo(bool useZ0 = true, bool useHiggs0 = true, bool useW = true);

  // Physics Variable Get Functions
  // Event Variables
  Long64_t GetMCID();
  int GetStrippedMCID();
  /* unsigned int RunNumber(){return Ntp->Event_RunNumber;} */
  /* unsigned int EventNumber(){ return Ntp->Event_EventNumber;} */
  /* int BunchCrossing(){ return Ntp->Event_bunchCrossing;} */
  /* int OrbitNumber(){ return Ntp->Event_orbitNumber;} */
  /* unsigned int LuminosityBlock(){return Ntp->Event_luminosityBlock;} */
  /* float           PileupInfo_TrueNumInteractions_nm1(){return Ntp->PileupInfo_TrueNumInteractions_nm1;} */
  /* float           PileupInfo_TrueNumInteractions_n0(){return Ntp->PileupInfo_TrueNumInteractions_n0;} */
  /* float           PileupInfo_TrueNumInteractions_np1(){return Ntp->PileupInfo_TrueNumInteractions_np1;} */
  /* float        PUWeight(){return Ntp->PUWeight;} */
  /* float        PUWeight_p5(){return Ntp->PUWeight_p5;} */
  /* float        PUWeight3D_m5(){return Ntp->PUWeight3D_m5;} */
  /* float        PUWeight3D(){return Ntp->PUWeight3D;} */
  /* float        PUWeight3D_p5(){return Ntp->PUWeight3D_p5;} */
  /* float        PUWeight_m5(){return Ntp->PUWeight_m5;} */
  /* float		   PUWeightFineBins(){return Ntp->PUWeightFineBins;} */

  // embedding
  /* float		Embedding_TauSpinnerWeight(){return Ntp->TauSpinnerWeight;}; */
  /* float 	Embedding_SelEffWeight(){return Ntp->SelEffWeight;} */
  /* float 	Embedding_MinVisPtFilter(){return Ntp->MinVisPtFilter;} */
  /* float 	Embedding_KinWeightPt(){return Ntp->KinWeightPt;} */
  /* float 	Embedding_KinWeightEta(){return Ntp->KinWeightEta;} */
  /* float 	Embedding_KinWeightMassPt(){return Ntp->KinWeightMassPt;} */
  // don't use this combined weight blindly. Check which weights you should apply.
  /* float		EmbeddedWeight(){ */
  /* 	  return Ntp->TauSpinnerWeight * Ntp->SelEffWeight * Ntp->MinVisPtFilter * Ntp->KinWeightPt * Ntp->KinWeightEta * Ntp->KinWeightMassPt; */
  /* } */

  /* TVectorT<double>      beamspot_par(){TVectorT<double> BS(NBS_par);for(unsigned int i=0;i<NBS_par;i++)BS(i)=Ntp->beamspot_par->at(i);return BS;} */

  /* TMatrixTSym<double>   beamspot_cov(){ */
  /*   TMatrixTSym<double> BS_cov(NBS_par); */
  /*   unsigned int l=0; */
  /*   for(unsigned int i=0;i<NBS_par;i++){ */
  /*     for(unsigned int j=i;j<NBS_par;j++){ */
  /* 	BS_cov(i,j)=Ntp->beamspot_cov->at(l); */
  /*     } */
  /*   } */
  /*   return BS_cov; */
  /* } */
  
  /* double  beamspot_emittanceX(){return Ntp->beamspot_emittanceX;} */
  /* double  beamspot_emittanceY(){return Ntp->beamspot_emittanceY;} */
  /* double  beamspot_betaStar(){return Ntp->beamspot_betaStar;} */



  /* // Vertex Information */
   unsigned int NVtx(){return Ntp->npv;}
   TVector3       PVtx(){return TVector3(Ntp->pv_x,Ntp->pv_y,Ntp->pv_z);}
   double pv_z(){return Ntp->pv_z;}
   TMatrixTSym<float> PFTau_TIP_primaryVertex_cov();
   bool isPVCovAvailable();
   /* bool           isPVtxRefit(){return Ntp->isRefitPV;} */
   /* TVector3       PVRefitNoBSOld(){return TVector3(Ntp->pvRefitNoBS_x,Ntp->pvRefitNoBS_y,Ntp->pvRefitNoBS_z);} */
   /* TVector3       PVRefitBSOld(){return TVector3(Ntp->pvRefitBS_x,Ntp->pvRefitBS_y,Ntp->pvRefitBS_z);} */
   /* TVector3       PVRefitNoBSZNominalOld(){return TVector3(Ntp->pvRefitNoBS_x,Ntp->pvRefitNoBS_y,Ntp->pv_z);} */
   /* TVector3       PVRefitBSZNominalOld(){return TVector3(Ntp->pvRefitBS_x,Ntp->pvRefitBS_y,Ntp->pv_z);} */
   /* unsigned int   NPVRefitNoBSOld(){return Ntp->pvRefitNoBS_x->size();} */
   /* unsigned int   NPVRefitBSOld(){return Ntp->pvRefitBS_x->size();} */
   /* bool           isPVRefitNoBSOld(){return (Ntp->pvRefitNoBS_x->size()>0);} */
   /* bool           isPVRefitBSOld(){return (Ntp->pvRefitBS_x->size()>0);} */
   
   /* TVector3       PVRefitNoBSOneTrackRemovedOld(unsigned int i){return TVector3(Ntp->pvRefitNoBSTracksRemoved_x->at(i),Ntp->pvRefitNoBSTracksRemoved_y->at(i),Ntp->pvRefitNoBSTracksRemoved_z->at(i));} */
   /* TVector3       PVRefitBSOneTrackRemovedOld(unsigned int i){return TVector3(Ntp->pvRefitBSTracksRemoved_x->at(i),Ntp->pvRefitBSTracksRemoved_y->at(i),Ntp->pvRefitBSTracksRemoved_z->at(i));} */
   /* TVector3       PVRefitNoBSOneTrackRemovedZNominalOld(unsigned int i){return TVector3(Ntp->pvRefitNoBSTracksRemoved_x->at(i),Ntp->pvRefitNoBSTracksRemoved_y->at(i),Ntp->pv_z);} */
   /* TVector3       PVRefitBSOneTrackRemovedZNominalOld(unsigned int i){return TVector3(Ntp->pvRefitBSTracksRemoved_x->at(i),Ntp->pvRefitBSTracksRemoved_y->at(i),Ntp->pv_z);} */
   /* unsigned int   NPVRefitNoBSOneTrackRemovedOld(){return Ntp->pvRefitNoBSTracksRemoved_x->size();} */
   /* unsigned int   NPVRefitBSOneTrackRemovedOld(){return Ntp->pvRefitBSTracksRemoved_x->size();} */
   /* bool           isPVRefitNoBSOneTrackRemovedOld(){return (Ntp->pvRefitNoBSTracksRemoved_x->size()>0);} */
   /* bool           isPVRefitBSOneTrackRemovedOld(){return (Ntp->pvRefitBSTracksRemoved_x->size()>0);} */

   /* TVector3       PVRefitNoBSTracksRemovedOld(unsigned int k){return TVector3(Ntp->RefitPVNoBSNew_x->at(k),Ntp->RefitPVNoBSNew_y->at(k),Ntp->RefitPVNoBSNew_z->at(k));} */
   /* TVector3       PVRefitBSTracksRemovedOld(unsigned int k){return TVector3(Ntp->RefitPVBSNew_x->at(k),Ntp->RefitPVBSNew_y->at(k),Ntp->RefitPVBSNew_z->at(k));} */
   /* TVector3       PVRefitNoBSTracksRemovedZNominalOld(unsigned int k){return TVector3(Ntp->RefitPVNoBSNew_x->at(k),Ntp->RefitPVNoBSNew_y->at(k),Ntp->pv_z);} */
   /* TVector3       PVRefitBSTracksRemovedZNominalOld(unsigned int k){return TVector3(Ntp->RefitPVBSNew_x->at(k),Ntp->RefitPVBSNew_y->at(k),Ntp->pv_z);} */
   /* double RefitPVNoBSTracksRemovedOld_x(unsigned int i){return Ntp->RefitPVNoBSNew_x->at(i);} */
   /* double RefitPVNoBSTracksRemovedOld_y(unsigned int i){return Ntp->RefitPVNoBSNew_y->at(i);} */
   /* double RefitPVNoBSTracksRemovedOld_z(unsigned int i){return Ntp->RefitPVNoBSNew_z->at(i);} */
   /* double RefitPVBSTracksRemovedOld_x(unsigned int i){return Ntp->RefitPVBSNew_x->at(i);} */
   /* double RefitPVBSTracksRemovedOld_y(unsigned int i){return Ntp->RefitPVBSNew_y->at(i);} */
   /* double RefitPVBSTracksRemovedOld_z(unsigned int i){return Ntp->RefitPVBSNew_z->at(i);} */
   /* unsigned int   NPVRefitNoBSTracksRemovedOld(){return Ntp->RefitPVNoBSNew_x->size();} */
   /* unsigned int   NPVRefitBSTracksRemovedOld(){return Ntp->RefitPVBSNew_x->size();} */
   /* bool           isPVRefitNoBSTracksRemovedOld(){return (Ntp->RefitPVNoBSNew_x->size()>0);} */
   /* bool           isPVRefitBSTracksRemovedOld(){return (Ntp->RefitPVBSNew_x->size()>0);} */
   
   /* TVector3       PVRefitNoBSNew(unsigned int k){return TVector3(Ntp->RefitPVNoBS_x->at(k),Ntp->RefitPVNoBS_y->at(k),Ntp->RefitPVNoBS_z->at(k));} */
   /* TVector3       PVRefitBSNew(unsigned int k){return TVector3(Ntp->RefitPVBS_x->at(k),Ntp->RefitPVBS_y->at(k),Ntp->RefitPVBS_z->at(k));} */
   /* TVector3       PVRefitNoBSZNominalNew(unsigned int k){return TVector3(Ntp->RefitPVNoBS_x->at(k),Ntp->RefitPVNoBS_y->at(k),Ntp->pv_z);} */
   /* TVector3       PVRefitBSZNominalNew(unsigned int k){return TVector3(Ntp->RefitPVBS_x->at(k),Ntp->RefitPVBS_y->at(k),Ntp->pv_z);} */
   double RefitPVNoBS_x(unsigned int k){return Ntp->RefitPVNoBS_x->at(k);}
   double RefitPVNoBS_y(unsigned int k){return Ntp->RefitPVNoBS_y->at(k);}
   double RefitPVNoBS_z(unsigned int k){return Ntp->RefitPVNoBS_z->at(k);}
   double RefitPVBS_x(unsigned int k){return Ntp->RefitPVBS_x->at(k);}
   double RefitPVBS_y(unsigned int k){return Ntp->RefitPVBS_y->at(k);}
   double RefitPVBS_z(unsigned int k){return Ntp->RefitPVBS_z->at(k);}

   double RefitPVWithTracksBS_x(){return Ntp->RefitPVWithTracksBS_x;}
   double RefitPVWithTracksBS_y(){return Ntp->RefitPVWithTracksBS_y;}
   double RefitPVWithTracksBS_z(){return Ntp->RefitPVWithTracksBS_z;}
   double RefitPVWithTracksBS_xError(){return Ntp->RefitPVWithTracksBS_xError;}
   double RefitPVWithTracksBS_yError(){return Ntp->RefitPVWithTracksBS_yError;}
   double RefitPVWithTracksBS_zError(){return Ntp->RefitPVWithTracksBS_zError;}   
   
   unsigned int   NPVRefitNoBS(){return Ntp->RefitPVNoBS_x->size();}
   unsigned int   NPVRefitBS(){return Ntp->RefitPVBS_x->size();}
   bool           isPVRefitNoBS(){return (Ntp->RefitPVNoBS_x->size()>0);}
   bool           isPVRefitBS(){return (Ntp->RefitPVBS_x->size()>0);}

   /* TVector3       PVRefitNoBSTIP(unsigned int i){return TVector3(Ntp->PFTau_TIP_PVPosNoBS->at(i).at(0),Ntp->PFTau_TIP_PVPosNoBS->at(i).at(1),Ntp->PFTau_TIP_PVPosNoBS->at(i).at(2));} */
   /* TVector3       PVRefitBSTIP(unsigned int i){return TVector3(Ntp->PFTau_TIP_PVPosBS->at(i).at(0),Ntp->PFTau_TIP_PVPosBS->at(i).at(1),Ntp->PFTau_TIP_PVPosBS->at(i).at(2));} */
   /* TVector3       PVRefitNoBSZNominalTIP(unsigned int i){return TVector3(Ntp->PFTau_TIP_PVPosNoBS->at(i).at(0),Ntp->PFTau_TIP_PVPosNoBS->at(i).at(1),Ntp->pv_z);} */
   /* TVector3       PVRefitBSZNominalTIP(unsigned int i){return TVector3(Ntp->PFTau_TIP_PVPosBS->at(i).at(0),Ntp->PFTau_TIP_PVPosBS->at(i).at(1),Ntp->pv_z);} */
   /* double PVRefitTIPNoBS_Size(unsigned int i){return  Ntp->PFTau_TIP_PVPosNoBS->at(i).size();} */
   /* double PVRefitTIPBS_Size(unsigned int i){return  Ntp->PFTau_TIP_PVPosBS->at(i).size();} */
   /* bool           isPVRefitNoBSNew(unsigned int i){return (Ntp->PFTau_TIP_PVPosNoBS_x->at(i).size()>0);} */
   /* bool           isPVRefitBSNew(unsigned int i){return (Ntp->PFTau_TIP_PVPosBS_x->at(i).size()>0);} */
   
   TVector3       PVtx_Gen(){return TVector3(Ntp->pvGen_x,Ntp->pvGen_y,Ntp->pvGen_z);}

   unsigned int NRefitVtx(){return Ntp->RefitPVNoBS_x->size();}
   //bool         isPVtxRefit(){return (Ntp->RefitPVNoBS_x->size()>0);}
   //TVector3 PVRefitNoBS(unsigned int k){return TVector3(Ntp->RefitPVNoBS_x->at(k),Ntp->RefitPVNoBS_y->at(k),Ntp->RefitPVNoBS_z->at(k));}
   //TVector3 PVRefitBS(unsigned int k){return TVector3(Ntp->RefitPVBS_x->at(k),Ntp->RefitPVBS_y->at(k),Ntp->RefitPVBS_z->at(k));}
   //TVector3 PFTau_TIP_primaryVertex_pos(unsigned int i){return  TVector3(Ntp->PFTau_TIP_PVPos->at(i).at(0),Ntp->PFTau_TIP_PVPos->at(i).at(1),Ntp->PFTau_TIP_PVPos->at(i).at(2));}
   
   //double PFTau_TIP_primaryVertex_pos_Size(unsigned int i){return  Ntp->PFTau_TIP_PVPos->at(i).size();}

   size_t LeptonHash(unsigned int k){return Ntp->LeptonHash->at(k);}
   //vector<size_t> LeptonHash(){return Ntp->LeptonHash;}
   //unsigned int NLeptonHash(){return Ntp->LeptonHash->size();}
   size_t VertexHashNoBS1(unsigned int k){return Ntp->VertexHashNoBS1->at(k);}
   size_t VertexHashNoBS2(unsigned int k){return Ntp->VertexHashNoBS2->at(k);}
   unsigned int NVertexHashNoBS(){return Ntp->VertexHashNoBS1->size();}
   size_t VertexHashBS1(unsigned int k){return Ntp->VertexHashBS1->at(k);}
   size_t VertexHashBS2(unsigned int k){return Ntp->VertexHashBS2->at(k);}
   unsigned int NVertexHashBS(){return Ntp->VertexHashBS1->size();}
   
   /* size_t VertexHashNoBSTracksRemovedOld1(unsigned int k){return Ntp->VertexHashNoBSNew1->at(k);} */
   /* size_t VertexHashNoBSTracksRemovedOld2(unsigned int k){return Ntp->VertexHashNoBSNew2->at(k);} */
   /* unsigned int NVertexHashNoBSTracksRemovedOld(){return Ntp->VertexHashNoBSNew1->size();} */
   /* size_t VertexHashBSTracksRemovedOld1(unsigned int k){return Ntp->VertexHashBSNew1->at(k);} */
   /* size_t VertexHashBSTracksRemovedOld2(unsigned int k){return Ntp->VertexHashBSNew2->at(k);} */
   /* unsigned int NVertexHashBSTracksRemovedOld(){return Ntp->VertexHashBSNew1->size();} */

   unsigned int LeptonHashSize(){return Ntp->LeptonHash->size();}
   
   /* unsigned int VertexHashNoBS1Size(){return Ntp->VertexHashNoBS1->size();} */
   /* unsigned int VertexHashNoBS2Size(){return Ntp->VertexHashNoBS2->size();} */
   /* unsigned int VertexHashBS1Size(){return Ntp->VertexHashBS1->size();} */
   /* unsigned int VertexHashBS2Size(){return Ntp->VertexHashBS2->size();} */

   /* unsigned int VertexHashNoBSNew1Size(){return Ntp->VertexHashNoBSNew1->size();} */
   /* unsigned int VertexHashNoBSNew2Size(){return Ntp->VertexHashNoBSNew2->size();} */
   /* unsigned int VertexHashBSNew1Size(){return Ntp->VertexHashBSNew1->size();} */
   /* unsigned int VertexHashBSNew2Size(){return Ntp->VertexHashBSNew2->size();} */
 
  /* double       Vtx_chi2(unsigned int i){return Ntp->Vtx_chi2->at(i);} */
  /* unsigned     Vtx_nTrk(unsigned int i){return Ntp->Vtx_nTrk->at(i);} */
  /* float        Vtx_ndof(unsigned int i){return Ntp->Vtx_ndof->at(i);} */
  /* TMatrixF     Vtx_Cov(unsigned int i); */
  /* std::vector<int>  Vtx_Track_idx(unsigned int i){return Ntp->Vtx_Track_idx->at(i);} */
  /* bool Vtx_isFake(unsigned int i){return Ntp->Vtx_isFake->at(i);} */
  /* TLorentzVector Vtx_TracksP4(unsigned int i, unsigned int j){return TLorentzVector(Ntp->Vtx_TracksP4->at(i).at(j).at(1),Ntp->Vtx_TracksP4->at(i).at(j).at(2),Ntp->Vtx_TracksP4->at(i).at(j).at(3),Ntp->Vtx_TracksP4->at(i).at(j).at(0));} */

  /* bool isVtxGood(unsigned int i); */
  /* bool isGoodVtx(unsigned int i); */


   ULong64_t EventNumber(){return Ntp->EventNumber;}
   Int_t RunNumber(){return Ntp->RunNumber;}
   Int_t           LuminosityBlock(){return Ntp->lumi;}
   //   Int_t           NBadMu(){return Ntp->NBadMu;}
   Long64_t        triggerbit(){return Ntp->triggerbit;}
   Long64_t           metfilterbit(){return Ntp->metfilterbit;}
   //Float_t         MET(){return Ntp->met;}
   //Float_t         METphi(){return Ntp->metphi;}
   //Float_t         PUPPImet(){return Ntp->PUPPImet;}
   //Float_t         PUPPImetphi(){return Ntp->PUPPImetphi;}
   //Float_t         PFMETCov00(){return Ntp->PFMETCov00;}
   //Float_t         PFMETCov01(){return Ntp->PFMETCov01;}
   //Float_t         PFMETCov10(){return Ntp->PFMETCov10;}
   //Float_t         PFMETCov11(){return Ntp->PFMETCov11;}
   Float_t         PFMETsignif(){return Ntp->PFMETsignif;}
   Float_t         npu(){return Ntp->npu;}
   Int_t           npv(){return Ntp->npv;}
   //Float_t         PUReweight(){return Ntp->PUReweight;}
   double         PUReweight();
   double         prefiringweight(){return Ntp->prefiringweight;}
   double         prefiringweightup(){return Ntp->prefiringweightup;}
   double         prefiringweightdown(){return Ntp->prefiringweightdown;}
   Int_t           PUNumInteractions(){return Ntp->PUNumInteractions;}
   Float_t         rho(){return Ntp->rho;}


   int getBitOfGivenTrigger(TString tname);


   unsigned int  NMothers(){return Ntp->mothers_px->size();}

   TLorentzVector Mothers_P4(unsigned int i){return TLorentzVector(Ntp->mothers_px->at(i), Ntp->mothers_py->at(i), Ntp->mothers_pz->at(i),Ntp->mothers_e->at(i));}
   Long64_t mothers_trgSeparateMatch(unsigned int i){return Ntp->mothers_trgSeparateMatch->at(i);}


   int NTriggers(){return Ntp->trigger_accept->size();}
   bool TriggerAccept(unsigned int i){return Ntp->trigger_accept->at(i);}
   TString TriggerName(unsigned int i){return Ntp->trigger_name->at(i);}
   bool         GetTriggerIndex(TString n,  int &i);
   std::vector<int> GetVectorTriggers(TString n);
   std::vector<int> GetVectorTriggers(std::vector<TString> v);
   std::vector<int> GetVectorTriggersFullMatch(std::vector<TString> v);
   std::vector<int> GetVectorCrossTriggers(TString n1,TString n2,TString f1,TString f2);
   bool  CheckIfAnyPassed(  std::vector<int> list);
   bool  passecalBadCalibFilterUpdate(){return Ntp->passecalBadCalibFilterUpdate;}
 
   Float_t         MC_weight(){return Ntp->MC_weight;}
   Float_t         MC_weight_scale_muF0p5(){return Ntp->MC_weight_scale_muF0p5;}
   Float_t         MC_weight_scale_muF2(){return Ntp->MC_weight_scale_muF2;}
   Float_t         MC_weight_scale_muR0p5(){return Ntp->MC_weight_scale_muR0p5;}
   Float_t         MC_weight_scale_muR2(){return Ntp->MC_weight_scale_muR2;}
   Double_t        nominal_wt(){return Ntp->nominal_wt;}
   
   int             TheoreticalPSUncSize(){return Ntp->TheoreticalPSUnc->size();}
   Double_t        TheoreticalPSUnc(int t){return Ntp->TheoreticalPSUnc->at(t);}
   Double_t        TheoreticalScaleUnc1005(){return Ntp->TheoreticalScaleUnc1005;}
   Double_t        TheoreticalScaleUnc1009(){return Ntp->TheoreticalScaleUnc1009;}
   Double_t        TheoreticalScaleUnc5(){return Ntp->TheoreticalScaleUnc5;}
   Double_t        TheoreticalScaleUnc9(){return Ntp->TheoreticalScaleUnc9;}
   
   Float_t         lheHt(){return Ntp->lheHt;}
   Int_t            lheNOutPartons(){return Ntp->lheNOutPartons;}
   Int_t            lheNOutB(){return Ntp->lheNOutB;}
   Int_t            lheNOutC(){return Ntp->lheNOutC;}
   //Float_t         aMCatNLOweight(){return Ntp->aMCatNLOweight;}
   Int_t            DataMC_Type(){return Ntp->DataMC_Type_idx;}


   unsigned int    NGenParts(){return Ntp->genpart_px->size();}
   TLorentzVector Genpart_P4(unsigned int i){return TLorentzVector(Ntp->genpart_px->at(i), Ntp->genpart_py->at(i), Ntp->genpart_pz->at(i),Ntp->genpart_e->at(i));}
   int Genpart_pdg(unsigned int i){return Ntp->genpart_pdg->at(i);}
   int Genpart_status(unsigned int i){return Ntp->genpart_status->at(i);}
   int Genpart_HMothInd(unsigned int i){return Ntp->genpart_HMothInd->at(i);}
   int Genpart_MSSMHMothInd(unsigned int i){return Ntp->genpart_MSSMHMothInd->at(i);}
   int Genpart_TopMothInd(unsigned int i){return Ntp->genpart_TopMothInd->at(i);}
   int Genpart_TauMothInd(unsigned int i){return Ntp->genpart_TauMothInd->at(i);}
   int Genpart_ZMothInd(unsigned int i){return Ntp->genpart_ZMothInd->at(i);}
   int Genpart_WMothInd(unsigned int i){return Ntp->genpart_WMothInd->at(i);}
   int Genpart_bMothInd(unsigned int i){return Ntp->genpart_bMothInd->at(i);}
   int Genpart_HZDecayMode(unsigned int i){return Ntp->genpart_HZDecayMode->at(i);}
   int Genpart_TopDecayMode(unsigned int i){return Ntp->genpart_TopDecayMode->at(i);}
   int Genpart_WDecayMode(unsigned int i){return Ntp->genpart_WDecayMode->at(i);}
   int Genpart_TauGenDecayMode(unsigned int i){return Ntp->genpart_TauGenDecayMode->at(i);}
   int Genpart_TauGenDetailedDecayMode(unsigned int i){return Ntp->genpart_TauGenDetailedDecayMode->at(i);}
   int Genpart_flags(unsigned int i){return Ntp->genpart_flags->at(i);}


   int NGenJets(){return Ntp->genjet_px->size();}
   TLorentzVector GenJet_P4(unsigned int i){return TLorentzVector(Ntp->genjet_px->at(i), Ntp->genjet_py->at(i), Ntp->genjet_pz->at(i),Ntp->genjet_e->at(i));}
   int Genjet_partonFlavour(unsigned int i){return Ntp->genjet_partonFlavour->at(i);}
   int Genjet_hadronFlavour(unsigned int i){return Ntp->genjet_hadronFlavour->at(i);}
   TLorentzVector TauP4_Corrected(unsigned int i);
   TLorentzVector P4Corrected(unsigned int i,int genmatch,string Unc="Nom");
   int GetGenMatch(int tauindex);
 
   Int_t NUP(){return Ntp->NUP;}
   /* bool isSVFitInfoAvailable(){if(Ntp->SVfit_fitMETPhiTauUp->size()!=0) return true; return false;} */
   /* Int_t NPairCandidate(){return Ntp->SVfit_fitMETPhiTauUp->size();} */
   /* float  SVfit_fitMETPhiTauUp(unsigned int i){return Ntp->SVfit_fitMETPhiTauUp->at(i);} */
   /* float  SVfit_fitMETPhiTauDown(unsigned int i){return Ntp->SVfit_fitMETPhiTauDown->at(i);} */
   /* float  SVfit_fitMETRhoTauUp(unsigned int i){return Ntp->SVfit_fitMETRhoTauUp->at(i);} */
   /* float  SVfit_fitMETRhoTauDown(unsigned int i){return Ntp->SVfit_fitMETRhoTauDown->at(i);} */
   /* float  SVfit_phiUncTauUp(unsigned int i){return Ntp->SVfit_phiUncTauUp->at(i);} */
   /* float  SVfit_phiUncTauDown(unsigned int i){return Ntp->SVfit_phiUncTauDown->at(i);} */
   /* float  SVfit_phiTauUp(unsigned int i){return Ntp->SVfit_phiTauUp->at(i);} */
   /* float  SVfit_phiTauDown(unsigned int i){return Ntp->SVfit_phiTauDown->at(i);} */
   /* float  SVfit_etaUncTauUp(unsigned int i){return Ntp->SVfit_etaUncTauUp->at(i);} */
   /* float  SVfit_etaUncTauDown(unsigned int i){return Ntp->SVfit_etaUncTauDown->at(i);} */
   /* float  SVfit_etaTauUp(unsigned int i){return Ntp->SVfit_etaTauUp->at(i);} */
   /* float  SVfit_etaTauDown(unsigned int i){return Ntp->SVfit_etaTauDown->at(i);} */
   /* float  SVfit_ptUncTauUp(unsigned int i){return Ntp->SVfit_ptUncTauUp->at(i);} */
   /* float  SVfit_ptUncTauDown(unsigned int i){return Ntp->SVfit_ptUncTauDown->at(i);} */
   /* float  SVfit_ptTauUp(unsigned int i){return Ntp->SVfit_ptTauUp->at(i);} */
   /* float  SVfit_ptTauDown(unsigned int i){return Ntp->SVfit_ptTauDown->at(i);} */
   /* float  SVfitTransverseMassTauUp(unsigned int i){return Ntp->SVfitTransverseMassTauUp->at(i);} */
   /* float  SVfitTransverseMassTauDown(unsigned int i){return Ntp->SVfitTransverseMassTauDown->at(i);} */
   /* float  SVfitMassTauUp(unsigned int i){return Ntp->SVfitMassTauUp->at(i);} */
   /* float  SVfitMassTauDown(unsigned int i){return Ntp->SVfitMassTauDown->at(i);} */
   
   /* float  SVfitMass(unsigned int i){return Ntp->SVfitMass->at(i);} */
   /* float  SVfitTransverseMass(unsigned int i){return Ntp->SVfitTransverseMass->at(i);} */
   /* float  SVfit_pt(unsigned int i){return Ntp->SVfit_pt->at(i);} */
   /* float  SVfit_ptUnc(unsigned int i){return Ntp->SVfit_ptUnc->at(i);} */
   /* float  SVfit_eta(unsigned int i){return Ntp->SVfit_eta->at(i);} */
   /* float  SVfit_etaUnc(unsigned int i){return Ntp->SVfit_etaUnc->at(i);} */
   /* float  SVfit_phi(unsigned int i){return Ntp->SVfit_phi->at(i);} */
   /* float  SVfit_phiUnc(unsigned int i){return Ntp->SVfit_phiUnc->at(i);} */
   /* float  SVfit_fitMETRho(unsigned int i){return Ntp->SVfit_fitMETRho->at(i);} */
   /* float  SVfit_fitMETPhi(unsigned int i){return Ntp->SVfit_fitMETPhi->at(i);} */



   void RecoilCorr(TLorentzVector Gen,TLorentzVector Vis, int Index,float &PUPPImetCorr_px,float &PUPPImetCorr_py, string JER="Nom", string METScale="Nom", string METReso="Nom");
   float  METx(unsigned int i){return Ntp->METx->at(i);}   // index here is a pair
   float  METy(unsigned int i){return Ntp->METy->at(i);}
   //float  uncorrMETx(unsigned int i){return Ntp->uncorrMETx->at(i);}
   //float  uncorrMETy(unsigned int i){return Ntp->uncorrMETy->at(i);}
   /* float  MET_cov00(unsigned int i){return Ntp->MET_cov00->at(i);} */
   /* float  MET_cov01(unsigned int i){return Ntp->MET_cov01->at(i);} */
   /* float  MET_cov10(unsigned int i){return Ntp->MET_cov10->at(i);} */
   /* float  MET_cov11(unsigned int i){return Ntp->MET_cov11->at(i);} */
   /* float  MET_significance(unsigned int i){return Ntp->MET_significance->at(i);} */


   //int NPairCandidates(){return Ntp->isOSCand->size(); } // by construction the size  is always one, no idea why they fill it this way
   //bool isOSCand(unsigned int i){return Ntp->isOSCand->at(i);}
   /* float mT_Dau1(unsigned int i){return Ntp->mT_Dau1->at(i);} */
   /* float mT_Dau2(unsigned int i){return Ntp->mT_Dau2->at(i);} */
   int   indexDau1(unsigned int i){return Ntp->indexDau1->at(i);}
   int   indexDau2(unsigned int i){return Ntp->indexDau2->at(i);}
   int   getBestPairHTauTau (TString whatApply = "All", bool debug = false); // returns best pair formed by idx1, idx2, using HTauTau strategy - for studies
   int   getBestPairPtAndRawIsoOrd ( TString whatApply = "All", bool debug = false); // returns best pair formed by idx1, idx2, sorting them by pt in each pair, then by raw iso
   static bool pairSort (const tauPair_t& pA, const tauPair_t& pB);
   static bool pairSortRawIso (const tauPair_t& pA, const tauPair_t& pB);

   /* static bool pairSort (const tauPair_tuple& pA, const tauPair_tuple& pB); */
   /* static bool pairSortRawIso (const tauPair_tuple& pA, const tauPair_tuple& pB); */

   unsigned int    NDaughters() { return Ntp->daughters_px->size();}
   TLorentzVector Daughters_P4(unsigned int i){return TLorentzVector(Ntp->daughters_px->at(i), Ntp->daughters_py->at(i), Ntp->daughters_pz->at(i),Ntp->daughters_e->at(i));}
   int Daughters_P4_Size(){return Ntp->daughters_px->size();}
   int                   Daughters_charge(unsigned int i){return Ntp->daughters_charge->at(i);}
   //bool                 Daughters_TauUpExists(unsigned int i){return Ntp->daughters_TauUpExists->at(i);}
   //   TLorentzVector Daughters_px_TauUp_P4(unsigned int i){return TLorentzVector(Ntp->daughters_px_TauUp->at(i), Ntp->daughters_py_TauUp->at(i), Ntp->daughters_pz_TauUp->at(i),Ntp->daughters_e_TauUp->at(i));}
   //bool                 Daughters_TauDownExists(unsigned int i){return Ntp->daughters_TauDownExists->at(i);}
   //   TLorentzVector Daughters_px_TauDown_P4(unsigned int i){return TLorentzVector(Ntp->daughters_px_TauDown->at(i), Ntp->daughters_py_TauDown->at(i), Ntp->daughters_pz_TauDown->at(i),Ntp->daughters_e_TauDown->at(i));}
   int                   Daughters_genindex(unsigned int i){return Ntp->daughters_genindex->at(i);}


   unsigned int    NChargedDaughters(){return Ntp->daughters_charged_px->size();}
   TLorentzVector ChargedDaughters_P4(unsigned int i){return TLorentzVector(Ntp->daughters_charged_px->at(i), Ntp->daughters_charged_py->at(i), Ntp->daughters_charged_pz->at(i),Ntp->daughters_charged_e->at(i));}
   unsigned int    NNeutralDaughters(){return Ntp->daughters_neutral_px->size();}
   TLorentzVector NeutralDaughters_P4(unsigned int i){return TLorentzVector(Ntp->daughters_neutral_px->at(i), Ntp->daughters_neutral_py->at(i), Ntp->daughters_neutral_pz->at(i),Ntp->daughters_neutral_e->at(i));}


   //bool  Daughters_iseleBDT(unsigned int i){return Ntp->daughters_iseleBDT->at(i);}
bool  Daughters_iseleWP80(unsigned int i){return Ntp->daughters_iseleWP80->at(i);}
bool  Daughters_iseleWP90(unsigned int i){return Ntp->daughters_iseleWP90->at(i);}
bool  Daughters_iseleNoIsoWP90(unsigned int i){return Ntp->daughters_iseleNoIsoWP90->at(i);}
bool  Daughters_iseleNoIsoWPLoose(unsigned int i){return Ntp->daughters_iseleNoIsoWPLoose->at(i);}
float Daughters_eleMVAnt(unsigned int i){return Ntp->daughters_eleMVAnt->at(i);}
float Daughters_eleMVA_HZZ(unsigned int i){return Ntp->daughters_eleMVA_HZZ->at(i);}
bool  Daughters_passConversionVeto(unsigned int i){return Ntp->daughters_passConversionVeto->at(i);}
int   Daughters_eleMissingHits(unsigned int i){return Ntp->daughters_eleMissingHits->at(i);}
bool  Daughters_iseleChargeConsistent(unsigned int i){return Ntp->daughters_iseleChargeConsistent->at(i);}
//int   Daughters_eleCUTID(unsigned int i){return Ntp->daughters_eleCUTID->at(i);}
int   decayMode(unsigned int i){return Ntp->decayMode->at(i);}
/* 0: 1prong + 0 pi0 */
/* 1: 1prong + 1pi0 */
/* 2: 1prong + 2pi0s */
/* 5: 2prongs + 0pi0 */
/* 6: 2prongs + 1 pi0s */
/* 7: 2prongs + 2 pi0s */
/* 10: 3prongs + 0 pi0s */
/* 11: 3prongs + 1 pi0 */



 Long64_t  tauID(unsigned int i){return Ntp->tauID->at(i);}
 //int MVADM2016(unsigned int i){return Ntp->MVADM2016->at(i);}
 //int MVADM2016Size(){return Ntp->MVADM2016->size();}
 int MVADM2017(unsigned int i){return Ntp->MVADM2017->at(i);}
 int MVADM2017Size(){return Ntp->MVADM2017->size();}
 float  combreliso(unsigned int i){return Ntp->combreliso->at(i);}
 float  combreliso03(unsigned int i){return Ntp->combreliso03->at(i);}
 int    PDGIdDaughters(unsigned int i){return Ntp->PDGIdDaughters->at(i);}
 /* static const int ntauIds = 30;  */
 /* TString tauIDStrings[ntauIds] = { */
 /*     "byLoosePileupWeightedIsolation3Hits",  */
 /*     "byMediumPileupWeightedIsolation3Hits", */
 /*     "byTightPileupWeightedIsolation3Hits", */
 /*     "byLooseCombinedIsolationDeltaBetaCorr3Hits", */
 /*     "byMediumCombinedIsolationDeltaBetaCorr3Hits", */
 /*     "byTightCombinedIsolationDeltaBetaCorr3Hits",  */
 /*     "againstMuonLoose3",  */
 /*     "againstMuonTight3",  */
 /*     "againstElectronVLooseMVA6", */
 /*     "againstElectronLooseMVA6",  */
 /*     "againstElectronMediumMVA6", */
 /*     "againstElectronTightMVA6",  */
 /*     "againstElectronVTightMVA6", */
 /*     "byVLooseIsolationMVArun2v1DBoldDMwLT",  */
 /*     "byLooseIsolationMVArun2v1DBoldDMwLT",  */
 /*     "byMediumIsolationMVArun2v1DBoldDMwLT", */
 /*     "byTightIsolationMVArun2v1DBoldDMwLT",  */
 /*     "byVTightIsolationMVArun2v1DBoldDMwLT", */
 /*     "byVLooseIsolationMVArun2v1DBnewDMwLT", */
 /*     "byLooseIsolationMVArun2v1DBnewDMwLT",  */
 /*     "byMediumIsolationMVArun2v1DBnewDMwLT",  */
 /*     "byTightIsolationMVArun2v1DBnewDMwLT",  */
 /*     "byVTightIsolationMVArun2v1DBnewDMwLT",  */
 /*     "byLooseCombinedIsolationDeltaBetaCorr3HitsdR03", */
 /*     "byMediumCombinedIsolationDeltaBetaCorr3HitsdR03",  */
 /*     "byTightCombinedIsolationDeltaBetaCorr3HitsdR03", */
 /*     "byLooseIsolationMVArun2v1DBdR03oldDMwLT",  */
 /*     "byMediumIsolationMVArun2v1DBdR03oldDMwLT",  */
 /*     "byTightIsolationMVArun2v1DBdR03oldDMwLT",  */
 /*     "byVTightIsolationMVArun2v1DBdR03oldDMwLT"  */
 /*   };  */

float  Daughters_depositR03_tracker(unsigned int i){return Ntp->daughters_depositR03_tracker->at(i);}
float  Daughters_depositR03_ecal(unsigned int i){return Ntp->daughters_depositR03_ecal->at(i);}
float  Daughters_depositR03_hcal(unsigned int i){return Ntp->daughters_depositR03_hcal->at(i);}
int     Daughters_decayModeFindingOldDMs(unsigned int i){return Ntp->daughters_decayModeFindingOldDMs->at(i);}
int     Daughters_decayModeFindingNewDMs(unsigned int i){return Ntp->daughters_decayModeFindingNewDMs->at(i);}


//float  Daughters_SCeta(unsigned int i){return Ntp->daughters_SCeta->at(i);}
//float  againstElectronMVA5category(unsigned int i){return Ntp->againstElectronMVA5category->at(i);}
//float  againstElectronMVA5raw(unsigned int i){return Ntp->againstElectronMVA5raw->at(i);}
//float  byPileupWeightedIsolationRaw3Hits(unsigned int i){return Ntp->byPileupWeightedIsolationRaw3Hits->at(i);}
float  footprintCorrection(unsigned int i){return Ntp->footprintCorrection->at(i);}
float  neutralIsoPtSumWeight(unsigned int i){return Ntp->neutralIsoPtSumWeight->at(i);}
float  photonPtSumOutsideSignalCone(unsigned int i){return Ntp->photonPtSumOutsideSignalCone->at(i);}
//float  Daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits(unsigned int i){return Ntp->daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits->at(i);}
//float  Daughters_byIsolationMVA3oldDMwoLTraw(unsigned int i){return Ntp->daughters_byIsolationMVA3oldDMwoLTraw->at(i);}
//float  Daughters_byIsolationMVA3oldDMwLTraw(unsigned int i){return Ntp->daughters_byIsolationMVA3oldDMwLTraw->at(i);}
//float  Daughters_byIsolationMVA3newDMwoLTraw(unsigned int i){return Ntp->daughters_byIsolationMVA3newDMwoLTraw->at(i);}
//float  Daughters_byIsolationMVA3newDMwLTraw(unsigned int i){return Ntp->daughters_byIsolationMVA3newDMwLTraw->at(i);}
//float  Daughters_byIsolationMVArun2v1DBoldDMwLTraw(unsigned int i){return Ntp->daughters_byIsolationMVArun2v1DBoldDMwLTraw->at(i);}
float  Daughters_byDeepTau2017v2p1VSjetraw(unsigned int i){return Ntp->daughters_byDeepTau2017v2p1VSjetraw->at(i);}
float  Daughters_byDeepTau2017v2p1VSeraw(unsigned int i){return Ntp->daughters_byDeepTau2017v2p1VSeraw->at(i);}
float  Daughters_byDeepTau2017v2p1VSmuraw(unsigned int i){return Ntp->daughters_byDeepTau2017v2p1VSmuraw->at(i);}
float  Daughters_chargedIsoPtSum(unsigned int i){return Ntp->daughters_chargedIsoPtSum->at(i);}
float  Daughters_neutralIsoPtSum(unsigned int i){return Ntp->daughters_neutralIsoPtSum->at(i);}
float  Daughters_puCorrPtSum(unsigned int i){return Ntp->daughters_puCorrPtSum->at(i);}
int     Daughters_numChargedParticlesSignalCone(unsigned int i){return Ntp->daughters_numChargedParticlesSignalCone->at(i);}
int     Daughters_numNeutralHadronsSignalCone(unsigned int i){return Ntp->daughters_numNeutralHadronsSignalCone->at(i);}
int     Daughters_numPhotonsSignalCone(unsigned int i){return Ntp->daughters_numPhotonsSignalCone->at(i);}
int     Daughters_daughters_numParticlesSignalCone(unsigned int i){return Ntp->daughters_daughters_numParticlesSignalCone->at(i);}
int     Daughters_numChargedParticlesIsoCone(unsigned int i){return Ntp->daughters_numChargedParticlesIsoCone->at(i);}
int     Daughters_numNeutralHadronsIsoCone(unsigned int i){return Ntp->daughters_numNeutralHadronsIsoCone->at(i);}
int     Daughters_numPhotonsIsoCone(unsigned int i){return Ntp->daughters_numPhotonsIsoCone->at(i);}
int     Daughters_numParticlesIsoCone(unsigned int i){return Ntp->daughters_numParticlesIsoCone->at(i);}
float  Daughters_leadChargedParticlePt(unsigned int i){return Ntp->daughters_leadChargedParticlePt->at(i);}
float  Daughters_trackRefPt(unsigned int i){return Ntp->daughters_trackRefPt->at(i);}
//int     Daughters_isLastTriggerObjectforPath(unsigned int i){return Ntp->daughters_isLastTriggerObjectforPath->at(i);}
Long64_t Daughters_trgMatched(unsigned int i){return Ntp->daughters_trgMatched->at(i);}
Long64_t Daughters_FilterFired(unsigned int i){return Ntp->daughters_FilterFired->at(i);}
Long64_t Daughters_isGoodTriggerType(unsigned int i){return Ntp->daughters_isGoodTriggerType->at(i);}
Long64_t Daughters_L3FilterFired(unsigned int i){return Ntp->daughters_L3FilterFired->at(i);}
Long64_t Daughters_L3FilterFiredLast(unsigned int i){return Ntp->daughters_L3FilterFiredLast->at(i);}
float  Daughters_HLTpt(unsigned int i){return Ntp->daughters_HLTpt->at(i);}
//bool  Daughters_isL1IsoTau28Matched(unsigned int i){return Ntp->daughters_isL1IsoTau28Matched->at(i);}

int     Daughters_jetNDauChargedMVASel(unsigned int i){return Ntp->daughters_jetNDauChargedMVASel->at(i);}
float  Daughters_miniRelIsoCharged(unsigned int i){return Ntp->daughters_miniRelIsoCharged->at(i);}
float  Daughters_miniRelIsoNeutral(unsigned int i){return Ntp->daughters_miniRelIsoNeutral->at(i);}
float  Daughters_jetPtRel(unsigned int i){return Ntp->daughters_jetPtRel->at(i);}
float  Daughters_jetPtRatio(unsigned int i){return Ntp->daughters_jetPtRatio->at(i);}
float  Daughters_jetBTagCSV(unsigned int i){return Ntp->daughters_jetBTagCSV->at(i);}
//float  Daughters_lepMVA_mvaId(unsigned int i){return Ntp->daughters_lepMVA_mvaId->at(i);}


 TVector3 Daughters_pca(unsigned int i){return TVector3(Ntp->daughters_pca_x->at(i),Ntp->daughters_pca_y->at(i),Ntp->daughters_pca_z->at(i));}
 TVector3 Daughters_pcaRefitPV(unsigned int i){return TVector3(Ntp->daughters_pcaRefitPV_x->at(i),Ntp->daughters_pcaRefitPV_y->at(i),Ntp->daughters_pcaRefitPV_z->at(i));}
 TVector3 Daughters_pcaGenPV(unsigned int i){return TVector3(Ntp->daughters_pcaGenPV_x->at(i),Ntp->daughters_pcaGenPV_y->at(i),Ntp->daughters_pcaGenPV_z->at(i));}

 float  Daughters_IetaIeta(unsigned int i){return Ntp->daughters_IetaIeta->at(i);}
 float  Daughters_full5x5_IetaIeta(unsigned int i){return Ntp->daughters_full5x5_IetaIeta->at(i);}
 float  Daughters_hOverE(unsigned int i){return Ntp->daughters_hOverE->at(i);}
 float  Daughters_deltaEtaSuperClusterTrackAtVtx(unsigned int i){return Ntp->daughters_deltaEtaSuperClusterTrackAtVtx->at(i);}
 float  Daughters_deltaPhiSuperClusterTrackAtVtx(unsigned int i){return Ntp->daughters_deltaPhiSuperClusterTrackAtVtx->at(i);}
 float  Daughters_IoEmIoP(unsigned int i){return Ntp->daughters_IoEmIoP->at(i);}
 float  Daughters_IoEmIoP_ttH(unsigned int i){return Ntp->daughters_IoEmIoP_ttH->at(i);}

 float  dxy(unsigned int i){return Ntp->dxy->at(i);}
 float  dz(unsigned int i){return Ntp->dz->at(i);}
 float  dxy_innerTrack(unsigned int i){return Ntp->dxy_innerTrack->at(i);}
 //float  dz_innerTrack(unsigned int i){return Ntp->dz_innerTrack->at(i);}
 float  Daughters_rel_error_trackpt(unsigned int i){return Ntp->daughters_rel_error_trackpt->at(i);}
 float  SIP(unsigned int i){return Ntp->SIP->at(i);}
 int  Daughters_muonID(unsigned int i){return Ntp->daughters_muonID->at(i);} //bitwise (bit 0 loose, 1 soft , 2 medium, 3 tight, 4 highPT 5 tight_noVtx)
 int  Daughters_typeOfMuon(unsigned int i){return Ntp->daughters_typeOfMuon->at(i);}  // //bitwise, 0=PF, 1=Global, 2=Tracker
 int particleType(unsigned int i){return Ntp->particleType->at(i);} // 0 - muon, 1- electron, 2 - tau
 // float discriminator(unsigned int i){return Ntp->discriminator->at(i);}


 bool PFTau_hassecondaryVertex(unsigned int i){if(Ntp->PFTauSVPos->at(i).size()==3)return true; return false;}
 TVector3 PFTau_secondaryVertex_pos(unsigned int i){return  TVector3(Ntp->PFTauSVPos->at(i).at(0),Ntp->PFTauSVPos->at(i).at(1),Ntp->PFTauSVPos->at(i).at(2));}
 int PFTau_secondaryVertex_pos_Size(){return  Ntp->PFTauSVPos->size();}
 TMatrixTSym<double> PFTau_TIP_secondaryVertex_cov(unsigned int i);
 double PFTau_secondaryVertex_vtxchi2(unsigned int i){if(Ntp->PFTauSVChi2NDofMatchingQuality->at(i).size()==3) return  Ntp->PFTauSVChi2NDofMatchingQuality->at(i).at(0); return 0;}
 double PFTau_secondaryVertex_vtxndof(unsigned int i){if(Ntp->PFTauSVChi2NDofMatchingQuality->at(i).size()==3) return  Ntp->PFTauSVChi2NDofMatchingQuality->at(i).at(1);  return 0;}
 double PFTau_secondaryVertex_TracksMatchingQuality(unsigned int i){if(Ntp->PFTauSVChi2NDofMatchingQuality->at(i).size()==3) return  Ntp->PFTauSVChi2NDofMatchingQuality->at(i).at(2);  return 0;}



 bool PFtauHasPions(unsigned int i){if(Ntp->PFTauPionsP4->at(i).size()!=0){ return true;} return false;}
 bool PFtauHasThreePions(unsigned int i){if(Ntp->PFTauPionsP4->at(i).size()==3){ return true;} return false;}
 int NPions(unsigned int i){return Ntp->PFTauPionsP4->at(i).size();}
 TLorentzVector PFTau_PionsP4(unsigned int i, unsigned int j){return TLorentzVector(Ntp->PFTauPionsP4->at(i).at(j).at(1),Ntp->PFTauPionsP4->at(i).at(j).at(2),Ntp->PFTauPionsP4->at(i).at(j).at(3),Ntp->PFTauPionsP4->at(i).at(j).at(0) );}
 int PFTau_PionsP4_Size(){return Ntp->PFTauPionsP4->size();}
 int PFTau_PionsP4_SizePions(unsigned int i){return Ntp->PFTauPionsP4->at(i).size();}
 TLorentzVector PFTauRefit_PionsP4(unsigned int i, unsigned int j){return TLorentzVector(Ntp->PFTauRefitPionsP4->at(i).at(j).at(1),Ntp->PFTauRefitPionsP4->at(i).at(j).at(2),Ntp->PFTauRefitPionsP4->at(i).at(j).at(3),Ntp->PFTauRefitPionsP4->at(i).at(j).at(0) );}
 int PFTauRefit_PionsP4_Size(){return Ntp->PFTauRefitPionsP4->size();}
 int PFTauRefit_PionsP4_SizePions(unsigned int i){return Ntp->PFTauRefitPionsP4->at(i).size();}
 
 double PFTau_PionsCharge(unsigned int i, unsigned int j){return Ntp->PFTauPionsCharge->at(i).at(j);}
 double PFTauRefit_PionsCharge(unsigned int i, unsigned int j){return Ntp->PFTauRefitPionsCharge->at(i).at(j);}
 
 // MC Information
 // Signal particles (Z0,W+/-,H0,H+/-)
 unsigned int               NMCSignalParticles(){return Ntp->MCSignalParticle_p4->size();}
 TLorentzVector             MCSignalParticle_p4(unsigned int i){return TLorentzVector(Ntp->MCSignalParticle_p4->at(i).at(1),Ntp->MCSignalParticle_p4->at(i).at(2),Ntp->MCSignalParticle_p4->at(i).at(3),Ntp->MCSignalParticle_p4->at(i).at(0));}
 int                        MCSignalParticle_pdgid(unsigned int i){return Ntp->MCSignalParticle_pdgid->at(i);}
 int                        MCSignalParticle_charge(unsigned int i){return Ntp->MCSignalParticle_charge->at(i);}
 TVector3                   MCSignalParticle_Poca(unsigned int i){return TVector3(Ntp->MCSignalParticle_Poca->at(i).at(0),Ntp->MCSignalParticle_Poca->at(i).at(1),Ntp->MCSignalParticle_Poca->at(i).at(2));}
 std::vector<unsigned int>  MCSignalParticle_Tauidx(unsigned int i){return Ntp->MCSignalParticle_Tauidx->at(i);}
 // full MC chain
 unsigned int               NMCParticles(){return Ntp->MC_p4->size();}
 TLorentzVector             MCParticle_p4(unsigned int i){return TLorentzVector(Ntp->MC_p4->at(i).at(1),Ntp->MC_p4->at(i).at(2),Ntp->MC_p4->at(i).at(3),Ntp->MC_p4->at(i).at(0));}
 int                        MCParticle_pdgid(unsigned int i){return Ntp->MC_pdgid->at(i);}
 int                        MCParticle_charge(unsigned int i){return Ntp->MC_charge->at(i);}
 //int              		  MCParticle_midx(unsigned int i){return Ntp->MC_midx->at(i);}
 int              		  MCParticle_midx(unsigned int i){return Ntp->MC_midx->at(i);}
 std::vector<int>           MCParticle_childpdgid(unsigned int i){return Ntp->MC_childpdgid->at(i);}
 std::vector<int>           MCParticle_childidx(unsigned int i){return Ntp->MC_childidx->at(i);}
 int						  MCParticle_status(unsigned int i){return Ntp->MC_status->at(i);}
 int 						  getMatchTruthIndex(TLorentzVector tvector, int pid, double dr);
 int						  matchTruth(TLorentzVector tvector);
 bool						  matchTruth(TLorentzVector tvector, int pid, double dr);
 // decay tree functionality
 bool						  MCParticle_hasMother(unsigned int i){return Ntp->MC_midx->at(i) >= 0;}
 void						  printMCDecayChainOfMother(unsigned int i, bool printStatus = false, bool printPt = false, bool printEtaPhi = false, bool printQCD = false); // decay chain of object i
 void						  printMCDecayChainOfEvent(bool printStatus = false, bool printPt = false, bool printEtaPhi = false, bool printQCD = false); // full event decay chain
 std::string				  MCParticleToString(unsigned int par, bool printStatus = false, bool printPt = false, bool printEtaPhi = false);
 bool CheckDecayID( unsigned int jak1,unsigned  int jak2);
 TLorentzVector GetTruthTauLV(unsigned int jak, unsigned  int number);
 TLorentzVector GetTruthTauProductLV(unsigned int jak, int pdgID, unsigned  int number);
 std::vector<TLorentzVector> GetTruthPionsFromA1(unsigned int number=0);
 std::vector<TLorentzVector> GetTruthPionsFromRho(unsigned int number=0);
 


 // Tau decays (Tau is first element of vector)
 unsigned int NMCTaus(){return Ntp->MCTauandProd_p4->size();}
 TLorentzVector MCTau_p4(unsigned int i){return MCTauandProd_p4(i,0);}
 int MCTau_pdgid(unsigned int i){return MCTauandProd_pdgid(i,0);}
 int MCTau_charge(unsigned int i){return MCTauandProd_charge(i,0);}
 unsigned int MCTau_JAK(unsigned int i){return Ntp->MCTau_JAK->at(i);}
 unsigned int MCTau_DecayBitMask(unsigned int i){return Ntp->MCTau_DecayBitMask->at(i);}
 int MCTau_getDaughterOfType(unsigned int i_mcTau, int daughter_pdgid, bool ignoreCharge = true);
 int MCTau_true3prongAmbiguity(unsigned int i);
 int matchTauTruth(unsigned int i_hpsTau, bool onlyHadrDecays = false);
 TLorentzVector BoostToRestFrame(TLorentzVector TLV1, TLorentzVector TLV2);
 TLorentzVector MCTau_invisiblePart(unsigned int i);
 TLorentzVector MCTau_visiblePart(unsigned int i);
 
 //Tau and decay products
 int NMCTauDecayProducts(unsigned int i){if(0<=i && i<(unsigned int)NMCTaus()) return Ntp->MCTauandProd_p4->at(i).size(); return 0;}
 TLorentzVector MCTauandProd_p4(unsigned int i, unsigned int j){return TLorentzVector(Ntp->MCTauandProd_p4->at(i).at(j).at(1),Ntp->MCTauandProd_p4->at(i).at(j).at(2),Ntp->MCTauandProd_p4->at(i).at(j).at(3),Ntp->MCTauandProd_p4->at(i).at(j).at(0));}
 int MCTauandProd_pdgid(unsigned int i, unsigned int j){return Ntp->MCTauandProd_pdgid->at(i).at(j);}
 unsigned int MCTauandProd_midx(unsigned int i, unsigned int j){return Ntp->MCTauandProd_midx->at(i).at(j);}
 int MCTauandProd_charge(unsigned int i, unsigned int j){return Ntp->MCTauandProd_charge->at(i).at(j);}
 TVector3 MCTauandProd_Vertex(unsigned int i, unsigned int j){
   return TVector3(Ntp->MCTauandProd_Vertex->at(i).at(j).at(0),Ntp->MCTauandProd_Vertex->at(i).at(j).at(1),Ntp->MCTauandProd_Vertex->at(i).at(j).at(2));
 }
 int MCTauandProd_VertexSize(){return Ntp->MCTauandProd_Vertex->size();}
 bool hasSignalTauDecay(PDGInfo::PDGMCNumbering parent_pdgid,unsigned int &Boson_idx,TauDecay::JAK tau_jak, unsigned int &idx);
 bool hasSignalTauDecay(PDGInfo::PDGMCNumbering parent_pdgid,unsigned int &Boson_idx,unsigned int &tau1_idx, unsigned int &tau2_idx);


 
 //Int_t           JetsNumber(){return Ntp->JetsNumber;}
 bool tightJetID(unsigned int i){return Ntp->tightJetID->at(i);}
 unsigned int NJets(){return Ntp->jets_px->size();}
 unsigned int NJetsDown(){return Ntp->jetsDown_px->size();}
 unsigned int NJetsUp(){return Ntp->jetsUp_px->size();}
 
 TLorentzVector Jet_P4(unsigned int i){return TLorentzVector(Ntp->jets_px->at(i), Ntp->jets_py->at(i), Ntp->jets_pz->at(i),Ntp->jets_e->at(i));}
 TLorentzVector JetDown_P4(unsigned int i){return TLorentzVector(Ntp->jetsDown_px->at(i), Ntp->jetsDown_py->at(i), Ntp->jetsDown_pz->at(i),Ntp->jetsDown_e->at(i));}
 TLorentzVector JetUp_P4(unsigned int i){return TLorentzVector(Ntp->jetsUp_px->at(i), Ntp->jetsUp_py->at(i), Ntp->jetsUp_pz->at(i),Ntp->jetsUp_e->at(i));}
 
 // float jets_rawPt(unsigned int i){return Ntp->jets_rawPt->at(i);}
 float jets_area(unsigned int i){return Ntp->jets_area->at(i);}
 float jets_mT(unsigned int i){return Ntp->jets_mT->at(i);}
 int jets_Flavour(unsigned int i){return Ntp->jets_Flavour->at(i);}
 int jets_HadronFlavour(unsigned int i){return Ntp->jets_HadronFlavour->at(i);}
 int jets_genjetIndex(unsigned int i){return Ntp->jets_genjetIndex->at(i);}
 float jets_PUJetID(unsigned int i){return Ntp->jets_PUJetID->at(i);}
 float jets_PUJetIDupdated(unsigned int i){return Ntp->jets_PUJetIDupdated->at(i);}
 float jets_vtxPt(unsigned int i){return Ntp->jets_vtxPt->at(i);}
 float jets_vtxMass(unsigned int i){return Ntp->jets_vtxMass->at(i);}
 float jets_vtx3dL(unsigned int i){return Ntp->jets_vtx3dL->at(i);}
 float jets_vtxNtrk(unsigned int i){return Ntp->jets_vtxNtrk->at(i);}
 float jets_vtx3deL(unsigned int i){return Ntp->jets_vtx3deL->at(i);}
 float jets_leadTrackPt(unsigned int i){return Ntp->jets_leadTrackPt->at(i);}
 float jets_leptonPtRel(unsigned int i){return Ntp->jets_leptonPtRel->at(i);}
 float jets_leptonPt(unsigned int i){return Ntp->jets_leptonPt->at(i);}
 float jets_leptonDeltaR(unsigned int i){return Ntp->jets_leptonDeltaR->at(i);}
 float jets_chEmEF(unsigned int i){return Ntp->_jets_chEmEF->at(i);}
 float jets_chHEF(unsigned int i){return Ntp->_jets_chHEF->at(i);}
 float jets_nEmEF(unsigned int i){return Ntp->_jets_nEmEF->at(i);}
 float jets_nHEF(unsigned int i){return Ntp->_jets_nHEF->at(i);}
 float jets_MUF(unsigned int i){return Ntp->jets_MUF->at(i);}
 int   jets_neMult(unsigned int i){return Ntp->_jets_neMult->at(i);}
 int   jets_chMult(unsigned int i){return Ntp->_jets_chMult->at(i);}
 float jets_jecUnc(unsigned int i){return Ntp->jets_jecUnc->at(i);}
 int jets_jetUnc_AbsoluteSize_up(){return Ntp->jets_jetUnc_Absolute_up->size();}
 float jets_jetUnc_Absolute_up(int i){return Ntp->jets_jetUnc_Absolute_up->at(i);}
 float jets_jetUnc_FlavorQCD_up(int i){return Ntp->jets_jetUnc_FlavorQCD_up->at(i);}
 float jets_jetUnc_RelativeBal_up(int i){return Ntp->jets_jetUnc_RelativeBal_up->at(i);}
 float jets_jetUnc_HF_up(int i){return Ntp->jets_jetUnc_HF_up->at(i);}
 float jets_jetUnc_BBEC1_up(int i){return Ntp->jets_jetUnc_BBEC1_up->at(i);}
 float jets_jetUnc_EC2_up(int i){return Ntp->jets_jetUnc_EC2_up->at(i);}
 float jets_jetUnc_BBEC1_YEAR_up(int i){return Ntp->jets_jetUnc_BBEC1_YEAR_up->at(i);}
 float jets_jetUnc_EC2_YEAR_up(int i){return Ntp->jets_jetUnc_EC2_YEAR_up->at(i);}
 float jets_jetUnc_Absolute_YEAR_up(int i){return Ntp->jets_jetUnc_Absolute_YEAR_up->at(i);}
 float jets_jetUnc_HF_YEAR_up(int i){return Ntp->jets_jetUnc_HF_YEAR_up->at(i);}
 float jets_jetUnc_RelativeSample_YEAR_up(int i){return Ntp->jets_jetUnc_RelativeSample_YEAR_up->at(i);}
 float jets_jetUnc_Absolute_dw(int i){return Ntp->jets_jetUnc_Absolute_dw->at(i);}
 float jets_jetUnc_FlavorQCD_dw(int i){return Ntp->jets_jetUnc_FlavorQCD_dw->at(i);}
 float jets_jetUnc_RelativeBal_dw(int i){return Ntp->jets_jetUnc_RelativeBal_dw->at(i);}
 float jets_jetUnc_HF_dw(int i){return Ntp->jets_jetUnc_HF_dw->at(i);}
 float jets_jetUnc_BBEC1_dw(int i){return Ntp->jets_jetUnc_BBEC1_dw->at(i);}
 float jets_jetUnc_EC2_dw(int i){return Ntp->jets_jetUnc_EC2_dw->at(i);}
 float jets_jetUnc_BBEC1_YEAR_dw(int i){return Ntp->jets_jetUnc_BBEC1_YEAR_dw->at(i);}
 float jets_jetUnc_EC2_YEAR_dw(int i){return Ntp->jets_jetUnc_EC2_YEAR_dw->at(i);}
 float jets_jetUnc_Absolute_YEAR_dw(int i){return Ntp->jets_jetUnc_Absolute_YEAR_dw->at(i);}
 float jets_jetUnc_HF_YEAR_dw(int i){return Ntp->jets_jetUnc_HF_YEAR_dw->at(i);}
 float jets_jetUnc_RelativeSample_YEAR_dw(int i){return Ntp->jets_jetUnc_RelativeSample_YEAR_dw->at(i);}
 
 float  bDiscriminator(unsigned int i){return Ntp->bDiscriminator->at(i);}
 float  bCSVscore(unsigned int i){return Ntp->bCSVscore->at(i);}
 float  DeepCSV_probb(unsigned int i){return Ntp->bDeepCSV_probb->at(i);}
 float  DeepCSV_probbb(unsigned int i){return Ntp->bDeepCSV_probbb->at(i);}
 
 float  pfCombinedMVAV2BJetTags(unsigned int i){return Ntp->pfCombinedMVAV2BJetTags->at(i);}
 float  jetRawf(unsigned int i){return Ntp->jetRawf->at(i);}
 // int  PFjetID(unsigned int i){return Ntp->PFjetID->at(i);}
 float ttbarPtWeight();

 unsigned int NAK8Jets(){return Ntp->ak8jets_px->size();}
 TLorentzVector AK8Jet_P4(unsigned int i){return TLorentzVector(Ntp->ak8jets_px->at(i), Ntp->ak8jets_py->at(i), Ntp->ak8jets_pz->at(i),Ntp->ak8jets_e->at(i));}
 float ak8jets_SoftDropMass(unsigned int i){return Ntp->ak8jets_SoftDropMass->at(i);}
 /* float ak8jets_PrunedMass(unsigned int i){return Ntp->ak8jets_PrunedMass->at(i);} */
 /* float ak8jets_TrimmedMass(unsigned int i){return Ntp->ak8jets_TrimmedMass->at(i);} */
 /* float ak8jets_FilteredMass(unsigned int i){return Ntp->ak8jets_FilteredMass->at(i);} */
 float ak8jets_tau1(unsigned int i){return Ntp->ak8jets_tau1->at(i);}
 float ak8jets_tau2(unsigned int i){return Ntp->ak8jets_tau2->at(i);}
 float ak8jets_tau3(unsigned int i){return Ntp->ak8jets_tau3->at(i);}
 float ak8jets_CSV(unsigned int i){return Ntp->ak8jets_CSV->at(i);}
 int ak8jets_nsubjets(unsigned int i){return Ntp->ak8jets_nsubjets->at(i);}
 
 unsigned int NSubJets(){return Ntp->subjets_px->size();}
 TLorentzVector SUBJet_P4(unsigned int i){return TLorentzVector(Ntp->subjets_px->at(i), Ntp->subjets_py->at(i), Ntp->subjets_pz->at(i),Ntp->subjets_e->at(i));}
 
 float subjets_CSV(unsigned int i){return Ntp->subjets_CSV->at(i);}
 int subjets_ak8MotherIdx(unsigned int i){return Ntp->subjets_ak8MotherIdx->at(i);}

 bool CHECK_BIT(unsigned long long var, int pos){  
   unsigned long long CHECK1=(unsigned long long)(var & ((unsigned long long)(1) << pos));
   unsigned long long CHECK2=(unsigned long long)((unsigned long long)(1) << pos);
   return ( CHECK1==CHECK2 );
 }
 //
 double stitch_weight(bool isDY1050);
 
 Int_t year(){return Ntp->year;}
 Int_t SelectedPairs(){return Ntp->SelectedPairs;}
 Bool_t eleveto(){return Ntp->eleveto;}
 Bool_t muonveto(){return Ntp->muonveto;}
 Bool_t trg_doubletau(int k){return Ntp->trg_doubletau->at(k);}
  
 Int_t genmatch (int i){return Ntp->genmatch->at(i);}

 Float_t puppimt_1(int k){return Ntp->puppimt_1->at(k);}
 Float_t gen_match_1(int k){return Ntp->gen_match_1->at(k);}
  
 Float_t trigweight_1(int k){return Ntp->trigweight_1->at(k);}
 Float_t idisoweight_1(int k){return Ntp->idisoweight_1->at(k);}
 Float_t antieweight_1(int k){return Ntp->antieweight_1->at(k);}
 Float_t antimuweight_1(int k){return Ntp->antimuweight_1->at(k);}

 Float_t puppimt_2(int k){return Ntp->puppimt_2->at(k);}
 Float_t gen_match_2(int k){return Ntp->gen_match_2->at(k);}

 Float_t trigweight_2(int k){return Ntp->trigweight_2->at(k);}
 Float_t idisoweight_2(int k){return Ntp->idisoweight_2->at(k);}
 Float_t antieweight_2(int k){return Ntp->antieweight_2->at(k);}
 Float_t antimuweight_2(int k){return Ntp->antimuweight_2->at(k);}

 double IDSF(int i, int genmatch, string TES="Nom", string particle="tau", string Unc="Nom");
 double TriggerSF(int i,int genmatch,string TES="Nom", string Unc="Nom");
 double ZPtReweight(TLorentzVector GenP4);
 double EmbeddingSelectionSF(int imc1, int imc2);
 void FillHist(unsigned int t, int idx, bool isOS, bool GenMatchSelection, double value, std::pair<float, int> max_pair, double w, double wData, double wMC, std::vector<TH2D> *histHiggs=nullptr, std::vector<TH2D> *histJetFakes=nullptr, std::vector<TH2D> *histZTT=nullptr, std::vector<TH2D> *histWfakesHiggs=nullptr, std::vector<TH2D> *histWfakesJetFakes=nullptr, std::vector<TH2D> *histWfakesZTT=nullptr, std::vector<TH2D> *histHiggsQCDMC=nullptr, std::vector<TH2D> *histJetFakesQCDMC=nullptr, std::vector<TH2D> *histZTTQCDMC=nullptr);

 //Float_t pt_tt(int k){return Ntp->pt_tt->at(k);} 
 /* Float_t pt_vis(int k){return Ntp->pt_vis->at(k);} */
 Float_t mt_tot(int k){return Ntp->mt_tot->at(k);}
 // Float_t m_vis(int k){return Ntp->m_vis->at(k);}

 Float_t MET(){return Ntp->met;}
 Float_t METphi(){return Ntp->metphi;}
 Float_t PUPPImet(){return Ntp->PUPPImet;}
 Float_t puppimet_ex_UnclusteredEnUp(){return Ntp->puppimet_ex_UnclusteredEnUp;}
 Float_t puppimet_ey_UnclusteredEnUp(){return Ntp->puppimet_ey_UnclusteredEnUp;}
 Float_t puppimet_ex_UnclusteredEnDown(){return Ntp->puppimet_ex_UnclusteredEnDown;}
 Float_t puppimet_ey_UnclusteredEnDown(){return Ntp->puppimet_ey_UnclusteredEnDown;}
 Float_t PUPPImetphi(){return Ntp->PUPPImetphi;}
 Float_t PFMETCov00(){return Ntp->PFMETCov00;}
 Float_t PFMETCov01(){return Ntp->PFMETCov01;}
 Float_t PFMETCov10(){return Ntp->PFMETCov10;}
 Float_t PFMETCov11(){return Ntp->PFMETCov11;}
 Float_t PUPPIMETCov00(){return Ntp->PUPPIMETCov00;}
 Float_t PUPPIMETCov01(){return Ntp->PUPPIMETCov01;}
 Float_t PUPPIMETCov10(){return Ntp->PUPPIMETCov10;}
 Float_t PUPPIMETCov11(){return Ntp->PUPPIMETCov11;}

 Float_t mjj(int k){return Ntp->mjj->at(k);}
 Float_t jdeta(int k){return Ntp->jdeta->at(k);}
 Float_t mjjDown(int k){return Ntp->mjjDown->at(k);}
 Float_t jdetaDown(int k){return Ntp->jdetaDown->at(k);}
 Float_t mjjUp(int k){return Ntp->mjjUp->at(k);}
 Float_t jdetaUp(int k){return Ntp->jdetaUp->at(k);}
 
 Int_t njetingap(int k){return Ntp->njetingap->at(k);}
 Int_t njetingap20(int k){return Ntp->njetingap20->at(k);}
 Float_t jdphi(int k){return Ntp->jdphi->at(k);}
 Float_t dijetpt(int k){return Ntp->dijetpt->at(k);}
 Float_t dijetphi(int k){return Ntp->dijetphi->at(k);}
 //Float_t ptvis(int k){return Ntp->ptvis->at(k);}
 Int_t nbtag(int k){return Ntp->nbtag->at(k);}
 Int_t njetsSize(){return Ntp->njets->size();}
 Int_t njets(int k){return Ntp->njets->at(k);}
 Int_t njetspt20Size(){return Ntp->njetspt20->size();}
 Int_t njetspt20(int k){return Ntp->njetspt20->at(k);}
 Int_t njetsSizeUp(){return Ntp->njetsUp->size();}
 Int_t njetsUp(int k){return Ntp->njetsUp->at(k);}
 Int_t njetspt20Up(int k){return Ntp->njetspt20Up->at(k);}
 Int_t njetspt20SizeUp(){return Ntp->njetspt20Up->size();}
 Int_t njetsSizeDown(){return Ntp->njetsDown->size();}
 Int_t njetsDown(int k){return Ntp->njetsDown->at(k);}
 Int_t njetspt20Down(int k){return Ntp->njetspt20Down->at(k);}
 Int_t njetspt20SizeDown(){return Ntp->njetspt20Down->size();}
 int jptSize_1(){return Ntp->jpt_1->size();}
 Float_t jpt_1(int k){return Ntp->jpt_1->at(k);}
 Float_t jptDown_1(int k){return Ntp->jptDown_1->at(k);}
 Float_t jptUp_1(int k){return Ntp->jptUp_1->at(k);}
 
 Float_t jeta_1(int k){return Ntp->jeta_1->at(k);}
 Float_t jphi_1(int k){return Ntp->jphi_1->at(k);}
 Float_t jcsv_1(int k){return Ntp->jcsv_1->at(k);}
 Float_t jpt_2(int k){return Ntp->jpt_2->at(k);}
 Float_t jeta_2(int k){return Ntp->jeta_2->at(k);}
 Float_t jphi_2(int k){return Ntp->jphi_2->at(k);}
 Float_t jcsv_2(int k){return Ntp->jcsv_2->at(k);}
 Float_t bpt_1(int k){return Ntp->bpt_1->at(k);}
 Float_t beta_1(int k){return Ntp->beta_1->at(k);}
 Float_t bphi_1(int k){return Ntp->bphi_1->at(k);}
 Float_t bcsv_1(int k){return Ntp->bcsv_1->at(k);}
 Float_t bpt_2(int k){return Ntp->bpt_2->at(k);}
 Float_t beta_2(int k){return Ntp->beta_2->at(k);}
 Float_t bphi_2(int k){return Ntp->bphi_2->at(k);}
 Float_t bcsv_2(int k){return Ntp->bcsv_2->at(k);}

 Float_t puweight(){return Ntp->puweight;}

 Float_t weight(int k){return Ntp->weight->at(k);}

 /* Float_t jpfid_1(int k){return Ntp->jpfid_1->at(k);} */
 /* Float_t jpuid_1(int k){return Ntp->jpuid_1->at(k);} */
 /* Float_t jpfid_2(int k){return Ntp->jpfid_2->at(k);} */
 /* Float_t jpuid_2(int k){return Ntp->jpuid_2->at(k);} */
 /* Float_t bpfid_1(int k){return Ntp->bpfid_1->at(k);} */
 /* Float_t bpuid_1(int k){return Ntp->bpuid_1->at(k);} */
 /* Float_t bpfid_2(int k){return Ntp->bpfid_2->at(k);} */
 /* Float_t bpuid_2(int k){return Ntp->bpuid_2->at(k);} */

 // Float_t byDeepTau2017v2p1VSjetraw_1(int k){return Ntp->byDeepTau2017v2p1VSjetraw_1->at(k);}
 Float_t byVVVLooseDeepTau2017v2p1VSjet_1(int k){return Ntp->byVVVLooseDeepTau2017v2p1VSjet_1->at(k);}
 Float_t byVVLooseDeepTau2017v2p1VSjet_1(int k){return Ntp->byVVLooseDeepTau2017v2p1VSjet_1->at(k);} 
 Float_t byVLooseDeepTau2017v2p1VSjet_1(int k){return Ntp->byVLooseDeepTau2017v2p1VSjet_1->at(k);}  
 Float_t byLooseDeepTau2017v2p1VSjet_1(int k){return Ntp->byLooseDeepTau2017v2p1VSjet_1->at(k);}   
 Float_t byMediumDeepTau2017v2p1VSjet_1(int k){return Ntp->byMediumDeepTau2017v2p1VSjet_1->at(k);}  
 Float_t byTightDeepTau2017v2p1VSjet_1(int k){return Ntp->byTightDeepTau2017v2p1VSjet_1->at(k);}   
 Float_t byVTightDeepTau2017v2p1VSjet_1(int k){return Ntp->byVTightDeepTau2017v2p1VSjet_1->at(k);}  
 Float_t byVVTightDeepTau2017v2p1VSjet_1(int k){return Ntp->byVVTightDeepTau2017v2p1VSjet_1->at(k);} 
  
 // Float_t byDeepTau2017v2p1VSeleraw_1(int k){return Ntp->byDeepTau2017v2p1VSeleraw_1->at(k);}
 Float_t byVVVLooseDeepTau2017v2p1VSe_1(int k){return Ntp->byVVVLooseDeepTau2017v2p1VSe_1->at(k);}  
 Float_t byVVLooseDeepTau2017v2p1VSe_1(int k){return Ntp->byVVLooseDeepTau2017v2p1VSe_1->at(k);} 
 Float_t byVLooseDeepTau2017v2p1VSe_1(int k){return Ntp->byVLooseDeepTau2017v2p1VSe_1->at(k);}   
 Float_t byLooseDeepTau2017v2p1VSe_1(int k){return Ntp->byLooseDeepTau2017v2p1VSe_1->at(k);}	
 Float_t byMediumDeepTau2017v2p1VSe_1(int k){return Ntp->byMediumDeepTau2017v2p1VSe_1->at(k);}   
 Float_t byTightDeepTau2017v2p1VSe_1(int k){return Ntp->byTightDeepTau2017v2p1VSe_1->at(k);}	
 Float_t byVTightDeepTau2017v2p1VSe_1(int k){return Ntp->byVTightDeepTau2017v2p1VSe_1->at(k);}   
 Float_t byVVTightDeepTau2017v2p1VSe_1(int k){return Ntp->byVVTightDeepTau2017v2p1VSe_1->at(k);}
  
 // Float_t byDeepTau2017v2p1VSmuraw_1(int k){return Ntp->byDeepTau2017v2p1VSmuraw_1->at(k);}
 Float_t byVLooseDeepTau2017v2p1VSmu_1(int k){return Ntp->byVLooseDeepTau2017v2p1VSmu_1->at(k);} 
 Float_t byLooseDeepTau2017v2p1VSmu_1(int k){return Ntp->byLooseDeepTau2017v2p1VSmu_1->at(k);} 
 Float_t byMediumDeepTau2017v2p1VSmu_1(int k){return Ntp->byMediumDeepTau2017v2p1VSmu_1->at(k);}   
 Float_t byTightDeepTau2017v2p1VSmu_1(int k){return Ntp->byTightDeepTau2017v2p1VSmu_1->at(k);}

 // Float_t byDeepTau2017v2p1VSjetraw_2(int k){return Ntp->byDeepTau2017v2p1VSjetraw_2->at(k);}
 Float_t byVVVLooseDeepTau2017v2p1VSjet_2(int k){return Ntp->byVVVLooseDeepTau2017v2p1VSjet_2->at(k);}
 Float_t byVVLooseDeepTau2017v2p1VSjet_2(int k){return Ntp->byVVLooseDeepTau2017v2p1VSjet_2->at(k);} 
 Float_t byVLooseDeepTau2017v2p1VSjet_2(int k){return Ntp->byVLooseDeepTau2017v2p1VSjet_2->at(k);}  
 Float_t byLooseDeepTau2017v2p1VSjet_2(int k){return Ntp->byLooseDeepTau2017v2p1VSjet_2->at(k);}   
 Float_t byMediumDeepTau2017v2p1VSjet_2(int k){return Ntp->byMediumDeepTau2017v2p1VSjet_2->at(k);}  
 Float_t byTightDeepTau2017v2p1VSjet_2(int k){return Ntp->byTightDeepTau2017v2p1VSjet_2->at(k);}   
 Float_t byVTightDeepTau2017v2p1VSjet_2(int k){return Ntp->byVTightDeepTau2017v2p1VSjet_2->at(k);}  
 Float_t byVVTightDeepTau2017v2p1VSjet_2(int k){return Ntp->byVVTightDeepTau2017v2p1VSjet_2->at(k);} 
  
 // Float_t byDeepTau2017v2p1VSeleraw_2(int k){return Ntp->byDeepTau2017v2p1VSeleraw_2->at(k);}
 Float_t byVVVLooseDeepTau2017v2p1VSe_2(int k){return Ntp->byVVVLooseDeepTau2017v2p1VSe_2->at(k);}  
 Float_t byVVLooseDeepTau2017v2p1VSe_2(int k){return Ntp->byVVLooseDeepTau2017v2p1VSe_2->at(k);} 
 Float_t byVLooseDeepTau2017v2p1VSe_2(int k){return Ntp->byVLooseDeepTau2017v2p1VSe_2->at(k);}   
 Float_t byLooseDeepTau2017v2p1VSe_2(int k){return Ntp->byLooseDeepTau2017v2p1VSe_2->at(k);}	
 Float_t byMediumDeepTau2017v2p1VSe_2(int k){return Ntp->byMediumDeepTau2017v2p1VSe_2->at(k);}   
 Float_t byTightDeepTau2017v2p1VSe_2(int k){return Ntp->byTightDeepTau2017v2p1VSe_2->at(k);}	
 Float_t byVTightDeepTau2017v2p1VSe_2(int k){return Ntp->byVTightDeepTau2017v2p1VSe_2->at(k);}   
 Float_t byVVTightDeepTau2017v2p1VSe_2(int k){return Ntp->byVVTightDeepTau2017v2p1VSe_2->at(k);}
  
 // Float_t byDeepTau2017v2p1VSmuraw_2(int k){return Ntp->byDeepTau2017v2p1VSmuraw_2->at(k);}
 Float_t byVLooseDeepTau2017v2p1VSmu_2(int k){return Ntp->byVLooseDeepTau2017v2p1VSmu_2->at(k);} 
 Float_t byLooseDeepTau2017v2p1VSmu_2(int k){return Ntp->byLooseDeepTau2017v2p1VSmu_2->at(k);} 
 Float_t byMediumDeepTau2017v2p1VSmu_2(int k){return Ntp->byMediumDeepTau2017v2p1VSmu_2->at(k);}   
 Float_t byTightDeepTau2017v2p1VSmu_2(int k){return Ntp->byTightDeepTau2017v2p1VSmu_2->at(k);}

 /* Float_t pvx(int k){return Ntp->pvx->at(k);} */
 /* Float_t pvy(int k){return Ntp->pvy->at(k);} */
 /* Float_t pvz(int k){return Ntp->pvz->at(k);} */

 Float_t dm_1(int k){return Ntp->dm_1->at(k);}
 Float_t dmMVA_1(int k){return Ntp->dmMVA_1->at(k);}

 Float_t dm_2(int k){return Ntp->dm_2->at(k);}
 Float_t dmMVA_2(int k){return Ntp->dmMVA_2->at(k);}

 /* Float_t svx_1(int k){return Ntp->svx_1->at(k);} */
 /* Float_t svy_1(int k){return Ntp->svy_1->at(k);} */
 /* Float_t svz_1(int k){return Ntp->svz_1->at(k);} */

 /* Float_t svx_2(int k){return Ntp->svx_2->at(k);} */
 /* Float_t svy_2(int k){return Ntp->svy_2->at(k);} */
 /* Float_t svz_2(int k){return Ntp->svz_2->at(k);} */

 Int_t tau1IndexVect(int k){return Ntp->tau1IndexVect->at(k);}
 Int_t tau2IndexVect(int k){return Ntp->tau2IndexVect->at(k);}


 // bool res = word & (1 << bitpos);


 // CHECK_BIT(Ntp->triggerbit(), getBitOfGivenTrigger("HLT_IsoMu24"))






























  /* // Muon information */
  /* unsigned int   NMuons(){return Ntp->Muon_p4->size();} */
  /* TLorentzVector Muon_p4(unsigned int i, TString corr = "default"); */
  /* TVector3       Muon_Poca(unsigned int i){return TVector3(Ntp->Muon_Poca->at(i).at(0),Ntp->Muon_Poca->at(i).at(1),Ntp->Muon_Poca->at(i).at(2));} */
  // bool           Muon_isGlobalMuon(unsigned int i){return Ntp->Muon_isGlobalMuon->at(i);}
  /* bool           Muon_isStandAloneMuon(unsigned int i){return Ntp->Muon_isStandAloneMuon->at(i);} */
  // bool           Muon_isTrackerMuon(unsigned int i){return Ntp->Muon_isTrackerMuon->at(i);}
  /* bool           Muon_isCaloMuon(unsigned int i){return Ntp->Muon_isCaloMuon->at(i);} */
  /* bool           Muon_isIsolationValid(unsigned int i){return Ntp->Muon_isIsolationValid->at(i);} */
  /* bool           Muon_isQualityValid(unsigned int i){return Ntp->Muon_isQualityValid->at(i);} */
  /* bool           Muon_isTimeValid(unsigned int i){return Ntp->Muon_isTimeValid->at(i);} */
  /* float          Muon_emEt03(unsigned int i){return Ntp->Muon_emEt03->at(i);} */
  /* float          Muon_emVetoEt03(unsigned int i){return Ntp->Muon_emVetoEt03->at(i);} */
  /* float          Muon_hadEt03(unsigned int i){return Ntp->Muon_hadEt03->at(i);} */
  /* float          Muon_hadVetoEt03(unsigned int i){return Ntp->Muon_hadVetoEt03->at(i);} */
  /* int 	         Muon_nJets03(unsigned int i){return Ntp->Muon_nJets03->at(i);} */
  /* int            Muon_nTracks03(unsigned int i){return Ntp->Muon_nTracks03->at(i);} */
  /* float          Muon_sumPt03(unsigned int i){return Ntp->Muon_sumPt03->at(i);} */
  /* float          Muon_trackerVetoPt03(unsigned int i){return Ntp->Muon_trackerVetoPt03->at(i);} */
  /* float          Muon_emEt05(unsigned int i){return Ntp->Muon_emEt05->at(i);} */
  /* float          Muon_emVetoEt05(unsigned int i){return Ntp->Muon_emVetoEt05->at(i);} */
  /* float          Muon_hadEt05(unsigned int i){return Ntp->Muon_hadEt05->at(i);} */
  /* float          Muon_hadVetoEt05(unsigned int i){return Ntp->Muon_hadVetoEt05->at(i);} */
  /* int            Muon_nJets05(unsigned int i){return Ntp->Muon_nJets05->at(i);} */
  /* int            Muon_nTracks05(unsigned int i){return Ntp->Muon_nTracks05->at(i);} */
  /* float          Muon_sumPt05(unsigned int i){return Ntp->Muon_sumPt05->at(i);} */
  /* float          Muon_trackerVetoPt05(unsigned int i){return Ntp->Muon_trackerVetoPt05->at(i);} */
  /* unsigned int   Muon_Track_idx(unsigned int i){return Ntp->Muon_Track_idx->at(i);} */
  /* int            Muon_hitPattern_pixelLayerwithMeas(unsigned int i){return Ntp->Muon_hitPattern_pixelLayerwithMeas->at(i);} */
  /* int            Muon_numberOfMatchedStations(unsigned int i){return Ntp->Muon_numberOfMatchedStations->at(i);} */
  /* float          Muon_normChi2(unsigned int i){return Ntp->Muon_normChi2->at(i);} */
  /* int            Muon_hitPattern_numberOfValidMuonHits(unsigned int i){return Ntp->Muon_hitPattern_numberOfValidMuonHits->at(i);} */
  /* int            Muon_innerTrack_numberofValidHits(unsigned int i){return Ntp->Muon_innerTrack_numberofValidHits->at(i);} */
  /* int            Muon_numberOfMatches(unsigned int i){return Ntp->Muon_numberOfMatches->at(i);} */
  /* int            Muon_numberOfChambers(unsigned int i){return Ntp->Muon_numberOfChambers->at(i);} */
  /* int            Muon_Charge(unsigned int i){return Ntp->Muon_charge->at(i);} */
  /* int            Muon_trackCharge(unsigned int i){return Ntp->Muon_trackCharge->at(i);} */

  // bool           Muon_isPFMuon(unsigned int i){return Ntp->Muon_isPFMuon->at(i);}
  /* float          Muon_sumChargedHadronPt03(unsigned int i){return Ntp->Muon_sumChargedHadronPt03->at(i);}                             // sum-pt of charged Hadron					                                */
  /* float          Muon_sumChargedParticlePt03(unsigned int i){return Ntp->Muon_sumChargedParticlePt03->at(i);}			      // sum-pt of charged Particles(inludes e/mu)			     */
  /* float          Muon_sumNeutralHadronEt03(unsigned int i){return Ntp->Muon_sumNeutralHadronEt03->at(i);}			      // sum pt of neutral hadrons					     */
  /* float          Muon_sumNeutralHadronEtHighThreshold03(unsigned int i){return Ntp->Muon_sumNeutralHadronEtHighThreshold03->at(i);}   // sum pt of neutral hadrons with a higher threshold		     */
  /* float          Muon_sumPhotonEt03(unsigned int i){return Ntp->Muon_sumPhotonEt03->at(i);}					      // sum pt of PF photons					     */
  /* float          Muon_sumPhotonEtHighThreshold03(unsigned int i){return Ntp->Muon_sumPhotonEtHighThreshold03->at(i);}		      // sum pt of PF photons with a higher threshold		     */
  /* float          Muon_sumPUPt03(unsigned int i){return Ntp->Muon_sumPUPt03->at(i);}						      // sum pt of charged Particles not from PV (for Pu corrections)   */
																      								    
  /* float          Muon_sumChargedHadronPt04(unsigned int i){return Ntp->Muon_sumChargedHadronPt04->at(i);}			      // sum-pt of charged Hadron					     */
  /* float          Muon_sumChargedParticlePt04(unsigned int i){return Ntp->Muon_sumChargedParticlePt04->at(i);}			      // sum-pt of charged Particles(inludes e/mu)			     */
  /* float          Muon_sumNeutralHadronEt04(unsigned int i){return Ntp->Muon_sumNeutralHadronEt04->at(i);}			      // sum pt of neutral hadrons					     */
  /* float          Muon_sumNeutralHadronEtHighThreshold04(unsigned int i){return Ntp->Muon_sumNeutralHadronEtHighThreshold04->at(i);}   // sum pt of neutral hadrons with a higher threshold		     */
  /* float          Muon_sumPhotonEt04(unsigned int i){return Ntp->Muon_sumPhotonEt04->at(i);}					      // sum pt of PF photons					     */
  /* float          Muon_sumPhotonEtHighThreshold04(unsigned int i){return Ntp->Muon_sumPhotonEtHighThreshold04->at(i);}		      // sum pt of PF photons with a higher threshold		     */
  /* float          Muon_sumPUPt04(unsigned int i){return Ntp->Muon_sumPUPt04->at(i);}						      // sum pt of charged Particles not from PV (for Pu corrections)   */

  /* int            Muon_numberofValidPixelHits(unsigned int i){return Ntp->Muon_numberofValidPixelHits->at(i);} */
  /* int            Muon_trackerLayersWithMeasurement(unsigned int i){return Ntp->Muon_trackerLayersWithMeasurement->at(i);} */

 bool Muon_TrackParticleHasMomentum(unsigned int i){if(Ntp->Muon_par->at(i).size()!=0)return true; return false;}
   /* TrackParticle Muon_TrackParticle(unsigned int i){ */
   /*   TMatrixT<double>    mu_par(TrackParticle::NHelixPar,1); */
   /*   TMatrixTSym<double> mu_cov(TrackParticle::NHelixPar); */
   /*   unsigned int l=0; */
   /*   for(int k=0; k<TrackParticle::NHelixPar; k++){ */
   /*     mu_par(k,0)=Ntp->Muon_par->at(i).at(k); */
   /*     for(int j=k; j<TrackParticle::NHelixPar; j++){ */
   /*       mu_cov(k,j)=Ntp->Muon_cov->at(i).at(l); */
   /*       l++; */
   /*     } */
   /*   } */
   /*   return TrackParticle(mu_par,mu_cov,Ntp->Muon_pdgid->at(i),Ntp->Muon_M->at(i),Ntp->Muon_trackCharge->at(i),Ntp->Muon_B->at(i)); */
   /* } */

  /* bool           isGoodMuon(unsigned int i); */

   bool           isElectron( int i);
   bool           isMuon( int i);
   bool           isLooseGoodMuon( int i);
   bool           isSoftGoodMuon( int i);
   bool           isMediumGoodMuon( int i);
   bool           isTightGoodMuon( int i);

   bool           isTau( int i);
   bool           isLooseGoodTau( int i);
   bool           isMediumGoodTau( int i);
   bool           isTightGoodTau( int i);

   bool           isIsolatedTau(int i, TString isotype);


   bool           tauBaselineSelection(int i, int genmatch, double cutPt, double cutEta, int aele, int amu, int ajet, string TES="Nom");
   bool           muonBaselineSelection(int i, float ptMin, float etaMax, int muWP);
   bool           electronBaselineSelection( int i, double cutPt, double cutEta);

   bool           tauBaseline (int iDau, float ptMin, float etaMax, int againstEleWP, int againstMuWP, float isoRaw3Hits, TString whatApply, bool debug);
   bool           muBaseline (int iDau, float ptMin, float etaMax, float relIso, int muIDWP, TString whatApply, bool debug);
   bool           eleBaseline (int iDau, float ptMin, float etaMax, float relIso, int MVAIDflag, TString whatApply, bool debug);
   int            getPairType (int type1, int type2);

   bool           MuonVeto(unsigned int i);
   bool           ElectronVeto(unsigned int i);
   bool           DiMuonVeto();
   bool           DiEleVeto();
   
   
   static  bool ComparePairsbyPt(TLorentzVector i, TLorentzVector j);
   std::vector<int>  SortTauHTauHPair (std::vector<int>  PairIndices);
   std::vector<int>  SortPair(std::vector<int>  PairIndices,  std::vector<int>  PairsIndex1, std::vector<int>  PairsIndex2);

   bool           pairPassBaseline (int iPair, TString whatApply, bool debug);
   float          DeltaRDau(int dau1idx, int dau2idx);

  /* bool			 isTightMuon(unsigned int i); */
  /* bool			 isTightMuon(unsigned int i, unsigned int j, TString corr = "default"); */
  /* bool           isSelectedMuon(unsigned int i, unsigned int ja, double impact_xy, double impact_z, TString corr = "default"); */
  /* bool			 isLooseMuon(unsigned int i); */
  /* float          Muon_RelIso(unsigned int i, TString corr = "default"); */

  /* //Base Tau Information (PF) */
  /*  unsigned int      NPFTaus(){return Ntp->PFTau_p4->size();} */
  /*  TLorentzVector	 PFTau_p4(unsigned int i, TString corr = "default"); */
  /*  TVector3          PFTau_Poca(unsigned int i){return TVector3(Ntp->PFTau_Poca->at(i).at(0),Ntp->PFTau_Poca->at(i).at(1),Ntp->PFTau_Poca->at(i).at(2));} */
  /*  bool PFTau_isTightIsolation(unsigned int i){return Ntp->PFTau_isTightIsolation->at(i);} */
  /*  bool PFTau_isMediumIsolation(unsigned int i){return  Ntp->PFTau_isMediumIsolation->at(i);} */
  /*  bool PFTau_isLooseIsolation(unsigned int i){return  Ntp->PFTau_isLooseIsolation->at(i);} */
  /*  bool PFTau_isTightIsolationDBSumPtCorr(unsigned int i){return  Ntp->PFTau_isTightIsolationDBSumPtCorr->at(i);} */
  /*  bool PFTau_isMediumIsolationDBSumPtCorr(unsigned int i){return  Ntp->PFTau_isMediumIsolationDBSumPtCorr->at(i);} */
  /*  bool PFTau_isLooseIsolationDBSumPtCorr(unsigned int i){return  Ntp->PFTau_isLooseIsolationDBSumPtCorr->at(i);} */
  /*  bool PFTau_isVLooseIsolationDBSumPtCorr(unsigned int i){return  Ntp->PFTau_isVLooseIsolationDBSumPtCorr->at(i);} */
  /*  bool PFTau_isHPSAgainstElectronsLoose(unsigned int i){return  Ntp->PFTau_isHPSAgainstElectronsLoose->at(i);} */
  /*  bool PFTau_isHPSAgainstElectronsMedium(unsigned int i){return  Ntp->PFTau_isHPSAgainstElectronsMedium->at(i);} */
  /*  bool PFTau_isHPSAgainstElectronsTight(unsigned int i){return  Ntp->PFTau_isHPSAgainstElectronsTight->at(i);} */
  /*  bool PFTau_isHPSAgainstMuonLoose(unsigned int i){return  Ntp->PFTau_isHPSAgainstMuonLoose->at(i);} */
  /*  bool PFTau_isHPSAgainstMuonMedium(unsigned int i){return  Ntp->PFTau_isHPSAgainstMuonMedium->at(i);} */
  /*  bool PFTau_isHPSAgainstMuonTight(unsigned int i){return  Ntp->PFTau_isHPSAgainstMuonTight->at(i);} */
  /*  bool PFTau_isHPSAgainstMuonLoose2(unsigned int i){return  Ntp->PFTau_isHPSAgainstMuonLoose2->at(i);} */
  /*  bool PFTau_isHPSAgainstMuonMedium2(unsigned int i){return  Ntp->PFTau_isHPSAgainstMuonMedium2->at(i);} */
  /*  bool PFTau_isHPSAgainstMuonTight2(unsigned int i){return  Ntp->PFTau_isHPSAgainstMuonTight2->at(i);} */
  /*  bool PFTau_isHPSByDecayModeFinding(unsigned int i){return  Ntp->PFTau_isHPSByDecayModeFinding->at(i);} */
  /*  bool PFTau_HPSPFTauDiscriminationByMVA3LooseElectronRejection(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByMVA3LooseElectronRejection->at(i);} */
  /*  bool PFTau_HPSPFTauDiscriminationByMVA3MediumElectronRejection(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByMVA3MediumElectronRejection->at(i);} */
  /*  bool PFTau_HPSPFTauDiscriminationByMVA3TightElectronRejection(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByMVA3TightElectronRejection->at(i);} */
  /*  bool PFTau_HPSPFTauDiscriminationByMVA3VTightElectronRejection(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByMVA3VTightElectronRejection->at(i);} */
  /*  bool PFTau_HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits->at(i);} */
  /*  bool PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits->at(i);} */
  /*  bool PFTau_HPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits->at(i);} */
  /*  float PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits(unsigned int i){return Ntp->PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits->at(i);} */
  /*  bool PFTau_HPSPFTauDiscriminationByLooseIsolationMVA(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByLooseIsolationMVA->at(i);} */
  /*  bool PFTau_HPSPFTauDiscriminationByMediumIsolationMVA(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByMediumIsolationMVA->at(i);} */
  /*  bool PFTau_HPSPFTauDiscriminationByTightIsolationMVA(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByTightIsolationMVA->at(i);} */
  /*  bool PFTau_HPSPFTauDiscriminationByLooseIsolationMVA2(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByLooseIsolationMVA2->at(i);} */
  /*  bool PFTau_HPSPFTauDiscriminationByMediumIsolationMVA2(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByMediumIsolationMVA2->at(i);} */
  /*  bool PFTau_HPSPFTauDiscriminationByTightIsolationMVA2(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByTightIsolationMVA2->at(i);} */
  /*  int PFTau_hpsDecayMode(unsigned int i){return  Ntp->PFTau_hpsDecayMode->at(i);} */
  /*  int PFTau_Charge(unsigned int i){return  Ntp->PFTau_Charge->at(i);} */
  /*  std::vector<int> PFTau_Track_idx(unsigned int i){return  Ntp->PFTau_Track_idx->at(i);} */
    
  /*  TMatrixTSym<double> PFTau_TIP_primaryVertex_cov(unsigned int i); */
  /*  bool PFTau_TIP_hassecondaryVertex(unsigned int i){if(Ntp->PFTau_TIP_secondaryVertex_pos->at(i).size()==3)return true; return false;} */
  /*  TVector3 PFTau_TIP_secondaryVertex_pos(unsigned int i){return  TVector3(Ntp->PFTau_TIP_secondaryVertex_pos->at(i).at(0),Ntp->PFTau_TIP_secondaryVertex_pos->at(i).at(1),Ntp->PFTau_TIP_secondaryVertex_pos->at(i).at(2));} */
  /*  TMatrixTSym<double> PFTau_TIP_secondaryVertex_cov(unsigned int i); */
  /*  double PFTau_TIP_secondaryVertex_vtxchi2(unsigned int i){if(Ntp->PFTau_TIP_secondaryVertex_vtxchi2->at(i).size()==1) return  Ntp->PFTau_TIP_secondaryVertex_vtxchi2->at(i).at(0); return 0;} */
  /*  double PFTau_TIP_secondaryVertex_vtxndof(unsigned int i){if(Ntp->PFTau_TIP_secondaryVertex_vtxndof->at(i).size()==1) return  Ntp->PFTau_TIP_secondaryVertex_vtxndof->at(i).at(0);  return 0;} */
    bool PFTau_TIP_hasA1Momentum(unsigned int i){if(Ntp->PFTau_a1_lvp->at(i).size()==LorentzVectorParticle::NLorentzandVertexPar)return true; return false;}
    LorentzVectorParticle PFTau_a1_lvp(unsigned int i);
    TLorentzVector PFTau_3PS_A1_LV(unsigned int i){return PFTau_a1_lvp(i).LV();}
  /*  std::vector<TrackParticle> PFTau_daughterTracks(unsigned int i); */
  /*  std::vector<TVector3> PFTau_daughterTracks_poca(unsigned int i);    */
  /*  TMatrixTSym<double> PFTau_FlightLength3d_cov(unsigned int i){return  PFTau_TIP_secondaryVertex_cov(i)+PFTau_TIP_primaryVertex_cov(i);} */
  /*  TVector3 PFTau_FlightLength3d(unsigned int i){return PFTau_TIP_secondaryVertex_pos(i)-PFTau_TIP_primaryVertex_pos(i);} */
 
  /*  double	PFTau_FlightLength_significance(unsigned int i); */
  /*  double   PFTau_FlightLength(unsigned int i){return PFTau_FlightLength3d(i).Mag();} */

  /*     // Jet Information */
  /*  unsigned int       NPFJets(){return Ntp->PFJet_p4->size();} */
  /*  TLorentzVector     PFJet_p4(unsigned int i, TString corr = "default"); */
  /*  float              PFJet_chargedEmEnergy(unsigned int i){return Ntp->PFJet_chargedEmEnergy->at(i);} */
  /*  float              PFJet_chargedHadronEnergy(unsigned int i){return Ntp->PFJet_chargedHadronEnergy->at(i);} */
  /*  int	              PFJet_chargedHadronMultiplicity(unsigned int i){return Ntp->PFJet_chargedHadronMultiplicity->at(i);} */
  /*  float              PFJet_chargedMuEnergy(unsigned int i){return Ntp->PFJet_chargedMuEnergy->at(i);} */
  // float	              PFJet_chargedMultiplicity(unsigned int i){return Ntp->_jets_chMult->at(i);}
  // float	              PFJet_neutralMultiplicity(unsigned int i){return Ntp->_jets_neMult->at(i);}
    float                NumConst(unsigned int i){return (jets_chMult(i) +jets_neMult(i));}
  /*  float              PFJet_electronEnergy(unsigned int i){return Ntp->PFJet_electronEnergy->at(i);} */
  /*  int	              PFJet_electronMultiplicity(unsigned int i){return Ntp->PFJet_electronMultiplicity->at(i);} */
  /*  float              PFJet_HFEMEnergy(unsigned int i){return Ntp->PFJet_HFEMEnergy->at(i);} */
  /*  int	              PFJet_HFEMMultiplicity(unsigned int i){return Ntp->PFJet_HFEMMultiplicity->at(i);} */
  /*  float              PFJet_HFHadronEnergy(unsigned int i){return Ntp->PFJet_HFHadronEnergy->at(i);} */
  /*  int	              PFJet_HFHadronMultiplicity(unsigned int i){return Ntp->PFJet_HFHadronMultiplicity->at(i);} */
  /*  float              PFJet_muonEnergy(unsigned int i){return Ntp->PFJet_muonEnergy->at(i);} */
  /*  int	              PFJet_muonMultiplicity(unsigned int i){return Ntp->PFJet_muonMultiplicity->at(i);} */
  /*  float              PFJet_neutralEmEnergy(unsigned int i){return Ntp->PFJet_neutralEmEnergy->at(i);} */
  /*  float              PFJet_neutralHadronEnergy(unsigned int i){return Ntp->PFJet_neutralHadronEnergy->at(i);} */
  /*  int	              PFJet_neutralHadronMultiplicity(unsigned int i){return Ntp->PFJet_neutralHadronMultiplicity->at(i);} */
  /*  float              PFJet_photonEnergy(unsigned int i){return Ntp->PFJet_photonEnergy->at(i);} */
  /*  int	              PFJet_photonMultiplicity(unsigned int i){return Ntp->PFJet_photonMultiplicity->at(i);} */
  /*  float              PFJet_jetArea(unsigned int i){return Ntp->PFJet_jetArea->at(i);} */
  /*  float              PFJet_maxDistance(unsigned int i){return Ntp->PFJet_maxDistance->at(i);} */
  /*  int                PFJet_nConstituents(unsigned int i){return Ntp->PFJet_nConstituents->at(i);} */
  /*  float              PFJet_pileup(unsigned int i){return Ntp->PFJet_pileup->at(i);} */
  /*  float              PFJet_etaetaMoment(unsigned int i){return Ntp->PFJet_etaetaMoment->at(i);} */
  /*  float              PFJet_etaphiMoment(unsigned int i){return Ntp->PFJet_etaphiMoment->at(i);} */
  /*  std::vector<int>   PFJet_Track_idx(unsigned int i){return Ntp->PFJet_Track_idx->at(i);} */
  /*  int                PFJet_MatchedHPS_idx(unsigned int i){return Ntp->PFJet_MatchedHPS_idx->at(i);} */
  /*  int                PFJet_numberOfDaughters(unsigned int i){return Ntp->PFJet_numberOfDaughters->at(i);} */
  //  float              PFJet_chargedEmEnergyFraction(unsigned int i){return Ntp->_jets_chEmEF->at(i);}
  //  float              PFJet_chargedHadronEnergyFraction(unsigned int i){return Ntp->_jets_chHEF->at(i);}
  //  float              PFJet_neutralHadronEnergyFraction(unsigned int i){return Ntp->_jets_nHEF->at(i);}
  //  float              PFJet_neutralEmEnergyFraction(unsigned int i){return Ntp->_jets_nEmEF->at(i);}
  /*  bool               isGoodJet(unsigned int i); */
  /*  bool               isGoodJet_nooverlapremoval(unsigned int i); */
  /*  bool               isJetID(unsigned int i, TString corr = "default"); */
  /*  int	              PFJet_nTrk(unsigned int i){return Ntp->PFJet_nTrk->at(i);} */
  /*  TLorentzVector     PFJet_TracksP4(unsigned int i, unsigned int j){return TLorentzVector(Ntp->PFJet_TracksP4->at(i).at(j).at(1),Ntp->PFJet_TracksP4->at(i).at(j).at(2),Ntp->PFJet_TracksP4->at(i).at(j).at(3),Ntp->PFJet_TracksP4->at(i).at(j).at(0));} */
  /*  int                PFJet_nTracks(unsigned int i){return Ntp->PFJet_TracksP4->at(i).size();} */
  /*  float			  PFJet_JECuncertainty(unsigned int i){return Ntp->PFJet_JECuncertainty->at(i);} */
  /*  unsigned int       PFJet_NGenJets(){return Ntp->PFJet_GenJet_p4->size();} */
  /*  unsigned int       PFJet_NGenJetsNoNu(){return Ntp->PFJet_GenJetNoNu_p4->size();} */
  /*  TLorentzVector     PFJet_GenJet_p4(unsigned int i){return TLorentzVector(Ntp->PFJet_GenJet_p4->at(i).at(1),Ntp->PFJet_GenJet_p4->at(i).at(2),Ntp->PFJet_GenJet_p4->at(i).at(3),Ntp->PFJet_GenJet_p4->at(i).at(0));} */
  /*  TLorentzVector     PFJet_GenJet_Constituents_p4(unsigned int i, unsigned int j){return TLorentzVector(Ntp->PFJet_GenJet_Constituents_p4->at(i).at(j).at(1),Ntp->PFJet_GenJet_Constituents_p4->at(i).at(j).at(2),Ntp->PFJet_GenJet_Constituents_p4->at(i).at(j).at(3),Ntp->PFJet_GenJet_Constituents_p4->at(i).at(j).at(0));} */
  /*  TLorentzVector     PFJet_GenJetNoNu_p4(unsigned int i){return TLorentzVector(Ntp->PFJet_GenJetNoNu_p4->at(i).at(1),Ntp->PFJet_GenJetNoNu_p4->at(i).at(2),Ntp->PFJet_GenJetNoNu_p4->at(i).at(3),Ntp->PFJet_GenJetNoNu_p4->at(i).at(0));} */
  /*  TLorentzVector     PFJet_GenJetNoNu_Constituents_p4(unsigned int i, unsigned int j){return TLorentzVector(Ntp->PFJet_GenJetNoNu_Constituents_p4->at(i).at(j).at(1),Ntp->PFJet_GenJetNoNu_Constituents_p4->at(i).at(j).at(2),Ntp->PFJet_GenJetNoNu_Constituents_p4->at(i).at(j).at(3),Ntp->PFJet_GenJetNoNu_Constituents_p4->at(i).at(j).at(0));} */

  /*  float              PFJet_PUJetID_discr(unsigned int i){return Ntp->PFJet_PUJetID_discr->at(i);} */
  /*  bool	              PFJet_PUJetID_looseWP(unsigned int i){return Ntp->PFJet_PUJetID_looseWP->at(i);} */
  /*  bool	              PFJet_PUJetID_mediumWP(unsigned int i){return Ntp->PFJet_PUJetID_mediumWP->at(i);} */
  /*  bool	              PFJet_PUJetID_tightWP(unsigned int i){return Ntp->PFJet_PUJetID_tightWP->at(i);} */

  /*  int	              PFJet_partonFlavour(unsigned int i){return Ntp->PFJet_partonFlavour->at(i);} */
  /*  float              PFJet_bDiscriminator(unsigned int i){return Ntp->PFJet_bDiscriminator->at(i);} */
  /*  //float              PFJet_BTagWeight(unsigned int i){return Ntp->PFJet_BTagWeight->at(i);} // not implemented at the moment */

  /*  double 			  rundependentJetPtCorrection(double jeteta, int runnumber); */
  /*  double             JERCorrection(TLorentzVector jet, double dr=0.25, TString corr = "default"); // dr=0.25 from AN2013_416_v4 */
  /*  TLorentzVector     PFJet_matchGenJet(TLorentzVector jet, double dr); */
  /*  double             JetEnergyResolutionCorr(double jeteta); */
  /*  double             JetEnergyResolutionCorrErr(double jeteta); */


  /*  //  double             MET_CorrCaloT1T2_ey(){return Ntp->MET_CorrCaloT1T2_et*sin(Ntp->MET_CorrCaloT1T2_phi);} */

  /*  double             MET_CorrMVA_et(){return Ntp->MET_CorrMVA_et;} */
  /*  double             MET_CorrMVA_phi(){return Ntp->MET_CorrMVA_phi;} */
  /*  double             MET_CorrMVA_ex(){return Ntp->MET_CorrMVA_et*cos(Ntp->MET_CorrMVA_phi);} */
  /*  double             MET_CorrMVA_ey(){return Ntp->MET_CorrMVA_et*sin(Ntp->MET_CorrMVA_phi);} */
  /*  double             MET_CorrMVA_significance(){return Ntp->MET_CorrMVA_significance;} */
  /*  double             MET_CorrMVA_significance_xx(){return Ntp->MET_CorrMVA_significance_xx;} */
  /*  double             MET_CorrMVA_significance_xy(){return Ntp->MET_CorrMVA_significance_xy;} */
  /*  double             MET_CorrMVA_significance_yy(){return Ntp->MET_CorrMVA_significance_yy;} */
  /*  unsigned int   	  NMET_CorrMVA_srcMuons(){return Ntp->MET_CorrMVA_srcMuon_p4->size();} */
  /*  TLorentzVector     MET_CorrMVA_srcMuon_p4(unsigned int i){return TLorentzVector(Ntp->MET_CorrMVA_srcMuon_p4->at(i).at(1),Ntp->MET_CorrMVA_srcMuon_p4->at(i).at(2),Ntp->MET_CorrMVA_srcMuon_p4->at(i).at(3),Ntp->MET_CorrMVA_srcMuon_p4->at(i).at(0));} */
  /*  bool				  findCorrMVASrcMuon(unsigned int muon_idx, int &mvaSrcMuon_idx, float &dR ); */
  /*  unsigned int   	  NMET_CorrMVA_srcElectrons(){return Ntp->MET_CorrMVA_srcElectron_p4->size();} */
  /*  TLorentzVector     MET_CorrMVA_srcElectron_p4(unsigned int i){return TLorentzVector(Ntp->MET_CorrMVA_srcElectron_p4->at(i).at(1),Ntp->MET_CorrMVA_srcElectron_p4->at(i).at(2),Ntp->MET_CorrMVA_srcElectron_p4->at(i).at(3),Ntp->MET_CorrMVA_srcElectron_p4->at(i).at(0));} */
  /*  bool				  findCorrMVASrcElectron(unsigned int elec_idx, int &mvaSrcElectron_idx, float &dR ); */
  /*  unsigned int   	  NMET_CorrMVA_srcTaus(){return Ntp->MET_CorrMVA_srcTau_p4->size();} */
  /*  TLorentzVector     MET_CorrMVA_srcTau_p4(unsigned int i){return TLorentzVector(Ntp->MET_CorrMVA_srcTau_p4->at(i).at(1),Ntp->MET_CorrMVA_srcTau_p4->at(i).at(2),Ntp->MET_CorrMVA_srcTau_p4->at(i).at(3),Ntp->MET_CorrMVA_srcTau_p4->at(i).at(0));} */
  /*  bool				  findCorrMVASrcTau(unsigned int tau_idx, int &mvaSrcTau_idx, float &dR ); */

  /*  double             MET_CorrMVAMuTau_et(){return Ntp->MET_CorrMVAMuTau_et;} */
  /*  double             MET_CorrMVAMuTau_phi(){return Ntp->MET_CorrMVAMuTau_phi;} */
  /*  double             MET_CorrMVAMuTau_ex(){return Ntp->MET_CorrMVAMuTau_et*cos(Ntp->MET_CorrMVAMuTau_phi);} */
  /*  double             MET_CorrMVAMuTau_ey(){return Ntp->MET_CorrMVAMuTau_et*sin(Ntp->MET_CorrMVAMuTau_phi);} */
  /*  double             MET_CorrMVAMuTau_significance(){return Ntp->MET_CorrMVAMuTau_significance;} */
  /*  double             MET_CorrMVAMuTau_significance_xx(){return Ntp->MET_CorrMVAMuTau_significance_xx;} */
  /*  double             MET_CorrMVAMuTau_significance_xy(){return Ntp->MET_CorrMVAMuTau_significance_xy;} */
  /*  double             MET_CorrMVAMuTau_significance_yy(){return Ntp->MET_CorrMVAMuTau_significance_yy;} */
  /*  unsigned int   	  NMET_CorrMVAMuTau_srcMuons(){return Ntp->MET_CorrMVAMuTau_srcMuon_p4->size();} */
  /*  TLorentzVector     MET_CorrMVAMuTau_srcMuon_p4(unsigned int i){return TLorentzVector(Ntp->MET_CorrMVAMuTau_srcMuon_p4->at(i).at(1),Ntp->MET_CorrMVAMuTau_srcMuon_p4->at(i).at(2),Ntp->MET_CorrMVAMuTau_srcMuon_p4->at(i).at(3),Ntp->MET_CorrMVAMuTau_srcMuon_p4->at(i).at(0));} */
  /*  bool				  findCorrMVAMuTauSrcMuon(unsigned int muon_idx, int &mvaMuTauSrcMuon_idx, float &dR ); */
  /*  unsigned int   	  NMET_CorrMVAMuTau_srcTaus(){return Ntp->MET_CorrMVAMuTau_srcTau_p4->size();} */
  /*  TLorentzVector     MET_CorrMVAMuTau_srcTau_p4(unsigned int i){return TLorentzVector(Ntp->MET_CorrMVAMuTau_srcTau_p4->at(i).at(1),Ntp->MET_CorrMVAMuTau_srcTau_p4->at(i).at(2),Ntp->MET_CorrMVAMuTau_srcTau_p4->at(i).at(3),Ntp->MET_CorrMVAMuTau_srcTau_p4->at(i).at(0));} */
  /*  bool				  findCorrMVAMuTauSrcTau(unsigned int tau_idx, int &mvaMuTauSrcTau_idx, float &dR ); */

  /*  /////////////////////////////// */
  /*  // */
  /*  // MET uncertainties */
  /*  // */
  /*  // Type1      = T0 + T1 + Txy (recommendation by JetMET POG and Christian Veelken) */
  /*  // Type1Type2 = T0 + T1 + Txy + calibration for unclustered energy (better MET response, worse MET resolution) */
  /*  // */
  /*  double             MET_Type1CorrElectronUp_et(){ return Ntp->MET_Type1CorrElectronUp_et; } */
  /*  double             MET_Type1CorrElectronDown_et(){ return Ntp->MET_Type1CorrElectronDown_et; } */
  /*  double             MET_Type1CorrMuonUp_et(){ return Ntp->MET_Type1CorrMuonUp_et; } */
  /*  double             MET_Type1CorrMuonDown_et(){ return Ntp->MET_Type1CorrMuonDown_et; } */
  /*  double             MET_Type1CorrTauUp_et(){ return Ntp->MET_Type1CorrTauUp_et; } */
  /*  double             MET_Type1CorrTauDown_et(){ return Ntp->MET_Type1CorrTauDown_et; } */
  /*  double             MET_Type1CorrJetResUp_et(){ return Ntp->MET_Type1CorrJetResUp_et; } */
  /*  double             MET_Type1CorrJetResDown_et(){ return Ntp->MET_Type1CorrJetResDown_et; } */
  /*  double             MET_Type1CorrJetEnUp_et(){ return Ntp->MET_Type1CorrJetEnUp_et; } */
  /*  double             MET_Type1CorrJetEnDown_et(){ return Ntp->MET_Type1CorrJetEnDown_et; } */
  /*  double             MET_Type1CorrUnClusteredUp_et(){ return Ntp->MET_Type1CorrUnclusteredUp_et; } */
  /*  double             MET_Type1CorrUnClusteredDown_et(){ return Ntp->MET_Type1CorrUnclusteredDown_et; } */
  /*  double             MET_Type1p2CorrElectronUp_et(){ return Ntp->MET_Type1p2CorrElectronUp_et; } */
  /*  double             MET_Type1p2CorrElectronDown_et(){ return Ntp->MET_Type1p2CorrElectronDown_et; } */
  /*  double             MET_Type1p2CorrMuonUp_et(){ return Ntp->MET_Type1p2CorrMuonUp_et; } */
  /*  double             MET_Type1p2CorrMuonDown_et(){ return Ntp->MET_Type1p2CorrMuonDown_et; } */
  /*  double             MET_Type1p2CorrTauUp_et(){ return Ntp->MET_Type1p2CorrTauUp_et; } */
  /*  double             MET_Type1p2CorrTauDown_et(){ return Ntp->MET_Type1p2CorrTauDown_et; } */
  /*  double             MET_Type1p2CorrJetResUp_et(){ return Ntp->MET_Type1p2CorrJetResUp_et; } */
  /*  double             MET_Type1p2CorrJetResDown_et(){ return Ntp->MET_Type1p2CorrJetResDown_et; } */
  /*  double             MET_Type1p2CorrJetEnUp_et(){ return Ntp->MET_Type1p2CorrJetEnUp_et; } */
  /*  double             MET_Type1p2CorrJetEnDown_et(){ return Ntp->MET_Type1p2CorrJetEnDown_et; } */
  /*  double             MET_Type1p2CorrUnclusteredUp_et(){ return Ntp->MET_Type1p2CorrUnclusteredUp_et; } */
  /*  double             MET_Type1p2CorrUnclusteredDown_et(){ return Ntp->MET_Type1p2CorrUnclusteredDown_et; } */

  /*  //Track Information */
  /*  unsigned int      NTracks(){return Ntp->Track_p4->size();} */
  /*  TLorentzVector    Track_p4(unsigned int i){return TLorentzVector(Ntp->Track_p4->at(i).at(1),Ntp->Track_p4->at(i).at(2),Ntp->Track_p4->at(i).at(3),Ntp->Track_p4->at(i).at(0));} */
  /*  TVector3          Track_Poca(unsigned int i){return TVector3(Ntp->Track_Poca->at(i).at(0),Ntp->Track_Poca->at(i).at(1),Ntp->Track_Poca->at(i).at(2));} */
  /*  int               Track_charge(unsigned int i){return Ntp->Track_charge->at(i);} */
  /*  double            Track_chi2(unsigned int i){return Ntp->Track_chi2->at(i);} */
  /*  double            Track_ndof(unsigned int i){return Ntp->Track_ndof->at(i);} */
  /*  unsigned short    Track_numberOfLostHits(unsigned int i){return Ntp->Track_numberOfLostHits->at(i);} */
  /*  unsigned short    Track_numberOfValidHits(unsigned int i){return Ntp->Track_numberOfValidHits->at(i);} */
  /*  unsigned int      Track_qualityMask(unsigned int i){return Ntp->Track_qualityMask->at(i);} */

  /*  TrackParticle Track_TrackParticle(unsigned int i){ */
  /*    TMatrixT<double>    track_par(TrackParticle::NHelixPar,1); */
  /*    TMatrixTSym<double> track_cov(TrackParticle::NHelixPar); */
  /*    unsigned int l=0; */
  /*    for(int k=0; k<TrackParticle::NHelixPar; k++){ */
  /*      track_par(k,0)=Ntp->Track_par->at(i).at(k); */
  /*      for(int j=k; j<TrackParticle::NHelixPar; j++){ */
  /* 	 track_cov(k,j)=Ntp->Track_cov->at(i).at(l); */
  /* 	 l++; */
  /*      } */
  /*    } */
  /*    return TrackParticle(track_par,track_cov,Ntp->Track_pdgid->at(i),Ntp->Track_M->at(i),Ntp->Track_charge->at(i),Ntp->Track_B->at(i)); */
  /*  } */

   // MC Information
   // Signal particles (Z0,W+/-,H0,H+/-)

   /* bool hasSignalTauDecay(PDGInfo::PDGMCNumbering parent_pdgid,unsigned int &Boson_idx,TauDecay::JAK tau_jak, unsigned int &idx); */
   /* bool hasSignalTauDecay(PDGInfo::PDGMCNumbering parent_pdgid,unsigned int &Boson_idx,unsigned int &tau1_idx, unsigned int &tau2_idx); */


   /* // overlap of objects */
   /* bool jethasMuonOverlap(unsigned int jet_idx,unsigned int &muon_idx); */
   /* bool muonhasJetOverlap(unsigned int muon_idx,unsigned int &jet_idx); */
   /* bool muonhasJetMatch(unsigned int muon_idx,unsigned int &jet_idx); */


   /* // Electrons */
   /* unsigned int       NElectrons(){return Ntp->Electron_p4->size();} */
   /* TLorentzVector	  Electron_p4(unsigned int i, TString corr = "default"); */
   /* TVector3           Electron_Poca(unsigned int i){return TVector3(Ntp->Electron_Poca->at(i).at(0),Ntp->Electron_Poca->at(i).at(1),Ntp->Electron_Poca->at(i).at(2));} */
   /* int   Electron_Charge(unsigned int i){return Ntp->Electron_charge->at(i);} */
   /* float   Electron_Gsf_deltaEtaEleClusterTrackAtCalo(unsigned int i){return Ntp->Electron_Gsf_deltaEtaEleClusterTrackAtCalo->at(i);} */
   /* float   Electron_Gsf_deltaEtaSeedClusterTrackAtCalo(unsigned int i){return Ntp->Electron_Gsf_deltaEtaSeedClusterTrackAtCalo->at(i);} */
   /* float   Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(unsigned int i){return Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx->at(i);} */
   /* float   Electron_Gsf_deltaPhiEleClusterTrackAtCalo(unsigned int i){return Ntp->Electron_Gsf_deltaPhiEleClusterTrackAtCalo->at(i);} */
   /* float   Electron_Gsf_deltaPhiSeedClusterTrackAtCalo(unsigned int i){return Ntp->Electron_Gsf_deltaPhiSeedClusterTrackAtCalo->at(i);} */
   /* float   Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(unsigned int i){return Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx->at(i);} */
   /* float   Electron_Gsf_dr03EcalRecHitSumE(unsigned int i){return Ntp->Electron_Gsf_dr03EcalRecHitSumE->at(i);} */
   /* float   Electron_Gsf_dr03HcalDepth1TowerSumEt(unsigned int i){return Ntp->Electron_Gsf_dr03HcalDepth1TowerSumEt->at(i);} */
   /* float   Electron_Gsf_dr03HcalDepth1TowerSumEtBc(unsigned int i){return Ntp->Electron_Gsf_dr03HcalDepth1TowerSumEtBc->at(i);} */
   /* float   Electron_Gsf_dr03HcalDepth2TowerSumEt(unsigned int i){return Ntp->Electron_Gsf_dr03HcalDepth2TowerSumEt->at(i);} */
   /* float   Electron_Gsf_dr03HcalDepth2TowerSumEtBc(unsigned int i){return Ntp->Electron_Gsf_dr03HcalDepth2TowerSumEtBc->at(i);} */
   /* float   Electron_Gsf_dr03HcalTowerSumEt(unsigned int i){return Ntp->Electron_Gsf_dr03HcalTowerSumEt->at(i);} */
   /* float   Electron_Gsf_dr03HcalTowerSumEtBc(unsigned int i){return Ntp->Electron_Gsf_dr03HcalTowerSumEtBc->at(i);} */
   /* float   Electron_Gsf_dr03TkSumPt(unsigned int i){return Ntp->Electron_Gsf_dr03TkSumPt->at(i);} */
   /* bool    Electron_Gsf_passingCutBasedPreselection(unsigned int i){return Ntp->Electron_Gsf_passingCutBasedPreselection->at(i);} */
   /* bool    Electron_Gsf_passingMvaPreselection(unsigned int i){return Ntp->Electron_Gsf_passingMvaPreselection->at(i);} */
   /* int     Electron_gsftrack_trackerExpectedHitsInner_numberOfLostHits(unsigned int i){return Ntp->Electron_gsftrack_trackerExpectedHitsInner_numberOfLostHits->at(i);} */
   /* double  Electron_supercluster_e(unsigned int i){return Ntp->Electron_supercluster_e->at(i);} */
   /* double  Electron_supercluster_phi(unsigned int i){return Ntp->Electron_supercluster_phi->at(i);} */
   /* double  Electron_supercluster_eta(unsigned int i){return Ntp->Electron_supercluster_eta->at(i);} */
   /* float   Electron_supercluster_centroid_x(unsigned int i){return Ntp->Electron_supercluster_centroid_x->at(i);} */
   /* float   Electron_supercluster_centroid_y(unsigned int i){return Ntp->Electron_supercluster_centroid_y->at(i);} */
   /* float   Electron_supercluster_centroid_z(unsigned int i){return Ntp->Electron_supercluster_centroid_z->at(i);} */
   /* unsigned int Electron_Track_idx(unsigned int i){return Ntp->Electron_Track_idx->at(i);} */

   /* float    Electron_ecalRecHitSumEt03(unsigned int i){return Ntp->Electron_ecalRecHitSumEt03->at(i);} */
   /* float    Electron_hcalDepth1TowerSumEt03(unsigned int i){return Ntp->Electron_hcalDepth1TowerSumEt03->at(i);} */
   /* float    Electron_hcalDepth1TowerSumEtBc03(unsigned int i){return Ntp->Electron_hcalDepth1TowerSumEtBc03->at(i);} */
   /* float    Electron_hcalDepth2TowerSumEt03(unsigned int i){return Ntp->Electron_hcalDepth2TowerSumEt03->at(i);} */
   /* float    Electron_hcalDepth2TowerSumEtBc03(unsigned int i){return Ntp->Electron_hcalDepth2TowerSumEtBc03->at(i);} */
   /* float    Electron_tkSumPt03(unsigned int i){return Ntp->Electron_tkSumPt03->at(i);} */
   /* float    Electron_ecalRecHitSumEt04(unsigned int i){return Ntp->Electron_ecalRecHitSumEt04->at(i);} */
   /* float    Electron_hcalDepth1TowerSumEt04(unsigned int i){return Ntp->Electron_hcalDepth1TowerSumEt04->at(i);} */
   /* float    Electron_hcalDepth1TowerSumEtBc04(unsigned int i){return Ntp->Electron_hcalDepth1TowerSumEtBc04->at(i);} */
   /* float    Electron_hcalDepth2TowerSumEt04(unsigned int i){return Ntp->Electron_hcalDepth2TowerSumEt04->at(i);} */
   /* float    Electron_hcalDepth2TowerSumEtBc04(unsigned int i){return Ntp->Electron_hcalDepth2TowerSumEtBc04->at(i);} */
   /* float    Electron_tkSumPt04(unsigned int i){return Ntp->Electron_tkSumPt04->at(i);} */
   
   /* float    Electron_chargedHadronIso(unsigned int i){return Ntp->Electron_chargedHadronIso->at(i);} */
   /* float    Electron_neutralHadronIso(unsigned int i){return Ntp->Electron_neutralHadronIso->at(i);} */
   /* float    Electron_photonIso(unsigned int i){return Ntp->Electron_photonIso->at(i);} */
   
   /* double   Electron_isoDeposits_chargedHadronIso04(unsigned int i){return Ntp->Electron_isoDeposits_chargedHadronIso04->at(i);} */
   /* double   Electron_isoDeposits_neutralHadronIso04(unsigned int i){return Ntp->Electron_isoDeposits_neutralHadronIso04->at(i);} */
   /* double   Electron_isoDeposits_photonIso04(unsigned int i){return Ntp->Electron_isoDeposits_photonIso04->at(i);} */
   /* double   Electron_isoDeposits_chargedHadronIso03(unsigned int i){return Ntp->Electron_isoDeposits_chargedHadronIso03->at(i);} */
   /* double   Electron_isoDeposits_neutralHadronIso03(unsigned int i){return Ntp->Electron_isoDeposits_neutralHadronIso03->at(i);} */
   /* double   Electron_isoDeposits_photonIso03(unsigned int i){return Ntp->Electron_isoDeposits_photonIso03->at(i);} */

   /* float    Electron_sigmaIetaIeta(unsigned int i){return Ntp->Electron_sigmaIetaIeta->at(i);} */
   /* float    Electron_hadronicOverEm(unsigned int i){return Ntp->Electron_hadronicOverEm->at(i);} */
   /* float    Electron_fbrem(unsigned int i){return Ntp->Electron_fbrem->at(i);} */
   /* float    Electron_eSuperClusterOverP(unsigned int i){return Ntp->Electron_eSuperClusterOverP->at(i);} */
   /* float    Electron_ecalEnergy(unsigned int i){return Ntp->Electron_ecalEnergy->at(i);} */
   /* float    Electron_trackMomentumAtVtx(unsigned int i){return Ntp->Electron_trackMomentumAtVtx->at(i);} */
   /* int      Electron_numberOfMissedHits(unsigned int i){return Ntp->Electron_numberOfMissedHits->at(i);} */
   /* bool     Electron_HasMatchedConversions(unsigned int i){return Ntp->Electron_HasMatchedConversions->at(i);} */

   /* double   Electron_MVA_Trig_discriminator(unsigned int i){return Ntp->Electron_MVA_Trig_discriminator->at(i);} */
   /* double   Electron_MVA_TrigNoIP_discriminator(unsigned int i){return Ntp->Electron_MVA_TrigNoIP_discriminator->at(i);} */
   /* double   Electron_MVA_NonTrig_discriminator(unsigned int i){return Ntp->Electron_MVA_NonTrig_discriminator->at(i);} */
   /* double   RhoIsolationAllInputTags(){return Ntp->RhoIsolationAllInputTags;} */

   /* double   Electron_RegEnergy(unsigned int i){return Ntp->Electron_RegEnergy->at(i);} */
   /* double   Electron_RegEnergyError(unsigned int i){return Ntp->Electron_RegEnergyError->at(i);} */

   /* TrackParticle Electron_TrackParticle(unsigned int i){ */
   /*   TMatrixT<double>    e_par(TrackParticle::NHelixPar,1); */
   /*   TMatrixTSym<double> e_cov(TrackParticle::NHelixPar); */
   /*   unsigned int l=0; */
   /*   for(int k=0; k<TrackParticle::NHelixPar; k++){ */
   /*     e_par(k,0)=Ntp->Electron_par->at(i).at(k); */
   /*     for(int j=k; j<TrackParticle::NHelixPar; j++){ */
   /* 	 e_cov(k,j)=Ntp->Electron_cov->at(i).at(l); */
   /* 	 l++; */
   /*     } */
   /*   } */
   /*   return TrackParticle(e_par,e_cov,Ntp->Electron_pdgid->at(i),Ntp->Electron_M->at(i),Ntp->Electron_charge->at(i),Ntp->Electron_B->at(i)); */
   /* } */

   /* bool isTrigPreselElectron(unsigned int i); */
   /* bool isTrigNoIPPreselElectron(unsigned int i); */
   /* bool isMVATrigElectron(unsigned int i, TString corr = "default"); */
   /* bool isMVATrigNoIPElectron(unsigned int i, TString corr = "default"); */
   /* bool isMVANonTrigElectron(unsigned int i, unsigned int j, TString corr = "default"); */
   /* bool isTightElectron(unsigned int i, TString corr = "default"); */
   /* bool isTightElectron(unsigned int i, unsigned int j, TString corr = "default"); */
   /* float Electron_RelIso03(unsigned int i, TString corr = "default"); */
   /* double Electron_RelIsoDep03(unsigned int i, TString corr = "default");//todo: better name */
   /* float Electron_RelIso04(unsigned int i, TString corr = "default"); */
   /* double Electron_RelIsoDep04(unsigned int i, TString corr = "default");//todo: better name */
   /* double Electron_Aeff_R04(double Eta); */
   /* double Electron_Aeff_R03(double Eta); */
   /* bool isSelectedElectron(unsigned int i, unsigned int j, double impact_xy, double impact_z, TString corr = "default"); */

   // Trigger
   /* bool         TriggerAccept(TString n); */
   /* unsigned int HLTPrescale(TString n); */
   /* unsigned int L1SEEDPrescale(TString n); */
   /* bool         GetTriggerIndex(TString n, unsigned int &i); */
   /* double 		matchTrigger(TLorentzVector obj, std::vector<TString> trigger, std::string objectType); */
   /* bool 		matchTrigger(TLorentzVector obj, double dr_cut, std::vector<TString> trigger, std::string objectType); */
   /* bool			matchTrigger(TLorentzVector obj, double dr_cut, TString trigger, std::string objectType); */
   /* unsigned int NHLTTriggers(){return Ntp->HTLTriggerName->size();} */
   /* std::string  HTLTriggerName(unsigned int i){return Ntp->HTLTriggerName->at(i);} */
   /* bool         TriggerAccept(unsigned int i){return Ntp->TriggerAccept->at(i);} */
   /* bool         TriggerError(unsigned int i){return Ntp->TriggerError->at(i);} */
   /* bool         TriggerWasRun(unsigned int i){return Ntp->TriggerWasRun->at(i);} */
   /* unsigned int HLTPrescale(unsigned int i){return Ntp->HLTPrescale->at(i);} */
   /* unsigned int NHLTL1GTSeeds(unsigned int i){return Ntp->NHLTL1GTSeeds->at(i);} */
   /* unsigned int L1SEEDPrescale(unsigned int i){return Ntp->L1SEEDPrescale->at(i);} */
   /* bool         L1SEEDInvalidPrescale(unsigned int i){return Ntp->L1SEEDInvalidPrescale->at(i);} */
   /* unsigned int NHLTTriggerObject(unsigned int i){return Ntp->HLTTrigger_objs_Eta->at(i).size();} */
   /* TLorentzVector HLTTriggerObject_p4(unsigned int i, unsigned int j){ */
   /*   TLorentzVector L(0,0,0,0);  */
   /*   if(j<Ntp->HLTTrigger_objs_Eta->at(i).size())L.SetPtEtaPhiM(Ntp->HLTTrigger_objs_Pt->at(i).at(j),Ntp->HLTTrigger_objs_Eta->at(i).at(j), Ntp->HLTTrigger_objs_Phi->at(i).at(j),0.0); */
   /*   return L; */
   /* } */
   /* int          NHLTTrigger_objs(){return Ntp->HLTTrigger_objs_Pt->size();} */
   /* int          NHLTTrigger_objs(unsigned int i){return Ntp->HLTTrigger_objs_Pt->at(i).size();} */
   /* float        HLTTrigger_objs_Pt(unsigned int i, unsigned int j){return Ntp->HLTTrigger_objs_Pt->at(i).at(j);} */
   /* float        HLTTrigger_objs_Eta(unsigned int i, unsigned int j){return Ntp->HLTTrigger_objs_Eta->at(i).at(j);} */
   /* float        HLTTrigger_objs_Phi(unsigned int i, unsigned int j){return Ntp->HLTTrigger_objs_Phi->at(i).at(j);} */
   /* float        HLTTrigger_objs_E(unsigned int i,unsigned int j){return Ntp->HLTTrigger_objs_E->at(i).at(j);} */
   /* int          HLTTrigger_objs_Id(unsigned int i,unsigned int j){return Ntp->HLTTrigger_objs_Id->at(i).at(j);} */
   /* std::string  HLTTrigger_objs_trigger(unsigned int i){return Ntp->HLTTrigger_objs_trigger->at(i);} */

   // helper functions

   double	PFTau_FlightLength_significance(TVector3 pv,TMatrixTSym<double> PVcov, TVector3 sv, TMatrixTSym<double> SVcov );
   double       dxy(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);
   double	    dxySigned(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);
   double       dz(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);
   double       dzSigned(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);
   double       vertexSignificance(TVector3 vec, unsigned int vertex);
   double       transverseMass(double pt1, double phi1, double pt2, double phi2){return sqrt(2 * pt1 * pt2 * (1 - cos(phi1 - phi2)));}
   std::vector<int> sortObjects(std::vector<int> indices, std::vector<double> values);
   std::vector<int> sortPFJetsByPt();
   std::vector<int> sortDefaultObjectsByPt(TString objectType);

   void deb(int i){std::cout<<"  deb   "<<i<<std::endl;}

   double DeltaPhi(double angle1,double angle2);
     


};

#endif

