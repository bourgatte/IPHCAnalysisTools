#include "HCPTauTau.h"
#include "TLorentzVector.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "SVFitObject.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"
//#include "SVfitProvider.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "PDG_Var.h"
#include "SkimConfig.h"
#include "TauSpinerInterface.h"


#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
#include "TVector3.h"
#include "TMath.h"
#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/MultiProngTauSolver.h"
#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"
#include "SimpleFits/FitSoftware/interface/TauA1NuConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/DiTauConstrainedFitter.h"
#include "SimpleFits/FitSoftware/interface/GlobalEventFit.h"
#include "Objects.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
//#include "TauAnalysis/ClassicSVfit/interface/FastMTT.h"
#include "TauPolSoftware/TauDecaysInterface/interface/fonction_a1.h"
#include "TauPolSoftware/TauDecaysInterface/interface/SCalculator.h"



HCPTauTau::HCPTauTau(TString Name_, TString id_):
  Selection(Name_,id_)//,
  //DataMC_Corr(true,true,false),
  //tauTrgSF("tight")
{
  ChargeSumDummy = -999;
  selMuon_IsoDummy = 999.;
  
  WorkSpaceFF2016=TFile::Open(((std::string)std::getenv("workdir")+"Code/fake_factors_cpdecay_v2/fakefactors_ws_tt_lite_2016.root").c_str(), "READ");
  wFF2016= (RooWorkspace*)gDirectory->Get("w");
  WorkSpaceFF2016->Close();
}

HCPTauTau::~HCPTauTau(){
  for(unsigned int j=0; j<Npassed.size(); j++){
    Logger(Logger::Info) << "Selection Summary before: "
			 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
			 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
  delete wFF2016;
}

void  HCPTauTau::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    //if(i==Trigger)             cut.at(Trigger)=1;
    //if(i==Id_and_Kin)          cut.at(Id_and_Kin)=1;
    //if(i==NPairsFound)         cut.at(NPairsFound)=1;
    //if(i==GoodIndex)           cut.at(GoodIndex)=1.;
    if(i==METFiltersAndBTagVeto) cut.at(METFiltersAndBTagVeto)=1.;
    if(i==TausIsolation)       cut.at(TausIsolation)=1.;
    if(i==AgainstEleMu)       cut.at(AgainstEleMu)=1.;
    //if(i==Tau2Isolation)       cut.at(Tau2Isolation)=1.;
    if(i==LeptonVeto)          cut.at(LeptonVeto)=0.;
    if(i==PairCharge)          cut.at(PairCharge)=1.;
    //if(i==PairMass)            cut.at(PairMass)=85.;
    //if(i==MTM)                 cut.at(MTM)=40;
    
  }
  // Setup cut plots
  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
    // if(i==PrimeVtx){
    //   title.at(i)="Number of Prime Vertices $(N>$";
    //   title.at(i)+=cut.at(PrimeVtx);
    //   title.at(i)+=")";
    //   htitle=title.at(i);
    //   htitle.ReplaceAll("$","");
    //   htitle.ReplaceAll("\\","#");
    //   hlabel="Number of Prime Vertices";
    //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PrimeVtx_",htitle,51,-0.5,50.5,hlabel,"Events"));
    //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PrimeVtx_",htitle,51,-0.5,50.5,hlabel,"Events"));
    // }
    // if(i==Trigger){
    //   title.at(i)="Trigger+Matching";
    //   hlabel="At least 1 good pair with Trig+Matching";
    //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Trigger_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Trigger_",htitle,2,-0.5,1.5,hlabel,"Events"));
    // }
    // else if(i==Id_and_Kin){
    //   title.at(i)="Id and Kinematic";
    //   hlabel="Number of Event with good particles";
    //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Id_and_Kin_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Id_and_Kin_",htitle,2,-0.5,1.5,hlabel,"Events"));
    // }
    // else if(i==NPairsFound){
    //   title.at(i)="Pairs with good DeltaR";
    //   hlabel="Pairs with good DeltaR";
    //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NPairsFound_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NPairsFound_",htitle,2,-0.5,1.5,hlabel,"Events"));
    // }
    //if(i==GoodIndex){
      //title.at(i)="Valid Index";
      //hlabel="Valid Index";
      //Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_GoodIndex_",htitle,2,-0.5,1.5,hlabel,"Events"));
      //Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_GoodIndex_",htitle,2,-0.5,1.5,hlabel,"Events"));
      // }
    if(i==METFiltersAndBTagVeto){
      title.at(i)="MET Filters and BTag Veto";
      hlabel="MET Filters and BTag Veto";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_METFiltersAndBTagVeto_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_METFiltersAndBTagVeto_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==TausIsolation){
      title.at(i)="Taus Isolation";
      hlabel="Isolation of Taus";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TausIsolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TausIsolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==AgainstEleMu){
      title.at(i)="Against Electrons and Muons";
      hlabel="Against Electrons and Muons";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_AgainstEleMu_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_AgainstEleMu_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    // else if(i==Tau2Isolation){
    //   title.at(i)="Tau2 Isolation";
    //   hlabel="Isolation of Tau2";
    //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Tau2Isolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Tau2Isolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
    // }
    else if(i==LeptonVeto){
      title.at(i)="Third Lepton Veto";
      hlabel="Third Lepton Veto";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_LeptonVeto_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_LeptonVeto_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==PairCharge){
      title.at(i)="Pair Charge";
      hlabel="is pair OS";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PairCharge_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PairCharge_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
     // else if(i==PairMass){
     //   title.at(i)="Pair Visible Mass";
     //   hlabel="M(tau-tau)";
     //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PairMass_",htitle,30,0,150,hlabel,"Events"));
     //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PairMass_",htitle,30,0,150,hlabel,"Events"));
     // }

    /* else if(i==MTM){
       title.at(i)="Missing Transverse Mass";
       hlabel="Missing Transverse Mass";
       Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MTM_",htitle,30,0,100,hlabel,"Events"));
       Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MTM_",htitle,30,0,100,hlabel,"Events"));
       }*/
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  
  Tau1PT=HConfig.GetTH1D(Name+"_Tau1PT","Transverse momentum of selected #tau1 candidate",20,40,140," P_{T}(#tau1), GeV","Events");
  Tau1E=HConfig.GetTH1D(Name+"_Tau1E","Energy of selected #tau1 candidate",20,35.,140," E(#tau1), GeV","Events");
  Tau1Mass=HConfig.GetTH1D(Name+"_Tau1Mass","Mass of selected #tau1 candidate",20,0,2.," M(#tau1), GeV","Events");
  Tau1Phi=HConfig.GetTH1D(Name+"_Tau1Phi","Phi of selected #tau1 candidate",10,-3.5,3.5," #phi(#tau1)","Events");
  Tau1Eta=HConfig.GetTH1D(Name+"_Tau1Eta","Pseudorapidity tau1",15,-2.7,2.7," #eta(#tau1)","Events");
  Tau1dz=HConfig.GetTH1D(Name+"_Tau1dz","Tau1 dz",20,-0.04,0.04,"Tau1 dz","Events");
  Tau1HPSDecayMode=HConfig.GetTH1D(Name+"_Tau1HPSDecayMode","Decay mode of the selected #tau candidate",12,-0.5,11.5,"HPS Mode","Events");
  Tau1MVADecayMode=HConfig.GetTH1D(Name+"_Tau1MVADecayMode","MVA decay mode of the selected #tau1 candidate",13,-1.5,11.5,"MVA DM","Events");
  
  Tau2PT=HConfig.GetTH1D(Name+"_Tau2PT","Transverse momentum of selected #tau2 candidate",12,40,100," P_{T}(#tau2), GeV","Events");
  Tau2E=HConfig.GetTH1D(Name+"_Tau2E","Energy of selected #tau2 candidate",20,30.,140," E(#tau2), GeV","Events");
  Tau2Mass=HConfig.GetTH1D(Name+"_Tau2Mass","Mass of selected #tau2 candidate",20,0,2.," M(#tau2), GeV","Events");
  Tau2Phi=HConfig.GetTH1D(Name+"_Tau2Phi","Phi of selected #tau2 candidate",10,-3.5,3.5," #phi(#tau2)","Events");
  Tau2Eta=HConfig.GetTH1D(Name+"_Tau2Eta","Pseudorapidity Tau2",15,-2.7,2.7," #eta(#tau2)","Events");
  Tau2dz=HConfig.GetTH1D(Name+"_Tau2dz","Tau2dz",20,-0.04,0.04,"Tau2 dz","Events");
  Tau2HPSDecayMode=HConfig.GetTH1D(Name+"_Tau2HPSDecayMode","Decay mode of the selected #tau candidate",12,-0.5,11.5," HPS Mode ","Events");
  Tau2MVADecayMode=HConfig.GetTH1D(Name+"_Tau2MVADecayMode","MVA decay mode of the selected #tau2 candidate",13,-1.5,11.5,"MVA DM","Events");
  
  // Tau1isolation=HConfig.GetTH1D(Name+"_Tau1isolation","First Tau isolation 0- VVVLoose, 1- VVLoose, 2 VLoose, 3-Loose, 4-Medium, 5-Tight, 6-VTight, 7-VVTight",8,-0.5,7.5,"Discriminator","Events");
  // Tau2isolation=HConfig.GetTH1D(Name+"_Tau2isolation","Second Tau isolation 0- VVVLoose, 1- VVLoose, 2 VLoose, 3-Loose, 4-Medium, 5-Tight, 6-VTight, 7-VVTight",8,-0.5,7.5," Discriminator","Events");
  
    /*
    againstElectronVLooseMVA6_Tau1=HConfig.GetTH1D(Name+"_againstElectronVLooseMVA6_Tau1","againstElectronVLooseMVA6_Tau1",2,-0.5,1.5,"againstElectronVLooseMVA6_Tau1","Events");
    againstElectronLooseMVA6_Tau1=HConfig.GetTH1D(Name+"_againstElectronLooseMVA6_Tau1","againstElectronLooseMVA6_Tau1",2,-0.5,1.5,"againstElectronLooseMVA6_Tau1","Events");
    againstElectronMediumMVA6_Tau1=HConfig.GetTH1D(Name+"_againstElectronMediumMVA6_Tau1","againstElectronMediumMVA6_Tau1",2,-0.5,1.5,"againstElectronMediumMVA6_Tau1","Events");
    againstElectronTightMVA6_Tau1=HConfig.GetTH1D(Name+"_againstElectronTightMVA6_Tau1","againstElectronTightMVA6_Tau1",2,-0.5,1.5,"againstElectronTightMVA6_Tau1","Events");
    againstElectronVTightMVA6_Tau1=HConfig.GetTH1D(Name+"_againstElectronVTightMVA6_Tau1","againstElectronVTightMVA6_Tau1",2,-0.5,1.5,"againstElectronVTightMVA6_Tau1","Events");
    againstMuonLoose3_Tau1=HConfig.GetTH1D(Name+"_againstMuonLoose3_Tau1","againstMuonLoose3_Tau1",2,-0.5,1.5,"againstMuonLoose3_Tau1","Events");
    againstMuonTight3_Tau1=HConfig.GetTH1D(Name+"_againstMuonTight3_Tau1","againstMuonTight3_Tau1",2,-0.5,1.5,"againstMuonTight3_Tau1","Events");
    byCombinedIsolationDeltaBetaCorrRaw3Hits_Tau1=HConfig.GetTH1D(Name+"_byCombinedIsolationDeltaBetaCorrRaw3Hits_Tau1","byCombinedIsolationDeltaBetaCorrRaw3Hits_Tau1",10,0,20,"byCombinedIsolationDeltaBetaCorrRaw3Hits_Tau1","Events");

    againstElectronVLooseMVA6_Tau2=HConfig.GetTH1D(Name+"_againstElectronVLooseMVA6_Tau2","againstElectronVLooseMVA6_Tau2",2,-0.5,1.5,"againstElectronVLooseMVA6_Tau2","Events");
    againstElectronLooseMVA6_Tau2=HConfig.GetTH1D(Name+"_againstElectronLooseMVA6_Tau2","againstElectronLooseMVA6_Tau2",2,-0.5,1.5,"againstElectronLooseMVA6_Tau2","Events");
    againstElectronMediumMVA6_Tau2=HConfig.GetTH1D(Name+"_againstElectronMediumMVA6_Tau2","againstElectronMediumMVA6_Tau2",2,-0.5,1.5,"againstElectronMediumMVA6_Tau2","Events");
    againstElectronTightMVA6_Tau2=HConfig.GetTH1D(Name+"_againstElectronTightMVA6_Tau2","againstElectronTightMVA6_Tau2",2,-0.5,1.5,"againstElectronTightMVA6_Tau2","Events");
    againstElectronVTightMVA6_Tau2=HConfig.GetTH1D(Name+"_againstElectronVTightMVA6_Tau2","againstElectronVTightMVA6_Tau2",2,-0.5,1.5,"againstElectronVTightMVA6_Tau2","Events");
    againstMuonLoose3_Tau2=HConfig.GetTH1D(Name+"_againstMuonLoose3_Tau2","againstMuonLoose3_Tau2",2,-0.5,1.5,"againstMuonLoose3_Tau2","Events");
    againstMuonTight3_Tau2=HConfig.GetTH1D(Name+"_againstMuonTight3_Tau2","againstMuonTight3_Tau2",2,-0.5,1.5,"againstMuonTight3_Tau2","Events");
    byCombinedIsolationDeltaBetaCorrRaw3Hits_Tau2=HConfig.GetTH1D(Name+"_byCombinedIsolationDeltaBetaCorrRaw3Hits_Tau2","byCombinedIsolationDeltaBetaCorrRaw3Hits_Tau2",10,0,20,"byCombinedIsolationDeltaBetaCorrRaw3Hits_Tau2","Events");
  */
  ExtraLeptonVeto=HConfig.GetTH1D(Name+"_ExtraLeptonVeto","ExtraLeptonVeto",2,-0.5,1.5,"ExtraLeptonVeto","Events");
  
  QCDShape=HConfig.GetTH1D(Name+"_QCDShape","QCDShape",2,-0.5,1.5,"QCD Shape","");

  TauTauVisMass=HConfig.GetTH1D(Name+"_TauTauVisMass","Visible invariant mass of a tau pair",23,20,250," M(#tau#tau)_{vis}, GeV","Events");
  //  TauTauTruthMass=HConfig.GetTH1D(Name+"_TauTauTruthMass","Truth invariant mass of a tau pair",40,0,150," M(#tau#tau)_{truth}, GeV","Events");
  //TauTauFullMass=HConfig.GetTH1D(Name+"_TauTauFullMass","Full invariant mass of a tau pair",40,0,150," M(#tau#tau)_{full}, GeV","Events");

  NQCD=HConfig.GetTH1D(Name+"_NQCD","NQCD",4,0.5,4.5,"NQCD in ABCD","Events");
  // TauTauFullMass_B=HConfig.GetTH1D(Name+"_TauTauFullMass_B","TauTauFullMass_B",40,0,150,"#tau_h#tau_h SVFit Mass in B","Events");
  // TauTauFullMass_C=HConfig.GetTH1D(Name+"_TauTauFullMass_C","TauTauFullMass_C",40,0,150,"#tau_h#tau_h SVFit Mass in C","Events");
  // TauTauFullMass_D=HConfig.GetTH1D(Name+"_TauTauFullMass_D","TauTauFullMass_D",40,0,150,"#tau_h#tau_h SVFit Mass in D","Events");

  NFFData=HConfig.GetTH1D(Name+"_NFFData","NFFData",1,0.5,1.5,"NFFData","Events");
  NFFLeadMC=HConfig.GetTH1D(Name+"_NFFLeadMC","NFFLeadMC",1,0.5,1.5,"NFFLeadMC","Events");
  W_res=HConfig.GetTH1D(Name+"_W_res","W_res",1,0.5,1.5,"W_res","Events");
  dRTauTau=HConfig.GetTH1D(Name+"_dRTauTau","#Delta R",20,0.,3.5," #Delta R","Events");

  MET=HConfig.GetTH1D(Name+"_MET","MET",20,0,80,"MET, GeV","Events");
  METphi=HConfig.GetTH1D(Name+"_METphi","METphi",10,-3.5,3.5,"METphi","Events");
  PUPPImet=HConfig.GetTH1D(Name+"_PUPPImet","PUPPImet",10,0,75,"PUPPImet, GeV","Events");
  PUPPImetphi=HConfig.GetTH1D(Name+"_PUPPImetphi","PUPPImetphi",10,-3.5,3.5,"PUPPImetphi","Events");
  PUPPImetcorr=HConfig.GetTH1D(Name+"_PUPPImetcorr","PUPPImetcorr",20,0,200,"PUPPImetcorr, GeV","Events");
  PUPPImetcorrphi=HConfig.GetTH1D(Name+"_PUPPImetcorrphi","PUPPImetcorrphi",10,-3.5,3.5,"PUPPImetcorrphi","Events");
  TransverseMass=HConfig.GetTH1D(Name+"_TransverseMass","TransverseMass, GeV",40,0,110,"TransverseMass","Events");
  
  NPrimeVtx=HConfig.GetTH1D(Name+"_NPrimeVtx","NPrimeVtx",10,0,50,"N vtx","Events");
  NPU=HConfig.GetTH1D(Name+"_npu","npu",10,0,50,"N pu","Events");
  RHO=HConfig.GetTH1D(Name+"_rho","rho",10,0,30,"rho","Events");
  
  NbJets=HConfig.GetTH1D(Name+"_NbJets","NbJets",5,-0.5,4.5,"N_{jets}","Events");
																			    
  ZPtVis=HConfig.GetTH1D(Name+"_ZPtVis","Visible Pt_{Z}",40,0,100,"","Events");

  TauTauVisPT=HConfig.GetTH1D(Name+"_TauTauVisPT","Visible Pt_{#tau#tau}",40,0,100,"","Events");
  Mjj=HConfig.GetTH1D(Name+"_Mjj","m_{jj}",15,0,400,"m_{jj} GeV","Events");

  ResolPullTauTauFroma1a1MZMomentum=HConfig.GetTH1D(Name+"_ResolPullTauTauFroma1a1MZMomentum","ResolPullTauTauFroma1a1MZMomentum",30,-1,1,"","Events");


  ResolPullTauminusFroma1a1MZMomentum=HConfig.GetTH1D(Name+"_ResolPullTauminusFroma1a1MZMomentum","ResolPullTauminusFroma1a1MZMomentum",30,-1,1,"","Events");


  ResolPullTauplusFroma1a1MZMomentum=HConfig.GetTH1D(Name+"_ResolPullTauplusFroma1a1MZMomentum","ResolPullTauplusFroma1a1MZMomentum",30,-1,1,"","Events");



  ResolPullTauFroma1a1MZMomentum=HConfig.GetTH1D(Name+"_ResolPullTauFroma1a1MZMomentum","ResolPullTauFroma1a1MZMomentum",30,-1,1,"","Events");


  ResolPullXVtxIna1a1=HConfig.GetTH1D(Name+"_ResolPullXVtxIna1a1","ResolPullXVtxIna1a1",30,-1,1,"","Events");
  ResolPullYVtxIna1a1=HConfig.GetTH1D(Name+"_ResolPullYVtxIna1a1","ResolPullYVtxIna1a1",30,-1,1,"","Events");
  ResolPullZVtxIna1a1=HConfig.GetTH1D(Name+"_ResolPullZVtxIna1a1","ResolPullZVtxIna1a1",30,-1,1,"","Events");

  tauminusa1a1MomentumPairConstraint=HConfig.GetTH1D(Name+"_tauminusa1a1MomentumPairConstraint","tauminusa1a1MomentumPairConstraint",15,0,200,"","Events");       
  tauplusa1a1MomentumPairConstraint=HConfig.GetTH1D(Name+"_tauplusa1a1MomentumPairConstraint","tauplusa1a1MomentumPairConstraint",15,0,200,"","Events");

  polarimetricAcopAngle=HConfig.GetTH1D(Name+"_polarimetricAcopAngle","GEF",60,0.,2*TMath::Pi(),"GEF","Events");
  
  

  polarimetricAcopAnglePVRefitNoBS=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBS","GEF PVRefit no BS constraint new",60,0.,2*TMath::Pi(),"GEF PVRefit no BS constraint new","Events");
  polarimetricAcopAnglePVRefitBS=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBS","GEF PVRefit BS constraint new",60,0.,2*TMath::Pi(),"GEF PVRefit BS constraint new","Events");
  polarimetricAcopAnglePVRefitNoBSZNominal=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBSZNominal","GEF PVRefit no BS constraint nominal Z new",60,0.,2*TMath::Pi(),"GEF PVRefit no BS constraint nominal Z new","Events");
  polarimetricAcopAnglePVRefitBSZNominal=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSZNominal","GEF PVRefit BS constraint nominal Z new",60,0.,2*TMath::Pi(),"GEF PVRefit BS constraint nominal Z new","Events");

  polarimetricAcopAngleMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAngleMVA","GEF MVA",60,0.,2*TMath::Pi(),"GEF MVA","Events");
  polarimetricAcopAnglePVRefitNoBSNewMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBSNewMVA","GEF PVRefit no BS constraint new MVA",60,0.,2*TMath::Pi(),"GEF PVRefit no BS constraint new MVA","Events");
  polarimetricAcopAnglePVRefitBSNewMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSNewMVA","GEF PVRefit BS constraint new MVA",60,0.,2*TMath::Pi(),"GEF PVRefit BS constraint new MVA","Events");
  polarimetricAcopAnglePVRefitNoBSZNominalNewMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBSZNominalNewMVA","GEF PVRefit no BS constraint nominal Z new MVA",60,0.,2*TMath::Pi(),"GEF PVRefit no BS constraint nominal Z new MVA","Events");
  polarimetricAcopAnglePVRefitBSZNominalNewMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSZNominalNewMVA","GEF PVRefit BS constraint nominal Z new MVA",60,0.,2*TMath::Pi(),"GEF PVRefit BS constraint nominal Z new MVA","Events");


  polarimetricAcopAnglePVRefitOnlyNoBS=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitOnlyNoBS","GEF PVRefit Only no BS constraint",60,0.,2*TMath::Pi(),"GEF PVRefit Only no BS constraint","Events");
  polarimetricAcopAnglePVRefitOnlyBS=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitOnlyBS","GEF PVRefit Only BS constraint",60,0.,2*TMath::Pi(),"GEF PVRefit Only BS constraint","Events");
  polarimetricAcopAnglePVRefitOnlyNoBSZNominal=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitOnlyNoBSZNominal","GEF PVRefit Only no BS constraint nominal Z",60,0.,2*TMath::Pi(),"GEF PVRefit Only no BS constraint nominal Z","Events");
  polarimetricAcopAnglePVRefitOnlyBSZNominal=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitOnlyBSZNominal","GEF PVRefit Only BS constraint nominal Z",60,0.,2*TMath::Pi(),"GEF PVRefit Only BS constraint nominal Z","Events");


  polarimetricAcopAnglePVRefitOnlyNoBSNewMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVOnlyRefitNoBSNewMVA","GEF PVRefit Only no BS constraint MVA",60,0.,2*TMath::Pi(),"GEF PVRefit Only no BS constraint MVA","Events");
  polarimetricAcopAnglePVRefitOnlyBSNewMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVOnlyRefitBSNewMVA","GEF PVRefit Only BS constraint MVA",60,0.,2*TMath::Pi(),"GEF PVRefit Only BS constraint MVA","Events");
  polarimetricAcopAnglePVRefitOnlyNoBSZNominalNewMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVOnlyRefitNoBSZNominalNewMVA","GEF PVRefit Only no BS constraint nominal Z MVA",60,0.,2*TMath::Pi(),"GEF PVRefit Only no BS constraint nominal Z MVA","Events");
  polarimetricAcopAnglePVRefitOnlyBSZNominalNewMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVOnlyRefitBSZNominalNewMVA","GEF PVRefit Only BS constraint nominal Z MVA",60,0.,2*TMath::Pi(),"GEF PVRefit Only BS constraint nominal Z MVA","Events");

  

  PVXResol=HConfig.GetTH1D(Name+"_PVXResol","PV_{X} pull",50,-0.1,0.1,"PV_{X} pull","Events");
 
  PVXNoBSResol=HConfig.GetTH1D(Name+"_PVXNoBSResol","No BS PV_{X}Refit pull",50,-0.1,0.1,"No BS PV_{X}Refit pull","Events");
  
  PVXBSResol=HConfig.GetTH1D(Name+"_PVXBSResol","BS PV_{X}Refit pull",50,-0.1,0.1,"BS PV_{X}Refit pull","Events");
  
  PVYResol=HConfig.GetTH1D(Name+"_PVYResol","PV_{Y} pull",50,-0.1,0.1,"PV_{Y} pull","Events");
  
  PVYNoBSResol=HConfig.GetTH1D(Name+"_PVYNoBSResol","No BS PV_{Y}Refit pull",50,-0.1,0.1,"No BS PV_{Y}Refit pull","Events");
  
  PVYBSResol=HConfig.GetTH1D(Name+"_PVYBSResol","BS PV_{Y}Refit pull",50,-0.1,0.1,"BS PV_{Y}Refit pull","Events");

  PVZResol=HConfig.GetTH1D(Name+"_PVZResol","PV_{Z} pull",50,-0.01,0.01,"PV_{Z} pull","Events");

  PVZNoBSResol=HConfig.GetTH1D(Name+"_PVZNoBSResol","No BS PV_{Z}Refit pull",50,-0.01,0.01,"No BS PV_{Z}Refit pull","Events");

  PVZBSResol=HConfig.GetTH1D(Name+"_PVZBSResol","BS PV_{Z}Refit pull",50,-0.01,0.01,"BS PV_{Z}Refit pull","Events");
  

  PVXNoBSOnlyResol=HConfig.GetTH1D(Name+"_PVXNoBSOnlyResol","No BS PV_{X}Refit Only pull",50,-0.1,0.1,"No BS PV_{X}Refit Only pull","Events");
 
  PVXBSOnlyResol=HConfig.GetTH1D(Name+"_PVXBSOnlyResol","BS PV_{X}Refit Only pull",50,-0.1,0.1,"BS PV_{X}Refit Only pull","Events");

  PVYNoBSOnlyResol=HConfig.GetTH1D(Name+"_PVYNoBSOnlyResol","No BS PV_{Y}Refit Only pull",50,-0.1,0.1,"No BS PV_{Y}Refit Only pull","Events");
  PVYBSOnlyResol=HConfig.GetTH1D(Name+"_PVYBSOnlyResol","BS PV_{Y}Refit Only pull",50,-0.1,0.1,"BS PV_{Y}Refit Only pull","Events");

  PVZNoBSOnlyResol=HConfig.GetTH1D(Name+"_PVZNoBSOnlyResol","No BS PV_{Z}Refit Only pull",50,-0.01,0.01,"No BS PV_{Z}Refit Only pull","Events");

  PVZBSOnlyResol=HConfig.GetTH1D(Name+"_PVZBSOnlyResol","BS PV_{Z}Refit Only pull",50,-0.01,0.01,"BS PV_{Z}Refit Only pull","Events");
  
  polarimetricAcopAngleTruthA1=HConfig.GetTH1D(Name+"_polarimetricAcopAngleTruthA1","polarimetricAcopAngleTruthA1",60,0,2*TMath::Pi(),"","Events");

  polarimetricAcopAngleDecayPlane=HConfig.GetTH1D(Name+"_polarimetricAcopAngleDecayPlane","polarimetricAcopAngleDecayPlane",60,0.,2*TMath::Pi(),"","Events");
  
  test=HConfig.GetTH1D(Name+"_test","test",60,0.,2*TMath::Pi(),"","Events"); 

  PurityDM=HConfig.GetTH1D(Name+"_PurityDM","PurityDM",29,0.,29,"","Events");
  PurityNewMVA=HConfig.GetTH1D(Name+"_PurityNewMVA","PurityNewMVA",29,0.,29,"","Events");

  TauPxResPull=HConfig.GetTH1D(Name+"_TauPxResPull","TauPxResPull",50,-2,2,"P_{X} pull of #tau","Events");
  TauPyResPull=HConfig.GetTH1D(Name+"_TauPyResPull","TauPyResPull",50,-2,2,"P_{Y} pull of #tau","Events");
  TauPzResPull=HConfig.GetTH1D(Name+"_TauPzResPull","TauPzResPull",50,-2,2,"P_{Z} pull of #tau","Events");

  TauPxResPullMVA=HConfig.GetTH1D(Name+"_TauPxResPullMVA","TauPxResPullMVA",50,-2,2,"P_{X} pull of #tau","Events");
  TauPyResPullMVA=HConfig.GetTH1D(Name+"_TauPyResPullMVA","TauPyResPullMVA",50,-2,2,"P_{Y} pull of #tau","Events");
  TauPzResPullMVA=HConfig.GetTH1D(Name+"_TauPzResPullMVA","TauPzResPullMVA",50,-2,2,"P_{Z} pull of #tau","Events");

  PUPPImetcorrQCDMC=HConfig.GetTH1D(Name+"_PUPPImetcorrQCDMC","PUPPImetcorrQCDMC",20,0,200,"PUPPImetcorr, GeV","Events");

  Tau1PTQCDMC=HConfig.GetTH1D(Name+"_Tau1PTQCDMC","Transverse momentum of selected #tau1 candidate",20,40,140," P_{T}(#tau1), GeV","Events");
  Tau1EQCDMC=HConfig.GetTH1D(Name+"_Tau1EQCDMC","Energy of selected #tau1 candidate",20,35.,140," E(#tau1), GeV","Events");
  Tau1MassQCDMC=HConfig.GetTH1D(Name+"_Tau1MassQCDMC","Mass of selected #tau1 candidate",20,0,2.," M(#tau1), GeV","Events");
  Tau1PhiQCDMC=HConfig.GetTH1D(Name+"_Tau1PhiQCDMC","Phi of selected #tau1 candidate",10,-3.5,3.5," #phi(#tau1)","Events");
  Tau1EtaQCDMC=HConfig.GetTH1D(Name+"_Tau1EtaQCDMC","Pseudorapidity tau1",15,-2.7,2.7," #eta(#tau1)","Events");
  Tau1dzQCDMC=HConfig.GetTH1D(Name+"_Tau1dzQCDMC","Tau1 dz",20,-0.04,0.04,"Tau1 dz","Events");
  Tau1HPSDecayModeQCDMC=HConfig.GetTH1D(Name+"_Tau1HPSDecayModeQCDMC","Decay mode of the selected #tau candidate",12,-0.5,11.5,"HPS Mode","Events");
  Tau1MVADecayModeQCDMC=HConfig.GetTH1D(Name+"_Tau1MVADecayModeQCDMC","MVA decay mode of the selected #tau1 candidate",13,-1.5,11.5,"MVA DM","Events");
  
  Tau2PTQCDMC=HConfig.GetTH1D(Name+"_Tau2PTQCDMC","Transverse momentum of selected #tau2 candidate",12,40,100," P_{T}(#tau2), GeV","Events");
  Tau2EQCDMC=HConfig.GetTH1D(Name+"_Tau2EQCDMC","Energy of selected #tau2 candidate",20,30.,140," E(#tau2), GeV","Events");
  Tau2MassQCDMC=HConfig.GetTH1D(Name+"_Tau2MassQCDMC","Mass of selected #tau2 candidate",20,0,2.," M(#tau2), GeV","Events");
  Tau2PhiQCDMC=HConfig.GetTH1D(Name+"_Tau2PhiQCDMC","Phi of selected #tau2 candidate",10,-3.5,3.5," #phi(#tau2)","Events");
  Tau2EtaQCDMC=HConfig.GetTH1D(Name+"_Tau2EtaQCDMC","Pseudorapidity Tau2",15,-2.7,2.7," #eta(#tau2)","Events");
  Tau2dzQCDMC=HConfig.GetTH1D(Name+"_Tau2dzQCDMC","Tau2dz",20,-0.04,0.04,"Tau2 dz","Events");
  Tau2HPSDecayModeQCDMC=HConfig.GetTH1D(Name+"_Tau2HPSDecayModeQCDMC","Decay mode of the selected #tau candidate",12,-0.5,11.5," HPS Mode ","Events");
  Tau2MVADecayModeQCDMC=HConfig.GetTH1D(Name+"_Tau2MVADecayModeQCDMC","MVA decay mode of the selected #tau2 candidate",13,-1.5,11.5,"MVA DM","Events");
  
  NbJetsQCDMC=HConfig.GetTH1D(Name+"_NbJetsQCDMC","NbJetsQCDMC",5,-0.5,4.5,"N_{jets}","Events");
  TauTauVisMassQCDMC=HConfig.GetTH1D(Name+"_TauTauVisMassQCDMC","Visible invariant mass of a tau pair",23,20,250," M(#tau#tau)_{vis}, GeV","Events");
  TauTauVisPTQCDMC=HConfig.GetTH1D(Name+"_TauTauVisPTQCDMC","Visible Pt_{#tau#tau}",40,0,100,"","Events");
  MjjQCDMC=HConfig.GetTH1D(Name+"_MjjQCDMC","m_{jj}",15,0,400,"m_{jj} GeV","Events");
  
  Selection::ConfigureHistograms();   //   do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);  // do not remove
  
}

void  HCPTauTau::Store_ExtraDist(){

  //every new histo should be addedd to Extradist1d vector, just push it back;
  Extradist1d.push_back(&Tau1PT);
  Extradist1d.push_back(&Tau1E);
  Extradist1d.push_back(&Tau1Mass);
  Extradist1d.push_back(&Tau1Phi);
  Extradist1d.push_back(&Tau1Eta);
  Extradist1d.push_back(&Tau1dz);
  Extradist1d.push_back(&Tau1HPSDecayMode);
  Extradist1d.push_back(&Tau1MVADecayMode);
  Extradist1d.push_back(&Tau2PT);
  Extradist1d.push_back(&Tau2E);
  Extradist1d.push_back(&Tau2Mass);
  Extradist1d.push_back(&Tau2Phi);
  Extradist1d.push_back(&Tau2Eta);
  Extradist1d.push_back(&Tau2dz);
  Extradist1d.push_back(&Tau2HPSDecayMode);
  Extradist1d.push_back(&Tau2MVADecayMode);
  // Extradist1d.push_back(&Tau1isolation);
  // Extradist1d.push_back(&Tau2isolation);
    /*
    Extradist1d.push_back(&againstElectronVLooseMVA6_Tau1);
    Extradist1d.push_back(&againstElectronLooseMVA6_Tau1);
    Extradist1d.push_back(&againstElectronMediumMVA6_Tau1);
    Extradist1d.push_back(&againstElectronTightMVA6_Tau1);
    Extradist1d.push_back(&againstElectronVTightMVA6_Tau1);
    Extradist1d.push_back(&againstMuonLoose3_Tau1);
    Extradist1d.push_back(&againstMuonTight3_Tau1);
    Extradist1d.push_back(&byCombinedIsolationDeltaBetaCorrRaw3Hits_Tau1);

    Extradist1d.push_back(&againstElectronVLooseMVA6_Tau2);
    Extradist1d.push_back(&againstElectronLooseMVA6_Tau2);
    Extradist1d.push_back(&againstElectronMediumMVA6_Tau2);
    Extradist1d.push_back(&againstElectronTightMVA6_Tau2);
    Extradist1d.push_back(&againstElectronVTightMVA6_Tau2);
    Extradist1d.push_back(&againstMuonLoose3_Tau2);
    Extradist1d.push_back(&againstMuonTight3_Tau2);
    Extradist1d.push_back(&byCombinedIsolationDeltaBetaCorrRaw3Hits_Tau2);
  */
  Extradist1d.push_back(&ExtraLeptonVeto);


  Extradist1d.push_back(&dRTauTau);
  Extradist1d.push_back(&TauTauVisMass);
  //  Extradist1d.push_back(&TauTauTruthMass);
  //Extradist1d.push_back(&TauTauFullMass);
  
  Extradist1d.push_back(&QCDShape);
  Extradist1d.push_back(&NQCD);
  // Extradist1d.push_back(&TauTauFullMass_B);
  // Extradist1d.push_back(&TauTauFullMass_C);
  // Extradist1d.push_back(&TauTauFullMass_D);

  Extradist1d.push_back(&NFFData);
  Extradist1d.push_back(&NFFLeadMC);
  Extradist1d.push_back(&W_res);

  Extradist1d.push_back(&MET);
  Extradist1d.push_back(&METphi);
  Extradist1d.push_back(&PUPPImet);
  Extradist1d.push_back(&PUPPImetphi);
  Extradist1d.push_back(&PUPPImetcorr);
  Extradist1d.push_back(&PUPPImetcorrphi);
  Extradist1d.push_back(&TransverseMass);

  Extradist1d.push_back(&NPrimeVtx);
  Extradist1d.push_back(&NPU);
  Extradist1d.push_back(&RHO);


  Extradist1d.push_back(&NbJets);

  Extradist1d.push_back(&ZPtVis);
  Extradist1d.push_back(&TauTauVisPT);
  Extradist1d.push_back(&Mjj);

  Extradist1d.push_back(&ResolPullTauTauFroma1a1MZMomentum);
  Extradist1d.push_back(&ResolPullTauminusFroma1a1MZMomentum);

  Extradist1d.push_back(&ResolPullTauplusFroma1a1MZMomentum);
  Extradist1d.push_back(&ResolPullTauFroma1a1MZMomentum);

  Extradist1d.push_back(&ResolPullXVtxIna1a1);
  Extradist1d.push_back(&ResolPullYVtxIna1a1);
  Extradist1d.push_back(&ResolPullZVtxIna1a1);

  Extradist1d.push_back(&tauminusa1a1MomentumPairConstraint);
  Extradist1d.push_back(&tauplusa1a1MomentumPairConstraint);
  
  Extradist1d.push_back(&polarimetricAcopAngle);
 
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBS);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitBS);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBSZNominal);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSZNominal);

  Extradist1d.push_back(&polarimetricAcopAnglePVRefitOnlyNoBS);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitOnlyBS);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitOnlyNoBSZNominal);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitOnlyBSZNominal);

  Extradist1d.push_back(&polarimetricAcopAngleMVA);

  Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBSNewMVA);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSNewMVA);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBSZNominalNewMVA);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSZNominalNewMVA);

  Extradist1d.push_back(&polarimetricAcopAnglePVRefitOnlyNoBSNewMVA);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitOnlyBSNewMVA);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitOnlyNoBSZNominalNewMVA);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitOnlyBSZNominalNewMVA);

  
  Extradist1d.push_back(&PVXResol);
  Extradist1d.push_back(&PVXNoBSResol);
  Extradist1d.push_back(&PVXBSResol);

  Extradist1d.push_back(&PVYResol);
  Extradist1d.push_back(&PVYNoBSResol);
  Extradist1d.push_back(&PVYBSResol);
  
  Extradist1d.push_back(&PVZResol);
  Extradist1d.push_back(&PVZNoBSResol);
  Extradist1d.push_back(&PVZBSResol); 
 
  Extradist1d.push_back(&PVXNoBSOnlyResol);
  Extradist1d.push_back(&PVXBSOnlyResol);
  
  Extradist1d.push_back(&PVYNoBSOnlyResol);
  Extradist1d.push_back(&PVYBSOnlyResol);
  
  Extradist1d.push_back(&PVZNoBSOnlyResol);
  Extradist1d.push_back(&PVZBSOnlyResol);

  Extradist1d.push_back(&polarimetricAcopAngleTruthA1);
  
  Extradist1d.push_back(&test);  
  
  Extradist1d.push_back(&PurityDM);
  Extradist1d.push_back(&PurityNewMVA);

  Extradist1d.push_back(&TauPxResPull);
  Extradist1d.push_back(&TauPyResPull);
  Extradist1d.push_back(&TauPzResPull);

  Extradist1d.push_back(&TauPxResPullMVA);
  Extradist1d.push_back(&TauPyResPullMVA);
  Extradist1d.push_back(&TauPzResPullMVA);
  
}

void  HCPTauTau::doEvent()  { //  Method called on every event
  unsigned int t;                // sample type, you may manage in your further analysis, if needed
  int id(Ntp->GetMCID());  //read event ID of a sample
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}  //  gives a warning if list of samples in Histo.txt  and SkimSummary.log do not coincide 
  //  std::cout<<"------------------ New Event -----------------------"<<std::endl;
  Charge = ChargeSumDummy;

  std::vector<int> Tau1IndexVect;
  std::vector<int> Tau2IndexVect;
  std::vector<int> IndexSelected;
  bool METFilters=true;
  if(Ntp->year()==2017 || Ntp->year()==2018){
    if(!Ntp->passecalBadCalibFilterUpdate()){
      METFilters=false;
    }
  }
  if(Ntp->metfilterbit()!=255)METFilters=false; // 255 bits because 8 filters but will be 127 when new production launched !!!
  bool passbtagveto=true;
  int nbtagloose=0;
  for(int ijet=0; ijet< Ntp->NJets(); ijet++) {
    if(Ntp->Jet_P4(ijet).Pt()>25 && Ntp->Jet_P4(ijet).Eta()<2.4 && Ntp->year()==2016 && (Ntp->DeepCSV_probb(ijet)+Ntp->DeepCSV_probbb(ijet))> 0.6321)passbtagveto=false;
    if(Ntp->Jet_P4(ijet).Pt()>25 && Ntp->Jet_P4(ijet).Eta()<2.4 && Ntp->year()==2017 && (Ntp->DeepCSV_probb(ijet)+Ntp->DeepCSV_probbb(ijet))> 0.4941)passbtagveto=false;
    if(Ntp->Jet_P4(ijet).Pt()>25 && Ntp->Jet_P4(ijet).Eta()<2.4 && Ntp->year()==2018 && (Ntp->DeepCSV_probb(ijet)+Ntp->DeepCSV_probbb(ijet))> 0.4184)passbtagveto=false;
    
    if(Ntp->Jet_P4(ijet).Pt()>25 && Ntp->Jet_P4(ijet).Eta()<2.4 && Ntp->year()==2016 && (Ntp->DeepCSV_probb(ijet)+Ntp->DeepCSV_probbb(ijet))> 0.2217)nbtagloose++;
    if(Ntp->Jet_P4(ijet).Pt()>25 && Ntp->Jet_P4(ijet).Eta()<2.4 && Ntp->year()==2017 && (Ntp->DeepCSV_probb(ijet)+Ntp->DeepCSV_probbb(ijet))> 0.1522)nbtagloose++;
    if(Ntp->Jet_P4(ijet).Pt()>25 && Ntp->Jet_P4(ijet).Eta()<2.4 && Ntp->year()==2018 && (Ntp->DeepCSV_probb(ijet)+Ntp->DeepCSV_probbb(ijet))> 0.1241)nbtagloose++;
  }
  if(nbtagloose>1 || !passbtagveto)passbtagveto=false;
  
  value.at(METFiltersAndBTagVeto)=0;
  value.at(METFiltersAndBTagVeto)=(METFilters && passbtagveto);
  pass.at(METFiltersAndBTagVeto)=value.at(METFiltersAndBTagVeto);

  for(int iDaughter=0;   iDaughter  <  Ntp->SelectedPairs() ;iDaughter++ ) {
    if(Ntp->trg_doubletau(iDaughter))
      {
  	Tau1IndexVect.push_back(Ntp->tau1IndexVect(iDaughter));
  	Tau2IndexVect.push_back(Ntp->tau2IndexVect(iDaughter));
  	IndexSelected.push_back(iDaughter);
      }
  }
  std::vector<int> goodTaus1Index;
  std::vector<int> goodTaus2Index;
  std::vector<int> IndexSelectedTemp;
  std::vector<int>  PairsIndexTemp;
  int j=0;
  for(unsigned int iDaughter=0;   iDaughter  <IndexSelected.size() ;iDaughter++ ) {
    if(Ntp->tauBaselineSelection(Tau1IndexVect.at(iDaughter),40., 2.1, 1,0,4) && Ntp->tauBaselineSelection(Tau2IndexVect.at(iDaughter),40., 2.1, 1,0,4))
      {
  	goodTaus1Index.push_back(Tau1IndexVect.at(iDaughter));
  	goodTaus2Index.push_back(Tau2IndexVect.at(iDaughter));
  	IndexSelectedTemp.push_back(IndexSelected.at(iDaughter));
  	PairsIndexTemp.push_back(j);
  	j++;
      }
  }

  int Tau1= -1;
  int Tau2= -1;
  std::vector<int>  Sorted;
  if(PairsIndexTemp.size()>0)
    {
      //Sorted = Ntp->SortPair(PairsIndexTemp,goodTaus1Index,goodTaus2Index);
      //Tau1=goodTaus1Index.at(Sorted.back());
      //Tau2=goodTaus2Index.at(Sorted.back());
      Tau1=goodTaus1Index.at(0);
      Tau2=goodTaus2Index.at(0);
      //value.at(GoodIndex)=0;
      //value.at(GoodIndex)=(Tau1!=-99 && Tau2!=-99);
      //pass.at(GoodIndex) = value.at(GoodIndex);
      value.at(TausIsolation)=0;
      value.at(TausIsolation) = (Ntp->byMediumDeepTau2017v2p1VSjet_1(IndexSelectedTemp.at(0)) && Ntp->byMediumDeepTau2017v2p1VSjet_2(IndexSelectedTemp.at(0)));
      pass.at(TausIsolation) = value.at(TausIsolation);
      value.at(AgainstEleMu)=0;
      value.at(AgainstEleMu) = (Ntp->byVVLooseDeepTau2017v2p1VSe_1(IndexSelectedTemp.at(0)) && Ntp->byVVLooseDeepTau2017v2p1VSe_2(IndexSelectedTemp.at(0)) && Ntp->byVLooseDeepTau2017v2p1VSmu_1(IndexSelectedTemp.at(0)) && Ntp->byVLooseDeepTau2017v2p1VSmu_2(IndexSelectedTemp.at(0)) );
      pass.at(AgainstEleMu) = value.at(AgainstEleMu);
      value.at(LeptonVeto)=0;
      int thirdlepton=0;
      for(int iDaughter=0;   iDaughter  <  Ntp->NDaughters() ;iDaughter++ ) {
	if((iDaughter!=Tau1)&&(iDaughter!=Tau2)){
	  if(Ntp->ElectronVeto(iDaughter) || Ntp->MuonVeto(iDaughter))thirdlepton++;
	}
      }
      //if(Ntp->eleveto() || Ntp->muonveto())thirdlepton++;
      value.at(LeptonVeto) = thirdlepton>0;
      pass.at(LeptonVeto) = (value.at(LeptonVeto)==cut.at(LeptonVeto));
      value.at(PairCharge)=0;
      bool isOS=false;
      isOS=((Ntp->Daughters_charge(Tau1)/abs(Ntp->Daughters_charge(Tau1))) != (Ntp->Daughters_charge(Tau2)/abs(Ntp->Daughters_charge(Tau2))));
      if(isOS)value.at(PairCharge) = 1;
      pass.at(PairCharge) = value.at(PairCharge);
       // value.at(PairMass) = 999.;
      // // //value.at(MTM) = 999.;
      // // //value.at(MTM) = .;
      //  value.at(PairMass)=(Ntp->P4Corrected(Tau1) + Ntp->P4Corrected(Tau2)).M();
      //  pass.at(PairMass) = ((value.at(PairMass) > 85.) && (value.at(PairMass) < 115.));
      // // //pass.at(MTM) = (value.at(MTM) <= cut.at(MTM));
    }

  // Here you can defined different type of weights you want to apply to events.
  double wobs=1;
  double w=1;
  if(!Ntp->isData() && id!=DataMCType::QCD) {
    //    w *= reweight.weight(2016,26,Ntp->PUNumInteractions());
    w*=Ntp->PUReweight();
    // if(Ntp->year()==2017) w*=Ntp->PUReweight(Ntp->year(),filePUdistribution2017_data,filePUdistribution2017_MC);
    // if(Ntp->year()==2018) w*=Ntp->PUReweight(Ntp->year(),filePUdistribution2018_data,filePUdistribution2018_MC);
    //w *= reweight.PUweightHTT(Ntp->npu());
    //std::cout<<" pu weigh HTT  "<< reweight.PUweightHTT(Ntp->npu())<<std::endl;
    if(!Ntp->isData() && PairsIndexTemp.size()>0/* pass.at(NPairsFound)*/ ){
      double w1=1.,w2=1.;

      if(Ntp->year()==2016)
	{
	  w1=Ntp->TriggerSF(Tau1);
	  w2=Ntp->TriggerSF(Tau2);
	}
     
      // else if(Ntp->year()==2017)
      // 	{
      // 	  w1=Ntp->TriggerSF(Tau1,w2017_);
      // 	  w2=Ntp->TriggerSF(Tau2,w2017_);
      // 	}
      // else if(Ntp->year()==2018)
      // 	{
      // 	  w1=Ntp->TriggerSF(Tau1,w2018_);
      // 	  w2=Ntp->TriggerSF(Tau2,w2018_);
      // 	}

      //double w1 = tauTrgSF.getSF(Ntp->P4Corrected(Tau1).Pt(),  Ntp->decayMode(Tau1)) ;  //from Luca
      //double w2 = tauTrgSF.getSF(Ntp->P4Corrected(Tau2).Pt(),  Ntp->decayMode(Tau2)) ;
      w*=w1;
      w*=w2;

      double wIDSF1=1.,wIDSF2=1.;
      wIDSF1=Ntp->IDSF(Tau1);
      wIDSF2=Ntp->IDSF(Tau2);
      w*=wIDSF1*wIDSF2;
    }
    // if(!Ntp->isData() && PairsIndexTemp.size()>0/*pass.at(NPairsFound)*/ && (id==33 || id == 10110333 || id == 10110433|| id == 10130533|| id ==10210333|| id == 10210433|| id == 10230533|| id ==10310333 || id ==10330533 || id ==10410433 || id == 10410333|| id == 10430533|| id == 30530533 || id ==11 || id ==12)){
    //   w *= 0.95*0.95;
    // }
    
  }
  TLorentzVector genMomentum(0,0,0,0);
  if(id==33 || id == 10110333 || id == 10110433|| id == 10130533|| id ==10210333|| id == 10210433|| id == 10230533|| id ==10310333 || id ==10330533 || id ==10410433 || id == 10410333|| id == 10430533|| id == 30530533) {
    for(unsigned int imc=0; imc < Ntp->NGenParts(); imc++){
      if(fabs(Ntp->Genpart_pdg(imc)) ==15   &&  Ntp->CHECK_BIT(Ntp->Genpart_flags(imc),0)&& Ntp->Genpart_status(imc) ==2) {
	if(Ntp->Genpart_P4(imc).Pt() > 8){
	  genMomentum+=Ntp->Genpart_P4(imc);
	}
      }
    }
  }
  if( id == 30) {
    for(unsigned int imc=0; imc < Ntp->NGenParts(); imc++){
      if((fabs(Ntp->Genpart_pdg(imc)) ==11 || fabs(Ntp->Genpart_pdg(imc)) ==13) && Ntp->Genpart_status(imc) ==1  ){
	if(Ntp->Genpart_P4(imc).Pt() > 8){
	  genMomentum+=Ntp->Genpart_P4(imc);
	}
      }
    }
  }
  float zptw(1.);
  if(!Ntp->isData() && genMomentum.Pt()!=0 && genMomentum.M() > 75 && genMomentum.M() < 120) {
    if(Ntp->year()==2016)zptw = Ntp->ZPtReweight(genMomentum);
    //if(Ntp->year()==2017)zptw = Ntp->ZPtReweight(genMomentum,w2017_);
    //if(Ntp->year()==2018)zptw = Ntp->ZPtReweight(genMomentum,w2018_);
    //zptw = DataMC_Corr.ZPTWeight(genMomentum.M(),genMomentum.Pt());
  }
  w*=zptw;
  float PUPPImetCorr_px=Ntp->PUPPImet()*cos(Ntp->PUPPImetphi());
  float PUPPImetCorr_py=Ntp->PUPPImet()*sin(Ntp->PUPPImetphi());
  if(id==33 || id == 10110333 || id == 10110433|| id == 10130533|| id ==10210333|| id == 10210433|| id == 10230533|| id ==10310333 || id ==10330533 || id ==10410433 || id == 10410333|| id == 10430533|| id == 30530533 || id==30 || id==11 || id==12 || id==20 || id==21 || id==22 ||id==23 ||id==45  ||id==46)
    {
      if(PairsIndexTemp.size()>0)
	{
	  TLorentzVector Vis,Gen;
	  for(unsigned i = 0; i < Ntp->NGenParts(); ++i) {
	    unsigned pdgid = abs(Ntp->Genpart_pdg(i));
	    unsigned status = abs(Ntp->Genpart_status(i));

	    if ( (pdgid >= 11 && pdgid <= 16 && Ntp->CHECK_BIT(Ntp->Genpart_flags(i),8) && status==1) || Ntp->CHECK_BIT(Ntp->Genpart_flags(i),10)) Gen+=Ntp->Genpart_P4(i);

	    if ( ( (pdgid == 11 || pdgid == 13 || pdgid == 15) && Ntp->CHECK_BIT(Ntp->Genpart_flags(i),8) && status==1) || (Ntp->CHECK_BIT(Ntp->Genpart_flags(i),10) && !(pdgid==12||pdgid==14||pdgid==16))) Vis+=Ntp->Genpart_P4(i) ;
	  }
	  Ntp->RecoilCorr(Gen, Vis,IndexSelectedTemp.at(0),PUPPImetCorr_px,PUPPImetCorr_py);
	}
    }
  if(!Ntp->isData() && id!=DataMCType::QCD)w*=Ntp->MC_weight();//cout<<"id: "<<id<<"MC_weight(): "<< Ntp->MC_weight()<<endl; //generator weight because negative weights for this samples

  if(id==11)w*=0.2455;
  else if(id==12)w*=0.2727;
  else if(id==45)w*=0.2546;
  //else if(Ntp->GetInputDatasetName().Contains("WminusHToTauTauUncorrelatedDecay_Filtered_M125_TuneCUETP8M1_13TeV-powheg-pythia8"))w*=0.2596;
  //else if(Ntp->GetInputDatasetName().Contains("WplusHToTauTauUncorrelatedDecay_Filtered_M125_TuneCUETP8M1_13TeV-powheg-pythia8"))w*=0.2425;
  w*=Ntp->stitch_weight(isDY1050);
  //cout<<"isDY1050: "<<isDY1050<<endl;

  if(!Ntp->isData() && (Ntp->year()==2016||Ntp->year()==2017))w*=Ntp->prefiringweight();

  std::vector<unsigned int> exclude_cuts;
  exclude_cuts.push_back(TausIsolation);
  // exclude_cuts.push_back(Tau2Isolation);
  //exclude_cuts.push_back(PairCharge);
  classic_svFit::LorentzVector tau1P4;
  classic_svFit::LorentzVector tau2P4;

  TLorentzVector Tau1P4;
  TLorentzVector Tau2P4;
  if(passAllBut(exclude_cuts)){
    Tau1P4 = Ntp->P4Corrected(Tau1);
    Tau2P4 = Ntp->P4Corrected(Tau2);
  }
  // QCD ABCD BG Method
  /*******************
   *        |        *
   *    C   |    D   *  SS
   *        |        *       S
   * ----------------*------ i
   *        |        *       g
   *    A   |    B   *  OS   n
   *        |        *
   *******************
   *  Iso   | AntiIso
   *
   *     TauIsolation
   */

  // std::cout<<" before  " << pass.at(TriggerOk) << "    " <<   pass.at(PrimeVtx) << "    " <<  pass.at(NPairsFound)<< "    " <<   pass.at(FirstTauIsolation) << "    " <<  pass.at(SecondTauIsolation) << "    " <<  pass.at(nGoodMuons) << "    " <<  pass.at(PairCharge) << "  passAllBut  " << passAllBut(exclude_cuts) <<std::endl;
  
  int pt20=-1;
  if(PairsIndexTemp.size()>0){
    for(int k=0;k<=IndexSelectedTemp.at(0);k++)
      {
	if(Ntp->njetspt20(k)>=2)pt20++;
	if(k==IndexSelectedTemp.at(0) && Ntp->njetspt20(k)<2)pt20=-99;
      }
  }
  if(passAllBut(exclude_cuts)) { 
    //    for(unsigned int ia=0; ia<pass.size(); ia++){         std::cout<<" ia  "<< ia <<  "   pass  " <<pass.at(ia) << std::endl;    }
    // if(pass.at(FirstTauIsolation) && pass.at(SecondTauIsolation)){
    //   if(pass.at(PairCharge)){
    // 	NQCD.at(t).Fill(1.,w); //A
    //   }  
    //   if(!pass.at(PairCharge)){
    // 	NQCD.at(t).Fill(2.,w); //B
    //   }
    //   if(Ntp->isIsolatedTau(TauIndex_1,"Medium") && Ntp->isIsolatedTau(TauIndex_2,"Loose")){
    // 	if(pass.at(PairCharge)){
    // 	  NQCD.at(t).Fill(3.,w); //С
    // 	}
    // 	if(!pass.at(PairCharge)){
    // 	  NQCD.at(t).Fill(4.,w); //В
    // 	}
    //   }
    // }
    if(pass.at(PairCharge)) {
      if(pass.at(TausIsolation)){
    	NQCD.at(t).Fill(1.,w); //A
      }
      if(((Ntp->byMediumDeepTau2017v2p1VSjet_1(IndexSelectedTemp.at(0))) && !Ntp->isIsolatedTau(Tau2,"Medium") && Ntp->isIsolatedTau(Tau2,"VLoose")) || (Ntp->isIsolatedTau(Tau2,"Medium") && !Ntp->isIsolatedTau(Tau1,"Medium") && Ntp->isIsolatedTau(Tau1,"VLoose"))){
      	NQCD.at(t).Fill(2.,w); //B
      	//TauTauFullMass_B.at(t).Fill((tau1P4+tau2P4).M(),w);
      }
    }
    if(!pass.at(PairCharge)) {
      if(pass.at(TausIsolation)){
    	NQCD.at(t).Fill(3.,w); //C
    	//TauTauFullMass_C.at(t).Fill((tau1P4+tau2P4).M(),w);
      }
      if((Ntp->isIsolatedTau(Tau1,"Medium") && !Ntp->isIsolatedTau(Tau2,"Medium") && Ntp->isIsolatedTau(Tau2,"VLoose")) || (Ntp->isIsolatedTau(Tau2,"Medium") && !Ntp->isIsolatedTau(Tau1,"Medium") && Ntp->isIsolatedTau(Tau1,"VLoose"))){
    	NQCD.at(t).Fill(4.,w); //D
    	//TauTauFullMass_D.at(t).Fill((tau1P4+tau2P4).M(),w);
      }
    }
    
    double PUPPIMETCorr=sqrt(PUPPImetCorr_px*PUPPImetCorr_px+PUPPImetCorr_py*PUPPImetCorr_py);
    if(Ntp->isData() && !Ntp->byMediumDeepTau2017v2p1VSjet_1(IndexSelectedTemp.at(0)) && Ntp->byVVVLooseDeepTau2017v2p1VSjet_1(IndexSelectedTemp.at(0)) && Ntp->byMediumDeepTau2017v2p1VSjet_2(IndexSelectedTemp.at(0))){
      double isOS=(Ntp->Daughters_charge(Tau1)/abs(Ntp->Daughters_charge(Tau1))) != (Ntp->Daughters_charge(Tau2)/abs(Ntp->Daughters_charge(Tau2)));
      auto args = std::vector<double>{Tau1P4.Pt(),(double)Ntp->MVADM2017(Tau1),(double)Ntp->njets(IndexSelectedTemp.at(0)),Tau2P4.Pt(),isOS,PUPPIMETCorr*cos(TLorentzVector(Tau1P4.X(),Tau1P4.Y(),Tau1P4.Z(),0.).DeltaPhi(TLorentzVector(PUPPImetCorr_px,PUPPImetCorr_py,0.,0.)))/Tau1P4.Pt()};
      double FFData2016 = std::shared_ptr<RooFunctor>(wFF2016->function("ff_tt_medium_mvadmbins_nosig")->functor(wFF2016->argSet("pt,mvadm,njets,pt_2,os,met_var_qcd")))->eval(args.data());
      Tau1PT.at(1).Fill(Tau1P4.Pt(),FFData2016);
      Tau1E.at(1).Fill(Tau1P4.E(),FFData2016);
      Tau1Mass.at(1).Fill(Tau1P4.M(),FFData2016);
      Tau1Phi.at(1).Fill(Tau1P4.Phi(),FFData2016);
      Tau1Eta.at(1).Fill(Tau1P4.Eta(),FFData2016);
      Tau1dz.at(1).Fill(Ntp->dz(Tau1),FFData2016);
      Tau1HPSDecayMode.at(1).Fill(Ntp->decayMode(Tau1),FFData2016);
      Tau1MVADecayMode.at(1).Fill(Ntp->MVADM2017(Tau1),FFData2016);
      Tau2PT.at(1).Fill(Tau2P4.Pt(),FFData2016);
      Tau2E.at(1).Fill(Tau2P4.E(),FFData2016);
      Tau2Mass.at(1).Fill(Tau2P4.M(),FFData2016);
      Tau2Phi.at(1).Fill(Tau2P4.Phi(),FFData2016);
      Tau2Eta.at(1).Fill(Tau2P4.Eta(),FFData2016);
      Tau2dz.at(1).Fill(Ntp->dz(Tau2),FFData2016);
      Tau2HPSDecayMode.at(1).Fill(Ntp->decayMode(Tau2),FFData2016);
      Tau2MVADecayMode.at(1).Fill(Ntp->MVADM2017(Tau2),FFData2016);
      PUPPImetcorr.at(1).Fill(sqrt(PUPPImetCorr_px*PUPPImetCorr_px+PUPPImetCorr_py*PUPPImetCorr_py),FFData2016);

      NbJets.at(1).Fill(Ntp->njets(IndexSelectedTemp.at(0)),FFData2016);
      TauTauVisMass.at(1).Fill((Tau1P4+Tau2P4).M(),FFData2016);
      TauTauVisPT.at(1).Fill((Tau1P4+Tau2P4).Pt(),FFData2016);
      if(pt20>=0)Mjj.at(1).Fill(Ntp->dijetpt(pt20),FFData2016);
      //cout<<"FF2016: "<<FF2016<<endl;
      // cout<<"w: "<<w<<endl;
      NFFData.at(t).Fill(1,FFData2016);
      //NFFData.at(t).Fill(1,w);
    }
    if(!Ntp->isData() && !Ntp->byMediumDeepTau2017v2p1VSjet_1(IndexSelectedTemp.at(0)) && Ntp->byVVVLooseDeepTau2017v2p1VSjet_1(IndexSelectedTemp.at(0))  && Ntp->byMediumDeepTau2017v2p1VSjet_2(IndexSelectedTemp.at(0)) && Ntp->genmatch(Tau1)!=6){
      double isOS=(Ntp->Daughters_charge(Tau1)/abs(Ntp->Daughters_charge(Tau1))) != (Ntp->Daughters_charge(Tau2)/abs(Ntp->Daughters_charge(Tau2)));
      auto args = std::vector<double>{Tau1P4.Pt(),(double)Ntp->MVADM2017(Tau1),(double)Ntp->njets(IndexSelectedTemp.at(0)),Tau2P4.Pt(),isOS,PUPPIMETCorr*cos(TLorentzVector(Tau1P4.X(),Tau1P4.Y(),0.,0.).DeltaPhi(TLorentzVector(PUPPImetCorr_px,PUPPImetCorr_py,0.,0.)))/Tau1P4.Pt()};
      double FFMC2016 = std::shared_ptr<RooFunctor>(wFF2016->function("ff_tt_medium_mvadmbins_nosig")->functor(wFF2016->argSet("pt,mvadm,njets,pt_2,os,met_var_qcd")))->eval(args.data());
      Tau1PTQCDMC.at(t).Fill(Tau1P4.Pt(),FFMC2016*w);
      Tau1EQCDMC.at(t).Fill(Tau1P4.E(),FFMC2016*w);
      Tau1MassQCDMC.at(t).Fill(Tau1P4.M(),FFMC2016*w);
      Tau1PhiQCDMC.at(t).Fill(Tau1P4.Phi(),FFMC2016*w);
      Tau1EtaQCDMC.at(t).Fill(Tau1P4.Eta(),FFMC2016*w);
      Tau1dzQCDMC.at(t).Fill(Ntp->dz(Tau1),FFMC2016*w);
      Tau1HPSDecayModeQCDMC.at(t).Fill(Ntp->decayMode(Tau1),FFMC2016*w);
      Tau1MVADecayModeQCDMC.at(t).Fill(Ntp->MVADM2017(Tau1),FFMC2016*w);
      Tau2PTQCDMC.at(t).Fill(Tau2P4.Pt(),FFMC2016*w);
      Tau2EQCDMC.at(t).Fill(Tau2P4.E(),FFMC2016*w);
      Tau2MassQCDMC.at(t).Fill(Tau2P4.M(),FFMC2016*w);
      Tau2PhiQCDMC.at(t).Fill(Tau2P4.Phi(),FFMC2016*w);
      Tau2EtaQCDMC.at(t).Fill(Tau2P4.Eta(),FFMC2016*w);
      Tau2dzQCDMC.at(t).Fill(Ntp->dz(Tau2),FFMC2016*w);
      Tau2HPSDecayModeQCDMC.at(t).Fill(Ntp->decayMode(Tau2),FFMC2016*w);
      Tau2MVADecayModeQCDMC.at(t).Fill(Ntp->MVADM2017(Tau2),FFMC2016*w);
      PUPPImetcorrQCDMC.at(t).Fill(sqrt(PUPPImetCorr_px*PUPPImetCorr_px+PUPPImetCorr_py*PUPPImetCorr_py),FFMC2016*w);
      NbJetsQCDMC.at(t).Fill(Ntp->njets(IndexSelectedTemp.at(0)),FFMC2016*w);
      TauTauVisMassQCDMC.at(t).Fill((Tau1P4+Tau2P4).M(),FFMC2016*w);
      TauTauVisPTQCDMC.at(t).Fill((Tau1P4+Tau2P4).Pt(),FFMC2016*w);
      if(pt20>=0)MjjQCDMC.at(t).Fill(Ntp->dijetpt(pt20),FFMC2016*w);
      NFFLeadMC.at(t).Fill(1,FFMC2016*w);
      
      //NFFLeadMC.at(t).Fill(1,w);
    }
    if(!Ntp->isData() && pass.at(TausIsolation) && Ntp->genmatch(Tau1)!=6 && Ntp->genmatch(Tau2)==6){
      W_res.at(t).Fill(1,w);
    }
  }
  // bool IsQCDEvent = false;
  // if(passAllBut(exclude_cuts)){
  //   if(pass.at(PairCharge)){
  //     if(!Ntp->isIsolatedTau(Tau1,"Medium") && Ntp->isIsolatedTau(Tau1,"VVVLoose")){
  //     //if((Ntp->isIsolatedTau(Tau1,"Medium") && !Ntp->isIsolatedTau(Tau2,"Medium") && Ntp->isIsolatedTau(Tau2,"VVVLoose")) || (Ntp->isIsolatedTau(Tau2,"Medium") && !Ntp->isIsolatedTau(Tau1,"Medium") && Ntp->isIsolatedTau(Tau1,"VVVLoose"))){
  // 	if(id == DataMCType::Data){
  // 	  QCDShape.at(t).Fill(1,w);
  // 	  t=HConfig.GetType(DataMCType::QCD);
  // 	  IsQCDEvent = true;
  // 	}
  //     }
  //     //}
  //   }
  // }

  //if(IsQCDEvent){ pass.at(PairCharge)= true;pass.at(TausIsolation)= true;}
  
    // std::vector<unsigned int> exclude_cuts_ForTauIso;
    // exclude_cuts_ForTauIso.push_back(TausIsolation);
    // //  exclude_cuts_ForTauIso.push_back(Tau2Isolation);
    // if(passAllBut(exclude_cuts_ForTauIso)) {
    // if(Ntp->isIsolatedTau(Tau1,"VVVLoose"))Tau1isolation.at(t).Fill(0.);
    // if(Ntp->isIsolatedTau(Tau1,"VVLoose"))Tau1isolation.at(t).Fill(1.);
    // if(Ntp->isIsolatedTau(Tau1,"VLoose"))Tau1isolation.at(t).Fill(2.);
    // if(Ntp->isIsolatedTau(Tau1,"Loose"))Tau1isolation.at(t).Fill(3.);
    // if(Ntp->isIsolatedTau(Tau1,"Medium"))Tau1isolation.at(t).Fill(4.);
    // if(Ntp->isIsolatedTau(Tau1,"Tight"))Tau1isolation.at(t).Fill(5.);
    // if(Ntp->isIsolatedTau(Tau1,"VTight"))Tau1isolation.at(t).Fill(6.);
    // if(Ntp->isIsolatedTau(Tau1,"VVTight"))Tau1isolation.at(t).Fill(7.);
    // if(Ntp->isIsolatedTau(Tau2,"VVVLoose"))Tau2isolation.at(t).Fill(0.);
    // if(Ntp->isIsolatedTau(Tau2,"VVLoose"))Tau2isolation.at(t).Fill(1.);
    // if(Ntp->isIsolatedTau(Tau2,"VLoose"))Tau2isolation.at(t).Fill(2.);
    // if(Ntp->isIsolatedTau(Tau2,"Loose"))Tau2isolation.at(t).Fill(3.);
    // if(Ntp->isIsolatedTau(Tau2,"Medium"))Tau2isolation.at(t).Fill(4.);
    // if(Ntp->isIsolatedTau(Tau2,"Tight"))Tau2isolation.at(t).Fill(5.);
    // if(Ntp->isIsolatedTau(Tau2,"VTight"))Tau2isolation.at(t).Fill(6.);
    // if(Ntp->isIsolatedTau(Tau2,"VVTight"))Tau2isolation.at(t).Fill(7.);
    // }
  
  bool status=AnalysisCuts(t,w,wobs);  // boolean that say whether your event passed critera defined in pass vector. The whole vector must be true for status = true
  ///////////////////////////////////////////////////////////
  // Analyse events which passed selection
  if(status) {
    double pvx(0);
    pvx =  Ntp->npv();
    // if(id == DataMCType::Data) pvx =  Ntp->npv();
    if(id !=DataMCType::Data && id !=DataMCType::QCD)	  pvx = Ntp->PUNumInteractions();
    NPrimeVtx.at(t).Fill(pvx,w);
    NPU.at(t).Fill(Ntp->npu(),w);
    RHO.at(t).Fill(Ntp->rho(),w);
    

    std::vector<int> thirdLepton;


    //    TLorentzVector taunew(tau1P4.Px(), );

    //svfTau1E.at(t).Fill(tau1P4.E(),w);
    //svfTau2E.at(t).Fill(tau2P4.E(),w);

    Tau1PT.at(t).Fill(Tau1P4.Pt(),w);
    Tau1E.at(t).Fill(Tau1P4.E(),w);
    Tau1Mass.at(t).Fill(Tau1P4.M(),w);
    Tau1Phi.at(t).Fill(Tau1P4.Phi(),w);
    Tau1Eta.at(t).Fill(Tau1P4.Eta(),w);
    Tau1dz.at(t).Fill(Ntp->dz(Tau1),w);
    Tau1HPSDecayMode.at(t).Fill(Ntp->decayMode(Tau1),w);
    Tau1MVADecayMode.at(t).Fill(Ntp->MVADM2017(Tau1),w);
    Tau2PT.at(t).Fill(Tau2P4.Pt(),w);
    Tau2E.at(t).Fill(Tau2P4.E(),w);
    Tau2Mass.at(t).Fill(Tau2P4.M(),w);
    Tau2Phi.at(t).Fill(Tau2P4.Phi(),w);
    Tau2Eta.at(t).Fill(Tau2P4.Eta(),w);
    Tau2dz.at(t).Fill(Ntp->dz(Tau2),w);
    Tau2HPSDecayMode.at(t).Fill(Ntp->decayMode(Tau2),w);
    Tau2MVADecayMode.at(t).Fill(Ntp->MVADM2017(Tau2),w);
    /*
      againstElectronVLooseMVA6_Tau1.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau1),Ntp->Bit_againstElectronVLooseMVA6),w);
      againstElectronLooseMVA6_Tau1.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau1),Ntp->Bit_againstElectronLooseMVA6),w);
      againstElectronMediumMVA6_Tau1.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau1),Ntp->Bit_againstElectronMediumMVA6),w);
      againstElectronTightMVA6_Tau1.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau1),Ntp->Bit_againstElectronTightMVA6),w);
      againstElectronVTightMVA6_Tau1.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau1),Ntp->Bit_againstElectronVTightMVA6),w);
      againstMuonLoose3_Tau1.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau1),Ntp->Bit_againstMuonLoose3),w);
      againstMuonTight3_Tau1.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau1),Ntp->Bit_againstMuonTight3),w);
      byCombinedIsolationDeltaBetaCorrRaw3Hits_Tau1.at(t).Fill(Ntp->Daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits(Tau1),w);

      againstElectronVLooseMVA6_Tau2.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau2),Ntp->Bit_againstElectronVLooseMVA6),w);
      againstElectronLooseMVA6_Tau2.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau2),Ntp->Bit_againstElectronLooseMVA6),w);
      againstElectronMediumMVA6_Tau2.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau2),Ntp->Bit_againstElectronMediumMVA6),w);
      againstElectronTightMVA6_Tau2.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau2),Ntp->Bit_againstElectronTightMVA6),w);
      againstElectronVTightMVA6_Tau2.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau2),Ntp->Bit_againstElectronVTightMVA6),w);
      againstMuonLoose3_Tau2.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau2),Ntp->Bit_againstMuonLoose3),w);
      againstMuonTight3_Tau2.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau2),Ntp->Bit_againstMuonTight3),w);
      byCombinedIsolationDeltaBetaCorrRaw3Hits_Tau2.at(t).Fill(Ntp->Daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits(Tau2),w);
    */
    for(int iDaughter=0;   iDaughter  <  Ntp->NDaughters() ;iDaughter++ ) {
      if((iDaughter!=Tau1)&&(iDaughter!=Tau2)){
	if(Ntp->ElectronVeto(iDaughter) || Ntp->MuonVeto(iDaughter))thirdLepton.push_back(iDaughter);
      }
    }
    if(thirdLepton.size()>0)ExtraLeptonVeto.at(t).Fill(1.,w);
    else ExtraLeptonVeto.at(t).Fill(0.,w);

    TauTauVisMass.at(t).Fill((Tau1P4+Tau2P4).M(),w);
    //TauTauFullMass.at(t).Fill((tau1P4+tau2P4).M(),w);
    dRTauTau.at(t).Fill(Tau1P4.DeltaR(Tau2P4),w);
    
    MET.at(t).Fill(Ntp->MET(),w);
    METphi.at(t).Fill(Ntp->METphi(),w);
    PUPPImet.at(t).Fill(Ntp->PUPPImet(),w);
    PUPPImetphi.at(t).Fill(Ntp->PUPPImetphi(),w);
    
    PUPPImetcorr.at(t).Fill(sqrt(PUPPImetCorr_px*PUPPImetCorr_px+PUPPImetCorr_py*PUPPImetCorr_py),w);
    PUPPImetcorrphi.at(t).Fill(TMath::ATan2(PUPPImetCorr_py,PUPPImetCorr_px),w);
    
    TransverseMass.at(t).Fill(Ntp->transverseMass(Tau1P4.Pt(), Tau1P4.Phi(), Tau2P4.Pt(), Tau2P4.Phi()),w);
    NbJets.at(t).Fill(Ntp->njets(IndexSelectedTemp.at(0)),w);
    
    if(pt20>=0)Mjj.at(t).Fill(Ntp->dijetpt(pt20),w);
    TLorentzVector Tauplusvis;
    TLorentzVector Tauminusvis;
    TLorentzVector Pi0RECO;
    TLorentzVector Tauplustruth;
    TLorentzVector Tauminustruth;
    unsigned int Tauplus=0;
    unsigned int Tauminus=0;

    if(Ntp->Daughters_charge(Tau1)>0)
      {
	//Tauplussvfit=tau1P4;
	//Tauminussvfit=tau2P4;
	Tauplusvis=Tau1P4;
	Tauminusvis=Tau2P4;
	Tauplus=Tau1;
	Tauminus=Tau2;
      }
    else
      {
	//Tauplussvfit=tau2P4;
        //Tauminussvfit=tau1P4;
	Tauplusvis=Tau2P4;
	Tauminusvis=Tau1P4;
	Tauplus=Tau2;
	Tauminus=Tau1;
      }

    ZPtVis.at(t).Fill((Tau1P4+Tau2P4).Pt(),w);
    TauTauVisPT.at(t).Fill((Tau1P4+Tau2P4).Pt(),w);

    TLorentzVector Tau1Truth; 
    TLorentzVector Tau2Truth;
    TLorentzVector TruthDecayFromTau1;
    TLorentzVector TruthDecayFromTau2; 
    std::vector<TLorentzVector> Pions1;
    std::vector<TLorentzVector> Pions2;
    std::vector<double> Pions1Charge;
    std::vector<double> Pions2Charge;
    bool decay=0;
	
    if(Ntp->CheckDecayID(1,3)){
	
      Tau1Truth=Ntp->GetTruthTauLV(1,0);
      Tau2Truth=Ntp->GetTruthTauLV(3,1);
      TruthDecayFromTau1=Ntp->GetTruthTauProductLV(1,11,0);
      TruthDecayFromTau2=Ntp->GetTruthTauProductLV(3,211,1);decay=1;
      if(TruthDecayFromTau1==TruthDecayFromTau2)cout<<"Same decay particle 1-3"<<endl;
    }
    if(Ntp->CheckDecayID(1,4)){
	
      Tau1Truth=Ntp->GetTruthTauLV(1,0);
      Tau2Truth=Ntp->GetTruthTauLV(4,1);
      TruthDecayFromTau1=Ntp->GetTruthTauProductLV(1,11,0);
      Pions2=Ntp->GetTruthPionsFromRho(1);
      TruthDecayFromTau2=Pions2.at(0)+Pions2.at(1);decay=1;
      if(TruthDecayFromTau1==TruthDecayFromTau2)cout<<"Same decay particle 1-4"<<endl;
    }
    if(Ntp->CheckDecayID(1,5)){
      Tau1Truth=Ntp->GetTruthTauLV(1,0);
      Tau2Truth=Ntp->GetTruthTauLV(5,1);
      TruthDecayFromTau1=Ntp->GetTruthTauProductLV(1,11,0);
      Pions2=Ntp->GetTruthPionsFromA1(1);
      TruthDecayFromTau2=Pions2.at(0)+Pions2.at(1)+Pions2.at(2);decay=1;
      if(TruthDecayFromTau1==TruthDecayFromTau2)cout<<"Same decay particle 1-5"<<endl;
    }


    if(Ntp->CheckDecayID(2,3)){
      Tau1Truth=Ntp->GetTruthTauLV(2,0);
      Tau2Truth=Ntp->GetTruthTauLV(3,1);
      TruthDecayFromTau1=Ntp->GetTruthTauProductLV(2,13,0);
      TruthDecayFromTau2=Ntp->GetTruthTauProductLV(3,211,1);decay=1;
      if(TruthDecayFromTau1==TruthDecayFromTau2)cout<<"Same decay particle 2-3"<<endl;
    }
    if(Ntp->CheckDecayID(2,4)){
      Tau1Truth=Ntp->GetTruthTauLV(2,0);
      Tau2Truth=Ntp->GetTruthTauLV(4,1);
      TruthDecayFromTau1=Ntp->GetTruthTauProductLV(2,13,0);
      Pions2=Ntp->GetTruthPionsFromRho(1);
      TruthDecayFromTau2=Pions2.at(0)+Pions2.at(1);decay=1;
      if(TruthDecayFromTau1==TruthDecayFromTau2)cout<<"Same decay particle 2-4"<<endl;
    }
    if(Ntp->CheckDecayID(2,5)){
      Tau1Truth=Ntp->GetTruthTauLV(2,0);
      Tau2Truth=Ntp->GetTruthTauLV(5,1);
      TruthDecayFromTau1=Ntp->GetTruthTauProductLV(2,13,0);
      Pions2=Ntp->GetTruthPionsFromA1(1);
      TruthDecayFromTau2=Pions2.at(0)+Pions2.at(1)+Pions2.at(2);decay=1;
      if(TruthDecayFromTau1==TruthDecayFromTau2)cout<<"Same decay particle 2-5"<<endl;
    }

    if(Ntp->CheckDecayID(3,3)){
      Tau1Truth=Ntp->GetTruthTauLV(3,0);
      Tau2Truth=Ntp->GetTruthTauLV(3,1);
      TruthDecayFromTau1=Ntp->GetTruthTauProductLV(3,211,0);
      TruthDecayFromTau2=Ntp->GetTruthTauProductLV(3,211,1);decay=1;
      if(TruthDecayFromTau1==TruthDecayFromTau2)cout<<"Same decay particle 3-3"<<endl;

    }
    if(Ntp->CheckDecayID(3,5)){
      Tau1Truth=Ntp->GetTruthTauLV(3,0);
      Tau2Truth=Ntp->GetTruthTauLV(5,1);
      TruthDecayFromTau1=Ntp->GetTruthTauProductLV(3,211,0);
      Pions2=Ntp->GetTruthPionsFromA1(1);
      TruthDecayFromTau2=Pions2.at(0)+Pions2.at(1)+Pions2.at(2);decay=1;
      if(TruthDecayFromTau1==TruthDecayFromTau2)cout<<"Same decay particle 3-5"<<endl;
    }


    if(Ntp->CheckDecayID(4,4)){
      Tau1Truth=Ntp->GetTruthTauLV(4,0);
      Tau2Truth=Ntp->GetTruthTauLV(4,1);
      Pions1=Ntp->GetTruthPionsFromRho(0);
      TruthDecayFromTau1=Pions1.at(0)+Pions1.at(1);
      Pions2=Ntp->GetTruthPionsFromRho(1);
      TruthDecayFromTau2=Pions2.at(0)+Pions2.at(1);decay=1;
      if(TruthDecayFromTau1==TruthDecayFromTau2)cout<<"Same decay particle 4-4"<<endl;
    }
    if(Ntp->CheckDecayID(4,3)){
      Tau1Truth=Ntp->GetTruthTauLV(4,0);
      Tau2Truth=Ntp->GetTruthTauLV(3,1);
      Pions1=Ntp->GetTruthPionsFromRho(0);
      TruthDecayFromTau1=Pions1.at(0)+Pions1.at(1);
      TruthDecayFromTau2=Ntp->GetTruthTauProductLV(3,211,1);decay=1;
      if(TruthDecayFromTau1==TruthDecayFromTau2)cout<<"Same decay particle 4-3"<<endl;

    }
    if(Ntp->CheckDecayID(4,5)){
      Tau1Truth=Ntp->GetTruthTauLV(4,0);
      Tau2Truth=Ntp->GetTruthTauLV(5,1);
      Pions1=Ntp->GetTruthPionsFromRho(0);
      TruthDecayFromTau1=Pions1.at(0)+Pions1.at(1);
      Pions2=Ntp->GetTruthPionsFromA1(1);
      TruthDecayFromTau2=Pions2.at(0)+Pions2.at(1)+Pions2.at(2);decay=1;
      if(TruthDecayFromTau1==TruthDecayFromTau2)cout<<"Same decay particle 4-5"<<endl;

    }


    if(Ntp->CheckDecayID(5,5) ){
      Tau1Truth=Ntp->GetTruthTauLV(5,0);
      Tau2Truth=Ntp->GetTruthTauLV(5,1);
      Pions1=Ntp->GetTruthPionsFromA1(0);
      TruthDecayFromTau1=Pions1.at(0)+Pions1.at(1)+Pions1.at(2);
      Pions2=Ntp->GetTruthPionsFromA1(1);
      TruthDecayFromTau2=Pions2.at(0)+Pions2.at(1)+Pions2.at(2);decay=1;
    }
    if (decay==1)
      {
	if(Tau1Truth==Tau2Truth)cout<<"Same Taus"<<endl;

	TLorentzVector TruthZ    = Tau1Truth+Tau2Truth;
	TLorentzVector VisDecayfromTauplus;
	TLorentzVector VisDecayfromTauminus;
	  
	//TauTauTruthMass.at(t).Fill((Tau1Truth+Tau2Truth).M(),w);
	Tauplustruth=Tau2Truth;
	Tauminustruth=Tau1Truth;
	VisDecayfromTauplus=TruthDecayFromTau2;
	VisDecayfromTauminus=TruthDecayFromTau1;
	    
      }

    TVector3 tauPrimaryVertex , tauNoBSPrimaryVertex, tauBSPrimaryVertex, tauNoBSZNominalPrimaryVertex, tauBSZNominalPrimaryVertex, TauminusSecondaryVertex , TauplusSecondaryVertex;

    TVector3 TauminusDirection , TauplusDirection, TauplusDirectionNoBS, TauplusDirectionBS, TauplusDirectionNoBSZNominal, TauplusDirectionBSZNominal,TauminusDirectionNoBS, TauminusDirectionBS, TauminusDirectionNoBSZNominal, TauminusDirectionBSZNominal;

    double thetaGJ_Tauminus , thetaGJ_Tauplus;
    TLorentzVector a1LV_Tauminus , a1LV_Tauplus, a1LVRefit_Tauminus , a1LVRefit_Tauplus;
    TLorentzVector TauminusPairConstraint, TauplusPairConstraintNoBS, TauplusPairConstraintBS, TauplusPairConstraintNoBSZNominal, TauplusPairConstraintBSZNominal, TauminusPairConstraintNoBS, TauminusPairConstraintBS, TauminusPairConstraintNoBSZNominal, TauminusPairConstraintBSZNominal;
    TLorentzVector TauplusSmall, TauplusLarge, TauplusPairConstraint,TauplusPairConstraintMVA,TauplusPairConstraintNoBSNewMVA,TauplusPairConstraintBSNewMVA;
    bool isPlusReal=true, isMinusReal=true,  a1a1=false,a1a1MVA=false, a1a1TruthSVFit=false,  a1a1TruthSVFitMVA =false, a1rhoTruthSVFit  =false,  a1rhoTruthSVFitMVA =false, a1piTruthSVFit =false,  a1piTruthSVFitMVA, piminus=false, piplus=false, rhominus=false, rhoplus=false, a1minus=false, a1plus=false,a1minusMVA=false,a1plusMVA=false, piminusMVA=false, piplusMVA=false, rhominusMVA=false, rhoplusMVA=false, a1minusTruthSVFitMVA=false,a1plusTruthSVFitMVA=false, a1minusTruthSVFit=false,a1plusTruthSVFit=false, piminusTruthSVFitMVA=false,piplusTruthSVFitMVA=false, piminusTruthSVFit=false,piplusTruthSVFit=false, rhominusTruthSVFitMVA=false,rhoplusTruthSVFitMVA=false, rhominusTruthSVFit=false,rhoplusTruthSVFit=false;
    std::vector<TLorentzVector> solutions, solutionsNoBS, solutionsBS, solutionsNoBSZNominal, solutionsBSZNominal;
    TLorentzVector ZPairConstraint;

    std::vector<size_t> hashes;
    size_t hash = 0;


    bool hasNoBSNewMVA=false, hasBSNewMVA=false;
    bool isRefitNoBS=true;
    bool isRefitBS=true;
    bool isRefitNoBSZNominal=true;
    bool isRefitBSZNominal=true;

    double Spin_WT=Ntp->TauSpinerGet(TauSpinerInterface::Spin);
    double UnSpin_WT=Ntp->TauSpinerGet(TauSpinerInterface::UnSpin);
    double FlipSpin_WT=Ntp->TauSpinerGet(TauSpinerInterface::FlipSpin);
    double hplus=Ntp->TauSpinerGet(TauSpinerInterface::hplus);
    double hminus=Ntp->TauSpinerGet(TauSpinerInterface::hminus);//1-hplus;
		 
    double Wspin=w*Spin_WT;
    
    if(Ntp->decayMode(Tauminus) == 10 &&  Ntp->PFTau_hassecondaryVertex(Tauminus) && Ntp->PFtauHasThreePions(Tauminus))a1minus=true;
    if(Ntp->decayMode(Tauplus) == 10 && Ntp->PFTau_hassecondaryVertex(Tauplus) && Ntp->PFtauHasThreePions(Tauplus))a1plus=true;
    
    if(Ntp->MVADM2017(Tauminus) == 10 &&  Ntp->PFTau_hassecondaryVertex(Tauminus) && Ntp->PFtauHasThreePions(Tauminus))a1minusMVA=true;
    if(Ntp->MVADM2017(Tauplus) == 10 && Ntp->PFTau_hassecondaryVertex(Tauplus) && Ntp->PFtauHasThreePions(Tauplus))a1plusMVA=true;
    
    if(a1minus && a1plus ) a1a1=true;  //a1-a1
    if(a1minusMVA && a1plusMVA ) a1a1MVA=true;  //a1-a1

    //TLorentzVector Tauplussvfit;
    //TLorentzVector Tauminussvfit;
    
    // ClassicSVfit svfitAlgo1;
    // //FastMTT FastMTTAlgo;
    // if(a1a1TruthSVFit || a1a1TruthSVFitMVA /*|| a1rhoTruthSVFit || a1rhoTruthSVFitMVA ||a1piTruthSVFit || a1piTruthSVFitMVA*/)
    //   {
    // 	// // //---------  svfit ---------------------
    // 	std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;
    // 	classic_svFit::MeasuredTauLepton lep1(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Tau1P4.Pt(), Tau1P4.Eta(),  Tau1P4.Phi(), Tau1P4.M(),10);
    // 	classic_svFit::MeasuredTauLepton lep2(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Tau2P4.Pt(), Tau2P4.Eta(),  Tau2P4.Phi(), Tau2P4.M(),10);
	
    // 	measuredTauLeptons.push_back(lep1);
    // 	measuredTauLeptons.push_back(lep2);
    // 	TMatrixD metcov(2,2);
    // 	double metx = Ntp->MET()*cos(Ntp->METphi());
    // 	double mety = Ntp->MET()*sin(Ntp->METphi());
	
    // 	metcov[0][0] = Ntp->PFMETCov00();
    // 	metcov[1][0] = Ntp->PFMETCov01();
    // 	metcov[0][1] = Ntp->PFMETCov10();
    // 	metcov[1][1] = Ntp->PFMETCov11();
	
    // 	svfitAlgo1.setHistogramAdapter(new classic_svFit::TauTauHistogramAdapter());
    // 	svfitAlgo1.addLogM_fixed(true,5.0);
    // 	svfitAlgo1.setDiTauMassConstraint(125.10);
    // 	//FastMTTAlgo.run(measuredTauLeptons, metx, mety, metcov);
    // 	svfitAlgo1.integrate(measuredTauLeptons,metx,mety, metcov );
    // 	if(svfitAlgo1.isValidSolution()){
	  
    // 	  double higgsmass  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svfitAlgo1.getHistogramAdapter())->getMass();
    // 	  h_SVFitMass.at(t).Fill(higgsmass,w); 
    // 	  //tau1P4 = FastMTTAlgo.getTau1P4();
    // 	  //tau2P4 = FastMTTAlgo.getTau2P4();
    // 	  tau1P4 = static_cast<classic_svFit::TauTauHistogramAdapter*>(svfitAlgo1.getHistogramAdapter())->GetFittedTau1LV();
	  
    // 	  tau2P4 = static_cast<classic_svFit::TauTauHistogramAdapter*>(svfitAlgo1.getHistogramAdapter())->GetFittedTau2LV();
	  
    	
    // 	  // ClassicSVfit svfitAlgo2;
    // 	  // svfitAlgo2.setHistogramAdapter(new classic_svFit::TauTauHistogramAdapter());
    // 	  // svfitAlgo2.addLogM_fixed(true, 5.);
    // 	  // svfitAlgo2.integrate(measuredTauLeptons,metx,mety, metcov );
    // 	  // tau1P4 = static_cast<classic_svFit::TauTauHistogramAdapter*>(svfitAlgo2.getHistogramAdapter())->GetFittedTau1LV();
    // 	  // tau2P4 = static_cast<classic_svFit::TauTauHistogramAdapter*>(svfitAlgo2.getHistogramAdapter())->GetFittedTau2LV();
	
    // 	  // //---------  svfit ---------------------
    // 	  if(svfitAlgo1.isValidSolution()){
    // 	    if(Ntp->Daughters_charge(Tau1)>0)
    // 	      {
    // 	 	Tauplussvfit.SetPxPyPzE(tau1P4.x(),tau1P4.y(),tau1P4.z(),tau1P4.t());
    // 	 	Tauminussvfit.SetPxPyPzE(tau2P4.x(),tau2P4.y(),tau2P4.z(),tau2P4.t());
	
    // 	      }
    // 	    else
    // 	      {
    // 	 	Tauplussvfit.SetPxPyPzE(tau2P4.x(),tau2P4.y(),tau2P4.z(),tau2P4.t());
    // 	 	Tauminussvfit.SetPxPyPzE(tau1P4.x(),tau1P4.y(),tau1P4.z(),tau1P4.t());
	
    // 	      }
    // 	  }
    // 	}
    //}

    if((a1a1 ||a1a1MVA) && std::isnan(Wspin)!=true && Ntp->PFTauRefit_PionsP4_SizePions(Tauminus)==3 && Ntp->PFTauRefit_PionsP4_SizePions(Tauplus))
      {
	std::vector<double> PVRefitNoBS_X_temp, PVRefitNoBS_Y_temp, PVRefitNoBS_Z_temp;
	for(unsigned int i=0;i<Ntp->NPVRefitNoBS();i++)
	  {
	    PVRefitNoBS_X_temp.push_back(Ntp->RefitPVNoBS_x(i));
	    PVRefitNoBS_Y_temp.push_back(Ntp->RefitPVNoBS_y(i));
	    PVRefitNoBS_Z_temp.push_back(Ntp->RefitPVNoBS_z(i));
	  }

	std::vector<double> PVRefitBS_X_temp, PVRefitBS_Y_temp, PVRefitBS_Z_temp;
	for(unsigned int i=0;i<Ntp->NPVRefitBS();i++)
	  {
	    PVRefitBS_X_temp.push_back(Ntp->RefitPVBS_x(i));
	    PVRefitBS_Y_temp.push_back(Ntp->RefitPVBS_y(i));
	    PVRefitBS_Z_temp.push_back(Ntp->RefitPVBS_z(i));
	  }
	

	vector<size_t> VertexHashNoBS1_temp, VertexHashNoBS2_temp;
	for(unsigned int i=0;i<Ntp->NVertexHashNoBS();i++)
	  {
	    VertexHashNoBS1_temp.push_back(Ntp->VertexHashNoBS1(i));
	    VertexHashNoBS2_temp.push_back(Ntp->VertexHashNoBS2(i));
	    
	  }
	vector<size_t> VertexHashBS1_temp, VertexHashBS2_temp;
	for(unsigned int i=0;i<Ntp->NVertexHashBS();i++)
	  {
	    VertexHashBS1_temp.push_back(Ntp->VertexHashBS1(i));
	    VertexHashBS2_temp.push_back(Ntp->VertexHashBS2(i));
	    
	  }

	boost::hash_combine(hash, Ntp->LeptonHash(Tauminus));
	boost::hash_combine(hash, Ntp->LeptonHash(Tauplus));
	hashes.push_back(hash);
	hash = 0;
	boost::hash_combine(hash, Ntp->LeptonHash(Tauplus));
	boost::hash_combine(hash, Ntp->LeptonHash(Tauminus));
	hashes.push_back(hash);

	TauminusSecondaryVertex = Ntp->PFTau_secondaryVertex_pos(Tauminus);
	TauplusSecondaryVertex = Ntp->PFTau_secondaryVertex_pos(Tauplus);

	a1LV_Tauminus = Ntp->PFTau_PionsP4(Tauminus,0) + Ntp->PFTau_PionsP4(Tauminus,1) + Ntp->PFTau_PionsP4(Tauminus,2);
	a1LV_Tauplus = Ntp->PFTau_PionsP4(Tauplus,0) + Ntp->PFTau_PionsP4(Tauplus,1) + Ntp->PFTau_PionsP4(Tauplus,2);

	a1LVRefit_Tauminus = Ntp->PFTauRefit_PionsP4(Tauminus,0) + Ntp->PFTauRefit_PionsP4(Tauminus,1) + Ntp->PFTauRefit_PionsP4(Tauminus,2);
	a1LVRefit_Tauplus = Ntp->PFTauRefit_PionsP4(Tauplus,0) + Ntp->PFTauRefit_PionsP4(Tauplus,1) + Ntp->PFTauRefit_PionsP4(Tauplus,2);

	tauPrimaryVertex=Ntp->PVtx();

	tauNoBSPrimaryVertex = GetRefittedPV(hashes, tauPrimaryVertex,PVRefitNoBS_X_temp ,PVRefitNoBS_Y_temp,PVRefitNoBS_Z_temp ,VertexHashNoBS1_temp, VertexHashNoBS2_temp,isRefitNoBS);
	tauBSPrimaryVertex = GetRefittedPV(hashes, tauPrimaryVertex, PVRefitBS_X_temp ,PVRefitBS_Y_temp,PVRefitBS_Z_temp ,VertexHashBS1_temp, VertexHashBS2_temp,isRefitBS);
	tauNoBSZNominalPrimaryVertex = GetRefittedPV(hashes, tauPrimaryVertex, PVRefitNoBS_X_temp ,PVRefitNoBS_Y_temp, Ntp->pv_z() ,VertexHashNoBS1_temp, VertexHashNoBS2_temp,isRefitNoBSZNominal);
	tauBSZNominalPrimaryVertex = GetRefittedPV(hashes, tauPrimaryVertex,PVRefitBS_X_temp ,PVRefitBS_Y_temp , Ntp->pv_z() ,VertexHashBS1_temp, VertexHashBS2_temp,isRefitBSZNominal);

	TauminusDirection = TauminusSecondaryVertex - tauPrimaryVertex;
	TauplusDirection = TauplusSecondaryVertex - tauPrimaryVertex;

	solutionsNoBS=tauPairMomentumSolutions(TauminusSecondaryVertex-tauNoBSPrimaryVertex, a1LVRefit_Tauminus, a1LV_Tauminus, isMinusReal, TauplusSecondaryVertex-tauNoBSPrimaryVertex, a1LVRefit_Tauplus, a1LVRefit_Tauplus, isPlusReal,isRefitNoBS); 
	solutionsBS=tauPairMomentumSolutions(TauminusSecondaryVertex-tauBSPrimaryVertex, a1LVRefit_Tauminus, a1LV_Tauminus, isMinusReal, TauplusSecondaryVertex-tauBSPrimaryVertex, a1LVRefit_Tauplus, a1LVRefit_Tauplus, isPlusReal,isRefitBS); 
	solutionsNoBSZNominal=tauPairMomentumSolutions(TauminusSecondaryVertex-tauNoBSZNominalPrimaryVertex, a1LVRefit_Tauminus, a1LV_Tauminus, isMinusReal, TauplusSecondaryVertex-tauNoBSZNominalPrimaryVertex, a1LVRefit_Tauplus, a1LVRefit_Tauplus, isPlusReal,isRefitNoBSZNominal); 
	solutionsBSZNominal=tauPairMomentumSolutions(TauminusSecondaryVertex-tauBSZNominalPrimaryVertex, a1LVRefit_Tauminus, a1LV_Tauminus, isMinusReal, TauplusSecondaryVertex-tauBSZNominalPrimaryVertex, a1LVRefit_Tauplus, a1LVRefit_Tauplus, isPlusReal,isRefitBSZNominal); 

	solutions=tauPairMomentumSolutions(TauminusDirection, a1LVRefit_Tauminus, a1LV_Tauminus, isMinusReal, TauplusDirection, a1LVRefit_Tauplus, a1LV_Tauplus, isPlusReal,false);

	TauminusPairConstraint = solutions.at(3);
	TauplusPairConstraint = solutions.at(7);

	TauminusPairConstraintNoBS=solutionsNoBS.at(3);
	TauminusPairConstraintBS=solutionsBS.at(3);
	TauminusPairConstraintNoBSZNominal=solutionsNoBSZNominal.at(3);
	TauminusPairConstraintBSZNominal=solutionsBSZNominal.at(3);

	TauplusPairConstraintNoBS=solutionsNoBS.at(7);
	TauplusPairConstraintBS=solutionsBS.at(7);
	TauplusPairConstraintNoBSZNominal=solutionsNoBSZNominal.at(7);
	TauplusPairConstraintBSZNominal=solutionsBSZNominal.at(7);

	ZPairConstraint= TauplusPairConstraint+TauminusPairConstraint;
      }
   
    std::vector<int> decay0;
    std::vector<int> decay1;
    bool a10=false;
    bool a11=false;
    bool a1pi00=false;
    bool a1pi01=false;
    bool rho0=false;
    bool rho1=false;
    bool pi0=false;
    bool pi1=false;
    bool pi23pi00=false;
    bool pi23pi01=false;
    bool e0=false;
    bool e1=false;
    bool mu0=false;
    bool mu1=false;

    bool purityDM=false;
    bool purityNewMVA=false;

   
    if(!Ntp->isData()){
      for(int i=0;i<Ntp->NMCTauDecayProducts(0);i++)
	{
	  decay0.push_back(Ntp->MCTauandProd_pdgid(0,i));
	}
      for(int i=0;i<Ntp->NMCTauDecayProducts(1);i++)
	{
	  decay1.push_back(Ntp->MCTauandProd_pdgid(1,i));
	}
    
      if(((((count(decay0.begin(), decay0.end(), 211)==2) && (count(decay0.begin(), decay0.end(), -211)==1)) ||((count(decay0.begin(), decay0.end(), -211)==2) && (count(decay0.begin(), decay0.end(), 211)==1))) && (count(decay0.begin(), decay0.end(), 111)==0))==true) a10=true;
      if((((count(decay0.begin(), decay0.end(), 111)==1) && (count(decay0.begin(), decay0.end(), 211)==1) && (count(decay0.begin(), decay0.end(), -211)==0)) ||((count(decay0.begin(), decay0.end(), 111)==1) && (count(decay0.begin(), decay0.end(), -211)==1) && (count(decay0.begin(), decay0.end(), 211)==0)))==true)rho0=true;
      if((((count(decay0.begin(), decay0.end(), 211)==1 && count(decay0.begin(), decay0.end(), -211)==0) ||(count(decay0.begin(), decay0.end(), -211)==1 && count(decay0.begin(), decay0.end(), 211)==0)) && count(decay0.begin(), decay0.end(), 111)==0)==true)pi0=true;
      if((((count(decay0.begin(), decay0.end(), 11)==1 && count(decay0.begin(), decay0.end(), -12)==1) ||(count(decay0.begin(), decay0.end(), -11)==1 && count(decay0.begin(), decay0.end(), 12)==1)) && count(decay0.begin(), decay0.end(), 111)==0) ==true)e0=true;
      if((((count(decay0.begin(), decay0.end(), 13)==1 && count(decay0.begin(), decay0.end(), -14)==1) ||(count(decay0.begin(), decay0.end(), -13)==1 && count(decay0.begin(), decay0.end(), 14)==1)) && count(decay0.begin(), decay0.end(), 111)==0) ==true)mu0=true;
    
      if(((((count(decay0.begin(), decay0.end(), 211)==2) && (count(decay0.begin(), decay0.end(), -211)==1)) ||((count(decay0.begin(), decay0.end(), -211)==2) && (count(decay0.begin(), decay0.end(), 211)==1))) && (count(decay0.begin(), decay0.end(), 111)==1))==true) a1pi00=true;
      if((((count(decay0.begin(), decay0.end(), 211)==1 && count(decay0.begin(), decay0.end(), -211)==0) ||(count(decay0.begin(), decay0.end(), -211)==1 && count(decay0.begin(), decay0.end(), 211)==0)) && ((count(decay0.begin(), decay0.end(), 111)==2) || count(decay0.begin(), decay0.end(), 111)==3))==true)pi23pi00=true;

      if(((((count(decay1.begin(), decay1.end(), 211)==2) && (count(decay1.begin(), decay1.end(), -211)==1)) ||((count(decay1.begin(), decay1.end(), -211)==2) && (count(decay1.begin(), decay1.end(), 211)==1))) && count(decay1.begin(), decay1.end(), 111)==0)==true) a11=true;
      if((((count(decay1.begin(), decay1.end(), 111)==1) && (count(decay1.begin(), decay1.end(), 211)==1) && (count(decay1.begin(), decay1.end(), -211)==0)) ||((count(decay1.begin(), decay1.end(), 111)==1) && (count(decay1.begin(), decay1.end(), -211)==1) && (count(decay1.begin(), decay1.end(), 211)==0)))==true)rho1=true;
      if((((count(decay1.begin(), decay1.end(), 211)==1 && count(decay1.begin(), decay1.end(), -211)==0) ||(count(decay1.begin(), decay1.end(), -211)==1 && count(decay1.begin(), decay1.end(), 211)==0)) && count(decay1.begin(), decay1.end(), 111)==0)==true)pi1=true;
      if((((count(decay1.begin(), decay1.end(), 11)==1 && count(decay1.begin(), decay1.end(), -12)==1) ||(count(decay1.begin(), decay1.end(), -11)==1 && count(decay1.begin(), decay1.end(), 12)==1)) && count(decay1.begin(), decay1.end(), 111)==0) ==true)e1=true;
      if((((count(decay1.begin(), decay1.end(), 13)==1 && count(decay1.begin(), decay1.end(), -14)==1) ||(count(decay1.begin(), decay1.end(), -13)==1 && count(decay1.begin(), decay1.end(), 14)==1)) && count(decay1.begin(), decay1.end(), 111)==0) ==true)mu1=true;
     
      if(((((count(decay1.begin(), decay1.end(), 211)==2) && (count(decay1.begin(), decay1.end(), -211)==1)) ||((count(decay1.begin(), decay1.end(), -211)==2) && (count(decay1.begin(), decay1.end(), 211)==1))) && (count(decay1.begin(), decay1.end(), 111)==1))==true) a1pi01=true;
      if((((count(decay1.begin(), decay1.end(), 211)==1 && count(decay1.begin(), decay1.end(), -211)==0) ||(count(decay1.begin(), decay1.end(), -211)==1 && count(decay1.begin(), decay1.end(), 211)==0)) && ((count(decay1.begin(), decay1.end(), 111)==2) || count(decay1.begin(), decay1.end(), 111))==3)==true)pi23pi01=true;
    }
    bool a1a1_=(a10 && a11);
    bool a1rho=((a10 && rho1) ||(a11 && rho0));
    bool a1pi=((a10 && pi1)||(a11 && pi0));
    bool a1e=((a10 && e1)||(a11 && e0));
    bool a1mu=((a10 && mu1)||(a11 && mu0));
    bool rhorho=(rho0 && rho1);
    bool rhopi=((rho0 && pi1)||(rho1 && pi0));
    bool rhoe=((rho0 && e1)||(rho1 && e0));
    bool rhomu=((rho0 && mu1)||(rho1 && mu0));
    bool pipi=(pi0 && pi1);
    bool pie=((pi0 && e1)||(pi1 && e0));
    bool pimu=((pi0 && mu1)||(pi1 && mu0));
    bool ee=(e0 && e1);
    bool emu=((e0 && mu1)||(e1 && mu0));
    bool mumu=(mu0 && mu1);

    bool a1pi0a1pi0=(a1pi00 && a1pi01);
    bool pi23pi0pi23pi0=(pi23pi00 && pi23pi01);
    
    bool a1a1pi0=((a10 && a1pi01) || (a11 && a1pi00));
    bool a1pi23pi0=((a10 && pi23pi01) || (a11 && pi23pi00));

    bool a1pi0pi23pi0=((a1pi00 && pi23pi01) || (a1pi01 && pi23pi00));

    bool a1pi0rho=((rho0 && a1pi01) || (rho1 && a1pi00));
    bool rhopi23pi0=((rho0 && pi23pi01) || (rho1 && pi23pi00));

    bool a1pi0pi=((pi0 && a1pi01) || (pi1 && a1pi00));
    bool pipi23pi0=((pi0 && pi23pi01) || (pi1 && pi23pi00));

    bool a1pi0e=((e0 && a1pi01) || (e1 && a1pi00));
    bool pi23pi0e=((e0 && pi23pi01) || (e1 && pi23pi00));

    bool a1pi0mu=((mu0 && a1pi01) || (mu1 && a1pi00));
    bool pi23pi0mu=((mu0 && pi23pi01) || (mu1 && pi23pi00));
    
    TLorentzVector testTruth(0,0,0,0);
    Pions1.clear();
    Pions2.clear();
    Pions1Charge.clear();
    Pions2Charge.clear();

    vector<TLorentzVector> HadPionsTruth_minus;
    vector<TLorentzVector> HadPionsTruth_plus;
    vector<double> HadPionsChargeTruth_minus;
    vector<double> HadPionsChargeTruth_plus;
    vector<TLorentzVector> tauandprodTruthminus;
    vector<TLorentzVector> tauandprodTruthplus;
		    
    vector<TLorentzVector> HadPions_minus;   
    vector<double> HadPionsCharge_minus;
    vector<TLorentzVector> HadPions_plus;    
    vector<double> HadPionsCharge_plus;
    
    vector<TLorentzVector> HadRefitPions_minus;   
    vector<double> HadRefitPionsCharge_minus;
    vector<TLorentzVector> HadRefitPions_plus;    
    vector<double> HadRefitPionsCharge_plus;
    
    if(a1a1 || a1a1MVA) 
      {
	HadPions_minus.push_back(Ntp->PFTau_PionsP4(Tauminus,0));
	HadPions_minus.push_back(Ntp->PFTau_PionsP4(Tauminus,1));
	HadPions_minus.push_back(Ntp->PFTau_PionsP4(Tauminus,2));
	
	HadPionsCharge_minus.push_back(Ntp->PFTau_PionsCharge(Tauminus, 0));
	HadPionsCharge_minus.push_back(Ntp->PFTau_PionsCharge(Tauminus, 1));
	HadPionsCharge_minus.push_back(Ntp->PFTau_PionsCharge(Tauminus, 2));

	SCalculator Scalcminusa1("a1");
	Scalcminusa1.SortPions(HadPions_minus,HadPionsCharge_minus);

	HadPions_plus.push_back(Ntp->PFTau_PionsP4(Tauplus,0));
	HadPions_plus.push_back(Ntp->PFTau_PionsP4(Tauplus,1));
	HadPions_plus.push_back(Ntp->PFTau_PionsP4(Tauplus,2));
	
	HadPionsCharge_plus.push_back(Ntp->PFTau_PionsCharge(Tauplus, 0));
	HadPionsCharge_plus.push_back(Ntp->PFTau_PionsCharge(Tauplus, 1));
	HadPionsCharge_plus.push_back(Ntp->PFTau_PionsCharge(Tauplus, 2));

	SCalculator Scalcplusa1("a1");
	Scalcplusa1.SortPions(HadPions_plus,HadPionsCharge_plus);

	

	HadRefitPions_minus.push_back(Ntp->PFTauRefit_PionsP4(Tauminus,0));
	HadRefitPions_minus.push_back(Ntp->PFTauRefit_PionsP4(Tauminus,1));
	HadRefitPions_minus.push_back(Ntp->PFTauRefit_PionsP4(Tauminus,2));
	
	HadRefitPionsCharge_minus.push_back(Ntp->PFTauRefit_PionsCharge(Tauminus, 0));
	HadRefitPionsCharge_minus.push_back(Ntp->PFTauRefit_PionsCharge(Tauminus, 1));
	HadRefitPionsCharge_minus.push_back(Ntp->PFTauRefit_PionsCharge(Tauminus, 2));

	SCalculator ScalcRefitminusa1("a1");
	ScalcRefitminusa1.SortPions(HadRefitPions_minus,HadRefitPionsCharge_minus);

	HadRefitPions_plus.push_back(Ntp->PFTauRefit_PionsP4(Tauplus,0));
	HadRefitPions_plus.push_back(Ntp->PFTauRefit_PionsP4(Tauplus,1));
	HadRefitPions_plus.push_back(Ntp->PFTauRefit_PionsP4(Tauplus,2));
	
	HadRefitPionsCharge_plus.push_back(Ntp->PFTauRefit_PionsCharge(Tauplus, 0));
	HadRefitPionsCharge_plus.push_back(Ntp->PFTauRefit_PionsCharge(Tauplus, 1));
	HadRefitPionsCharge_plus.push_back(Ntp->PFTauRefit_PionsCharge(Tauplus, 2));

	SCalculator ScalcRefitplusa1("a1");
	ScalcRefitplusa1.SortPions(HadRefitPions_plus,HadRefitPionsCharge_plus);
      }

    SCalculator Scalc("a1");
    SCalculator ScalcPVRefitNoBS("a1");
    SCalculator ScalcPVRefitBS("a1");
    SCalculator ScalcPVRefitNoBSZNominal("a1");
    SCalculator ScalcPVRefitBSZNominal("a1");

    TLorentzVector zeroLV(0,0,0,0);
    std::vector<TLorentzVector> VectZeroLV;
    VectZeroLV.push_back(zeroLV);
    VectZeroLV.push_back(zeroLV);
    VectZeroLV.push_back(zeroLV);

    
    if (std::isnan(Wspin)!=true)
      {
	if ((HadPions_minus!=HadPions_plus) && (HadPions_minus!=VectZeroLV) && (HadPions_plus!=VectZeroLV) && (HadRefitPions_minus!=HadRefitPions_plus) && (HadRefitPions_minus!=VectZeroLV) && (HadRefitPions_plus!=VectZeroLV)/* && (std::isnan(h1Norm)!=true) && (std::isnan(h2Norm)!=true) && (TauminusPairConstraint!=zeroLV) && ( TauplusPairConstraint!=zeroLV) && (h1.Cross(tauminus_HRF.Vect().Unit())).Mag()!=0 && (h2.Cross(tauplus_HRF.Vect().Unit())).Mag()!=0*/)
	  {
	    if(a1a1)
	      {
		if(TauminusPairConstraint!=TauplusPairConstraint && TauminusPairConstraint!=zeroLV && TauplusPairConstraint!=zeroLV )
		  {
		
		    SCalculator Scalc1test("a1");
		    SCalculator Scalc2test("a1");

		    vector<TLorentzVector> tauandprodminustest;
		    vector<TLorentzVector> tauandprodplustest;
		    
		    tauandprodminustest.push_back(TauminusPairConstraint);
		    tauandprodminustest.push_back(HadPions_minus.at(0));
		    tauandprodminustest.push_back(HadPions_minus.at(1));
		    tauandprodminustest.push_back(HadPions_minus.at(2));

		    tauandprodplustest.push_back(TauplusPairConstraint); 
		    tauandprodplustest.push_back(HadPions_plus.at(0));  
		    tauandprodplustest.push_back(HadPions_plus.at(1)); 
		    tauandprodplustest.push_back(HadPions_plus.at(2)); 

		    Scalc1test.Configure(tauandprodminustest,tauandprodminustest.at(0)+tauandprodplustest.at(0), -1);
		    TVector3 h1test=Scalc1test.pv();

		    Scalc2test.Configure(tauandprodplustest,tauandprodminustest.at(0)+tauandprodplustest.at(0), +1);
		    TVector3 h2test=Scalc2test.pv();
		    
		    
		    TLorentzVector tauminustest_HRF = Scalc1test.Boost(tauandprodminustest.at(0),tauandprodminustest.at(0)+tauandprodplustest.at(0));
		    TLorentzVector tauplustest_HRF  = Scalc2test.Boost(tauandprodplustest.at(0),tauandprodminustest.at(0)+tauandprodplustest.at(0));

		    if(tauminustest_HRF==zeroLV){cout<<endl;tauminustest_HRF.Print();cout<<endl;}
		    if(tauplustest_HRF==zeroLV){cout<<endl;tauplustest_HRF.Print();cout<<endl;}
		    
		    
		    double h1Normtest=1./h1test.Mag();
		    double h2Normtest=1./h2test.Mag();

		    if(std::isnan(h1Normtest)!=true && std::isnan(h2Normtest)!=true)
		      {
			h1test=h1test*h1Normtest;
			h2test=h2test*h2Normtest;
		    
			double k1Normtest=1./((h1test.Cross(tauminustest_HRF.Vect().Unit())).Mag());
			double k2Normtest=1./((h2test.Cross(tauplustest_HRF.Vect().Unit())).Mag());
			TVector3 k1test = (h1test.Cross(tauminustest_HRF.Vect().Unit()))*k1Normtest;
			TVector3 k2test = (h2test.Cross(tauplustest_HRF.Vect().Unit()))*k2Normtest;

		
			if(((h1test.Cross(h2test))*(tauminustest_HRF.Vect().Unit()))<=0) test.at(t).Fill(TMath::ATan2((k1test.Cross(k2test)).Mag(),k1test*k2test),Wspin);
			else test.at(t).Fill(2.*TMath::Pi()-TMath::ATan2((k1test.Cross(k2test)).Mag(),k1test*k2test),Wspin);

			
			polarimetricAcopAngle.at(t).Fill(Scalc.AcopAngle("a1", "a1", TauminusPairConstraint, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraint , HadPions_plus, HadPionsCharge_plus),Wspin);

			if(a1a1)
			  {

			    if(a1a1_ ||a1rho  ||a1pi  || a1e || a1mu || rhorho || rhopi || rhoe || rhomu ||pipi || pie||pimu || ee||emu ||mumu || a1pi0a1pi0||pi23pi0pi23pi0 ||a1a1pi0||a1pi0pi23pi0 ||a1pi23pi0 ||a1pi0rho ||rhopi23pi0 ||a1pi0pi ||pipi23pi0 || a1pi0e||pi23pi0e ||a1pi0mu ||pi23pi0mu)purityDM=true;
	    
			    if(!purityDM) { PurityDM.at(t).Fill(28.,w); //other
			      // cout<<endl;
			      // cout<<" pdgid0: ";
			      // for(int i=0;i<Ntp->NMCTauDecayProducts(0);i++)
			      //   {
			      //     cout<<Ntp->MCTauandProd_pdgid(0,i)<<"  ";
		  
			      //   }
			      // cout<<endl;
			      // cout<<" pdgid1: ";
			      // for(int i=0;i<Ntp->NMCTauDecayProducts(1);i++)
			      //   {
			      //     cout<<Ntp->MCTauandProd_pdgid(1,i)<<"  ";
		  
			      //   }
			      // cout<<endl;
			      //cout<<a1a1_ <<a1rho  <<a1pi  << a1e << a1mu << rhorho << rhopi << rhoe << rhomu <<pipi << pie<<pimu << ee<<emu <<mumu<<endl;
			      //cout<<endl;
			    }
			    if(purityDM)
			      {
				if (a1a1_){PurityDM.at(t).Fill(0.,w);
				  // cout<<endl;
				  // cout<<" a1a1 pdgid0: ";
				  // for(int i=0;i<Ntp->NMCTauDecayProducts(0);i++)
				  // 	{
				  // 	  cout<<Ntp->MCTauandProd_pdgid(0,i)<<"  ";
		  
				  // 	}
				  // cout<<endl;
				  // cout<<" a1a1 pdgid1: ";
				  // for(int i=0;i<Ntp->NMCTauDecayProducts(1);i++)
				  // 	{
				  // 	  cout<<Ntp->MCTauandProd_pdgid(1,i)<<"  ";
		  
				  // 	}
				  // cout<<endl;
				}
				    
				else if(a1a1pi0)PurityDM.at(t).Fill(1.,w);
				else if (a1rho){PurityDM.at(t).Fill(2.,w);
				  // cout<<endl;
				  // for(int i=0;i<Ntp->NMCTauDecayProducts(0);i++)
				  // 	{
				  // 	  cout<<" a1rho pdgid0: "<<Ntp->MCTauandProd_pdgid(0,i);
		  
				  // 	}
				  // cout<<endl;
				  // for(int i=0;i<Ntp->NMCTauDecayProducts(1);i++)
				  // 	{
				  // 	  cout<<" a1rho pdgid1: "<<Ntp->MCTauandProd_pdgid(1,i);
		  
				  // 	}
				  // cout<<endl;
				}
				else if (a1pi)PurityDM.at(t).Fill(3.,w);
				else if (a1pi23pi0)PurityDM.at(t).Fill(4.,w);
				else if (a1e)PurityDM.at(t).Fill(5.,w);
				else if (a1mu)PurityDM.at(t).Fill(6.,w);

				else if (a1pi0a1pi0)PurityDM.at(t).Fill(7.,w);
				else if (a1pi0rho)PurityDM.at(t).Fill(8.,w);
				else if (a1pi0pi)PurityDM.at(t).Fill(9.,w);
				else if (a1pi0pi23pi0)PurityDM.at(t).Fill(10.,w);
				else if (a1pi0e)PurityDM.at(t).Fill(11.,w);
				else if (a1pi0mu)PurityDM.at(t).Fill(12.,w);
				    
				else if (rhorho)PurityDM.at(t).Fill(13.,w);
				else if (rhopi)PurityDM.at(t).Fill(14.,w);
				else if (rhopi23pi0)PurityDM.at(t).Fill(15.,w);
				else if (rhoe)PurityDM.at(t).Fill(16.,w);
				else if (rhomu)PurityDM.at(t).Fill(17.,w);
				    
				else if (pipi)PurityDM.at(t).Fill(18.,w);
				else if (pipi23pi0)PurityDM.at(t).Fill(19.,w);
				else if (pie)PurityDM.at(t).Fill(20.,w);
				else if (pimu)PurityDM.at(t).Fill(21.,w);
				else if (pi23pi0pi23pi0)PurityDM.at(t).Fill(22.,w);
				else if (pi23pi0e)PurityDM.at(t).Fill(23.,w);
				else if (pi23pi0mu)PurityDM.at(t).Fill(24.,w);
				else if (ee)PurityDM.at(t).Fill(25.,w);
				else if (emu)PurityDM.at(t).Fill(26.,w);
				else if (mumu)PurityDM.at(t).Fill(27.,w);
			      }
			    //}
			  }

		      }
		  }
		if(TauminusPairConstraintNoBS!=zeroLV && TauplusPairConstraintNoBS!=zeroLV && ScalcPVRefitNoBS.isOk("a1", "a1", TauminusPairConstraintNoBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBS, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitNoBS.at(t).Fill(ScalcPVRefitNoBS.AcopAngle("a1", "a1", TauminusPairConstraintNoBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBS, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(TauminusPairConstraintBS!=zeroLV && TauplusPairConstraintBS!=zeroLV && ScalcPVRefitBS.isOk("a1", "a1", TauminusPairConstraintBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBS, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitBS.at(t).Fill(ScalcPVRefitBS.AcopAngle("a1", "a1", TauminusPairConstraintBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBS, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(TauminusPairConstraintNoBSZNominal!=zeroLV && TauplusPairConstraintNoBSZNominal!=zeroLV && ScalcPVRefitNoBSZNominal.isOk("a1", "a1", TauminusPairConstraintNoBSZNominal, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominal, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitNoBSZNominal.at(t).Fill(ScalcPVRefitNoBSZNominal.AcopAngle("a1", "a1", TauminusPairConstraintNoBSZNominal, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominal, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(TauminusPairConstraintBSZNominal!=zeroLV && TauplusPairConstraintBSZNominal!=zeroLV && ScalcPVRefitBSZNominal.isOk("a1", "a1", TauminusPairConstraintBSZNominal, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominal, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitBSZNominal.at(t).Fill(ScalcPVRefitBSZNominal.AcopAngle("a1", "a1", TauminusPairConstraintBSZNominal, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominal, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);
		
		
		if(isRefitNoBS && TauminusPairConstraintNoBS!=zeroLV && TauplusPairConstraintNoBS!=zeroLV && ScalcPVRefitNoBS.isOk("a1", "a1", TauminusPairConstraintNoBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBS, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitOnlyNoBS.at(t).Fill(ScalcPVRefitNoBS.AcopAngle("a1", "a1", TauminusPairConstraintNoBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBS, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(isRefitBS && TauminusPairConstraintBS!=zeroLV && TauplusPairConstraintBS!=zeroLV && ScalcPVRefitBS.isOk("a1", "a1", TauminusPairConstraintBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBS, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitOnlyBS.at(t).Fill(ScalcPVRefitBS.AcopAngle("a1", "a1", TauminusPairConstraintBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBS, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(isRefitNoBSZNominal && TauminusPairConstraintNoBSZNominal!=zeroLV && TauplusPairConstraintNoBSZNominal!=zeroLV && ScalcPVRefitNoBSZNominal.isOk("a1", "a1", TauminusPairConstraintNoBSZNominal, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominal, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitOnlyNoBSZNominal.at(t).Fill(ScalcPVRefitNoBSZNominal.AcopAngle("a1", "a1", TauminusPairConstraintNoBSZNominal, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominal, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(isRefitBSZNominal && TauminusPairConstraintBSZNominal!=zeroLV && TauplusPairConstraintBSZNominal!=zeroLV && ScalcPVRefitBSZNominal.isOk("a1", "a1", TauminusPairConstraintBSZNominal, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominal, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitOnlyBSZNominal.at(t).Fill(ScalcPVRefitBSZNominal.AcopAngle("a1", "a1", TauminusPairConstraintBSZNominal, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominal, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

	      }
	    if(a1a1MVA)
	      {	      
	    	if(TauminusPairConstraint!=zeroLV && TauplusPairConstraint!=zeroLV && Scalc.isOk("a1", "a1", TauminusPairConstraint, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraint , HadPions_plus, HadPionsCharge_plus)==true) polarimetricAcopAngleMVA.at(t).Fill(Scalc.AcopAngle("a1", "a1", TauminusPairConstraint, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraint, HadPions_plus, HadPionsCharge_plus),Wspin);
		
		if(TauminusPairConstraintNoBS!=zeroLV && TauplusPairConstraintNoBS!=zeroLV && ScalcPVRefitNoBS.isOk("a1", "a1", TauminusPairConstraintNoBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBS, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitNoBSNewMVA.at(t).Fill(ScalcPVRefitNoBS.AcopAngle("a1", "a1", TauminusPairConstraintNoBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBS, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(TauminusPairConstraintBS!=zeroLV && TauplusPairConstraintBS!=zeroLV && ScalcPVRefitBS.isOk("a1", "a1", TauminusPairConstraintBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBS, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitBSNewMVA.at(t).Fill(ScalcPVRefitBS.AcopAngle("a1", "a1", TauminusPairConstraintBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBS, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(TauminusPairConstraintNoBSZNominal!=zeroLV && TauplusPairConstraintNoBSZNominal!=zeroLV && ScalcPVRefitNoBSZNominal.isOk("a1", "a1", TauminusPairConstraintNoBSZNominal, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominal, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitNoBSZNominalNewMVA.at(t).Fill(ScalcPVRefitNoBSZNominal.AcopAngle("a1", "a1", TauminusPairConstraintNoBSZNominal, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominal, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(TauminusPairConstraintBSZNominal!=zeroLV && TauplusPairConstraintBSZNominal!=zeroLV && ScalcPVRefitBSZNominal.isOk("a1", "a1", TauminusPairConstraintBSZNominal, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominal, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitBSZNominalNewMVA.at(t).Fill(ScalcPVRefitBSZNominal.AcopAngle("a1", "a1", TauminusPairConstraintBSZNominal, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominal, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);
		
		if(isRefitNoBS && TauminusPairConstraintNoBS!=zeroLV && TauplusPairConstraintNoBS!=zeroLV && ScalcPVRefitNoBS.isOk("a1", "a1", TauminusPairConstraintNoBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBS, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitOnlyNoBSNewMVA.at(t).Fill(ScalcPVRefitNoBS.AcopAngle("a1", "a1", TauminusPairConstraintNoBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBS, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(isRefitBS && TauminusPairConstraintBS!=zeroLV && TauplusPairConstraintBS!=zeroLV && ScalcPVRefitBS.isOk("a1", "a1", TauminusPairConstraintBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBS, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitOnlyBSNewMVA.at(t).Fill(ScalcPVRefitBS.AcopAngle("a1", "a1", TauminusPairConstraintBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBS, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(isRefitNoBSZNominal && TauminusPairConstraintNoBSZNominal!=zeroLV && TauplusPairConstraintNoBSZNominal!=zeroLV && ScalcPVRefitNoBSZNominal.isOk("a1", "a1", TauminusPairConstraintNoBSZNominal, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominal, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitOnlyNoBSZNominalNewMVA.at(t).Fill(ScalcPVRefitNoBSZNominal.AcopAngle("a1", "a1", TauminusPairConstraintNoBSZNominal, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominal, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);
		if(isRefitBSZNominal && TauminusPairConstraintBSZNominal!=zeroLV && TauplusPairConstraintBSZNominal!=zeroLV && ScalcPVRefitBSZNominal.isOk("a1", "a1", TauminusPairConstraintBSZNominal, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominal, HadRefitPions_plus, HadRefitPionsCharge_plus)==true){ polarimetricAcopAnglePVRefitOnlyBSZNominalNewMVA.at(t).Fill(ScalcPVRefitBSZNominal.AcopAngle("a1", "a1", TauminusPairConstraintBSZNominal, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominal, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);}
		if(!Ntp->isData()){
		  if(a1a1_ ||a1rho  ||a1pi  || a1e || a1mu || rhorho || rhopi || rhoe || rhomu ||pipi || pie||pimu || ee||emu ||mumu || a1pi0a1pi0||pi23pi0pi23pi0 ||a1a1pi0||a1pi0pi23pi0 ||a1pi23pi0 ||a1pi0rho ||rhopi23pi0 ||a1pi0pi ||pipi23pi0 || a1pi0e||pi23pi0e ||a1pi0mu ||pi23pi0mu)purityNewMVA=true;

		  if(!purityNewMVA) { PurityNewMVA.at(t).Fill(28.,w); //other
		    // cout<<endl;
		    // cout<<" pdgid0: ";
		    // for(int i=0;i<Ntp->NMCTauDecayProducts(0);i++)
		    //   {
		    //     cout<<Ntp->MCTauandProd_pdgid(0,i)<<"  ";
		  
		    //   }
		    // cout<<endl;
		    // cout<<" pdgid1: ";
		    // for(int i=0;i<Ntp->NMCTauDecayProducts(1);i++)
		    //   {
		    //     cout<<Ntp->MCTauandProd_pdgid(1,i)<<"  ";
		  
		    //   }
		    // cout<<endl;
		    //cout<<a1a1_ <<a1rho  <<a1pi  << a1e << a1mu << rhorho << rhopi << rhoe << rhomu <<pipi << pie<<pimu << ee<<emu <<mumu<<endl;
		    //cout<<endl;
		  }
		  if(purityNewMVA)
		    {
		      if (a1a1_){PurityNewMVA.at(t).Fill(0.,w);
			// cout<<endl;
			// cout<<" a1a1 pdgid0: ";
			// for(int i=0;i<Ntp->NMCTauDecayProducts(0);i++)
			// 	{
			// 	  cout<<Ntp->MCTauandProd_pdgid(0,i)<<"  ";
		  
			// 	}
			// cout<<endl;
			// cout<<" a1a1 pdgid1: ";
			// for(int i=0;i<Ntp->NMCTauDecayProducts(1);i++)
			// 	{
			// 	  cout<<Ntp->MCTauandProd_pdgid(1,i)<<"  ";
		  
			// 	}
			// cout<<endl;
		      }
				    
		      else if(a1a1pi0)PurityNewMVA.at(t).Fill(1.,w);
		      else if (a1rho){PurityNewMVA.at(t).Fill(2.,w);
			// cout<<endl;
			// for(int i=0;i<Ntp->NMCTauDecayProducts(0);i++)
			// 	{
			// 	  cout<<" a1rho pdgid0: "<<Ntp->MCTauandProd_pdgid(0,i);
		  
			// 	}
			// cout<<endl;
			// for(int i=0;i<Ntp->NMCTauDecayProducts(1);i++)
			// 	{
			// 	  cout<<" a1rho pdgid1: "<<Ntp->MCTauandProd_pdgid(1,i);
		  
			// 	}
			// cout<<endl;
		      }
		      else if (a1pi)PurityNewMVA.at(t).Fill(3.,w);
		      else if (a1pi23pi0)PurityNewMVA.at(t).Fill(4.,w);
		      else if (a1e)PurityNewMVA.at(t).Fill(5.,w);
		      else if (a1mu)PurityNewMVA.at(t).Fill(6.,w);

		      else if (a1pi0a1pi0)PurityNewMVA.at(t).Fill(7.,w);
		      else if (a1pi0rho)PurityNewMVA.at(t).Fill(8.,w);
		      else if (a1pi0pi)PurityNewMVA.at(t).Fill(9.,w);
		      else if (a1pi0pi23pi0)PurityNewMVA.at(t).Fill(10.,w);
		      else if (a1pi0e)PurityNewMVA.at(t).Fill(11.,w);
		      else if (a1pi0mu)PurityNewMVA.at(t).Fill(12.,w);
				    
		      else if (rhorho)PurityNewMVA.at(t).Fill(13.,w);
		      else if (rhopi)PurityNewMVA.at(t).Fill(14.,w);
		      else if (rhopi23pi0)PurityNewMVA.at(t).Fill(15.,w);
		      else if (rhoe)PurityNewMVA.at(t).Fill(16.,w);
		      else if (rhomu)PurityNewMVA.at(t).Fill(17.,w);
				    
		      else if (pipi)PurityNewMVA.at(t).Fill(18.,w);
		      else if (pipi23pi0)PurityNewMVA.at(t).Fill(19.,w);
		      else if (pie)PurityNewMVA.at(t).Fill(20.,w);
		      else if (pimu)PurityNewMVA.at(t).Fill(21.,w);
		      else if (pi23pi0pi23pi0)PurityNewMVA.at(t).Fill(22.,w);
		      else if (pi23pi0e)PurityNewMVA.at(t).Fill(23.,w);
		      else if (pi23pi0mu)PurityNewMVA.at(t).Fill(24.,w);
		      else if (ee)PurityNewMVA.at(t).Fill(25.,w);
		      else if (emu)PurityNewMVA.at(t).Fill(26.,w);
		      else if (mumu)PurityNewMVA.at(t).Fill(27.,w);
		    }
		}
		//}
	      }

	  }
	if(!Ntp->isData())
	  {
	    if(a1a1_ && Ntp->CheckDecayID(5,5))
	      {
		Pions1=Ntp->GetTruthPionsFromA1(0);
		Pions1Charge.push_back(1);
		Pions1Charge.push_back(-1);
		Pions1Charge.push_back(-1);
		Pions2=Ntp->GetTruthPionsFromA1(1);
		Pions2Charge.push_back(-1);
		Pions2Charge.push_back(1);
		Pions2Charge.push_back(1);
		SCalculator Scalc1Truth("a1");
		SCalculator Scalc2Truth("a1");
		Scalc1Truth.SortPions(Pions1, Pions1Charge);
		Scalc2Truth.SortPions(Pions2, Pions2Charge);
		if(Ntp->MCTau_pdgid(0)==15)
		  {
		    tauandprodTruthminus.push_back(Ntp->GetTruthTauLV(5,0));
		    tauandprodTruthminus.push_back(Pions1.at(0));
		    tauandprodTruthminus.push_back(Pions1.at(1));
		    tauandprodTruthminus.push_back(Pions1.at(2));
		    tauandprodTruthplus.push_back(Ntp->GetTruthTauLV(5,1));   
		    tauandprodTruthplus.push_back(Pions2.at(0));   
		    tauandprodTruthplus.push_back(Pions2.at(1));   
		    tauandprodTruthplus.push_back(Pions2.at(2));   
		  }
		else if (Ntp->MCTau_pdgid(0)==-15)
		  {
		    tauandprodTruthminus.push_back(Ntp->GetTruthTauLV(5,1));
		    tauandprodTruthminus.push_back(Pions2.at(0));
		    tauandprodTruthminus.push_back(Pions2.at(1));
		    tauandprodTruthminus.push_back(Pions2.at(2));
		    tauandprodTruthplus.push_back(Ntp->GetTruthTauLV(5,0));
		    tauandprodTruthplus.push_back(Pions1.at(0));   
		    tauandprodTruthplus.push_back(Pions1.at(1));  
		    tauandprodTruthplus.push_back(Pions1.at(2));   
		  }	
	      }
	    if((Pions1!=Pions2) && (Pions1!=VectZeroLV) && (Pions2!=VectZeroLV) && tauandprodTruthminus.at(0)!=zeroLV &&tauandprodTruthplus.at(0)!=zeroLV &&tauandprodTruthminus.at(0)!=tauandprodTruthplus.at(0))
	      {
		if(a1a1_  && Ntp->CheckDecayID(5,5))
		  {
		    
		    SCalculator Scalc1Truth("a1");
		    SCalculator Scalc2Truth("a1");
		    
		    Scalc1Truth.Configure(tauandprodTruthminus,tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0), -1);
		    TVector3 h1Truth=Scalc1Truth.pv();
			
		    Scalc2Truth.Configure(tauandprodTruthplus,tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0), +1);
		    TVector3 h2Truth=Scalc2Truth.pv();
			
		    double h1TruthNorm=1./h1Truth.Mag();
		    double h2TruthNorm=1./h2Truth.Mag();
		    
		    if(std::isnan(h1TruthNorm)!=true && std::isnan(h2TruthNorm)!=true)
		      { 
			TLorentzVector tauminusTruth_HRF = Scalc1Truth.Boost(tauandprodTruthminus.at(0),tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0));
			TLorentzVector tauplusTruth_HRF  = Scalc2Truth.Boost(tauandprodTruthplus.at(0),tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0));
			
			double norm1Truth=1./(((h1Truth*h1TruthNorm).Cross(tauminusTruth_HRF.Vect().Unit())).Mag());
			double norm2Truth=1./(((h2Truth*h2TruthNorm).Cross(tauplusTruth_HRF.Vect().Unit())).Mag());
			TVector3 k1Truth = ((h1Truth*h1TruthNorm).Cross(tauminusTruth_HRF.Vect().Unit()))*norm1Truth;
			TVector3 k2Truth = ((h2Truth*h2TruthNorm).Cross(tauplusTruth_HRF.Vect().Unit()))*norm2Truth;
			
			if(((h1Truth*h1TruthNorm).Cross(h2Truth*h2TruthNorm))*tauminusTruth_HRF.Vect().Unit()<=0) {polarimetricAcopAngleTruthA1.at(t).Fill(TMath::ATan2((k1Truth.Cross(k2Truth)).Mag(),k1Truth*k2Truth),Wspin);/*cout<<"Angle Truth: "<<TMath::ATan2((k1Truth.Cross(k2Truth)).Mag(),k1Truth*k2Truth)<<endl;*/}
			else{ polarimetricAcopAngleTruthA1.at(t).Fill(2*TMath::Pi()-TMath::ATan2((k1Truth.Cross(k2Truth)).Mag(),k1Truth*k2Truth),Wspin);/*cout<<"Angle Truth: "<<2*TMath::Pi()-TMath::ATan2((k1Truth.Cross(k2Truth)).Mag(),k1Truth*k2Truth)<<endl;*/}
				
		      }
		  }
	      }
	    if(a1a1_  && Ntp->CheckDecayID(5,5))
	      {
		PVXResol.at(t).Fill((Ntp->PVtx().X()-Ntp->PVtx_Gen().X())/Ntp->PVtx_Gen().X(),w);
		PVXNoBSResol.at(t).Fill((tauNoBSPrimaryVertex.X()-Ntp->PVtx_Gen().X())/Ntp->PVtx_Gen().X(),w);
		PVXBSResol.at(t).Fill((tauBSPrimaryVertex.X()-Ntp->PVtx_Gen().X())/Ntp->PVtx_Gen().X(),w);
		PVYResol.at(t).Fill((Ntp->PVtx().Y()-Ntp->PVtx_Gen().Y())/Ntp->PVtx_Gen().Y(),w);
		PVYNoBSResol.at(t).Fill((tauNoBSPrimaryVertex.Y()-Ntp->PVtx_Gen().Y())/Ntp->PVtx_Gen().Y(),w);
		PVYBSResol.at(t).Fill((tauBSPrimaryVertex.Y()-Ntp->PVtx_Gen().Y())/Ntp->PVtx_Gen().Y(),w);
		PVZResol.at(t).Fill((Ntp->PVtx().Z()-Ntp->PVtx_Gen().Z())/Ntp->PVtx_Gen().Z(),w);
		PVZNoBSResol.at(t).Fill((tauNoBSPrimaryVertex.Z()-Ntp->PVtx_Gen().Z())/Ntp->PVtx_Gen().Z(),w);
		PVZBSResol.at(t).Fill((tauBSPrimaryVertex.Z()-Ntp->PVtx_Gen().Z())/Ntp->PVtx_Gen().Z(),w);
		if(isRefitNoBS)PVXNoBSOnlyResol.at(t).Fill((tauNoBSPrimaryVertex.X()-Ntp->PVtx_Gen().X())/Ntp->PVtx_Gen().X(),w);
		if(isRefitBS)PVXBSOnlyResol.at(t).Fill((tauBSPrimaryVertex.X()-Ntp->PVtx_Gen().X())/Ntp->PVtx_Gen().X(),w);
		if(isRefitNoBS)PVYNoBSOnlyResol.at(t).Fill((tauNoBSPrimaryVertex.Y()-Ntp->PVtx_Gen().Y())/Ntp->PVtx_Gen().Y(),w);
		if(isRefitBS)PVYBSOnlyResol.at(t).Fill((tauBSPrimaryVertex.Y()-Ntp->PVtx_Gen().Y())/Ntp->PVtx_Gen().Y(),w);
		if(isRefitNoBS)PVZNoBSOnlyResol.at(t).Fill((tauNoBSPrimaryVertex.Z()-Ntp->PVtx_Gen().Z())/Ntp->PVtx_Gen().Z(),w);
		if(isRefitBS)PVZBSOnlyResol.at(t).Fill((tauBSPrimaryVertex.Z()-Ntp->PVtx_Gen().Z())/Ntp->PVtx_Gen().Z(),w);
		if(Ntp->MCTau_pdgid(0)==15 && Ntp->PFTau_secondaryVertex_pos_Size()>1 && Ntp->MCTauandProd_VertexSize()==2)//nplus==1 && nminus==2)
		  {
		    ResolPullXVtxIna1a1.at(t).Fill((TauminusSecondaryVertex.X()-Ntp->MCTauandProd_Vertex(0,1).X())/Ntp->MCTauandProd_Vertex(0,1).X());
		    ResolPullXVtxIna1a1.at(t).Fill((TauplusSecondaryVertex.X()-Ntp->MCTauandProd_Vertex(1,1).X())/Ntp->MCTauandProd_Vertex(1,1).X());
		    ResolPullYVtxIna1a1.at(t).Fill((TauminusSecondaryVertex.Y()-Ntp->MCTauandProd_Vertex(0,1).Y())/Ntp->MCTauandProd_Vertex(0,1).Y());
		    ResolPullYVtxIna1a1.at(t).Fill((TauplusSecondaryVertex.Y()-Ntp->MCTauandProd_Vertex(1,1).Y())/Ntp->MCTauandProd_Vertex(1,1).Y());
		    ResolPullZVtxIna1a1.at(t).Fill((TauminusSecondaryVertex.Z()-Ntp->MCTauandProd_Vertex(0,1).Z())/Ntp->MCTauandProd_Vertex(0,1).Z());
		    ResolPullZVtxIna1a1.at(t).Fill((TauplusSecondaryVertex.Z()-Ntp->MCTauandProd_Vertex(1,1).Z())/Ntp->MCTauandProd_Vertex(1,1).Z());
		  }
		else if(Ntp->MCTau_pdgid(0)==-15 && Ntp->PFTau_secondaryVertex_pos_Size()>1 && Ntp->MCTauandProd_VertexSize()==2)
		  {
		    ResolPullXVtxIna1a1.at(t).Fill((TauminusSecondaryVertex.X()-Ntp->MCTauandProd_Vertex(1,1).X())/Ntp->MCTauandProd_Vertex(1,1).X());
		    ResolPullXVtxIna1a1.at(t).Fill((TauplusSecondaryVertex.X()-Ntp->MCTauandProd_Vertex(0,1).X())/Ntp->MCTauandProd_Vertex(0,1).X());
		    ResolPullYVtxIna1a1.at(t).Fill((TauminusSecondaryVertex.Y()-Ntp->MCTauandProd_Vertex(1,1).Y())/Ntp->MCTauandProd_Vertex(1,1).Y());
		    ResolPullYVtxIna1a1.at(t).Fill((TauplusSecondaryVertex.Y()-Ntp->MCTauandProd_Vertex(0,1).Y())/Ntp->MCTauandProd_Vertex(0,1).Y());
		    ResolPullZVtxIna1a1.at(t).Fill((TauminusSecondaryVertex.Z()-Ntp->MCTauandProd_Vertex(1,1).Z())/Ntp->MCTauandProd_Vertex(1,1).Z());
		    ResolPullZVtxIna1a1.at(t).Fill((TauplusSecondaryVertex.Z()-Ntp->MCTauandProd_Vertex(0,1).Z())/Ntp->MCTauandProd_Vertex(0,1).Z());
		    
		  }
	      }
	  }
      }
  }
}
//}
//}
//  This is a function if you want to do something after the event loop
void HCPTauTau::Finish() {

  if(mode == RECONSTRUCT) {
    // std::cout<<" Starting Finish!  " <<std::endl;
    // std::cout<<"A  Data  "<< NQCD.at(0).GetBinContent(1) << std::endl;
    // std::cout<<"B  Data  "<< NQCD.at(0).GetBinContent(2) << std::endl;
    // std::cout<<"C  Data  "<< NQCD.at(0).GetBinContent(3) << std::endl;
    // std::cout<<"D  Data  "<< NQCD.at(0).GetBinContent(4) << std::endl;
    SkimConfig SC;
    SC.ApplySkimEfficiency(types,Npassed, Npassed_noweight);
    std::vector<double> QCD_Integral_B;
    double QCD_IntegralMC_B;
    double QCD_Integral_B_Data_minus_MC = 0;
    
    std::vector<double> QCD_Integral_C;
    double QCD_IntegralMC_C;
    double QCD_Integral_C_Data_minus_MC = 0;
    
    std::vector<double> QCD_Integral_D;
    double QCD_IntegralMC_D;
    double QCD_Integral_D_Data_minus_MC = 0;
    //Get Yields in ABCD for QCD Scalefactor      
    double NQCDShape=0;
    double NQCDData=0;
    double NQCDLead=0;
    //double NQCDSub=0;
    std::vector<double> QCDLead;
    //std::vector<double> QCDSub;

    for(unsigned i=0;i<CrossSectionandAcceptance.size();i++){
      QCD_Integral_B.push_back(NQCD.at(i).GetBinContent(2));
      QCD_Integral_C.push_back(NQCD.at(i).GetBinContent(3));
      QCD_Integral_D.push_back(NQCD.at(i).GetBinContent(4));

      QCDLead.push_back(NFFLeadMC.at(i).GetBinContent(1));
      //cout<<"QCDLead en i: "<<NFFLeadMC.at(i).GetBinContent(1)<<endl;

      //      QCDSub.push_back(NFFSubMC.at(i).GetBinContent(1));

	
      if(CrossSectionandAcceptance.at(i)>0){
	
    	QCD_Integral_B.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
    	QCD_Integral_C.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
    	QCD_Integral_D.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
	QCDLead.at(i) *=CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
	//cout<<"QCDLead en i normalisé: "<<QCDLead.at(i)<<endl;

	//	QCDSub.at(i) *=CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
      }
    }

    for(unsigned i=0;i<CrossSectionandAcceptance.size();i++){
      if(HConfig.GetID(i) == DataMCType::Data){
	NQCDData=NFFData.at(i).GetBinContent(1);
    	QCD_Integral_B_Data_minus_MC  += QCD_Integral_B.at(i);
    	QCD_Integral_C_Data_minus_MC += QCD_Integral_C.at(i);
    	QCD_Integral_D_Data_minus_MC += QCD_Integral_D.at(i);
      }
      if(HConfig.GetID(i) == DataMCType::QCD){
	//NQCDShape+=QCDShape.at(i).GetBinContent(2);
      }
      if(CrossSectionandAcceptance.at(i)>0){
    	QCD_IntegralMC_B  += QCD_Integral_B.at(i);
    	QCD_IntegralMC_C  += QCD_Integral_C.at(i);
    	QCD_IntegralMC_D  += QCD_Integral_D.at(i);
	NQCDLead += QCDLead.at(i);
	//cout<<"NQCDLead additioné: "<<NQCDLead<<endl;
	//	NQCDSub += QCDSub.at(i);
      }
    }
    double QCD=NQCDData-NQCDLead;
    //PUPPImetcorr.at(1).Scale(NQCDData / Npassed.at(0).GetBinContent(0));
    //NQCDData=1000000.;
    //+NQCDSub;
    //double FF=NQCDShape/QCD;
    cout<<"Data: "<<NQCDData<<endl;
    cout<<"Lead: "<<NQCDLead<<endl;
    //cout<<"Sub: "<<NQCDSub<<endl;
    cout<<"QCD: "<<QCD<<endl;
    //cout<<"NQCDShape: "<<NQCDShape<<endl;
    //cout<<"FF: "<<FF<<endl;
    //double QCD_Signal=
    //double CDFactor = (QCD_Integral_C_Data_minus_MC  - QCD_IntegralMC_C )/ (QCD_Integral_D_Data_minus_MC - QCD_IntegralMC_D);
    //double QCD_Signal = QCD_Integral_B_Data_minus_MC *CDFactor;

    // std::cout << "Factor: " << CDFactor << std::endl;
    // std::cout << "QCD_Signal: " << QCD_Signal << std::endl;
    // std::cout << "QCD in B region "<<  QCD_Integral_B_Data_minus_MC <<std::endl;
    // std::cout << "QCD_Integral_B_Data_minus_MC is: " << QCD_Integral_B_Data_minus_MC << std::endl;
    // std::cout << "QCD_Integral_C_Data_minus_MC is: " << QCD_Integral_C_Data_minus_MC << std::endl;
    // std::cout << "QCD_Integral_D_Data_minus_MC is: " << QCD_Integral_D_Data_minus_MC << std::endl;
    // std::cout << "QCD_IntegralMC_B is: " << QCD_IntegralMC_B << std::endl;
    // std::cout << "QCD_IntegralMC_C is: " << QCD_IntegralMC_C << std::endl;
    // std::cout << "QCD_IntegralMC_D is: " << QCD_IntegralMC_D << std::endl;

    // ScaleAllHistOfType(HConfig.GetType(DataMCType::QCD),QCD_Signal/Nminus0.at(0).at(HConfig.GetType(DataMCType::QCD)).Integral());
    // ScaleAllHistOfType(HConfig.GetType(DataMCType::QCD),(FF*(QCD_Integral_B_Data_minus_MC  - QCD_IntegralMC_B))/Nminus0.at(0).at(HConfig.GetType(DataMCType::QCD)).Integral());
    // ScaleAllHistOfType(HConfig.GetType(DataMCType::QCD),QCD/Nminus0.at(0).at(HConfig.GetType(DataMCType::QCD)).Integral());
    double norm=1.;
    for(unsigned i=0;i<CrossSectionandAcceptance.size();i++){
      if(CrossSectionandAcceptance.at(i)>0){
    	norm= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
    	PUPPImetcorr.at(1).Add(&PUPPImetcorrQCDMC.at(i),-norm);
    	Tau1PT.at(1).Add(&Tau1PTQCDMC.at(i),-norm);
    	Tau1E.at(1).Add(&Tau1EQCDMC.at(i),-norm);
    	Tau1Mass.at(1).Add(&Tau1MassQCDMC.at(i),-norm);
    	Tau1Phi.at(1).Add(&Tau1PhiQCDMC.at(i),-norm);
    	Tau1Eta.at(1).Add(&Tau1EtaQCDMC.at(i),-norm);
    	Tau1dz.at(1).Add(&Tau1dzQCDMC.at(i),-norm);
    	Tau1HPSDecayMode.at(1).Add(&Tau1HPSDecayModeQCDMC.at(i),-norm);
    	Tau1MVADecayMode.at(1).Add(&Tau1MVADecayModeQCDMC.at(i),-norm);
	
    	Tau2PT.at(1).Add(&Tau2PTQCDMC.at(i),-norm);
    	Tau2E.at(1).Add(&Tau2EQCDMC.at(i),-norm);
    	Tau2Mass.at(1).Add(&Tau2MassQCDMC.at(i),-norm);
    	Tau2Phi.at(1).Add(&Tau2PhiQCDMC.at(i),-norm);
    	Tau2Eta.at(1).Add(&Tau2EtaQCDMC.at(i),-norm);
    	Tau2dz.at(1).Add(&Tau2dzQCDMC.at(i),-norm);
    	Tau2HPSDecayMode.at(1).Add(&Tau2HPSDecayModeQCDMC.at(i),-norm);
    	Tau2MVADecayMode.at(1).Add(&Tau2MVADecayModeQCDMC.at(i),-norm);

	NbJets.at(1).Add(&NbJetsQCDMC.at(i),-norm);
	TauTauVisMass.at(1).Add(&TauTauVisMassQCDMC.at(i),-norm);
	TauTauVisPT.at(1).Add(&TauTauVisPTQCDMC.at(i),-norm);
	Mjj.at(1).Add(&MjjQCDMC.at(i),-norm);
      }
    }
  }
  Selection::Finish();
}
