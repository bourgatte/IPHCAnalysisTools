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
#include "TauAnalysis/ClassicSVfit/interface/FastMTT.h"
#include "TauPolSoftware/TauDecaysInterface/interface/fonction_a1.h"
#include "TauPolSoftware/TauDecaysInterface/interface/SCalculator.h"
#include <algorithm>

HCPTauTau::HCPTauTau(TString Name_, TString id_):
  Selection(Name_,id_)
  //DataMC_Corr(true,true,false),
  //tauTrgSF("tight")
{
  ChargeSumDummy = -999;
  selMuon_IsoDummy = 999.;
  WorkSpaceFF2016=TFile::Open(((std::string)std::getenv("workdir")+"Code/fake_factors_tt_dRcorr/fakefactors_ws_tt_lite_2016_dR_corr.root").c_str(), "READ");
  wFF2016= (RooWorkspace*)gDirectory->Get("w");
  WorkSpaceFF2016->Close();
  WorkSpaceFF2017=TFile::Open(((std::string)std::getenv("workdir")+"Code/fake_factors_tt_dRcorr/fakefactors_ws_tt_lite_2017_dR_corr.root").c_str(), "READ");
  wFF2017= (RooWorkspace*)gDirectory->Get("w");
  WorkSpaceFF2017->Close();
  WorkSpaceFF2018=TFile::Open(((std::string)std::getenv("workdir")+"Code/fake_factors_tt_dRcorr/fakefactors_ws_tt_lite_2018_dR_corr.root").c_str(), "READ");
  wFF2018= (RooWorkspace*)gDirectory->Get("w");
  WorkSpaceFF2018->Close();
  BDT=new BDTClassification();
  BDT->PreAnalysis();
}

HCPTauTau::~HCPTauTau(){
  for(unsigned int j=0; j<Npassed.size(); j++){
    Logger(Logger::Info) << "Selection Summary before: "
			 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
			 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
  delete wFF2016;
  delete wFF2017;
  delete wFF2018;
}

void  HCPTauTau::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==Trigger)             cut.at(Trigger)=1;
    if(i==Id_and_Kin)            cut.at(Id_and_Kin)=1;
    //if(i==NPairsFound)         cut.at(NPairsFound)=1;
    //if(i==GoodIndex)           cut.at(GoodIndex)=1.;
    //if(i==ZTTMC)                 cut.at(ZTTMC)=1.;
    //if(i==METFilters)            cut.at(METFilters)=1.;
    //if(i==genmatch)              cut.at(genmatch)=1;
    if(i==TausIsolation)         cut.at(TausIsolation)=1;
    if(i==AgainstEleMu)          cut.at(AgainstEleMu)=1;
    //if(i==Tau2Isolation)       cut.at(Tau2Isolation)=1.;
    if(i==LeptonVeto)            cut.at(LeptonVeto)=0;
    if(i==PairCharge)            cut.at(PairCharge)=1.;
    if(i==PairMass)              cut.at(PairMass)=40.;
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
    if(i==Trigger){
      title.at(i)="Trigger Matching";
      hlabel="At least 1 good pair with Trig+Matching";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Trigger_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Trigger_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
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
    // if(i==ZTTMC){
    //   title.at(i)="ZTT MC";
    //   hlabel="ZTT MC";
    //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ZTTMC_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ZTTMC_",htitle,2,-0.5,1.5,hlabel,"Events"));
    // }
    // else if(i==METFilters){
    //   title.at(i)="MET Filters";
    //   hlabel="MET Filters";
    //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_METFilters_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_METFilters_",htitle,2,-0.5,1.5,hlabel,"Events"));
    // }
    else if(i==Id_and_Kin){
      title.at(i)="Id and Kinematic";
      hlabel="Number of Event with good particles";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Id_and_Kin_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Id_and_Kin_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    // else if(i==genmatch){
    //   title.at(i)="genmatch";
    //   hlabel="genmatch";
    //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_genmatch_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_genmatch_",htitle,2,-0.5,1.5,hlabel,"Events"));
    // }
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
    else if(i==PairMass){
      title.at(i)="Pair Visible Mass";
      hlabel="M(tau-tau)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PairMass_",htitle,30,0,150,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PairMass_",htitle,30,0,150,hlabel,"Events"));
    }

    /* else if(i==MTM){
       title.at(i)="Missing Transverse Mass";
       hlabel="Missing Transverse Mass";
       Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MTM_",htitle,30,0,100,hlabel,"Events"));
       Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MTM_",htitle,30,0,100,hlabel,"Events"));
       }*/
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  
  // Tau1PT=HConfig.GetTH1D(Name+"_Tau1PT","Transverse momentum of selected #tau1 candidate",20,40,140," P_{T}(#tau1), GeV","Events");
  // Tau1E=HConfig.GetTH1D(Name+"_Tau1E","Energy of selected #tau1 candidate",20,35.,140," E(#tau1), GeV","Events");
  // Tau1Mass=HConfig.GetTH1D(Name+"_Tau1Mass","Mass of selected #tau1 candidate",20,0,2.," M(#tau1), GeV","Events");
  // Tau1Phi=HConfig.GetTH1D(Name+"_Tau1Phi","Phi of selected #tau1 candidate",10,-3.5,3.5," #phi(#tau1)","Events");
  // Tau1Eta=HConfig.GetTH1D(Name+"_Tau1Eta","Pseudorapidity tau1",15,-2.7,2.7," #eta(#tau1)","Events");
  // Tau1dz=HConfig.GetTH1D(Name+"_Tau1dz","Tau1 dz",20,-0.02,0.02,"Tau1 dz","Events");
  // Tau1HPSDecayMode=HConfig.GetTH1D(Name+"_Tau1HPSDecayMode","Decay mode of the selected #tau candidate",12,-0.5,11.5,"HPS Mode","Events");
  // Tau1MVADecayMode=HConfig.GetTH1D(Name+"_Tau1MVADecayMode","MVA decay mode of the selected #tau1 candidate",13,-1.5,11.5,"MVA DM","Events");
  // Tau1GenMatch=HConfig.GetTH1D(Name+"_Tau1GenMatch","GenMatch of the selected #tau1 candidate",7,0,7,"GenMatch","Events");
  
  // Tau2PT=HConfig.GetTH1D(Name+"_Tau2PT","Transverse momentum of selected #tau2 candidate",12,40,100," P_{T}(#tau2), GeV","Events");
  // Tau2E=HConfig.GetTH1D(Name+"_Tau2E","Energy of selected #tau2 candidate",20,30.,140," E(#tau2), GeV","Events");
  // Tau2Mass=HConfig.GetTH1D(Name+"_Tau2Mass","Mass of selected #tau2 candidate",20,0,2.," M(#tau2), GeV","Events");
  // Tau2Phi=HConfig.GetTH1D(Name+"_Tau2Phi","Phi of selected #tau2 candidate",10,-3.5,3.5," #phi(#tau2)","Events");
  // Tau2Eta=HConfig.GetTH1D(Name+"_Tau2Eta","Pseudorapidity Tau2",15,-2.7,2.7," #eta(#tau2)","Events");
  // Tau2dz=HConfig.GetTH1D(Name+"_Tau2dz","Tau2dz",20,-0.02,0.02,"Tau2 dz","Events");
  // Tau2HPSDecayMode=HConfig.GetTH1D(Name+"_Tau2HPSDecayMode","Decay mode of the selected #tau candidate",12,-0.5,11.5," HPS Mode ","Events");
  // Tau2MVADecayMode=HConfig.GetTH1D(Name+"_Tau2MVADecayMode","MVA decay mode of the selected #tau2 candidate",13,-1.5,11.5,"MVA DM","Events");
  // Tau2GenMatch=HConfig.GetTH1D(Name+"_Tau2GenMatch","GenMatch of the selected #tau2 candidate",7,0,7,"GenMatch","Events");
  
  //Tau1isolation=HConfig.GetTH1D(Name+"_Tau1isolation","First Tau isolation 0- VVVLoose, 1- VVLoose, 2 VLoose, 3-Loose, 4-Medium, 5-Tight, 6-VTight, 7-VVTight",8,-0.5,7.5,"Discriminator","Events");
  //Tau2isolation=HConfig.GetTH1D(Name+"_Tau2isolation","Second Tau isolation 0- VVVLoose, 1- VVLoose, 2 VLoose, 3-Loose, 4-Medium, 5-Tight, 6-VTight, 7-VVTight",8,-0.5,7.5," Discriminator","Events");
  
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
  // ExtraLeptonVeto=HConfig.GetTH1D(Name+"_ExtraLeptonVeto","ExtraLeptonVeto",2,-0.5,1.5,"ExtraLeptonVeto","Events");
  
  // QCDShape=HConfig.GetTH1D(Name+"_QCDShape","QCDShape",2,-0.5,1.5,"QCD Shape","");

  // TauTauVisMass=HConfig.GetTH1D(Name+"_TauTauVisMass","Visible invariant mass of a tau pair",23,20,250," M(#tau#tau)_{vis}, GeV","Events");
  // TauTauFullMass=HConfig.GetTH1D(Name+"_TauTauFullMass","Full invariant mass of a tau pair",15,50,300," M(#tau#tau)_{full}, GeV","Events");

  // NQCD=HConfig.GetTH1D(Name+"_NQCD","NQCD",4,0.5,4.5,"NQCD in ABCD","Events");
  // // TauTauFullMass_B=HConfig.GetTH1D(Name+"_TauTauFullMass_B","TauTauFullMass_B",40,0,150,"#tau_h#tau_h SVFit Mass in B","Events");
  // // TauTauFullMass_C=HConfig.GetTH1D(Name+"_TauTauFullMass_C","TauTauFullMass_C",40,0,150,"#tau_h#tau_h SVFit Mass in C","Events");
  // // TauTauFullMass_D=HConfig.GetTH1D(Name+"_TauTauFullMass_D","TauTauFullMass_D",40,0,150,"#tau_h#tau_h SVFit Mass in D","Events");

  // NFFData=HConfig.GetTH1D(Name+"_NFFData","NFFData",1,0.5,1.5,"NFFData","Events");
  //NFFLeadMC=HConfig.GetTH1D(Name+"_NFFLeadMC","NFFLeadMC",1,0.5,1.5,"NFFLeadMC","Events");
  WfakesHiggs=HConfig.GetTH2D(Name+"_WfakesHiggs","WfakesHiggs",60,0.,2*TMath::Pi(),7,0.3,1,"Wfakes","Events");
  WfakesJetFakes=HConfig.GetTH2D(Name+"_WfakesJetFakes","WfakesJetFakes",60,0.,2*TMath::Pi(),7,0.3,1,"Wfakes","Events");
  WfakesZTT=HConfig.GetTH2D(Name+"_WfakesZTT","WfakesZTT",60,0.,2*TMath::Pi(),7,0.3,1,"Wfakes","Events");

  WfakesHiggs_DP=HConfig.GetTH2D(Name+"_WfakesHiggs_DP","WfakesHiggs_DP",60,0.,2*TMath::Pi(),7,0.3,1,"Wfakes","Events");
  WfakesJetFakes_DP=HConfig.GetTH2D(Name+"_WfakesJetFakes_DP","WfakesJetFakes_DP",60,0.,2*TMath::Pi(),7,0.3,1,"Wfakes","Events");
  WfakesZTT_DP=HConfig.GetTH2D(Name+"_WfakesZTT_DP","WfakesZTT_DP",60,0.,2*TMath::Pi(),7,0.3,1,"Wfakes","Events");
  
  // dRTauTau=HConfig.GetTH1D(Name+"_dRTauTau","#Delta R",20,0.,3.5," #Delta R","Events");

  // MET=HConfig.GetTH1D(Name+"_MET","MET",20,0,80,"MET, GeV","Events");
  // METphi=HConfig.GetTH1D(Name+"_METphi","METphi",10,-3.5,3.5,"METphi","Events");
  // PUPPImet=HConfig.GetTH1D(Name+"_PUPPImet","PUPPImet",10,0,75,"PUPPImet, GeV","Events");
  // PUPPImetphi=HConfig.GetTH1D(Name+"_PUPPImetphi","PUPPImetphi",10,-3.5,3.5,"PUPPImetphi","Events");
  //  PUPPImetcorr=HConfig.GetTH1D(Name+"_PUPPImetcorr","PUPPImetcorr",10,0,100,"PUPPImetcorr, GeV","Events");
  // PUPPImetcorrphi=HConfig.GetTH1D(Name+"_PUPPImetcorrphi","PUPPImetcorrphi",10,-3.5,3.5,"PUPPImetcorrphi","Events");
  // TransverseMass=HConfig.GetTH1D(Name+"_TransverseMass","TransverseMass, GeV",40,0,110,"TransverseMass","Events");
  
  // NPrimeVtx=HConfig.GetTH1D(Name+"_NPrimeVtx","NPrimeVtx",10,0,50,"N vtx","Events");
  // NPU=HConfig.GetTH1D(Name+"_npu","npu",10,0,50,"N pu","Events");
  // RHO=HConfig.GetTH1D(Name+"_rho","rho",10,0,30,"rho","Events");
  
  //  NbJets=HConfig.GetTH1D(Name+"_NbJets","NbJets",4,-0.5,3.5,"N_{jets}","Events");
																			    
  // HPtVis=HConfig.GetTH1D(Name+"_HPtVis","Visible Pt_{H}",40,0,200,"","Events");

  //TauTauVisPT=HConfig.GetTH1D(Name+"_TauTauVisPT","Visible Pt_{#tau#tau}",20,0,100,"Visible Pt_{#tau#tau}","Events");
  //TauTauFullPT=HConfig.GetTH1D(Name+"_TauTauFullPT","Full Pt_{#tau#tau}",20,0,150,"Full Pt_{#tau#tau}","Events");
  // Mjj=HConfig.GetTH1D(Name+"_Mjj","m_{jj}",15,0,700,"m_{jj} GeV","Events");
  // dijetpt=HConfig.GetTH1D(Name+"_dijetpt","dijetpt",15,0,300,"dijetpt GeV","Events");
  // dijetphi=HConfig.GetTH1D(Name+"_dijetphi","dijetphi",10,-TMath::Pi(),TMath::Pi(),"dijetphi","Events");
  // jdeta=HConfig.GetTH1D(Name+"_jdeta","jdeta",10,0,9,"jdeta","Events");
  // jdphi=HConfig.GetTH1D(Name+"_jdphi","jdphi",10,-TMath::Pi(),TMath::Pi(),"jdphi","Events");
  // jpt_2=HConfig.GetTH1D(Name+"_jpt_2","jpt_2",20,0,200,"jpt_2 GeV","Events");
  // jeta_2=HConfig.GetTH1D(Name+"_jeta_2","jeta_2",20,-4.7,4.7,"jeta_2","Events");
  // jphi_2=HConfig.GetTH1D(Name+"_jphi_2","jphi_2",10,-TMath::Pi(),TMath::Pi(),"jphi_2","Events");
  // jpt_1=HConfig.GetTH1D(Name+"_jpt_1","jpt_1",10,0,150,"jpt_1 GeV","Events");
  // jeta_1=HConfig.GetTH1D(Name+"_jeta_1","jeta_1",20,-4.7,4.7,"jeta_1","Events");
  // jphi_1=HConfig.GetTH1D(Name+"_jphi_1","jphi_1",10,-TMath::Pi(),TMath::Pi(),"jphi_1","Events");

  // ResolPullTauTauFroma1a1MZMomentum=HConfig.GetTH1D(Name+"_ResolPullTauTauFroma1a1MZMomentum","ResolPullTauTauFroma1a1MZMomentum",30,-1,1,"","Events");


  // ResolPullTauminusFroma1a1MZMomentum=HConfig.GetTH1D(Name+"_ResolPullTauminusFroma1a1MZMomentum","ResolPullTauminusFroma1a1MZMomentum",30,-1,1,"","Events");


  // ResolPullTauplusFroma1a1MZMomentum=HConfig.GetTH1D(Name+"_ResolPullTauplusFroma1a1MZMomentum","ResolPullTauplusFroma1a1MZMomentum",30,-1,1,"","Events");



  // ResolPullTauFroma1a1MZMomentum=HConfig.GetTH1D(Name+"_ResolPullTauFroma1a1MZMomentum","ResolPullTauFroma1a1MZMomentum",30,-1,1,"","Events");


  // ResolPullXVtxIna1a1=HConfig.GetTH1D(Name+"_ResolPullXVtxIna1a1","ResolPullXVtxIna1a1",30,-1,1,"","Events");
  // ResolPullYVtxIna1a1=HConfig.GetTH1D(Name+"_ResolPullYVtxIna1a1","ResolPullYVtxIna1a1",30,-1,1,"","Events");
  // ResolPullZVtxIna1a1=HConfig.GetTH1D(Name+"_ResolPullZVtxIna1a1","ResolPullZVtxIna1a1",30,-1,1,"","Events");

  // tauminusa1a1MomentumPairConstraint=HConfig.GetTH1D(Name+"_tauminusa1a1MomentumPairConstraint","tauminusa1a1MomentumPairConstraint",15,0,200,"","Events");       
  // tauplusa1a1MomentumPairConstraint=HConfig.GetTH1D(Name+"_tauplusa1a1MomentumPairConstraint","tauplusa1a1MomentumPairConstraint",15,0,200,"","Events");

  // polarimetricAcopAngle=HConfig.GetTH1D(Name+"_polarimetricAcopAngle","GEF",11,0.,2*TMath::Pi(),"GEF","Events");

  // polarimetricAcopAnglePVRefitNoBS=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBS","GEF PVRefit no BS constraint Tracks Removed",11,0.,2*TMath::Pi(),"GEF PVRefit no BS constraint Tracks Removed","Events");
  // polarimetricAcopAnglePVRefitBS=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBS","GEF PVRefit BS constraint Tracks Removed",11,0.,2*TMath::Pi(),"GEF PVRefit BS constraint Tracks Removed","Events");
  // polarimetricAcopAnglePVRefitNoBSZNominal=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBSZNominal","GEF PVRefit no BS constraint nominal Z Tracks Removed",11,0.,2*TMath::Pi(),"GEF PVRefit no BS constraint nominal Z Tracks Removed","Events");
  // polarimetricAcopAnglePVRefitBSZNominal=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSZNominal","GEF PVRefit BS constraint nominal Z Tracks Removed",11,0.,2*TMath::Pi(),"GEF PVRefit BS constraint nominal Z Tracks Removed","Events");


  // polarimetricAcopAngleMVADM=HConfig.GetTH1D(Name+"_polarimetricAcopAngleMVADM","PV nominal",11,0.,2*TMath::Pi(),"PV nominal","Events");
  // polarimetricAcopAnglePVRefitBSMVADM=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSMVADM","PV refitted",11,0.,2*TMath::Pi(),"PV refitted","Events");
  // polarimetricAcopAnglePVRefitBSZNominalMVADM=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSZNominalMVADM","PV refitted Z nominal",11,0.,2*TMath::Pi(),"PV refitted Z nominal","Events");
  //   polarimetricAcopAnglePVRefitWithTracksBSMVADM=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADM","PV refitted with taus tracks",11,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");
  // polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADM=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADM","PV refitted with taus tracks Z nominal",11,0.,2*TMath::Pi(),"PV refitted with taus tracks Z nominal","Events");
  
  // polarimetricAcopAngleMVADMHiggs=HConfig.GetTH1D(Name+"_polarimetricAcopAngleMVADMHiggs","PV nominal",11,0.,2*TMath::Pi(),"PV nominal","Events");
  // polarimetricAcopAnglePVRefitBSMVADMHiggs=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSMVADMHiggs","PV refitted",11,0.,2*TMath::Pi(),"PV refitted","Events");
  // polarimetricAcopAnglePVRefitBSZNominalMVADMHiggs=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSZNominalMVADMHiggs","PV refitted Z nominal",11,0.,2*TMath::Pi(),"PV refitted Z nominal","Events");
  // polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggs=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggs","PV refitted with taus tracks",11,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");
  // polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMHiggs=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMHiggs","PV refitted with taus tracks Z nominal",11,0.,2*TMath::Pi(),"PV refitted with taus tracks Z nominal","Events");

  // polarimetricAcopAngleMVADMJetFakes=HConfig.GetTH1D(Name+"_polarimetricAcopAngleMVADMJetFakes","PV nominal",11,0.,2*TMath::Pi(),"PV nominal","Events");
  // polarimetricAcopAnglePVRefitBSMVADMJetFakes=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSMVADMJetFakes","PV refitted",11,0.,2*TMath::Pi(),"PV refitted","Events");
  // polarimetricAcopAnglePVRefitBSZNominalMVADMJetFakes=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSZNominalMVADMJetFakes","PV refitted Z nominal",11,0.,2*TMath::Pi(),"PV refitted Z nominal","Events");
  // polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakes=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakes","PV refitted with taus tracks",11,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");
  // polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMJetFakes=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMJetFakes","PV refitted with taus tracks Z nominal",11,0.,2*TMath::Pi(),"PV refitted with taus tracks Z nominal","Events");

  // polarimetricAcopAngleMVADMZTTEmbed=HConfig.GetTH1D(Name+"_polarimetricAcopAngleMVADMZTTEmbed","PV nominal",11,0.,2*TMath::Pi(),"PV nominal","Events");
  // polarimetricAcopAnglePVRefitBSMVADMZTTEmbed=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSMVADMZTTEmbed","PV refitted",11,0.,2*TMath::Pi(),"PV refitted","Events");
  // polarimetricAcopAnglePVRefitBSZNominalMVADMZTTEmbed=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSZNominalMVADMZTTEmbed","PV refitted Z nominal",11,0.,2*TMath::Pi(),"PV refitted Z nominal","Events");
  // polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTEmbed=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTEmbed","PV refitted with taus tracks",11,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");
  // polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMZTTEmbed=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMZTTEmbed","PV refitted with taus tracks Z nominal",11,0.,2*TMath::Pi(),"PV refitted with taus tracks Z nominal","Events");
  

  // polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled0005=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled0005","PV refitted with taus tracks",5,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");
  // polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled0506=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled0506","PV refitted with taus tracks",5,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");
  // polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled0607=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled0607","PV refitted with taus tracks",5,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");
  // polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled0710=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled0710","PV refitted with taus tracks",5,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");

  polarimetricAcopAnglePVRefitWithTracksBSMVADM=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADM","PV refitted with taus tracks",60,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");

  polarimetricAcopAnglePVRefitWithTracksBSMVADM_DP=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADM_DP","PV refitted with taus tracks",10,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");
  
  polarimetricAcopAnglePVRefitWithTracksBSMVADMQCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMQCDMC","PV refitted with taus tracks",60,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");

  
  polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled","PV refitted with taus tracks",180,0.,3*2*TMath::Pi(),"PV refitted with taus tracks","Events");
  
  polarimetricAcopAnglePVRefitWithTracksBSMVADMQCDMC_DP=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMQCDMC_DP","PV refitted with taus tracks",60,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");
  polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled_DP=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled_DP","PV refitted with taus tracks",180,0.,3*2*TMath::Pi(),"PV refitted with taus tracks","Events");
  
  // polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled0005=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled0005","PV refitted with taus tracks",5,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");
  // polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled0506=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled0506","PV refitted with taus tracks",5,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");
  // polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled0607=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled0607","PV refitted with taus tracks",5,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");
  // polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled0710=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled0710","PV refitted with taus tracks",5,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");
  
  polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled","PV refitted with taus tracks",180,0.,3*2*TMath::Pi(),"PV refitted with taus tracks","Events");

  polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled_DP=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled_DP","PV refitted with taus tracks",180,0.,3*2*TMath::Pi(),"PV refitted with taus tracks","Events");
  
  // polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled0005=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled0005","PV refitted with taus tracks",5,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");
  // polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled0506=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled0506","PV refitted with taus tracks",5,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");
  // polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled0607=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled0607","PV refitted with taus tracks",5,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");
  // polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled0710=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled0710","PV refitted with taus tracks",5,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");

  polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled","PV refitted with taus tracks",180,0.,3*2*TMath::Pi(),"PV refitted with taus tracks","Events");
  
  polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled_DP=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled_DP","PV refitted with taus tracks",180,0.,3*2*TMath::Pi(),"PV refitted with taus tracks","Events");
  
  // polarimetricAcopAngleMVADMHiggs=HConfig.GetTH2D(Name+"_polarimetricAcopAngleMVADMHiggs","Higgs MVA score VS acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA score");
  // polarimetricAcopAnglePVRefitBSMVADMHiggs=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitBSMVADMHiggs","Higgs MVA Score VS acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA score");
  // polarimetricAcopAnglePVRefitBSZNominalMVADMHiggs=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitBSZNominalMVADMHiggs","Higgs MVA Score VS acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA score");
  polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggs=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggs","Higgs MVA Score VS acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA core");

  polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggs_DP=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggs_DP","Higgs MVA Score VS acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA core");

  //  polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMHiggs=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMHiggs","Higgs MVA Score VS acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA score");

  // polarimetricAcopAngleMVADMJetFakes=HConfig.GetTH2D(Name+"_polarimetricAcopAngleMVADMJetFakes","JetFakes MVA score VS acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA score");
  // polarimetricAcopAnglePVRefitBSMVADMJetFakes=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitBSMVADMJetFakes","JetFakes MVA Score VS acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA score");
  // polarimetricAcopAnglePVRefitBSZNominalMVADMJetFakes=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitBSZNominalMVADMJetFakes","JetFakes MVA Score VS acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA score");
  polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakes=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakes","JetFakes MVA Score VS acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA core");

  polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakes_DP=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakes_DP","JetFakes MVA Score VS acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA core");
  
  //  polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMJetFakes=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMJetFakes","JetFakes MVA Score VS acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA score");

  // polarimetricAcopAngleMVADMZTT=HConfig.GetTH2D(Name+"_polarimetricAcopAngleMVADMZTT","ZTT MVA score VS acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA score");
  // polarimetricAcopAnglePVRefitBSMVADMZTT=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitBSMVADMZTT","ZTT MVA Score VS acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA score");
  // polarimetricAcopAnglePVRefitBSZNominalMVADMZTT=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitBSZNominalMVADMZTT","ZTT MVA Score VS acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA score");
  polarimetricAcopAnglePVRefitWithTracksBSMVADMZTT=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMZTT","ZTT MVA Score VS acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA core");
  
  polarimetricAcopAnglePVRefitWithTracksBSMVADMZTT_DP=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMZTT_DP","ZTT MVA Score VS acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA core");
  
  //  polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMZTT=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMZTT","ZTT MVA Score VS acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA score");


  // polarimetricAcopAnglePVRefitNoBSMVADM=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBSMVADM","GEF PVRefit no BS constraint MVADM",11,0.,2*TMath::Pi(),"GEF PVRefit no BS constraint MVADM","Events");
  // polarimetricAcopAnglePVRefitNoBSZNominalMVADM=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBSZNominalMVADM","GEF PVRefit no BS constraint nominal Z MVADM",11,0.,2*TMath::Pi(),"GEF PVRefit no BS constraint nominal Z MVADM","Events");

  // PVXResol=HConfig.GetTH1D(Name+"_PVXResol","PV_{X} pull",50,-0.1,0.1,"PV_{X} pull","Events");
 
  // PVXNoBSResol=HConfig.GetTH1D(Name+"_PVXNoBSResol","No BS PV_{X}Refit pull",50,-0.1,0.1,"No BS PV_{X}Refit pull","Events");
  
  // PVXBSResol=HConfig.GetTH1D(Name+"_PVXBSResol","BS PV_{X}Refit pull",50,-0.1,0.1,"BS PV_{X}Refit pull","Events");
  
  // PVYResol=HConfig.GetTH1D(Name+"_PVYResol","PV_{Y} pull",50,-0.1,0.1,"PV_{Y} pull","Events");
  
  // PVYNoBSResol=HConfig.GetTH1D(Name+"_PVYNoBSResol","No BS PV_{Y}Refit pull",50,-0.1,0.1,"No BS PV_{Y}Refit pull","Events");
  
  // PVYBSResol=HConfig.GetTH1D(Name+"_PVYBSResol","BS PV_{Y}Refit pull",50,-0.1,0.1,"BS PV_{Y}Refit pull","Events");

  // PVZResol=HConfig.GetTH1D(Name+"_PVZResol","PV_{Z} pull",50,-0.01,0.01,"PV_{Z} pull","Events");

  // PVZNoBSResol=HConfig.GetTH1D(Name+"_PVZNoBSResol","No BS PV_{Z}Refit pull",50,-0.01,0.01,"No BS PV_{Z}Refit pull","Events");

  // PVZBSResol=HConfig.GetTH1D(Name+"_PVZBSResol","BS PV_{Z}Refit pull",50,-0.01,0.01,"BS PV_{Z}Refit pull","Events");
  

  // PVXNoBSOnlyResol=HConfig.GetTH1D(Name+"_PVXNoBSOnlyResol","No BS PV_{X}Refit Only pull",50,-0.1,0.1,"No BS PV_{X}Refit Only pull","Events");
 
  // PVXBSOnlyResol=HConfig.GetTH1D(Name+"_PVXBSOnlyResol","BS PV_{X}Refit Only pull",50,-0.1,0.1,"BS PV_{X}Refit Only pull","Events");

  // PVYNoBSOnlyResol=HConfig.GetTH1D(Name+"_PVYNoBSOnlyResol","No BS PV_{Y}Refit Only pull",50,-0.1,0.1,"No BS PV_{Y}Refit Only pull","Events");
  // PVYBSOnlyResol=HConfig.GetTH1D(Name+"_PVYBSOnlyResol","BS PV_{Y}Refit Only pull",50,-0.1,0.1,"BS PV_{Y}Refit Only pull","Events");

  // PVZNoBSOnlyResol=HConfig.GetTH1D(Name+"_PVZNoBSOnlyResol","No BS PV_{Z}Refit Only pull",50,-0.01,0.01,"No BS PV_{Z}Refit Only pull","Events");

  // PVZBSOnlyResol=HConfig.GetTH1D(Name+"_PVZBSOnlyResol","BS PV_{Z}Refit Only pull",50,-0.01,0.01,"BS PV_{Z}Refit Only pull","Events");
  
  // polarimetricAcopAngleTruthA1=HConfig.GetTH1D(Name+"_polarimetricAcopAngleTruthA1","polarimetricAcopAngleTruthA1",11,0,2*TMath::Pi(),"","Events");

  // polarimetricAcopAngleDecayPlane=HConfig.GetTH1D(Name+"_polarimetricAcopAngleDecayPlane","polarimetricAcopAngleDecayPlane",11,0.,2*TMath::Pi(),"","Events");
  
  // test=HConfig.GetTH1D(Name+"_test","test",11,0.,2*TMath::Pi(),"","Events"); 

  // PurityDM=HConfig.GetTH1D(Name+"_PurityDM","PurityDM",29,0.,29,"","Events");
  // PurityNewMVA=HConfig.GetTH1D(Name+"_PurityNewMVA","PurityNewMVA",29,0.,29,"","Events");

  // TauPxResPull=HConfig.GetTH1D(Name+"_TauPxResPull","TauPxResPull",50,-2,2,"P_{X} pull of #tau","Events");
  // TauPyResPull=HConfig.GetTH1D(Name+"_TauPyResPull","TauPyResPull",50,-2,2,"P_{Y} pull of #tau","Events");
  // TauPzResPull=HConfig.GetTH1D(Name+"_TauPzResPull","TauPzResPull",50,-2,2,"P_{Z} pull of #tau","Events");

  // TauPxResPullMVA=HConfig.GetTH1D(Name+"_TauPxResPullMVA","TauPxResPullMVA",50,-2,2,"P_{X} pull of #tau","Events");
  // TauPyResPullMVA=HConfig.GetTH1D(Name+"_TauPyResPullMVA","TauPyResPullMVA",50,-2,2,"P_{Y} pull of #tau","Events");
  // TauPzResPullMVA=HConfig.GetTH1D(Name+"_TauPzResPullMVA","TauPzResPullMVA",50,-2,2,"P_{Z} pull of #tau","Events");

  //HiggsBDTScore=HConfig.GetTH1D(Name+"_HiggsBDTScore","HiggsBDTScore",7,0.3,1,"MVA Score","Events");
  //JetFakesBDTScore=HConfig.GetTH1D(Name+"_JetFakesBDTScore","JetFakesBDTScore",7,0.3,1,"MVA Score","Events");
  //ZTTBDTScore=HConfig.GetTH1D(Name+"_ZTTBDTScore","ZTTBDTScore",7,0.3,1,"MVA Score","Events");


  //PUPPImetcorrQCDMC=HConfig.GetTH1D(Name+"_PUPPImetcorrQCDMC","PUPPImetcorrQCDMC",20,0,200,"PUPPImetcorr, GeV","Events");
  // PUPPImetcorrphiQCDMC=HConfig.GetTH1D(Name+"_PUPPImetcorrphiQCDMC","PUPPImetcorrphi",10,-3.5,3.5,"PUPPImetcorrphi","Events");
  
  //Tau1PTQCDMC=HConfig.GetTH1D(Name+"_Tau1PTQCDMC","Transverse momentum of selected #tau1 candidate",20,40,140," P_{T}(#tau1), GeV","Events");
  // Tau1EQCDMC=HConfig.GetTH1D(Name+"_Tau1EQCDMC","Energy of selected #tau1 candidate",20,35.,140," E(#tau1), GeV","Events");
  // Tau1MassQCDMC=HConfig.GetTH1D(Name+"_Tau1MassQCDMC","Mass of selected #tau1 candidate",20,0,2.," M(#tau1), GeV","Events");
  // Tau1PhiQCDMC=HConfig.GetTH1D(Name+"_Tau1PhiQCDMC","Phi of selected #tau1 candidate",10,-3.5,3.5," #phi(#tau1)","Events");
  // Tau1EtaQCDMC=HConfig.GetTH1D(Name+"_Tau1EtaQCDMC","Pseudorapidity tau1",15,-2.7,2.7," #eta(#tau1)","Events");
  // Tau1dzQCDMC=HConfig.GetTH1D(Name+"_Tau1dzQCDMC","Tau1 dz",20,-0.02,0.02,"Tau1 dz","Events");
  // Tau1HPSDecayModeQCDMC=HConfig.GetTH1D(Name+"_Tau1HPSDecayModeQCDMC","Decay mode of the selected #tau candidate",12,-0.5,11.5,"HPS Mode","Events");
  // Tau1MVADecayModeQCDMC=HConfig.GetTH1D(Name+"_Tau1MVADecayModeQCDMC","MVA decay mode of the selected #tau1 candidate",13,-1.5,11.5,"MVA DM","Events");
  // Tau1GenMatchQCDMC=HConfig.GetTH1D(Name+"_Tau1GenMatchQCDMC","GenMatch of the selected #tau1 candidate",7,0,7,"GenMatch","Events");
  
  // Tau2PTQCDMC=HConfig.GetTH1D(Name+"_Tau2PTQCDMC","Transverse momentum of selected #tau2 candidate",12,40,100," P_{T}(#tau2), GeV","Events");
  // Tau2EQCDMC=HConfig.GetTH1D(Name+"_Tau2EQCDMC","Energy of selected #tau2 candidate",20,30.,140," E(#tau2), GeV","Events");
  // Tau2MassQCDMC=HConfig.GetTH1D(Name+"_Tau2MassQCDMC","Mass of selected #tau2 candidate",20,0,2.," M(#tau2), GeV","Events");
  // Tau2PhiQCDMC=HConfig.GetTH1D(Name+"_Tau2PhiQCDMC","Phi of selected #tau2 candidate",10,-3.5,3.5," #phi(#tau2)","Events");
  // Tau2EtaQCDMC=HConfig.GetTH1D(Name+"_Tau2EtaQCDMC","Pseudorapidity Tau2",15,-2.7,2.7," #eta(#tau2)","Events");
  // Tau2dzQCDMC=HConfig.GetTH1D(Name+"_Tau2dzQCDMC","Tau2dz",20,-0.02,0.02,"Tau2 dz","Events");
  // Tau2HPSDecayModeQCDMC=HConfig.GetTH1D(Name+"_Tau2HPSDecayModeQCDMC","Decay mode of the selected #tau candidate",12,-0.5,11.5," HPS Mode ","Events");
  // Tau2MVADecayModeQCDMC=HConfig.GetTH1D(Name+"_Tau2MVADecayModeQCDMC","MVA decay mode of the selected #tau2 candidate",13,-1.5,11.5,"MVA DM","Events");
  // Tau2GenMatchQCDMC=HConfig.GetTH1D(Name+"_Tau2GenMatchQCDMC","GenMatch of the selected #tau2 candidate",7,0,7,"GenMatch","Events");
  
  // NbJetsQCDMC=HConfig.GetTH1D(Name+"_NbJetsQCDMC","NbJetsQCDMC",5,-0.5,4.5,"N_{jets}","Events");
  // TauTauVisMassQCDMC=HConfig.GetTH1D(Name+"_TauTauVisMassQCDMC","Visible invariant mass of a tau pair",23,20,250," M(#tau#tau)_{vis}, GeV","Events");
  // TauTauFullMassQCDMC=HConfig.GetTH1D(Name+"_TauTauFullMassQCDMC","Full invariant mass of a tau pair",15,50,300," M(#tau#tau)_{full}, GeV","Events");
  // TauTauVisPTQCDMC=HConfig.GetTH1D(Name+"_TauTauVisPTQCDMC","Visible Pt_{#tau#tau}",20,0,100,"Visible Pt_{#tau#tau}","Events");
  // TauTauFullPTQCDMC=HConfig.GetTH1D(Name+"_TauTauFullPTQCDMC","Full Pt_{#tau#tau}",20,0,150,"Full Pt_{#tau#tau}","Events");
  // MjjQCDMC=HConfig.GetTH1D(Name+"_MjjQCDMC","m_{jj}",15,0,1500,"m_{jj} GeV","Events");
  // dijetptQCDMC=HConfig.GetTH1D(Name+"_dijetptQCDMC","dijetptQCDMC",15,0,300,"dijetpt GeV","Events");
  // dijetphiQCDMC=HConfig.GetTH1D(Name+"_dijetphiQCDMC","dijetphiQCDMC",10,-TMath::Pi(),TMath::Pi(),"dijetphi","Events");
  //jdetaQCDMC=HConfig.GetTH1D(Name+"_jdetaQCDMC","jdetaQCDMC",10,0,9,"jdeta","Events");
  // jdphiQCDMC=HConfig.GetTH1D(Name+"_jdphiQCDMC","jdphiQCDMC",10,-TMath::Pi(),TMath::Pi(),"jdphi","Events");
  // jpt_2QCDMC=HConfig.GetTH1D(Name+"_jpt_2QCDMC","jpt_2QCDMC",20,0,200,"jpt_2 GeV","Events");
  // jeta_2QCDMC=HConfig.GetTH1D(Name+"_jeta_2QCDMC","jeta_2QCDMC",20,-4.7,4.7,"jeta_2","Events");
  // jphi_2QCDMC=HConfig.GetTH1D(Name+"_jphi_2QCDMC","jphi_2QCDMC",10,-TMath::Pi(),TMath::Pi(),"jphi_2","Events");
  //jpt_1QCDMC=HConfig.GetTH1D(Name+"_jpt_1QCDMC","jpt_1QCDMC",20,0,500,"jpt_1 GeV","Events");
  // jeta_1QCDMC=HConfig.GetTH1D(Name+"_jeta_1QCDMC","jeta_1QCDMC",20,-4.7,4.7,"jeta_1","Events");
  // jphi_1QCDMC=HConfig.GetTH1D(Name+"_jphi_1QCDMC","jphi_1QCDMC",10,-TMath::Pi(),TMath::Pi(),"jphi_1","Events");

  //HiggsBDTScoreQCDMC=HConfig.GetTH1D(Name+"_HiggsBDTScoreQCDMC","HiggsBDTScoreQCDMC",7,0.3,1,"MVA Score","Events");
  //JetFakesBDTScoreQCDMC=HConfig.GetTH1D(Name+"_JetFakesBDTScoreQCDMC","JetFakesBDTScoreQCDMC",7,0.3,1,"MVA Score","Events");
  //ZTTBDTScoreQCDMC=HConfig.GetTH1D(Name+"_ZTTBDTScoreQCDMC","ZTTBDTScoreQCDMC",7,0.3,1,"MVA Score","Events");

  Tau1PTa1a1=HConfig.GetTH1D(Name+"_Tau1PTa1a1","Transverse momentum of selected #tau1 candidate",10,40,140," P_{T}(#tau1), GeV","Events");
  // Tau1Ea1a1=HConfig.GetTH1D(Name+"_Tau1Ea1a1","Energy of selected #tau1 candidate",10,35.,140," E(#tau1), GeV","Events");
  // Tau1Massa1a1=HConfig.GetTH1D(Name+"_Tau1Massa1a1","Mass of selected #tau1 candidate",20,0,2.," M(#tau1), GeV","Events");
  // Tau1Phia1a1=HConfig.GetTH1D(Name+"_Tau1Phia1a1","Phi of selected #tau1 candidate",10,-3.5,3.5," #phi(#tau1)","Events");
  // Tau1Etaa1a1=HConfig.GetTH1D(Name+"_Tau1Etaa1a1","Pseudorapidity tau1",15,-2.7,2.7," #eta(#tau1)","Events");
  // Tau1dza1a1=HConfig.GetTH1D(Name+"_Tau1dza1a1","Tau1 dz",20,-0.02,0.02,"Tau1 dz","Events");
  
  // Tau2PTa1a1=HConfig.GetTH1D(Name+"_Tau2PTa1a1","Transverse momentum of selected #tau2 candidate",12,40,100," P_{T}(#tau2), GeV","Events");
  // Tau2Ea1a1=HConfig.GetTH1D(Name+"_Tau2Ea1a1","Energy of selected #tau2 candidate",10,30.,140," E(#tau2), GeV","Events");
  // Tau2Massa1a1=HConfig.GetTH1D(Name+"_Tau2Massa1a1","Mass of selected #tau2 candidate",20,0,2.," M(#tau2), GeV","Events");
  // Tau2Phia1a1=HConfig.GetTH1D(Name+"_Tau2Phia1a1","Phi of selected #tau2 candidate",10,-3.5,3.5," #phi(#tau2)","Events");
  // Tau2Etaa1a1=HConfig.GetTH1D(Name+"_Tau2Etaa1a1","Pseudorapidity Tau2",15,-2.7,2.7," #eta(#tau2)","Events");
  // Tau2dza1a1=HConfig.GetTH1D(Name+"_Tau2dza1a1","Tau2dz",20,-0.02,0.02,"Tau2 dz","Events");

  TauTauVisMassa1a1=HConfig.GetTH1D(Name+"_TauTauVisMassa1a1","Visible invariant mass of a tau pair",15,20,250," M(#tau#tau)_{vis}, GeV","Events");
  TauTauFullMassa1a1=HConfig.GetTH1D(Name+"_TauTauFullMassa1a1","Full invariant mass of a tau pair",15,50,300," M(#tau#tau)_{full}, GeV","Events");
  PUPPImetcorra1a1=HConfig.GetTH1D(Name+"_PUPPImetcorra1a1","PUPPImetcorr",7,0.3,100,"PUPPImetcorr, GeV","Events");
  // PUPPImetcorrphia1a1=HConfig.GetTH1D(Name+"_PUPPImetcorrphia1a1","PUPPImetcorrphi",10,-3.5,3.5,"PUPPImetcorrphi","Events");
  NbJetsa1a1=HConfig.GetTH1D(Name+"_NbJetsa1a1","NbJets",5,-0.5,4.5,"N_{jets}","Events");
  TauTauVisPTa1a1=HConfig.GetTH1D(Name+"_TauTauVisPTa1a1","Visible Pt_{#tau#tau}",20,0,100,"Visible Pt_{#tau#tau}","Events");
  TauTauFullPTa1a1=HConfig.GetTH1D(Name+"_TauTauFullPTa1a1","Full Pt_{#tau#tau}",20,0,150,"Full Pt_{#tau#tau}","Events");
  Mjja1a1=HConfig.GetTH1D(Name+"_Mjja1a1","m_{jj}",10,0,800,"m_{jj} GeV","Events");
  // dijetpta1a1=HConfig.GetTH1D(Name+"_dijetpta1a1","dijetpt",10,0,200,"dijetpt GeV","Events");
  // dijetphia1a1=HConfig.GetTH1D(Name+"_dijetphia1a1","dijetphi",10,-TMath::Pi(),TMath::Pi(),"dijetphi","Events");
  jdetaa1a1=HConfig.GetTH1D(Name+"_jdetaa1a1","jdeta",7,0,7,"jdeta","Events");
  // jdphia1a1=HConfig.GetTH1D(Name+"_jdphia1a1","jdphi",10,-TMath::Pi(),TMath::Pi(),"jdphi","Events");
  // jpt_2a1a1=HConfig.GetTH1D(Name+"_jpt_2a1a1","jpt_2",7,0.3,100,"jpt_2 GeV","Events");
  // jeta_2a1a1=HConfig.GetTH1D(Name+"_jeta_2a1a1","jeta_2",10,-4.7,4.7,"jeta_2","Events");
  // jphi_2a1a1=HConfig.GetTH1D(Name+"_jphi_2a1a1","jphi_2",10,-TMath::Pi(),TMath::Pi(),"jphi_2","Events");
  jpt_1a1a1=HConfig.GetTH1D(Name+"_jpt_1a1a1","jpt_1",10,0,250,"jpt_1 GeV","Events");
  // jeta_1a1a1=HConfig.GetTH1D(Name+"_jeta_1a1a1","jeta_1",10,-4.7,4.7,"jeta_1","Events");
  // jphi_1a1a1=HConfig.GetTH1D(Name+"_jphi_1a1a1","jphi_1",10,-TMath::Pi(),TMath::Pi(),"jphi_1","Events");
  
  HiggsBDTScorea1a1=HConfig.GetTH1D(Name+"_HiggsBDTScorea1a1","HiggsBDTScore",7,0.3,1,"MVA Score","Events");
  JetFakesBDTScorea1a1=HConfig.GetTH1D(Name+"_JetFakesBDTScorea1a1","JetFakesBDTScore",7,0.3,1,"MVA Score","Events");
  ZTTBDTScorea1a1=HConfig.GetTH1D(Name+"_ZTTBDTScorea1a1","ZTTBDTScore",7,0.3,1,"MVA Score","Events");

  HiggsBDTScorea1a1_DP=HConfig.GetTH1D(Name+"_HiggsBDTScorea1a1_DP","HiggsBDTScore",7,0.3,1,"MVA Score","Events");
  JetFakesBDTScorea1a1_DP=HConfig.GetTH1D(Name+"_JetFakesBDTScorea1a1_DP","JetFakesBDTScore",7,0.3,1,"MVA Score","Events");
  ZTTBDTScorea1a1_DP=HConfig.GetTH1D(Name+"_ZTTBDTScorea1a1_DP","ZTTBDTScore",7,0.3,1,"MVA Score","Events");

  Tau1PTa1a1QCDMC=HConfig.GetTH1D(Name+"_Tau1PTa1a1QCDMC","Transverse momentum of selected #tau1 candidate",10,40,140," P_{T}(#tau1), GeV","Events");
  // Tau1Ea1a1QCDMC=HConfig.GetTH1D(Name+"_Tau1Ea1a1QCDMC","Energy of selected #tau1 candidate",10,35.,140," E(#tau1), GeV","Events");
  // Tau1Massa1a1QCDMC=HConfig.GetTH1D(Name+"_Tau1Massa1a1QCDMC","Mass of selected #tau1 candidate",20,0,2.," M(#tau1), GeV","Events");
  // Tau1Phia1a1QCDMC=HConfig.GetTH1D(Name+"_Tau1Phia1a1QCDMC","Phi of selected #tau1 candidate",10,-3.5,3.5," #phi(#tau1)","Events");
  // Tau1Etaa1a1QCDMC=HConfig.GetTH1D(Name+"_Tau1Etaa1a1QCDMC","Pseudorapidity tau1",15,-2.7,2.7," #eta(#tau1)","Events");
  // Tau1dza1a1QCDMC=HConfig.GetTH1D(Name+"_Tau1dza1a1QCDMC","Tau1 dz",20,-0.02,0.02,"Tau1 dz","Events");
  
  // Tau2PTa1a1QCDMC=HConfig.GetTH1D(Name+"_Tau2PTa1a1QCDMC","Transverse momentum of selected #tau2 candidate",12,40,100," P_{T}(#tau2), GeV","Events");
  // Tau2Ea1a1QCDMC=HConfig.GetTH1D(Name+"_Tau2Ea1a1QCDMC","Energy of selected #tau2 candidate",10,30.,140," E(#tau2), GeV","Events");
  // Tau2Massa1a1QCDMC=HConfig.GetTH1D(Name+"_Tau2Massa1a1QCDMC","Mass of selected #tau2 candidate",20,0,2.," M(#tau2), GeV","Events");
  // Tau2Phia1a1QCDMC=HConfig.GetTH1D(Name+"_Tau2Phia1a1QCDMC","Phi of selected #tau2 candidate",10,-3.5,3.5," #phi(#tau2)","Events");
  // Tau2Etaa1a1QCDMC=HConfig.GetTH1D(Name+"_Tau2Etaa1a1QCDMC","Pseudorapidity Tau2",15,-2.7,2.7," #eta(#tau2)","Events");
  // Tau2dza1a1QCDMC=HConfig.GetTH1D(Name+"_Tau2dza1a1QCDMC","Tau2dz",20,-0.02,0.02,"Tau2 dz","Events");

  TauTauVisMassa1a1QCDMC=HConfig.GetTH1D(Name+"_TauTauVisMassa1a1QCDMC","Visible invariant mass of a tau pair",15,20,250," M(#tau#tau)_{vis}, GeV","Events");
  TauTauFullMassa1a1QCDMC=HConfig.GetTH1D(Name+"_TauTauFullMassa1a1QCDMC","Full invariant mass of a tau pair",15,50,300," M(#tau#tau)_{full}, GeV","Events");
  PUPPImetcorra1a1QCDMC=HConfig.GetTH1D(Name+"_PUPPImetcorra1a1QCDMC","PUPPImetcorr",7,0.3,100,"PUPPImetcorr, GeV","Events");
  // PUPPImetcorrphia1a1QCDMC=HConfig.GetTH1D(Name+"_PUPPImetcorrphia1a1QCDMC","PUPPImetcorrphi",10,-3.5,3.5,"PUPPImetcorrphi","Events");
  NbJetsa1a1QCDMC=HConfig.GetTH1D(Name+"_NbJetsa1a1QCDMC","NbJets",5,-0.5,4.5,"N_{jets}","Events");
  TauTauVisPTa1a1QCDMC=HConfig.GetTH1D(Name+"_TauTauVisPTa1a1QCDMC","Visible Pt_{#tau#tau}",20,0,100,"Visible Pt_{#tau#tau}","Events");
  TauTauFullPTa1a1QCDMC=HConfig.GetTH1D(Name+"_TauTauFullPTa1a1QCDMC","Full Pt_{#tau#tau}",20,0,150,"Full Pt_{#tau#tau}","Events");
  Mjja1a1QCDMC=HConfig.GetTH1D(Name+"_Mjja1a1QCDMC","m_{jj}",10,0,800,"m_{jj} GeV","Events");
  // dijetpta1a1QCDMC=HConfig.GetTH1D(Name+"_dijetpta1a1QCDMC","dijetpt",10,0,200,"dijetpt GeV","Events");
  // dijetphia1a1QCDMC=HConfig.GetTH1D(Name+"_dijetphia1a1QCDMC","dijetphi",10,-TMath::Pi(),TMath::Pi(),"dijetphi","Events");
  jdetaa1a1QCDMC=HConfig.GetTH1D(Name+"_jdetaa1a1QCDMC","jdeta",7,0,7,"jdeta","Events");
  // jdphia1a1QCDMC=HConfig.GetTH1D(Name+"_jdphia1a1QCDMC","jdphi",10,-TMath::Pi(),TMath::Pi(),"jdphi","Events");
  // jpt_2a1a1QCDMC=HConfig.GetTH1D(Name+"_jpt_2a1a1QCDMC","jpt_2",7,0.3,100,"jpt_2 GeV","Events");
  // jeta_2a1a1QCDMC=HConfig.GetTH1D(Name+"_jeta_2a1a1QCDMC","jeta_2",10,-4.7,4.7,"jeta_2","Events");
  // jphi_2a1a1QCDMC=HConfig.GetTH1D(Name+"_jphi_2a1a1QCDMC","jphi_2",10,-TMath::Pi(),TMath::Pi(),"jphi_2","Events");
  jpt_1a1a1QCDMC=HConfig.GetTH1D(Name+"_jpt_1a1a1QCDMC","jpt_1",10,0,250,"jpt_1 GeV","Events");
  // jeta_1a1a1QCDMC=HConfig.GetTH1D(Name+"_jeta_1a1a1QCDMC","jeta_1",10,-4.7,4.7,"jeta_1","Events");
  // jphi_1a1a1QCDMC=HConfig.GetTH1D(Name+"_jphi_1a1a1QCDMC","jphi_1",10,-TMath::Pi(),TMath::Pi(),"jphi_1","Events");
  
  HiggsBDTScorea1a1QCDMC=HConfig.GetTH1D(Name+"_HiggsBDTScorea1a1QCDMC","HiggsBDTScore",7,0.3,1,"MVA Score","Events");
  JetFakesBDTScorea1a1QCDMC=HConfig.GetTH1D(Name+"_JetFakesBDTScorea1a1QCDMC","JetFakesBDTScore",7,0.3,1,"MVA Score","Events");
  ZTTBDTScorea1a1QCDMC=HConfig.GetTH1D(Name+"_ZTTBDTScorea1a1QCDMC","ZTTBDTScore",7,0.3,1,"MVA Score","Events");
  
  HiggsBDTScorea1a1QCDMC_DP=HConfig.GetTH1D(Name+"_HiggsBDTScorea1a1QCDMC_DP","HiggsBDTScore",7,0.3,1,"MVA Score","Events");
  JetFakesBDTScorea1a1QCDMC_DP=HConfig.GetTH1D(Name+"_JetFakesBDTScorea1a1QCDMC_DP","JetFakesBDTScore",7,0.3,1,"MVA Score","Events");
  ZTTBDTScorea1a1QCDMC_DP=HConfig.GetTH1D(Name+"_ZTTBDTScorea1a1QCDMC_DP","ZTTBDTScore",7,0.3,1,"MVA Score","Events");

  // polarimetricAcopAngleMVADMQCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAngleMVADMQCDMC","PV nominal",11,0.,2*TMath::Pi(),"PV nominal","Events");
  // polarimetricAcopAnglePVRefitBSMVADMQCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSMVADMQCDMC","PV refitted",11,0.,2*TMath::Pi(),"PV refitted","Events");
  // polarimetricAcopAnglePVRefitBSZNominalMVADMQCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSZNominalMVADMQCDMC","PV refitted Z nominal",11,0.,2*TMath::Pi(),"PV refitted Z nominal","Events");
  //polarimetricAcopAnglePVRefitWithTracksBSMVADMQCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMQCDMC","PV refitted with taus tracks",11,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");
  //polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMQCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMQCDMC","PV refitted with taus tracks Z nominal",11,0.,2*TMath::Pi(),"PV refitted with taus tracks Z nominal","Events");


  // polarimetricAcopAngleMVADMHiggsQCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAngleMVADMHiggsQCDMC","PV nominal",11,0.,2*TMath::Pi(),"PV nominal","Events");
  // polarimetricAcopAnglePVRefitBSMVADMHiggsQCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSMVADMHiggsQCDMC","PV refitted",11,0.,2*TMath::Pi(),"PV refitted","Events");
  // polarimetricAcopAnglePVRefitBSZNominalMVADMHiggsQCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSZNominalMVADMHiggsQCDMC","PV refitted Z nominal",11,0.,2*TMath::Pi(),"PV refitted Z nominal","Events");
  // polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsQCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsQCDMC","PV refitted with taus tracks",11,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");
  // polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMHiggsQCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMHiggsQCDMC","PV refitted with taus tracks Z nominal",11,0.,2*TMath::Pi(),"PV refitted with taus tracks Z nominal","Events");

  // polarimetricAcopAngleMVADMJetFakesQCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAngleMVADMJetFakesQCDMC","PV nominal",11,0.,2*TMath::Pi(),"PV nominal","Events");
  // polarimetricAcopAnglePVRefitBSMVADMJetFakesQCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSMVADMJetFakesQCDMC","PV refitted",11,0.,2*TMath::Pi(),"PV refitted","Events");
  // polarimetricAcopAnglePVRefitBSZNominalMVADMJetFakesQCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSZNominalMVADMJetFakesQCDMC","PV refitted Z nominal",11,0.,2*TMath::Pi(),"PV refitted Z nominal","Events");
  // polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesQCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesQCDMC","PV refitted with taus tracks",11,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");
  // polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMJetFakesQCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMJetFakesQCDMC","PV refitted with taus tracks Z nominal",11,0.,2*TMath::Pi(),"PV refitted with taus tracks Z nominal","Events");

  // polarimetricAcopAngleMVADMZTTEmbedQCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAngleMVADMZTTEmbedQCDMC","PV nominal",11,0.,2*TMath::Pi(),"PV nominal","Events");
  // polarimetricAcopAnglePVRefitBSMVADMZTTEmbedQCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSMVADMZTTEmbedQCDMC","PV refitted",11,0.,2*TMath::Pi(),"PV refitted","Events");
  // polarimetricAcopAnglePVRefitBSZNominalMVADMZTTEmbedQCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSZNominalMVADMZTTEmbedQCDMC","PV refitted Z nominal",11,0.,2*TMath::Pi(),"PV refitted Z nominal","Events");
  // polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTEmbedQCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTEmbedQCDMC","PV refitted with taus tracks",11,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");
  // polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMZTTEmbedQCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMZTTEmbedQCDMC","PV refitted with taus tracks Z nominal",11,0.,2*TMath::Pi(),"PV refitted with taus tracks Z nominal","Events");

  
  // polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled0005QCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolledQCDMC0005","PV refitted with taus tracks",5,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");
  // polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled0506QCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolledQCDMC0506","PV refitted with taus tracks",5,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");
  // polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled0607QCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolledQCDMC0607","PV refitted with taus tracks",5,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");
  // polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled0710QCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolledQCDMC0710","PV refitted with taus tracks",5,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");

  polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolledQCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolledQCDMC","PV refitted with taus tracks",180,0.,3*2*TMath::Pi(),"PV refitted with taus tracks","Events");
  
  polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolledQCDMC_DP=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolledQCDMC_DP","PV refitted with taus tracks",180,0.,3*2*TMath::Pi(),"PV refitted with taus tracks","Events");

  // polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled0005QCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolledQCDMC0005","PV refitted with taus tracks",5,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");
  // polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled0506QCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolledQCDMC0506","PV refitted with taus tracks",5,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");
  // polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled0607QCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolledQCDMC0607","PV refitted with taus tracks",5,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");
  // polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled0710QCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolledQCDMC0710","PV refitted with taus tracks",5,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");

  polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolledQCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolledQCDMC","PV refitted with taus tracks",180,0.,3*2*TMath::Pi(),"PV refitted with taus tracks","Events");
  
  polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolledQCDMC_DP=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolledQCDMC_DP","PV refitted with taus tracks",180,0.,3*2*TMath::Pi(),"PV refitted with taus tracks","Events");

  // polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled0005QCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolledQCDMC0005","PV refitted with taus tracks",5,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");
  // polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled0506QCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolledQCDMC0506","PV refitted with taus tracks",5,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");
  // polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled0607QCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolledQCDMC0607","PV refitted with taus tracks",5,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");
  // polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled0710QCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolledQCDMC0710","PV refitted with taus tracks",5,0.,2*TMath::Pi(),"PV refitted with taus tracks","Events");
  
  polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolledQCDMC=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolledQCDMC","PV refitted with taus tracks",180,0.,3*2*TMath::Pi(),"PV refitted with taus tracks","Events");

  polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolledQCDMC_DP=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolledQCDMC_DP","PV refitted with taus tracks",180,0.,3*2*TMath::Pi(),"PV refitted with taus tracks","Events");

  // polarimetricAcopAngleMVADMHiggsQCDMC=HConfig.GetTH2D(Name+"_polarimetricAcopAngleMVADMHiggsQCDMC","Higgs MVA score VS acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA score");
  // polarimetricAcopAnglePVRefitBSMVADMHiggsQCDMC=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitBSMVADMHiggsQCDMC","Higgs MVA Score VS acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA score");
  // polarimetricAcopAnglePVRefitBSZNominalMVADMHiggsQCDMC=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitBSZNominalMVADMHiggsQCDMC","Higgs MVA Score VS acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA score");
  polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsQCDMC=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsQCDMC","Higgs MVA Score VS acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA Score");

  polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsQCDMC_DP=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsQCDMC_DP","Higgs MVA Score VS acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA Score");

  //polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMHiggsQCDMC=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMHiggsQCDMC","Higgs MVA Score VS acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA score");

  // polarimetricAcopAngleMVADMJetFakesQCDMC=HConfig.GetTH2D(Name+"_polarimetricAcopAngleMVADMJetFakesQCDMC","JetFakes MVA score VS acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA score");
  // polarimetricAcopAnglePVRefitBSMVADMJetFakesQCDMC=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitBSMVADMJetFakesQCDMC","JetFakes MVA Score VS acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA score");
  // polarimetricAcopAnglePVRefitBSZNominalMVADMJetFakesQCDMC=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitBSZNominalMVADMJetFakesQCDMC","JetFakes MVA Score VS acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA score");
  polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesQCDMC=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesQCDMC","JetFakes MVA Score VS acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA Score");
  polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesQCDMC_DP=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesQCDMC_DP","JetFakes MVA Score VS acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA Score");
  //polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMJetFakesQCDMC=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMJetFakesQCDMC","JetFakes MVA Score VS acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA score");

  // polarimetricAcopAngleMVADMZTTQCDMC=HConfig.GetTH2D(Name+"_polarimetricAcopAngleMVADMZTTQCDMC","ZTT MVA score VS acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA score");
  // polarimetricAcopAnglePVRefitBSMVADMZTTQCDMC=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitBSMVADMZTTQCDMC","ZTT MVA Score VS acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA score");
  // polarimetricAcopAnglePVRefitBSZNominalMVADMZTTQCDMC=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitBSZNominalMVADMZTTQCDMC","ZTT MVA Score VS acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA score");
  polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTQCDMC=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTQCDMC","ZTT MVA Score VS acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA Score");
  polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTQCDMC_DP=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTQCDMC_DP","ZTT MVA Score VS acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA Score");
  //polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMZTTQCDMC=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMZTTQCDMC","ZTT MVA Score VS acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","MVA score");

  ShapeSystPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_ShapeSystPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");

  CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_res_j_13TeVUpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_res_j_13TeVUpPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_res_j_13TeVDownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_res_j_13TeVDownPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_scale_gg_13TeVUpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_scale_gg_13TeVUpPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_scale_gg_13TeVDownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_scale_gg_13TeVDownPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
	  
  jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");


  PrefiringUpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_PrefiringUpPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  PrefiringDownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_PrefiringDownPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  
  ShapeSystPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_ShapeSystPVRefitWithTracksBSWfakesHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");

  polarimetricAcopAnglePVRefitWithTracksBSMVADMWfakesHiggs=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMWfakesHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSWfakesHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSWfakesHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSWfakesHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSWfakesHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSWfakesHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSWfakesHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_res_j_13TeVUpPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_res_j_13TeVUpPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_res_j_13TeVDownPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_res_j_13TeVDownPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSWfakesHiggs","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSWfakesHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_scale_gg_13TeVUpPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_scale_gg_13TeVUpPVRefitWithTracksBSWfakesHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSWfakesHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSWfakesHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSWfakesHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_scale_gg_13TeVDownPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_scale_gg_13TeVDownPVRefitWithTracksBSWfakesHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSWfakesHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSWfakesHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  PrefiringUpPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_PrefiringUpPVRefitWithTracksBSWfakesHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  PrefiringDownPVRefitWithTracksBSWfakesHiggs=HConfig.GetTH2D(Name+"_PrefiringDownPVRefitWithTracksBSWfakesHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");

  
  ShapeSystPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_ShapeSystPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");

  CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_res_j_13TeVUpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_res_j_13TeVUpPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_res_j_13TeVDownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_res_j_13TeVDownPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_scale_gg_13TeVUpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_gg_13TeVUpPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_scale_gg_13TeVDownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_gg_13TeVDownPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");

  jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  
  PrefiringUpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_PrefiringUpPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  PrefiringDownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_PrefiringDownPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  
  
  ShapeSystPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_ShapeSystPVRefitWithTracksBSWfakesJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");

  polarimetricAcopAnglePVRefitWithTracksBSMVADMWfakesJetFakes=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMWfakesJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSWfakesJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSWfakesJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSWfakesJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSWfakesJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSWfakesJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSWfakesJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_res_j_13TeVUpPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_res_j_13TeVUpPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_res_j_13TeVDownPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_res_j_13TeVDownPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSWfakesJetFakes","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSWfakesJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_scale_gg_13TeVUpPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_gg_13TeVUpPVRefitWithTracksBSWfakesJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSWfakesJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSWfakesJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSWfakesJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_scale_gg_13TeVDownPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_scale_gg_13TeVDownPVRefitWithTracksBSWfakesJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSWfakesJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSWfakesJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  PrefiringUpPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_PrefiringUpPVRefitWithTracksBSWfakesJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  PrefiringDownPVRefitWithTracksBSWfakesJetFakes=HConfig.GetTH2D(Name+"_PrefiringDownPVRefitWithTracksBSWfakesJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");


  ShapeSystPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_ShapeSystPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");

  CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_res_j_13TeVUpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_res_j_13TeVUpPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_res_j_13TeVDownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_res_j_13TeVDownPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_scale_gg_13TeVUpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_scale_gg_13TeVUpPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_scale_gg_13TeVDownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_scale_gg_13TeVDownPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");

  PrefiringUpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_PrefiringUpPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  PrefiringDownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_PrefiringDownPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  
  ShapeSystPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_ShapeSystPVRefitWithTracksBSWfakesZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");

  polarimetricAcopAnglePVRefitWithTracksBSMVADMWfakesZTT=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMWfakesZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSWfakesZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSWfakesZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSWfakesZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSWfakesZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSWfakesZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSWfakesZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_res_j_13TeVUpPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_res_j_13TeVUpPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_res_j_13TeVDownPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_res_j_13TeVDownPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSWfakesZTT","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSWfakesZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_scale_gg_13TeVUpPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_scale_gg_13TeVUpPVRefitWithTracksBSWfakesZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSWfakesZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSWfakesZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSWfakesZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_scale_gg_13TeVDownPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_scale_gg_13TeVDownPVRefitWithTracksBSWfakesZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSWfakesZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSWfakesZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  PrefiringUpPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_PrefiringUpPVRefitWithTracksBSWfakesZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  PrefiringDownPVRefitWithTracksBSWfakesZTT=HConfig.GetTH2D(Name+"_PrefiringDownPVRefitWithTracksBSWfakesZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");

  ttbarcontaminationHiggs=HConfig.GetTH2D(Name+"_ttbarcontaminationHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  ttbarcontaminationJetFakes=HConfig.GetTH2D(Name+"_ttbarcontaminationJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  ttbarcontaminationZTT=HConfig.GetTH2D(Name+"_ttbarcontaminationZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
 
  ttbarcontaminationHiggs_DP=HConfig.GetTH2D(Name+"_ttbarcontaminationHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  ttbarcontaminationJetFakes_DP=HConfig.GetTH2D(Name+"_ttbarcontaminationJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  ttbarcontaminationZTT_DP=HConfig.GetTH2D(Name+"_ttbarcontaminationZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");

  ttbarcontaminationWfakesHiggs=HConfig.GetTH2D(Name+"_ttbarcontaminationWfakesHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  ttbarcontaminationWfakesJetFakes=HConfig.GetTH2D(Name+"_ttbarcontaminationWfakesJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  ttbarcontaminationWfakesZTT=HConfig.GetTH2D(Name+"_ttbarcontaminationWfakesZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
 
  ttbarcontaminationWfakesHiggs_DP=HConfig.GetTH2D(Name+"_ttbarcontaminationWfakesHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  ttbarcontaminationWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_ttbarcontaminationWfakesJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  ttbarcontaminationWfakesZTT_DP=HConfig.GetTH2D(Name+"_ttbarcontaminationWfakesZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");

  jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSHiggsQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSHiggsQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSHiggsQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSHiggsQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSHiggsQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSHiggsQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSHiggsQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSHiggsQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSHiggsQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSHiggsQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSHiggsQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSHiggsQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");

  jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSJetFakesQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSJetFakesQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSJetFakesQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSJetFakesQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSJetFakesQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSJetFakesQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSJetFakesQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSJetFakesQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSJetFakesQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSJetFakesQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSJetFakesQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSJetFakesQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");

  jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSZTTQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSZTTQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSZTTQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSZTTQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSZTTQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSZTTQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSZTTQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSZTTQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSZTTQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSZTTQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSZTTQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSZTTQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");

  //-------------------------DP
   
  ShapeSystPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_ShapeSystPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");

  CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
   
  CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
   
  CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
   
  CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
   
  // CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
   
  // CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
   
  // CMS_res_j_13TeVUpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_res_j_13TeVUpPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
   
  // CMS_res_j_13TeVDownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_res_j_13TeVDownPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
   
  // CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
   
  // CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
   
  // CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
   
  // CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
   
  // CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
   
  // CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events"); CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events"); CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
   
  // CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
   
  // CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_scale_gg_13TeVUpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_gg_13TeVUpPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_scale_gg_13TeVDownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_gg_13TeVDownPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
 	  
  jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");

  PrefiringUpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_PrefiringUpPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  PrefiringDownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_PrefiringDownPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
   
  ShapeSystPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_ShapeSystPVRefitWithTracksBSWfakesHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");

  polarimetricAcopAnglePVRefitWithTracksBSMVADMWfakesHiggs_DP=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMWfakesHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_res_j_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_res_j_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_res_j_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_res_j_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_scale_gg_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_gg_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_scale_gg_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_scale_gg_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  PrefiringUpPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_PrefiringUpPVRefitWithTracksBSWfakesHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  PrefiringDownPVRefitWithTracksBSWfakesHiggs_DP=HConfig.GetTH2D(Name+"_PrefiringDownPVRefitWithTracksBSWfakesHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
 
  ShapeSystPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_ShapeSystPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");

  CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_res_j_13TeVUpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_res_j_13TeVUpPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_res_j_13TeVDownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_res_j_13TeVDownPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_scale_gg_13TeVUpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_gg_13TeVUpPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_scale_gg_13TeVDownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_gg_13TeVDownPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
 
  jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");

  PrefiringUpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_PrefiringUpPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  PrefiringDownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_PrefiringDownPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
   
  ShapeSystPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_ShapeSystPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");

  polarimetricAcopAnglePVRefitWithTracksBSMVADMWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMWfakesJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_res_j_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_res_j_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_res_j_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_res_j_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_scale_gg_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_gg_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_scale_gg_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_scale_gg_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  PrefiringUpPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_PrefiringUpPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  PrefiringDownPVRefitWithTracksBSWfakesJetFakes_DP=HConfig.GetTH2D(Name+"_PrefiringDownPVRefitWithTracksBSWfakesJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
 
  ShapeSystPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_ShapeSystPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");

  CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_res_j_13TeVUpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_res_j_13TeVUpPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_res_j_13TeVDownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_res_j_13TeVDownPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_scale_gg_13TeVUpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_gg_13TeVUpPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_scale_gg_13TeVDownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_gg_13TeVDownPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  PrefiringUpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_PrefiringUpPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  PrefiringDownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_PrefiringDownPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");

  ShapeSystPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_ShapeSystPVRefitWithTracksBSWfakesZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");

  polarimetricAcopAnglePVRefitWithTracksBSMVADMWfakesZTT_DP=HConfig.GetTH2D(Name+"_polarimetricAcopAnglePVRefitWithTracksBSMVADMWfakesZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSWfakesZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSWfakesZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSWfakesZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSWfakesZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSWfakesZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSWfakesZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_res_j_13TeVUpPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_res_j_13TeVUpPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_res_j_13TeVDownPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_res_j_13TeVDownPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSWfakesZTT_DP","acop angle",5,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSWfakesZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_scale_gg_13TeVUpPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_gg_13TeVUpPVRefitWithTracksBSWfakesZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSWfakesZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSWfakesZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSWfakesZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_scale_gg_13TeVDownPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_scale_gg_13TeVDownPVRefitWithTracksBSWfakesZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSWfakesZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSWfakesZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");

  PrefiringUpPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_PrefiringUpPVRefitWithTracksBSWfakesZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  PrefiringDownPVRefitWithTracksBSWfakesZTT_DP=HConfig.GetTH2D(Name+"_PrefiringDownPVRefitWithTracksBSWfakesZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  
  jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  
  //-- 
  jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
 
  jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSHiggsQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSHiggsQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSHiggsQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSHiggsQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSHiggsQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSHiggsQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSHiggsQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSHiggsQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSHiggsQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSHiggsQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSHiggsQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSHiggsQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
 
  jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSJetFakesQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSJetFakesQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
 
  // jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSJetFakesQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSJetFakesQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSJetFakesQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSJetFakesQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSJetFakesQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSJetFakesQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSJetFakesQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSJetFakesQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSJetFakesQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSJetFakesQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
 
  jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
 
  //BUG_AC 
  jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSZTTQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSZTTQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSZTTQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSZTTQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
   
  // jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSZTTQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSZTTQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  // jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSZTTQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSZTTQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSZTTQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSZTTQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSZTTQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSZTTQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");


  jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSHiggs=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSHiggs","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");

  jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSJetFakes=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSJetFakes","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  
  jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSZTT=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSZTT","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  

  jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSHiggsQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSHiggsQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSHiggsQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSHiggsQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSHiggsQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSHiggsQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSHiggsQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSHiggsQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSHiggsQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSHiggsQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSHiggsQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSHiggsQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSHiggsQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSHiggsQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSHiggsQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSHiggsQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSHiggsQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSHiggsQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSHiggsQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSHiggsQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSHiggsQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSHiggsQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSHiggsQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSHiggsQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");

  jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSJetFakesQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSJetFakesQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSJetFakesQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSJetFakesQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSJetFakesQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSJetFakesQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSJetFakesQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSJetFakesQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSJetFakesQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSJetFakesQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSJetFakesQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSJetFakesQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSJetFakesQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSJetFakesQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSJetFakesQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSJetFakesQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSJetFakesQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSJetFakesQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSJetFakesQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSJetFakesQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSJetFakesQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSJetFakesQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSJetFakesQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSJetFakesQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  
  jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSZTTQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSZTTQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSZTTQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSZTTQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSZTTQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSZTTQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSZTTQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSZTTQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSZTTQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSZTTQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSZTTQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSZTTQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSZTTQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSZTTQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSZTTQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSZTTQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSZTTQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSZTTQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSZTTQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSZTTQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSZTTQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSZTTQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSZTTQCDMC=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSZTTQCDMC","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  
  

  jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSHiggs_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSHiggs_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");

  jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSJetFakes_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSJetFakes_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  
  jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSZTT_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSZTT_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  

  jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSHiggsQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSHiggsQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSHiggsQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSHiggsQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSHiggsQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSHiggsQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSHiggsQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSHiggsQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSHiggsQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSHiggsQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSHiggsQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSHiggsQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSHiggsQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSHiggsQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSHiggsQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSHiggsQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSHiggsQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSHiggsQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSHiggsQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSHiggsQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSHiggsQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSHiggsQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSHiggsQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSHiggsQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");

  jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSJetFakesQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSJetFakesQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSJetFakesQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSJetFakesQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSJetFakesQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSJetFakesQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSJetFakesQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSJetFakesQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSJetFakesQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSJetFakesQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSJetFakesQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSJetFakesQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSJetFakesQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSJetFakesQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSJetFakesQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSJetFakesQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSJetFakesQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSJetFakesQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSJetFakesQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSJetFakesQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSJetFakesQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSJetFakesQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSJetFakesQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSJetFakesQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  
  jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSZTTQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSZTTQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSZTTQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSZTTQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSZTTQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSZTTQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSZTTQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSZTTQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSZTTQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSZTTQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSZTTQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSZTTQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSZTTQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSZTTQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSZTTQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSZTTQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSZTTQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSZTTQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSZTTQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSZTTQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSZTTQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSZTTQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSZTTQCDMC_DP=HConfig.GetTH2D(Name+"_jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSZTTQCDMC_DP","acop angle",60,0.,2*TMath::Pi(),7,0.3,1,"Acop angle","Events");
  

  Selection::ConfigureHistograms();   //   do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);  // do not remove
  
}

void  HCPTauTau::Store_ExtraDist(){

  //every new histo should be addedd to Extradist1d vector, just push it back;
  // Extradist1d.push_back(&Tau1PT);
  // Extradist1d.push_back(&Tau1E);
  // Extradist1d.push_back(&Tau1Mass);
  // Extradist1d.push_back(&Tau1Phi);
  // Extradist1d.push_back(&Tau1Eta);
  // Extradist1d.push_back(&Tau1dz);
  // Extradist1d.push_back(&Tau1HPSDecayMode);
  // Extradist1d.push_back(&Tau1MVADecayMode);
  // Extradist1d.push_back(&Tau1GenMatch);
  // Extradist1d.push_back(&Tau2PT);
  // Extradist1d.push_back(&Tau2E);
  // Extradist1d.push_back(&Tau2Mass);
  // Extradist1d.push_back(&Tau2Phi);
  // Extradist1d.push_back(&Tau2Eta);
  // Extradist1d.push_back(&Tau2dz);
  // Extradist1d.push_back(&Tau2HPSDecayMode);
  // Extradist1d.push_back(&Tau2MVADecayMode);
  // Extradist1d.push_back(&Tau2GenMatch);
  
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
  // Extradist1d.push_back(&ExtraLeptonVeto);

  // Extradist1d.push_back(&dRTauTau);
  // Extradist1d.push_back(&TauTauVisMass);
  // //  Extradist1d.push_back(&TauTauTruthMass);
  // Extradist1d.push_back(&TauTauFullMass);
  
  // Extradist1d.push_back(&QCDShape);
  // Extradist1d.push_back(&NQCD);
  // // Extradist1d.push_back(&TauTauFullMass_B);
  // // Extradist1d.push_back(&TauTauFullMass_C);
  // // Extradist1d.push_back(&TauTauFullMass_D);

  // Extradist1d.push_back(&NFFData);
  // Extradist1d.push_back(&NFFLeadMC);
  Extradist2d.push_back(&WfakesHiggs);
  Extradist2d.push_back(&WfakesJetFakes);
  Extradist2d.push_back(&WfakesZTT);

  Extradist2d.push_back(&WfakesHiggs_DP);
  Extradist2d.push_back(&WfakesJetFakes_DP);
  Extradist2d.push_back(&WfakesZTT_DP);

  // Extradist1d.push_back(&MET);
  // Extradist1d.push_back(&METphi);
  // Extradist1d.push_back(&PUPPImet);
  // Extradist1d.push_back(&PUPPImetphi);
  // Extradist1d.push_back(&PUPPImetcorr);
  // Extradist1d.push_back(&PUPPImetcorrphi);
  // Extradist1d.push_back(&TransverseMass);

  // Extradist1d.push_back(&NPrimeVtx);
  // Extradist1d.push_back(&NPU);
  // Extradist1d.push_back(&RHO);


  //Extradist1d.push_back(&NbJets);

  // Extradist1d.push_back(&HPtVis);
  // Extradist1d.push_back(&TauTauVisPT);
  //Extradist1d.push_back(&TauTauFullPT);
  //Extradist1d.push_back(&Mjj);
  // Extradist1d.push_back(&dijetpt);
  // Extradist1d.push_back(&dijetphi);
  //Extradist1d.push_back(&jdeta);
  // Extradist1d.push_back(&jdphi);
  // Extradist1d.push_back(&jpt_2);
  // Extradist1d.push_back(&jeta_2);
  // Extradist1d.push_back(&jphi_2);
  //Extradist1d.push_back(&jpt_1);
  // Extradist1d.push_back(&jeta_1);
  // Extradist1d.push_back(&jphi_1);

  // Extradist1d.push_back(&ResolPullTauTauFroma1a1MZMomentum);
  // Extradist1d.push_back(&ResolPullTauminusFroma1a1MZMomentum);

  // Extradist1d.push_back(&ResolPullTauplusFroma1a1MZMomentum);
  // Extradist1d.push_back(&ResolPullTauFroma1a1MZMomentum);

  // Extradist1d.push_back(&ResolPullXVtxIna1a1);
  // Extradist1d.push_back(&ResolPullYVtxIna1a1);
  // Extradist1d.push_back(&ResolPullZVtxIna1a1);

  // Extradist1d.push_back(&tauminusa1a1MomentumPairConstraint);
  // Extradist1d.push_back(&tauplusa1a1MomentumPairConstraint);
  
  // Extradist1d.push_back(&polarimetricAcopAngle);
 
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBS);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitBS);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBSZNominal);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSZNominal);

  // Extradist1d.push_back(&polarimetricAcopAngleMVADM);

  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBSMVADM);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSMVADM);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBSZNominalMVADM);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSZNominalMVADM);

  //  Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADM);
  //Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADM);

  // Extradist1d.push_back(&PVXResol);
  // Extradist1d.push_back(&PVXNoBSResol);
  // Extradist1d.push_back(&PVXBSResol);

  // Extradist1d.push_back(&PVYResol);
  // Extradist1d.push_back(&PVYNoBSResol);
  // Extradist1d.push_back(&PVYBSResol);
  
  // Extradist1d.push_back(&PVZResol);
  // Extradist1d.push_back(&PVZNoBSResol);
  // Extradist1d.push_back(&PVZBSResol); 
 
  // Extradist1d.push_back(&PVXNoBSOnlyResol);
  // Extradist1d.push_back(&PVXBSOnlyResol);
  
  // Extradist1d.push_back(&PVYNoBSOnlyResol);
  // Extradist1d.push_back(&PVYBSOnlyResol);
  
  // Extradist1d.push_back(&PVZNoBSOnlyResol);
  // Extradist1d.push_back(&PVZBSOnlyResol);

  // Extradist1d.push_back(&polarimetricAcopAngleTruthA1);
  
  // Extradist1d.push_back(&test);  
  
  // Extradist1d.push_back(&PurityDM);
  // Extradist1d.push_back(&PurityNewMVA);

  // Extradist1d.push_back(&TauPxResPull);
  // Extradist1d.push_back(&TauPyResPull);
  // Extradist1d.push_back(&TauPzResPull);

  // Extradist1d.push_back(&TauPxResPullMVA);
  // Extradist1d.push_back(&TauPyResPullMVA);
  // Extradist1d.push_back(&TauPzResPullMVA);

  // Extradist1d.push_back(&HiggsBDTScore);
  // Extradist1d.push_back(&JetFakesBDTScore);
  // Extradist1d.push_back(&ZTTBDTScore);


  Extradist1d.push_back(&Tau1PTa1a1);
  // Extradist1d.push_back(&Tau1Ea1a1);
  // Extradist1d.push_back(&Tau1Massa1a1);
  // Extradist1d.push_back(&Tau1Phia1a1);
  // Extradist1d.push_back(&Tau1Etaa1a1);
  // Extradist1d.push_back(&Tau1dza1a1);

  // Extradist1d.push_back(&Tau2PTa1a1);
  // Extradist1d.push_back(&Tau2Ea1a1);
  // Extradist1d.push_back(&Tau2Massa1a1);
  // Extradist1d.push_back(&Tau2Phia1a1);
  // Extradist1d.push_back(&Tau2Etaa1a1);
  // Extradist1d.push_back(&Tau2dza1a1);
  
  Extradist1d.push_back(&TauTauVisMassa1a1);
  Extradist1d.push_back(&TauTauFullMassa1a1);
  Extradist1d.push_back(&TauTauVisPTa1a1);
  Extradist1d.push_back(&TauTauFullPTa1a1);
  Extradist1d.push_back(&Mjja1a1);
  // Extradist1d.push_back(&dijetpta1a1);
  // Extradist1d.push_back(&dijetphia1a1);
  Extradist1d.push_back(&jdetaa1a1);
  // Extradist1d.push_back(&jdphia1a1);
  // Extradist1d.push_back(&jpt_2a1a1);
  // Extradist1d.push_back(&jeta_2a1a1);
  // Extradist1d.push_back(&jphi_2a1a1);
  Extradist1d.push_back(&jpt_1a1a1);
  // Extradist1d.push_back(&jeta_1a1a1);
  // Extradist1d.push_back(&jphi_1a1a1);
  
  Extradist1d.push_back(&PUPPImetcorra1a1);
  // Extradist1d.push_back(&PUPPImetcorrphia1a1);

  Extradist1d.push_back(&NbJetsa1a1);

  Extradist1d.push_back(&HiggsBDTScorea1a1);
  Extradist1d.push_back(&JetFakesBDTScorea1a1);
  Extradist1d.push_back(&ZTTBDTScorea1a1);
   
  Extradist1d.push_back(&HiggsBDTScorea1a1_DP);
  Extradist1d.push_back(&JetFakesBDTScorea1a1_DP);
  Extradist1d.push_back(&ZTTBDTScorea1a1_DP);
   
  Extradist1d.push_back(&Tau1PTa1a1QCDMC);
  Extradist1d.push_back(&TauTauVisMassa1a1QCDMC);
  Extradist1d.push_back(&TauTauFullMassa1a1QCDMC);
  Extradist1d.push_back(&TauTauVisPTa1a1QCDMC);
  Extradist1d.push_back(&TauTauFullPTa1a1QCDMC);
  Extradist1d.push_back(&Mjja1a1QCDMC);
  Extradist1d.push_back(&jdetaa1a1QCDMC);
  Extradist1d.push_back(&jpt_1a1a1QCDMC);
  Extradist1d.push_back(&PUPPImetcorra1a1QCDMC);
  Extradist1d.push_back(&NbJetsa1a1QCDMC);
  // Extradist1d.push_back(&polarimetricAcopAngleMVADMHiggs);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSMVADMHiggs);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSZNominalMVADMHiggs);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggs);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMHiggs);

  // Extradist1d.push_back(&polarimetricAcopAngleMVADMJetFakes);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSMVADMJetFakes);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSZNominalMVADMJetFakes);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakes);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMJetFakes);

  // Extradist1d.push_back(&polarimetricAcopAngleMVADMZTTEmbed);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSMVADMZTTEmbed);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSZNominalMVADMZTTEmbed);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTEmbed);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMZTTEmbed);

  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled0005);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled0506);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled0607);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled0710);
  
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled0005);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled0506);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled0607);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled0710);

  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled0005);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled0506);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled0607);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled0710);

  Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADM);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMQCDMC);

  Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled);

  Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADM_DP);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMQCDMC_DP);

  Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled_DP);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled_DP);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled_DP);

  // Extradist2d.push_back(&polarimetricAcopAngleMVADMHiggs);
  // Extradist2d.push_back(&polarimetricAcopAnglePVRefitBSMVADMHiggs);
  // Extradist2d.push_back(&polarimetricAcopAnglePVRefitBSZNominalMVADMHiggs);
  Extradist2d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggs);
  Extradist2d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggs_DP);
  //  Extradist2d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMHiggs);

  //Extradist2d.push_back(&polarimetricAcopAngleMVADMJetFakes);
  //Extradist2d.push_back(&polarimetricAcopAnglePVRefitBSMVADMJetFakes);
  //Extradist2d.push_back(&polarimetricAcopAnglePVRefitBSZNominalMVADMJetFakes);
  Extradist2d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakes);
  Extradist2d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakes_DP);
  //  Extradist2d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMJetFakes);

  //Extradist2d.push_back(&polarimetricAcopAngleMVADMZTT);
  //Extradist2d.push_back(&polarimetricAcopAnglePVRefitBSMVADMZTT);
  //Extradist2d.push_back(&polarimetricAcopAnglePVRefitBSZNominalMVADMZTT);
  Extradist2d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMZTT);
  Extradist2d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMZTT_DP);
  
  //  Extradist2d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMZTT);

  Extradist2d.push_back(&ShapeSystPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_res_j_13TeVUpPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_res_j_13TeVDownPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&CMS_scale_gg_13TeVUpPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&CMS_scale_gg_13TeVDownPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSHiggs);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSHiggs);
  
  Extradist2d.push_back(&PrefiringUpPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&PrefiringDownPVRefitWithTracksBSHiggs);
  
  Extradist2d.push_back(&ShapeSystPVRefitWithTracksBSWfakesHiggs);
  Extradist2d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMWfakesHiggs);
  Extradist2d.push_back(&CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSWfakesHiggs);
  Extradist2d.push_back(&CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSWfakesHiggs);
  Extradist2d.push_back(&CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSWfakesHiggs);
  Extradist2d.push_back(&CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSWfakesHiggs);
  Extradist2d.push_back(&CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSWfakesHiggs);
  Extradist2d.push_back(&CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_res_j_13TeVUpPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_res_j_13TeVDownPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSWfakesHiggs);
  // Extradist2d.push_back(&CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSWfakesHiggs);
  Extradist2d.push_back(&CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesHiggs);
  Extradist2d.push_back(&CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSWfakesHiggs);
  Extradist2d.push_back(&CMS_scale_gg_13TeVUpPVRefitWithTracksBSWfakesHiggs);
  Extradist2d.push_back(&CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSWfakesHiggs);
  Extradist2d.push_back(&CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSWfakesHiggs);
  Extradist2d.push_back(&CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesHiggs);
  Extradist2d.push_back(&CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSWfakesHiggs);
  Extradist2d.push_back(&CMS_scale_gg_13TeVDownPVRefitWithTracksBSWfakesHiggs);
  Extradist2d.push_back(&CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSWfakesHiggs);
  Extradist2d.push_back(&CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSWfakesHiggs);
  Extradist2d.push_back(&PrefiringUpPVRefitWithTracksBSWfakesHiggs);
  Extradist2d.push_back(&PrefiringDownPVRefitWithTracksBSWfakesHiggs);
  
  Extradist2d.push_back(&ShapeSystPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_res_j_13TeVUpPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_res_j_13TeVDownPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&CMS_scale_gg_13TeVUpPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&CMS_scale_gg_13TeVDownPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSJetFakes);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSJetFakes);
  
  Extradist2d.push_back(&PrefiringUpPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&PrefiringDownPVRefitWithTracksBSJetFakes);
  
  Extradist2d.push_back(&ShapeSystPVRefitWithTracksBSWfakesJetFakes);
  Extradist2d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMWfakesJetFakes);
  Extradist2d.push_back(&CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSWfakesJetFakes);
  Extradist2d.push_back(&CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSWfakesJetFakes);
  Extradist2d.push_back(&CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSWfakesJetFakes);
  Extradist2d.push_back(&CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSWfakesJetFakes);
  Extradist2d.push_back(&CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSWfakesJetFakes);
  Extradist2d.push_back(&CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_res_j_13TeVUpPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_res_j_13TeVDownPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSWfakesJetFakes);
  // Extradist2d.push_back(&CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSWfakesJetFakes);
  Extradist2d.push_back(&CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesJetFakes);
  Extradist2d.push_back(&CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSWfakesJetFakes);
  Extradist2d.push_back(&CMS_scale_gg_13TeVUpPVRefitWithTracksBSWfakesJetFakes);
  Extradist2d.push_back(&CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSWfakesJetFakes);
  Extradist2d.push_back(&CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSWfakesJetFakes);
  Extradist2d.push_back(&CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesJetFakes);
  Extradist2d.push_back(&CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSWfakesJetFakes);
  Extradist2d.push_back(&CMS_scale_gg_13TeVDownPVRefitWithTracksBSWfakesJetFakes);
  Extradist2d.push_back(&CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSWfakesJetFakes);
  Extradist2d.push_back(&CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSWfakesJetFakes);
  Extradist2d.push_back(&PrefiringUpPVRefitWithTracksBSWfakesJetFakes);
  Extradist2d.push_back(&PrefiringDownPVRefitWithTracksBSWfakesJetFakes);
  
  Extradist2d.push_back(&ShapeSystPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_res_j_13TeVUpPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_res_j_13TeVDownPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&CMS_scale_gg_13TeVUpPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&CMS_scale_gg_13TeVDownPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSZTT);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSZTT);
  
  Extradist2d.push_back(&PrefiringUpPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&PrefiringDownPVRefitWithTracksBSZTT);
  
  Extradist2d.push_back(&ShapeSystPVRefitWithTracksBSWfakesZTT);
  Extradist2d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMWfakesZTT);
  Extradist2d.push_back(&CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSWfakesZTT);
  Extradist2d.push_back(&CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSWfakesZTT);
  Extradist2d.push_back(&CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSWfakesZTT);
  Extradist2d.push_back(&CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSWfakesZTT);
  Extradist2d.push_back(&CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSWfakesZTT);
  Extradist2d.push_back(&CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_res_j_13TeVUpPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_res_j_13TeVDownPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSWfakesZTT);
  // Extradist2d.push_back(&CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSWfakesZTT);
  Extradist2d.push_back(&CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesZTT);
  Extradist2d.push_back(&CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSWfakesZTT);
  Extradist2d.push_back(&CMS_scale_gg_13TeVUpPVRefitWithTracksBSWfakesZTT);
  Extradist2d.push_back(&CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSWfakesZTT);
  Extradist2d.push_back(&CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSWfakesZTT);
  Extradist2d.push_back(&CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesZTT);
  Extradist2d.push_back(&CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSWfakesZTT);
  Extradist2d.push_back(&CMS_scale_gg_13TeVDownPVRefitWithTracksBSWfakesZTT);
  Extradist2d.push_back(&CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSWfakesZTT);
  Extradist2d.push_back(&CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSWfakesZTT);
  Extradist2d.push_back(&PrefiringUpPVRefitWithTracksBSWfakesZTT);
  Extradist2d.push_back(&PrefiringDownPVRefitWithTracksBSWfakesZTT);


  Extradist2d.push_back(&ttbarcontaminationHiggs);
  Extradist2d.push_back(&ttbarcontaminationJetFakes);
  Extradist2d.push_back(&ttbarcontaminationZTT);

  Extradist2d.push_back(&ttbarcontaminationHiggs_DP);
  Extradist2d.push_back(&ttbarcontaminationJetFakes_DP);
  Extradist2d.push_back(&ttbarcontaminationZTT_DP);

  Extradist2d.push_back(&ttbarcontaminationWfakesHiggs);
  Extradist2d.push_back(&ttbarcontaminationWfakesJetFakes);
  Extradist2d.push_back(&ttbarcontaminationWfakesZTT);

  Extradist2d.push_back(&ttbarcontaminationWfakesHiggs_DP);
  Extradist2d.push_back(&ttbarcontaminationWfakesJetFakes_DP);
  Extradist2d.push_back(&ttbarcontaminationWfakesZTT_DP);

  Extradist1d.push_back(&HiggsBDTScorea1a1QCDMC);
  Extradist1d.push_back(&JetFakesBDTScorea1a1QCDMC);
  Extradist1d.push_back(&ZTTBDTScorea1a1QCDMC);
  
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolledQCDMC);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolledQCDMC);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolledQCDMC);
  Extradist2d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsQCDMC);
  Extradist2d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesQCDMC);
  Extradist2d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTQCDMC);

  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSHiggsQCDMC);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSHiggsQCDMC);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSHiggsQCDMC);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSHiggsQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSHiggsQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSHiggsQCDMC);

  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSJetFakesQCDMC);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSJetFakesQCDMC);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSJetFakesQCDMC);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSJetFakesQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSJetFakesQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSJetFakesQCDMC);

  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSZTTQCDMC);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSZTTQCDMC);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSZTTQCDMC);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSZTTQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSZTTQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSZTTQCDMC);

  //DP-------------------
   
  Extradist2d.push_back(&ShapeSystPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_res_j_13TeVUpPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_res_j_13TeVDownPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&CMS_scale_gg_13TeVUpPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&CMS_scale_gg_13TeVDownPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSHiggs_DP);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSHiggs_DP);
  
  Extradist2d.push_back(&PrefiringUpPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&PrefiringDownPVRefitWithTracksBSHiggs_DP);
   
  Extradist2d.push_back(&ShapeSystPVRefitWithTracksBSWfakesHiggs_DP);
  Extradist2d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMWfakesHiggs_DP);
  Extradist2d.push_back(&CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP);
  Extradist2d.push_back(&CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP);
  Extradist2d.push_back(&CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP);
  Extradist2d.push_back(&CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP);
  Extradist2d.push_back(&CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP);
  Extradist2d.push_back(&CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_res_j_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_res_j_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP);
  // Extradist2d.push_back(&CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP);
  Extradist2d.push_back(&CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP);
  Extradist2d.push_back(&CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP);
  Extradist2d.push_back(&CMS_scale_gg_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP);
  Extradist2d.push_back(&CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP);
  Extradist2d.push_back(&CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP);
  Extradist2d.push_back(&CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP);
  Extradist2d.push_back(&CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP);
  Extradist2d.push_back(&CMS_scale_gg_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP);
  Extradist2d.push_back(&CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP);
  Extradist2d.push_back(&CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP);
  Extradist2d.push_back(&PrefiringUpPVRefitWithTracksBSWfakesHiggs_DP);
  Extradist2d.push_back(&PrefiringDownPVRefitWithTracksBSWfakesHiggs_DP);
   
  Extradist2d.push_back(&ShapeSystPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_res_j_13TeVUpPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_res_j_13TeVDownPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&CMS_scale_gg_13TeVUpPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&CMS_scale_gg_13TeVDownPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSJetFakes_DP);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSJetFakes_DP);

  Extradist2d.push_back(&PrefiringUpPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&PrefiringDownPVRefitWithTracksBSJetFakes_DP);

  Extradist2d.push_back(&ShapeSystPVRefitWithTracksBSWfakesJetFakes_DP);
  Extradist2d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMWfakesJetFakes_DP);
  Extradist2d.push_back(&CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP);
  Extradist2d.push_back(&CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP);
  Extradist2d.push_back(&CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP);
  Extradist2d.push_back(&CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP);
  Extradist2d.push_back(&CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP);
  Extradist2d.push_back(&CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_res_j_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_res_j_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP);
  // Extradist2d.push_back(&CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP);
  Extradist2d.push_back(&CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP);
  Extradist2d.push_back(&CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP);
  Extradist2d.push_back(&CMS_scale_gg_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP);
  Extradist2d.push_back(&CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP);
  Extradist2d.push_back(&CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP);
  Extradist2d.push_back(&CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP);
  Extradist2d.push_back(&CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP);
  Extradist2d.push_back(&CMS_scale_gg_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP);
  Extradist2d.push_back(&CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP);
  Extradist2d.push_back(&CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP);
  Extradist2d.push_back(&PrefiringUpPVRefitWithTracksBSWfakesJetFakes_DP);
  Extradist2d.push_back(&PrefiringDownPVRefitWithTracksBSWfakesJetFakes_DP);
   
  Extradist2d.push_back(&ShapeSystPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_res_j_13TeVUpPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_res_j_13TeVDownPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&CMS_scale_gg_13TeVUpPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&CMS_scale_gg_13TeVDownPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSZTT_DP);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSZTT_DP);
   
  Extradist2d.push_back(&PrefiringUpPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&PrefiringDownPVRefitWithTracksBSZTT_DP);
   
  Extradist2d.push_back(&ShapeSystPVRefitWithTracksBSWfakesZTT_DP);
  Extradist2d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMWfakesZTT_DP);
  Extradist2d.push_back(&CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSWfakesZTT_DP);
  Extradist2d.push_back(&CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSWfakesZTT_DP);
  Extradist2d.push_back(&CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSWfakesZTT_DP);
  Extradist2d.push_back(&CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSWfakesZTT_DP);
  Extradist2d.push_back(&CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSWfakesZTT_DP);
  Extradist2d.push_back(&CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_res_j_13TeVUpPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_res_j_13TeVDownPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSWfakesZTT_DP);
  // Extradist2d.push_back(&CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSWfakesZTT_DP);
  Extradist2d.push_back(&CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesZTT_DP);
  Extradist2d.push_back(&CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSWfakesZTT_DP);
  Extradist2d.push_back(&CMS_scale_gg_13TeVUpPVRefitWithTracksBSWfakesZTT_DP);
  Extradist2d.push_back(&CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSWfakesZTT_DP);
  Extradist2d.push_back(&CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSWfakesZTT_DP);
  Extradist2d.push_back(&CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesZTT_DP);
  Extradist2d.push_back(&CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSWfakesZTT_DP);
  Extradist2d.push_back(&CMS_scale_gg_13TeVDownPVRefitWithTracksBSWfakesZTT_DP);
  Extradist2d.push_back(&CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSWfakesZTT_DP);
  Extradist2d.push_back(&CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSWfakesZTT_DP);
  Extradist2d.push_back(&PrefiringUpPVRefitWithTracksBSWfakesZTT_DP);
  Extradist2d.push_back(&PrefiringDownPVRefitWithTracksBSWfakesZTT_DP);
 
  Extradist1d.push_back(&HiggsBDTScorea1a1QCDMC_DP);
  Extradist1d.push_back(&JetFakesBDTScorea1a1QCDMC_DP);
  Extradist1d.push_back(&ZTTBDTScorea1a1QCDMC_DP);
   
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolledQCDMC_DP);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolledQCDMC_DP);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolledQCDMC_DP);
  Extradist2d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsQCDMC_DP);
  Extradist2d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesQCDMC_DP);
  Extradist2d.push_back(&polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTQCDMC_DP);
 
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSHiggsQCDMC_DP);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSHiggsQCDMC_DP);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSHiggsQCDMC_DP);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSHiggsQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSHiggsQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSHiggsQCDMC_DP);
 
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSJetFakesQCDMC_DP);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSJetFakesQCDMC_DP);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSJetFakesQCDMC_DP);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSJetFakesQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSJetFakesQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSJetFakesQCDMC_DP);
 

  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSZTTQCDMC_DP);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSZTTQCDMC_DP);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSZTTQCDMC_DP);
  // Extradist2d.push_back(&jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSZTTQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSZTTQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSZTTQCDMC_DP);


  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSHiggs);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSHiggs);
 
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSJetFakes);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSJetFakes);
   
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSZTT);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSZTT);

  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSHiggsQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSHiggsQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSHiggsQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSHiggsQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSHiggsQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSHiggsQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSHiggsQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSHiggsQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSHiggsQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSHiggsQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSHiggsQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSHiggsQCDMC);

  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSJetFakesQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSJetFakesQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSJetFakesQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSJetFakesQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSJetFakesQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSJetFakesQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSJetFakesQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSJetFakesQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSJetFakesQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSJetFakesQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSJetFakesQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSJetFakesQCDMC);

  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSZTTQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSZTTQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSZTTQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSZTTQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSZTTQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSZTTQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSZTTQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSZTTQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSZTTQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSZTTQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSZTTQCDMC);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSZTTQCDMC);

  
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSHiggs_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSHiggs_DP);
 
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSJetFakes_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSJetFakes_DP);
   
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSZTT_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSZTT_DP);

  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSHiggsQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSHiggsQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSHiggsQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSHiggsQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSHiggsQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSHiggsQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSHiggsQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSHiggsQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSHiggsQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSHiggsQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSHiggsQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSHiggsQCDMC_DP);

  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSJetFakesQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSJetFakesQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSJetFakesQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSJetFakesQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSJetFakesQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSJetFakesQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSJetFakesQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSJetFakesQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSJetFakesQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSJetFakesQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSJetFakesQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSJetFakesQCDMC_DP);

  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSZTTQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSZTTQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSZTTQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSZTTQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSZTTQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSZTTQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSZTTQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSZTTQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSZTTQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSZTTQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSZTTQCDMC_DP);
  Extradist2d.push_back(&jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSZTTQCDMC_DP);

}

void  HCPTauTau::doEvent()  { //  Method called on every event

  unsigned int t;                // sample type, you may manage in your further analysis, if needed
  int id(Ntp->GetMCID());  //read event ID of a sample

  //Because some problems with the prod
  // if (InputNtuplePath.Contains("EWKWMinus2Jets"))id=202;
  // else if (InputNtuplePath.Contains("EWKWPlus2Jets"))id=201;
  // else if (InputNtuplePath.Contains("EWKZ2Jets_Z_ll"))id=203;
  // else if (InputNtuplePath.Contains("EWKZ2Jets_Z_nunu"))id=27;
  //if (InputNtuplePath.Contains("WplusHToTauTau"))id=460;
  //else if (InputNtuplePath.Contains("WminusHToTauTau"))id=461;
  // else if (InputNtuplePath.Contains("TTTo2L2Nu"))id=701;
  // else if (InputNtuplePath.Contains("TTToHadronic"))id=702;
  // else if (InputNtuplePath.Contains("TTToSemiLeptonic"))id=703;
  // else if (InputNtuplePath.Contains("Embed"))id=36;
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}  //  gives a warning if list of samples in Histo.txt  and SkimSummary.log do not coincide 

  bool TESUp=false;
  bool TESDown=false;
  bool JERUp=false;
  bool JERDown=false;
  bool jetAbsoluteUp=false;
  bool jetBBEC1Up=false;
  bool jetEC2Up=false;
  bool jetFlavorQCDUp=false;
  bool jetHFUp=false;
  bool jetRelativeBalUp=false;
  bool jetAbsoluteYearUp=false;
  bool jetBBEC1YearUp=false;
  bool jetEC2YearUp=false;
  bool jetHFYearUp=false;
  bool jetRelativeSampleYearUp=false;
  bool jetAbsoluteDown=false;
  bool jetBBEC1Down=false;
  bool jetEC2Down=false;
  bool jetFlavorQCDDown=false;
  bool jetHFDown=false;
  bool jetRelativeBalDown=false;
  bool jetAbsoluteYearDown=false;
  bool jetBBEC1YearDown=false;
  bool jetEC2YearDown=false;
  bool jetHFYearDown=false;
  bool jetRelativeSampleYearDown=false;
  bool METScaleUp=false;
  bool METResoUp=false;
  bool METUnclusteredScaleUp=false;
  bool METScaleDown=false;
  bool METResoDown=false;
  bool METUnclusteredScaleDown=false;
  // bool ScaleggUp=false;
  // bool ISRUp=false;
  // bool FSRUp=false;
  // bool ScaleggDown=false;
  // bool ISRDown=false;
  // bool FSRDown=false;


  string TES="Nom";
  if(TESUp) TES="Up";
  else if(TESDown) TES="Down";

  string JER="Nom";
  if(JERUp) JER="Up";
  else if(JERDown) JER="Down";

  string METScale="Nom";
  if(METScaleUp) METScale="Up";
  else if(METScaleDown) METScale="Down";

  string METReso="Nom";
  if(METResoUp) METReso="Up";
  else if(METResoDown) METReso="Down";

  string METUnclusteredScale="Nom";
  if(METUnclusteredScaleUp) METReso="Up";
  else if(METUnclusteredScaleDown) METReso="Down";
  
  bool JESUp=false,JESDown=false;
  if(jetAbsoluteUp||jetBBEC1Up ||jetEC2Up ||jetFlavorQCDUp ||jetHFUp ||jetRelativeBalUp || jetAbsoluteYearUp||jetBBEC1YearUp ||jetEC2YearUp ||jetHFYearUp ||jetRelativeSampleYearUp)JESUp=true;
  else if(jetAbsoluteDown||jetBBEC1Down ||jetEC2Down ||jetFlavorQCDDown ||jetHFDown ||jetRelativeBalDown || jetAbsoluteYearDown||jetBBEC1YearDown ||jetEC2YearDown ||jetHFYearDown ||jetRelativeSampleYearDown)JESDown=true;

  bool isEmbed=(id==36);
  //cout<<isEmbed<<endl;
  // value.at(ZTTMC)=1;
  // if(id==33)value.at(ZTTMC)=false;
  // else value.at(ZTTMC)=true;
  // pass.at(ZTTMC)=value.at(ZTTMC);
  //if(!pass.at(ZTTMC))cout<<"Passe pas !!!!!!!!!!!!!!!!"<<endl;

  bool trig=false;
  std::vector<int> TauIndex ;
  std::vector<int> TriggerIndexVector ;
  std::vector<TString>  MatchedTriggerNames;
  value.at(Trigger)=0;

  std::vector<int> Tau1IndexVect;
  std::vector<int> Tau2IndexVect;
  std::vector<int> IndexSelected;
  for(int iDaughter=0;   iDaughter  <  Ntp->SelectedPairs() ;iDaughter++ ) {
    Tau1IndexVect.push_back(Ntp->tau1IndexVect(iDaughter));
    Tau2IndexVect.push_back(Ntp->tau2IndexVect(iDaughter));
  }      

  if(Ntp->year()==2018 && isEmbed)
    {
      MatchedTriggerNames.push_back("HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v");
      TriggerIndexVector=Ntp->GetVectorTriggers(MatchedTriggerNames);
      for(unsigned int iDaughter=0;   iDaughter  <  Ntp->SelectedPairs() ;iDaughter++ ) {
	//cout<<Ntp->CHECK_BIT(Ntp->Daughters_trgMatched(Tau1IndexVect.at(iDaughter)),31)<<" "<<Ntp->CHECK_BIT(Ntp->Daughters_trgMatched(Tau2IndexVect.at(iDaughter)),31)<<endl;
	if(Ntp->CHECK_BIT(Ntp->Daughters_trgMatched(Tau1IndexVect.at(iDaughter)),31) && Ntp->CHECK_BIT(Ntp->Daughters_trgMatched(Tau2IndexVect.at(iDaughter)),31))
	  {
	    IndexSelected.push_back(iDaughter);
	  }
      }
    }
  else if(Ntp->year()==2018 && (!Ntp->isData() || (Ntp->isData() && Ntp->RunNumber() >= 317509)))
    {
      MatchedTriggerNames.push_back("HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v");
      TriggerIndexVector=Ntp->GetVectorTriggers(MatchedTriggerNames);
      for(unsigned int iDaughter=0;   iDaughter  <  Ntp->SelectedPairs() ;iDaughter++ ) {
	//cout<<Ntp->CHECK_BIT(Ntp->Daughters_trgMatched(Tau1IndexVect.at(iDaughter)),31)<<" "<<Ntp->CHECK_BIT(Ntp->Daughters_trgMatched(Tau2IndexVect.at(iDaughter)),31)<<endl;
	if(Ntp->CHECK_BIT(Ntp->Daughters_trgMatched(Tau1IndexVect.at(iDaughter)),31) && Ntp->CHECK_BIT(Ntp->Daughters_trgMatched(Tau2IndexVect.at(iDaughter)),31))
	  {
	    IndexSelected.push_back(iDaughter);
	  }
      }
    }
  else if(Ntp->year()==2018 && (Ntp->isData() && Ntp->RunNumber() < 317509))
    {
      MatchedTriggerNames.push_back("HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v");
      MatchedTriggerNames.push_back("HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v");
      MatchedTriggerNames.push_back("HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v");
      TriggerIndexVector=Ntp->GetVectorTriggers(MatchedTriggerNames);
      for(unsigned int iDaughter=0;   iDaughter  <  Ntp->SelectedPairs() ;iDaughter++ ) {
	if((Ntp->CHECK_BIT(Ntp->Daughters_trgMatched(Tau1IndexVect.at(iDaughter)),42) || Ntp->CHECK_BIT(Ntp->Daughters_trgMatched(Tau1IndexVect.at(iDaughter)),44)|| Ntp->CHECK_BIT(Ntp->Daughters_trgMatched(Tau1IndexVect.at(iDaughter)),45)) && (Ntp->CHECK_BIT(Ntp->Daughters_trgMatched(Tau2IndexVect.at(iDaughter)),42) || Ntp->CHECK_BIT(Ntp->Daughters_trgMatched(Tau2IndexVect.at(iDaughter)),44)|| Ntp->CHECK_BIT(Ntp->Daughters_trgMatched(Tau2IndexVect.at(iDaughter)),45)))
	  {
	    IndexSelected.push_back(iDaughter);
	  }
      }
    }
  else if(Ntp->year()==2017){
    MatchedTriggerNames.push_back("HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v");
    MatchedTriggerNames.push_back("HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v");
    MatchedTriggerNames.push_back("HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v");
    TriggerIndexVector=Ntp->GetVectorTriggers(MatchedTriggerNames);
    for(unsigned int iDaughter=0;   iDaughter  <  Ntp->SelectedPairs() ;iDaughter++ ) {
      IndexSelected.push_back(iDaughter);
    }
  }
  else if(Ntp->year()==2016){
    MatchedTriggerNames.push_back("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v");
    MatchedTriggerNames.push_back("HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v");
    TriggerIndexVector=Ntp->GetVectorTriggers(MatchedTriggerNames);
    for(unsigned int iDaughter=0;   iDaughter  <  Ntp->SelectedPairs() ;iDaughter++ ) {
      IndexSelected.push_back(iDaughter);
    }
  }
  else
    {
      for(unsigned int iDaughter=0;   iDaughter  <  Ntp->SelectedPairs() ;iDaughter++ ) {
	IndexSelected.push_back(iDaughter);
      }
    }
  
  std::vector<int> out;
  for(unsigned int i=0; i<Ntp->NTriggers();i++){
    TString name=Ntp->TriggerName(i);
    //if(Ntp->TriggerAccept(i))cout<<Ntp->TriggerName(i)<<endl;
    for(unsigned int j=0; j<MatchedTriggerNames.size(); j++){
      
      if(name.Contains(MatchedTriggerNames.at(j))) out.push_back(i) ;
    }
  }

  TriggerIndexVector=out;
  for(unsigned int itrig = 0; itrig < TriggerIndexVector.size(); itrig++){
    //cout<<"indice: "<<TriggerIndexVector.at(itrig)<<endl;
    //cout<<"Accept: "<<Ntp->TriggerAccept(TriggerIndexVector.at(itrig))<<endl;
    if(Ntp->TriggerAccept(TriggerIndexVector.at(itrig))){
      trig=true;
    }
  }
  
  if(isEmbed)trig=true;

  if(trig && IndexSelected.size()>0)value.at(Trigger)=1;
  else value.at(Trigger)=0;
  //if(Ntp->year()!=2018)value.at(Trigger)=1;
  pass.at(Trigger)=(value.at(Trigger)==cut.at(Trigger));

  bool METFilters=true;

  if(Ntp->year()==2017 || Ntp->year()==2018){
    if(!Ntp->passecalBadCalibFilterUpdate()){
      METFilters=false;
    }
  }
  if(Ntp->metfilterbit()!=127)METFilters=false;
  /*
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
  */

  std::vector<int> goodTaus1Index;
  std::vector<int> goodTaus2Index;
  std::vector<int> IndexSelectedTemp;
  std::vector<int>  PairsIndexTemp;
  int j=0;
  int GenMatch1=6, GenMatch2=6;

  for(unsigned int iDaughter=0;   iDaughter  <IndexSelected.size() ;iDaughter++ ) {
    if(!Ntp->isData()|| isEmbed){GenMatch1=Ntp->gen_match_1(IndexSelected.at(iDaughter));GenMatch2=Ntp->gen_match_2(IndexSelected.at(iDaughter));}
    else{GenMatch1=6;GenMatch2=6;}

    if(Ntp->tauBaselineSelection(Tau1IndexVect.at(iDaughter),GenMatch1, 40., 2.1, 1,0,4,TES) && Ntp->tauBaselineSelection(Tau2IndexVect.at(iDaughter),GenMatch2,40., 2.1, 1,0,4,TES) && METFilters)
      {
  	goodTaus1Index.push_back(Tau1IndexVect.at(iDaughter));
  	goodTaus2Index.push_back(Tau2IndexVect.at(iDaughter));
  	IndexSelectedTemp.push_back(IndexSelected.at(iDaughter));
  	PairsIndexTemp.push_back(j);
  	j++;
      }
  }

  value.at(Id_and_Kin)=0;
  value.at(Id_and_Kin)=(IndexSelectedTemp.size()>0);
  pass.at(Id_and_Kin)=value.at(Id_and_Kin);

  int Tau1= -1;
  int Tau2= -1;
  std::vector<int>  Sorted;
  bool dileptonveto;
  bool GenMatchSelection=false;
  if(PairsIndexTemp.size()>0)
    {

      Sorted = Ntp->SortPair(PairsIndexTemp,goodTaus1Index,goodTaus2Index);
      Tau1=goodTaus1Index.at(Sorted.back());
      Tau2=goodTaus2Index.at(Sorted.back());

      //Tau1=goodTaus1Index.at(0);
      //Tau2=goodTaus2Index.at(0);
      if(!Ntp->isData()|| isEmbed){GenMatch1=Ntp->gen_match_1(IndexSelectedTemp.at(Sorted.back()));GenMatch2=Ntp->gen_match_2(IndexSelectedTemp.at(Sorted.back()));}
      else{GenMatch1=6;GenMatch2=6;}

      //      cout<<IndexSelectedTemp.at(Sorted.back())<<"  "<<IndexSelectedTemp.at(0)<<endl;
      
      //value.at(genmatch)=false;
      //cout<<id<<"  "<<GenMatch1<<"  "<<GenMatch2<<endl;
      //cout<<id<<" "<<(!(GenMatch1==6 || GenMatch2==6) && !(GenMatch1==5 && GenMatch2==5))<<endl;
      if(id==30 || id==33 ||id==203)GenMatchSelection=(!(GenMatch1==6 || GenMatch2==6) && !(GenMatch1==5 && GenMatch2==5)); //ZL
      else if(id==36)GenMatchSelection=(GenMatch1==5 && GenMatch2==5);
      else if(id==23 && id==20)GenMatchSelection=false; //WTaunu
      else if(id==201|| id==202)GenMatchSelection=!((GenMatch1==5 && GenMatch2==6)||(GenMatch2==5 && GenMatch1==6));
      //else if((Ntp->year()==2016 && id==70) || ((Ntp->year()==2017 || Ntp->year()==2018 ) && (id==701 ||id==702 ||id==703)))GenMatchSelection = !((GenMatch1==6 || GenMatch2==6) || ((GenMatch1 == 5) && (GenMatch2 == 5))); //ttbar veto due to embedded ttbar contamination
      else if(id==71||id==72||id==73||id==74||id==47||id==48||id==49||id==50||id==51||id==52||id==53||id==54||id==55||id==56||id==57||id==58)GenMatchSelection=(!(GenMatch1==6 || GenMatch2==6) &&  !(GenMatch1==5 && GenMatch2==5)); // VVT
      else if(id==70 ||id==701 ||id==702 ||id==703 )GenMatchSelection=(!(GenMatch1==6 || GenMatch2==6) &&  !(GenMatch1==5 && GenMatch2==5));
      else GenMatchSelection=true;

      //value.at(genmatch)=GenMatchSelection;
      //pass.at(genmatch)=GenMatchSelection;
      //value.at(GoodIndex)=0;
      //value.at(GoodIndex)=(Tau1!=-99 && Tau2!=-99);
      //pass.at(GoodIndex) = value.at(GoodIndex);
      value.at(TausIsolation)=0;
      value.at(TausIsolation) = (Ntp->byMediumDeepTau2017v2p1VSjet_1(IndexSelectedTemp.at(Sorted.back())) && Ntp->byMediumDeepTau2017v2p1VSjet_2(IndexSelectedTemp.at(Sorted.back())));
      pass.at(TausIsolation) = value.at(TausIsolation);
      value.at(AgainstEleMu)=0;
      value.at(AgainstEleMu) = (Ntp->byVVLooseDeepTau2017v2p1VSe_1(IndexSelectedTemp.at(Sorted.back())) && Ntp->byVVLooseDeepTau2017v2p1VSe_2(IndexSelectedTemp.at(Sorted.back())) && Ntp->byVLooseDeepTau2017v2p1VSmu_1(IndexSelectedTemp.at(Sorted.back())) && Ntp->byVLooseDeepTau2017v2p1VSmu_2(IndexSelectedTemp.at(Sorted.back())) );
      pass.at(AgainstEleMu) = value.at(AgainstEleMu);
      value.at(LeptonVeto)=0;
      int thirdlepton=0;
      for(int iDaughter=0;   iDaughter  <  Ntp->NDaughters() ;iDaughter++ ) {
	if((iDaughter!=Tau1)&&(iDaughter!=Tau2)){
	  if(Ntp->ElectronVeto(iDaughter) || Ntp->MuonVeto(iDaughter))thirdlepton++;
	}
      }

      dileptonveto=(Ntp->DiEleVeto() || Ntp->DiMuonVeto());
      //if(id==30 && value.at(AgainstEleMu)==1 && value.at(TausIsolation)==1)cout<<Ntp->DiEleVeto() <<endl;
      //if(Ntp->eleveto() || Ntp->muonveto())thirdlepton++;
      value.at(LeptonVeto) = (thirdlepton>0 || dileptonveto==true);
      pass.at(LeptonVeto) = (value.at(LeptonVeto)==cut.at(LeptonVeto));
      value.at(PairCharge)=false;
      bool isOS=false;
      isOS=((Ntp->Daughters_charge(Tau1)/abs(Ntp->Daughters_charge(Tau1))) != (Ntp->Daughters_charge(Tau2)/abs(Ntp->Daughters_charge(Tau2))));
      if(isOS)value.at(PairCharge) = true;
      pass.at(PairCharge) = value.at(PairCharge);
      value.at(PairMass) = 999.;
      // // //value.at(MTM) = 999.;
      // // //value.at(MTM) = .;
      value.at(PairMass)=(Ntp->P4Corrected(Tau1,GenMatch1,TES) + Ntp->P4Corrected(Tau2,GenMatch2,TES)).M();
      pass.at(PairMass) = (value.at(PairMass) > cut.at(PairMass));
      // // //pass.at(MTM) = (value.at(MTM) <= cut.at(MTM));
    }

  // Here you can defined different type of weights you want to apply to events.
  double wobs=1;
  double w=1;
  double wIDSF1=1.,wIDSF2=1.;
  double wIDSFAntiE1=1.,wIDSFAntiE2=1.;
  double wIDSFAntiMu1=1.,wIDSFAntiMu2=1.;
  double wTrgSF1=1.,wTrgSF2=1.;

  if(!Ntp->isData() && id!=DataMCType::QCD && !isEmbed) w*=Ntp->PUReweight();

  if((!Ntp->isData() && id!=DataMCType::QCD) || isEmbed) {
    if(PairsIndexTemp.size()>0){
      //double w1=1.,w2=1.;
      
      wTrgSF1=Ntp->TriggerSF(Tau1,GenMatch1,TES,"Nom");
      wTrgSF2=Ntp->TriggerSF(Tau2,GenMatch2,TES,"Nom");

      w*=wTrgSF1;
      w*=wTrgSF2;
      
      wIDSF1=Ntp->IDSF(Tau1,GenMatch1,TES);
      wIDSF2=Ntp->IDSF(Tau2,GenMatch2,TES);
      if(!isEmbed){
	wIDSFAntiE1=Ntp->IDSF(Tau1,GenMatch1,TES,"ele");
	wIDSFAntiE2=Ntp->IDSF(Tau2,GenMatch2,TES,"ele");
	wIDSFAntiMu1=Ntp->IDSF(Tau1,GenMatch1,TES,"mu");
	wIDSFAntiMu2=Ntp->IDSF(Tau2,GenMatch2,TES,"mu");
      }
      w*=wIDSF1*wIDSF2*wIDSFAntiE1*wIDSFAntiE2*wIDSFAntiMu1*wIDSFAntiMu2;
      
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
    zptw = Ntp->ZPtReweight(genMomentum);
    //zptw = DataMC_Corr.ZPTWeight(genMomentum.M(),genMomentum.Pt());
  }  
  w*=zptw;
  
  double top_wt = 1.0;
  if(id==70 ||id==701||id==702||id==703){
    for (unsigned i = 0; i < Ntp->NGenParts(); i++){
      unsigned pdgid = abs(Ntp->Genpart_pdg(i));
      if(pdgid == 6 && Ntp->Genpart_flags(i)==8 && Ntp->Genpart_flags(i)==13){ // Flag From Hard Process and isLastCopy
	double pt = Ntp->Genpart_P4(i).Pt();
	pt = std::min(pt, 472.);
	double a = 0.088, b = -0.00087, c = 9.2e-07;
	top_wt *= std::exp(a + b * pt + c * pt*pt); 
      }
    }
    top_wt = std::sqrt(top_wt);
    w*=top_wt;
  }
  int imc1=-1,imc2=-1;
  if(isEmbed){
    for(unsigned int imc=0; imc < Ntp->NGenParts(); imc++){
      if(fabs(Ntp->Genpart_pdg(imc))==15 && imc1==-1)imc1=imc;
      else if(fabs(Ntp->Genpart_pdg(imc))==15 && imc1!=-1)imc2=imc;
    }
    w*=Ntp->EmbeddingSelectionSF(imc1,imc2);
  }

  float PUPPImetCorr_px=Ntp->PUPPImet()*cos(Ntp->PUPPImetphi());
  float PUPPImetCorr_py=Ntp->PUPPImet()*sin(Ntp->PUPPImetphi());

  double jdeta_=-9999, jpt_1_=-9999, m_sv_=-9999, m_vis_=-9999, met_=-9999, mjj_=-9999, n_jets_=-9999, pt_1_=-9999,  pt_tt_=-9999, pt_vis_=-9999;
  int njets20=-99;

  TLorentzVector NullTLV(0,0,0,0);
  int maxjets=0;
  int maxjetsUp=0;
  int maxjetsDown=0;
  vector<TLorentzVector> JetsVectNominal;
  
  int indice=0;
  int njetsJES=0;
  // int njetsJESUp;
  // int njetsJERDown;
  int njetsJERUp=0;
  int njetsJERDown=0;
  int m=0;

  if(PairsIndexTemp.size()>0){
    if(JERUp) {n_jets_=Ntp->njetsUp(IndexSelectedTemp.at(Sorted.back()));njets20=Ntp->njetspt20Up(IndexSelectedTemp.at(Sorted.back()));}
    else if(JERDown) {n_jets_=Ntp->njetsDown(IndexSelectedTemp.at(Sorted.back()));njets20=Ntp->njetspt20Down(IndexSelectedTemp.at(Sorted.back()));}
    else {n_jets_=Ntp->njets(IndexSelectedTemp.at(Sorted.back()));njets20=Ntp->njetspt20(IndexSelectedTemp.at(Sorted.back()));}
  }


  // for(int i=0;i<Ntp->njetspt20Size();i++)if(maxjets<Ntp->njetspt20(i))maxjets=Ntp->njetspt20(i);
  // vector<vector<TLorentzVector> >JetsVect(Ntp->njetspt20Size(),vector<TLorentzVector>(maxjets,NullTLV));

  // for(int i=0;i<Ntp->njetspt20SizeUp();i++)if(maxjetsUp<Ntp->njetspt20Up(i))maxjetsUp=Ntp->njetspt20Up(i);
  // vector<vector<TLorentzVector> >JetsVectJERUp(Ntp->njetspt20SizeUp(),vector<TLorentzVector>(maxjetsUp,NullTLV));

  // for(int i=0;i<Ntp->njetspt20SizeDown();i++)if(maxjetsDown<Ntp->njetspt20Down(i))maxjetsDown=Ntp->njetspt20Down(i);
  // vector<vector<TLorentzVector> >JetsVectJERDown(Ntp->njetspt20SizeDown(),vector<TLorentzVector>(maxjetsDown,NullTLV));

  vector<TLorentzVector> JetsVect;
  vector<TLorentzVector> JetsVectJERUp;
  vector<TLorentzVector> JetsVectJERDown;

  // if(Ntp->njetspt20Size()>0)njetsJES=Ntp->njetspt20(indice);
  // if(Ntp->njetspt20SizeUp()>0)njetsJERUp=Ntp->njetspt20Up(indice);
  // if(Ntp->njetspt20SizeDown()>0)njetsJERDown=Ntp->njetspt20Down(indice);

  if(Ntp->njetsSize()>0)njetsJES=Ntp->njetspt20(indice);
  if(Ntp->njetsSizeUp()>0)njetsJERUp=Ntp->njetsUp(indice);
  if(Ntp->njetsSizeDown()>0)njetsJERDown=Ntp->njetsDown(indice);


  //   cout<<"PairsIndexTemp.size: "<<PairsIndexTemp.size()<<endl;
  // int njetpercand[Ntp->njetspt20Size()];
  // for(int i=0;i<Ntp->njetspt20Size();i++)
  //    {
  //      njetpercand[i]=Ntp->njetspt20(i);
  //    }
  int njetpercand[Ntp->njetspt20Size()];
  for(int i=0;i<Ntp->njetspt20Size();i++)
    {
      njetpercand[i]=Ntp->njetspt20(i);
    }
  

  if(PairsIndexTemp.size()>0){
    // cout<<"NJets: "<<Ntp->NJets()<<endl;
    // cout<<"Liste: ";
    // for(int i=0;i<Ntp->njetspt20Size();i++)
    //  {
    //    cout<<Ntp->njetspt20(i)<<"  ";
    //  }
    // cout<<endl;
    // cout<<"Candidat: "<<IndexSelectedTemp.at(0)<<endl;
    // // cout<<1.4<<endl;
    // for(int i=0;i<Ntp->NJets();i++)
    //   {
    // 	cout<<"i: "<<i<<endl;
    // 	cout<<"njetsJES: "<<njetsJES<<endl;
    // 	cout<<"indice0: "<<indice<<endl;
    // 	if(Ntp->njetspt20(indice)==0){indice++;}
    // 	if(i>=njetsJES){indice++;njetsJES+=Ntp->njetspt20(indice);}
    // 	//if(Ntp->njetspt20(indice)==0)indice++;

    // 	if(indice==IndexSelectedTemp.at(Sorted.back())){JetsVect.push_back(Ntp->Jet_P4(i));}
    // 	cout<<"njetsJES1: "<<njetsJES<<endl;
    // 	cout<<"indice1: "<<indice<<endl;
    //   }
    int j=0;
    m=0;
    // for(int i=0;i<Ntp->njetspt20Size();i++)
    //   {
    // 	m+=Ntp->njetspt20(i);
    // 	while(j<m){
    // 	  if(i==IndexSelectedTemp.at(Sorted.back()))JetsVect.push_back(Ntp->Jet_P4(j));
    // 	  j++;
    // 	}
    //   }
    for(int i=0;i<Ntp->njetspt20Size();i++)
      {
	m+=Ntp->njetspt20(i);
	while(j<m){
	  if(i==IndexSelectedTemp.at(Sorted.back()))JetsVect.push_back(Ntp->Jet_P4(j));
	  j++;
	}
      }

    indice=0;
    j=0;
    m=0;
 
    // for(int i=0;i<Ntp->njetspt20SizeUp();i++)
    //   {
    // 	m+=Ntp->njetspt20Up(i);
    // 	while(j<m){
    // 	  if(i==IndexSelectedTemp.at(Sorted.back()))JetsVectJERUp.push_back(Ntp->JetUp_P4(j));
    // 	  j++;
    // 	}
    //   }

    for(int i=0;i<Ntp->njetsSizeUp();i++)
      {
	m+=Ntp->njetsUp(i);
	while(j<m){
	  if(i==IndexSelectedTemp.at(Sorted.back()))JetsVectJERUp.push_back(Ntp->JetUp_P4(j));
	  j++;
	}
      }
    j=0;
    indice=0;
    m=0;

    // for(int i=0;i<Ntp->njetspt20SizeDown();i++)
    //   {
    // 	m+=Ntp->njetspt20Down(i);
    // 	while(j<m){
    // 	  if(i==IndexSelectedTemp.at(Sorted.back()))JetsVectJERDown.push_back(Ntp->JetDown_P4(j));
    // 	  j++;
    // 	}
    //   }

    for(int i=0;i<Ntp->njetsSizeDown();i++)
      {
	m+=Ntp->njetsDown(i);
	while(j<m){
	  if(i==IndexSelectedTemp.at(Sorted.back()))JetsVectJERDown.push_back(Ntp->JetDown_P4(j));
	  j++;
	}
      }
  }

  indice=0;
  m=0;

  int decalage=0;
  double SumPxNom=0., SumPyNom=0.,SumPxShift=0., SumPyShift=0.;
  int njets20Nom=0,njetsNom=0;
  vector<TLorentzVector> JetsVectJES;
  if(PairsIndexTemp.size()>0){

    if(njets20>0)
      {

	for(int i=0;i<Ntp->njetspt20Size();i++)
	  {

	    if(i<IndexSelectedTemp.at(Sorted.back()))decalage+=Ntp->njetspt20(i);
	  }

	//njets20Nom=Ntp->njetspt20(IndexSelectedTemp.at(Sorted.back()));
	njetsNom=Ntp->njetspt20(IndexSelectedTemp.at(Sorted.back()));
	//float JESshift[njets20]={0.};
	float JESshift[njets20]={0.};
	//for(int i=0;i<njets20;i++)

	for(int i=0;i<njets20;i++)
	  {

	    JESshift[i]=0;
	    if(jetAbsoluteUp)JESshift[i]=Ntp->jets_jetUnc_Absolute_up(decalage+i);
	    else if(jetBBEC1Up)JESshift[i]=Ntp->jets_jetUnc_BBEC1_up(decalage+i);
	    else if(jetEC2Up)JESshift[i]=Ntp->jets_jetUnc_EC2_up(decalage+i);
	    else if(jetFlavorQCDUp)JESshift[i]=Ntp->jets_jetUnc_FlavorQCD_up(decalage+i);
	    else if(jetHFUp)JESshift[i]=Ntp->jets_jetUnc_HF_up(decalage+i);
	    else if(jetRelativeBalUp)JESshift[i]=Ntp->jets_jetUnc_RelativeBal_up(decalage+i);
	    else if(jetAbsoluteYearUp)JESshift[i]=Ntp->jets_jetUnc_Absolute_YEAR_up(decalage+i);
	    else if(jetBBEC1YearUp)JESshift[i]=Ntp->jets_jetUnc_BBEC1_YEAR_up(decalage+i);
	    else if(jetEC2YearUp)JESshift[i]=Ntp->jets_jetUnc_EC2_YEAR_up(decalage+i);
	    else if(jetHFYearUp)JESshift[i]=Ntp->jets_jetUnc_HF_YEAR_up(decalage+i);
	    else if(jetRelativeSampleYearUp)JESshift[i]=Ntp->jets_jetUnc_RelativeSample_YEAR_up(decalage+i);
	    else if(jetAbsoluteDown)JESshift[i]=-Ntp->jets_jetUnc_Absolute_dw(decalage+i);
	    else if(jetBBEC1Down)JESshift[i]=-Ntp->jets_jetUnc_BBEC1_dw(decalage+i);
	    else if(jetEC2Down)JESshift[i]=-Ntp->jets_jetUnc_EC2_dw(decalage+i);
	    else if(jetFlavorQCDDown)JESshift[i]=-Ntp->jets_jetUnc_FlavorQCD_dw(decalage+i);
	    else if(jetHFDown)JESshift[i]=-Ntp->jets_jetUnc_HF_dw(decalage+i);
	    else if(jetRelativeBalDown)JESshift[i]=-Ntp->jets_jetUnc_RelativeBal_dw(decalage+i);
	    else if(jetAbsoluteYearDown)JESshift[i]=-Ntp->jets_jetUnc_Absolute_YEAR_dw(decalage+i);
	    else if(jetBBEC1YearDown)JESshift[i]=-Ntp->jets_jetUnc_BBEC1_YEAR_dw(decalage+i);
	    else if(jetEC2YearDown)JESshift[i]=-Ntp->jets_jetUnc_EC2_YEAR_dw(decalage+i);
	    else if(jetHFYearDown)JESshift[i]=-Ntp->jets_jetUnc_HF_YEAR_dw(decalage+i);
	    else if(jetRelativeSampleYearDown)JESshift[i]=-Ntp->jets_jetUnc_RelativeSample_YEAR_dw(decalage+i);

	  }

	JetsVectNominal=JetsVect;

	//for(int i=0;i<njets20;i++)
	for(int i=0;i<JetsVect.size();i++)
	  {
	    //    cout<<"i: "<<i<<endl;
	    //cout<<"JetsVect: ";JetsVect[i].Print();
	    if((JetsVect[i]+JESshift[i]*JetsVect[i]).Pt()>30.)JetsVectJES.push_back(JetsVect[i]+JESshift[i]*JetsVect[i]);
	  }

	if(JetsVectJES.size()>1)std::sort(JetsVectJES.begin(),JetsVectJES.end(),Ntp->ComparePairsbyPt);
	if(JetsVectJERUp.size()>1)std::sort(JetsVectJERUp.begin(),JetsVectJERUp.end(),Ntp->ComparePairsbyPt);
	if(JetsVectJERDown.size()>1)std::sort(JetsVectJERDown.begin(),JetsVectJERDown.end(),Ntp->ComparePairsbyPt);

	for(int i=0;i<JetsVectJES.size()/*njets20*/;i++)
	  {
	    if(JESUp || JESDown){
	      SumPxShift+=JetsVectJES[i].Px();
	      SumPyShift+=JetsVectJES[i].Py();
	    }
	  }

	for(int i=0;i<JetsVectNominal.size()/*njets20Nom*/;i++)
	  {
	    SumPxNom+=JetsVectNominal[i].Px();
	    SumPyNom+=JetsVectNominal[i].Py();
	  }

	for(int i=0;i<JetsVectJERUp.size()/*njets20*/;i++)
	  {
	    if(JERUp){
	      SumPxShift+=JetsVectJERUp[i].Px();
	      SumPyShift+=JetsVectJERUp[i].Py();
	    }
	  }

	for(int i=0;i<JetsVectJERDown.size()/*njets20*/;i++){
	  if(JERDown){
	    SumPxShift+=JetsVectJERDown[i].Px();
	    SumPyShift+=JetsVectJERDown[i].Py();
	  }
	}

      }

    if(TESUp || TESDown){
      PUPPImetCorr_px=PUPPImetCorr_px + Ntp->Daughters_P4(Tau1).Px()+Ntp->Daughters_P4(Tau2).Px()-(Ntp->P4Corrected(Tau1,GenMatch1,TES).Px()+Ntp->P4Corrected(Tau2,GenMatch2,TES).Px());
      PUPPImetCorr_py=PUPPImetCorr_py + Ntp->Daughters_P4(Tau1).Py()+Ntp->Daughters_P4(Tau2).Py()-(Ntp->P4Corrected(Tau1,GenMatch1,TES).Py()+Ntp->P4Corrected(Tau2,GenMatch2,TES).Py());
    }

  }

  

  if(id==33 || id == 10110333 || id == 10110433|| id == 10130533|| id ==10210333|| id == 10210433|| id == 10230533|| id ==10310333 || id ==10330533 || id ==10410433 || id == 10410333|| id == 10430533|| id == 30530533 || id==30 || id==11 || id==12 || id==20 || id==23 || id==21 || id==22 ||id==23 ||id==45  ||id==460||id == 461)
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
	  Ntp->RecoilCorr(Gen, Vis, IndexSelectedTemp.at(Sorted.back()),PUPPImetCorr_px,PUPPImetCorr_py,JER,METScale,METReso);
	}
    }
  else if(PairsIndexTemp.size()>0 && !Ntp->isData() && !isEmbed){
    if(METUnclusteredScale=="Up"){
      PUPPImetCorr_px=Ntp->puppimet_ex_UnclusteredEnUp();
      PUPPImetCorr_py=Ntp->puppimet_ey_UnclusteredEnUp();
    }
    else if(METUnclusteredScale=="Down"){
      PUPPImetCorr_px=Ntp->puppimet_ex_UnclusteredEnDown();
      PUPPImetCorr_py=Ntp->puppimet_ey_UnclusteredEnDown();
    }

    if(JERUp || JERDown){
      PUPPImetCorr_px=PUPPImetCorr_px + SumPxNom-SumPxShift;
      PUPPImetCorr_py=PUPPImetCorr_py + SumPyNom-SumPyShift;
    }

    if(JESUp || JESDown){
      PUPPImetCorr_px=PUPPImetCorr_px + SumPxNom-SumPxShift;
      PUPPImetCorr_py=PUPPImetCorr_py + SumPyNom-SumPyShift;
    }
  }

  if((!Ntp->isData() && id!=DataMCType::QCD) || isEmbed){
    
    //if(Ntp->MC_weight()>0.4)
    //cout<<"Ntp->MC_weight()"<<Ntp->MC_weight()<<endl;
    if (isEmbed && Ntp->MC_weight()>10000.)w*=Ntp->MC_weight()*0.000000001; //problem with 2016
    //if(Ntp->MC_weight()<0)w*=-Ntp->MC_weight();
    else w*=Ntp->MC_weight();
    
  }  //generator weight because negative weights for this samples
  if(Ntp->year()==2016)
    {
      if(id==11)w*=0.2455;
      else if(id==12)w*=0.2727;
      else if(id==45)w*=0.2546;
      else if(id==461)w*=0.2596;
      else if(id==460)w*=0.2425;
    }
  else if(Ntp->year()==2017)
    {
      if(id==11)w*=0.2447;
      else if(id==12)w*=0.2697;
      else if(id==45)w*=0.2514;
      else if(id==461)w*=0.2567;
      else if(id==460)w*=0.2394;
    }
  else if(Ntp->year()==2018)
    {
      if(id==11)w*=0.2446;
      else if(id==12)w*=0.2695;
      else if(id==45)w*=0.2513;
      else if(id==461)w*=0.2563;
      else if(id==460)w*=0.2397;
    }
  w*=Ntp->stitch_weight(isDY1050);

  //cout<<"w: "<<Ntp->stitch_weight(isDY1050)<<endl;
  //cout<<"isDY1050: "<<isDY1050<<endl;

  if(!Ntp->isData() && !isEmbed && (Ntp->year()==2016||Ntp->year()==2017))w*=Ntp->prefiringweight();
  double prefup;
  double prefdown;
  if((Ntp->prefiringweightup()/Ntp->prefiringweight())>1.2)prefup=Ntp->prefiringweightup();
  else prefup=Ntp->prefiringweight()*1.2;
  if((Ntp->prefiringweightdown()/Ntp->prefiringweight())<0.8)prefdown=Ntp->prefiringweightdown();
  else prefdown=Ntp->prefiringweight()*0.8;
  // if(isEmbed && PairsIndexTemp.size()>0)
  //    {
  //      if(GenMatch1==5)
  // 	 {
  // 	   if(Ntp->decayMode(Tau1)==1)w*=0.975;
  // 	   if(Ntp->decayMode(Tau1)==2)w*=0.975*1.051;
  // 	   if(Ntp->decayMode(Tau1)==5 || Ntp->decayMode(Tau1)==6)w*=pow(0.975,2);
  // 	   if(Ntp->decayMode(Tau1)==10)w*=pow(0.975,3);
  // 	   if(Ntp->decayMode(Tau1)==11)w*=pow(0.975,3)*1.051;
  // 	 }
  //      if(GenMatch2==5){
  // 	 if(Ntp->decayMode(Tau2)==1)w*=0.975;
  // 	 if(Ntp->decayMode(Tau2)==2)w*=0.975*1.051;
  // 	 if(Ntp->decayMode(Tau2)==5 || Ntp->decayMode(Tau2)==6)w*=pow(0.975,2);
  // 	 if(Ntp->decayMode(Tau2)==10)w*=pow(0.975,3);
  // 	 if(Ntp->decayMode(Tau2)==11)w*=pow(0.975,3)*1.051;
  //      }
  //    }
  
  std::vector<unsigned int> exclude_cuts;
  exclude_cuts.push_back(TausIsolation);
  classic_svFit::LorentzVector tau1P4;
  classic_svFit::LorentzVector tau2P4;

  TLorentzVector Tau1P4;
  TLorentzVector Tau2P4;
  if(passAllBut(exclude_cuts)){
    Tau1P4 = Ntp->P4Corrected(Tau1,GenMatch1,TES);
    Tau2P4 = Ntp->P4Corrected(Tau2,GenMatch2,TES);
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
  
  // int Has2jets=-1;
  // int Has1jet=-1;
  // if(PairsIndexTemp.size()>0){
  //   for(int k=0;k<=IndexSelectedTemp.at(Sorted.back());k++)
  //     {
	
  // if (JERUp){
  //   if(Ntp->njetspt20Up(k)>=2)Has2jets++;
  //   if(Ntp->njetspt20Up(k)>=1)Has1jet++;
  // }
  // else if (JERDown){
  //   if(Ntp->njetspt20Down(k)>=2)Has2jets++;
  //   if(Ntp->njetspt20Down(k)>=1)Has1jet++;
  // }
  // else {
  //   if(Ntp->njetspt20(k)>=2)Has2jets++;
  //   if(Ntp->njetspt20(k)>=1)Has1jet++;
  // }
  //     }
  // }

  unsigned long long evt_ = Ntp->EventNumber();
  std::pair<float, int> max_pair;
  std::vector<float> scores = {};
  double PUPPIMET=sqrt(PUPPImetCorr_px*PUPPImetCorr_px+PUPPImetCorr_py*PUPPImetCorr_py);
  double PUPPIMETPhi=(TVector3(PUPPImetCorr_px,PUPPImetCorr_py,0)).Phi();
  std::vector<TLorentzVector> Pions1;
  std::vector<TLorentzVector> Pions2;
  std::vector<double> Pions1Charge;
  std::vector<double> Pions2Charge;

  TVector3 tauPrimaryVertex , tauNoBSPrimaryVertex, tauBSPrimaryVertex, tauNoBSZNominalPrimaryVertex, tauBSZNominalPrimaryVertex, TauminusSecondaryVertex , TauplusSecondaryVertex,tauWithTracksBSZNominalPrimaryVertex,tauWithTracksBSPrimaryVertex;

  TVector3 TauminusDirection , TauplusDirection, TauplusDirectionNoBS, TauplusDirectionBS, TauplusDirectionNoBSZNominal, TauplusDirectionBSZNominal,TauminusDirectionNoBS, TauminusDirectionBS, TauminusDirectionNoBSZNominal, TauminusDirectionBSZNominal;

  double thetaGJ_Tauminus , thetaGJ_Tauplus;
  TLorentzVector a1LV_Tauminus , a1LV_Tauplus, a1LVRefit_Tauminus , a1LVRefit_Tauplus;
  TLorentzVector TauminusPairConstraint, TauplusPairConstraintNoBS, TauplusPairConstraintBS, TauplusPairConstraintNoBSZNominal, TauplusPairConstraintBSZNominal, TauminusPairConstraintNoBS, TauminusPairConstraintBS, TauminusPairConstraintNoBSZNominal, TauminusPairConstraintBSZNominal, TauminusPairConstraintWithTracksBSZNominal,TauminusPairConstraintWithTracksBS, TauplusPairConstraintWithTracksBSZNominal,TauplusPairConstraintWithTracksBS;
  TLorentzVector TauplusSmall, TauplusLarge, TauplusPairConstraint,TauplusPairConstraintMVA,TauplusPairConstraintNoBSNewMVA,TauplusPairConstraintBSNewMVA;
  bool isPlusReal=true, isMinusReal=true,  a1a1=false,a1a1MVA=false, a1minus=false, a1plus=false,a1minusMVA=false,a1plusMVA=false;
  std::vector<TLorentzVector> solutions, solutionsNoBS, solutionsBS, solutionsNoBSZNominal, solutionsBSZNominal, solutionsWithTracksBSZNominal, solutionsWithTracksBS;
  TLorentzVector HPairConstraint;

  std::vector<size_t> hashes;
  size_t hash = 0;


  bool isRefitNoBS=true;
  bool isRefitBS=true;
  bool isRefitNoBSZNominal=true;
  bool isRefitBSZNominal=true;

  double Spin_WT=Ntp->TauSpinerGet(TauSpinerInterface::Spin);
  // double UnSpin_WT=Ntp->TauSpinerGet(TauSpinerInterface::UnSpin);
  double FlipSpin_WT=Ntp->TauSpinerGet(TauSpinerInterface::FlipSpin);
  // double hplus=Ntp->TauSpinerGet(TauSpinerInterface::hplus);
  // double hminus=Ntp->TauSpinerGet(TauSpinerInterface::hminus);//1-hplus;
		 

  double Wspin;
  double Wflipspin;
  if(id == 11 || id == 12|| id == 45 || id ==461||id == 460)Wspin=w*Spin_WT;
  else Wspin=w;
  if(id == 11 || id == 12|| id == 45 ||id ==461 ||id ==460)Wflipspin=w*FlipSpin_WT;
  else Wflipspin=w;
  TLorentzVector testTruth(0,0,0,0);

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

  SCalculator Scalc("a1");
  SCalculator ScalcPVRefitNoBS("a1");
  SCalculator ScalcPVRefitBS("a1");
  SCalculator ScalcPVRefitNoBSZNominal("a1");
  SCalculator ScalcPVRefitBSZNominal("a1");
    
  SCalculator ScalcPVRefitWithTracksBS("a1");
  SCalculator ScalcPVRefitWithTracksBSZNominal("a1");
  double acop = -99.;
  
  TLorentzVector zeroLV(0,0,0,0);
  std::vector<TLorentzVector> VectZeroLV;
  VectZeroLV.push_back(zeroLV);
  VectZeroLV.push_back(zeroLV);
  VectZeroLV.push_back(zeroLV);


  TLorentzVector Tauplusvis;
  TLorentzVector Tauminusvis;
  TLorentzVector Pi0RECO;
  TLorentzVector Tauplustruth;
  TLorentzVector Tauminustruth;
  unsigned int Tauplus=0;
  unsigned int Tauminus=0;

  if(passAllBut(exclude_cuts)) { 
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

    if(Ntp->decayMode(Tauminus) == 10 &&  Ntp->PFTau_hassecondaryVertex(Tauminus) && Ntp->PFtauHasThreePions(Tauminus))a1minus=true;
    if(Ntp->decayMode(Tauplus) == 10 && Ntp->PFTau_hassecondaryVertex(Tauplus) && Ntp->PFtauHasThreePions(Tauplus))a1plus=true;
    
    if(Ntp->MVADM2017(Tauminus) == 10 &&  Ntp->PFTau_hassecondaryVertex(Tauminus) && Ntp->PFtauHasThreePions(Tauminus))a1minusMVA=true;
    if(Ntp->MVADM2017(Tauplus) == 10 && Ntp->PFTau_hassecondaryVertex(Tauplus) && Ntp->PFtauHasThreePions(Tauplus))a1plusMVA=true;
    
    if(a1minus && a1plus ) a1a1=true;  //a1-a1
    if(a1minusMVA && a1plusMVA ) a1a1MVA=true;  //a1-a1

    TLorentzVector Tauplussvfit;
    TLorentzVector Tauminussvfit;

    ClassicSVfit svfitAlgo1(0);
    double higgsmass;
    if(a1a1MVA) {

      //FastMTT FastMTTAlgo;
      //if(a1a1TruthSVFit || a1a1TruthSVFitMVA /*|| a1rhoTruthSVFit || a1rhoTruthSVFitMVA ||a1piTruthSVFit || a1piTruthSVFitMVA*/)
      //{
      // // //---------  svfit ---------------------
      std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;
      classic_svFit::MeasuredTauLepton lep1(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Tau1P4.Pt(), Tau1P4.Eta(),  Tau1P4.Phi(), Tau1P4.M(),Ntp->MVADM2017(Tau1));
      classic_svFit::MeasuredTauLepton lep2(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Tau2P4.Pt(), Tau2P4.Eta(),  Tau2P4.Phi(), Tau2P4.M(),Ntp->MVADM2017(Tau2));

      measuredTauLeptons.push_back(lep1);
      measuredTauLeptons.push_back(lep2);
      TMatrixD metcov(2,2);
      double metx = PUPPIMET*cos(PUPPIMETPhi);
      double mety = PUPPIMET*sin(PUPPIMETPhi);

      metcov[0][0] = Ntp->PUPPIMETCov00();
      metcov[1][0] = Ntp->PUPPIMETCov01();
      metcov[0][1] = Ntp->PUPPIMETCov10();
      metcov[1][1] = Ntp->PUPPIMETCov11();

      svfitAlgo1.setHistogramAdapter(new classic_svFit::TauTauHistogramAdapter());
      svfitAlgo1.addLogM_fixed(true,5.);
      // svfitAlgo1.setDiTauMassConstraint(125.10);
      //FastMTTAlgo.run(measuredTauLeptons, metx, mety, metcov);
      svfitAlgo1.integrate(measuredTauLeptons,metx,mety, metcov );
      if(svfitAlgo1.isValidSolution()) {

	higgsmass  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svfitAlgo1.getHistogramAdapter())->getMass();
	//h_SVFitMass.at(t).Fill(higgsmass,w); 
	//tau1P4 = FastMTTAlgo.getTau1P4();
	//tau2P4 = FastMTTAlgo.getTau2P4();
    
	tau1P4 = static_cast<classic_svFit::TauTauHistogramAdapter*>(svfitAlgo1.getHistogramAdapter())->GetFittedTau1LV();	  
	tau2P4 = static_cast<classic_svFit::TauTauHistogramAdapter*>(svfitAlgo1.getHistogramAdapter())->GetFittedTau2LV();
	  

	// ClassicSVfit svfitAlgo2;
	// svfitAlgo2.setHistogramAdapter(new classic_svFit::TauTauHistogramAdapter());
	// svfitAlgo2.addLogM_fixed(true, 5.);
	// svfitAlgo2.integrate(measuredTauLeptons,metx,mety, metcov );
	// tau1P4 = static_cast<classic_svFit::TauTauHistogramAdapter*>(svfitAlgo2.getHistogramAdapter())->GetFittedTau1LV();
	// tau2P4 = static_cast<classic_svFit::TauTauHistogramAdapter*>(svfitAlgo2.getHistogramAdapter())->GetFittedTau2LV();
	
	// //---------  svfit ---------------------
	//if(svfitAlgo1.isValidSolution()){
	if(Ntp->Daughters_charge(Tau1)>0)
	  {
	    Tauplussvfit.SetPxPyPzE(tau1P4.x(),tau1P4.y(),tau1P4.z(),tau1P4.t());
	    Tauminussvfit.SetPxPyPzE(tau2P4.x(),tau2P4.y(),tau2P4.z(),tau2P4.t());
	
	  }
	else
	  {
	    Tauplussvfit.SetPxPyPzE(tau2P4.x(),tau2P4.y(),tau2P4.z(),tau2P4.t());
	    Tauminussvfit.SetPxPyPzE(tau1P4.x(),tau1P4.y(),tau1P4.z(),tau1P4.t());
	  }
      }
      //    }
      // }

      if(JetsVectJERUp.size()>1){
	if (JERUp && (JetsVectJERUp[0].Pt()>30. && JetsVectJERUp[1].Pt()>30.)){
	  jdeta_=Ntp->jdetaUp(IndexSelectedTemp.at(Sorted.back()));
	  mjj_=Ntp->mjjUp(IndexSelectedTemp.at(Sorted.back()));
	}
      }
      if(JetsVectJERDown.size()>1){
	if (JERDown && (JetsVectJERDown[0].Pt()>30. && JetsVectJERDown[1].Pt()>30.)){
	  jdeta_=Ntp->jdetaDown(IndexSelectedTemp.at(Sorted.back()));
	  mjj_=Ntp->mjjDown(IndexSelectedTemp.at(Sorted.back()));
	}
      }
      if(JetsVectJES.size()>1){
	if((JESUp || JESDown) && (JetsVectJES[0].Pt()>30. && JetsVectJES[1].Pt()>30.)){
	  jdeta_=Ntp->jdeta(IndexSelectedTemp.at(Sorted.back()));
	  mjj_=(JetsVectJES[0]+JetsVectJES[1]).M();
	}
	else if(JetsVectJES[0].Pt()>30. && JetsVectJES[1].Pt()>30.){
	  jdeta_=Ntp->jdeta(IndexSelectedTemp.at(Sorted.back()));
	  mjj_=(JetsVectJES[0]+JetsVectJES[1]).M();
	}
      }

      //cout<<Ntp->jptSize_1()<<"  "<<IndexSelectedTemp.at(Sorted.back())<<"  "<<njets20<<"  "<<Ntp->jets_jetUnc_AbsoluteSize_up()<<endl;
      if(JetsVectJERUp.size()>0){
	if (JERUp && JetsVectJERUp[0].Pt()>30.){
	  jpt_1_=JetsVectJERUp[0].Pt();
	}
      }
      if(JetsVectJERDown.size()>0){
	if (JERDown && JetsVectJERDown[0].Pt()>30.){
	  jpt_1_=JetsVectJERDown[0].Pt();
	}
      }
      if(JetsVectJES.size()>0){
	if ((JESUp || JESDown) && JetsVectJES[0].Pt()>30.){ 
	  jpt_1_=JetsVectJES[0].Pt();
	}
	else if(JetsVectJES[0].Pt()>30.){
	  jpt_1_=JetsVectJES[0].Pt();
	}
      }
      
      if(svfitAlgo1.isValidSolution()) m_sv_=higgsmass;
      else m_sv_=-1;
      m_vis_=(Tau1P4+Tau2P4).M();
      met_=PUPPIMET;

      pt_1_=Tau1P4.Pt();
      //pt_2_=Tau2P4.Pt();
      pt_tt_=(Tau1P4+Tau2P4+TLorentzVector(PUPPImetCorr_px,PUPPImetCorr_py,0,0)).Pt();
      pt_vis_=(Tau1P4+Tau2P4).Pt();
      //cout<<jdeta_<<" "<<jpt_1_<<"  "<<m_vis_<<"  "<<met_<<"  "<<mjj_<<"  "<<n_jets_<<endl;

      JetsVect.clear();
      JetsVectJERUp.clear();
      JetsVectJERDown.clear();
      JetsVectJES.clear();
      JetsVectNominal.clear();

      BDT->Execute(jdeta_, jpt_1_, m_vis_, met_, mjj_, n_jets_, pt_1_, pt_tt_, pt_vis_, m_sv_, evt_,scores,max_pair); 

      // cout<<"--------------------"<<endl;
      // 	for(unsigned int i=0; i<3; i++){
      // 	  cout<<Ntp->PFTauRefit_PionsCharge(Tauplus,i)<<endl;
      // 	}
      // 	for(unsigned int i=0; i<3; i++){
      // 	cout<<Ntp->PFTauRefit_PionsCharge(Tauminus,i)<<endl;
      // 	}
      // cout<<svfitAlgo1.isValidSolution()<<endl;
      // Acoplanar Angle
    }

    if(a1a1MVA && std::isnan(Wspin)!=true && Ntp->PFTauRefit_PionsP4_SizePions(Tauminus)==3 && Ntp->PFTauRefit_PionsP4_SizePions(Tauplus)==3)
      {
	// 	if((!Ntp->isData() || isEmbed) && GenMatchSelection){
	// 	  if(Ntp->byVVVLooseDeepTau2017v2p1VSjet_1(0))Tau1isolation.at(t).Fill(0.,w);
	// 	  if(Ntp->byVVLooseDeepTau2017v2p1VSjet_1(0))Tau1isolation.at(t).Fill(1,w);
	// 	  if(Ntp->byVLooseDeepTau2017v2p1VSjet_1(0))Tau1isolation.at(t).Fill(2,w);
	// 	  if(Ntp->byLooseDeepTau2017v2p1VSjet_1(0))Tau1isolation.at(t).Fill(3,w);
	// 	  if(Ntp->byMediumDeepTau2017v2p1VSjet_1(0))Tau1isolation.at(t).Fill(4,w);
	// 	  if(Ntp->byTightDeepTau2017v2p1VSjet_1(0))Tau1isolation.at(t).Fill(5,w);
	// 	  if(Ntp->byVTightDeepTau2017v2p1VSjet_1(0))Tau1isolation.at(t).Fill(6,w);
	// 	  if(Ntp->byVVTightDeepTau2017v2p1VSjet_1(0))Tau1isolation.at(t).Fill(7,w);
	// 	  if(Ntp->byVVVLooseDeepTau2017v2p1VSjet_2(0))Tau2isolation.at(t).Fill(0.,w);
	// 	  if(Ntp->byVVLooseDeepTau2017v2p1VSjet_2(0))Tau2isolation.at(t).Fill(1,w);
	// 	  if(Ntp->byVLooseDeepTau2017v2p1VSjet_2(0))Tau2isolation.at(t).Fill(2,w);
	// 	  if(Ntp->byLooseDeepTau2017v2p1VSjet_2(0))Tau2isolation.at(t).Fill(3,w);
	// 	  if(Ntp->byMediumDeepTau2017v2p1VSjet_2(0))Tau2isolation.at(t).Fill(4,w);
	// 	  if(Ntp->byTightDeepTau2017v2p1VSjet_2(0))Tau2isolation.at(t).Fill(5,w);
	// 	  if(Ntp->byVTightDeepTau2017v2p1VSjet_2(0))Tau2isolation.at(t).Fill(6,w);
	// 	  if(Ntp->byVVTightDeepTau2017v2p1VSjet_2(0))Tau2isolation.at(t).Fill(7,w);
	// 	}
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

	tauWithTracksBSPrimaryVertex = TVector3(Ntp->RefitPVWithTracksBS_x(),Ntp->RefitPVWithTracksBS_y(),Ntp->RefitPVWithTracksBS_z());
	tauWithTracksBSZNominalPrimaryVertex = TVector3(Ntp->RefitPVWithTracksBS_x(),Ntp->RefitPVWithTracksBS_y(),Ntp->PVtx().Z());


	TauminusDirection = TauminusSecondaryVertex - tauPrimaryVertex;
	TauplusDirection = TauplusSecondaryVertex - tauPrimaryVertex;

	solutionsNoBS=tauPairMomentumSolutions(TauminusSecondaryVertex-tauNoBSPrimaryVertex, a1LVRefit_Tauminus, a1LV_Tauminus, isMinusReal, TauplusSecondaryVertex-tauNoBSPrimaryVertex, a1LVRefit_Tauplus, a1LVRefit_Tauplus, isPlusReal,isRefitNoBS); 
	solutionsBS=tauPairMomentumSolutions(TauminusSecondaryVertex-tauBSPrimaryVertex, a1LVRefit_Tauminus, a1LV_Tauminus, isMinusReal, TauplusSecondaryVertex-tauBSPrimaryVertex, a1LVRefit_Tauplus, a1LVRefit_Tauplus, isPlusReal,isRefitBS); 
	solutionsNoBSZNominal=tauPairMomentumSolutions(TauminusSecondaryVertex-tauNoBSZNominalPrimaryVertex, a1LVRefit_Tauminus, a1LV_Tauminus, isMinusReal, TauplusSecondaryVertex-tauNoBSZNominalPrimaryVertex, a1LVRefit_Tauplus, a1LVRefit_Tauplus, isPlusReal,isRefitNoBSZNominal); 
	solutionsBSZNominal=tauPairMomentumSolutions(TauminusSecondaryVertex-tauBSZNominalPrimaryVertex, a1LVRefit_Tauminus, a1LV_Tauminus, isMinusReal, TauplusSecondaryVertex-tauBSZNominalPrimaryVertex, a1LVRefit_Tauplus, a1LVRefit_Tauplus, isPlusReal,isRefitBSZNominal); 

	solutions=tauPairMomentumSolutions(TauminusDirection, a1LVRefit_Tauminus, a1LV_Tauminus, isMinusReal, TauplusDirection, a1LVRefit_Tauplus, a1LV_Tauplus, isPlusReal,false);
	
	solutionsWithTracksBS=tauPairMomentumSolutions(TauminusSecondaryVertex-tauWithTracksBSPrimaryVertex, a1LVRefit_Tauminus, a1LV_Tauminus, isMinusReal, TauplusSecondaryVertex-tauWithTracksBSPrimaryVertex, a1LVRefit_Tauplus, a1LVRefit_Tauplus, isPlusReal,true); 
	solutionsWithTracksBSZNominal=tauPairMomentumSolutions(TauminusSecondaryVertex-tauWithTracksBSZNominalPrimaryVertex, a1LVRefit_Tauminus, a1LV_Tauminus, isMinusReal, TauplusSecondaryVertex-tauBSZNominalPrimaryVertex, a1LVRefit_Tauplus, a1LVRefit_Tauplus, isPlusReal,true); 

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

	TauminusPairConstraintWithTracksBS=solutionsWithTracksBS.at(3);
	TauminusPairConstraintWithTracksBSZNominal=solutionsWithTracksBSZNominal.at(3);
	TauplusPairConstraintWithTracksBS=solutionsWithTracksBS.at(7);
	TauplusPairConstraintWithTracksBSZNominal=solutionsWithTracksBSZNominal.at(7);
	
	HPairConstraint= TauplusPairConstraint+TauminusPairConstraint;


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

	//DPM---------------------------------------------------
	//////////////////////
        //Decay Plane Method//
        //////////////////////
        
    	unsigned int idx_PiMinus=99, idx_PiPlus = 99;
    	unsigned int idx_PiPlusfromRho=99, idx_PiMinusfromRho = 99;
        unsigned int idx_PiPlusfromA1=99, idx_PiMinusfromA1 = 99;

    	Double_t y_Tau;
    	Double_t RhoZero_mass = 0.77526; //GeV

    	std::vector<Double_t> idx;
    	std::vector<Double_t> InvMass;

    	////////////////////////
    	///////TAU MINUS////////
    	////////////////////////
        //std::cout <<" DPM 0 "<<Tauminus << " " << Tauplus << std::endl;
 
    	//Look for the PionPlus index associated with TauMinus
    	for(unsigned int i=0; i<3; i++){
	  //std::cout << "aloa " << Ntp->PFTauRefit_PionsCharge(Tauminus,i) <<  " "<< Ntp->PFTau_PionsCharge(Tauminus,i) << std::endl;
    	  if(Ntp->PFTauRefit_PionsCharge(Tauminus,i)==+1){
	    if (Ntp->PFTauRefit_PionsCharge(Tauminus,i)==0) std::cout <<" BIG BUG "<<std::endl;
    	    idx_PiPlus = i;
    	  }
    	}

	if(idx_PiPlus != 99)
	  {
	    //Calculate inv. mass of the two possible pions pairs for TauMinus
	    for(unsigned int i=0; i<3; i++){
	      if(i==idx_PiPlus){continue;}
	      else{
		InvMass.push_back((Ntp->PFTauRefit_PionsP4(Tauminus,idx_PiPlus) + Ntp->PFTauRefit_PionsP4(Tauminus,i)).Mag());
		idx.push_back(i);
	      }
	    }
 
	    //Look for the closest inv. mass from rho 0 mass for TauMinus
	    if(fabs(RhoZero_mass-InvMass[0])<fabs(RhoZero_mass-InvMass[1])){idx_PiMinusfromRho = idx[0];idx_PiMinusfromA1 = idx[1];}
	    else{idx_PiMinusfromRho = idx[1];idx_PiMinusfromA1 = idx[0];}

	  }
    	//Clear vectors
    	idx.clear();
    	InvMass.clear();

    	////////////////////////
    	////////////////////////

    	////////////////////////
    	////////TAU PLUS////////
    	////////////////////////

    	//Look for the PionMinus index associated with TauPlus
    	for(unsigned int i=0; i<3; i++){
	  //std::cout << "aloa " << Ntp->PFTauRefit_PionsCharge(Tauplus,i) <<  " "<< Ntp->PFTau_PionsCharge(Tauplus,i) << std::endl;
    	  if(Ntp->PFTauRefit_PionsCharge(Tauplus,i)==-1 ){
	    if (Ntp->PFTauRefit_PionsCharge(Tauplus,i)==0) std::cout <<" BIG BUG "<<std::endl;
    	    idx_PiMinus = i;
    	    //break;
    	  }
    	}
	
	if(idx_PiMinus != 99)
	  {
	    //Calculate inv. mass of the two possible pions pairs for TauPlus
	    for(unsigned int i=0; i<3; i++){
	      if(i==idx_PiMinus){continue;}
	      else{
		InvMass.push_back((Ntp->PFTauRefit_PionsP4(Tauplus,idx_PiMinus) + Ntp->PFTauRefit_PionsP4(Tauplus,i)).Mag());
		idx.push_back(i);
	      }
	    }
	    //Look for the closest inv. mass from rho 0 mass for TauPlus
	    if(fabs(RhoZero_mass-InvMass[0])<fabs(RhoZero_mass-InvMass[1])){idx_PiPlusfromRho = idx[0];idx_PiPlusfromA1 = idx[1];}
	    else{idx_PiPlusfromRho = idx[1];idx_PiPlusfromA1 = idx[0];}
	  }

    	//Clear vectors
    	idx.clear();
    	InvMass.clear();
	
	if(idx_PiMinus != 99 && idx_PiPlus != 99)
	  {
	    TLorentzVector Pi1 = Ntp->PFTauRefit_PionsP4(Tauminus,idx_PiMinusfromRho);
	    TLorentzVector Pi2 = Ntp->PFTauRefit_PionsP4(Tauplus,idx_PiPlusfromRho);

	    TLorentzVector ref1 = Ntp->PFTauRefit_PionsP4(Tauminus,idx_PiPlus);
	    TLorentzVector ref2 = Ntp->PFTauRefit_PionsP4(Tauplus,idx_PiMinus);
      
	    TLorentzVector Pi1a1 = Ntp->PFTauRefit_PionsP4(Tauminus,idx_PiMinusfromA1);
	    TLorentzVector Pi2a1 = Ntp->PFTauRefit_PionsP4(Tauplus,idx_PiPlusfromA1);


	    double y1 = 1;
	    double y2 = 1;

	    y1 = (Pi1.E() - ref1.E())/(Pi1.E() + ref1.E());

	    y2 = (Pi2.E() - ref2.E())/(Pi2.E() + ref2.E());

	    y_Tau = y1*y2;

	    TLorentzVector Prongsum = TauminusPairConstraintWithTracksBS + TauplusPairConstraintWithTracksBS;
	    TVector3 boost = -Prongsum.BoostVector();

	    TLorentzVector Rho1 = Pi1+ref1;
	    TLorentzVector Rho2 = Pi2+ref2;

	    Pi1.Boost(boost);
	    Pi2.Boost(boost);
	    ref1.Boost(boost);
	    ref2.Boost(boost);
	    Pi1a1.Boost(boost);
	    Pi2a1.Boost(boost);
	    Rho1.Boost(boost);
	    Rho2.Boost(boost);

	    TVector3 vecPi1 = Pi1.Vect();
	    TVector3 vecPi2 = Pi2.Vect();
	    TVector3 vecRef1 = ref1.Vect();
	    TVector3 vecRef2 = ref2.Vect();
	    TVector3 vecPi1a1 = Pi1a1.Vect();
	    TVector3 vecPi2a1 = Pi2a1.Vect();
	    TVector3 vecRho1 = Rho1.Vect();
	    TVector3 vecRho2 = Rho2.Vect();
       
	    vecPi1 *= 1/vecPi1.Mag();
	    vecPi2 *= 1/vecPi2.Mag();
	    vecRef1 *= 1/vecRef1.Mag();
	    vecRef2 *= 1/vecRef2.Mag();
	    vecPi1a1 *= 1/vecPi1a1.Mag();
	    vecPi2a1 *= 1/vecPi2a1.Mag();
	    vecRho1 *= 1/vecRho1.Mag();
	    vecRho2 *= 1/vecRho2.Mag();
    
	    // transverse components
	    TVector3 vecRef1transv = vecRef1 - vecPi1*(vecPi1*vecRef1);
	    TVector3 vecRef2transv = vecRef2 - vecPi2*(vecPi2*vecRef2);

	    vecRef1transv *= 1/vecRef1transv.Mag();
	    vecRef2transv *= 1/vecRef2transv.Mag();

	    acop = TMath::ACos(vecRef1transv*vecRef2transv);
	    double acoporiginal=acop;
	    double sign = vecPi2 * vecRef1transv.Cross(vecRef2transv);
	    double psioriginal =sign;

	    if (sign<0) acop = 2.0*TMath::Pi() - acop;

	    if (y_Tau<0) {
	      acop = acop + TMath::Pi();
	      if (acop>2*TMath::Pi()) {
		acop = acop - 2*TMath::Pi();
	      }
	    }
	  }//if  idx_PiMinus != 99 && idx_PiPlus != 99)
	
      }
    RooWorkspace* wFF=NULL;
    double angle=-99;
    if (std::isnan(Wspin)!=true && a1a1MVA)
      {
	if(Ntp->year()==2016)wFF=wFF2016;
	else if(Ntp->year()==2017)wFF=wFF2017;
	else if(Ntp->year()==2018)wFF=wFF2018;

	if((HadRefitPions_minus!=HadRefitPions_plus) && (HadRefitPions_minus!=VectZeroLV) && (HadRefitPions_plus!=VectZeroLV) && TauminusPairConstraintWithTracksBS!=zeroLV && TauplusPairConstraintWithTracksBS!=zeroLV && ScalcPVRefitWithTracksBS.isOk("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus)){
	  
	  int idx=IndexSelectedTemp.at(Sorted.back());
	  angle= ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus);
	  bool isOS=(Ntp->Daughters_charge(Tau1)/abs(Ntp->Daughters_charge(Tau1))) != (Ntp->Daughters_charge(Tau2)/abs(Ntp->Daughters_charge(Tau2)));
	  auto args = std::vector<double>{Tau1P4.Pt(),(double)Ntp->MVADM2017(Tau1),(double)n_jets_,(double)isOS,PUPPIMET*cos(TLorentzVector(Tau1P4.X(),Tau1P4.Y(),Tau1P4.Z(),0.).DeltaPhi(TLorentzVector(PUPPImetCorr_px,PUPPImetCorr_py,0.,0.)))/Tau1P4.Pt(),Tau1P4.DeltaR(Tau2P4)};


	  double FFData= std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data());	
	  double FFMC= std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data());

	  std::vector<TH2D> emptyvector;
	  //cout<<"angle: "<<angle<<"  Wspin: "<<Wspin<<endl;
	  Ntp->FillHist(t,idx,isOS, GenMatchSelection,angle, max_pair ,Wspin, FFData,FFMC,&polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggs,&polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakes,&polarimetricAcopAnglePVRefitWithTracksBSMVADMZTT, &polarimetricAcopAnglePVRefitWithTracksBSMVADMWfakesHiggs,&polarimetricAcopAnglePVRefitWithTracksBSMVADMWfakesJetFakes,&polarimetricAcopAnglePVRefitWithTracksBSMVADMWfakesZTT, &polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsQCDMC,&polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesQCDMC,&polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTQCDMC);

	  Ntp->FillHist(t,idx,isOS, GenMatchSelection,angle, max_pair ,Wspin, FFData,FFMC,&ShapeSystPVRefitWithTracksBSHiggs,&ShapeSystPVRefitWithTracksBSJetFakes,&ShapeSystPVRefitWithTracksBSZTT,& ShapeSystPVRefitWithTracksBSWfakesHiggs,&ShapeSystPVRefitWithTracksBSWfakesJetFakes,&ShapeSystPVRefitWithTracksBSWfakesZTT,& emptyvector,&emptyvector,&emptyvector);
	  /*
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection,angle, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggs,& jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakes,& jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTT,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC,& jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC,& jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection,angle, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggs,&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakes ,& jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTT,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC,&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC ,& jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC);

	    Ntp->FillHist(t,idx,isOS, GenMatchSelection,angle, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggs,& jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakes,& jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTT,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC,& jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC,& jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection, angle, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggs,&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakes ,&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTT,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC,&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC ,&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection,angle, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggs,&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakes ,& jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTT,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC,&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC ,& jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection,angle, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggs,&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakes ,& jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTT,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC,&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC ,& jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection,angle, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggs,& jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakes,& jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTT,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC,& jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC,& jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection,angle, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggs,& jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakes,& jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTT,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC,& jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC,& jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection, angle, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggs,& jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakes,& jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTT,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC,& jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC,& jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection, angle, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggs,&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakes ,& jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTT,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC,&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC ,& jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection,angle, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggs,& jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakes,& jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTT,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC,& jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC,& jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection,angle, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggs,&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakes ,&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTT,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC,&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC ,&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection, angle, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSHiggs,& jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSJetFakes,&jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSZTT,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSHiggsQCDMC,& jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSJetFakesQCDMC,&jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSZTTQCDMC);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection,angle, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSHiggs,&jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSJetFakes ,& jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSZTT,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSHiggsQCDMC,&jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSJetFakesQCDMC ,& jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSZTTQCDMC);

	    Ntp->FillHist(t,idx,isOS, GenMatchSelection,angle, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSHiggs,& jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSJetFakes,& jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSZTT,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSHiggsQCDMC,& jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSJetFakesQCDMC,& jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSZTTQCDMC);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection,angle, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSHiggs,& jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSJetFakes,& jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSZTT,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSHiggsQCDMC,& jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSJetFakesQCDMC,& jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSZTTQCDMC);
	  */
	  Ntp->FillHist(t,idx,isOS, GenMatchSelection, angle, max_pair, (Wspin/(wIDSF1*wIDSF2))*Ntp->IDSF(Tau1,GenMatch1,"Nom","tau","Up")*Ntp->IDSF(Tau2,GenMatch2,"Nom","tau","Up"),FFData,FFMC,&CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSHiggs,&CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSJetFakes,& CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSZTT,&CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSWfakesHiggs,&CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSWfakesJetFakes,& CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSWfakesZTT,&emptyvector,&emptyvector,&emptyvector);
	  Ntp->FillHist(t,idx,isOS, GenMatchSelection,angle, max_pair,(Wspin/(wIDSF1*wIDSF2))*Ntp->IDSF(Tau1,GenMatch1,"Nom","tau","Down")*Ntp->IDSF(Tau2,GenMatch2,"Nom","tau","Down") , FFData,FFMC,&CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSHiggs,&CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSJetFakes,& CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSZTT,&CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSWfakesHiggs,&CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSWfakesJetFakes,& CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSWfakesZTT,&emptyvector,&emptyvector,&emptyvector);

	  Ntp->FillHist(t,idx,isOS, GenMatchSelection, angle, max_pair, (Wspin/(wTrgSF1*wTrgSF2))*Ntp->TriggerSF(Tau1,GenMatch1,"Nom","Up")*Ntp->TriggerSF(Tau2,GenMatch2,"Nom","Up"),FFData ,FFMC,&CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSHiggs,&CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSJetFakes,& CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSZTT,&CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSWfakesHiggs,&CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSWfakesJetFakes,& CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSWfakesZTT,&emptyvector,&emptyvector,&emptyvector);
	  Ntp->FillHist(t,idx,isOS, GenMatchSelection,angle, max_pair,(Wspin/(wTrgSF1*wTrgSF2))*Ntp->TriggerSF(Tau1,GenMatch1,"Nom","Down")*Ntp->TriggerSF(Tau2,GenMatch2,"Nom","Down"),FFData,FFMC,&CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSHiggs,&CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSJetFakes,& CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSZTT,&CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSWfakesHiggs,&CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSWfakesJetFakes,& CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSWfakesZTT,&emptyvector,&emptyvector,&emptyvector);

	  Ntp->FillHist(t,idx,isOS, GenMatchSelection, angle, max_pair ,zptw*zptw*(Wspin/zptw),FFData,FFMC,&CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSHiggs,&CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSJetFakes,& CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSZTT,&CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSWfakesHiggs,&CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSWfakesJetFakes,& CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSWfakesZTT,&emptyvector,&emptyvector,&emptyvector);
	  Ntp->FillHist(t,idx,isOS, GenMatchSelection, angle, max_pair, Wspin/zptw, FFData, FFMC,& CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSHiggs,&CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSJetFakes,& CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSZTT,&CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSWfakesHiggs,&CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSWfakesJetFakes,& CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSWfakesZTT,&emptyvector,&emptyvector,&emptyvector);

	  
	  if(id==11||id==12||id==45||id==460||id==461){
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection, angle, max_pair ,Wspin*Ntp->TheoreticalPSUnc(6)/Ntp->MC_weight(),FFData,FFMC,&CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSHiggs,&CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSJetFakes,& CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSZTT,&CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSWfakesHiggs,&CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSWfakesJetFakes,& CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSWfakesZTT,&emptyvector,&emptyvector,&emptyvector);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection, angle, max_pair, Wspin*Ntp->TheoreticalPSUnc(8)/Ntp->MC_weight(), FFData, FFMC,& CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSHiggs,&CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSJetFakes,& CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSZTT,&CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSWfakesHiggs,&CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSWfakesJetFakes,& CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSWfakesZTT,&emptyvector,&emptyvector,&emptyvector);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection, angle, max_pair ,Wspin*Ntp->TheoreticalPSUnc(7)/Ntp->MC_weight(),FFData,FFMC,&CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSHiggs,&CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSJetFakes,& CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSZTT,&CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSWfakesHiggs,&CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSWfakesJetFakes,& CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSWfakesZTT,&emptyvector,&emptyvector,&emptyvector);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection, angle, max_pair, Wspin*Ntp->TheoreticalPSUnc(9)/Ntp->MC_weight(), FFData, FFMC,& CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSHiggs,&CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSJetFakes,& CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSZTT,&CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSWfakesHiggs,&CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSWfakesJetFakes,& CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSWfakesZTT,&emptyvector,&emptyvector,&emptyvector);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection, angle, max_pair ,Wspin*Ntp->TheoreticalScaleUnc1005()/Ntp->nominal_wt(),FFData,FFMC,&CMS_scale_gg_13TeVUpPVRefitWithTracksBSHiggs,&CMS_scale_gg_13TeVUpPVRefitWithTracksBSJetFakes,& CMS_scale_gg_13TeVUpPVRefitWithTracksBSZTT,&CMS_scale_gg_13TeVUpPVRefitWithTracksBSWfakesHiggs,&CMS_scale_gg_13TeVUpPVRefitWithTracksBSWfakesJetFakes,& CMS_scale_gg_13TeVUpPVRefitWithTracksBSWfakesZTT,&emptyvector,&emptyvector,&emptyvector);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection, angle, max_pair, Wspin*Ntp->TheoreticalScaleUnc1009()/Ntp->nominal_wt(), FFData, FFMC,& CMS_scale_gg_13TeVDownPVRefitWithTracksBSHiggs,&CMS_scale_gg_13TeVDownPVRefitWithTracksBSJetFakes,& CMS_scale_gg_13TeVDownPVRefitWithTracksBSZTT,&CMS_scale_gg_13TeVDownPVRefitWithTracksBSWfakesHiggs,&CMS_scale_gg_13TeVDownPVRefitWithTracksBSWfakesJetFakes,& CMS_scale_gg_13TeVDownPVRefitWithTracksBSWfakesZTT,&emptyvector,&emptyvector,&emptyvector);
	  
	  }
	  Ntp->FillHist(t,idx,isOS, GenMatchSelection, angle, max_pair,top_wt*top_wt*(Wspin/top_wt) ,FFData,FFMC,&CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSHiggs,&CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSJetFakes,&CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSZTT,&CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSWfakesHiggs,&CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSWfakesJetFakes,&CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSWfakesZTT ,&emptyvector,&emptyvector,&emptyvector);
	  Ntp->FillHist(t,idx,isOS, GenMatchSelection, angle, max_pair,Wspin/top_wt ,FFData,FFMC,&CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSHiggs,&CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSJetFakes,&CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSZTT,&CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSWfakesHiggs,&CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSWfakesJetFakes,&CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSWfakesZTT ,&emptyvector,&emptyvector,&emptyvector);
	  

	  if(id==71||id==72||id==73||id==74 || id==47||id==48||id==49||id==50||id==51||id==52||id==53||id==54||id==55||id==56||id==57||id==58 || id==70 ||id==701 ||id==702 ||id==703){ // VV+TT
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection, angle, max_pair,Wspin ,FFData,FFMC,&CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSHiggs,&CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSJetFakes,&CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSZTT,&CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesHiggs,&CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesJetFakes,&CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesZTT ,&emptyvector,&emptyvector,&emptyvector);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection, angle, max_pair,Wspin ,FFData,FFMC,&CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSHiggs,&CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSJetFakes,&CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSZTT,&CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesHiggs,&CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesJetFakes,&CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesZTT ,&emptyvector,&emptyvector,&emptyvector);
	    GenMatchSelection=!(GenMatch1==6 || GenMatch2==6);
	    if (GenMatch1==5 && GenMatch2==5){
	      Ntp->FillHist(t,idx,isOS, GenMatchSelection, angle, max_pair,Wspin,FFData,FFMC,&ttbarcontaminationHiggs,&ttbarcontaminationJetFakes,&ttbarcontaminationZTT,&ttbarcontaminationWfakesHiggs,&ttbarcontaminationWfakesJetFakes,&ttbarcontaminationWfakesZTT ,&emptyvector,&emptyvector,&emptyvector);
	    }
	    GenMatchSelection=(!(GenMatch1==6 || GenMatch2==6) &&  !(GenMatch1==5 && GenMatch2==5));
	      
	  }
			      
	  Ntp->FillHist(t,idx,isOS, GenMatchSelection, angle, max_pair, prefup*(Wspin/Ntp->prefiringweight()),FFData, FFMC,&PrefiringUpPVRefitWithTracksBSHiggs,&PrefiringUpPVRefitWithTracksBSJetFakes,& PrefiringUpPVRefitWithTracksBSZTT,&PrefiringUpPVRefitWithTracksBSWfakesHiggs,&PrefiringUpPVRefitWithTracksBSWfakesJetFakes,& PrefiringUpPVRefitWithTracksBSWfakesZTT ,&emptyvector,&emptyvector,&emptyvector);
	  Ntp->FillHist(t,idx,isOS, GenMatchSelection, angle, max_pair, prefdown*(Wspin/Ntp->prefiringweight()),FFData,FFMC,& PrefiringDownPVRefitWithTracksBSHiggs,&PrefiringDownPVRefitWithTracksBSJetFakes,& PrefiringDownPVRefitWithTracksBSZTT,&PrefiringDownPVRefitWithTracksBSWfakesHiggs,&PrefiringDownPVRefitWithTracksBSWfakesJetFakes,& PrefiringDownPVRefitWithTracksBSWfakesZTT,&emptyvector,&emptyvector,&emptyvector);
	
	  
	  
	  //DPM//

	  

	  Ntp->FillHist(t,idx,isOS, GenMatchSelection,acop, max_pair ,Wspin, FFData,FFMC,&polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggs_DP,&polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakes_DP,&polarimetricAcopAnglePVRefitWithTracksBSMVADMZTT_DP,& polarimetricAcopAnglePVRefitWithTracksBSMVADMWfakesHiggs_DP,&polarimetricAcopAnglePVRefitWithTracksBSMVADMWfakesJetFakes_DP,&polarimetricAcopAnglePVRefitWithTracksBSMVADMWfakesZTT_DP,& polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsQCDMC_DP,&polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesQCDMC_DP,&polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTQCDMC_DP);

	  Ntp->FillHist(t,idx,isOS, GenMatchSelection,acop, max_pair ,Wspin, FFData,FFMC,&ShapeSystPVRefitWithTracksBSHiggs_DP,&ShapeSystPVRefitWithTracksBSJetFakes_DP,&ShapeSystPVRefitWithTracksBSZTT_DP,& ShapeSystPVRefitWithTracksBSWfakesHiggs_DP,&ShapeSystPVRefitWithTracksBSWfakesJetFakes_DP,&ShapeSystPVRefitWithTracksBSWfakesZTT_DP,& emptyvector,&emptyvector,&emptyvector);
	  /*
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection,acop, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggs_DP, &jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakes_DP, &jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTT_DP,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP, &jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP, &jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection,acop, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggs_DP,&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakes_DP , &jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTT_DP,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP,&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP , &jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP);

	    Ntp->FillHist(t,idx,isOS, GenMatchSelection,acop, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggs_DP, &jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakes_DP, &jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTT_DP,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP, &jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP, &jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggs_DP,&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakes_DP ,&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTT_DP,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP,&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP ,&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection,acop, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggs_DP,&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakes_DP , &jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTT_DP,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP,&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP , &jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection,acop, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggs_DP,&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakes_DP , &jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTT_DP,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP,&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP , &jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection,acop, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggs_DP, &jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakes_DP, &jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTT_DP,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP, &jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP, &jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection,acop, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggs_DP, &jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakes_DP, &jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTT_DP,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP, &jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP, &jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggs_DP, &jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakes_DP, &jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTT_DP,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP, &jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP, &jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggs_DP,&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakes_DP , &jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTT_DP,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP,&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP , &jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection,acop, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggs_DP, &jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakes_DP, &jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTT_DP,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP, &jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP, &jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection,acop, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggs_DP,&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakes_DP ,&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTT_DP,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP,&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP ,&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSHiggs_DP, &jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSJetFakes_DP,&jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSZTT_DP,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSHiggsQCDMC_DP, &jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSJetFakesQCDMC_DP,&jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSZTTQCDMC_DP);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection,acop, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSHiggs_DP,&jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSJetFakes_DP , &jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSZTT_DP,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSHiggsQCDMC_DP,&jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSJetFakesQCDMC_DP , &jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSZTTQCDMC_DP);

	    Ntp->FillHist(t,idx,isOS, GenMatchSelection,acop, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSHiggs_DP, &jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSJetFakes_DP, &jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSZTT_DP,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSHiggsQCDMC_DP, &jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSJetFakesQCDMC_DP, &jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSZTTQCDMC_DP);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection,acop, max_pair,Wspin , (std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin,&jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSHiggs_DP, &jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSJetFakes_DP, &jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSZTT_DP,&emptyvector,&emptyvector,&emptyvector,&jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSHiggsQCDMC_DP, &jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSJetFakesQCDMC_DP, &jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSZTTQCDMC_DP);
	  */
	  Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair, (Wspin/(wIDSF1*wIDSF2))*Ntp->IDSF(Tau1,GenMatch1,"Nom","tau","Up")*Ntp->IDSF(Tau2,GenMatch2,"Nom","tau","Up"),FFData,FFMC,&CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSHiggs_DP,&CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSJetFakes_DP, &CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSZTT_DP,&CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP,&CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_eff_t_pThigh_MVADM10_13TeVUpPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  Ntp->FillHist(t,idx,isOS, GenMatchSelection,acop, max_pair,(Wspin/(wIDSF1*wIDSF2))*Ntp->IDSF(Tau1,GenMatch1,"Nom","tau","Down")*Ntp->IDSF(Tau2,GenMatch2,"Nom","tau","Down") , FFData,FFMC,&CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSHiggs_DP,&CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSJetFakes_DP, &CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSZTT_DP,&CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP,&CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_eff_t_pThigh_MVADM10_13TeVDownPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);

	  Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair, (Wspin/(wTrgSF1*wTrgSF2))*Ntp->TriggerSF(Tau1,GenMatch1,"Nom","Up")*Ntp->TriggerSF(Tau2,GenMatch2,"Nom","Up"),FFData ,FFMC,&CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSHiggs_DP,&CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSJetFakes_DP, &CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSZTT_DP,&CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP,&CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_eff_t_trg_MVADM10_13TeVUpPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  Ntp->FillHist(t,idx,isOS, GenMatchSelection,acop, max_pair,(Wspin/(wTrgSF1*wTrgSF2))*Ntp->TriggerSF(Tau1,GenMatch1,"Nom","Down")*Ntp->TriggerSF(Tau2,GenMatch2,"Nom","Down"),FFData,FFMC,&CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSHiggs_DP,&CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSJetFakes_DP, &CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSZTT_DP,&CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP,&CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_eff_t_trg_MVADM10_13TeVDownPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);

	  Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair ,zptw*zptw*(Wspin/zptw),FFData,FFMC,&CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSHiggs_DP,&CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSJetFakes_DP, &CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSZTT_DP,&CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP,&CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_htt_dyShape_13TeVUpPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair, Wspin/zptw, FFData, FFMC, &CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSHiggs_DP,&CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSJetFakes_DP, &CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSZTT_DP,&CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP,&CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_htt_dyShape_13TeVDownPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);

	  // if(TESUp)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair ,Wspin,FFData,FFMC,&CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSHiggs_DP,&CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSJetFakes_DP, &CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSZTT_DP,&CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP,&CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_scale_t_3prong_13TeVUpPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  // if(TESDown)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair, Wspin, FFData, FFMC, &CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSHiggs_DP,&CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSJetFakes_DP, &CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSZTT_DP,&CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP,&CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_scale_t_3prong_13TeVDownPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);

	  // if(JERUp)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair ,Wspin,FFData,FFMC,&CMS_res_j_13TeVUpPVRefitWithTracksBSHiggs_DP,&CMS_res_j_13TeVUpPVRefitWithTracksBSJetFakes_DP, &CMS_res_j_13TeVUpPVRefitWithTracksBSZTT_DP,&CMS_res_j_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP,&CMS_res_j_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_res_j_13TeVUpPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  // if(JERDown)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair, Wspin, FFData, FFMC, &CMS_res_j_13TeVDownPVRefitWithTracksBSHiggs_DP,&CMS_res_j_13TeVDownPVRefitWithTracksBSJetFakes_DP, &CMS_res_j_13TeVDownPVRefitWithTracksBSZTT_DP,&CMS_res_j_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP,&CMS_res_j_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_res_j_13TeVDownPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  

	  
	  // if(jetAbsoluteUp)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair ,Wspin,FFData,FFMC,&CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSHiggs_DP,&CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSJetFakes_DP, &CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSZTT_DP,&CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP,&CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_scale_j_Absolute_13TeVUpPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  // if(jetBBEC1Up)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair ,Wspin,FFData,FFMC,&CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSHiggs_DP,&CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSJetFakes_DP, &CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSZTT_DP,&CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP,&CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_scale_j_BBEC1_13TeVUpPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  // if(jetEC2Up)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair ,Wspin,FFData,FFMC,&CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSHiggs_DP,&CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSJetFakes_DP, &CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSZTT_DP,&CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP,&CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_scale_j_EC2_13TeVUpPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  // if(jetFlavorQCDUp)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair ,Wspin,FFData,FFMC,&CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSHiggs_DP,&CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSJetFakes_DP, &CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSZTT_DP,&CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP,&CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_scale_j_FlavorQCD_13TeVUpPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  // if(jetHFUp)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair ,Wspin,FFData,FFMC,&CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSHiggs_DP,&CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSJetFakes_DP, &CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSZTT_DP,&CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP,&CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_scale_j_HF_13TeVUpPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  // if(jetRelativeBalUp)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair ,Wspin,FFData,FFMC,&CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSHiggs_DP,&CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSJetFakes_DP, &CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSZTT_DP,&CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP,&CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_scale_j_RelativeBal_13TeVUpPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  // if(jetAbsoluteYearUp)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair ,Wspin,FFData,FFMC,&CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSHiggs_DP,&CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSJetFakes_DP, &CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSZTT_DP,&CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP,&CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_scale_j_Absolute_Year_13TeVUpPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  // if(jetBBEC1YearUp)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair ,Wspin,FFData,FFMC,&CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSHiggs_DP,&CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSJetFakes_DP, &CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSZTT_DP,&CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP,&CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_scale_j_BBEC1_Year_13TeVUpPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  // if(jetEC2YearUp)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair ,Wspin,FFData,FFMC,&CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSHiggs_DP,&CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSJetFakes_DP, &CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSZTT_DP,&CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP,&CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_scale_j_EC2_Year_13TeVUpPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  // if(jetHFYearUp)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair ,Wspin,FFData,FFMC,&CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSHiggs_DP,&CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSJetFakes_DP, &CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSZTT_DP,&CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP,&CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_scale_j_HF_Year_13TeVUpPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  // if(jetRelativeSampleYearUp)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair ,Wspin,FFData,FFMC,&CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSHiggs_DP,&CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSJetFakes_DP, &CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSZTT_DP,&CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP,&CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_scale_j_RelativeSample_Year_13TeVUpPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	

	  // if(jetAbsoluteDown)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair ,Wspin,FFData,FFMC,&CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSHiggs_DP,&CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSJetFakes_DP, &CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSZTT_DP,&CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP,&CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_scale_j_Absolute_13TeVDownPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  // if(jetBBEC1Down)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair ,Wspin,FFData,FFMC,&CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSHiggs_DP,&CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSJetFakes_DP, &CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSZTT_DP,&CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP,&CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_scale_j_BBEC1_13TeVDownPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  // if(jetEC2Down)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair ,Wspin,FFData,FFMC,&CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSHiggs_DP,&CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSJetFakes_DP, &CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSZTT_DP,&CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP,&CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_scale_j_EC2_13TeVDownPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  // if(jetFlavorQCDDown)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair ,Wspin,FFData,FFMC,&CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSHiggs_DP,&CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSJetFakes_DP, &CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSZTT_DP,&CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP,&CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_scale_j_FlavorQCD_13TeVDownPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  // if(jetHFDown)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair ,Wspin,FFData,FFMC,&CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSHiggs_DP,&CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSJetFakes_DP, &CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSZTT_DP,&CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP,&CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_scale_j_HF_13TeVDownPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  // if(jetRelativeBalDown)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair ,Wspin,FFData,FFMC,&CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSHiggs_DP,&CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSJetFakes_DP, &CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSZTT_DP,&CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP,&CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_scale_j_RelativeBal_13TeVDownPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  // if(jetAbsoluteYearDown)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair ,Wspin,FFData,FFMC,&CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSHiggs_DP,&CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSJetFakes_DP, &CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSZTT_DP,&CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP,&CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_scale_j_Absolute_Year_13TeVDownPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  // if(jetBBEC1YearDown)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair ,Wspin,FFData,FFMC,&CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSHiggs_DP,&CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSJetFakes_DP, &CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSZTT_DP,&CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP,&CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_scale_j_BBEC1_Year_13TeVDownPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  // if(jetEC2YearDown)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair ,Wspin,FFData,FFMC,&CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSHiggs_DP,&CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSJetFakes_DP, &CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSZTT_DP,&CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP,&CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_scale_j_EC2_Year_13TeVDownPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  // if(jetHFYearDown)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair ,Wspin,FFData,FFMC,&CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSHiggs_DP,&CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSJetFakes_DP, &CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSZTT_DP,&CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP,&CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_scale_j_HF_Year_13TeVDownPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  // if(jetRelativeSampleYearDown)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair ,Wspin,FFData,FFMC,&CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSHiggs_DP,&CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSJetFakes_DP, &CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSZTT_DP,&CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP,&CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_scale_j_RelativeSample_Year_13TeVDownPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);


	  // if(METScaleUp)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair ,Wspin,FFData,FFMC,&CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSHiggs_DP,&CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSJetFakes_DP, &CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSZTT_DP,&CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP,&CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_htt_boson_scale_met_13TeVUpPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  // if(METScaleDown)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair, Wspin, FFData, FFMC, &CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSHiggs_DP,&CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSJetFakes_DP, &CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSZTT_DP,&CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP,&CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_htt_boson_scale_met_13TeVDownPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  // if(METResoUp)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair ,Wspin,FFData,FFMC,&CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSHiggs_DP,&CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSJetFakes_DP, &CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSZTT_DP,&CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP,&CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_htt_boson_reso_met_13TeVUpPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  // if(METResoDown)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair, Wspin, FFData, FFMC, &CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSHiggs_DP,&CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSJetFakes_DP, &CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSZTT_DP,&CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP,&CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_htt_boson_reso_met_13TeVDownPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  // if(METUnclusteredScaleUp)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair ,Wspin,FFData,FFMC,&CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSHiggs_DP,&CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSJetFakes_DP, &CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSZTT_DP,&CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP,&CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_scale_met_unclustered_13TeVUpPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  // if(METUnclusteredScaleDown)Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair, Wspin, FFData, FFMC, &CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSHiggs_DP,&CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSJetFakes_DP, &CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSZTT_DP,&CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP,&CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_scale_met_unclustered_13TeVDownPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  
	  if(id==11||id==12||id==45||id==460||id==461){
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair ,Wspin*Ntp->TheoreticalPSUnc(6)/Ntp->MC_weight(),FFData,FFMC,&CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSHiggs_DP,&CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSJetFakes_DP, &CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSZTT_DP,&CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP,&CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_PS_ISR_ggH_13TeVUpPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair, Wspin*Ntp->TheoreticalPSUnc(8)/Ntp->MC_weight(), FFData, FFMC, &CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSHiggs_DP,&CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSJetFakes_DP, &CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSZTT_DP,&CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP,&CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_PS_ISR_ggH_13TeVDownPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair ,Wspin*Ntp->TheoreticalPSUnc(7)/Ntp->MC_weight(),FFData,FFMC,&CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSHiggs_DP,&CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSJetFakes_DP, &CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSZTT_DP,&CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP,&CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_PS_FSR_ggH_13TeVUpPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair, Wspin*Ntp->TheoreticalPSUnc(9)/Ntp->MC_weight(), FFData, FFMC, &CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSHiggs_DP,&CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSJetFakes_DP, &CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSZTT_DP,&CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP,&CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_PS_FSR_ggH_13TeVDownPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);

	    Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair ,Wspin*Ntp->TheoreticalScaleUnc1005()/Ntp->nominal_wt(),FFData,FFMC,&CMS_scale_gg_13TeVUpPVRefitWithTracksBSHiggs_DP,&CMS_scale_gg_13TeVUpPVRefitWithTracksBSJetFakes_DP, &CMS_scale_gg_13TeVUpPVRefitWithTracksBSZTT_DP,&CMS_scale_gg_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP,&CMS_scale_gg_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_scale_gg_13TeVUpPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair, Wspin*Ntp->TheoreticalScaleUnc1009()/Ntp->nominal_wt(), FFData, FFMC, &CMS_scale_gg_13TeVDownPVRefitWithTracksBSHiggs_DP,&CMS_scale_gg_13TeVDownPVRefitWithTracksBSJetFakes_DP, &CMS_scale_gg_13TeVDownPVRefitWithTracksBSZTT_DP,&CMS_scale_gg_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP,&CMS_scale_gg_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP, &CMS_scale_gg_13TeVDownPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);
	  }
	  
	  Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair,top_wt*top_wt*(Wspin/top_wt) ,FFData,FFMC,&CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSHiggs_DP,&CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSJetFakes_DP,&CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSZTT_DP,&CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP,&CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP,&CMS_htt_ttbarShape_13TeVUpPVRefitWithTracksBSWfakesZTT ,&emptyvector,&emptyvector,&emptyvector);
	  Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair,Wspin/top_wt ,FFData,FFMC,&CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSHiggs_DP,&CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSJetFakes_DP,&CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSZTT_DP,&CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP,&CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP,&CMS_htt_ttbarShape_13TeVDownPVRefitWithTracksBSWfakesZTT ,&emptyvector,&emptyvector,&emptyvector);
	  
	  if(id==71||id==72||id==73||id==74 || id==47||id==48||id==49||id==50||id==51||id==52||id==53||id==54||id==55||id==56||id==57||id==58 || id==70 ||id==701 ||id==702 ||id==703){ // VV+TT
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair,Wspin ,FFData,FFMC,&CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSHiggs_DP,&CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSJetFakes_DP,&CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSZTT_DP,&CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesHiggs_DP,&CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesJetFakes_DP,&CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesZTT_DP ,&emptyvector,&emptyvector,&emptyvector);
	    Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair,Wspin ,FFData,FFMC,&CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSHiggs_DP,&CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSJetFakes_DP,&CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSZTT_DP,&CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesHiggs_DP,&CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesJetFakes_DP,&CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesZTT_DP ,&emptyvector,&emptyvector,&emptyvector);
	    GenMatchSelection=!(GenMatch1==6 || GenMatch2==6);
	    if (GenMatch1==5 && GenMatch2==5){
	      Ntp->FillHist(t,idx,isOS, GenMatchSelection, angle, max_pair,Wspin,FFData,FFMC,&ttbarcontaminationHiggs_DP,&ttbarcontaminationJetFakes_DP,&ttbarcontaminationZTT_DP,&ttbarcontaminationWfakesHiggs_DP,&ttbarcontaminationWfakesJetFakes_DP,&ttbarcontaminationWfakesZTT_DP ,&emptyvector,&emptyvector,&emptyvector);
	    }
	    GenMatchSelection=(!(GenMatch1==6 || GenMatch2==6) &&  !(GenMatch1==5 && GenMatch2==5));
	  }

			      
	  Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair, prefup*(Wspin/Ntp->prefiringweight()),FFData, FFMC,&PrefiringUpPVRefitWithTracksBSHiggs_DP,&PrefiringUpPVRefitWithTracksBSJetFakes_DP, &PrefiringUpPVRefitWithTracksBSZTT_DP,&PrefiringUpPVRefitWithTracksBSWfakesHiggs_DP,&PrefiringUpPVRefitWithTracksBSWfakesJetFakes_DP, &PrefiringUpPVRefitWithTracksBSWfakesZTT ,&emptyvector,&emptyvector,&emptyvector);
	  Ntp->FillHist(t,idx,isOS, GenMatchSelection, acop, max_pair, prefdown*(Wspin/Ntp->prefiringweight()),FFData,FFMC, &PrefiringDownPVRefitWithTracksBSHiggs_DP,&PrefiringDownPVRefitWithTracksBSJetFakes_DP, &PrefiringDownPVRefitWithTracksBSZTT_DP,&PrefiringDownPVRefitWithTracksBSWfakesHiggs_DP,&PrefiringDownPVRefitWithTracksBSWfakesJetFakes_DP, &PrefiringDownPVRefitWithTracksBSWfakesZTT_DP,&emptyvector,&emptyvector,&emptyvector);

	}
      }
    if(Ntp->isData() && !isEmbed && !Ntp->byMediumDeepTau2017v2p1VSjet_1(IndexSelectedTemp.at(Sorted.back())) && Ntp->byVVVLooseDeepTau2017v2p1VSjet_1(IndexSelectedTemp.at(Sorted.back())) && Ntp->byMediumDeepTau2017v2p1VSjet_2(IndexSelectedTemp.at(Sorted.back()))) {

      double isOS=(Ntp->Daughters_charge(Tau1)/abs(Ntp->Daughters_charge(Tau1))) != (Ntp->Daughters_charge(Tau2)/abs(Ntp->Daughters_charge(Tau2)));
      auto args = std::vector<double>{Tau1P4.Pt(),(double)Ntp->MVADM2017(Tau1),(double)Ntp->njets(IndexSelectedTemp.at(Sorted.back())),isOS,PUPPIMET*cos(TLorentzVector(Tau1P4.X(),Tau1P4.Y(),Tau1P4.Z(),0.).DeltaPhi(TLorentzVector(PUPPImetCorr_px,PUPPImetCorr_py,0.,0.)))/Tau1P4.Pt(),Tau1P4.DeltaR(Tau2P4)};
      double FFData;
      if(Ntp->year()==2016)FFData= std::shared_ptr<RooFunctor>(wFF2016->function("ff_tt_medium_mvadmbins")->functor(wFF2016->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data());
      if(Ntp->year()==2017)FFData= std::shared_ptr<RooFunctor>(wFF2017->function("ff_tt_medium_mvadmbins")->functor(wFF2017->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data());
      if(Ntp->year()==2018)FFData= std::shared_ptr<RooFunctor>(wFF2018->function("ff_tt_medium_mvadmbins")->functor(wFF2018->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data());

      if (std::isnan(Wspin)!=true)
	{
	  if(a1a1MVA)
	    {
	      
	      Tau1PTa1a1.at(1).Fill(pt_1_,FFData*w);
	      TauTauVisMassa1a1.at(1).Fill(m_vis_,FFData*w);
	      TauTauFullMassa1a1.at(1).Fill(m_sv_,FFData*w);
	      TauTauVisPTa1a1.at(1).Fill(pt_vis_,FFData*w);
	      TauTauFullPTa1a1.at(1).Fill(pt_tt_,FFData*w);
	      Mjja1a1.at(1).Fill(mjj_,FFData*w);
	      jdetaa1a1.at(1).Fill(jdeta_,FFData*w);
	      jpt_1a1a1.at(1).Fill(jpt_1_,FFData*w);
	      PUPPImetcorra1a1.at(1).Fill(met_,FFData*w);
	      NbJetsa1a1.at(1).Fill(n_jets_,FFData*w);
	      
	      if ((HadRefitPions_minus!=HadRefitPions_plus) && (HadRefitPions_minus!=VectZeroLV) && (HadRefitPions_plus!=VectZeroLV))
		{
		  if(TauminusPairConstraintWithTracksBS!=TauplusPairConstraintWithTracksBS && TauminusPairConstraintWithTracksBS!=zeroLV && TauplusPairConstraintWithTracksBS!=zeroLV )
		    {
		      if(ScalcPVRefitWithTracksBS.isOk("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus)==true){polarimetricAcopAnglePVRefitWithTracksBSMVADM.at(1).Fill(ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),FFData);
			
			polarimetricAcopAnglePVRefitWithTracksBSMVADM_DP.at(1).Fill(acop,FFData);
		      }
		      if(max_pair.second==0)
			{

			  if(ScalcPVRefitWithTracksBS.isOk("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus)==true){
			    
			    HiggsBDTScorea1a1.at(1).Fill(max_pair.first,FFData);
			    HiggsBDTScorea1a1_DP.at(1).Fill(max_pair.first,FFData);//AC?

			    if(max_pair.first>0 && max_pair.first<0.7){polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled.at(1).Fill(ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),FFData*Wspin);
			      polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled_DP.at(1).Fill(acop,FFData*Wspin);
			    }
			    if(max_pair.first>0.7 && max_pair.first<0.8){polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled.at(1).Fill(2*TMath::Pi()+ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),FFData*Wspin);
			      polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled_DP.at(1).Fill(2*TMath::Pi()+acop,FFData*Wspin);
			    }
			    if(max_pair.first>0.8 && max_pair.first<1.0){polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled.at(1).Fill(2*2*TMath::Pi()+ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),FFData*Wspin);
			      polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled_DP.at(1).Fill(2*2*TMath::Pi()+acop,FFData*Wspin);
			    }

			    if(Ntp->njets(IndexSelectedTemp.at(Sorted.back()))==0){
			      jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

			      jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			    }
			    else{
			      jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

			      jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			    }
			    if(Ntp->njets(IndexSelectedTemp.at(Sorted.back()))==1){
			      jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

			      jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			    }
			    else{
			      jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

			      jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			    }
			    if(Ntp->njets(IndexSelectedTemp.at(Sorted.back()))==2){
			      jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

			      jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			    }
			    else{
			      jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSHiggs.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

			      jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSHiggs_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			    }


                           
                            // jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSHiggs.at(1).Fill(ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			    // jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSHiggs.at(1).Fill(ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
                            // jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSHiggs.at(1).Fill(ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
                            // jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSHiggs.at(1).Fill(ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

			  }
			}
		      if(max_pair.second==1)
			{
			  
			  if(ScalcPVRefitWithTracksBS.isOk("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus)==true){
			    JetFakesBDTScorea1a1.at(1).Fill(max_pair.first,FFData);
			    JetFakesBDTScorea1a1_DP.at(1).Fill(max_pair.first,FFData);//AC?  

			    if(max_pair.first>0 && max_pair.first<0.7){polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled.at(1).Fill(ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),FFData);
			      polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled_DP.at(1).Fill(acop,FFData);
			    }
			    if(max_pair.first>0.7 && max_pair.first<0.8){polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled.at(1).Fill(2*TMath::Pi()+ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),FFData);
			      polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled_DP.at(1).Fill(2*TMath::Pi()+acop,FFData);
			    }
			    if(max_pair.first>0.8 && max_pair.first<1.0){polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled.at(1).Fill(2*2*TMath::Pi()+ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),FFData);
			      polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled_DP.at(1).Fill(2*2*TMath::Pi()+acop,FFData);
			    }

			    if(Ntp->njets(IndexSelectedTemp.at(Sorted.back()))==0){
			      jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

			      jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			    }
			    else{
			      jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

			      jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			    }
			    if(Ntp->njets(IndexSelectedTemp.at(Sorted.back()))==1){
			      jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

			      jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			    }
			    else{
			      jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

			      jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			    }
			    if(Ntp->njets(IndexSelectedTemp.at(Sorted.back()))==2){
			      jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

			      jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			    }
			    else{
			      jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSJetFakes.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

			      jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSJetFakes_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

			    }
			  }
			}
		      if(max_pair.second==2)
			{

			  if(ScalcPVRefitWithTracksBS.isOk("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus)==true){
			    ZTTBDTScorea1a1.at(1).Fill(max_pair.first,FFData);
			    ZTTBDTScorea1a1_DP.at(1).Fill(max_pair.first,FFData);

			    if(max_pair.first>0 && max_pair.first<0.7){polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled.at(1).Fill(ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),FFData);
			      polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled_DP.at(1).Fill(acop,FFData);
			    }
			    if(max_pair.first>0.7 && max_pair.first<0.8){polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled.at(1).Fill(2*TMath::Pi()+ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),FFData);
			      polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled_DP.at(1).Fill(2*TMath::Pi()+acop,FFData);
			    }
			    if(max_pair.first>0.8 && max_pair.first<1.0){polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled.at(1).Fill(2*2*TMath::Pi()+ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),FFData);
			    
			      polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled_DP.at(1).Fill(2*2*TMath::Pi()+acop,FFData);
			    }

			    if(Ntp->njets(IndexSelectedTemp.at(Sorted.back()))==0){

			      jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

			      jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

			      jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

			      jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

			      jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

			      jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

			      jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

			    }
			    else{

			      jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

			      jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

			    }
			    if(Ntp->njets(IndexSelectedTemp.at(Sorted.back()))==1){

			      jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

			      jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			    }
			    else{
			      jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

			      jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			    }
			    if(Ntp->njets(IndexSelectedTemp.at(Sorted.back()))==2){

			      jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

			      jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			    }
			    else{
			      jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSZTT.at(1).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

			      jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
			      jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
			      jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSZTT_DP.at(1).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

			    }
			    		    
			  }
			}
		    }
		}
	    }
	}
    }

    if((!Ntp->isData() /*|| isEmbed*/) && !Ntp->byMediumDeepTau2017v2p1VSjet_1(IndexSelectedTemp.at(Sorted.back())) && Ntp->byVVVLooseDeepTau2017v2p1VSjet_1(IndexSelectedTemp.at(Sorted.back()))  && Ntp->byMediumDeepTau2017v2p1VSjet_2(IndexSelectedTemp.at(Sorted.back()))) {

      if( GenMatch1!=6) {

	double isOS=(Ntp->Daughters_charge(Tau1)/abs(Ntp->Daughters_charge(Tau1))) != (Ntp->Daughters_charge(Tau2)/abs(Ntp->Daughters_charge(Tau2)));
	auto args = std::vector<double>{Tau1P4.Pt(),(double)Ntp->MVADM2017(Tau1),(double)Ntp->njets(IndexSelectedTemp.at(Sorted.back())),isOS,PUPPIMET*cos(TLorentzVector(Tau1P4.X(),Tau1P4.Y(),0.,0.).DeltaPhi(TLorentzVector(PUPPImetCorr_px,PUPPImetCorr_py,0.,0.)))/Tau1P4.Pt(),Tau1P4.DeltaR(Tau2P4)};
	double FFMC;
	if(Ntp->year()==2016)FFMC= std::shared_ptr<RooFunctor>(wFF2016->function("ff_tt_medium_mvadmbins")->functor(wFF2016->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data());
	if(Ntp->year()==2017)FFMC= std::shared_ptr<RooFunctor>(wFF2017->function("ff_tt_medium_mvadmbins")->functor(wFF2017->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data());
	if(Ntp->year()==2018)FFMC= std::shared_ptr<RooFunctor>(wFF2018->function("ff_tt_medium_mvadmbins")->functor(wFF2018->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data());
	
	if (std::isnan(Wspin)!=true)
	  {	  
	    if(a1a1MVA)
	      {	  
		Tau1PTa1a1QCDMC.at(t).Fill(pt_1_,FFMC*w);
		TauTauVisMassa1a1QCDMC.at(t).Fill(m_vis_,FFMC*w);
		TauTauFullMassa1a1QCDMC.at(t).Fill(m_sv_,FFMC*w);
		TauTauVisPTa1a1QCDMC.at(t).Fill(pt_vis_,FFMC*w);
		TauTauFullPTa1a1QCDMC.at(t).Fill(pt_tt_,FFMC*w);
		Mjja1a1QCDMC.at(t).Fill(mjj_,FFMC*w);
		jdetaa1a1QCDMC.at(t).Fill(jdeta_,FFMC*w);
		jpt_1a1a1QCDMC.at(t).Fill(jpt_1_,FFMC*w);
		PUPPImetcorra1a1QCDMC.at(t).Fill(met_,FFMC*w);
		NbJetsa1a1QCDMC.at(t).Fill(n_jets_,FFMC*w);
		if ((HadRefitPions_minus!=HadRefitPions_plus) && (HadRefitPions_minus!=VectZeroLV) && (HadRefitPions_plus!=VectZeroLV))
		  {
		    if(TauminusPairConstraintWithTracksBS!=TauplusPairConstraintWithTracksBS && TauminusPairConstraintWithTracksBS!=zeroLV && TauplusPairConstraintWithTracksBS!=zeroLV )
		      {
			
			if(ScalcPVRefitWithTracksBS.isOk("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus)==true){polarimetricAcopAnglePVRefitWithTracksBSMVADMQCDMC.at(t).Fill(ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),FFMC*w);

			  polarimetricAcopAnglePVRefitWithTracksBSMVADMQCDMC_DP.at(t).Fill(acop,FFMC*w);
			}
			if(max_pair.second==0)
			  {
			    
			    
			    if(ScalcPVRefitWithTracksBS.isOk("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus)==true){
			      HiggsBDTScorea1a1QCDMC.at(t).Fill(max_pair.first,FFMC*w);
			      HiggsBDTScorea1a1QCDMC_DP.at(t).Fill(max_pair.first,FFMC*w);//AC?
			      
			      if(max_pair.first>0 && max_pair.first<0.7){polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolledQCDMC.at(t).Fill(ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),FFMC*Wspin);
				polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolledQCDMC_DP.at(t).Fill(acop,FFMC*Wspin);
			      }
			      if(max_pair.first>0.7 && max_pair.first<0.8){polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolledQCDMC.at(t).Fill(2*TMath::Pi()+ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),FFMC*Wspin);
				polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolledQCDMC_DP.at(t).Fill(2*TMath::Pi()+acop,FFMC*Wspin);
			      }
			      if(max_pair.first>0.8 && max_pair.first<1.){polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolledQCDMC.at(t).Fill(2*2*TMath::Pi()+ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),FFMC*Wspin);
				polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolledQCDMC_DP.at(t).Fill(2*2*TMath::Pi()+acop,FFMC*Wspin);
			      }  

			      if(Ntp->njets(IndexSelectedTemp.at(Sorted.back()))==0){
				jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

				jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      }
			      else{
				jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

				jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      }
			      if(Ntp->njets(IndexSelectedTemp.at(Sorted.back()))==1){
				jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

				jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      }
			      else{
				jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

				jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      }
			      if(Ntp->njets(IndexSelectedTemp.at(Sorted.back()))==2){
				jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

				jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      }
			      else{
				jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSHiggsQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

				jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSHiggsQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      }
			    }
			    
			  }
			if(max_pair.second==1)
			  {
			    

			    if(ScalcPVRefitWithTracksBS.isOk("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus)==true){ 
			      JetFakesBDTScorea1a1QCDMC.at(t).Fill(max_pair.first,FFMC*w);
			      JetFakesBDTScorea1a1QCDMC_DP.at(t).Fill(max_pair.first,FFMC*w);//AC?

			      if(max_pair.first>0 && max_pair.first<0.7){polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolledQCDMC.at(t).Fill(ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),FFMC*w);
				polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolledQCDMC_DP.at(t).Fill(acop,FFMC*w);  
			      }
			      if(max_pair.first>0.7 && max_pair.first<0.8){polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolledQCDMC.at(t).Fill(2*TMath::Pi()+ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),FFMC*w);
				polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolledQCDMC_DP.at(t).Fill(2*TMath::Pi()+acop,FFMC*w);
			      }
			      if(max_pair.first>0.8 && max_pair.first<1.){polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolledQCDMC.at(t).Fill(2*2*TMath::Pi()+ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),FFMC*w);
				polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolledQCDMC_DP.at(t).Fill(2*2*TMath::Pi()+acop,FFMC*w);   
			      } 
			      if(Ntp->njets(IndexSelectedTemp.at(Sorted.back()))==0){
				jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

				jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      }
			      else{
				jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

				jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      }
			      if(Ntp->njets(IndexSelectedTemp.at(Sorted.back()))==1){
				jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

				jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      }
			      else{
				jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

				jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      }
			      if(Ntp->njets(IndexSelectedTemp.at(Sorted.back()))==2){
				jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

				jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      }
			      else{
				jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSJetFakesQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

				jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      }
			    }
			  }
			if(max_pair.second==2)
			  {
			    if(ScalcPVRefitWithTracksBS.isOk("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus)==true){
			      
			      ZTTBDTScorea1a1QCDMC.at(t).Fill(max_pair.first,FFMC*w);
			      ZTTBDTScorea1a1QCDMC_DP.at(t).Fill(max_pair.first,FFMC*w);//AC?

			      if(max_pair.first>0 && max_pair.first<0.7){polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolledQCDMC.at(t).Fill(ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),FFMC*w);
				polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolledQCDMC_DP.at(t).Fill(acop,FFMC*w);
			      }
			      if(max_pair.first>0.7 && max_pair.first<0.8){polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolledQCDMC.at(t).Fill(2*TMath::Pi()+ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),FFMC*w);
				polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolledQCDMC_DP.at(t).Fill(2*TMath::Pi()+acop,FFMC*w);
			      }
			      if(max_pair.first>0.8 && max_pair.first<1.){polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolledQCDMC.at(t).Fill(2*2*TMath::Pi()+ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),FFMC*w);
				polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolledQCDMC_DP.at(t).Fill(2*2*TMath::Pi()+acop,FFMC*w);
			      }

			      if(Ntp->njets(IndexSelectedTemp.at(Sorted.back()))==0){
				jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

				jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet0_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      }
			      else{
				jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

				jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      }
			      if(Ntp->njets(IndexSelectedTemp.at(Sorted.back()))==1){
				jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

				jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet1_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      }
			      else{
				jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

				jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      }
			      if(Ntp->njets(IndexSelectedTemp.at(Sorted.back()))==2){
				jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

				jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc1_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_stat_unc2_njet2_mvadm10_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_met_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_up")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins_qcd_syst_down")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      }
			      else{
				jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSZTTQCDMC.at(t).Fill(angle,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);

				jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);			    
				jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      
				jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
				jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSZTTQCDMC_DP.at(t).Fill(acop,max_pair.first,(std::shared_ptr<RooFunctor>(wFF->function("ff_tt_medium_mvadmbins")->functor(wFF->argSet("pt,mvadm,njets,os,met_var_qcd,dR")))->eval(args.data()))*Wspin);
			      }
			    }
			  }
		      }
		  }
	      
	      }
	  }
	
      }
    }

    if((!Ntp->isData()/*|| isEmbed*/) && pass.at(TausIsolation)){
      if(GenMatch1!=6 && GenMatch2==6){
	if (std::isnan(Wspin)!=true)
	  {

	    if(a1a1MVA)
	      {
		if ((HadPions_minus!=HadPions_plus) && (HadPions_minus!=VectZeroLV) && (HadPions_plus!=VectZeroLV) && (HadRefitPions_minus!=HadRefitPions_plus) && (HadRefitPions_minus!=VectZeroLV) && (HadRefitPions_plus!=VectZeroLV))
		  {
		    if(TauminusPairConstraintWithTracksBS!=TauplusPairConstraintWithTracksBS && TauminusPairConstraintWithTracksBS!=zeroLV && TauplusPairConstraintWithTracksBS!=zeroLV )
		      {
			if(ScalcPVRefitWithTracksBS.isOk("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus)==true){polarimetricAcopAnglePVRefitWithTracksBSMVADM.at(t).Fill(ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);
			  polarimetricAcopAnglePVRefitWithTracksBSMVADM_DP.at(t).Fill(acop,Wspin);
			}
			if(max_pair.second==0)
			  {

			    if(ScalcPVRefitWithTracksBS.isOk("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus)==true){
			      WfakesHiggs.at(t).Fill(ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),max_pair.first,Wspin);
			      
			      WfakesHiggs_DP.at(t).Fill(acop,max_pair.first,Wspin); 
			    }
			  }
			if(max_pair.second==1)
			  {

			    if(ScalcPVRefitWithTracksBS.isOk("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus)==true){
			      WfakesJetFakes.at(t).Fill(ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),max_pair.first,Wspin);
			      
			      WfakesJetFakes_DP.at(t).Fill(acop,max_pair.first,Wspin);
			    }
			  }
			if(max_pair.second==2)
			  {

			    if(ScalcPVRefitWithTracksBS.isOk("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus)==true){

			      WfakesZTT.at(t).Fill(ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),max_pair.first,Wspin);

			      WfakesZTT_DP.at(t).Fill(acop,max_pair.first,Wspin);
			    }
			  }
		      }
		  }
	      }
	  }
      }
    }
  }

  bool status=AnalysisCuts(t,w,wobs);  // boolean that say whether your event passed critera defined in pass vector. The whole vector must be true for status = true
  ///////////////////////////////////////////////////////////
  // Analyse events which passed selection
  if(status && GenMatchSelection) {

    //cout<<"genMatch1: "<<GenMatch1<<"  "<<"genMatch2: "<<GenMatch2<<endl;
    //cout<<w<<endl;
    double pvx(0);
    pvx =  Ntp->npv();
    // if(id == DataMCType::Data) pvx =  Ntp->npv();
    if(id !=DataMCType::Data && id !=DataMCType::QCD)	  pvx = Ntp->PUNumInteractions();
    //NPrimeVtx.at(t).Fill(pvx,w);
    //NPU.at(t).Fill(Ntp->npu(),w);
    //RHO.at(t).Fill(Ntp->rho(),w);

    std::vector<int> thirdLepton;

    //    TLorentzVector taunew(tau1P4.Px(), );

    //svfTau1E.at(t).Fill(tau1P4.E(),w);
    //svfTau2E.at(t).Fill(tau2P4.E(),w);

    // Tau1PT.at(t).Fill(Tau1P4.Pt(),w);
    // Tau1E.at(t).Fill(Tau1P4.E(),w);
    // Tau1Mass.at(t).Fill(Tau1P4.M(),w);
    // Tau1Phi.at(t).Fill(Tau1P4.Phi(),w);
    // Tau1Eta.at(t).Fill(Tau1P4.Eta(),w);
    // Tau1dz.at(t).Fill(Ntp->dz(Tau1),w);
    // Tau1HPSDecayMode.at(t).Fill(Ntp->decayMode(Tau1),w);
    // Tau1MVADecayMode.at(t).Fill(Ntp->MVADM2017(Tau1),w);
    // if(!Ntp->isData())Tau1GenMatch.at(t).Fill(GenMatch1,w);
    // //    cout<<"genmatch:"<<GenMatch1<<endl;
    // Tau2PT.at(t).Fill(Tau2P4.Pt(),w);
    // Tau2E.at(t).Fill(Tau2P4.E(),w);
    // Tau2Mass.at(t).Fill(Tau2P4.M(),w);
    // Tau2Phi.at(t).Fill(Tau2P4.Phi(),w);
    // Tau2Eta.at(t).Fill(Tau2P4.Eta(),w);
    // Tau2dz.at(t).Fill(Ntp->dz(Tau2),w);
    // Tau2HPSDecayMode.at(t).Fill(Ntp->decayMode(Tau2),w);
    // Tau2MVADecayMode.at(t).Fill(Ntp->MVADM2017(Tau2),w);
    //    if(!Ntp->isData())Tau2GenMatch.at(t).Fill(GenMatch2,w);
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
    //if(thirdLepton.size()>0)ExtraLeptonVeto.at(t).Fill(1.,w);
    //else ExtraLeptonVeto.at(t).Fill(0.,w);

    // TauTauVisMass.at(t).Fill((Tau1P4+Tau2P4).M(),w);
    // TauTauFullMass.at(t).Fill((tau1P4+tau2P4).M(),w);
    // TauTauFullPT.at(t).Fill((tau1P4+tau2P4).Pt(),w);
    // dRTauTau.at(t).Fill(Tau1P4.DeltaR(Tau2P4),w);
    
    // MET.at(t).Fill(Ntp->MET(),w);
    // METphi.at(t).Fill(Ntp->METphi(),w);
    // PUPPImet.at(t).Fill(Ntp->PUPPImet(),w);
    // PUPPImetphi.at(t).Fill(Ntp->PUPPImetphi(),w);
    
    // PUPPImetcorr.at(t).Fill(sqrt(PUPPImetCorr_px*PUPPImetCorr_px+PUPPImetCorr_py*PUPPImetCorr_py),w);
    // PUPPImetcorrphi.at(t).Fill(TMath::ATan2(PUPPImetCorr_py,PUPPImetCorr_px),w);
    
    // TransverseMass.at(t).Fill(Ntp->transverseMass(Tau1P4.Pt(), Tau1P4.Phi(), Tau2P4.Pt(), Tau2P4.Phi()),w);
    // NbJets.at(t).Fill(Ntp->njets(IndexSelectedTemp.at(Sorted.back())),w);

    // if(Has2jets>=0)
    //   {
    // 	Mjj.at(t).Fill(Ntp->mjj(Has2jets),w);
    // 	dijetpt.at(t).Fill(Ntp->dijetpt(Has2jets),w);
    // 	dijetphi.at(t).Fill(Ntp->dijetphi(Has2jets),w);
    // 	jdeta.at(t).Fill(Ntp->jdeta(Has2jets),w);
    // 	jdphi.at(t).Fill(Ntp->jdphi(Has2jets),w);
    // 	jpt_2.at(t).Fill(Ntp->jpt_2(Has2jets),w);
    // 	jeta_2.at(t).Fill(Ntp->jeta_2(Has2jets),w);
    // 	jphi_2.at(t).Fill(Ntp->jphi_2(Has2jets),w);
	  
    //   }

    // if(Has1jet>=0)
    //   {
    // 	jpt_1.at(t).Fill(Ntp->jpt_1(Has1jet),w);
    // 	jeta_1.at(t).Fill(Ntp->jeta_1(Has1jet),w);
    // 	jphi_1.at(t).Fill(Ntp->jphi_1(Has1jet),w);
    //   }


    if (std::isnan(Wspin)!=true)
      {
	if(a1a1MVA)
	  {
	   
	    Tau1PTa1a1.at(t).Fill(pt_1_,w);
	    TauTauVisMassa1a1.at(t).Fill(m_vis_,w);
	    TauTauFullMassa1a1.at(t).Fill(m_sv_,w);
	    TauTauVisPTa1a1.at(t).Fill(pt_vis_,w);
	    TauTauFullPTa1a1.at(t).Fill(pt_tt_,w);
	    Mjja1a1.at(t).Fill(mjj_,w);
	    jdetaa1a1.at(t).Fill(jdeta_,w);
	    jpt_1a1a1.at(t).Fill(jpt_1_,w);
	    PUPPImetcorra1a1.at(t).Fill(met_,w);
	    NbJetsa1a1.at(t).Fill(n_jets_,w);

	    if ((HadPions_minus!=HadPions_plus) && (HadPions_minus!=VectZeroLV) && (HadPions_plus!=VectZeroLV) && (HadRefitPions_minus!=HadRefitPions_plus) && (HadRefitPions_minus!=VectZeroLV) && (HadRefitPions_plus!=VectZeroLV))
	      {
		if(TauminusPairConstraintWithTracksBS!=TauplusPairConstraintWithTracksBS && TauminusPairConstraintWithTracksBS!=zeroLV && TauplusPairConstraintWithTracksBS!=zeroLV)
		  {
		    if(ScalcPVRefitWithTracksBS.isOk("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus)==true){polarimetricAcopAnglePVRefitWithTracksBSMVADM.at(t).Fill(ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);
		    
		      polarimetricAcopAnglePVRefitWithTracksBSMVADM_DP.at(t).Fill(acop,Wspin);
		    }
		    if(max_pair.second==0)
		      {
			

			if(ScalcPVRefitWithTracksBS.isOk("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus)==true){
			  HiggsBDTScorea1a1.at(t).Fill(max_pair.first,w);
			  HiggsBDTScorea1a1_DP.at(t).Fill(max_pair.first,w);//AC?
			  
			  if(max_pair.first>0 && max_pair.first<0.7){polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled.at(t).Fill(ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

			    polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled_DP.at(t).Fill(acop,Wspin);
			  }
			  if(max_pair.first>0.7 && max_pair.first<0.8){polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled.at(t).Fill(2*TMath::Pi()+ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

			    polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled_DP.at(t).Fill(2*TMath::Pi()+acop,Wspin);
			  }
			  if(max_pair.first>0.8 && max_pair.first<1.){polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled.at(t).Fill(2*2*TMath::Pi()+ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

			    polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled_DP.at(t).Fill(2*2*TMath::Pi()+acop,Wspin);
			  }	
			  if(max_pair.first>0.7 && max_pair.first<1.){polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled.at(t).Fill(3*2*TMath::Pi()+ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);
			    polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled_DP.at(t).Fill(3*2*TMath::Pi()+acop,Wspin);
			  }
			}
			
		      }
		    if(max_pair.second==1)
		      {
			if(ScalcPVRefitWithTracksBS.isOk("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus)==true){
			  JetFakesBDTScorea1a1.at(t).Fill(max_pair.first,w);
			  JetFakesBDTScorea1a1_DP.at(t).Fill(max_pair.first,w);//AC?
			  
			  if(max_pair.first>0 && max_pair.first<0.7){polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled.at(t).Fill(ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);
			    polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled_DP.at(t).Fill(acop,Wspin);
			  }
			  if(max_pair.first>0.7 && max_pair.first<0.8){polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled.at(t).Fill(2*TMath::Pi()+ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);
			    polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled_DP.at(t).Fill(2*TMath::Pi()+acop,Wspin);
			  }
			  if(max_pair.first>0.8 && max_pair.first<1.){polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled.at(t).Fill(2*2*TMath::Pi()+ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);
			    polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled_DP.at(t).Fill(2*2*TMath::Pi()+acop,Wspin);
			  }
			  
			}
	     
		      }
		    if(max_pair.second==2)
		      {
			if(ScalcPVRefitWithTracksBS.isOk("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus)==true){
			  
			  ZTTBDTScorea1a1.at(t).Fill(max_pair.first,w);
			  ZTTBDTScorea1a1_DP.at(t).Fill(max_pair.first,w);//AC?
			  if(max_pair.first>0 && max_pair.first<0.7){polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled.at(t).Fill(ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);
			    polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled_DP.at(t).Fill(acop,Wspin);
			  }
			  if(max_pair.first>0.7 && max_pair.first<0.8){polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled.at(t).Fill(2*TMath::Pi()+ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);
			    polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled_DP.at(t).Fill(2*TMath::Pi()+acop,Wspin);
			  }
			  if(max_pair.first>0.8 && max_pair.first<1.){polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled.at(t).Fill(2*2*TMath::Pi()+ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);
			    polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled_DP.at(t).Fill(2*2*TMath::Pi()+acop,Wspin);
			  
			  }
			}
		      }

		    // Tau1PTa1a1.at(t).Fill(Tau1P4.Pt(),w);
		    // Tau1Ea1a1.at(t).Fill(Tau1P4.E(),w);
		    // Tau1Massa1a1.at(t).Fill(Tau1P4.M(),w);
		    // Tau1Phia1a1.at(t).Fill(Tau1P4.Phi(),w);
		    // Tau1Etaa1a1.at(t).Fill(Tau1P4.Eta(),w);
		    // Tau1dza1a1.at(t).Fill(Ntp->dz(Tau1),w);
		    // Tau2PTa1a1.at(t).Fill(Tau2P4.Pt(),w);
		    // Tau2Ea1a1.at(t).Fill(Tau2P4.E(),w);
		    // Tau2Massa1a1.at(t).Fill(Tau2P4.M(),w);
		    // Tau2Phia1a1.at(t).Fill(Tau2P4.Phi(),w);
		    // Tau2Etaa1a1.at(t).Fill(Tau2P4.Eta(),w);
		    // Tau2dza1a1.at(t).Fill(Ntp->dz(Tau2),w);
		    // TauTauVisMassa1a1.at(t).Fill((Tau1P4+Tau2P4).M(),w);
		    // TauTauVisPTa1a1.at(t).Fill((Tau1P4+Tau2P4).Pt(),w);
		    // TauTauFullMassa1a1.at(t).Fill((tau1P4+tau2P4).M(),w);
		    // TauTauFullPTa1a1.at(t).Fill((tau1P4+tau2P4).Pt(),w);
		    // PUPPImetcorra1a1.at(t).Fill(sqrt(PUPPImetCorr_px*PUPPImetCorr_px+PUPPImetCorr_py*PUPPImetCorr_py),w);
		    // PUPPImetcorrphia1a1.at(t).Fill(TMath::ATan2(PUPPImetCorr_py,PUPPImetCorr_px),w);
		    // NbJetsa1a1.at(t).Fill(Ntp->njets(IndexSelectedTemp.at(Sorted.back())),w);
		    // if(Has2jets>=0)
		    //   {
		    // 	Mjja1a1.at(t).Fill(Ntp->mjj(Has2jets),w);
		    // 	dijetpta1a1.at(t).Fill(Ntp->dijetpt(Has2jets),w);
		    // 	dijetphia1a1.at(t).Fill(Ntp->dijetphi(Has2jets),w);
		    // 	jdetaa1a1.at(t).Fill(Ntp->jdeta(Has2jets),w);
		    // 	jdphia1a1.at(t).Fill(Ntp->jdphi(Has2jets),w);
		    // 	jpt_2a1a1.at(t).Fill(Ntp->jpt_2(Has2jets),w);
		    // 	jeta_2a1a1.at(t).Fill(Ntp->jeta_2(Has2jets),w);
		    // 	jphi_2a1a1.at(t).Fill(Ntp->jphi_2(Has2jets),w);
		    //   }
		    // if(Has1jet>=0)
		    //   {
		    // 	jpt_1a1a1.at(t).Fill(Ntp->jpt_1(Has1jet),w);
		    // 	jeta_1a1a1.at(t).Fill(Ntp->jeta_1(Has1jet),w);
		    // 	jphi_1a1a1.at(t).Fill(Ntp->jphi_1(Has1jet),w);
		    //   }
		  }
	      }
	  }
      }

    //HPtVis.at(t).Fill((Tau1P4+Tau2P4).Pt(),w);
    //TauTauVisPT.at(t).Fill((Tau1P4+Tau2P4).Pt(),w);
    
    // std::vector<int> decay0;
    // std::vector<int> decay1;
    // bool a10=false;
    // bool a11=false;
    // bool a1pi00=false;
    // bool a1pi01=false;
    // bool rho0=false;
    // bool rho1=false;
    // bool pi0=false;
    // bool pi1=false;
    // bool pi23pi00=false;
    // bool pi23pi01=false;
    // bool e0=false;
    // bool e1=false;
    // bool mu0=false;
    // bool mu1=false;

    // bool purityDM=false;
    // bool purityNewMVA=false;

   
    // if(!Ntp->isData()){
    //   for(int i=0;i<Ntp->NMCTauDecayProducts(0);i++)
    // 	{
    // 	  decay0.push_back(Ntp->MCTauandProd_pdgid(0,i));
    // 	}
    //   for(int i=0;i<Ntp->NMCTauDecayProducts(1);i++)
    // 	{
    // 	  decay1.push_back(Ntp->MCTauandProd_pdgid(1,i));
    // 	}
    
    //   if(((((count(decay0.begin(), decay0.end(), 211)==2) && (count(decay0.begin(), decay0.end(), -211)==1)) ||((count(decay0.begin(), decay0.end(), -211)==2) && (count(decay0.begin(), decay0.end(), 211)==1))) && (count(decay0.begin(), decay0.end(), 111)==0))==true) a10=true;
    //   if((((count(decay0.begin(), decay0.end(), 111)==1) && (count(decay0.begin(), decay0.end(), 211)==1) && (count(decay0.begin(), decay0.end(), -211)==0)) ||((count(decay0.begin(), decay0.end(), 111)==1) && (count(decay0.begin(), decay0.end(), -211)==1) && (count(decay0.begin(), decay0.end(), 211)==0)))==true)rho0=true;
    //   if((((count(decay0.begin(), decay0.end(), 211)==1 && count(decay0.begin(), decay0.end(), -211)==0) ||(count(decay0.begin(), decay0.end(), -211)==1 && count(decay0.begin(), decay0.end(), 211)==0)) && count(decay0.begin(), decay0.end(), 111)==0)==true)pi0=true;
    //   if((((count(decay0.begin(), decay0.end(), 11)==1 && count(decay0.begin(), decay0.end(), -12)==1) ||(count(decay0.begin(), decay0.end(), -11)==1 && count(decay0.begin(), decay0.end(), 12)==1)) && count(decay0.begin(), decay0.end(), 111)==0) ==true)e0=true;
    //   if((((count(decay0.begin(), decay0.end(), 13)==1 && count(decay0.begin(), decay0.end(), -14)==1) ||(count(decay0.begin(), decay0.end(), -13)==1 && count(decay0.begin(), decay0.end(), 14)==1)) && count(decay0.begin(), decay0.end(), 111)==0) ==true)mu0=true;
    
    //   if(((((count(decay0.begin(), decay0.end(), 211)==2) && (count(decay0.begin(), decay0.end(), -211)==1)) ||((count(decay0.begin(), decay0.end(), -211)==2) && (count(decay0.begin(), decay0.end(), 211)==1))) && (count(decay0.begin(), decay0.end(), 111)==1))==true) a1pi00=true;
    //   if((((count(decay0.begin(), decay0.end(), 211)==1 && count(decay0.begin(), decay0.end(), -211)==0) ||(count(decay0.begin(), decay0.end(), -211)==1 && count(decay0.begin(), decay0.end(), 211)==0)) && ((count(decay0.begin(), decay0.end(), 111)==2) || count(decay0.begin(), decay0.end(), 111)==3))==true)pi23pi00=true;

    //   if(((((count(decay1.begin(), decay1.end(), 211)==2) && (count(decay1.begin(), decay1.end(), -211)==1)) ||((count(decay1.begin(), decay1.end(), -211)==2) && (count(decay1.begin(), decay1.end(), 211)==1))) && count(decay1.begin(), decay1.end(), 111)==0)==true) a11=true;
    //   if((((count(decay1.begin(), decay1.end(), 111)==1) && (count(decay1.begin(), decay1.end(), 211)==1) && (count(decay1.begin(), decay1.end(), -211)==0)) ||((count(decay1.begin(), decay1.end(), 111)==1) && (count(decay1.begin(), decay1.end(), -211)==1) && (count(decay1.begin(), decay1.end(), 211)==0)))==true)rho1=true;
    //   if((((count(decay1.begin(), decay1.end(), 211)==1 && count(decay1.begin(), decay1.end(), -211)==0) ||(count(decay1.begin(), decay1.end(), -211)==1 && count(decay1.begin(), decay1.end(), 211)==0)) && count(decay1.begin(), decay1.end(), 111)==0)==true)pi1=true;
    //   if((((count(decay1.begin(), decay1.end(), 11)==1 && count(decay1.begin(), decay1.end(), -12)==1) ||(count(decay1.begin(), decay1.end(), -11)==1 && count(decay1.begin(), decay1.end(), 12)==1)) && count(decay1.begin(), decay1.end(), 111)==0) ==true)e1=true;
    //   if((((count(decay1.begin(), decay1.end(), 13)==1 && count(decay1.begin(), decay1.end(), -14)==1) ||(count(decay1.begin(), decay1.end(), -13)==1 && count(decay1.begin(), decay1.end(), 14)==1)) && count(decay1.begin(), decay1.end(), 111)==0) ==true)mu1=true;
     
    //   if(((((count(decay1.begin(), decay1.end(), 211)==2) && (count(decay1.begin(), decay1.end(), -211)==1)) ||((count(decay1.begin(), decay1.end(), -211)==2) && (count(decay1.begin(), decay1.end(), 211)==1))) && (count(decay1.begin(), decay1.end(), 111)==1))==true) a1pi01=true;
    //   if((((count(decay1.begin(), decay1.end(), 211)==1 && count(decay1.begin(), decay1.end(), -211)==0) ||(count(decay1.begin(), decay1.end(), -211)==1 && count(decay1.begin(), decay1.end(), 211)==0)) && ((count(decay1.begin(), decay1.end(), 111)==2) || count(decay1.begin(), decay1.end(), 111))==3)==true)pi23pi01=true;
    // }
    
    // bool a1a1_=(a10 && a11);
    // bool a1rho=((a10 && rho1) ||(a11 && rho0));
    // bool a1pi=((a10 && pi1)||(a11 && pi0));
    // bool a1e=((a10 && e1)||(a11 && e0));
    // bool a1mu=((a10 && mu1)||(a11 && mu0));
    // bool rhorho=(rho0 && rho1);
    // bool rhopi=((rho0 && pi1)||(rho1 && pi0));
    // bool rhoe=((rho0 && e1)||(rho1 && e0));
    // bool rhomu=((rho0 && mu1)||(rho1 && mu0));
    // bool pipi=(pi0 && pi1);
    // bool pie=((pi0 && e1)||(pi1 && e0));
    // bool pimu=((pi0 && mu1)||(pi1 && mu0));
    // bool ee=(e0 && e1);
    // bool emu=((e0 && mu1)||(e1 && mu0));
    // bool mumu=(mu0 && mu1);

    // bool a1pi0a1pi0=(a1pi00 && a1pi01);
    // bool pi23pi0pi23pi0=(pi23pi00 && pi23pi01);
    
    // bool a1a1pi0=((a10 && a1pi01) || (a11 && a1pi00));
    // bool a1pi23pi0=((a10 && pi23pi01) || (a11 && pi23pi00));

    // bool a1pi0pi23pi0=((a1pi00 && pi23pi01) || (a1pi01 && pi23pi00));

    // bool a1pi0rho=((rho0 && a1pi01) || (rho1 && a1pi00));
    // bool rhopi23pi0=((rho0 && pi23pi01) || (rho1 && pi23pi00));

    // bool a1pi0pi=((pi0 && a1pi01) || (pi1 && a1pi00));
    // bool pipi23pi0=((pi0 && pi23pi01) || (pi1 && pi23pi00));

    // bool a1pi0e=((e0 && a1pi01) || (e1 && a1pi00));
    // bool pi23pi0e=((e0 && pi23pi01) || (e1 && pi23pi00));

    // bool a1pi0mu=((mu0 && a1pi01) || (mu1 && a1pi00));
    // bool pi23pi0mu=((mu0 && pi23pi01) || (mu1 && pi23pi00));
    
    // if (std::isnan(Wspin)!=true)
    //   {
    // 	if(a1a1)
    // 	  {
    // 	    if ((HadPions_minus!=HadPions_plus) && (HadPions_minus!=VectZeroLV) && (HadPions_plus!=VectZeroLV) && (HadRefitPions_minus!=HadRefitPions_plus) && (HadRefitPions_minus!=VectZeroLV) && (HadRefitPions_plus!=VectZeroLV))
    // 	      {
    // 		if(TauminusPairConstraint!=TauplusPairConstraint && TauminusPairConstraint!=zeroLV && TauplusPairConstraint!=zeroLV )
    // 		  {

    // 		    polarimetricAcopAngle.at(t).Fill(Scalc.AcopAngle("a1", "a1", TauminusPairConstraint, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraint , HadPions_plus, HadPionsCharge_plus),Wspin);

    // 		    if(a1a1_ ||a1rho  ||a1pi  || a1e || a1mu || rhorho || rhopi || rhoe || rhomu ||pipi || pie||pimu || ee||emu ||mumu || a1pi0a1pi0||pi23pi0pi23pi0 ||a1a1pi0||a1pi0pi23pi0 ||a1pi23pi0 ||a1pi0rho ||rhopi23pi0 ||a1pi0pi ||pipi23pi0 || a1pi0e||pi23pi0e ||a1pi0mu ||pi23pi0mu)purityDM=true;
	    
    // 		    if(!purityDM) { PurityDM.at(t).Fill(28.,w); //other
    // 		      // cout<<endl;
    // 		      // cout<<" pdgid0: ";
    // 		      // for(int i=0;i<Ntp->NMCTauDecayProducts(0);i++)
    // 		      //   {
    // 		      //     cout<<Ntp->MCTauandProd_pdgid(0,i)<<"  ";
		  
    // 		      //   }
    // 		      // cout<<endl;
    // 		      // cout<<" pdgid1: ";
    // 		      // for(int i=0;i<Ntp->NMCTauDecayProducts(1);i++)
    // 		      //   {
    // 		      //     cout<<Ntp->MCTauandProd_pdgid(1,i)<<"  ";
		  
    // 		      //   }
    // 		      // cout<<endl;
    // 		      //cout<<a1a1_ <<a1rho  <<a1pi  << a1e << a1mu << rhorho << rhopi << rhoe << rhomu <<pipi << pie<<pimu << ee<<emu <<mumu<<endl;
    // 		      //cout<<endl;
    // 		    }
    // 		    if(purityDM)
    // 		      {
    // 			if (a1a1_){PurityDM.at(t).Fill(0.,w);
    // 			  // cout<<endl;
    // 			  // cout<<" a1a1 pdgid0: ";
    // 			  // for(int i=0;i<Ntp->NMCTauDecayProducts(0);i++)
    // 			  // 	{
    // 			  // 	  cout<<Ntp->MCTauandProd_pdgid(0,i)<<"  ";
		  
    // 			  // 	}
    // 			  // cout<<endl;
    // 			  // cout<<" a1a1 pdgid1: ";
    // 			  // for(int i=0;i<Ntp->NMCTauDecayProducts(1);i++)
    // 			  // 	{
    // 			  // 	  cout<<Ntp->MCTauandProd_pdgid(1,i)<<"  ";
		  
    // 			  // 	}
    // 			  // cout<<endl;
    // 			}
				    
    // 			else if(a1a1pi0)PurityDM.at(t).Fill(1.,w);
    // 			else if (a1rho){PurityDM.at(t).Fill(2.,w);
    // 			  // cout<<endl;
    // 			  // for(int i=0;i<Ntp->NMCTauDecayProducts(0);i++)
    // 			  // 	{
    // 			  // 	  cout<<" a1rho pdgid0: "<<Ntp->MCTauandProd_pdgid(0,i);
		  
    // 			  // 	}
    // 			  // cout<<endl;
    // 			  // for(int i=0;i<Ntp->NMCTauDecayProducts(1);i++)
    // 			  // 	{
    // 			  // 	  cout<<" a1rho pdgid1: "<<Ntp->MCTauandProd_pdgid(1,i);
		  
    // 			  // 	}
    // 			  // cout<<endl;
    // 			}
    // 			else if (a1pi)PurityDM.at(t).Fill(3.,w);
    // 			else if (a1pi23pi0)PurityDM.at(t).Fill(4.,w);
    // 			else if (a1e)PurityDM.at(t).Fill(5.,w);
    // 			else if (a1mu)PurityDM.at(t).Fill(6.,w);

    // 			else if (a1pi0a1pi0)PurityDM.at(t).Fill(7.,w);
    // 			else if (a1pi0rho)PurityDM.at(t).Fill(8.,w);
    // 			else if (a1pi0pi)PurityDM.at(t).Fill(9.,w);
    // 			else if (a1pi0pi23pi0)PurityDM.at(t).Fill(10.,w);
    // 			else if (a1pi0e)PurityDM.at(t).Fill(11.,w);
    // 			else if (a1pi0mu)PurityDM.at(t).Fill(12.,w);
				    
    // 			else if (rhorho)PurityDM.at(t).Fill(13.,w);
    // 			else if (rhopi)PurityDM.at(t).Fill(14.,w);
    // 			else if (rhopi23pi0)PurityDM.at(t).Fill(15.,w);
    // 			else if (rhoe)PurityDM.at(t).Fill(16.,w);
    // 			else if (rhomu)PurityDM.at(t).Fill(17.,w);
				    
    // 			else if (pipi)PurityDM.at(t).Fill(18.,w);
    // 			else if (pipi23pi0)PurityDM.at(t).Fill(19.,w);
    // 			else if (pie)PurityDM.at(t).Fill(20.,w);
    // 			else if (pimu)PurityDM.at(t).Fill(21.,w);
    // 			else if (pi23pi0pi23pi0)PurityDM.at(t).Fill(22.,w);
    // 			else if (pi23pi0e)PurityDM.at(t).Fill(23.,w);
    // 			else if (pi23pi0mu)PurityDM.at(t).Fill(24.,w);
    // 			else if (ee)PurityDM.at(t).Fill(25.,w);
    // 			else if (emu)PurityDM.at(t).Fill(26.,w);
    // 			else if (mumu)PurityDM.at(t).Fill(27.,w);
    // 		      }
    // 		    //}
    // 		  }
    // 	      }
	
    // 	    if(TauminusPairConstraint!=zeroLV && TauplusPairConstraint!=zeroLV && Scalc.isOk("a1", "a1", TauminusPairConstraint, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraint, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAngle.at(t).Fill(Scalc.AcopAngle("a1", "a1", TauminusPairConstraint, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraint, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);
	
    // 	    if(TauminusPairConstraintNoBS!=zeroLV && TauplusPairConstraintNoBS!=zeroLV && ScalcPVRefitNoBS.isOk("a1", "a1", TauminusPairConstraintNoBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBS, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitNoBS.at(t).Fill(ScalcPVRefitNoBS.AcopAngle("a1", "a1", TauminusPairConstraintNoBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBS, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

    // 	    if(TauminusPairConstraintBS!=zeroLV && TauplusPairConstraintBS!=zeroLV && ScalcPVRefitBS.isOk("a1", "a1", TauminusPairConstraintBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBS, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitBS.at(t).Fill(ScalcPVRefitBS.AcopAngle("a1", "a1", TauminusPairConstraintBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBS, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

    // 	    if(TauminusPairConstraintNoBSZNominal!=zeroLV && TauplusPairConstraintNoBSZNominal!=zeroLV && ScalcPVRefitNoBSZNominal.isOk("a1", "a1", TauminusPairConstraintNoBSZNominal, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominal, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitNoBSZNominal.at(t).Fill(ScalcPVRefitNoBSZNominal.AcopAngle("a1", "a1", TauminusPairConstraintNoBSZNominal, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominal, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

    // 	    if(TauminusPairConstraintBSZNominal!=zeroLV && TauplusPairConstraintBSZNominal!=zeroLV && ScalcPVRefitBSZNominal.isOk("a1", "a1", TauminusPairConstraintBSZNominal, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominal, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitBSZNominal.at(t).Fill(ScalcPVRefitBSZNominal.AcopAngle("a1", "a1", TauminusPairConstraintBSZNominal, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominal, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

    // 	  }
    // 	if(a1a1MVA)
    // 	  {	      
    // 	    if(TauminusPairConstraint!=zeroLV && TauplusPairConstraint!=zeroLV && Scalc.isOk("a1", "a1", TauminusPairConstraint, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraint , HadPions_plus, HadPionsCharge_plus)==true) polarimetricAcopAngleMVADM.at(t).Fill(Scalc.AcopAngle("a1", "a1", TauminusPairConstraint, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraint, HadPions_plus, HadPionsCharge_plus),Wspin);
		
    // 	    if(TauminusPairConstraintNoBS!=zeroLV && TauplusPairConstraintNoBS!=zeroLV && ScalcPVRefitNoBS.isOk("a1", "a1", TauminusPairConstraintNoBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBS, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitNoBSMVADM.at(t).Fill(ScalcPVRefitNoBS.AcopAngle("a1", "a1", TauminusPairConstraintNoBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBS, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

    // 	    if(TauminusPairConstraintBS!=zeroLV && TauplusPairConstraintBS!=zeroLV && ScalcPVRefitBS.isOk("a1", "a1", TauminusPairConstraintBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBS, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitBSMVADM.at(t).Fill(ScalcPVRefitBS.AcopAngle("a1", "a1", TauminusPairConstraintBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBS, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

    // 	    if(TauminusPairConstraintNoBSZNominal!=zeroLV && TauplusPairConstraintNoBSZNominal!=zeroLV && ScalcPVRefitNoBSZNominal.isOk("a1", "a1", TauminusPairConstraintNoBSZNominal, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominal, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitNoBSZNominalMVADM.at(t).Fill(ScalcPVRefitNoBSZNominal.AcopAngle("a1", "a1", TauminusPairConstraintNoBSZNominal, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominal, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

    // 	    if(TauminusPairConstraintBSZNominal!=zeroLV && TauplusPairConstraintBSZNominal!=zeroLV && ScalcPVRefitBSZNominal.isOk("a1", "a1", TauminusPairConstraintBSZNominal, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominal, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitBSZNominalMVADM.at(t).Fill(ScalcPVRefitBSZNominal.AcopAngle("a1", "a1", TauminusPairConstraintBSZNominal, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominal, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

    // 	    if(TauminusPairConstraintWithTracksBS!=zeroLV && TauplusPairConstraintWithTracksBS!=zeroLV && ScalcPVRefitWithTracksBS.isOk("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitWithTracksBSMVADM.at(t).Fill(ScalcPVRefitWithTracksBS.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBS, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBS, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

    // 	    if(TauminusPairConstraintWithTracksBSZNominal!=zeroLV && TauplusPairConstraintWithTracksBSZNominal!=zeroLV && ScalcPVRefitWithTracksBSZNominal.isOk("a1", "a1", TauminusPairConstraintWithTracksBSZNominal, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBSZNominal, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADM.at(t).Fill(ScalcPVRefitWithTracksBSZNominal.AcopAngle("a1", "a1", TauminusPairConstraintWithTracksBSZNominal, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintWithTracksBSZNominal, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);
		

    // 	    if(!Ntp->isData()){
    // 	      if(a1a1_ ||a1rho  ||a1pi  || a1e || a1mu || rhorho || rhopi || rhoe || rhomu ||pipi || pie||pimu || ee||emu ||mumu || a1pi0a1pi0||pi23pi0pi23pi0 ||a1a1pi0||a1pi0pi23pi0 ||a1pi23pi0 ||a1pi0rho ||rhopi23pi0 ||a1pi0pi ||pipi23pi0 || a1pi0e||pi23pi0e ||a1pi0mu ||pi23pi0mu)purityNewMVA=true;

    // 	      if(!purityNewMVA) { PurityNewMVA.at(t).Fill(28.,w); //other
    // 		// cout<<endl;
    // 		// cout<<" pdgid0: ";
    // 		// for(int i=0;i<Ntp->NMCTauDecayProducts(0);i++)
    // 		//   {
    // 		//     cout<<Ntp->MCTauandProd_pdgid(0,i)<<"  ";
		  
    // 		//   }
    // 		// cout<<endl;
    // 		// cout<<" pdgid1: ";
    // 		// for(int i=0;i<Ntp->NMCTauDecayProducts(1);i++)
    // 		//   {
    // 		//     cout<<Ntp->MCTauandProd_pdgid(1,i)<<"  ";
		  
    // 		//   }
    // 		// cout<<endl;
    // 		//cout<<a1a1_ <<a1rho  <<a1pi  << a1e << a1mu << rhorho << rhopi << rhoe << rhomu <<pipi << pie<<pimu << ee<<emu <<mumu<<endl;
    // 		//cout<<endl;
    // 	      }
    // 	      if(purityNewMVA)
    // 		{
    // 		  if (a1a1_){PurityNewMVA.at(t).Fill(0.,w);
    // 		    // cout<<endl;
    // 		    // cout<<" a1a1 pdgid0: ";
    // 		    // for(int i=0;i<Ntp->NMCTauDecayProducts(0);i++)
    // 		    // 	{
    // 		    // 	  cout<<Ntp->MCTauandProd_pdgid(0,i)<<"  ";
		  
    // 		    // 	}
    // 		    // cout<<endl;
    // 		    // cout<<" a1a1 pdgid1: ";
    // 		    // for(int i=0;i<Ntp->NMCTauDecayProducts(1);i++)
    // 		    // 	{
    // 		    // 	  cout<<Ntp->MCTauandProd_pdgid(1,i)<<"  ";
		  
    // 		    // 	}
    // 		    // cout<<endl;
    // 		  }
				    
    // 		  else if(a1a1pi0)PurityNewMVA.at(t).Fill(1.,w);
    // 		  else if (a1rho){PurityNewMVA.at(t).Fill(2.,w);
    // 		    // cout<<endl;
    // 		    // for(int i=0;i<Ntp->NMCTauDecayProducts(0);i++)
    // 		    // 	{
    // 		    // 	  cout<<" a1rho pdgid0: "<<Ntp->MCTauandProd_pdgid(0,i);
		  
    // 		    // 	}
    // 		    // cout<<endl;
    // 		    // for(int i=0;i<Ntp->NMCTauDecayProducts(1);i++)
    // 		    // 	{
    // 		    // 	  cout<<" a1rho pdgid1: "<<Ntp->MCTauandProd_pdgid(1,i);
		  
    // 		    // 	}
    // 		    // cout<<endl;
    // 		  }
    // 		  else if (a1pi)PurityNewMVA.at(t).Fill(3.,w);
    // 		  else if (a1pi23pi0)PurityNewMVA.at(t).Fill(4.,w);
    // 		  else if (a1e)PurityNewMVA.at(t).Fill(5.,w);
    // 		  else if (a1mu)PurityNewMVA.at(t).Fill(6.,w);

    // 		  else if (a1pi0a1pi0)PurityNewMVA.at(t).Fill(7.,w);
    // 		  else if (a1pi0rho)PurityNewMVA.at(t).Fill(8.,w);
    // 		  else if (a1pi0pi)PurityNewMVA.at(t).Fill(9.,w);
    // 		  else if (a1pi0pi23pi0)PurityNewMVA.at(t).Fill(10.,w);
    // 		  else if (a1pi0e)PurityNewMVA.at(t).Fill(11.,w);
    // 		  else if (a1pi0mu)PurityNewMVA.at(t).Fill(12.,w);
				    
    // 		  else if (rhorho)PurityNewMVA.at(t).Fill(13.,w);
    // 		  else if (rhopi)PurityNewMVA.at(t).Fill(14.,w);
    // 		  else if (rhopi23pi0)PurityNewMVA.at(t).Fill(15.,w);
    // 		  else if (rhoe)PurityNewMVA.at(t).Fill(16.,w);
    // 		  else if (rhomu)PurityNewMVA.at(t).Fill(17.,w);
				    
    // 		  else if (pipi)PurityNewMVA.at(t).Fill(18.,w);
    // 		  else if (pipi23pi0)PurityNewMVA.at(t).Fill(19.,w);
    // 		  else if (pie)PurityNewMVA.at(t).Fill(20.,w);
    // 		  else if (pimu)PurityNewMVA.at(t).Fill(21.,w);
    // 		  else if (pi23pi0pi23pi0)PurityNewMVA.at(t).Fill(22.,w);
    // 		  else if (pi23pi0e)PurityNewMVA.at(t).Fill(23.,w);
    // 		  else if (pi23pi0mu)PurityNewMVA.at(t).Fill(24.,w);
    // 		  else if (ee)PurityNewMVA.at(t).Fill(25.,w);
    // 		  else if (emu)PurityNewMVA.at(t).Fill(26.,w);
    // 		  else if (mumu)PurityNewMVA.at(t).Fill(27.,w);
    // 		}
    // 	    }
    // 	    //}
    // 	  }

    //   }



    // if(!Ntp->isData())
    //   {
    // 	if(a1a1_ && Ntp->CheckDecayID(5,5))
    // 	  {
    // 	    Pions1=Ntp->GetTruthPionsFromA1(0);
    // 	    Pions1Charge.push_back(1);
    // 	    Pions1Charge.push_back(-1);
    // 	    Pions1Charge.push_back(-1);
    // 	    Pions2=Ntp->GetTruthPionsFromA1(1);
    // 	    Pions2Charge.push_back(-1);
    // 	    Pions2Charge.push_back(1);
    // 	    Pions2Charge.push_back(1);
    // 	    SCalculator Scalc1Truth("a1");
    // 	    SCalculator Scalc2Truth("a1");
    // 	    Scalc1Truth.SortPions(Pions1, Pions1Charge);
    // 	    Scalc2Truth.SortPions(Pions2, Pions2Charge);
    // 	    if(Ntp->MCTau_pdgid(0)==15)
    // 	      {
    // 		tauandprodTruthminus.push_back(Ntp->GetTruthTauLV(5,0));
    // 		tauandprodTruthminus.push_back(Pions1.at(0));
    // 		tauandprodTruthminus.push_back(Pions1.at(1));
    // 		tauandprodTruthminus.push_back(Pions1.at(2));
    // 		tauandprodTruthplus.push_back(Ntp->GetTruthTauLV(5,1));   
    // 		tauandprodTruthplus.push_back(Pions2.at(0));   
    // 		tauandprodTruthplus.push_back(Pions2.at(1));   
    // 		tauandprodTruthplus.push_back(Pions2.at(2));   
    // 	      }
    // 	    else if (Ntp->MCTau_pdgid(0)==-15)
    // 	      {
    // 		tauandprodTruthminus.push_back(Ntp->GetTruthTauLV(5,1));
    // 		tauandprodTruthminus.push_back(Pions2.at(0));
    // 		tauandprodTruthminus.push_back(Pions2.at(1));
    // 		tauandprodTruthminus.push_back(Pions2.at(2));
    // 		tauandprodTruthplus.push_back(Ntp->GetTruthTauLV(5,0));
    // 		tauandprodTruthplus.push_back(Pions1.at(0));   
    // 		tauandprodTruthplus.push_back(Pions1.at(1));  
    // 		tauandprodTruthplus.push_back(Pions1.at(2));   
    // 	      }	
    // 	  }
    // 	if((Pions1!=Pions2) && (Pions1!=VectZeroLV) && (Pions2!=VectZeroLV) && tauandprodTruthminus.at(0)!=zeroLV &&tauandprodTruthplus.at(0)!=zeroLV &&tauandprodTruthminus.at(0)!=tauandprodTruthplus.at(0))
    // 	  {
    // 	    if(a1a1_  && Ntp->CheckDecayID(5,5))
    // 	      {
		    
    // 		SCalculator Scalc1Truth("a1");
    // 		SCalculator Scalc2Truth("a1");
		    
    // 		Scalc1Truth.Configure(tauandprodTruthminus,tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0), -1);
    // 		TVector3 h1Truth=Scalc1Truth.pv();
			
    // 		Scalc2Truth.Configure(tauandprodTruthplus,tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0), +1);
    // 		TVector3 h2Truth=Scalc2Truth.pv();
			
    // 		double h1TruthNorm=1./h1Truth.Mag();
    // 		double h2TruthNorm=1./h2Truth.Mag();
		    
    // 		if(std::isnan(h1TruthNorm)!=true && std::isnan(h2TruthNorm)!=true)
    // 		  { 
    // 		    TLorentzVector tauminusTruth_HRF = Scalc1Truth.Boost(tauandprodTruthminus.at(0),tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0));
    // 		    TLorentzVector tauplusTruth_HRF  = Scalc2Truth.Boost(tauandprodTruthplus.at(0),tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0));
			
    // 		    double norm1Truth=1./(((h1Truth*h1TruthNorm).Cross(tauminusTruth_HRF.Vect().Unit())).Mag());
    // 		    double norm2Truth=1./(((h2Truth*h2TruthNorm).Cross(tauplusTruth_HRF.Vect().Unit())).Mag());
    // 		    TVector3 k1Truth = ((h1Truth*h1TruthNorm).Cross(tauminusTruth_HRF.Vect().Unit()))*norm1Truth;
    // 		    TVector3 k2Truth = ((h2Truth*h2TruthNorm).Cross(tauplusTruth_HRF.Vect().Unit()))*norm2Truth;
			
    // 		    if(((h1Truth*h1TruthNorm).Cross(h2Truth*h2TruthNorm))*tauminusTruth_HRF.Vect().Unit()<=0) {polarimetricAcopAngleTruthA1.at(t).Fill(TMath::ATan2((k1Truth.Cross(k2Truth)).Mag(),k1Truth*k2Truth),Wspin);/*cout<<"Angle Truth: "<<TMath::ATan2((k1Truth.Cross(k2Truth)).Mag(),k1Truth*k2Truth)<<endl;*/}
    // 		    else{ polarimetricAcopAngleTruthA1.at(t).Fill(2*TMath::Pi()-TMath::ATan2((k1Truth.Cross(k2Truth)).Mag(),k1Truth*k2Truth),Wspin);/*cout<<"Angle Truth: "<<2*TMath::Pi()-TMath::ATan2((k1Truth.Cross(k2Truth)).Mag(),k1Truth*k2Truth)<<endl;*/}
				
    // 		  }
    // 	      }
    // 	  }
    // 	if(a1a1_  && Ntp->CheckDecayID(5,5))
    // 	  {
    // 	    PVXResol.at(t).Fill((Ntp->PVtx().X()-Ntp->PVtx_Gen().X())/Ntp->PVtx_Gen().X(),w);
    // 	    PVXNoBSResol.at(t).Fill((tauNoBSPrimaryVertex.X()-Ntp->PVtx_Gen().X())/Ntp->PVtx_Gen().X(),w);
    // 	    PVXBSResol.at(t).Fill((tauBSPrimaryVertex.X()-Ntp->PVtx_Gen().X())/Ntp->PVtx_Gen().X(),w);
    // 	    PVYResol.at(t).Fill((Ntp->PVtx().Y()-Ntp->PVtx_Gen().Y())/Ntp->PVtx_Gen().Y(),w);
    // 	    PVYNoBSResol.at(t).Fill((tauNoBSPrimaryVertex.Y()-Ntp->PVtx_Gen().Y())/Ntp->PVtx_Gen().Y(),w);
    // 	    PVYBSResol.at(t).Fill((tauBSPrimaryVertex.Y()-Ntp->PVtx_Gen().Y())/Ntp->PVtx_Gen().Y(),w);
    // 	    PVZResol.at(t).Fill((Ntp->PVtx().Z()-Ntp->PVtx_Gen().Z())/Ntp->PVtx_Gen().Z(),w);
    // 	    PVZNoBSResol.at(t).Fill((tauNoBSPrimaryVertex.Z()-Ntp->PVtx_Gen().Z())/Ntp->PVtx_Gen().Z(),w);
    // 	    PVZBSResol.at(t).Fill((tauBSPrimaryVertex.Z()-Ntp->PVtx_Gen().Z())/Ntp->PVtx_Gen().Z(),w);
    // 	    if(isRefitNoBS)PVXNoBSOnlyResol.at(t).Fill((tauNoBSPrimaryVertex.X()-Ntp->PVtx_Gen().X())/Ntp->PVtx_Gen().X(),w);
    // 	    if(isRefitBS)PVXBSOnlyResol.at(t).Fill((tauBSPrimaryVertex.X()-Ntp->PVtx_Gen().X())/Ntp->PVtx_Gen().X(),w);
    // 	    if(isRefitNoBS)PVYNoBSOnlyResol.at(t).Fill((tauNoBSPrimaryVertex.Y()-Ntp->PVtx_Gen().Y())/Ntp->PVtx_Gen().Y(),w);
    // 	    if(isRefitBS)PVYBSOnlyResol.at(t).Fill((tauBSPrimaryVertex.Y()-Ntp->PVtx_Gen().Y())/Ntp->PVtx_Gen().Y(),w);
    // 	    if(isRefitNoBS)PVZNoBSOnlyResol.at(t).Fill((tauNoBSPrimaryVertex.Z()-Ntp->PVtx_Gen().Z())/Ntp->PVtx_Gen().Z(),w);
    // 	    if(isRefitBS)PVZBSOnlyResol.at(t).Fill((tauBSPrimaryVertex.Z()-Ntp->PVtx_Gen().Z())/Ntp->PVtx_Gen().Z(),w);
    // 	    if(Ntp->MCTau_pdgid(0)==15 && Ntp->PFTau_secondaryVertex_pos_Size()>1 && Ntp->MCTauandProd_VertexSize()==2)//nplus==1 && nminus==2)
    // 	      {
    // 		ResolPullXVtxIna1a1.at(t).Fill((TauminusSecondaryVertex.X()-Ntp->MCTauandProd_Vertex(0,1).X())/Ntp->MCTauandProd_Vertex(0,1).X());
    // 		ResolPullXVtxIna1a1.at(t).Fill((TauplusSecondaryVertex.X()-Ntp->MCTauandProd_Vertex(1,1).X())/Ntp->MCTauandProd_Vertex(1,1).X());
    // 		ResolPullYVtxIna1a1.at(t).Fill((TauminusSecondaryVertex.Y()-Ntp->MCTauandProd_Vertex(0,1).Y())/Ntp->MCTauandProd_Vertex(0,1).Y());
    // 		ResolPullYVtxIna1a1.at(t).Fill((TauplusSecondaryVertex.Y()-Ntp->MCTauandProd_Vertex(1,1).Y())/Ntp->MCTauandProd_Vertex(1,1).Y());
    // 		ResolPullZVtxIna1a1.at(t).Fill((TauminusSecondaryVertex.Z()-Ntp->MCTauandProd_Vertex(0,1).Z())/Ntp->MCTauandProd_Vertex(0,1).Z());
    // 		ResolPullZVtxIna1a1.at(t).Fill((TauplusSecondaryVertex.Z()-Ntp->MCTauandProd_Vertex(1,1).Z())/Ntp->MCTauandProd_Vertex(1,1).Z());
    // 	      }
    // 	    else if(Ntp->MCTau_pdgid(0)==-15 && Ntp->PFTau_secondaryVertex_pos_Size()>1 && Ntp->MCTauandProd_VertexSize()==2)
    // 	      {
    // 		ResolPullXVtxIna1a1.at(t).Fill((TauminusSecondaryVertex.X()-Ntp->MCTauandProd_Vertex(1,1).X())/Ntp->MCTauandProd_Vertex(1,1).X());
    // 		ResolPullXVtxIna1a1.at(t).Fill((TauplusSecondaryVertex.X()-Ntp->MCTauandProd_Vertex(0,1).X())/Ntp->MCTauandProd_Vertex(0,1).X());
    // 		ResolPullYVtxIna1a1.at(t).Fill((TauminusSecondaryVertex.Y()-Ntp->MCTauandProd_Vertex(1,1).Y())/Ntp->MCTauandProd_Vertex(1,1).Y());
    // 		ResolPullYVtxIna1a1.at(t).Fill((TauplusSecondaryVertex.Y()-Ntp->MCTauandProd_Vertex(0,1).Y())/Ntp->MCTauandProd_Vertex(0,1).Y());
    // 		ResolPullZVtxIna1a1.at(t).Fill((TauminusSecondaryVertex.Z()-Ntp->MCTauandProd_Vertex(1,1).Z())/Ntp->MCTauandProd_Vertex(1,1).Z());
    // 		ResolPullZVtxIna1a1.at(t).Fill((TauplusSecondaryVertex.Z()-Ntp->MCTauandProd_Vertex(0,1).Z())/Ntp->MCTauandProd_Vertex(0,1).Z());
		    
    // 	      }
    //	  }
    //   }
  }
}
//}
//  This is a function if you want to do something after the event loop
void HCPTauTau::Finish() {

  if(mode == RECONSTRUCT) {
    SkimConfig SC;
    SC.ApplySkimEfficiency(types,Npassed, Npassed_noweight);

    double norm=1.;

    for(unsigned i=0;i<CrossSectionandAcceptance.size();i++){
      if(CrossSectionandAcceptance.at(i)>0 || HConfig.GetID(i)==36 || HConfig.GetID(i)==20 || HConfig.GetID(i)==23 || HConfig.GetID(i)==30 || HConfig.GetID(i)==33){
	if(CrossSectionandAcceptance.at(i)>0)norm= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
	else norm=1.;
	// PUPPImetcorr.at(1).Add(&PUPPImetcorrQCDMC.at(i),-norm);
	// PUPPImetcorrphi.at(1).Add(&PUPPImetcorrphiQCDMC.at(i),-norm);
	// Tau1PT.at(1).Add(&Tau1PTQCDMC.at(i),-norm);
	// Tau1E.at(1).Add(&Tau1EQCDMC.at(i),-norm);
	// Tau1Mass.at(1).Add(&Tau1MassQCDMC.at(i),-norm);
	// Tau1Phi.at(1).Add(&Tau1PhiQCDMC.at(i),-norm);
	// Tau1Eta.at(1).Add(&Tau1EtaQCDMC.at(i),-norm);
	// Tau1dz.at(1).Add(&Tau1dzQCDMC.at(i),-norm);
	// Tau1HPSDecayMode.at(1).Add(&Tau1HPSDecayModeQCDMC.at(i),-norm);
	// Tau1MVADecayMode.at(1).Add(&Tau1MVADecayModeQCDMC.at(i),-norm);
	// //	Tau1GenMatch.at(1).Add(&Tau1GenMatchQCDMC.at(i),-norm);
	// Tau2PT.at(1).Add(&Tau2PTQCDMC.at(i),-norm);
	// Tau2E.at(1).Add(&Tau2EQCDMC.at(i),-norm);
	// Tau2Mass.at(1).Add(&Tau2MassQCDMC.at(i),-norm);
	// Tau2Phi.at(1).Add(&Tau2PhiQCDMC.at(i),-norm);
	// Tau2Eta.at(1).Add(&Tau2EtaQCDMC.at(i),-norm);
	// Tau2dz.at(1).Add(&Tau2dzQCDMC.at(i),-norm);
	// Tau2HPSDecayMode.at(1).Add(&Tau2HPSDecayModeQCDMC.at(i),-norm);
	// Tau2MVADecayMode.at(1).Add(&Tau2MVADecayModeQCDMC.at(i),-norm);
	// //	Tau2GenMatch.at(1).Add(&Tau2GenMatchQCDMC.at(i),-norm);
	// NbJets.at(1).Add(&NbJetsQCDMC.at(i),-norm);
	// TauTauVisMass.at(1).Add(&TauTauVisMassQCDMC.at(i),-norm);
	// TauTauFullMass.at(1).Add(&TauTauFullMassQCDMC.at(i),-norm);
	// TauTauVisPT.at(1).Add(&TauTauVisPTQCDMC.at(i),-norm);
	// TauTauFullPT.at(1).Add(&TauTauFullPTQCDMC.at(i),-norm);
	// Mjj.at(1).Add(&MjjQCDMC.at(i),-norm);
	// dijetpt.at(1).Add(&dijetptQCDMC.at(i),-norm);
	// dijetphi.at(1).Add(&dijetphiQCDMC.at(i),-norm);
	// jdeta.at(1).Add(&jdetaQCDMC.at(i),-norm);
	// jdphi.at(1).Add(&jdphiQCDMC.at(i),-norm);
	// jpt_2.at(1).Add(&jpt_2QCDMC.at(i),-norm);
	// jeta_2.at(1).Add(&jeta_2QCDMC.at(i),-norm);
	// jphi_2.at(1).Add(&jphi_2QCDMC.at(i),-norm);
	// jpt_1.at(1).Add(&jpt_1QCDMC.at(i),-norm);
	// jeta_1.at(1).Add(&jeta_1QCDMC.at(i),-norm);
	// jphi_1.at(1).Add(&jphi_1QCDMC.at(i),-norm);
	
	// HiggsBDTScore.at(1).Add(&HiggsBDTScoreQCDMC.at(i),-norm);
	// JetFakesBDTScore.at(1).Add(&JetFakesBDTScoreQCDMC.at(i),-norm);
	// ZTTBDTScore.at(1).Add(&ZTTBDTScoreQCDMC.at(i),-norm);
	
	

	//polarimetricAcopAngleMVADM.at(1).Add(&polarimetricAcopAngleMVADMQCDMC.at(i),-norm);
	//polarimetricAcopAnglePVRefitBSMVADM.at(1).Add(&polarimetricAcopAnglePVRefitBSMVADMQCDMC.at(i),-norm);
	//polarimetricAcopAnglePVRefitBSZNominalMVADM.at(1).Add(&polarimetricAcopAnglePVRefitBSZNominalMVADMQCDMC.at(i),-norm);
	//polarimetricAcopAnglePVRefitWithTracksBSMVADM.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMQCDMC.at(i),-norm);
	//polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADM.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMQCDMC.at(i),-norm);

	//polarimetricAcopAngleMVADMHiggs.at(1).Add(&polarimetricAcopAngleMVADMHiggsQCDMC.at(i),-norm);
	//polarimetricAcopAnglePVRefitBSMVADMHiggs.at(1).Add(&polarimetricAcopAnglePVRefitBSMVADMHiggsQCDMC.at(i),-norm);
	//polarimetricAcopAnglePVRefitBSZNominalMVADMHiggs.at(1).Add(&polarimetricAcopAnglePVRefitBSZNominalMVADMHiggsQCDMC.at(i),-norm);

	polarimetricAcopAnglePVRefitWithTracksBSMVADM.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMQCDMC.at(i),-norm);
	polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggs.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsQCDMC.at(i),-norm);
	polarimetricAcopAnglePVRefitWithTracksBSMVADM_DP.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMQCDMC_DP.at(i),-norm);
	polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggs_DP.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsQCDMC_DP.at(i),-norm);
	//polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMHiggs.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMHiggsQCDMC.at(i),-norm);

	//polarimetricAcopAngleMVADMJetFakes.at(1).Add(&polarimetricAcopAngleMVADMJetFakesQCDMC.at(i),-norm);
	//polarimetricAcopAnglePVRefitBSMVADMJetFakes.at(1).Add(&polarimetricAcopAnglePVRefitBSMVADMJetFakesQCDMC.at(i),-norm);
	//polarimetricAcopAnglePVRefitBSZNominalMVADMJetFakes.at(1).Add(&polarimetricAcopAnglePVRefitBSZNominalMVADMJetFakesQCDMC.at(i),-norm);
	polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakes.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesQCDMC.at(i),-norm);
	polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakes_DP.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesQCDMC_DP.at(i),-norm);

	//polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMJetFakes.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMJetFakesQCDMC.at(i),-norm);

	//polarimetricAcopAngleMVADMZTT.at(1).Add(&polarimetricAcopAngleMVADMZTTQCDMC.at(i),-norm);
	//polarimetricAcopAnglePVRefitBSMVADMZTT.at(1).Add(&polarimetricAcopAnglePVRefitBSMVADMZTTQCDMC.at(i),-norm);
	//polarimetricAcopAnglePVRefitBSZNominalMVADMZTT.at(1).Add(&polarimetricAcopAnglePVRefitBSZNominalMVADMZTTQCDMC.at(i),-norm);
	polarimetricAcopAnglePVRefitWithTracksBSMVADMZTT.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTQCDMC.at(i),-norm);
	polarimetricAcopAnglePVRefitWithTracksBSMVADMZTT_DP.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTQCDMC_DP.at(i),-norm);

	//polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMZTT.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSZNominalMVADMZTTQCDMC.at(i),-norm);
	  
	
	PUPPImetcorra1a1.at(1).Add(&PUPPImetcorra1a1QCDMC.at(i),-norm);
	// PUPPImetcorrphia1a1.at(1).Add(&PUPPImetcorrphia1a1QCDMC.at(i),-norm);
	Tau1PTa1a1.at(1).Add(&Tau1PTa1a1QCDMC.at(i),-norm);
	// Tau1Ea1a1.at(1).Add(&Tau1Ea1a1QCDMC.at(i),-norm);
	// Tau1Massa1a1.at(1).Add(&Tau1Massa1a1QCDMC.at(i),-norm);
	// Tau1Phia1a1.at(1).Add(&Tau1Phia1a1QCDMC.at(i),-norm);
	// Tau1Etaa1a1.at(1).Add(&Tau1Etaa1a1QCDMC.at(i),-norm);
	// Tau1dza1a1.at(1).Add(&Tau1dza1a1QCDMC.at(i),-norm);
	
	// Tau2PTa1a1.at(1).Add(&Tau2PTa1a1QCDMC.at(i),-norm);
	// Tau2Ea1a1.at(1).Add(&Tau2Ea1a1QCDMC.at(i),-norm);
	// Tau2Massa1a1.at(1).Add(&Tau2Massa1a1QCDMC.at(i),-norm);
	// Tau2Phia1a1.at(1).Add(&Tau2Phia1a1QCDMC.at(i),-norm);
	// Tau2Etaa1a1.at(1).Add(&Tau2Etaa1a1QCDMC.at(i),-norm);
	// Tau2dza1a1.at(1).Add(&Tau2dza1a1QCDMC.at(i),-norm);

	NbJetsa1a1.at(1).Add(&NbJetsa1a1QCDMC.at(i),-norm);
	TauTauVisMassa1a1.at(1).Add(&TauTauVisMassa1a1QCDMC.at(i),-norm);
	TauTauFullMassa1a1.at(1).Add(&TauTauFullMassa1a1QCDMC.at(i),-norm);
	TauTauVisPTa1a1.at(1).Add(&TauTauVisPTa1a1QCDMC.at(i),-norm);
	TauTauFullPTa1a1.at(1).Add(&TauTauFullPTa1a1QCDMC.at(i),-norm);
	Mjja1a1.at(1).Add(&Mjja1a1QCDMC.at(i),-norm);
	// dijetpta1a1.at(1).Add(&dijetpta1a1QCDMC.at(i),-norm);
	// dijetphia1a1.at(1).Add(&dijetphia1a1QCDMC.at(i),-norm);
	jdetaa1a1.at(1).Add(&jdetaa1a1QCDMC.at(i),-norm);
	// jdphia1a1.at(1).Add(&jdphia1a1QCDMC.at(i),-norm);
	// jpt_2a1a1.at(1).Add(&jpt_2a1a1QCDMC.at(i),-norm);
	// jeta_2a1a1.at(1).Add(&jeta_2a1a1QCDMC.at(i),-norm);
	// jphi_2a1a1.at(1).Add(&jphi_2a1a1QCDMC.at(i),-norm);
	jpt_1a1a1.at(1).Add(&jpt_1a1a1QCDMC.at(i),-norm);
	// jeta_1a1a1.at(1).Add(&jeta_1a1a1QCDMC.at(i),-norm);
	// jphi_1a1a1.at(1).Add(&jphi_1a1a1QCDMC.at(i),-norm);
	cout<<"Soustraction du QCD: "<<endl;
	//JetFakesBDTScorea1a1.at(1).Fill(20);
	//JetFakesBDTScorea1a1.at(i).Fill(90);
	//cout<<JetFakesBDTScorea1a1QCDMC.at(1).Integral()<<endl;
	cout<<JetFakesBDTScorea1a1QCDMC.at(i).Integral()*norm<<endl;
	cout<<HiggsBDTScorea1a1QCDMC.at(i).Integral()*norm<<endl;
	cout<<ZTTBDTScorea1a1QCDMC.at(i).Integral()*norm<<endl;
	HiggsBDTScorea1a1.at(1).Add(&HiggsBDTScorea1a1QCDMC.at(i),-norm);
	JetFakesBDTScorea1a1.at(1).Add(&JetFakesBDTScorea1a1QCDMC.at(i),-norm);
	ZTTBDTScorea1a1.at(1).Add(&ZTTBDTScorea1a1QCDMC.at(i),-norm);


	polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolledQCDMC.at(i),-norm);
	//polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolledQCDMC.at(i),-norm);
	//polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolledQCDMC.at(i),-norm);
	//polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolledQCDMC.at(i),-norm);

	polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolledQCDMC.at(i),-norm);
	//polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolledQCDMC.at(i),-norm);
	//polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolledQCDMC.at(i),-norm);
	//polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolledQCDMC.at(i),-norm);

	polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolledQCDMC.at(i),-norm);
	//polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolledQCDMC.at(i),-norm);
	//polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolledQCDMC.at(i),-norm);
	//polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolledQCDMC.at(i),-norm);

	jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggs.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggs.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggs.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggs.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggs.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggs.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggs.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggs.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggs.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggs.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggs.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggs.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC.at(i),-norm);
	// jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSHiggs.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSHiggsQCDMC.at(i),-norm);
	// jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSHiggs.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSHiggsQCDMC.at(i),-norm);
	// jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSHiggs.at(1).Add(&jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSHiggsQCDMC.at(i),-norm);
	// jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSHiggs.at(1).Add(&jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSHiggsQCDMC.at(i),-norm);
	jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSHiggs.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsQCDMC.at(i),-norm*1.1);
	jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSHiggs.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsQCDMC.at(i),-norm*0.9);

	jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakes.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakes.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakes.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakes.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakes.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakes.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakes.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakes.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakes.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakes.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakes.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakes.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC.at(i),-norm);
	// jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSJetFakes.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSJetFakesQCDMC.at(i),-norm);
	// jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSJetFakes.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSJetFakesQCDMC.at(i),-norm);
	// jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSJetFakes.at(1).Add(&jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSJetFakesQCDMC.at(i),-norm);
	// jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSJetFakes.at(1).Add(&jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSJetFakesQCDMC.at(i),-norm);
	jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSJetFakes.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesQCDMC.at(i),-norm*1.1);
	jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSJetFakes.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesQCDMC.at(i),-norm*0.9);

	jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTT.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTT.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTT.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTT.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTT.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTT.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTT.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTT.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTT.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTT.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTT.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTT.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC.at(i),-norm);
	// jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSZTT.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSZTTQCDMC.at(i),-norm);
	// jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSZTT.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSZTTQCDMC.at(i),-norm);
	// jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSZTT.at(1).Add(&jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSZTTQCDMC.at(i),-norm);
	// jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSZTT.at(1).Add(&jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSZTTQCDMC.at(i),-norm);
	jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSZTT.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTQCDMC.at(i),-norm*1.1);
	jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSZTT.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTQCDMC.at(i),-norm*0.9);
	
	//std::cout <<"soustraction du QCD -DP" << std::endl;
  
	//----------------DP---------------------
	
	HiggsBDTScorea1a1_DP.at(1).Add(&HiggsBDTScorea1a1QCDMC_DP.at(i),-norm);
	JetFakesBDTScorea1a1_DP.at(1).Add(&JetFakesBDTScorea1a1QCDMC_DP.at(i),-norm);
	ZTTBDTScorea1a1_DP.at(1).Add(&ZTTBDTScorea1a1QCDMC_DP.at(i),-norm);

	polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled_DP.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolledQCDMC_DP.at(i),-norm);
	//polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled_DP.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolledQCDMC_DP.at(i),-norm);
	//polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled_DP.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolledQCDMC_DP.at(i),-norm);
	//polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolled_DP.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsUnrolledQCDMC_DP.at(i),-norm);

	polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled_DP.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolledQCDMC_DP.at(i),-norm);
	//polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled_DP.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolledQCDMC_DP.at(i),-norm);
	//polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled_DP.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolledQCDMC_DP.at(i),-norm);
	//polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolled_DP.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesUnrolledQCDMC_DP.at(i),-norm);

	polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled_DP.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolledQCDMC_DP.at(i),-norm);
	//polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled_DP.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolledQCDMC_DP.at(i),-norm);
	//polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled_DP.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolledQCDMC_DP.at(i),-norm);
	//polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolled_DP.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTUnrolledQCDMC_DP.at(i),-norm);

	jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggs_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggs_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggs_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggs_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggs_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggs_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggs_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggs_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggs_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggs_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggs_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSHiggsQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggs_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSHiggsQCDMC_DP.at(i),-norm);
	// jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSHiggs_DP.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSHiggsQCDMC_DP.at(i),-norm);
	// jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSHiggs_DP.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSHiggsQCDMC_DP.at(i),-norm);
	// jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSHiggs_DP.at(1).Add(&jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSHiggsQCDMC_DP.at(i),-norm);
	// jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSHiggs_DP.at(1).Add(&jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSHiggsQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSHiggs_DP.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsQCDMC_DP.at(i),-norm*1.1);
	jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSHiggs_DP.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMHiggsQCDMC_DP.at(i),-norm*0.9);

	jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakes_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakes_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakes_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakes_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakes_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakes_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakes_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakes_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakes_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakes_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakes_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakes_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(i),-norm);
	// jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSJetFakes_DP.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSJetFakesQCDMC_DP.at(i),-norm);
	// jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSJetFakes_DP.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSJetFakesQCDMC_DP.at(i),-norm);
	// jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSJetFakes_DP.at(1).Add(&jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSJetFakesQCDMC_DP.at(i),-norm);
	// jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSJetFakes_DP.at(1).Add(&jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSJetFakesQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSJetFakes_DP.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesQCDMC_DP.at(i),-norm*1.1);
	jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSJetFakes_DP.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMJetFakesQCDMC_DP.at(i),-norm*0.9);

	jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTT_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTT_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTT_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTT_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTT_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTT_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc1_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTT_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTT_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets0_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTT_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTT_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets1_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTT_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10UpPVRefitWithTracksBSZTTQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTT_DP.at(1).Add(&jetFakes_ff_tt_qcd_stat_unc2_njets2_mvadm10DownPVRefitWithTracksBSZTTQCDMC_DP.at(i),-norm);
	// jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSZTT_DP.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_systUpPVRefitWithTracksBSZTTQCDMC_DP.at(i),-norm);
	// jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSZTT_DP.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_systDownPVRefitWithTracksBSZTTQCDMC_DP.at(i),-norm);
	// jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSZTT_DP.at(1).Add(&jetFakes_ff_tt_qcd_systUpPVRefitWithTracksBSZTTQCDMC_DP.at(i),-norm);
	// jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSZTT_DP.at(1).Add(&jetFakes_ff_tt_qcd_systDownPVRefitWithTracksBSZTTQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_sub_systUpPVRefitWithTracksBSZTT_DP.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTQCDMC_DP.at(i),-norm*1.1);
	jetFakes_ff_tt_sub_systDownPVRefitWithTracksBSZTT_DP.at(1).Add(&polarimetricAcopAnglePVRefitWithTracksBSMVADMZTTQCDMC_DP.at(i),-norm*0.9);

	CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSHiggs.at(i).Add(&ttbarcontaminationHiggs.at(i),0.1*norm);
	CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSJetFakes.at(i).Add(&ttbarcontaminationJetFakes.at(i),0.1*norm);
	CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSZTT.at(i).Add(&ttbarcontaminationZTT.at(i),0.1*norm);
	CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSHiggs.at(i).Add(&ttbarcontaminationHiggs.at(i),-0.1*norm);
	CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSJetFakes.at(i).Add(&ttbarcontaminationJetFakes.at(i),-0.1*norm);
	CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSZTT.at(i).Add(&ttbarcontaminationZTT.at(i),-0.1*norm);

	CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesHiggs.at(i).Add(&ttbarcontaminationWfakesHiggs.at(i),0.1*norm);
	CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesJetFakes.at(i).Add(&ttbarcontaminationWfakesJetFakes.at(i),0.1*norm);
	CMS_ttbar_embeded_13TeVUpPVRefitWithTracksBSWfakesZTT.at(i).Add(&ttbarcontaminationWfakesZTT.at(i),0.1*norm);
	CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesHiggs.at(i).Add(&ttbarcontaminationWfakesHiggs.at(i),-0.1*norm);
	CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesJetFakes.at(i).Add(&ttbarcontaminationWfakesJetFakes.at(i),-0.1*norm);
	CMS_ttbar_embeded_13TeVDownPVRefitWithTracksBSWfakesZTT.at(i).Add(&ttbarcontaminationWfakesZTT.at(i),-0.1*norm);


	
	jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSHiggs.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSHiggsQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSHiggs.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSHiggsQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSHiggs.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSHiggsQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSHiggs.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSHiggsQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSHiggs.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSHiggsQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSHiggs.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSHiggsQCDMC.at(i),-norm);

	jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSJetFakes.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSJetFakesQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSJetFakes.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSJetFakesQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSJetFakes.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSJetFakesQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSJetFakes.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSJetFakesQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSJetFakes.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSJetFakesQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSJetFakes.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSJetFakesQCDMC.at(i),-norm);

	jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSZTT.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSZTTQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSZTT.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSZTTQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSZTT.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSZTTQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSZTT.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSZTTQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSZTT.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSZTTQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSZTT.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSZTTQCDMC.at(i),-norm);

	
	jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSHiggs_DP.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSHiggsQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSHiggs_DP.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSHiggsQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSHiggs_DP.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSHiggsQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSHiggs_DP.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSHiggsQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSHiggs_DP.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSHiggsQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSHiggs_DP.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSHiggsQCDMC_DP.at(i),-norm);

	jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSJetFakes_DP.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSJetFakes_DP.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSJetFakes_DP.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSJetFakes_DP.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSJetFakes_DP.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSJetFakes_DP.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSJetFakesQCDMC_DP.at(i),-norm);

	jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSZTT_DP.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets0UpPVRefitWithTracksBSZTTQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSZTT_DP.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets0UpPVRefitWithTracksBSZTTQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSZTT_DP.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets1UpPVRefitWithTracksBSZTTQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSZTT_DP.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets1UpPVRefitWithTracksBSZTTQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSZTT_DP.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets2UpPVRefitWithTracksBSZTTQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSZTT_DP.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets2UpPVRefitWithTracksBSZTTQCDMC_DP.at(i),-norm);

	

	jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSHiggs.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSHiggsQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSHiggs.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSHiggsQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSHiggs.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSHiggsQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSHiggs.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSHiggsQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSHiggs.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSHiggsQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSHiggs.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSHiggsQCDMC.at(i),-norm);

	jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSJetFakes.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSJetFakesQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSJetFakes.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSJetFakesQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSJetFakes.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSJetFakesQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSJetFakes.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSJetFakesQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSJetFakes.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSJetFakesQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSJetFakes.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSJetFakesQCDMC.at(i),-norm);

	jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSZTT.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSZTTQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSZTT.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSZTTQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSZTT.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSZTTQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSZTT.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSZTTQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSZTT.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSZTTQCDMC.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSZTT.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSZTTQCDMC.at(i),-norm);

	
	jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSHiggs_DP.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSHiggsQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSHiggs_DP.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSHiggsQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSHiggs_DP.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSHiggsQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSHiggs_DP.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSHiggsQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSHiggs_DP.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSHiggsQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSHiggs_DP.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSHiggsQCDMC_DP.at(i),-norm);

	jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSJetFakes_DP.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSJetFakes_DP.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSJetFakes_DP.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSJetFakes_DP.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSJetFakes_DP.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSJetFakes_DP.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSJetFakesQCDMC_DP.at(i),-norm);

	jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSZTT_DP.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets0DownPVRefitWithTracksBSZTTQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSZTT_DP.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets0DownPVRefitWithTracksBSZTTQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSZTT_DP.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets1DownPVRefitWithTracksBSZTTQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSZTT_DP.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets1DownPVRefitWithTracksBSZTTQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSZTT_DP.at(1).Add(&jetFakes_ff_tt_qcd_met_closure_syst_njets2DownPVRefitWithTracksBSZTTQCDMC_DP.at(i),-norm);
	jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSZTT_DP.at(1).Add(&jetFakes_ff_tt_qcd_syst_njets2DownPVRefitWithTracksBSZTTQCDMC_DP.at(i),-norm);

	std::cout << "end " << std::endl;
      }
    }
  }
  Selection::Finish();
}
