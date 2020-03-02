#include "HTauTau.h"
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



HTauTau::HTauTau(TString Name_, TString id_):
  Selection(Name_,id_),
  DataMC_Corr(true,true,false),
  tauTrgSF("tight")
{
  ChargeSumDummy = -999;
  selMuon_IsoDummy = 999.;
}

HTauTau::~HTauTau(){
  for(unsigned int j=0; j<Npassed.size(); j++){
    Logger(Logger::Info) << "Selection Summary before: "
			 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
			 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  HTauTau::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    //if(i==Trigger)             cut.at(Trigger)=1;
    if(i==Id_and_Kin)          cut.at(Id_and_Kin)=1;
    //if(i==NPairsFound)         cut.at(NPairsFound)=1;
    //if(i==Tau1Isolation)       cut.at(Tau1Isolation)=1.;
    //if(i==Tau2Isolation)       cut.at(Tau2Isolation)=1.;
    //if(i==LeptonVeto)          cut.at(LeptonVeto)=0.;
    if(i==PairCharge)          cut.at(PairCharge)=1.;
    //if(i==PairMass)            cut.at(PairMass)=1000.;
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
    //   title.at(i)="At least 1 good pair with Trig+Matching";
    //   hlabel="At least 1 good pair with Trig+Matching";
    //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Trigger_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Trigger_",htitle,2,-0.5,1.5,hlabel,"Events"));
    // }
    /* else */if(i==Id_and_Kin){
      title.at(i)="Id and Kinematic";
      hlabel="Number of Event with good particles";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Id_and_Kin_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Id_and_Kin_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    // else if(i==NPairsFound){
    //   title.at(i)="Pairs with good DeltaR";
    //   hlabel="Pairs with good DeltaR";
    //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NPairsFound_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NPairsFound_",htitle,2,-0.5,1.5,hlabel,"Events"));
    // }
    // else if(i==Tau1Isolation){
    //   title.at(i)="Tau1 Isolation";
    //   hlabel="Isolation of Tau1";
    //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Tau1Isolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Tau1Isolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
    // }
    // else if(i==Tau2Isolation){
    //   title.at(i)="Tau2 Isolation";
    //   hlabel="Isolation of Tau2";
    //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Tau2Isolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Tau2Isolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
    // }
    // else if(i==LeptonVeto){
    //   title.at(i)="Lepton Veto";
    //   hlabel="Third Lepton Veto  ";
    //   Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_LeptonVeto_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //   Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_LeptonVeto_",htitle,2,-0.5,1.5,hlabel,"Events"));
    //}
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
  
  Tau1PT=HConfig.GetTH1D(Name+"_Tau1PT","Transverse momentum of selected #tau1 candidate",20,35,80," P_{T}(#tau1), GeV","Events");
  Tau1E=HConfig.GetTH1D(Name+"_Tau1E","Energy of selected #tau1 candidate",20,35.,140," E(#tau1), GeV","Events");
  Tau1Mass=HConfig.GetTH1D(Name+"_Tau1Mass","Mass of selected #tau1 candidate",20,0,2.," M(#tau1), GeV","Events");
  Tau1Phi=HConfig.GetTH1D(Name+"_Tau1Phi","Phi of selected #tau1 candidate",10,-3.5,3.5," #phi(#tau1)","Events");
  Tau1Eta=HConfig.GetTH1D(Name+"_Tau1Eta","Pseudorapidity tau1",15,-2.7,2.7," #eta(#tau1)","Events");
  Tau1dz=HConfig.GetTH1D(Name+"_Tau1dz","Tau1 dz",20,-0.04,0.04,"Tau1 dz","Events");
  Tau1HPSDecayMode=HConfig.GetTH1D(Name+"_Tau1HPSDecayMode","Decay mode of the selected #tau candidate",11,-0.5,10.5," HPS Mode ","Events");
  
  Tau2PT=HConfig.GetTH1D(Name+"_Tau2PT","Transverse momentum of selected #tau2 candidate",20,30,55," P_{T}(#tau2), GeV","Events");
  Tau2E=HConfig.GetTH1D(Name+"_Tau2E","Energy of selected #tau2 candidate",20,30.,140," E(#tau2), GeV","Events");
  Tau2Mass=HConfig.GetTH1D(Name+"_Tau2Mass","Mass of selected #tau2 candidate",20,0,2.," M(#tau2), GeV","Events");
  Tau2Phi=HConfig.GetTH1D(Name+"_Tau2Phi","Phi of selected #tau2 candidate",10,-3.5,3.5," #phi(#tau2)","Events");
  Tau2Eta=HConfig.GetTH1D(Name+"_Tau2Eta","Pseudorapidity Tau2",15,-2.7,2.7," #eta(#tau2)","Events");
  Tau2dz=HConfig.GetTH1D(Name+"_Tau2dz","Tau2dz",20,-0.04,0.04,"Tau2 dz","Events");
  Tau2HPSDecayMode=HConfig.GetTH1D(Name+"_Tau2HPSDecayMode","Decay mode of the selected #tau candidate",11,-0.5,10.5," HPS Mode ","Events");
  /*
    Tau1isolation=HConfig.GetTH1D(Name+"_Tau1isolation","First Tau isolation 1- Loose, 2- Medium, 3 Tight, 4-VTight",4,0.5,4.5," Discriminator","Events");
    Tau2isolation=HConfig.GetTH1D(Name+"_Tau2isolation","First Tau isolation 1- Loose, 2- Medium, 3 Tight, 4-VTight",4,0.5,4.5," Discriminator","Events");
  
  
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

  TauTauVisMass=HConfig.GetTH1D(Name+"_TauTauVisMass","Visible invariant mass of a tau pair",40,0,200," M(#tau#tau)_{vis}, GeV","Events");
  TauTauTruthMass=HConfig.GetTH1D(Name+"_TauTauTruthMass","Truth invariant mass of a tau pair",40,0,150," M(#tau#tau)_{truth}, GeV","Events");
  //TauTauFullMass=HConfig.GetTH1D(Name+"_TauTauFullMass","Full invariant mass of a tau pair",40,0,150," M(#tau#tau)_{full}, GeV","Events");

  NQCD=HConfig.GetTH1D(Name+"_NQCD","NQCD",4,0.5,4.5,"NQCD in ABCD","Events");
  // TauTauFullMass_B=HConfig.GetTH1D(Name+"_TauTauFullMass_B","TauTauFullMass_B",40,0,150,"#tau_h#tau_h SVFit Mass in B","Events");
  // TauTauFullMass_C=HConfig.GetTH1D(Name+"_TauTauFullMass_C","TauTauFullMass_C",40,0,150,"#tau_h#tau_h SVFit Mass in C","Events");
  // TauTauFullMass_D=HConfig.GetTH1D(Name+"_TauTauFullMass_D","TauTauFullMass_D",40,0,150,"#tau_h#tau_h SVFit Mass in D","Events");

  dRTauTau=HConfig.GetTH1D(Name+"_dRTauTau","#Delta R",20,0.,3.5," #Delta R","Events");

  MET=HConfig.GetTH1D(Name+"_MET","MET",20,0,80,"MET, GeV","Events");
  METphi=HConfig.GetTH1D(Name+"_METphi","METphi",10,-3.5,3.5,"METphi","Events");
  PUPPImet=HConfig.GetTH1D(Name+"_PUPPImet","PUPPImet",10,0,75,"PUPPImet, GeV","Events");
  PUPPImetphi=HConfig.GetTH1D(Name+"_PUPPImetphi","PUPPImetphi",10,-3.5,3.5,"PUPPImetphi","Events");
  TransverseMass=HConfig.GetTH1D(Name+"_TransverseMass","TransverseMass, GeV",40,0,110,"TransverseMass","Events");
  
  NPrimeVtx=HConfig.GetTH1D(Name+"_NPrimeVtx","NPrimeVtx",10,0,50,"N vtx","Events");
  NPU=HConfig.GetTH1D(Name+"_npu","npu",10,0,50,"N pu","Events");
  RHO=HConfig.GetTH1D(Name+"_rho","rho",10,0,30,"rho","Events");
  
  NbJets=HConfig.GetTH1D(Name+"_NbJets","NbJets",12,0,12,"Number of jets","Events");
  
  h_SVFitMass = HConfig.GetTH1D(Name+"_SVFitMass","SVFitMass",100,0.,200.,"m_{SVfit}(#tau_{h},#tau_{h})/GeV");
  h_SVFitStatus = HConfig.GetTH1D(Name+"_SVFitStatus", "SVFitStatus", 5, -0.5, 4.5, "Status of SVFit calculation");

  // svfTau1E = HConfig.GetTH1D(Name+"_svfTau1E","svFitTau1E",40,20.,120.,"E_{SVfit}(#tau_{h}1)/GeV");
  // svfTau2E = HConfig.GetTH1D(Name+"_svfTau2E","svFitTau2E",40,20.,120.,"E_{SVfit}(#tau_{h}2)/GeV");

  
  // PhiDatasvfitpipi = HConfig.GetTH1D(Name+"_PhiDatasvfitpipi","Phisvfit #pi#pi",30,-TMath::Pi(),TMath::Pi(),"Angle with svfit between Tau+ and an initial proton in the new xy plan for #pi#pi channel");
  // PhiDatasvfitpirho = HConfig.GetTH1D(Name+"_PhiDatasvfitpirho","Phisvfit #pi#rho",30,-TMath::Pi(),TMath::Pi(),"Angle with svfit between Tau+ and an initial proton in the new xy plan for #pi#rho channel");
  // PhiDatasvfitpia1 = HConfig.GetTH1D(Name+"_PhiDatasvfitpia1","Phisvfit a1#pi",30,-TMath::Pi(),TMath::Pi(),"Angle with svfit between Tau+ and an initial proton in the new xy plan for a1#pi channel");
  // PhiDatasvfitrhorho = HConfig.GetTH1D(Name+"_PhiDatasvfitrhorho","Phisvfit #rho#rho",30,-TMath::Pi(),TMath::Pi(),"Angle with svfit between Tau+ and an initial proton in the new xy plan for #rho#rho channel");
  // PhiDatasvfitrhoa1 = HConfig.GetTH1D(Name+"_PhiDatasvfitrhoa1","Phisvfit a1#rho",30,-TMath::Pi(),TMath::Pi(),"Angle with svfit between Tau+ and an initial proton in the new xy plan for a1#rho channel");
  // PhiDatasvfita1a1 = HConfig.GetTH1D(Name+"_PhiDatasvfita1a1","Phisvfit a1a1",30,-TMath::Pi(),TMath::Pi(),"Angle with svfit between Tau+ and an initial proton in the new xy plan for a1a1 channel");
  /*
    PhiDatavispipi = HConfig.GetTH1D(Name+"_PhiDatavispipi","Phivis #pi#pi",30,-TMath::Pi(),TMath::Pi(),"Angle with vis between Tau+ and an initial proton in the new xy plan for #pi#pi channel");
    PhiDatavispirho = HConfig.GetTH1D(Name+"_PhiDatavispirho","Phivis #pi#rho",30,-TMath::Pi(),TMath::Pi(),"Angle with vis between Tau+ and an initial proton in the new xy plan for #pi#rho channel");
    PhiDatavispia1 = HConfig.GetTH1D(Name+"_PhiDatavispia1","Phivis a1#pi",30,-TMath::Pi(),TMath::Pi(),"Angle with vis between Tau+ and an initial proton in the new xy plan for a1#pi channel");
    PhiDatavisrhorho = HConfig.GetTH1D(Name+"_PhiDatavisrhorho","Phivis #rho#rho",30,-TMath::Pi(),TMath::Pi(),"Angle with vis between Tau+ and an initial proton in the new xy plan for #rho#rho channel");
    PhiDatavisrhoa1 = HConfig.GetTH1D(Name+"_PhiDatavisrhoa1","Phivis a1#rho",30,-TMath::Pi(),TMath::Pi(),"Angle with vis between Tau+ and an initial proton in the new xy plan for #rho#pi channel");
    PhiDatavisa1a1 = HConfig.GetTH1D(Name+"_PhiDatavisa1a1","Phivis a1a1",30,-TMath::Pi(),TMath::Pi(),"Angle with vis between Tau+ and an initial proton in the new xy plan for a1a1 channel");
  */
  
  // Etasvfit = HConfig.GetTH1D(Name+"_Etasvfit", "Etasvfit", 30, -10, 10, "Eta  with svfit between Tau+ and Tau- in the new xy plan");

  // Phisvfitpipi = HConfig.GetTH1D(Name+"_Phisvfitpipi","Phisvfit #pi#pi",30,-TMath::Pi(),TMath::Pi(),"Angle with svfit between Tau+ and an initial proton in the new xy plan for #pi#pi channel");
  // Phisvfitpirho = HConfig.GetTH1D(Name+"_Phisvfitpirho","Phisvfit #pi#rho",30,-TMath::Pi(),TMath::Pi(),"Angle with svfit between Tau+ and an initial proton in the new xy plan for #pi#rho channel");
  // // Phisvfitlpi = HConfig.GetTH1D(Name+"_Phisvfitlpi","Phisvfit l#pi",30,-TMath::Pi(),TMath::Pi(),"Angle with svfit between Tau+ and an initial proton in the new xy plan for #pi channel");
  // // Phisvfitlrho = HConfig.GetTH1D(Name+"_Phisvfitlrho","Phisvfit l#rho",30,-TMath::Pi(),TMath::Pi(),"Angle with svfit between Tau+ and an initial proton in the new xy plan for l#rho channel");
  // Phisvfitpia1 = HConfig.GetTH1D(Name+"_Phisvfitpia1","Phisvfit a1#pi",30,-TMath::Pi(),TMath::Pi(),"Angle with svfit between Tau+ and an initial proton in the new xy plan for a1#pi channel");
  //  Phisvfitrhorho = HConfig.GetTH1D(Name+"_Phisvfitrhorho","Phisvfit #rho#rho",30,-TMath::Pi(),TMath::Pi(),"Angle with svfit between Tau+ and an initial proton in the new xy plan for #rho#rho channel");
  // Phisvfitrhoa1 = HConfig.GetTH1D(Name+"_Phisvfitrhoa1","Phisvfit a1#rho",30,-TMath::Pi(),TMath::Pi(),"Angle with svfit between Tau+ and an initial proton in the new xy plan for a1#rho channel");
  // // Phisvfitla1 = HConfig.GetTH1D(Name+"_Phisvfitla1","Phisvfit la1",30,-TMath::Pi(),TMath::Pi(),"Angle with svfit between Tau+ and an initial proton in the new xy plan for la1 channel");
  // Phisvfita1a1 = HConfig.GetTH1D(Name+"_Phisvfita1a1","Phisvfit a1a1",30,-TMath::Pi(),TMath::Pi(),"Angle with svfit between Tau+ and an initial proton in the new xy plan for a1a1 channel");

  // Thetasvfit = HConfig.GetTH1D(Name+"_Thetasvfit","Thetasvfit",30,0.,3.,"Original theta of Tau- with svfit");
  

  Etavis = HConfig.GetTH1D(Name+"_Etavis", "Etavis", 30, -10, 10, "Eta between Tau+ and Tau- in the new xy plan with visible particles");

  Phivispipi = HConfig.GetTH1D(Name+"_Phivispipi","Phivis #pi#pi",30,-TMath::Pi(),TMath::Pi(),"Visible angle between Tau+ and an initial proton in the new xy plan for #pi#pi channel");
  Phivispirho = HConfig.GetTH1D(Name+"_Phivispirho","Phivis #pi#rho",30,-TMath::Pi(),TMath::Pi(),"Visible angle between Tau+ and an initial proton in the new xy plan for #pi#rho channel");
  // Phivislpi = HConfig.GetTH1D(Name+"_Phivislpi","Phivis l#pi",30,-TMath::Pi(),TMath::Pi(),"Visible angle between Tau+ and an initial proton in the new xy plan for #pi channel");
  // Phivislrho = HConfig.GetTH1D(Name+"_Phivislrho","Phivis l#rho",30,-TMath::Pi(),TMath::Pi(),"Visible angle between Tau+ and an initial proton in the new xy plan for l#rho channel");
  Phivispia1 = HConfig.GetTH1D(Name+"_Phivispia1","Phivis a1#pi",30,-TMath::Pi(),TMath::Pi(),"Visible angle between Tau+ and an initial proton in the new xy plan for a1#pi channel");
  Phivisrhorho = HConfig.GetTH1D(Name+"_Phivisrhorho","Phivis #rho#rho",30,-TMath::Pi(),TMath::Pi(),"Visible angle between Tau+ and an initial proton in the new xy plan for #rho#rho channel");
  Phivisrhoa1 = HConfig.GetTH1D(Name+"_Phivisrhoa1","Phivis a1#rho",30,-TMath::Pi(),TMath::Pi(),"Visible angle between Tau+ and an initial proton in the new xy plan for a1#rho channel");
  // Phivisla1 = HConfig.GetTH1D(Name+"_Phivisla1","Phivis la1",30,-TMath::Pi(),TMath::Pi(),"Visible angle between Tau+ and an initial proton in the new xy plan for la1 channel");
  Phivisa1a1 = HConfig.GetTH1D(Name+"_Phivisa1a1","Phivis a1a1",30,-TMath::Pi(),TMath::Pi(),"Visible angle between Tau+ and an initial proton in the new xy plan for a1a1 channel");

  Thetavis = HConfig.GetTH1D(Name+"_Thetavis","Thetavis",30,0.,3.,"Original theta of Tau- with visible particles");
  
  
  Etatruth = HConfig.GetTH1D(Name+"_Etatruth", "Etatruth", 30, -10, 10, "Real Eta between Tau+ and Tau- in the new xy plan");

  Phitruthpipi = HConfig.GetTH1D(Name+"_Phitruthpipi","Phitruth #pi#pi",30,-TMath::Pi(),TMath::Pi(),"Real angle between Tau+ and an initial proton in the new xy plan for #pi#pi channel");
  Phitruthpirho = HConfig.GetTH1D(Name+"_Phitruthpirho","Phitruth #pi#rho",30,-TMath::Pi(),TMath::Pi(),"Real angle between Tau+ and an initial proton in the new xy plan for #pi#rho channel");
  //Phitruthlpi = HConfig.GetTH1D(Name+"_Phitruthlpi","Phitruth l#pi",30,-TMath::Pi(),TMath::Pi(),"Real angle between Tau+ and an initial proton in the new xy plan for #pi channel");
  //Phitruthlrho = HConfig.GetTH1D(Name+"_Phitruthlrho","Phitruth l#rho",30,-TMath::Pi(),TMath::Pi(),"Real angle between Tau+ and an initial proton in the new xy plan for l#rho channel");
  Phitruthpia1 = HConfig.GetTH1D(Name+"_Phitruthpia1","Phitruth a1#pi",30,-TMath::Pi(),TMath::Pi(),"Real angle between Tau+ and an initial proton in the new xy plan for a1#pi channel");
  Phitruthrhorho = HConfig.GetTH1D(Name+"_Phitruthrhorho","Phitruth #rho#rho",30,-TMath::Pi(),TMath::Pi(),"Real angle between Tau+ and an initial proton in the new xy plan for #rho#rho channel");
  Phitruthrhoa1 = HConfig.GetTH1D(Name+"_Phitruthrhoa1","Phitruth a1#rho",30,-TMath::Pi(),TMath::Pi(),"Real angle between Tau+ and an initial proton in the new xy plan for a1#rho channel");
  // Phitruthla1 = HConfig.GetTH1D(Name+"_Phitruthla1","Phitruth la1",30,-TMath::Pi(),TMath::Pi(),"Real angle between Tau+ and an initial proton in the new xy plan for la1 channel");
  Phitrutha1a1 = HConfig.GetTH1D(Name+"_Phitrutha1a1","Phitruth a1a1",30,-TMath::Pi(),TMath::Pi(),"Real angle between Tau+ and an initial proton in the new xy plan for a1a1 channel");

  Thetatruth = HConfig.GetTH1D(Name+"_Thetatruth","Thetatruth",30,0.,3.,"Real original theta of Tau-");
  
  /*
    PhiSvFitRespipi=HConfig.GetTH1D(Name+"_PhiSvFitRespipi","Resolution of the full Phi",100,-TMath::Pi(),TMath::Pi(),"Full Phi resol for #pi#pi channel, GeV","Events");
    PhiSvFitRespirho=HConfig.GetTH1D(Name+"_PhiSvFitRespirho","Resolution of the full Phi",100,-TMath::Pi(),TMath::Pi(),"Full Phi resol for #pi#rho channel, GeV","Events");
    PhiSvFitReslpi=HConfig.GetTH1D(Name+"_PhiSvFitReslpi","Resolution of the full Phi",100,-TMath::Pi(),TMath::Pi(),"Full Phi resol for l#pi channel, GeV","Events");
    PhiSvFitRespia1=HConfig.GetTH1D(Name+"_PhiSvFitRespia1","Resolution of the full Phi",100,-TMath::Pi(),TMath::Pi(),"Full Phi resol for a1#pi channel, GeV","Events");
    PhiSvFitResrhorho=HConfig.GetTH1D(Name+"_PhiSvFitResrhorho","Resolution of the full Phi",100,-TMath::Pi(),TMath::Pi(),"Full Phi resol for #rho#rho channel, GeV","Events");
    PhiSvFitResrhoa1=HConfig.GetTH1D(Name+"_PhiSvFitResrhoa1","Resolution of the full Phi",100,-TMath::Pi(),TMath::Pi(),"Full Phi resol for a1#rho channel, GeV","Events");
    PhiSvFitReslrho=HConfig.GetTH1D(Name+"_PhiSvFitReslrho","Resolution of the full Phi",100,-TMath::Pi(),TMath::Pi(),"Full Phi resol for l#rho channel, GeV","Events");
    PhiSvFitResla1=HConfig.GetTH1D(Name+"_PhiSvFitResla1","Resolution of the full Phi",100,-TMath::Pi(),TMath::Pi(),"Full Phi resol for la1 channel, GeV","Events");
    PhiSvFitResa1a1=HConfig.GetTH1D(Name+"_PhiSvFitResa1a1","Resolution of the full Phi",100,-TMath::Pi(),TMath::Pi(),"Full Phi resol for a1a1 channel, GeV","Events");
  
    PhiVisRespipi=HConfig.GetTH1D(Name+"_PhiVisRespipi","Resolution of the full Phi",100,-TMath::Pi(),TMath::Pi(),"Visible Phi resol for #pi#pi channel, GeV","Events");
    PhiVisRespirho=HConfig.GetTH1D(Name+"_PhiVisRespirho","Resolution of the full Phi",100,-TMath::Pi(),TMath::Pi(),"Visible Phi resol for #pi#rho channel, GeV","Events");
    PhiVisReslpi=HConfig.GetTH1D(Name+"_PhiVisReslpi","Resolution of the full Phi",100,-TMath::Pi(),TMath::Pi(),"Visible Phi resol for l#pi channel, GeV","Events");
    PhiVisRespia1=HConfig.GetTH1D(Name+"_PhiVisRespia1","Resolution of the full Phi",100,-TMath::Pi(),TMath::Pi(),"Visible Phi resol for a1#pi channel, GeV","Events");
    PhiVisResrhorho=HConfig.GetTH1D(Name+"_PhiVisResrhorho","Resolution of the full Phi",100,-TMath::Pi(),TMath::Pi(),"Visible Phi resol for #rho#rho channel, GeV","Events");
    PhiVisResrhoa1=HConfig.GetTH1D(Name+"_PhiVisResrhoa1","Resolution of the full Phi",100,-TMath::Pi(),TMath::Pi(),"Visible Phi resol for a1#rho channel, GeV","Events");
    PhiVisReslrho=HConfig.GetTH1D(Name+"_PhiVisReslrho","Resolution of the full Phi",100,-TMath::Pi(),TMath::Pi(),"Visible Phi resol for l#rho channel, GeV","Events");
    PhiVisResla1=HConfig.GetTH1D(Name+"_PhiVisResla1","Resolution of the full Phi",100,-TMath::Pi(),TMath::Pi(),"Visible Phi resol for la1 channel, GeV","Events");
    PhiVisResa1a1=HConfig.GetTH1D(Name+"_PhiVisResa1a1","Resolution of the full Phi",100,-TMath::Pi(),TMath::Pi(),"Visible Phi resol for a1a1 channel, GeV","Events");
  

    TauTauFullPtRes=HConfig.GetTH1D(Name+"_TauTauFullPtRes","Resolution of the full TauTau Pt",100,-50.,50.,"Full TauTau Pt resol, GeV","Events");
    TauTauFullEtaRes=HConfig.GetTH1D(Name+"_TauTauFullEtaRes","Resolution of the full TauTau Eta",100,-5.,5.,"Full TauTau Eta resol, GeV","Events");
    TauTauFullPhiRes=HConfig.GetTH1D(Name+"_TauTauFullPhiRes","Resolution of the full TauTau Phi",100,-3.,3.,"Full TauTau Phi resol, GeV","Events");

    TauplusFullPtRes=HConfig.GetTH1D(Name+"_TauplusFullPtRes","Resolution of the full Tau+ Pt",100,-50.,50.,"Full Tau+ Pt resol, GeV","Events");
    TauplusFullEtaRes=HConfig.GetTH1D(Name+"_TauplusFullEtaRes","Resolution of the full Tau+ Eta",100,-5.,5.,"Full Tau+ Eta resol, GeV","Events");
    TauplusFullPhiRes=HConfig.GetTH1D(Name+"_TauplusFullPhiRes","Resolution of the full Tau+ Phi",100,-3.,3.,"Full Tau+ Phi resol, GeV","Events");

    TauminusFullPtRes=HConfig.GetTH1D(Name+"_TauminusFullPtRes","Resolution of the full Tau- Pt",100,-50.,50.,"Full Tau- Pt resol, GeV","Events");
    TauminusFullEtaRes=HConfig.GetTH1D(Name+"_TauminusFullEtaRes","Resolution of the full Tau- Eta",50,-2.,2.,"Full Tau- Eta resol, GeV","Events");
    TauminusFullPhiRes=HConfig.GetTH1D(Name+"_TauminusFullPhiRes","Resolution of the full Tau- Phi",100,-3.,3.,"Full Tau- Phi resol, GeV","Events");
 
    TauTauVisPtRes=HConfig.GetTH1D(Name+"_TauTauVisPtRes","Resolution of the visible TauTau Pt",100,-50.,50.,"Visible TauTau Pt resol, GeV","Events");
    TauTauVisEtaRes=HConfig.GetTH1D(Name+"_TauTauVisEtaRes","Resolution of the visible TauTau Eta",100,-5.,5.,"Visible TauTau Eta resol, GeV","Events");
    TauTauVisPhiRes=HConfig.GetTH1D(Name+"_TauTauVisPhiRes","Resolution of the visible TauTau Phi",100,-3.,3.,"Visible TauTau Phi resol, GeV","Events");

    TauplusVisPtRes=HConfig.GetTH1D(Name+"_TauplusVisPtRes","Resolution of the visible Tau+ Pt",100,-50.,50.,"Visible Tau+ Pt resol, GeV","Events");
    TauplusVisEtaRes=HConfig.GetTH1D(Name+"_TauplusVisEtaRes","Resolution of the visible Tau+ Eta",100,-5.,5.,"Visible Tau+ Eta resol, GeV","Events");
    TauplusVisPhiRes=HConfig.GetTH1D(Name+"_TauplusVisPhiRes","Resolution of the visible Tau+ Phi",100,-3.,3.,"Visible Tau+ Phi resol, GeV","Events");
    
    TauminusVisPtRes=HConfig.GetTH1D(Name+"_TauminusVisPtRes","Resolution of the visible Tau- Pt",100,-50.,50.,"Visible Tau- Pt resol, GeV","Events");
    TauminusVisEtaRes=HConfig.GetTH1D(Name+"_TauminusVisEtaRes","Resolution of the visible Tau- Eta",50,-2.,5.,"Visible Tau- Eta resol, GeV","Events");
    TauminusVisPhiRes=HConfig.GetTH1D(Name+"_TauminusVisPhiRes","Resolution of the visible Tau- Phi",100,-3.,3.,"Visible Tau- Phi resol, GeV","Events");



    TauTauFullPtResPull=HConfig.GetTH1D(Name+"_TauTauFullPtResPull","Pull Plot of the full TauTau Pt",50,-10.,10.,"Full TauTau Pt pull plot, GeV","Events");
    TauTauFullEtaResPull=HConfig.GetTH1D(Name+"_TauTauFullEtaResPull","Pull Plot of the full TauTau Eta",100,-5.,5.,"Full TauTau Eta pull plot, GeV","Events");
    TauTauFullPhiResPull=HConfig.GetTH1D(Name+"_TauTauFullPhiResPull","Pull Plot of the full TauTau Phi",100,-3.,3.,"Full TauTau Phi pull plot, GeV","Events");

    TauplusFullPtResPull=HConfig.GetTH1D(Name+"_TauplusFullPtResPull","Pull Plot of the full Tau+ Pt",50,-10.,10.,"Full Tau+ Pt pull plot, GeV","Events");
    TauplusFullEtaResPull=HConfig.GetTH1D(Name+"_TauplusFullEtaResPull","Pull Plot of the full Tau+ Eta",100,-5.,5.,"Full Tau+ Eta pull plot, GeV","Events");
    TauplusFullPhiResPull=HConfig.GetTH1D(Name+"_TauplusFullPhiResPull","Pull Plot of the full Tau+ Phi",100,-3.,3.,"Full Tau+ Phi pull plot, GeV","Events");

    TauminusFullPtResPull=HConfig.GetTH1D(Name+"_TauminusFullPtResPull","Pull Plot of the full Tau- Pt",50,-2.,2.,"Full Tau- Pt pull plot, GeV","Events");
    TauminusFullEtaResPull=HConfig.GetTH1D(Name+"_TauminusFullEtaResPull","Pull Plot of the full Tau- Eta",50,-2.,2.,"Full Tau- Eta pull plot, GeV","Events");
    TauminusFullPhiResPull=HConfig.GetTH1D(Name+"_TauminusFullPhiResPull","Pull Plot of the full Tau- Phi",100,-3.,3.,"Full Tau- Phi pull plot, GeV","Events");


    TauTauVisPtResPull=HConfig.GetTH1D(Name+"_TauTauVisPtResPull","Pull Plot of the visible TauTau Pt",50,-10.,10.,"Visible TauTau Pt pull plot, GeV","Events");
    TauTauVisEtaResPull=HConfig.GetTH1D(Name+"_TauTauVisEtaResPull","Pull Plot of the visible TauTau Eta",100,-5.,5.,"Visible TauTau Eta pull plot, GeV","Events");
    TauTauVisPhiResPull=HConfig.GetTH1D(Name+"_TauTauVisPhiResPull","Pull Plot of the visible TauTau Phi",100,-3.,3.,"Visible TauTau Phi pull plot, GeV","Events");

    TauplusVisPtResPull=HConfig.GetTH1D(Name+"_TauplusVisPtResPull","Pull Plot of the visible Tau+ Pt",50,-10.,10.,"Visible Tau+ Pt pull plot, GeV","Events");
    TauplusVisEtaResPull=HConfig.GetTH1D(Name+"_TauplusVisEtaResPull","Pull Plot of the visible Tau+ Eta",100,-5.,5.,"Visible Tau+ Eta pull plot, GeV","Events");
    TauplusVisPhiResPull=HConfig.GetTH1D(Name+"_TauplusVisPhiResPull","Pull Plot of the visible Tau+ Phi",100,-3.,3.,"Visible Tau+ Phi pull plot, GeV","Events");
    
    TauminusVisPtResPull=HConfig.GetTH1D(Name+"_TauminusVisPtResPull","Pull Plot of the visible Tau- Pt",50,-10.,10.,"Visible Tau- Pt pull plot, GeV","Events");
    TauminusVisEtaResPull=HConfig.GetTH1D(Name+"_TauminusVisEtaResPull","Pull Plot of the visible Tau- Eta",50,-2.,2.,"Visible Tau- Eta pull plot, GeV","Events");
    TauminusVisPhiResPull=HConfig.GetTH1D(Name+"_TauminusVisPhiResPull","Pull Plot of the visible Tau- Phi",100,-3.,3.,"Visible Tau- Phi pull plot, GeV","Events");
  */
  //DRTruth=HConfig.GetTH1D(Name+"_DRTruth","Delta R of Truth Tau",100,0.,5.,"Delta R of Full Tau","Events");
  // DRFull=HConfig.GetTH1D(Name+"_DRFull","Delta R of Full Tau",100,0.,5.,"Delta R of Full Tau","Events");
  //DRFullTruth=HConfig.GetTH1D(Name+"_DRTruth","Delta R of Truth Tau",100,0.,5.,"Delta R of Full Tau","Events");
  // DRVisTruth=HConfig.GetTH1D(Name+"_DRFull","Delta R of Full Tau",100,0.,5.,"Delta R of Full Tau","Events");
  
  Pi0EnergyRes=HConfig.GetTH1D(Name+"_Pi0EnergyRes","Energy resolution of Pi0",100,-50.,50.,"Energy resolution of Pi0, GeV","Events");
  Pi0EnergyResPull=HConfig.GetTH1D(Name+"_Pi0EnergyResPull","Energy Pull Plot of Pi0",100,-50.,50.,"Energy Pull Plot of Pi0, GeV","Events");
 

  NewPhivsDeltaPhi=HConfig.GetTH2D(Name+"_NewPhivsDeltaPhi","New #phi VS old #Delta#phi",30,-TMath::Pi(),TMath::Pi(),30,-TMath::Pi(),TMath::Pi(),"","Events");
  NewPhivsDeltaEta=HConfig.GetTH2D(Name+"_NewPhivsDeltaEta","New #phi VS old #Delta#eta",30,-TMath::Pi(),TMath::Pi(),20,-2.7,2.7,"","Events");
  NewPhivsPhiproton=HConfig.GetTH2D(Name+"_NewPhivsPhiproton","New #phi VS old #phi of the proton",30,-TMath::Pi(),TMath::Pi(),30,-TMath::Pi(),TMath::Pi(),"","Events");
  NewPhivsPhiTauplus=HConfig.GetTH2D(Name+"_NewPhivsPhiTauplus","New #phi VS old #phi of the Tau-",30,-TMath::Pi(),TMath::Pi(),30,-TMath::Pi(),TMath::Pi(),"","Events");
  NewPhivsEtaproton=HConfig.GetTH2D(Name+"_NewPhivsEtaproton","New #phi VS old #eta of the proton",30,-TMath::Pi(),TMath::Pi(),20,-2.7,2.7,"","Events");
  NewPhivsEtaTauplus=HConfig.GetTH2D(Name+"_NewPhivsEtaTauplus","New #phi VS old #eta of the Tau-",30,-TMath::Pi(),TMath::Pi(),20,-2.7,2.7,"","Events");
  NewPhivsZPt=HConfig.GetTH2D(Name+"_NewPhivsZPt","New #phi VS Pt_{Z}",30,-TMath::Pi(),TMath::Pi(),40,0,100,"","Events");

  NewPhiSignal=HConfig.GetTH1D(Name+"_NewPhiSignal","New Phi for all Signal",30,-TMath::Pi(),TMath::Pi(),"","Events");
  NewPhiQCD=HConfig.GetTH1D(Name+"_NewPhiQCD","New Phi for QCD",30,-TMath::Pi(),TMath::Pi(),"","Events");
																			      
  ZPtVis=HConfig.GetTH1D(Name+"_ZPtVis","Visible Pt_{Z}",40,0,100,"","Events");
  
  IstauminusvisPhysical=HConfig.GetTH1D(Name+"_IstauminusvisPhysical","IstauminusvisPhysical",2,-0.5,1.5,"","Events");
  IstauplusvisPhysical=HConfig.GetTH1D(Name+"_IstauplusvisPhysical","IstauplusvisPhysical",2,-0.5,1.5,"","Events");
  IsPairPhysical=HConfig.GetTH1D(Name+"_IsPairPhysical","IsPairPhysical",2,-0.5,1.5,"","Events");
  
  //ResolPullTauTauFroma1a1MeanEnergy=HConfig.GetTH1D(Name+"_ResolPullTauTauFroma1a1MeanEnergy","ResolPullTauTauFroma1a1MeanEnergy",30,-1,1,"","Events");
  //ResolPullTauTauFroma1a1MZEnergy=HConfig.GetTH1D(Name+"_ResolPullTauTauFroma1a1MZEnergy","ResolPullTauTauFroma1a1MZEnergy",30,-1,1,"","Events");
  ResolPullTauTauFroma1a1MeanMomentum=HConfig.GetTH1D(Name+"_ResolPullTauTauFroma1a1MeanMomentum","ResolPullTauTauFroma1a1MeanMomentum",30,-1,1,"","Events");
  ResolPullTauTauFroma1a1MZMomentum=HConfig.GetTH1D(Name+"_ResolPullTauTauFroma1a1MZMomentum","ResolPullTauTauFroma1a1MZMomentum",30,-1,1,"","Events");

  //ResolPullTauminusFroma1a1MeanEnergy=HConfig.GetTH1D(Name+"_ResolPullTauminusFroma1a1MeanEnergy","ResolPullTauminusFroma1a1MeanEnergy",30,-1,1,"","Events");
  //ResolPullTauminusFroma1a1MZEnergy=HConfig.GetTH1D(Name+"_ResolPullTauminusFroma1a1MZEnergy","ResolPullTauminusFroma1a1MZEnergy",30,-1,1,"","Events");
  ResolPullTauminusFroma1a1MeanMomentum=HConfig.GetTH1D(Name+"_ResolPullTauminusFroma1a1MeanMomentum","ResolPullTauminusFroma1a1MeanMomentum",30,-1,1,"","Events");
  ResolPullTauminusFroma1a1MZMomentum=HConfig.GetTH1D(Name+"_ResolPullTauminusFroma1a1MZMomentum","ResolPullTauminusFroma1a1MZMomentum",30,-1,1,"","Events");


  //ResolPullTauplusFroma1a1MeanEnergy=HConfig.GetTH1D(Name+"_ResolPullTauplusFroma1a1MeanEnergy","ResolPullTauplusFroma1a1MeanEnergy",30,-1,1,"","Events");
  //ResolPullTauplusFroma1a1MZEnergy=HConfig.GetTH1D(Name+"_ResolPullTauplusFroma1a1MZEnergy","ResolPullTauplusFroma1a1MZEnergy",30,-1,1,"","Events");
  ResolPullTauplusFroma1a1MeanMomentum=HConfig.GetTH1D(Name+"_ResolPullTauplusFroma1a1MeanMomentum","ResolPullTauplusFroma1a1MeanMomentum",30,-1,1,"","Events");
  ResolPullTauplusFroma1a1MZMomentum=HConfig.GetTH1D(Name+"_ResolPullTauplusFroma1a1MZMomentum","ResolPullTauplusFroma1a1MZMomentum",30,-1,1,"","Events");


  //ResolPullTauFroma1a1MeanEnergy=HConfig.GetTH1D(Name+"_ResolPullTauFroma1a1MeanEnergy","ResolPullTauFroma1a1MeanEnergy",30,-1,1,"","Events");
  //ResolPullTauFroma1a1MZEnergy=HConfig.GetTH1D(Name+"_ResolPullTauFroma1a1MZEnergy","ResolPullTauFroma1a1MZEnergy",30,-1,1,"","Events");
  ResolPullTauFroma1a1MeanMomentum=HConfig.GetTH1D(Name+"_ResolPullTauFroma1a1MeanMomentum","ResolPullTauFroma1a1MeanMomentum",30,-1,1,"","Events");
  ResolPullTauFroma1a1MZMomentum=HConfig.GetTH1D(Name+"_ResolPullTauFroma1a1MZMomentum","ResolPullTauFroma1a1MZMomentum",30,-1,1,"","Events");


  //ResolPullTauminusFroma1XMeanEnergy=HConfig.GetTH1D(Name+"_ResolPullTauminusFroma1XMeanEnergy","ResolPullTauminusFroma1XMeanEnergy",30,-1,1,"","Events");
  ResolPullTauminusFroma1XMeanMomentum=HConfig.GetTH1D(Name+"_ResolPullTauminusFroma1XMeanMomentum","ResolPullTauminusFroma1XMeanMomentum",30,-1,1,"","Events");
  //ResolPullXminusFroma1XMeanEnergy=HConfig.GetTH1D(Name+"_ResolPullXminusFroma1XMeanEnergy","ResolPullXminusFroma1XMeanEnergy",30,-1,1,"","Events");
  ResolPullXminusFroma1XMeanMomentum=HConfig.GetTH1D(Name+"_ResolPullXminusFroma1XMeanMomentum","ResolPullXminusFroma1XMeanMomentum",30,-1,1,"","Events");


  //ResolPullTauplusFroma1XMeanEnergy=HConfig.GetTH1D(Name+"_ResolPullTauplusFroma1XMeanEnergy","ResolPullTauplusFroma1XMeanEnergy",30,-1,1,"","Events");
  ResolPullTauplusFroma1XMeanMomentum=HConfig.GetTH1D(Name+"_ResolPullTauplusFroma1XMeanMomentum","ResolPullTauplusFroma1XMeanMomentum",30,-1,1,"","Events");
  //ResolPullXplusFroma1XMeanEnergy=HConfig.GetTH1D(Name+"_ResolPullXplusFroma1XMeanEnergy","ResolPullXplusFroma1XMeanEnergy",30,-1,1,"","Events");
  ResolPullXplusFroma1XMeanMomentum=HConfig.GetTH1D(Name+"_ResolPullXplusFroma1XMeanMomentum","ResolPullXplusFroma1XMeanMomentum",30,-1,1,"","Events");


  ResolPullXVtxIna1a1=HConfig.GetTH1D(Name+"_ResolPullXVtxIna1a1","ResolPullXVtxIna1a1",30,-1,1,"","Events");
  ResolPullYVtxIna1a1=HConfig.GetTH1D(Name+"_ResolPullYVtxIna1a1","ResolPullYVtxIna1a1",30,-1,1,"","Events");
  ResolPullZVtxIna1a1=HConfig.GetTH1D(Name+"_ResolPullZVtxIna1a1","ResolPullZVtxIna1a1",30,-1,1,"","Events");
  
  // tauminusa1a1MomentumVis=HConfig.GetTH1D(Name+"_tauminusa1a1MomentumVis","tauminusa1a1MomentumVis",15,20,200,"","Events");                             
  // tauplusa1a1MomentumVis=HConfig.GetTH1D(Name+"_tauplusa1a1MomentumVis","tauplusa1a1MomentumVis",15,20,200,"","Events");
  // InvariantMasstausa1a1Vis=HConfig.GetTH1D(Name+"_InvariantMasstausa1a1Vis","InvariantMasstausa1a1Vis",15,20,110,"","Events");
		    
  // tauminusa1a1MomentumMean=HConfig.GetTH1D(Name+"_tauminusa1a1MomentumMean","tauminusa1a1MomentumMean",15,20,200,"","Events");                             
  // tauplusa1a1MomentumMean=HConfig.GetTH1D(Name+"_tauplusa1a1MomentumMean","tauplusa1a1MomentumMean",15,20,200,"","Events");
  // InvariantMasstausa1a1Mean=HConfig.GetTH1D(Name+"_InvariantMasstausa1a1Mean","InvariantMasstausa1a1Mean",15,0,200,"","Events");

  tauminusa1a1MomentumPairConstraint=HConfig.GetTH1D(Name+"_tauminusa1a1MomentumPairConstraint","tauminusa1a1MomentumPairConstraint",15,0,200,"","Events");       
  tauplusa1a1MomentumPairConstraint=HConfig.GetTH1D(Name+"_tauplusa1a1MomentumPairConstraint","tauplusa1a1MomentumPairConstraint",15,0,200,"","Events");
  InvariantMasstausa1a1PairConstraint=HConfig.GetTH1D(Name+"_InvariantMasstausa1a1PairConstraint","InvariantMasstausa1a1PairConstraint",15,20,160,"","Events");
	
  tauminusa1XMomentumVis=HConfig.GetTH1D(Name+"_tauminusa1XMomentumVis","tauminusa1XMomentumVis",15,0,200,"","Events");                             
  tauplusa1XMomentumVis=HConfig.GetTH1D(Name+"_tauplusa1XMomentumVis","tauplusa1XMomentumVis",15,0,200,"","Events");
  InvariantMasstausa1XVis=HConfig.GetTH1D(Name+"_InvariantMasstausa1XVis","InvariantMasstausa1XVis",15,20,110,"","Events");
		    
  tauminusa1XMomentumMean=HConfig.GetTH1D(Name+"_tauminusa1XMomentumMean","tauminusa1XMomentumMean",15,0,200,"","Events");                             
  tauplusa1XMomentumMean=HConfig.GetTH1D(Name+"_tauplusa1XMomentumMean","tauplusa1XMomentumMean",15,0,200,"","Events");
  InvariantMasstausa1XMean=HConfig.GetTH1D(Name+"_InvariantMasstausa1XMean","InvariantMasstausa1XMean",15,20,160,"","Events");


  polarimetricAcopAngle=HConfig.GetTH1D(Name+"_polarimetricAcopAngle","GEF",60,0.,2*TMath::Pi(),"GEF","Events");
  
  polarimetricAcopAnglePVRefitNoBSOld=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBSOld","GEF with old refitted PV no BS",60,0.,2*TMath::Pi(),"GEF with old refitted PV no BS","Events");
  polarimetricAcopAnglePVRefitBSOld=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSOld","GEF with old refitted PV BS",60,0.,2*TMath::Pi(),"GEF with old refitted PV BS","Events");
  polarimetricAcopAnglePVRefitNoBSZNominalOld=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBSZNominalOld","GEF with old refitted PV no BS nominal Z",60,0.,2*TMath::Pi(),"GEF with old refitted PV no BS nominal Z","Events");
  polarimetricAcopAnglePVRefitBSZNominalOld=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSZNominalOld","GEF with old refitted PV BS nominal Z",60,0.,2*TMath::Pi(),"GEF with old refitted PV BS nominal Z","Events");
   
  // polarimetricAcopAnglePVRefitNoBSOneTrackRemovedOld=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBSOneTrackRemovedOld","GEF with old refitted PV no BS one track removed",60,0.,2*TMath::Pi(),"","Events");
  // polarimetricAcopAnglePVRefitBSOneTrackRemovedOld=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSOneTrackRemovedOld","GEF with old refitted PV BS one track removed",60,0.,2*TMath::Pi(),"","Events");
  // polarimetricAcopAnglePVRefitNoBSOneTrackRemovedZNominalOld=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBSOneTrackRemovedZNominalOld","GEF with old refitted PV no BS nominal Z one track removed",60,0.,2*TMath::Pi(),"","Events");
  // polarimetricAcopAnglePVRefitBSOneTrackRemovedZNominalOld=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSOneTrackRemovedZNominalOld","GEF with old refitted PV BS nominal Z one track removed",60,0.,2*TMath::Pi(),"","Events");

  polarimetricAcopAnglePVRefitNoBSTracksRemovedOld=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBSTracksRemovedOld","GEF personal PVRefit no BS constraint",60,0.,2*TMath::Pi(),"GEF personal PVRefit no BS constraint","Events");
  polarimetricAcopAnglePVRefitBSTracksRemovedOld=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSTracksRemovedOld","GEF personal PVRefit BS constraint",60,0.,2*TMath::Pi(),"GEF personal PVRefit BS constraint","Events");
  polarimetricAcopAnglePVRefitNoBSTracksRemovedZNominalOld=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBSTracksRemovedZNominalOld","GEF personal PVRefit no BS constraint nominal Z",60,0.,2*TMath::Pi(),"GEF personal PVRefit no BS constraint nominal Z","Events");
  polarimetricAcopAnglePVRefitBSTracksRemovedZNominalOld=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSTracksRemovedZNominalOld","GEF personal PVRefit BS constraint nominal Z",60,0.,2*TMath::Pi(),"GEF personal PVRefit BS constraint nominal Z","Events");

  polarimetricAcopAnglePVRefitNoBSNew=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBSNew","GEF PVRefit no BS constraint new",60,0.,2*TMath::Pi(),"GEF PVRefit no BS constraint new","Events");
  polarimetricAcopAnglePVRefitBSNew=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSNew","GEF PVRefit BS constraint new",60,0.,2*TMath::Pi(),"GEF PVRefit BS constraint new","Events");
  polarimetricAcopAnglePVRefitNoBSZNominalNew=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBSZNominalNew","GEF PVRefit no BS constraint nominal Z new",60,0.,2*TMath::Pi(),"GEF PVRefit no BS constraint nominal Z new","Events");
  polarimetricAcopAnglePVRefitBSZNominalNew=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSZNominalNew","GEF PVRefit BS constraint nominal Z new",60,0.,2*TMath::Pi(),"GEF PVRefit BS constraint nominal Z new","Events");
   
  // polarimetricAcopAnglePVRefitNoBSTIP=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBSTIP","GEF with old refitted PV no BS TIP",60,0.,2*TMath::Pi(),"","Events");
  // polarimetricAcopAnglePVRefitBSTIP=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSTIP","GEF with old refitted PV BS TIP",60,0.,2*TMath::Pi(),"","Events");
  // polarimetricAcopAnglePVRefitNoBSZNominalTIP=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBSZNominalTIP","GEF with old refitted PV no BS nominal Z TIP",60,0.,2*TMath::Pi(),"","Events");
  // polarimetricAcopAnglePVRefitBSZNominalTIP=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSZNominalTIP","GEF with old refitted PV BS nominal Z TIP",60,0.,2*TMath::Pi(),"","Events");


  polarimetricAcopAngleMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAngleMVA","GEF MVA",60,0.,2*TMath::Pi(),"GEF MVA","Events");
  
  polarimetricAcopAnglePVRefitNoBSOldMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBSOldMVA","GEF with old refitted PV no BS constraint MVA",60,0.,2*TMath::Pi(),"GEF with old refitted PV no BS constraint MVA","Events");
  polarimetricAcopAnglePVRefitBSOldMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSOldMVA","GEF with old refitted PV BS constraint MVA",60,0.,2*TMath::Pi(),"GEF with old refitted PV BS constraint MVA","Events");
  polarimetricAcopAnglePVRefitNoBSZNominalOldMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBSZNominalOldMVA","GEF with old refitted PV no BS constraint nominal Z MVA",60,0.,2*TMath::Pi(),"GEF with old refitted PV no BS constraint nominal Z MVA","Events");
  polarimetricAcopAnglePVRefitBSZNominalOldMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSZNominalOldMVA","GEF with old refitted PV BS constraint nominal Z MVA",60,0.,2*TMath::Pi(),"GEF with old refitted PV BS constraint nominal Z MVA","Events");
   
  // polarimetricAcopAnglePVRefitNoBSOneTrackRemovedOldMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBSOneTrackRemovedOldMVA","GEF with old refitted PV no BS one track removed MVA",60,0.,2*TMath::Pi(),"","Events");
  // polarimetricAcopAnglePVRefitBSOneTrackRemovedOldMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSOneTrackRemovedOldMVA","GEF with old refitted PV BS one track removed MVA",60,0.,2*TMath::Pi(),"","Events");
  // polarimetricAcopAnglePVRefitNoBSOneTrackRemovedZNominalOldMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBSOneTrackRemovedZNominalOldMVA","GEF with old refitted PV no BS nominal Z one track removed MVA",60,0.,2*TMath::Pi(),"","Events");
  // polarimetricAcopAnglePVRefitBSOneTrackRemovedZNominalOldMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSOneTrackRemovedZNominalOldMVA","GEF with old refitted PV BS nominal Z one track removed MVA",60,0.,2*TMath::Pi(),"","Events");

  polarimetricAcopAnglePVRefitNoBSTracksRemovedOldMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBSTracksRemovedOldMVA","GEF personal PVRefit no BS constraint MVA",60,0.,2*TMath::Pi(),"GEF personal PVRefit no BS constraint MVA","Events");
  polarimetricAcopAnglePVRefitBSTracksRemovedOldMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSTracksRemovedOldMVA","GEF personal PVRefit BS constraint MVA",60,0.,2*TMath::Pi(),"GEF personal PVRefit BS constraint MVA","Events");
  polarimetricAcopAnglePVRefitNoBSTracksRemovedZNominalOldMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBSTracksRemovedZNominalOldMVA","GEF personal PVRefit no BS constraint nominal Z MVA",60,0.,2*TMath::Pi(),"GEF personal PVRefit no BS constraint nominal Z MVA","Events");
  polarimetricAcopAnglePVRefitBSTracksRemovedZNominalOldMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSTracksRemovedZNominalOldMVA","GEF personal PVRefit BS constraint nominal Z MVA",60,0.,2*TMath::Pi(),"GEF personal PVRefit BS constraint nominal Z MVA","Events");

  polarimetricAcopAnglePVRefitNoBSNewMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBSNewMVA","GEF PVRefit no BS constraint new MVA",60,0.,2*TMath::Pi(),"GEF PVRefit no BS constraint new MVA","Events");
  polarimetricAcopAnglePVRefitBSNewMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSNewMVA","GEF PVRefit BS constraint new MVA",60,0.,2*TMath::Pi(),"GEF PVRefit BS constraint new MVA","Events");
  polarimetricAcopAnglePVRefitNoBSZNominalNewMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBSZNominalNewMVA","GEF PVRefit no BS constraint nominal Z new MVA",60,0.,2*TMath::Pi(),"GEF PVRefit no BS constraint nominal Z new MVA","Events");
  polarimetricAcopAnglePVRefitBSZNominalNewMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSZNominalNewMVA","GEF PVRefit BS constraint nominal Z new MVA",60,0.,2*TMath::Pi(),"GEF PVRefit BS constraint nominal Z new MVA","Events");

  polarimetricAcopAnglePVRefitOnlyNoBSTracksRemovedOld=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitOnlyNoBSTracksRemovedOld","GEF personal PVRefit Only no BS constraint",60,0.,2*TMath::Pi(),"GEF personal PVRefit no BS constraint","Events");
  polarimetricAcopAnglePVRefitOnlyBSTracksRemovedOld=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitOnlyBSTracksRemovedOld","GEF personal PVRefit Only BS constraint",60,0.,2*TMath::Pi(),"GEF personal PVRefit BS constraint","Events");
  polarimetricAcopAnglePVRefitOnlyNoBSTracksRemovedZNominalOld=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitOnlyNoBSTracksRemovedZNominalOld","GEF personal PVRefit Only no BS constraint nominal Z",60,0.,2*TMath::Pi(),"GEF personal PVRefit Only no BS constraint nominal Z","Events");
  polarimetricAcopAnglePVRefitOnlyBSTracksRemovedZNominalOld=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitOnlyBSTracksRemovedZNominalOld","GEF personal PVRefit Only BS constraint nominal Z",60,0.,2*TMath::Pi(),"GEF personal PVRefit Only BS constraint nominal Z","Events");

  polarimetricAcopAnglePVRefitOnlyNoBSNew=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitOnlyNoBSNew","GEF PVRefit Only no BS constraint",60,0.,2*TMath::Pi(),"GEF PVRefit Only no BS constraint","Events");
  polarimetricAcopAnglePVRefitOnlyBSNew=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitOnlyBSNew","GEF PVRefit Only BS constraint",60,0.,2*TMath::Pi(),"GEF PVRefit Only BS constraint","Events");
  polarimetricAcopAnglePVRefitOnlyNoBSZNominalNew=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitOnlyNoBSZNominalNew","GEF PVRefit Only no BS constraint nominal Z",60,0.,2*TMath::Pi(),"GEF PVRefit Only no BS constraint nominal Z","Events");
  polarimetricAcopAnglePVRefitOnlyBSZNominalNew=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitOnlyBSZNominalNew","GEF PVRefit Only BS constraint nominal Z",60,0.,2*TMath::Pi(),"GEF PVRefit Only BS constraint nominal Z","Events");



  polarimetricAcopAnglePVRefitOnlyNoBSTracksRemovedOldMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitOnlyNoBSTracksRemovedOldMVA","GEF personal PVRefit Only no BS constraint MVA",60,0.,2*TMath::Pi(),"GEF personal PVRefit Only no BS constraint MVA","Events");
  polarimetricAcopAnglePVRefitOnlyBSTracksRemovedOldMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitOnlyBSTracksRemovedOldMVA","GEF personal PVRefit Only BS constraint MVA",60,0.,2*TMath::Pi(),"GEF personal PVRefit Only BS constraint MVA","Events");
  polarimetricAcopAnglePVRefitOnlyNoBSTracksRemovedZNominalOldMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitOnlyNoBSTracksRemovedZNominalOldMVA","GEF personal PVRefit Only no BS constraint nominal Z MVA",60,0.,2*TMath::Pi(),"GEF personal PVRefit Only no BS constraint nominal Z MVA","Events");
  polarimetricAcopAnglePVRefitOnlyBSTracksRemovedZNominalOldMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitOnlyBSTracksRemovedZNominalOldMVA","GEF personal PVRefit Only BS constraint nominal Z MVA",60,0.,2*TMath::Pi(),"GEF personal PVRefit Only BS constraint nominal Z MVA","Events");

  polarimetricAcopAnglePVRefitOnlyNoBSNewMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVOnlyRefitNoBSNewMVA","GEF PVRefit Only no BS constraint MVA",60,0.,2*TMath::Pi(),"GEF PVRefit Only no BS constraint MVA","Events");
  polarimetricAcopAnglePVRefitOnlyBSNewMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVOnlyRefitBSNewMVA","GEF PVRefit Only BS constraint MVA",60,0.,2*TMath::Pi(),"GEF PVRefit Only BS constraint MVA","Events");
  polarimetricAcopAnglePVRefitOnlyNoBSZNominalNewMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVOnlyRefitNoBSZNominalNewMVA","GEF PVRefit Only no BS constraint nominal Z MVA",60,0.,2*TMath::Pi(),"GEF PVRefit Only no BS constraint nominal Z MVA","Events");
  polarimetricAcopAnglePVRefitOnlyBSZNominalNewMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVOnlyRefitBSZNominalNewMVA","GEF PVRefit Only BS constraint nominal Z MVA",60,0.,2*TMath::Pi(),"GEF PVRefit Only BS constraint nominal Z MVA","Events");

  

  PVXResol=HConfig.GetTH1D(Name+"_PVXResol","PV_{X} pull",50,-0.1,0.1,"PV_{X} pull","Events");
  PVXNoBSOldResol=HConfig.GetTH1D(Name+"_PVXNoBSOldResol","no BS PV_{X}Refit Old pull",50,-0.1,0.1,"no BS PV_{X}Refit Old pull","Events");
  PVXNoBSTracksRemovedOldResol=HConfig.GetTH1D(Name+"_PVXNoBSTracksRemovedOldResol","No BS personal PV_{X}Refit pull",50,-0.1,0.1,"No BS personal PV_{X}Refit pull","Events");
  PVXNoBSNewResol=HConfig.GetTH1D(Name+"_PVXNoBSNewResol","No BS PV_{X}Refit pull",50,-0.1,0.1,"No BS PV_{X}Refit pull","Events");
  PVXBSOldResol=HConfig.GetTH1D(Name+"_PVXBSOldResol","BS PV_{X}Refit Old pull",50,-0.1,0.1,"BS PV_{X}Refit Old pull","Events");
  PVXBSTracksRemovedOldResol=HConfig.GetTH1D(Name+"_PVXBSTracksRemovedOldResol","BS personal PV_{X}Refit pull",50,-0.1,0.1,"BS personal PV_{X}Refit pull","Events");
  PVXBSNewResol=HConfig.GetTH1D(Name+"_PVXBSNewResol","BS PV_{X}Refit pull",50,-0.1,0.1,"BS PV_{X}Refit pull","Events");
  
  PVYResol=HConfig.GetTH1D(Name+"_PVYResol","PV_{Y} pull",50,-0.1,0.1,"PV_{Y} pull","Events");
  PVYNoBSOldResol=HConfig.GetTH1D(Name+"_PVYNoBSOldResol","no BS PV_{Y}Refit Old pull",50,-0.1,0.1,"no BS PV_{Y}Refit Old pull","Events");
  PVYNoBSTracksRemovedOldResol=HConfig.GetTH1D(Name+"_PVYNoBSTracksRemovedOldResol","No BS personal PV_{Y}Refit pull",50,-0.1,0.1,"No BS personal PV_{Y}Refit pull","Events");
  PVYNoBSNewResol=HConfig.GetTH1D(Name+"_PVYNoBSNewResol","No BS PV_{Y}Refit pull",50,-0.1,0.1,"No BS PV_{Y}Refit pull","Events");
  PVYBSOldResol=HConfig.GetTH1D(Name+"_PVYBSOldResol","BS PV_{Y}Refit Old pull",50,-0.1,0.1,"BS PV_{Y}Refit Old pull","Events");
  PVYBSTracksRemovedOldResol=HConfig.GetTH1D(Name+"_PVYBSTracksRemovedOldResol","BS personal PV_{Y}Refit pull",50,-0.1,0.1,"BS personal PV_{Y}Refit pull","Events");
  PVYBSNewResol=HConfig.GetTH1D(Name+"_PVYBSNewResol","BS PV_{Y}Refit pull",50,-0.1,0.1,"BS PV_{Y}Refit pull","Events");

  PVZResol=HConfig.GetTH1D(Name+"_PVZResol","PV_{Z} pull",50,-0.01,0.01,"PV_{Z} pull","Events");
  PVZNoBSOldResol=HConfig.GetTH1D(Name+"_PVZNoBSOldResol","no BS PV_{Z}Refit Old pull",50,-0.01,0.01,"no BS PV_{Z}Refit Old pull","Events");
  PVZNoBSTracksRemovedOldResol=HConfig.GetTH1D(Name+"_PVZNoBSTracksRemovedOldResol","No BS personal PV_{Z}Refit pull",50,-0.01,0.01,"No BS personal PV_{Z}Refit pull","Events");
  PVZNoBSNewResol=HConfig.GetTH1D(Name+"_PVZNoBSNewResol","No BS PV_{Z}Refit pull",50,-0.01,0.01,"No BS PV_{Z}Refit pull","Events");
  PVZBSOldResol=HConfig.GetTH1D(Name+"_PVZBSOldResol","BS PV_{Z}Refit Old pull",50,-0.01,0.01,"BS PV_{Z}Refit Old pull","Events");
  PVZBSTracksRemovedOldResol=HConfig.GetTH1D(Name+"_PVZBSTracksRemovedOldResol","BS personal PV_{Z}Refit pull",50,-0.01,0.01,"BS personal PV_{Z}Refit pull","Events");
  PVZBSNewResol=HConfig.GetTH1D(Name+"_PVZBSNewResol","BS PV_{Z}Refit pull",50,-0.01,0.01,"BS PV_{Z}Refit pull","Events");
  
  
  PVXNoBSTracksRemovedOldOnlyResol=HConfig.GetTH1D(Name+"_PVXNoBSTracksRemovedOldOnlyResol","No BS personal PV_{X}Refit Only pull",50,-0.1,0.1,"No BS personal PV_{X}Refit Only pull","Events");
  PVXNoBSNewOnlyResol=HConfig.GetTH1D(Name+"_PVXNoBSNewOnlyResol","No BS PV_{X}Refit Only pull",50,-0.1,0.1,"No BS PV_{X}Refit Only pull","Events");
  PVXBSTracksRemovedOldOnlyResol=HConfig.GetTH1D(Name+"_PVXBSTracksRemovedOldOnlyResol","BS personal PV_{X}Refit Only pull",50,-0.1,0.1,"BS personal PV_{X}Refit Only pull","Events");
  PVXBSNewOnlyResol=HConfig.GetTH1D(Name+"_PVXBSNewOnlyResol","BS PV_{X}Refit Only pull",50,-0.1,0.1,"BS PV_{X}Refit Only pull","Events");

  PVYNoBSTracksRemovedOldOnlyResol=HConfig.GetTH1D(Name+"_PVYNoBSTracksRemovedOldOnlyResol","No BS personal PV_{Y}Refit Only pull",50,-0.1,0.1,"No BS personal PV_{Y}Refit Only pull","Events");
  PVYNoBSNewOnlyResol=HConfig.GetTH1D(Name+"_PVYNoBSNewOnlyResol","No BS PV_{Y}Refit Only pull",50,-0.1,0.1,"No BS PV_{Y}Refit Only pull","Events");
  PVYBSTracksRemovedOldOnlyResol=HConfig.GetTH1D(Name+"_PVYBSTracksRemovedOldOnlyResol","BS personal PV_{Y}Refit Only pull",50,-0.1,0.1,"BS personal PV_{Y}Refit Only pull","Events");
  PVYBSNewOnlyResol=HConfig.GetTH1D(Name+"_PVYBSNewOnlyResol","BS PV_{Y}Refit Only pull",50,-0.1,0.1,"BS PV_{Y}Refit Only pull","Events");

  PVZNoBSTracksRemovedOldOnlyResol=HConfig.GetTH1D(Name+"_PVZNoBSTracksRemovedOldOnlyResol","No BS personal PV_{Z}Refit Only pull",50,-0.01,0.01,"No BS personal PV_{Z}Refit Only pull","Events");
  PVZNoBSNewOnlyResol=HConfig.GetTH1D(Name+"_PVZNoBSNewOnlyResol","No BS PV_{Z}Refit Only pull",50,-0.01,0.01,"No BS PV_{Z}Refit Only pull","Events");
  PVZBSTracksRemovedOldOnlyResol=HConfig.GetTH1D(Name+"_PVZBSTracksRemovedOldOnlyResol","BS personal PV_{Z}Refit Only pull",50,-0.01,0.01,"BS personal PV_{Z}Refit Only pull","Events");
  PVZBSNewOnlyResol=HConfig.GetTH1D(Name+"_PVZBSNewOnlyResol","BS PV_{Z}Refit Only pull",50,-0.01,0.01,"BS PV_{Z}Refit Only pull","Events");
  
  // polarimetricAcopAnglePVRefitNoBSTIPMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBSTIPMVA","GEF with old refitted PV no BS TIP MVA",60,0.,2*TMath::Pi(),"","Events");
  // polarimetricAcopAnglePVRefitBSTIPMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSTIPMVA","GEF with old refitted PV BS TIP MVA",60,0.,2*TMath::Pi(),"","Events");
  // polarimetricAcopAnglePVRefitNoBSZNominalTIPMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBSZNominalTIPMVA","GEF with old refitted PV no BS nominal Z TIP MVA",60,0.,2*TMath::Pi(),"","Events");
  // polarimetricAcopAnglePVRefitBSZNominalTIPMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSZNominalTIPMVA","GEF with old refitted PV BS nominal Z TIP MVA",60,0.,2*TMath::Pi(),"","Events");
  
  polarimetricAcopAngleTruthA1=HConfig.GetTH1D(Name+"_polarimetricAcopAngleTruthA1","polarimetricAcopAngleTruthA1",60,0,2*TMath::Pi(),"","Events");
  
  // polarimetricAcopAngleTruthRho=HConfig.GetTH1D(Name+"_polarimetricAcopAngleTruthRho","polarimetricAcopAngleTruthRho",60,0,2*TMath::Pi(),"","Events");
  // polarimetricAcopAngleTruthPi=HConfig.GetTH1D(Name+"_polarimetricAcopAngleTruthPi","polarimetricAcopAngleTruthPi",60,0,2*TMath::Pi(),"","Events");

  // polarimetricAcopAngle30=HConfig.GetTH1D(Name+"_polarimetricAcopAngle30","polarimetricAcopAngle30",60,0.,2*TMath::Pi(),"","Events");
  // polarimetricAcopAngle25=HConfig.GetTH1D(Name+"_polarimetricAcopAngle25","polarimetricAcopAngle25",60,0.,2*TMath::Pi(),"","Events");
  // polarimetricAcopAngle20=HConfig.GetTH1D(Name+"_polarimetricAcopAngle20","polarimetricAcopAngle20",60,0.,2*TMath::Pi(),"","Events");
  // polarimetricAcopAngle15=HConfig.GetTH1D(Name+"_polarimetricAcopAngle15","polarimetricAcopAngle15",60,0.,2*TMath::Pi(),"","Events");
  // polarimetricAcopAngle10=HConfig.GetTH1D(Name+"_polarimetricAcopAngle10","polarimetricAcopAngle10",60,0.,2*TMath::Pi(),"","Events");
  // polarimetricAcopAngle5=HConfig.GetTH1D(Name+"_polarimetricAcopAngle5","polarimetricAcopAngle5",60,0.,2*TMath::Pi(),"","Events");
  
  
  // AcolAngle=HConfig.GetTH1D(Name+"_AcolAngle","AcolAngle",60,2.5,3.16,"","Events");
  // AcolAngleSVFit=HConfig.GetTH1D(Name+"_AcolAngleSVFit","AcolAngleSVFit",60,2.5,3.16,"","Events");
  // AcolAngleTruth=HConfig.GetTH1D(Name+"_AcolAngleTruth","AcolAngleTruth",60,2.5,3.16,"","Events");

  // polarimetricAcopAnglePVRefitNoBS=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBS","polarimetricAcopAnglePVRefitNoBS",60,0.,2*TMath::Pi(),"","Events");
  // polarimetricAcopAnglePVRefitBS=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBS","polarimetricAcopAnglePVRefitBS",60,0.,2*TMath::Pi(),"","Events");
  // polarimetricAcopAngleMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAngleMVA","polarimetricAcopAngleMVA",60,0.,2*TMath::Pi(),"","Events");
  // polarimetricAcopAnglePVRefitNoBSMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBSMVA","polarimetricAcopAnglePVRefitNoBSMVA",60,0.,2*TMath::Pi(),"","Events");
  // polarimetricAcopAnglePVRefitBSMVA=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSMVA","polarimetricAcopAnglePVRefitBSMVA",60,0.,2*TMath::Pi(),"","Events");
  
  polarimetricAcopAngleDecayPlane=HConfig.GetTH1D(Name+"_polarimetricAcopAngleDecayPlane","polarimetricAcopAngleDecayPlane",60,0.,2*TMath::Pi(),"","Events");
  
  // polarimetricAcopAnglePtTruthA1=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePtTruthA1","polarimetricAcopAnglePtTruthA1",60,0.,2*TMath::Pi(),"","Events");
  // polarimetricAcopAngleMVAPtTruthA1=HConfig.GetTH1D(Name+"_polarimetricAcopAngleMVAPtTruthA1","polarimetricAcopAngleMVAPtTruthA1",60,0.,2*TMath::Pi(),"","Events");


  // polarimetricAcopAnglePtTruthRho=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePtTruthRho","polarimetricAcopAnglePtTruthRho",60,0.,2*TMath::Pi(),"","Events");
  // polarimetricAcopAngleMVAPtTruthRho=HConfig.GetTH1D(Name+"_polarimetricAcopAngleMVAPtTruthRho","polarimetricAcopAngleMVAPtTruthRho",60,0.,2*TMath::Pi(),"","Events");


  // polarimetricAcopAnglePtTruthPi=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePtTruthPi","polarimetricAcopAnglePtTruthPi",60,0.,2*TMath::Pi(),"","Events");
  // polarimetricAcopAngleMVAPtTruthPi=HConfig.GetTH1D(Name+"_polarimetricAcopAngleMVAPtTruthPi","polarimetricAcopAngleMVAPtTruthPi",60,0.,2*TMath::Pi(),"","Events");


   polarimetricAcopAngleSVFit=HConfig.GetTH1D(Name+"_polarimetricAcopAngleSVFit","polarimetricAcopAngleSVFit",60,0.,2*TMath::Pi(),"SVfit","Events");
   polarimetricAcopAngleMVASVFit=HConfig.GetTH1D(Name+"_polarimetricAcopAngleMVASVFit","polarimetricAcopAngleMVASVFit",60,0.,2*TMath::Pi(),"SVFit MVA","Events");


  // polarimetricAcopAngleSVFitRho=HConfig.GetTH1D(Name+"_polarimetricAcopAngleSVFitRho","polarimetricAcopAngleSVFitRho",60,0.,2*TMath::Pi(),"","Events");
  // polarimetricAcopAngleMVASVFitRho=HConfig.GetTH1D(Name+"_polarimetricAcopAngleMVASVFitRho","polarimetricAcopAngleMVASVFitRho",60,0.,2*TMath::Pi(),"","Events");
  
  // polarimetricAcopAngleSVFitPi=HConfig.GetTH1D(Name+"_polarimetricAcopAngleSVFitPi","polarimetricAcopAngleSVFitPi",60,0.,2*TMath::Pi(),"","Events");
  // polarimetricAcopAngleMVASVFitPi=HConfig.GetTH1D(Name+"_polarimetricAcopAngleMVASVFitPi","polarimetricAcopAngleMVASVFitPi",60,0.,2*TMath::Pi(),"","Events");
  
  
  // polarimetricAcopAngleBackground=HConfig.GetTH1D(Name+"_polarimetricAcopAngleBackground","polarimetricAcopAngleBackground",60,0.,2*TMath::Pi(),"","Events");
  // polarimetricAcopAnglePVRefitNoBSBackground=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBSBackground","polarimetricAcopAnglePVRefitNoBSBackground",60,0.,2*TMath::Pi(),"","Events");
  // polarimetricAcopAnglePVRefitBSBackground=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSBackground","polarimetricAcopAnglePVRefitBSBackground",60,0.,2*TMath::Pi(),"","Events");
  // polarimetricAcopAngleMVABackground=HConfig.GetTH1D(Name+"_polarimetricAcopAngleMVABackground","polarimetricAcopAngleMVABackground",60,0.,2*TMath::Pi(),"","Events");
  // polarimetricAcopAnglePVRefitNoBSMVABackground=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBSMVABackground","polarimetricAcopAnglePVRefitNoBSMVABackground",60,0.,2*TMath::Pi(),"","Events");
  // polarimetricAcopAnglePVRefitBSMVABackground=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitBSMVABackground","polarimetricAcopAnglePVRefitBSMVABackground",60,0.,2*TMath::Pi(),"","Events");

  // polarimetricAcopAnglePVRefitNoBSWithoutWSpin=HConfig.GetTH1D(Name+"_polarimetricAcopAnglePVRefitNoBSWithoutWSpin","polarimetricAcopAnglePVRefitNoBSWithoutWSpin",60,0.,2*TMath::Pi(),"","Events"); 

  
   test=HConfig.GetTH1D(Name+"_test","test",60,0.,2*TMath::Pi(),"","Events"); 

   PurityDM=HConfig.GetTH1D(Name+"_PurityDM","PurityDM",29,0.,29,"","Events");
   PurityNewMVA=HConfig.GetTH1D(Name+"_PurityNewMVA","PurityNewMVA",29,0.,29,"","Events");


   TauSVFitPxResPull=HConfig.GetTH1D(Name+"_TauSVFitPxResPull","TauSVFitPxResPull",50,-2,2,"P_{X} pull of #tau with SVFit","Events");
   TauSVFitPyResPull=HConfig.GetTH1D(Name+"_TauSVFitPyResPull","TauSVFitPyResPull",50,-2,2,"P_{Y} pull of #tau with SVFit","Events");
   TauSVFitPzResPull=HConfig.GetTH1D(Name+"_TauSVFitPzResPull","TauSVFitPzResPull",50,-2,2,"P_{Z} pull of #tau with SVFit","Events");

   TauPxResPull=HConfig.GetTH1D(Name+"_TauPxResPull","TauPxResPull",50,-2,2,"P_{X} pull of #tau","Events");
   TauPyResPull=HConfig.GetTH1D(Name+"_TauPyResPull","TauPyResPull",50,-2,2,"P_{Y} pull of #tau","Events");
   TauPzResPull=HConfig.GetTH1D(Name+"_TauPzResPull","TauPzResPull",50,-2,2,"P_{Z} pull of #tau","Events");

   TauSVFitPxResPullMVA=HConfig.GetTH1D(Name+"_TauSVFitPxResPullMVA","TauSVFitPxResPullMVA",50,-2,2,"P_{X} pull of #tau with SVFit","Events");
   TauSVFitPyResPullMVA=HConfig.GetTH1D(Name+"_TauSVFitPyResPullMVA","TauSVFitPyResPullMVA",50,-2,2,"P_{Y} pull of #tau with SVFit","Events");
   TauSVFitPzResPullMVA=HConfig.GetTH1D(Name+"_TauSVFitPzResPullMVA","TauSVFitPzResPullMVA",50,-2,2,"P_{Z} pull of #tau with SVFit","Events");

   TauPxResPullMVA=HConfig.GetTH1D(Name+"_TauPxResPullMVA","TauPxResPullMVA",50,-2,2,"P_{X} pull of #tau","Events");
   TauPyResPullMVA=HConfig.GetTH1D(Name+"_TauPyResPullMVA","TauPyResPullMVA",50,-2,2,"P_{Y} pull of #tau","Events");
   TauPzResPullMVA=HConfig.GetTH1D(Name+"_TauPzResPullMVA","TauPzResPullMVA",50,-2,2,"P_{Z} pull of #tau","Events");

  
   Selection::ConfigureHistograms();   //   do not remove
   HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);  // do not remove
}

 

void  HTauTau::Store_ExtraDist(){

  //every new histo should be addedd to Extradist1d vector, just push it back;
  Extradist1d.push_back(&Tau1PT);
  Extradist1d.push_back(&Tau1E);
  Extradist1d.push_back(&Tau1Mass);
  Extradist1d.push_back(&Tau1Phi);
  Extradist1d.push_back(&Tau1Eta);
  Extradist1d.push_back(&Tau1dz);
  Extradist1d.push_back(&Tau1HPSDecayMode);
  
  Extradist1d.push_back(&Tau2PT);
  Extradist1d.push_back(&Tau2E);
  Extradist1d.push_back(&Tau2Mass);
  Extradist1d.push_back(&Tau2Phi);
  Extradist1d.push_back(&Tau2Eta);
  Extradist1d.push_back(&Tau2dz);
  Extradist1d.push_back(&Tau2HPSDecayMode);
  /*
    Extradist1d.push_back(&Tau1isolation);
    Extradist1d.push_back(&Tau2isolation);
  
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
  Extradist1d.push_back(&TauTauTruthMass);
  //Extradist1d.push_back(&TauTauFullMass);
  
  Extradist1d.push_back(&QCDShape);
  Extradist1d.push_back(&NQCD);
  // Extradist1d.push_back(&TauTauFullMass_B);
  // Extradist1d.push_back(&TauTauFullMass_C);
  // Extradist1d.push_back(&TauTauFullMass_D);

  Extradist1d.push_back(&MET);
  Extradist1d.push_back(&METphi);
  Extradist1d.push_back(&PUPPImet);
  Extradist1d.push_back(&PUPPImetphi);
  Extradist1d.push_back(&TransverseMass);

  Extradist1d.push_back(&NPrimeVtx);
  Extradist1d.push_back(&NPU);
  Extradist1d.push_back(&RHO);


  Extradist1d.push_back(&NbJets);

  Extradist1d.push_back(&h_SVFitMass);
  Extradist1d.push_back(&h_SVFitStatus);
  // Extradist1d.push_back(&svfTau1E);
  // Extradist1d.push_back(&svfTau2E);
  
  // Extradist1d.push_back(&PhiDatasvfitpipi);
  // Extradist1d.push_back(&PhiDatasvfitpirho);
  // Extradist1d.push_back(&PhiDatasvfitpia1);
  // Extradist1d.push_back(&PhiDatasvfitrhorho);
  // Extradist1d.push_back(&PhiDatasvfitrhoa1);
  // Extradist1d.push_back(&PhiDatasvfita1a1);
  /*
    Extradist1d.push_back(&PhiDatavispipi);
    Extradist1d.push_back(&PhiDatavispirho);
    Extradist1d.push_back(&PhiDatavispia1);
    Extradist1d.push_back(&PhiDatavisrhorho);
    Extradist1d.push_back(&PhiDatavisrhoa1);
    Extradist1d.push_back(&PhiDatavisa1a1);
  */
  
  // Extradist1d.push_back(&Etasvfit);
  // Extradist1d.push_back(&Phisvfitpipi);
  // Extradist1d.push_back(&Phisvfitpirho);
  // //Extradist1d.push_back(&Phisvfitlpi);
  // //Extradist1d.push_back(&Phisvfitlrho);
  // Extradist1d.push_back(&Phisvfitpia1);
  // Extradist1d.push_back(&Phisvfitrhorho);
  // Extradist1d.push_back(&Phisvfitrhoa1);
  // //Extradist1d.push_back(&Phisvfitla1);
  // Extradist1d.push_back(&Phisvfita1a1);
  // Extradist1d.push_back(&Thetasvfit);
  
  Extradist1d.push_back(&Etavis);
  Extradist1d.push_back(&Phivispipi);
  Extradist1d.push_back(&Phivispirho);
  //Extradist1d.push_back(&Phivislpi);
  // Extradist1d.push_back(&Phivislrho);
  Extradist1d.push_back(&Phivispia1);
  Extradist1d.push_back(&Phivisrhorho);
  Extradist1d.push_back(&Phivisrhoa1);
  //Extradist1d.push_back(&Phivisla1);
  Extradist1d.push_back(&Phivisa1a1);
  Extradist1d.push_back(&Thetavis);
  
  Extradist1d.push_back(&Etatruth);
  Extradist1d.push_back(&Phitruthpipi);
  Extradist1d.push_back(&Phitruthpirho);
  //Extradist1d.push_back(&Phitruthlpi);
  // Extradist1d.push_back(&Phitruthlrho);
  Extradist1d.push_back(&Phitruthpia1);
  Extradist1d.push_back(&Phitruthrhorho);
  Extradist1d.push_back(&Phitruthrhoa1);
  //Extradist1d.push_back(&Phitruthla1);
  Extradist1d.push_back(&Phitrutha1a1);
  Extradist1d.push_back(&Thetatruth);
  /*
    Extradist1d.push_back(&PhiSvFitRespipi);
    Extradist1d.push_back(&PhiSvFitRespirho);
    Extradist1d.push_back(&PhiSvFitReslpi);
    Extradist1d.push_back(&PhiSvFitReslrho);
    Extradist1d.push_back(&PhiSvFitRespia1);
    Extradist1d.push_back(&PhiSvFitResrhorho);
    Extradist1d.push_back(&PhiSvFitResrhoa1);
    Extradist1d.push_back(&PhiSvFitResla1);
    Extradist1d.push_back(&PhiSvFitResa1a1);
  
    Extradist1d.push_back(&PhiVisRespipi);
    Extradist1d.push_back(&PhiVisRespirho);
    Extradist1d.push_back(&PhiVisReslpi);
    Extradist1d.push_back(&PhiVisReslrho);
    Extradist1d.push_back(&PhiVisRespia1);
    Extradist1d.push_back(&PhiVisResrhorho);
    Extradist1d.push_back(&PhiVisResrhoa1);
    Extradist1d.push_back(&PhiVisResla1);
    Extradist1d.push_back(&PhiVisResa1a1);
  

    Extradist1d.push_back(&TauTauFullPtRes);
    Extradist1d.push_back(&TauTauFullEtaRes);
    Extradist1d.push_back(&TauTauFullPhiRes);

    Extradist1d.push_back(&TauplusFullPtRes);
    Extradist1d.push_back(&TauplusFullEtaRes);
    Extradist1d.push_back(&TauplusFullPhiRes);

    Extradist1d.push_back(&TauminusFullPtRes);
    Extradist1d.push_back(&TauminusFullEtaRes);
    Extradist1d.push_back(&TauminusFullPhiRes);
  
    Extradist1d.push_back(&TauTauVisPtRes);
    Extradist1d.push_back(&TauTauVisEtaRes);
    Extradist1d.push_back(&TauTauVisPhiRes);

    Extradist1d.push_back(&TauplusVisPtRes);
    Extradist1d.push_back(&TauplusVisEtaRes);
    Extradist1d.push_back(&TauplusVisPhiRes);
    
    Extradist1d.push_back(&TauminusVisPtRes);
    Extradist1d.push_back(&TauminusVisEtaRes);
    Extradist1d.push_back(&TauminusVisPhiRes);

    Extradist1d.push_back(&TauTauFullPtResPull);
    Extradist1d.push_back(&TauTauFullEtaResPull);
    Extradist1d.push_back(&TauTauFullPhiResPull);

    Extradist1d.push_back(&TauplusFullPtResPull);
    Extradist1d.push_back(&TauplusFullEtaResPull);
    Extradist1d.push_back(&TauplusFullPhiResPull);

    Extradist1d.push_back(&TauminusFullPtResPull);
    Extradist1d.push_back(&TauminusFullEtaResPull);
    Extradist1d.push_back(&TauminusFullPhiResPull);

    Extradist1d.push_back(&TauTauVisPtResPull);
    Extradist1d.push_back(&TauTauVisEtaResPull);
    Extradist1d.push_back(&TauTauVisPhiResPull);

    Extradist1d.push_back(&TauplusVisPtResPull);
    Extradist1d.push_back(&TauplusVisEtaResPull);
    Extradist1d.push_back(&TauplusVisPhiResPull);
    
    Extradist1d.push_back(&TauminusVisPtResPull);
    Extradist1d.push_back(&TauminusVisEtaResPull);
    Extradist1d.push_back(&TauminusVisPhiResPull);
  */
  //Extradist1d.push_back(&DRTruth);
  //Extradist1d.push_back(&DRFull);
  //Extradist1d.push_back(&DRFullTruth);
  //Extradist1d.push_back(&DRVisTruth);
  Extradist1d.push_back(&Pi0EnergyRes);
  Extradist1d.push_back(&Pi0EnergyResPull);

  Extradist1d.push_back(&ZPtVis);

  Extradist2d.push_back(&NewPhivsDeltaPhi);
  Extradist2d.push_back(&NewPhivsDeltaEta);
  Extradist2d.push_back(&NewPhivsPhiproton);
  Extradist2d.push_back(&NewPhivsPhiTauplus);
  Extradist2d.push_back(&NewPhivsEtaproton);
  Extradist2d.push_back(&NewPhivsEtaTauplus);
  Extradist2d.push_back(&NewPhivsZPt);
  Extradist1d.push_back(&NewPhiSignal);
  Extradist1d.push_back(&NewPhiQCD);


  Extradist1d.push_back(&IstauminusvisPhysical);
  Extradist1d.push_back(&IstauplusvisPhysical);
  Extradist1d.push_back(&IsPairPhysical);

  //Extradist1d.push_back(&ResolPullTauTauFroma1a1MeanEnergy);
  //Extradist1d.push_back(&ResolPullTauTauFroma1a1MZEnergy);
  Extradist1d.push_back(&ResolPullTauTauFroma1a1MeanMomentum);
  Extradist1d.push_back(&ResolPullTauTauFroma1a1MZMomentum);

  //Extradist1d.push_back(&ResolPullTauminusFroma1a1MeanEnergy);
  //Extradist1d.push_back(&ResolPullTauminusFroma1a1MZEnergy);
  Extradist1d.push_back(&ResolPullTauminusFroma1a1MeanMomentum);
  Extradist1d.push_back(&ResolPullTauminusFroma1a1MZMomentum);

  //Extradist1d.push_back(&ResolPullTauplusFroma1a1MeanEnergy);
  //Extradist1d.push_back(&ResolPullTauplusFroma1a1MZEnergy);
  Extradist1d.push_back(&ResolPullTauplusFroma1a1MeanMomentum);
  Extradist1d.push_back(&ResolPullTauplusFroma1a1MZMomentum);

  //Extradist1d.push_back(&ResolPullTauFroma1a1MeanEnergy);
  //Extradist1d.push_back(&ResolPullTauFroma1a1MZEnergy);
  Extradist1d.push_back(&ResolPullTauFroma1a1MeanMomentum);
  Extradist1d.push_back(&ResolPullTauFroma1a1MZMomentum);

  //Extradist1d.push_back(&ResolPullTauminusFroma1XMeanEnergy);
  Extradist1d.push_back(&ResolPullTauminusFroma1XMeanMomentum);
  //Extradist1d.push_back(&ResolPullXminusFroma1XMeanEnergy);
  Extradist1d.push_back(&ResolPullXminusFroma1XMeanMomentum);

  //Extradist1d.push_back(&ResolPullTauplusFroma1XMeanEnergy);
  Extradist1d.push_back(&ResolPullTauplusFroma1XMeanMomentum); 
  //Extradist1d.push_back(&ResolPullXplusFroma1XMeanEnergy);
  Extradist1d.push_back(&ResolPullXplusFroma1XMeanMomentum);


  Extradist1d.push_back(&ResolPullXVtxIna1a1);
  Extradist1d.push_back(&ResolPullYVtxIna1a1);
  Extradist1d.push_back(&ResolPullZVtxIna1a1);


  // Extradist1d.push_back(&tauminusa1a1MomentumVis);                        
  // Extradist1d.push_back(&tauplusa1a1MomentumVis);
  // Extradist1d.push_back(&InvariantMasstausa1a1Vis);
		    
  // Extradist1d.push_back(&tauminusa1a1MomentumMean);
  // Extradist1d.push_back(&tauplusa1a1MomentumMean);
  // Extradist1d.push_back(&InvariantMasstausa1a1Mean);

  Extradist1d.push_back(&tauminusa1a1MomentumPairConstraint);
  Extradist1d.push_back(&tauplusa1a1MomentumPairConstraint);
  Extradist1d.push_back(&InvariantMasstausa1a1PairConstraint);
	
  Extradist1d.push_back(&tauminusa1XMomentumVis);
  Extradist1d.push_back(&tauplusa1XMomentumVis);
  Extradist1d.push_back(&InvariantMasstausa1XVis);
		    
  Extradist1d.push_back(&tauminusa1XMomentumMean);
  Extradist1d.push_back(&tauplusa1XMomentumMean);
  Extradist1d.push_back(&InvariantMasstausa1XMean);

  
  Extradist1d.push_back(&polarimetricAcopAngle);
  
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBSOld);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSOld);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBSZNominalOld);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSZNominalOld);
   
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBSOneTrackRemovedOld);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSOneTrackRemovedOld);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBSOneTrackRemovedZNominalOld);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSOneTrackRemovedZNominalOld);

  Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBSTracksRemovedOld);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSTracksRemovedOld);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBSTracksRemovedZNominalOld);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSTracksRemovedZNominalOld);

  Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBSNew);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSNew);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBSZNominalNew);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSZNominalNew);
   
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitOnlyNoBSTracksRemovedOld);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitOnlyBSTracksRemovedOld);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitOnlyNoBSTracksRemovedZNominalOld);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitOnlyBSTracksRemovedZNominalOld);

  Extradist1d.push_back(&polarimetricAcopAnglePVRefitOnlyNoBSNew);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitOnlyBSNew);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitOnlyNoBSZNominalNew);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitOnlyBSZNominalNew);
   
  
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBSTIP);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSTIP);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBSZNominalTIP);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSZNominalTIP);


  Extradist1d.push_back(&polarimetricAcopAngleMVA);
  
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBSOldMVA);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSOldMVA);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBSZNominalOldMVA);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSZNominalOldMVA);
   
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBSOneTrackRemovedOldMVA);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSOneTrackRemovedOldMVA);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBSOneTrackRemovedZNominalOldMVA);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSOneTrackRemovedZNominalOldMVA);

  Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBSTracksRemovedOldMVA);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSTracksRemovedOldMVA);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBSTracksRemovedZNominalOldMVA);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSTracksRemovedZNominalOldMVA);

  Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBSNewMVA);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSNewMVA);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBSZNominalNewMVA);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSZNominalNewMVA);


  Extradist1d.push_back(&polarimetricAcopAnglePVRefitOnlyNoBSTracksRemovedOldMVA);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitOnlyBSTracksRemovedOldMVA);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitOnlyNoBSTracksRemovedZNominalOldMVA);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitOnlyBSTracksRemovedZNominalOldMVA);

  Extradist1d.push_back(&polarimetricAcopAnglePVRefitOnlyNoBSNewMVA);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitOnlyBSNewMVA);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitOnlyNoBSZNominalNewMVA);
  Extradist1d.push_back(&polarimetricAcopAnglePVRefitOnlyBSZNominalNewMVA);

  
  Extradist1d.push_back(&PVXResol);
  Extradist1d.push_back(&PVXNoBSOldResol);
  Extradist1d.push_back(&PVXNoBSTracksRemovedOldResol);
  Extradist1d.push_back(&PVXNoBSNewResol);
  Extradist1d.push_back(&PVXBSOldResol);
  Extradist1d.push_back(&PVXBSTracksRemovedOldResol);
  Extradist1d.push_back(&PVXBSNewResol);

  Extradist1d.push_back(&PVYResol);
  Extradist1d.push_back(&PVYNoBSOldResol);
  Extradist1d.push_back(&PVYNoBSTracksRemovedOldResol);
  Extradist1d.push_back(&PVYNoBSNewResol);
  Extradist1d.push_back(&PVYBSOldResol);
  Extradist1d.push_back(&PVYBSTracksRemovedOldResol);
  Extradist1d.push_back(&PVYBSNewResol);
  
  Extradist1d.push_back(&PVZResol);
  Extradist1d.push_back(&PVZNoBSOldResol);
  Extradist1d.push_back(&PVZNoBSTracksRemovedOldResol);
  Extradist1d.push_back(&PVZNoBSNewResol);
  Extradist1d.push_back(&PVZBSOldResol);
  Extradist1d.push_back(&PVZBSTracksRemovedOldResol);
  Extradist1d.push_back(&PVZBSNewResol);
  
  Extradist1d.push_back(&PVXNoBSTracksRemovedOldOnlyResol);
  Extradist1d.push_back(&PVXNoBSNewOnlyResol);
  Extradist1d.push_back(&PVXBSTracksRemovedOldOnlyResol);
  Extradist1d.push_back(&PVXBSNewOnlyResol);
  
  Extradist1d.push_back(&PVYNoBSTracksRemovedOldOnlyResol);
  Extradist1d.push_back(&PVYNoBSNewOnlyResol);
  Extradist1d.push_back(&PVYBSTracksRemovedOldOnlyResol);
  Extradist1d.push_back(&PVYBSNewOnlyResol);
  
  Extradist1d.push_back(&PVZNoBSTracksRemovedOldOnlyResol);
  Extradist1d.push_back(&PVZNoBSNewOnlyResol);
  Extradist1d.push_back(&PVZBSTracksRemovedOldOnlyResol);
  Extradist1d.push_back(&PVZBSNewOnlyResol);
   
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBSTIPMVA);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSTIPMVA);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBSZNominalTIPMVA);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSZNominalTIPMVA);

  Extradist1d.push_back(&polarimetricAcopAngleTruthA1);
  
  // Extradist1d.push_back(&polarimetricAcopAngleTruthRho);
  // Extradist1d.push_back(&polarimetricAcopAngleTruthPi);
  
  // Extradist1d.push_back(&polarimetricAcopAngle30);
  // Extradist1d.push_back(&polarimetricAcopAngle25);
  // Extradist1d.push_back(&polarimetricAcopAngle20);
  // Extradist1d.push_back(&polarimetricAcopAngle15);
  // Extradist1d.push_back(&polarimetricAcopAngle10);
  // Extradist1d.push_back(&polarimetricAcopAngle5);

  //Extradist1d.push_back(&AcolAngle);
  // Extradist1d.push_back(&AcolAngleSVFit);
  // Extradist1d.push_back(&AcolAngleTruth);

  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBS);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitBS);
  // Extradist1d.push_back(&polarimetricAcopAngleMVA);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBSMVA);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSMVA);

  // Extradist1d.push_back(&polarimetricAcopAngleDecayPlane);

  // Extradist1d.push_back(&polarimetricAcopAnglePtTruthA1);
  // Extradist1d.push_back(&polarimetricAcopAngleMVAPtTruthA1);

  // Extradist1d.push_back(&polarimetricAcopAnglePtTruthRho);
  
  // Extradist1d.push_back(&polarimetricAcopAngleMVAPtTruthRho);

  // Extradist1d.push_back(&polarimetricAcopAnglePtTruthPi);
  // Extradist1d.push_back(&polarimetricAcopAngleMVAPtTruthPi);

   Extradist1d.push_back(&polarimetricAcopAngleSVFit);
   Extradist1d.push_back(&polarimetricAcopAngleMVASVFit);

  // Extradist1d.push_back(&polarimetricAcopAngleSVFitRho);
  // Extradist1d.push_back(&polarimetricAcopAngleMVASVFitRho);
  
  // Extradist1d.push_back(&polarimetricAcopAngleSVFitPi);
  // Extradist1d.push_back(&polarimetricAcopAngleMVASVFitPi);
  
  
  // Extradist1d.push_back(&polarimetricAcopAngleBackground);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBSBackground);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSBackground);
  // Extradist1d.push_back(&polarimetricAcopAngleMVABackground);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBSMVABackground);
  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitBSMVABackground);

  // Extradist1d.push_back(&polarimetricAcopAnglePVRefitNoBSWithoutWSpin);  
  
  Extradist1d.push_back(&test);  
  
  
  Extradist1d.push_back(&PurityDM);
  Extradist1d.push_back(&PurityNewMVA);


  Extradist1d.push_back(&TauSVFitPxResPull);
  Extradist1d.push_back(&TauSVFitPyResPull);
  Extradist1d.push_back(&TauSVFitPzResPull);

  Extradist1d.push_back(&TauPxResPull);
  Extradist1d.push_back(&TauPyResPull);
  Extradist1d.push_back(&TauPzResPull);

  Extradist1d.push_back(&TauSVFitPxResPullMVA);
  Extradist1d.push_back(&TauSVFitPyResPullMVA);
  Extradist1d.push_back(&TauSVFitPzResPullMVA);

  Extradist1d.push_back(&TauPxResPullMVA);
  Extradist1d.push_back(&TauPyResPullMVA);
  Extradist1d.push_back(&TauPzResPullMVA);
  
}

void  HTauTau::doEvent()  { //  Method called on every event
  //cout<<"----------------------"<<endl;
  unsigned int t;                // sample type, you may manage in your further analysis, if needed
  int id(Ntp->GetMCID());  //read event ID of a sample
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}  //  gives a warning if list of samples in Histo.txt  and SkimSummary.log do not coincide 
  //  std::cout<<"------------------ New Event -----------------------"<<std::endl;
  Charge = ChargeSumDummy;
  bool trig=0;
  std::vector<int> TauIndex ;
  std::vector<int> TriggerIndexVector ;
  std::vector<TString>  MatchedTriggerNames;
  // value.at(Trigger)=0;
  // MatchedTriggerNames.push_back("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v");
  // MatchedTriggerNames.push_back("HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v");
  // TriggerIndexVector=Ntp->GetVectorTriggers(MatchedTriggerNames);
  /*for(int jj=0;jj<Ntp->NTriggers();jj++)
    {
    cout<<jj<<" "<<Ntp->TriggerName(jj)<<endl;
    }
    if(Ntp->TriggerAccept(TriggerIndexVector.at(0)))
    {
    if (Ntp->TriggerName(TriggerIndexVector.at(0)).Contains("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v"))
    {
    trig0=1;
    }
    else if (Ntp->TriggerName(TriggerIndexVector.at(0)).Contains("HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v"))
    {
    trig1=1;
    }
    }
    if(Ntp->TriggerAccept(TriggerIndexVector.at(1)))
    {
    if (Ntp->TriggerName(TriggerIndexVector.at(1)).Contains("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v"))
    {
    trig0=1;
    }
    else if (Ntp->TriggerName(TriggerIndexVector.at(1)).Contains("HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v"))
    {
    trig1=1;
    }
    }*/
  // for(unsigned int itrig = 0; itrig < TriggerIndexVector.size(); itrig++){
  //   if(Ntp->TriggerAccept(TriggerIndexVector.at(itrig))){
  //     trig=1;
  //     }
  // }
  // for(unsigned int iDaughter=0;   iDaughter  <  Ntp->NDaughters() ;iDaughter++ ) {
  //   if(Ntp->CHECK_BIT(Ntp->Daughters_trgMatched(iDaughter),29) || Ntp->CHECK_BIT(Ntp->Daughters_trgMatched(iDaughter),27))
  //   {
  // 	TauIndex.push_back(iDaughter);
  // 	}
  // }
  //if(trig && TauIndex.size()>1)value.at(Trigger)=1;
  //pass.at(Trigger)=(value.at(Trigger)==cut.at(Trigger));
  value.at(Id_and_Kin)=0;
  int goodTau_counter=0;
  std::vector<int> thirdLeptonCounter;
  std::vector<int> goodTausIndex;
  for(unsigned int iDaughter=0;   iDaughter  <Ntp->NDaughters()/* TauIndex.size() */ ;iDaughter++ ) {
    if(Ntp->tauBaselineSelection(iDaughter,20., 2.4, 0,0,0)){
      goodTausIndex.push_back(/*TauIndex.at(*/iDaughter/*)*/);
      goodTau_counter++;
    }
  }

  unsigned int Tau1= -1;
  unsigned int Tau2= -1;
  std::vector<int>  Sorted;
  std::vector<int>  PairsIndexTemp;
  std::vector<int>  PairsIndexTau1Temp;
  std::vector<int>  PairsIndexTau2Temp;
  int j=0;
  if(goodTau_counter>1)value.at(Id_and_Kin)=1;
  pass.at(Id_and_Kin)=(value.at(Id_and_Kin) == cut.at(Id_and_Kin));
  //value.at(NPairsFound)=0;

  for(int unsigned ipair =0; ipair <goodTausIndex.size(); ipair++)
    {
      for(int unsigned jpair =1; jpair <goodTausIndex.size(); jpair++)
	{
	  if(jpair>ipair)
	    {
	      // if((Ntp->Daughters_P4(goodTausIndex.at(ipair)).DeltaR(Ntp->Daughters_P4(goodTausIndex.at(jpair))))>0.5){
	      PairsIndexTemp.push_back(j);
	      PairsIndexTau1Temp.push_back(goodTausIndex.at(ipair));
	      PairsIndexTau2Temp.push_back(goodTausIndex.at(jpair));
	      j++;
	      // }
	    }
	}
    }
  //if(PairsIndexTemp.size()>0)value.at(NPairsFound)=1;
  //pass.at(NPairsFound)=(value.at(NPairsFound)==cut.at(NPairsFound));
  if(PairsIndexTemp.size()>0/*pass.at(NPairsFound)*/)
    {
      Sorted = Ntp->SortPair(PairsIndexTemp,PairsIndexTau1Temp,PairsIndexTau2Temp);
      Tau1=PairsIndexTau1Temp.at(Sorted.back());
      Tau2=PairsIndexTau2Temp.at(Sorted.back());
      // value.at(Tau1Isolation)=0;
      // value.at(Tau1Isolation) = (Ntp->isIsolatedTau(Tau1,"Tight"));
      // pass.at(Tau1Isolation) = value.at(Tau1Isolation);
      // value.at(Tau2Isolation)=0;
      // value.at(Tau2Isolation) = (Ntp->isIsolatedTau(Tau2,"Tight"));
      // pass.at(Tau2Isolation) = value.at(Tau2Isolation);
      // value.at(LeptonVeto)=0;
      // for(unsigned int iDaughter=0;   iDaughter  <  Ntp->NDaughters() ;iDaughter++ ) {  // loop over all daughters in the event
      // 	 if((iDaughter!=Tau1)&&(iDaughter!=Tau2)){
      // 	   if(Ntp->ElectronVeto(iDaughter) || Ntp->MuonVeto(iDaughter))thirdLeptonCounter.push_back(iDaughter);
      // 	 }
      // }
      // value.at(LeptonVeto) = thirdLeptonCounter.size()>0;
      // pass.at(LeptonVeto) = (value.at(LeptonVeto)==cut.at(LeptonVeto));

      value.at(PairCharge)=0;
      bool isOS=false;
      isOS=((Ntp->Daughters_charge(Tau1)/abs(Ntp->Daughters_charge(Tau1))) != (Ntp->Daughters_charge(Tau2)/abs(Ntp->Daughters_charge(Tau2))));
      if(isOS)value.at(PairCharge) = 1;
      pass.at(PairCharge) = value.at(PairCharge);
      //     value.at(PairMass) = 999.;
      //     //value.at(MTM) = 999.;
      //     //value.at(MTM) = .;
      //     value.at(PairMass)=(Ntp->Daughters_P4(Tau1)+Ntp->Daughters_P4(Tau2)).M();
      //     pass.at(PairMass) = (value.at(PairMass) < cut.at(PairMass));
      //     //pass.at(MTM) = (value.at(MTM) <= cut.at(MTM));
    }
  // Here you can defined different type of weights you want to apply to events.
  double wobs=1;
  double w=1;
  if(!Ntp->isData() && id!=DataMCType::QCD) {
    //    w *= reweight.weight(2016,26,Ntp->PUNumInteractions());
    w *= reweight.PUweightHTT(Ntp->npu());
    //std::cout<<" pu weigh HTT  "<< reweight.PUweightHTT(Ntp->npu())<<std::endl;
    if(!Ntp->isData() &&PairsIndexTemp.size()>0/* pass.at(NPairsFound)*/ ){
      double w1 = tauTrgSF.getSF(Ntp->TauP4_Corrected(Tau1).Pt(),  Ntp->decayMode(Tau1)) ;  //from Luca
      double w2 = tauTrgSF.getSF(Ntp->TauP4_Corrected(Tau2).Pt(),  Ntp->decayMode(Tau2)) ;
      w*=w1;
      w*=w2;
    }
    if(!Ntp->isData() && PairsIndexTemp.size()>0/*pass.at(NPairsFound)*/ && (id==33 || id == 10110333 || id == 10110433|| id == 10130533|| id ==10210333|| id == 10210433|| id == 10230533|| id ==10310333 || id ==10330533 || id ==10410433 || id == 10410333|| id == 10430533|| id == 30530533 || id ==11 || id ==12)){
      w *= 0.95*0.95;
    }
  }

  TLorentzVector genMomentum(0,0,0,0);
  if(id==33 || id == 10110333 || id == 10110433|| id == 10130533|| id ==10210333|| id == 10210433|| id == 10230533|| id ==10310333 || id ==10330533 || id ==10410433 || id == 10410333|| id == 10430533|| id == 30530533){
    for(unsigned int imc=0; imc < Ntp->NGenParts(); imc++){
      // if((fabs(Ntp->Genpart_pdg(imc)) ==11 || fabs(Ntp->Genpart_pdg(imc)) ==13)   &&  Ntp->CHECK_BIT(Ntp->Genpart_flags(imc),5)  && Ntp->Genpart_status(imc) ==2){
      //   if(Ntp->Genpart_P4(imc).Pt() > 8){
      //     // std::cout<<" pdgid   "<< Ntp->Genpart_pdg(imc)<< "  px  " << Ntp->Genpart_P4(imc).Px() << " status   "<<Ntp->Genpart_status(imc) <<" index  "  << imc<<std::endl;
      //     // if(Ntp->Genpart_ZMothInd(imc)!=-1) 	std::cout<<"Mother    pdgid   "<< Ntp->Genpart_pdg(Ntp->Genpart_ZMothInd(imc))<< "  px  " << Ntp->Genpart_P4(Ntp->Genpart_ZMothInd(imc)).Px() <<std::endl;
      //   }
      // }
      if(fabs(Ntp->Genpart_pdg(imc)) ==15   &&  Ntp->CHECK_BIT(Ntp->Genpart_flags(imc),0)&& Ntp->Genpart_status(imc) ==2) {
	if(Ntp->Genpart_P4(imc).Pt() > 8){
	  genMomentum+=Ntp->Genpart_P4(imc);
	}
      }
    }
  }
  if( id == 30){
    for(unsigned int imc=0; imc < Ntp->NGenParts(); imc++){
      if((fabs(Ntp->Genpart_pdg(imc)) ==11 || fabs(Ntp->Genpart_pdg(imc)) ==13) && Ntp->Genpart_status(imc) ==1  ){
	if(Ntp->Genpart_P4(imc).Pt() > 8){
	  genMomentum+=Ntp->Genpart_P4(imc);
	}
      }
    }
  }
  float zptw(1);
  if(genMomentum.Pt()!=0 && genMomentum.M() > 75 && genMomentum.M() < 120){
    zptw = DataMC_Corr.ZPTWeight(genMomentum.M(),genMomentum.Pt());
  }
  w*=zptw;
  
  //-------------------------  mu/e tau fake rate weights 
  double wAgainstMuon1(1);
  double wAgainstElectron1(1);
  double wAgainstMuon2(1);
  double wAgainstElectron2(1);
  if(id == 33  || id == 10110333 || id == 10110433|| id == 10130533|| id ==10210333|| id == 10210433|| id == 10230533|| id ==10310333 || id ==10330533 || id ==10410433 || id == 10410333|| id == 10430533|| id == 30530533 ||id==11 || id==12){
    if(PairsIndexTemp.size()>0/*pass.at(NPairsFound)*/){
      int matchedIndex1(-1);
      int matchedIndex2(-1);
      double dR1(999);
      double dR2(999);
      for(unsigned int imc=0; imc < Ntp->NGenParts(); imc++){
	if(fabs(Ntp->Genpart_pdg(imc)) ==11 || fabs(Ntp->Genpart_pdg(imc)) ==13)
	  {
	    if(sqrt(pow(Ntp->TauP4_Corrected(Tau1).Phi() - Ntp->Genpart_P4(imc).Phi(),2) + 
		    pow(Ntp->TauP4_Corrected(Tau2).Eta() - Ntp->Genpart_P4(imc).Eta(),2)) < dR1)
	      {
		dR1 = sqrt(pow(Ntp->TauP4_Corrected(Tau1).Phi() - Ntp->Genpart_P4(imc).Phi(),2) + 
			   pow(Ntp->TauP4_Corrected(Tau1).Eta() - Ntp->Genpart_P4(imc).Eta(),2));
		matchedIndex1=imc;
	      }
	    if(sqrt(pow(Ntp->TauP4_Corrected(Tau2).Phi() - Ntp->Genpart_P4(imc).Phi(),2) + 
		    pow(Ntp->TauP4_Corrected(Tau2).Eta() - Ntp->Genpart_P4(imc).Eta(),2)) < dR2)
	      {
		dR2 = sqrt(pow(Ntp->TauP4_Corrected(Tau2).Phi() - Ntp->Genpart_P4(imc).Phi(),2) + 
			   pow(Ntp->TauP4_Corrected(Tau2).Eta() - Ntp->Genpart_P4(imc).Eta(),2));
		matchedIndex2=imc;
	      }
	  }
      }
      if(dR1 < 0.2  && matchedIndex1!=-1 ){
	if(fabs(Ntp->Genpart_pdg(matchedIndex1)) ==13 &&  ( Ntp->CHECK_BIT(Ntp->Genpart_flags(matchedIndex1),0) || Ntp->CHECK_BIT(Ntp->Genpart_flags(matchedIndex1),5)   ) )
	  {
	    wAgainstMuon1 = DataMC_Corr.AgainstMuonDataMCCorrection(Ntp->TauP4_Corrected(Tau1),"AgainstMuonMVATight3");
	  }
	if(fabs(Ntp->Genpart_pdg(matchedIndex1)) ==11 &&  ( Ntp->CHECK_BIT(Ntp->Genpart_flags(matchedIndex1),0) || Ntp->CHECK_BIT(Ntp->Genpart_flags(matchedIndex1),5)   ) )
	  {
	    wAgainstElectron1 = DataMC_Corr.AgainstElectronDataMCCorrection(Ntp->TauP4_Corrected(Tau1),"AgainstElectronMVATight");
	  }
      }
      if(dR2 < 0.2  && matchedIndex2!=-1 ){
	if(fabs(Ntp->Genpart_pdg(matchedIndex2)) ==13 &&  ( Ntp->CHECK_BIT(Ntp->Genpart_flags(matchedIndex2),0) || Ntp->CHECK_BIT(Ntp->Genpart_flags(matchedIndex2),5)   ) )
	  {
	    wAgainstMuon2 = DataMC_Corr.AgainstMuonDataMCCorrection(Ntp->TauP4_Corrected(Tau2),"AgainstMuonMVATight3");
	  }
	if(fabs(Ntp->Genpart_pdg(matchedIndex2)) ==11 &&  ( Ntp->CHECK_BIT(Ntp->Genpart_flags(matchedIndex2),0) || Ntp->CHECK_BIT(Ntp->Genpart_flags(matchedIndex2),5)   ) )
	  {
	    wAgainstElectron2 = DataMC_Corr.AgainstElectronDataMCCorrection(Ntp->TauP4_Corrected(Tau2),"AgainstElectronMVATight");
	  }
      }
    }
  }
  w*=wAgainstMuon1;//cout<<"againstmu1: "<<w<<endl;
  w*=wAgainstElectron1;//cout<<"againstele1: "<<w<<endl;
  w*=wAgainstMuon2;//cout<<"againstmu2: "<<w<<endl;
  w*=wAgainstElectron2;//cout<<"againstele2: "<<w<<endl;
  if(!Ntp->isData() && id!=DataMCType::QCD)w*=Ntp->MC_weight();//cout<<"id: "<<id<<"MC_weight(): "<< Ntp->MC_weight()<<endl; //generator weight because negative weights for this samples
  
  if(id==11)w*=0.2514;
  if(id==12)w*=0.2825;

  w*=Ntp->stitch_weight();


  std::vector<unsigned int> exclude_cuts;
  // exclude_cuts.push_back(Tau1Isolation);
  // exclude_cuts.push_back(Tau2Isolation);
  //exclude_cuts.push_back(PairCharge);
  classic_svFit::LorentzVector tau1P4;
  classic_svFit::LorentzVector tau2P4;

  TLorentzVector Tau1P4;
  TLorentzVector Tau2P4;
  if(passAllBut(exclude_cuts)){
    Tau1P4 = Ntp->TauP4_Corrected(Tau1);
    Tau2P4 = Ntp->TauP4_Corrected(Tau2);
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

  if(passAllBut(exclude_cuts)){
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
    // if(pass.at(PairCharge)) {
    //   if(pass.at(Tau1Isolation) && pass.at(Tau2Isolation) ){
    // 	NQCD.at(t).Fill(1.,w); //A
    //   }
    //   if((Ntp->isIsolatedTau(Tau1,"Medium") && !Ntp->isIsolatedTau(Tau2,"Tight") && Ntp->isIsolatedTau(Tau2,"Loose")) || (Ntp->isIsolatedTau(Tau2,"Medium") && !Ntp->isIsolatedTau(Tau1,"Tight") && Ntp->isIsolatedTau(Tau1,"Loose"))){
    // 	NQCD.at(t).Fill(2.,w); //B
    // 	//TauTauFullMass_B.at(t).Fill((tau1P4+tau2P4).M(),w);
    //   }
    // }
    // if(!pass.at(PairCharge)) {
    //   if(pass.at(Tau1Isolation) && pass.at(Tau2Isolation)){
    // 	NQCD.at(t).Fill(3.,w); //C
    // 	//TauTauFullMass_C.at(t).Fill((tau1P4+tau2P4).M(),w);
    //   }
    //   if((Ntp->isIsolatedTau(Tau1,"Medium") && !Ntp->isIsolatedTau(Tau2,"Tight") && Ntp->isIsolatedTau(Tau2,"Loose")) || (Ntp->isIsolatedTau(Tau2,"Medium") && !Ntp->isIsolatedTau(Tau1,"Tight") && Ntp->isIsolatedTau(Tau1,"Loose"))){
    // 	NQCD.at(t).Fill(4.,w); //D
    // 	//TauTauFullMass_D.at(t).Fill((tau1P4+tau2P4).M(),w);
    //   }
    // }
  }

  bool IsQCDEvent = false;
  // if(passAllBut(exclude_cuts)){
  //   if(pass.at(PairCharge)){
  //     if((Ntp->isIsolatedTau(Tau1,"Medium") && !Ntp->isIsolatedTau(Tau2,"Tight") && Ntp->isIsolatedTau(Tau2,"Loose")) || (Ntp->isIsolatedTau(Tau2,"Medium") && !Ntp->isIsolatedTau(Tau1,"Tight") && Ntp->isIsolatedTau(Tau1,"Loose"))){
  // 	if(id == DataMCType::Data){
  // 	  QCDShape.at(t).Fill(1,w);
  // 	  t=HConfig.GetType(DataMCType::QCD);
  // 	  IsQCDEvent = true;
  // 	}
  //     }
  //   }
  // }

  //if(IsQCDEvent){ pass.at(PairCharge)= true;pass.at(Tau2Isolation)= true;pass.at(Tau1Isolation)=true;}
  /*
    std::vector<unsigned int> exclude_cuts_ForTauIso;
    exclude_cuts_ForTauIso.push_back(Tau1Isolation);
    exclude_cuts_ForTauIso.push_back(Tau2Isolation);
    if(passAllBut(exclude_cuts_ForTauIso)) {
    if(Ntp->isIsolatedTau(Tau1,"Loose"))Tau1isolation.at(t).Fill(1.);
    if(Ntp->isIsolatedTau(Tau1,"Medium"))Tau1isolation.at(t).Fill(2.);
    if(Ntp->isIsolatedTau(Tau1,"Tight"))Tau1isolation.at(t).Fill(3.);
    if(Ntp->isIsolatedTau(Tau1,"VTight"))Tau1isolation.at(t).Fill(4.);
    if(Ntp->isIsolatedTau(Tau2,"Loose"))Tau2isolation.at(t).Fill(1.);
    if(Ntp->isIsolatedTau(Tau2,"Medium"))Tau2isolation.at(t).Fill(2.);
    if(Ntp->isIsolatedTau(Tau2,"Tight"))Tau2isolation.at(t).Fill(3.);
    if(Ntp->isIsolatedTau(Tau2,"VTight"))Tau2isolation.at(t).Fill(4.);
    }
  */
  bool status=AnalysisCuts(t,w,wobs);  // boolean that say whether your event passed critera defined in pass vector. The whole vector must be true for status = true
  ///////////////////////////////////////////////////////////
  // Analyse events which passed selection
  if(status) {

    //cout<<"Refit Size: "<<Ntp->PFTauRefit_PionsP4_Size()<<endl;
    // if(Ntp->PFTauRefit_PionsP4_Size()>0)cout<<"Refit Size Pions: "<<Ntp->PFTauRefit_PionsP4_SizePions(0)<<endl;
    // cout<<"Lepton Hash size: "<<Ntp->NLeptonHash()<<endl;
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
    
    Tau2PT.at(t).Fill(Tau2P4.Pt(),w);
    Tau2E.at(t).Fill(Tau2P4.E(),w);
    Tau2Mass.at(t).Fill(Tau2P4.M(),w);
    Tau2Phi.at(t).Fill(Tau2P4.Phi(),w);
    Tau2Eta.at(t).Fill(Tau2P4.Eta(),w);
    Tau2dz.at(t).Fill(Ntp->dz(Tau2),w);
    Tau2HPSDecayMode.at(t).Fill(Ntp->decayMode(Tau2),w);
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
    for(unsigned int iDaughter=0;   iDaughter  <  Ntp->NDaughters() ;iDaughter++ ) {
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
    TransverseMass.at(t).Fill(Ntp->transverseMass(Tau1P4.Pt(), Tau1P4.Phi(), Tau2P4.Pt(), Tau2P4.Phi()),w);


    int jets_counter=0;
    for(int ijet=0; ijet< Ntp->JetsNumber(); ijet++) {
       
      if((((Ntp->Jet_P4(ijet)).Pt())>20) && (fabs((Ntp->Jet_P4(ijet).Eta())<4.7))) {

	if((((Ntp->Jet_P4(ijet)).DeltaR(Ntp->Daughters_P4(Tau1)))>0.5) && (((Ntp->Jet_P4(ijet)).DeltaR(Ntp->Daughters_P4(Tau2)))>0.5))jets_counter++; {

	  //if(((Ntp->Jet_P4(ijet).Pt())>20) && (fabs((Ntp->Jet_P4(ijet).Eta())<2.4)) && (((Ntp->Jet_P4(ijet))->(Ntp->bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags'))) > 0.800)) {
	  if((((Ntp->jets_nHEF(ijet))<0.99) && ((Ntp->jets_nEmEF(ijet)<0.99)) && ((Ntp->NumConst(ijet))>1)) && (((fabs((Ntp->Jet_P4(ijet).Eta())<=2.4)) && ((Ntp->jets_chHEF(ijet))>0) && ((Ntp->jets_chMult(ijet))>0) && ((Ntp->jets_chEmEF(ijet))<0.99)) || (fabs((Ntp->Jet_P4(ijet).Eta())>2.4))) && (fabs((Ntp->Jet_P4(ijet)).Eta()<=2.7)))jets_counter++;
	  else if(((Ntp->jets_nHEF(ijet))<0.98) && ((Ntp->jets_nEmEF(ijet))>0.01) && ((Ntp->jets_neMult(ijet))>2) && (fabs((Ntp->Jet_P4(ijet).Eta())>2.7)) && (fabs((Ntp->Jet_P4(ijet).Eta())<=3.0)))jets_counter++;
	  else if(((Ntp->jets_nEmEF(ijet))<0.90) && ((Ntp->jets_neMult(ijet))>10) && (fabs((Ntp->Jet_P4(ijet).Eta())>3.0)))jets_counter++;
	}
      }
    }

    NbJets.at(t).Fill(jets_counter,w);
    


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

  
    // TVector3 Tauplus3Dsvfit(Tauplussvfit.X(),Tauplussvfit.Y(),Tauplussvfit.Z());
    // TVector3 Tauminus3Dsvfit(Tauminussvfit.X(),Tauminussvfit.Y(),Tauminussvfit.Z());
    // TVector3 Zsvfit=(Tauminus3Dsvfit).Unit();
    // TVector3 Ysvfit=(Tauminus3Dsvfit.Cross(Tauplus3Dsvfit)).Unit();
    // TVector3 Xsvfit=(Ysvfit.Cross(Tauminus3Dsvfit)).Unit();
    // TVector3 protonsvfit(0,0,1);
    // TVector3 xyprotonsvfit(protonsvfit.Dot(Xsvfit),protonsvfit.Dot(Ysvfit),0);
    // TVector3 xytauplussvfit(Tauplus3Dsvfit.Dot(Xsvfit),Tauplus3Dsvfit.Dot(Ysvfit),0);
    
    
    //   Etasvfit.at(t).Fill(TMath::ATanH((Tauplus3Dsvfit*Zsvfit)/(sqrt((Tauplus3Dsvfit*Zsvfit)*(Tauplus3Dsvfit*Zsvfit)+(Tauplus3Dsvfit*Ysvfit)*(Tauplus3Dsvfit*Ysvfit)+(Tauplus3Dsvfit*Xsvfit)*(Tauplus3Dsvfit*Xsvfit)))),w);
    //   if(id==10310333){	    Phisvfitpipi.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi()),w);}
    //   if(id==10410333){	    Phisvfitpirho.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi()),w);}
    //   if(id==10410433){	    Phisvfitrhorho.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi()),w);}
    //   if(id==10330533){	    Phisvfitpia1.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi()),w);}
    //   if(id==10430533){	    Phisvfitrhoa1.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi()),w);}
    //   if(id==30530533){	    Phisvfita1a1.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi()),w);}
    //   //if(id==10110333 || id==10210333){	    Phisvfitlpi.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi()),w);}
    //   //if(id==10110433 || id==10210433){	    Phisvfitlrho.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi()),w);}
    //   //if(id==10130533 || id==10230533){	    Phisvfitla1.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi()),w);}
    //   Thetasvfit.at(t).Fill(Tauminussvfit.Theta(),w);
    
    // if(Ntp->decayMode(Tau1)==0 && Ntp->decayMode(Tau2)==0 ){PhiDatasvfitpipi.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi()),w);}
    // if((Ntp->decayMode(Tau1)==0 && Ntp->decayMode(Tau2)==1) || (Ntp->decayMode(Tau2)==0 && Ntp->decayMode(Tau1)==1)){PhiDatasvfitpirho.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi()),w);}
    // if((Ntp->decayMode(Tau1)==0 && Ntp->decayMode(Tau2)==10) || (Ntp->decayMode(Tau2)==0 && Ntp->decayMode(Tau1)==10)){PhiDatasvfitpia1.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi()),w);}
    // if(Ntp->decayMode(Tau1)==1 && Ntp->decayMode(Tau2)==1 ){PhiDatasvfitrhorho.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi()),w);}
    // if((Ntp->decayMode(Tau1)==1 && Ntp->decayMode(Tau2)==10) || (Ntp->decayMode(Tau2)==1 && Ntp->decayMode(Tau1)==10)){PhiDatasvfitrhoa1.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi()),w);}
    // if(Ntp->decayMode(Tau1)==10 && Ntp->decayMode(Tau2)==10 ){PhiDatasvfita1a1.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi()),w);}
    

    TVector3 Tauplus3Dvis(Tauplusvis.X(),Tauplusvis.Y(),Tauplusvis.Z());
    TVector3 Tauminus3Dvis(Tauminusvis.X(),Tauminusvis.Y(),Tauminusvis.Z());
    TVector3 Zvis=(Tauminus3Dvis).Unit();
    TVector3 Yvis=(Tauminus3Dvis.Cross(Tauplus3Dvis)).Unit();
    TVector3 Xvis=(Yvis.Cross(Tauminus3Dvis)).Unit();
    TVector3 protonvis(0,0,1);
    TVector3 xyprotonvis(protonvis.Dot(Xvis),protonvis.Dot(Yvis),0);
    TVector3 xytauplusvis(Tauplus3Dvis.Dot(Xvis),Tauplus3Dvis.Dot(Yvis),0);
    

    Etavis.at(t).Fill(TMath::ATanH((Tauplus3Dvis*Zvis)/(sqrt((Tauplus3Dvis*Zvis)*(Tauplus3Dvis*Zvis)+(Tauplus3Dvis*Yvis)*(Tauplus3Dvis*Yvis)+(Tauplus3Dvis*Xvis)*(Tauplus3Dvis*Xvis)))),w);
    if(id==10310333){	    Phivispipi.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),w);}
    if(id==10410333){	    Phivispirho.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),w);}
    if(id==10410433){	    Phivisrhorho.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),w);}
    if(id==10330533){	    Phivispia1.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),w);}
    if(id==10430533){	    Phivisrhoa1.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),w);}
    if(id==30530533){	    Phivisa1a1.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),w);}
    //if(id==10110333 || id==10210333){	    Phivislpi.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),w);}
    // if(id==10110433 || id==10210433){	    Phivislrho.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),w);}
    // if(id==10130533 || id==10230533){	    Phivisla1.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),w);}
    Thetavis.at(t).Fill(Tauminusvis.Theta(),w);

    ZPtVis.at(t).Fill((Tau1P4+Tau2P4).Pt(),w);

    if(id==33 || id == 10110333 || id == 10110433|| id == 10130533|| id ==10210333|| id == 10210433|| id == 10230533|| id ==10310333 || id ==10330533 || id ==10410433 || id == 10410333 || id == 10430533 || id == 30530533) {
      NewPhivsDeltaPhi.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),Ntp->DeltaPhi(protonvis.Phi(),Tauplus3Dvis.Phi()),w);
      NewPhivsDeltaEta.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),protonvis.Eta()-Tauplus3Dvis.Eta(),w);
      NewPhivsPhiproton.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),protonvis.Phi(),w);
      NewPhivsPhiTauplus.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),Tauplus3Dvis.Phi(),w);
      NewPhivsEtaproton.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),protonvis.Eta(),w);
      NewPhivsEtaTauplus.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),Tauplus3Dvis.Eta(),w);

      NewPhiSignal.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),w);
    }
    if(id==DataMCType::QCD || id==DataMCType::Data)NewPhiQCD.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),w);				

    
    // if(id ==11 || id == 12)
    //   {
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
      if((sqrt((Tau1P4.Eta()-Tau1Truth.Eta())*(Tau1P4.Eta()-Tau1Truth.Eta())+(Tau1P4.Phi()-Tau1Truth.Phi())*(Tau1P4.Phi()-Tau1Truth.Phi())))<0.5 && (sqrt((Tau2P4.Eta()-Tau2Truth.Eta())*(Tau2P4.Eta()-Tau2Truth.Eta())+(Tau2P4.Phi()-Tau2Truth.Phi())*(Tau2P4.Phi()-Tau2Truth.Phi())))<0.5)
	{
	  Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau2).E()-Pions2.at(1).E(),w);
	  Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau2).E()-Pions2.at(1).E())/Pions2.at(1).E(),w);
	}
      else{ Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau1).E()-Pions2.at(1).E(),w);Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau1).E()-Pions2.at(1).E())/Pions2.at(1).E(),w);}
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
      if((sqrt((Tau1P4.Eta()-Tau1Truth.Eta())*(Tau1P4.Eta()-Tau1Truth.Eta())+(Tau1P4.Phi()-Tau1Truth.Phi())*(Tau1P4.Phi()-Tau1Truth.Phi())))<0.5 && (sqrt((Tau2P4.Eta()-Tau2Truth.Eta())*(Tau2P4.Eta()-Tau2Truth.Eta())+(Tau2P4.Phi()-Tau2Truth.Phi())*(Tau2P4.Phi()-Tau2Truth.Phi())))<0.5)
	{
	  Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau2).E()-Pions2.at(1).E(),w);
	  Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau2).E()-Pions2.at(1).E())/Pions2.at(1).E(),w);
	}
      else{ Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau1).E()-Pions2.at(1).E(),w);Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau1).E()-Pions2.at(1).E())/Pions2.at(1).E(),w);}
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
      if((sqrt((Tau1P4.Eta()-Pions1.at(1).Eta())*(Tau1P4.Eta()-Pions1.at(1).Eta())+(Tau1P4.Phi()-Pions1.at(1).Phi())*(Tau1P4.Phi()-Pions1.at(1).Phi())))<0.5 && (sqrt((Tau2P4.Eta()-Pions2.at(1).Eta())*(Tau2P4.Eta()-Pions2.at(1).Eta())+(Tau2P4.Phi()-Pions2.at(1).Phi())*(Tau2P4.Phi()-Pions2.at(1).Phi())))<0.5)
	{
	  Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau1).E()-Pions1.at(1).E(),w);
	  Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau2).E()-Pions2.at(1).E(),w);
	  Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau1).E()-Pions1.at(1).E())/Pions1.at(1).E(),w);
	  Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau2).E()-Pions2.at(1).E())/Pions2.at(1).E(),w);
	}
      else
	{
	  Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau1).E()-Pions2.at(1).E(),w);
	  Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau2).E()-Pions1.at(1).E(),w);
	  Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau1).E()-Pions2.at(1).E())/Pions2.at(1).E(),w);
	  Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau2).E()-Pions1.at(1).E())/Pions1.at(1).E(),w);
	}

      //TLorentzVector pi0reco_tau1 =  Ntp->NeutralDaughters_P4(Tau1);
      //TLorentzVector pi0reco_tau2 =  Ntp->NeutralDaughters_P4(Tau2);


      // std::cout<<" pi0reco_tau1   "<< pi0reco_tau1.Phi() << " pi0reco_tau2    " << pi0reco_tau2.Phi() <<std::endl;
      // std::cout<<" pi0mc_tau1   "<< Pions1.at(1).Phi() << " pi0mc_tau2    " << Pions2.at(1).Phi() <<std::endl;
      if(TruthDecayFromTau1==TruthDecayFromTau2)cout<<"Same decay particle 4-4"<<endl;
    }
    if(Ntp->CheckDecayID(4,3)){
      Tau1Truth=Ntp->GetTruthTauLV(4,0);
      Tau2Truth=Ntp->GetTruthTauLV(3,1);
      Pions1=Ntp->GetTruthPionsFromRho(0);
      TruthDecayFromTau1=Pions1.at(0)+Pions1.at(1);
      TruthDecayFromTau2=Ntp->GetTruthTauProductLV(3,211,1);decay=1;
      if((sqrt((Tau1P4.Eta()-Tau1Truth.Eta())*(Tau1P4.Eta()-Tau1Truth.Eta())+(Tau1P4.Phi()-Tau1Truth.Phi())*(Tau1P4.Phi()-Tau1Truth.Phi())))<0.5 && (sqrt((Tau2P4.Eta()-Tau2Truth.Eta())*(Tau2P4.Eta()-Tau2Truth.Eta())+(Tau2P4.Phi()-Tau2Truth.Phi())*(Tau2P4.Phi()-Tau2Truth.Phi())))<0.5)
	{
	  Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau1).E()-Pions1.at(1).E(),w);
	  Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau1).E()-Pions1.at(1).E())/Pions1.at(1).E(),w);
	}
      else{ Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau2).E()-Pions1.at(1).E(),w);Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau2).E()-Pions1.at(1).E())/Pions1.at(1).E(),w);}
      if(TruthDecayFromTau1==TruthDecayFromTau2)cout<<"Same decay particle 4-3"<<endl;

    }
    if(Ntp->CheckDecayID(4,5)){
      Tau1Truth=Ntp->GetTruthTauLV(4,0);
      Tau2Truth=Ntp->GetTruthTauLV(5,1);
      Pions1=Ntp->GetTruthPionsFromRho(0);
      TruthDecayFromTau1=Pions1.at(0)+Pions1.at(1);
      Pions2=Ntp->GetTruthPionsFromA1(1);
      TruthDecayFromTau2=Pions2.at(0)+Pions2.at(1)+Pions2.at(2);decay=1;
      if((sqrt((Tau1P4.Eta()-Tau1Truth.Eta())*(Tau1P4.Eta()-Tau1Truth.Eta())+(Tau1P4.Phi()-Tau1Truth.Phi())*(Tau1P4.Phi()-Tau1Truth.Phi())))<0.5 && (sqrt((Tau2P4.Eta()-Tau2Truth.Eta())*(Tau2P4.Eta()-Tau2Truth.Eta())+(Tau2P4.Phi()-Tau2Truth.Phi())*(Tau2P4.Phi()-Tau2Truth.Phi())))<0.5)
	{
	  Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau1).E()-Pions1.at(1).E(),w);
	  Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau1).E()-Pions1.at(1).E())/Pions1.at(1).E(),w);
	}
      else{ Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau2).E()-Pions1.at(1).E(),w);Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau2).E()-Pions1.at(1).E())/Pions1.at(1).E(),w);}
      if(TruthDecayFromTau1==TruthDecayFromTau2)cout<<"Same decay particle 4-5"<<endl;

    }


    if(Ntp->CheckDecayID(5,5) ){
	  
      Tau1Truth=Ntp->GetTruthTauLV(5,0);
      Tau2Truth=Ntp->GetTruthTauLV(5,1);
      Pions1=Ntp->GetTruthPionsFromA1(0);
      TruthDecayFromTau1=Pions1.at(0)+Pions1.at(1)+Pions1.at(2);
      Pions2=Ntp->GetTruthPionsFromA1(1);
      TruthDecayFromTau2=Pions2.at(0)+Pions2.at(1)+Pions2.at(2);decay=1;
      //if(TruthDecayFromTau1==TruthDecayFromTau2)cout<<"Same decay particle 5-5"<<endl;
    }

    if (decay==1)
      {
	if(Tau1Truth==Tau2Truth)cout<<"Same Taus"<<endl;

	TLorentzVector TruthZ    = Tau1Truth+Tau2Truth;
	// TLorentzVector TruthrhoFromTau2 = (Ntp->GetTruthTauProductLV(3,211) + Ntp->GetTruthTauProductLV(3,111));
	    

	//double visiblePtTruth = (TruthDecayFromTau1 + TruthDecayFromTau2).Pt();

	TLorentzVector VisDecayfromTauplus;
	TLorentzVector VisDecayfromTauminus;
	  
	TauTauTruthMass.at(t).Fill((Tau1Truth+Tau2Truth).M(),w);
	    

	Tauplustruth=Tau2Truth;
	Tauminustruth=Tau1Truth;
	VisDecayfromTauplus=TruthDecayFromTau2;
	VisDecayfromTauminus=TruthDecayFromTau1;
	    
	// if((sqrt((Tauplusvis.Eta()-Tau1Truth.Eta())*(Tauplusvis.Eta()-Tau1Truth.Eta())+(Tauplusvis.Phi()-Tau1Truth.Phi())*(Tauplusvis.Phi()-Tau1Truth.Phi())))<0.5 && (sqrt((Tauminusvis.Eta()-Tau2Truth.Eta())*(Tauminusvis.Eta()-Tau2Truth.Eta())+(Tauminusvis.Phi()-Tau2Truth.Phi())*(Tauminusvis.Phi()-Tau2Truth.Phi())))<0.5)
	//   {
		
		
	
	// 	Tauplustruth=Tau1Truth;
	// 	Tauminustruth=Tau2Truth;
	// 	VisDecayfromTauplus=TruthDecayFromTau1;
	// 	VisDecayfromTauminus=TruthDecayFromTau2;
	//   }
	// else
	//   {
	// 	Tauplustruth=Tau2Truth;
	// 	Tauminustruth=Tau1Truth;
	// 	VisDecayfromTauplus=TruthDecayFromTau2;
	// 	VisDecayfromTauminus=TruthDecayFromTau1;
	//   }
	    
	TVector3 Tauplus3Dtruth(Tauplustruth.X(),Tauplustruth.Y(),Tauplustruth.Z());
	TVector3 Tauminus3Dtruth(Tauminustruth.X(),Tauminustruth.Y(),Tauminustruth.Z());
	TVector3 Ztruth=(Tauminus3Dtruth).Unit();//(Tauminus3Dtruth);//cout<<"Z: "<<Ztruth.X()<<" "<<Ztruth.Y()<<" "<<Ztruth.Z()<<" "<<endl;
	//   Ztruth=Ntp->Rotate(TVector3(0,0,1),Tauminus3Dtruth);
	//  cout<<"Z: "<<Ztruth.X()<<" "<<Ztruth.Y()<<" "<<Ztruth.Z()<<" "<<endl;
	//   Ztruth=Ntp->Rotate(Tauminus3Dtruth,TVector3(0,0,1));
	//   cout<<"Z: "<<Ztruth.X()<<" "<<Ztruth.Y()<<" "<<Ztruth.Z()<<" "<<endl;
	TVector3 Ytruth=(Tauminus3Dtruth.Cross(Tauplus3Dtruth)).Unit();//Ntp->Rotate(TVector3(0,1,0),Ztruth.Cross(Tauplus3Dtruth));
	TVector3 Xtruth=(Ytruth.Cross(Tauminus3Dtruth)).Unit();//Ntp->Rotate(TVector3(1,0,0),Ytruth.Cross(Ztruth));
	//TVector3 newtauminustruth=; //=Ntp->Rotate1(Tauminus3Dtruth,Tauminus3Dtruth,Tauplus3Dtruth);cout<<"newtauminustruth: "<<newtauminustruth.X()<<" "<<newtauminustruth.Y()<<" "<<newtauminustruth.Z()<<endl;//(Tauminus3Dtruth*Xtruth,Tauminus3Dtruth*Ytruth,Tauminus3Dtruth*Ztruth);
	// newtauminustruth =Ntp->Rotate2(newtauminustruth,Tauplus3Dtruth);cout<<"newtauminustruth: "<<newtauminustruth.X()<<" "<<newtauminustruth.Y()<<" "<<newtauminustruth.Z()<<endl;
	//TVector3 newtauplustruth=;//=Ntp->Rotate11(Tauplus3Dtruth,Tauminus3Dtruth);cout<<"newtauplustruth: "<<newtauplustruth.X()<<" "<<newtauplustruth.Y()<<" "<<newtauplustruth.Z()<<endl;//((Tauplus3Dtruth).Dot(Xtruth),(Tauplus3Dtruth).Dot(Ytruth),(Tauplus3Dtruth).Dot(Ztruth));cout<<"finish"<<endl;
	// newtauplustruth=Ntp->Rotate2(newtauplustruth,Tauplus3Dtruth);cout<<"newtauplustruth: "<<newtauplustruth.X()<<" "<<newtauplustruth.Y()<<" "<<newtauplustruth.Z()<<endl;
	TVector3 protontruth(0,0,1);
	TVector3 xyprotontruth(protontruth.Dot(Xtruth),protontruth.Dot(Ytruth),0);
	TVector3 xytauplustruth(Tauplus3Dtruth.Dot(Xtruth),Tauplus3Dtruth.Dot(Ytruth),0);
	    
	Etatruth.at(t).Fill(TMath::ATanH((Tauplus3Dtruth*Ztruth)/(sqrt((Tauplus3Dtruth*Ztruth)*(Tauplus3Dtruth*Ztruth)+(Tauplus3Dtruth*Ytruth)*(Tauplus3Dtruth*Ytruth)+(Tauplus3Dtruth*Xtruth)*(Tauplus3Dtruth*Xtruth)))),w);
	if(id==10310333) {	   
	  Phitruthpipi.at(t).Fill(Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	  //PhiSvFitRespipi.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	  //PhiVisRespipi.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	}
	if(id==10410333) {	  
	  Phitruthpirho.at(t).Fill(Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	  // PhiSvFitRespirho.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	  //PhiVisRespirho.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	}
	if(id==10410433) {	  
	  Phitruthrhorho.at(t).Fill(Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	  //PhiSvFitResrhorho.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	  // PhiVisResrhorho.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	}
	if(id==10330533) {	  
	  Phitruthpia1.at(t).Fill(Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	  // PhiSvFitRespia1.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	  // PhiVisRespia1.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	}
	if(id==10430533) {	  
	  Phitruthrhoa1.at(t).Fill(Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	  // PhiSvFitResrhoa1.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	  // PhiVisResrhoa1.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	}

      }
    //cout<<"Début Pola"<<endl;
    TVector3 tauPrimaryVertex , tauNoBSOldPrimaryVertex, tauBSOldPrimaryVertex, tauNoBSZNominalOldPrimaryVertex, tauBSZNominalOldPrimaryVertex, tauNoBSOneTrackRemovedOldPrimaryVertex, tauBSOneTrackRemovedOldPrimaryVertex, tauNoBSZNominalOneTrackRemovedOldPrimaryVertex, tauBSZNominalOneTrackRemovedOld,  tauNoBSTracksRemovedOldPrimaryVertex, tauBSTracksRemovedOldPrimaryVertex, tauNoBSZNominalTracksRemovedOldPrimaryVertex, tauBSZNominalTracksRemovedOldPrimaryVertex,PrimaryVertex,tauNoBSNewPrimaryVertex, tauBSNewPrimaryVertex, tauNoBSZNominalNewPrimaryVertex, tauBSZNominalNewPrimaryVertex, tauNoBSTIPPrimaryVertex, tauBSTIPPrimaryVertex, tauNoBSZNominalTIPPrimaryVertex, tauBSZNominalTIPPrimaryVertex, TauminusSecondaryVertex , TauplusSecondaryVertex;
    
    TVector3 TauminusDirection , TauplusDirection, TauplusDirectionNoBSOld, TauplusDirectionBSOld, TauplusDirectionNoBSZNominalOld, TauplusDirectionBSZNominalOld, TauplusDirectionNoBSOneTrackRemovedOld, TauplusDirectionBSOneTrackRemovedOld, TauplusDirectionNoBSZNominalOneTrackRemovedOld, TauplusDirectionBSZNominalOneTrackRemovedOld,  TauplusDirectionNoBSTracksRemovedOld, TauplusDirectionBSTracksRemovedOld, TauplusDirectionNoBSZNominalTracksRemovedOld, TauplusDirectionBSZNominalTracksRemovedOld,TauplusDirectionNoBSNew, TauplusDirectionBSNew, TauplusDirectionNoBSZNominalNew, TauplusDirectionBSZNominalNew, TauplusDirectionNoBSTIP, TauplusDirectionBSTIP, TauplusDirectionNoBSZNominalTIP, TauplusDirectionBSZNominalTIP, TauminusDirectionNoBSOld, TauminusDirectionBSOld, TauminusDirectionNoBSZNominalOld, TauminusDirectionBSZNominalOld, TauminusDirectionNoBSOneTrackRemovedOld, TauminusDirectionBSOneTrackRemovedOld, TauminusDirectionNoBSZNominalOneTrackRemovedOld, TauminusDirectionBSZNominalOneTrackRemovedOld, TauminusDirectionNoBSTracksRemovedOld, TauminusDirectionBSTracksRemovedOld, TauminusDirectionNoBSZNominalTracksRemovedOld, TauminusDirectionBSZNominalTracksRemovedOld,TauminusDirectionNoBSNew, TauminusDirectionBSNew, TauminusDirectionNoBSZNominalNew, TauminusDirectionBSZNominalNew, TauminusDirectionNoBSTIP, TauminusDirectionBSTIP, TauminusDirectionNoBSZNominalTIP, TauminusDirectionBSZNominalTIP;

    double thetaGJ_Tauminus , thetaGJ_Tauplus;
    TLorentzVector a1LV_Tauminus , a1LV_Tauplus, a1LVRefit_Tauminus , a1LVRefit_Tauplus,piLV_Tauminus , piLV_Tauplus, rhoLV_Tauminus , rhoLV_Tauplus;
    TLorentzVector TauminusPairConstraint,TauplusPairConstraintNoBSOld, TauplusPairConstraintBSOld, TauplusPairConstraintNoBSZNominalOld, TauplusPairConstraintBSZNominalOld, TauplusPairConstraintNoBSOneTrackRemovedOld, TauplusPairConstraintBSOneTrackRemovedOld, TauplusPairConstraintNoBSZNominalOneTrackRemovedOld, TauplusPairConstraintBSZNominalOneTrackRemovedOld,  TauplusPairConstraintNoBSTracksRemovedOld, TauplusPairConstraintBSTracksRemovedOld, TauplusPairConstraintNoBSZNominalTracksRemovedOld, TauplusPairConstraintBSZNominalTracksRemovedOld,TauplusPairConstraintNoBSNew, TauplusPairConstraintBSNew, TauplusPairConstraintNoBSZNominalNew, TauplusPairConstraintBSZNominalNew, TauplusPairConstraintNoBSTIP, TauplusPairConstraintBSTIP, TauplusPairConstraintNoBSZNominalTIP, TauplusPairConstraintBSZNominalTIP, TauminusPairConstraintNoBSOld, TauminusPairConstraintBSOld, TauminusPairConstraintNoBSZNominalOld, TauminusPairConstraintBSZNominalOld, TauminusPairConstraintNoBSOneTrackRemovedOld, TauminusPairConstraintBSOneTrackRemovedOld, TauminusPairConstraintNoBSZNominalOneTrackRemovedOld, TauminusPairConstraintBSZNominalOneTrackRemovedOld, TauminusPairConstraintNoBSTracksRemovedOld, TauminusPairConstraintBSTracksRemovedOld, TauminusPairConstraintNoBSZNominalTracksRemovedOld, TauminusPairConstraintBSZNominalTracksRemovedOld,TauminusPairConstraintNoBSNew, TauminusPairConstraintBSNew, TauminusPairConstraintNoBSZNominalNew, TauminusPairConstraintBSZNominalNew, TauminusPairConstraintNoBSTIP, TauminusPairConstraintBSTIP, TauminusPairConstraintNoBSZNominalTIP, TauminusPairConstraintBSZNominalTIP;
    TLorentzVector TauplusSmall, TauplusLarge, TauplusMean, TauplusPairConstraint,TauplusPairConstraintMVA,TauplusPairConstraintNoBS,TauplusPairConstraintNoBSMVA,TauplusPairConstraintBS,TauplusPairConstraintBSMVA;
    bool isPlusReal=true, isMinusReal=true, a1X=false, a1a1=false,a1a1MVA=false, a1a1TruthSVFit=false,  a1a1TruthSVFitMVA =false, a1rhoTruthSVFit  =false,  a1rhoTruthSVFitMVA =false, a1piTruthSVFit =false,  a1piTruthSVFitMVA, piminus=false, piplus=false, rhominus=false, rhoplus=false, a1minus=false, a1plus=false,a1minusMVA=false,a1plusMVA=false, piminusMVA=false, piplusMVA=false, rhominusMVA=false, rhoplusMVA=false, a1minusTruthSVFitMVA=false,a1plusTruthSVFitMVA=false, a1minusTruthSVFit=false,a1plusTruthSVFit=false, piminusTruthSVFitMVA=false,piplusTruthSVFitMVA=false, piminusTruthSVFit=false,piplusTruthSVFit=false, rhominusTruthSVFitMVA=false,rhoplusTruthSVFitMVA=false, rhominusTruthSVFit=false,rhoplusTruthSVFit=false;
    std::vector<TLorentzVector> solutions, solutionsNoBSOld, solutionsBSOld, solutionsNoBSZNominalOld, solutionsBSZNominalOld, solutionsNoBSOneTrackRemovedOld, solutionsBSOneTrackRemovedOld, solutionsNoBSOneTrackRemovedZNominalOld, solutionsBSOneTrackRemovedZNominalOld,  solutionsNoBSTracksRemovedOld, solutionsBSTracksRemovedOld, solutionsNoBSZNominalTracksRemovedOld, solutionsBSZNominalTracksRemovedOld,solutionsNoBSNew, solutionsBSNew, solutionsNoBSZNominalNew, solutionsBSZNominalNew, solutionsNoBSTIP, solutionsBSTIP, solutionsNoBSZNominalTIP, solutionsBSZNominalTIP;
    TLorentzVector ZMean,ZMeana1X;
    TLorentzVector ZPairConstraint;

    std::vector<size_t> hashes;
    size_t hash = 0;


    bool hasNoBSMVA=false, hasBSMVA=false;
    bool isRefitNoBSTracksRemovedOld=true;
    bool isRefitBSTracksRemovedOld=true;
    bool isRefitNoBSZNominalTracksRemovedOld=true;
    bool isRefitBSZNominalTracksRemovedOld=true;
    bool isRefitNoBSNew=true;
    bool isRefitBSNew=true;
    bool isRefitNoBSZNominalNew=true;
    bool isRefitBSZNominalNew=true;

    //bool PionsRefitGoodSize=false;

    double Spin_WT=Ntp->TauSpinerGet(TauSpinerInterface::Spin);
    double UnSpin_WT=Ntp->TauSpinerGet(TauSpinerInterface::UnSpin);
    double FlipSpin_WT=Ntp->TauSpinerGet(TauSpinerInterface::FlipSpin);
    double hplus=Ntp->TauSpinerGet(TauSpinerInterface::hplus);
    double hminus=Ntp->TauSpinerGet(TauSpinerInterface::hminus);//1-hplus;
		 
    double Wspin=w*Spin_WT;
    
    //for(int i=0;i<Ntp->LeptonHashSize();i++){cout<<Ntp->LeptonHash(i)<<endl;}
    //   cout<<"P4 size"<<Ntp->Daughters_P4_Size()<<endl;
    //cout<<"Refit Size: "<<Ntp->PFTauRefit_PionsP4_Size()<<endl;

    // cout<<"MVA size: "<<Ntp->MVADM2016Size()<<endl;
    if(Ntp->decayMode(Tauminus) == 0)piminusTruthSVFit=true;
    if(Ntp->decayMode(Tauplus) == 0)piplusTruthSVFit=true;
    if(Ntp->MVADM2016(Tauminus) == 0)piminusTruthSVFitMVA=true;
    if(Ntp->MVADM2016(Tauplus) == 0)piplusTruthSVFitMVA=true;
    if(Ntp->decayMode(Tauminus) == 1)rhominusTruthSVFit=true;
    if(Ntp->decayMode(Tauplus) == 1)rhoplusTruthSVFit=true;
    
    if(Ntp->MVADM2016(Tauminus) == 1)rhominusTruthSVFitMVA=true;
    if(Ntp->MVADM2016(Tauplus) == 1)rhoplusTruthSVFitMVA=true;
    if(Ntp->decayMode(Tauminus) == 10 &&  Ntp->PFTau_hassecondaryVertex(Tauminus) && Ntp->PFtauHasThreePions(Tauminus))a1minus=true;
    if(Ntp->decayMode(Tauplus) == 10 && Ntp->PFTau_hassecondaryVertex(Tauplus) && Ntp->PFtauHasThreePions(Tauplus))a1plus=true;
    //cout<<"MVA size: "<<Ntp->MVADM2016Size()<<endl;
    //cout<<"Tauminus: "<<Tauminus<<"  Tauplus: "<<Tauplus<<endl;
    if(Ntp->decayMode(Tauminus) == 10 &&  Ntp->PFtauHasThreePions(Tauminus))a1minusTruthSVFit=true;
    if(Ntp->decayMode(Tauplus) == 10 && Ntp->PFtauHasThreePions(Tauplus))a1plusTruthSVFit=true;
    if(Ntp->MVADM2016(Tauminus) == 10 &&  Ntp->PFtauHasThreePions(Tauminus))a1minusTruthSVFitMVA=true;
    if(Ntp->MVADM2016(Tauplus) == 10 &&  Ntp->PFtauHasThreePions(Tauplus))a1plusTruthSVFitMVA=true;
    if(Ntp->MVADM2016(Tauminus) == 10 &&  Ntp->PFTau_hassecondaryVertex(Tauminus) && Ntp->PFtauHasThreePions(Tauminus))a1minusMVA=true;
    if(Ntp->MVADM2016(Tauplus) == 10 && Ntp->PFTau_hassecondaryVertex(Tauplus) && Ntp->PFtauHasThreePions(Tauplus))a1plusMVA=true;
    
    
    if((a1minus && !a1plus ) ||( !a1minus && a1plus)) a1X=true;  //a1-X
    if(a1minus && a1plus ) a1a1=true;  //a1-a1
    if(a1minusMVA && a1plusMVA ) a1a1MVA=true;  //a1-a1

    if(a1plusTruthSVFit && a1minusTruthSVFit) a1a1TruthSVFit=true;
    if(a1plusTruthSVFitMVA && a1minusTruthSVFitMVA) a1a1TruthSVFitMVA=true;
    if((a1plusTruthSVFit && rhominusTruthSVFit) || (a1minusTruthSVFit && rhoplusTruthSVFit)) a1rhoTruthSVFit=true;
    if((a1plusTruthSVFitMVA && rhominusTruthSVFitMVA) || (a1minusTruthSVFitMVA && rhoplusTruthSVFitMVA)) a1rhoTruthSVFitMVA=true;
    if((a1plusTruthSVFit && piminusTruthSVFit) || (a1minusTruthSVFit && piplusTruthSVFit)) a1piTruthSVFit=true;
    if((a1plusTruthSVFitMVA && piminusTruthSVFitMVA) || (a1minusTruthSVFitMVA && piplusTruthSVFitMVA)) a1piTruthSVFitMVA=true;
   

    TLorentzVector Tauplussvfit;
    TLorentzVector Tauminussvfit;
    
    ClassicSVfit svfitAlgo1;
    //FastMTT FastMTTAlgo;
    if(a1a1TruthSVFit || a1a1TruthSVFitMVA /*|| a1rhoTruthSVFit || a1rhoTruthSVFitMVA ||a1piTruthSVFit || a1piTruthSVFitMVA*/)
      {
	// // //---------  svfit ---------------------
	std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;
	classic_svFit::MeasuredTauLepton lep1(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Tau1P4.Pt(), Tau1P4.Eta(),  Tau1P4.Phi(), Tau1P4.M(),10);
	classic_svFit::MeasuredTauLepton lep2(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Tau2P4.Pt(), Tau2P4.Eta(),  Tau2P4.Phi(), Tau2P4.M(),10);
	
	measuredTauLeptons.push_back(lep1);
	measuredTauLeptons.push_back(lep2);
	TMatrixD metcov(2,2);
	double metx = Ntp->MET()*cos(Ntp->METphi());
	double mety = Ntp->MET()*sin(Ntp->METphi());
	
	metcov[0][0] = Ntp->PFMETCov00();
	metcov[1][0] = Ntp->PFMETCov01();
	metcov[0][1] = Ntp->PFMETCov10();
	metcov[1][1] = Ntp->PFMETCov11();
	
	svfitAlgo1.setHistogramAdapter(new classic_svFit::TauTauHistogramAdapter());
	svfitAlgo1.addLogM_fixed(true,5.0);
	svfitAlgo1.setDiTauMassConstraint(125.10);
	//FastMTTAlgo.run(measuredTauLeptons, metx, mety, metcov);
	svfitAlgo1.integrate(measuredTauLeptons,metx,mety, metcov );
	if(svfitAlgo1.isValidSolution()){
	  
	  double higgsmass  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svfitAlgo1.getHistogramAdapter())->getMass();
	  h_SVFitMass.at(t).Fill(higgsmass,w); 
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
	  if(svfitAlgo1.isValidSolution()){
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
	}
	//}

	if((a1a1 ||a1a1MVA) && std::isnan(Wspin)!=true && Ntp->PFTauRefit_PionsP4_SizePions(Tauminus)==3 && Ntp->PFTauRefit_PionsP4_SizePions(Tauplus)/* && (id==11 || id==12)*/)
      {
	std::vector<double> PVRefitNoBSTracksRemovedOld_X_temp, PVRefitNoBSTracksRemovedOld_Y_temp, PVRefitNoBSTracksRemovedOld_Z_temp;
	for(unsigned int i=0;i<Ntp->NPVRefitNoBSTracksRemovedOld();i++)
	  {
	    PVRefitNoBSTracksRemovedOld_X_temp.push_back(Ntp->RefitPVNoBSTracksRemovedOld_x(i));
	    PVRefitNoBSTracksRemovedOld_Y_temp.push_back(Ntp->RefitPVNoBSTracksRemovedOld_y(i));
	    PVRefitNoBSTracksRemovedOld_Z_temp.push_back(Ntp->RefitPVNoBSTracksRemovedOld_z(i));
	  }
	
	std::vector<double> PVRefitBSTracksRemovedOld_X_temp, PVRefitBSTracksRemovedOld_Y_temp, PVRefitBSTracksRemovedOld_Z_temp;
	for(unsigned int i=0;i<Ntp->NPVRefitBSTracksRemovedOld();i++)
	  {
	    PVRefitBSTracksRemovedOld_X_temp.push_back(Ntp->RefitPVBSTracksRemovedOld_x(i));
	    PVRefitBSTracksRemovedOld_Y_temp.push_back(Ntp->RefitPVBSTracksRemovedOld_y(i));
	    PVRefitBSTracksRemovedOld_Z_temp.push_back(Ntp->RefitPVBSTracksRemovedOld_z(i));
	  }
	//cout<<"X: "<<Ntp->NPVRefitNoBSNew()<<endl;
	//cout<<"Z: "<<Ntp->NZPVRefitNoBSNew()<<endl;
	std::vector<double> PVRefitNoBSNew_X_temp, PVRefitNoBSNew_Y_temp, PVRefitNoBSNew_Z_temp;
	for(unsigned int i=0;i<Ntp->NPVRefitNoBSNew();i++)
	  {
	    PVRefitNoBSNew_X_temp.push_back(Ntp->RefitPVNoBSNew_x(i));
	    PVRefitNoBSNew_Y_temp.push_back(Ntp->RefitPVNoBSNew_y(i));
	    PVRefitNoBSNew_Z_temp.push_back(Ntp->RefitPVNoBSNew_z(i*2));
	  }
	
	std::vector<double> PVRefitBSNew_X_temp, PVRefitBSNew_Y_temp, PVRefitBSNew_Z_temp;
	for(unsigned int i=0;i<Ntp->NPVRefitBSNew();i++)
	  {
	    PVRefitBSNew_X_temp.push_back(Ntp->RefitPVBSNew_x(i));
	    PVRefitBSNew_Y_temp.push_back(Ntp->RefitPVBSNew_y(i));
	    PVRefitBSNew_Z_temp.push_back(Ntp->RefitPVBSNew_z(i));
	  }
	
	vector<size_t> VertexHashNoBSTracksRemovedOld1_temp, VertexHashNoBSTracksRemovedOld2_temp;
	for(unsigned int i=0;i<Ntp->NVertexHashNoBSTracksRemovedOld();i++)
	  {
	    VertexHashNoBSTracksRemovedOld1_temp.push_back(Ntp->VertexHashNoBSTracksRemovedOld1(i));
	    VertexHashNoBSTracksRemovedOld2_temp.push_back(Ntp->VertexHashNoBSTracksRemovedOld2(i));
	    
	  }
	vector<size_t> VertexHashBSTracksRemovedOld1_temp, VertexHashBSTracksRemovedOld2_temp;
	for(unsigned int i=0;i<Ntp->NVertexHashBSTracksRemovedOld();i++)
	  {
	    VertexHashBSTracksRemovedOld1_temp.push_back(Ntp->VertexHashBSTracksRemovedOld1(i));
	    VertexHashBSTracksRemovedOld2_temp.push_back(Ntp->VertexHashBSTracksRemovedOld2(i));
	    
	  }
	vector<size_t> VertexHashNoBSNew1_temp, VertexHashNoBSNew2_temp;
	for(unsigned int i=0;i<Ntp->NVertexHashNoBSNew();i++)
	  {
	    VertexHashNoBSNew1_temp.push_back(Ntp->VertexHashNoBSNew1(i));
	    VertexHashNoBSNew2_temp.push_back(Ntp->VertexHashNoBSNew2(i));
	    
	  }
	vector<size_t> VertexHashBSNew1_temp, VertexHashBSNew2_temp;
	for(unsigned int i=0;i<Ntp->NVertexHashBSNew();i++)
	  {
	    VertexHashBSNew1_temp.push_back(Ntp->VertexHashBSNew1(i));
	    VertexHashBSNew2_temp.push_back(Ntp->VertexHashBSNew2(i));
	    
	  }

	//cout<<Ntp->LeptonHashSize()<<endl;
	//cout<<Tauminus<<"  "<<Tauplus<<endl;
	//size_t hash = 0;
	//vector<size_t> hashes;
	//cout<<"Taille LeptonHash:  "<<Ntp->LeptonHashSize()<<endl;
	//cout<<Ntp->LeptonHash(Tauminus)<<endl;
	//cout<<Ntp->LeptonHash(Tauplus)<<endl;
	boost::hash_combine(hash, Ntp->LeptonHash(Tauminus));
	//cout<<1<<endl;
	boost::hash_combine(hash, Ntp->LeptonHash(Tauplus));
	
	hashes.push_back(hash);
	//cout<<"HashLepton1: "<<hash<<endl;
	hash = 0;
	
	boost::hash_combine(hash, Ntp->LeptonHash(Tauplus));
	boost::hash_combine(hash, Ntp->LeptonHash(Tauminus));
	//cout<<"HashLepton2: "<<hash<<endl;
	hashes.push_back(hash);
	 

	TauminusSecondaryVertex = Ntp->PFTau_secondaryVertex_pos(Tauminus);
	TauplusSecondaryVertex = Ntp->PFTau_secondaryVertex_pos(Tauplus);

	a1LV_Tauminus = Ntp->PFTau_PionsP4(Tauminus,0) + Ntp->PFTau_PionsP4(Tauminus,1) + Ntp->PFTau_PionsP4(Tauminus,2);
	a1LV_Tauplus = Ntp->PFTau_PionsP4(Tauplus,0) + Ntp->PFTau_PionsP4(Tauplus,1) + Ntp->PFTau_PionsP4(Tauplus,2);

	// cout<<"Refit size: "<<Ntp->PFTauRefit_PionsP4_Size()<<endl;
	// cout<<"Normal size: "<<Ntp->PFTau_PionsP4_Size()<<endl;
	// cout<<"Refit size: "<<Ntp->PFTauRefit_PionsP4_SizePions(Tauplus)<<endl;
	// cout<<"Normal size: "<<Ntp->PFTau_PionsP4_SizePions(Tauplus)<<endl;
	// cout<<"Refit size: "<<Ntp->PFTauRefit_PionsP4_SizePions(Tauminus)<<endl;
	// cout<<"Normal size: "<<Ntp->PFTau_PionsP4_SizePions(Tauminus)<<endl;

	a1LVRefit_Tauminus = Ntp->PFTauRefit_PionsP4(Tauminus,0) + Ntp->PFTauRefit_PionsP4(Tauminus,1) + Ntp->PFTauRefit_PionsP4(Tauminus,2);
	a1LVRefit_Tauplus = Ntp->PFTauRefit_PionsP4(Tauplus,0) + Ntp->PFTauRefit_PionsP4(Tauplus,1) + Ntp->PFTauRefit_PionsP4(Tauplus,2);
	  
	  

	tauPrimaryVertex=Ntp->PVtx();

	tauNoBSTracksRemovedOldPrimaryVertex = GetRefittedPV(hashes, tauPrimaryVertex,PVRefitNoBSTracksRemovedOld_X_temp ,PVRefitNoBSTracksRemovedOld_Y_temp,PVRefitNoBSTracksRemovedOld_Z_temp ,VertexHashNoBSTracksRemovedOld1_temp, VertexHashNoBSTracksRemovedOld2_temp,isRefitNoBSTracksRemovedOld);
	tauBSTracksRemovedOldPrimaryVertex = GetRefittedPV(hashes, tauPrimaryVertex, PVRefitBSTracksRemovedOld_X_temp ,PVRefitBSTracksRemovedOld_Y_temp,PVRefitBSTracksRemovedOld_Z_temp ,VertexHashBSTracksRemovedOld1_temp, VertexHashBSTracksRemovedOld2_temp,isRefitBSTracksRemovedOld);
	tauNoBSZNominalTracksRemovedOldPrimaryVertex = GetRefittedPV(hashes, tauPrimaryVertex, PVRefitNoBSTracksRemovedOld_X_temp ,PVRefitNoBSTracksRemovedOld_Y_temp, Ntp->pv_z() ,VertexHashNoBSTracksRemovedOld1_temp, VertexHashNoBSTracksRemovedOld2_temp,isRefitNoBSZNominalTracksRemovedOld);
	tauBSZNominalTracksRemovedOldPrimaryVertex = GetRefittedPV(hashes, tauPrimaryVertex,PVRefitBSTracksRemovedOld_X_temp ,PVRefitBSTracksRemovedOld_Y_temp , Ntp->pv_z() ,VertexHashBSTracksRemovedOld1_temp, VertexHashBSTracksRemovedOld2_temp,isRefitBSZNominalTracksRemovedOld);

	tauNoBSNewPrimaryVertex = GetRefittedPV(hashes, tauPrimaryVertex,PVRefitNoBSNew_X_temp ,PVRefitNoBSNew_Y_temp,PVRefitNoBSNew_Z_temp ,VertexHashNoBSNew1_temp, VertexHashNoBSNew2_temp,isRefitNoBSNew);
	tauBSNewPrimaryVertex = GetRefittedPV(hashes, tauPrimaryVertex, PVRefitBSNew_X_temp ,PVRefitBSNew_Y_temp,PVRefitBSNew_Z_temp ,VertexHashBSNew1_temp, VertexHashBSNew2_temp,isRefitBSNew);
	tauNoBSZNominalNewPrimaryVertex = GetRefittedPV(hashes, tauPrimaryVertex, PVRefitNoBSNew_X_temp ,PVRefitNoBSNew_Y_temp, Ntp->pv_z() ,VertexHashNoBSNew1_temp, VertexHashNoBSNew2_temp,isRefitNoBSZNominalNew);
	tauBSZNominalNewPrimaryVertex = GetRefittedPV(hashes, tauPrimaryVertex,PVRefitBSNew_X_temp ,PVRefitBSNew_Y_temp , Ntp->pv_z() ,VertexHashBSNew1_temp, VertexHashBSNew2_temp,isRefitBSZNominalNew);


	TauminusDirection = TauminusSecondaryVertex - tauPrimaryVertex;
	TauplusDirection = TauplusSecondaryVertex - tauPrimaryVertex;
	    

	solutionsNoBSOld=tauPairMomentumSolutions(TauminusSecondaryVertex-Ntp->PVRefitNoBSOld(), a1LVRefit_Tauminus, a1LV_Tauminus,isMinusReal, TauplusSecondaryVertex-Ntp->PVRefitNoBSOld(), a1LVRefit_Tauplus,a1LV_Tauplus, isPlusReal,false);
	solutionsBSOld=tauPairMomentumSolutions(TauminusSecondaryVertex-Ntp->PVRefitBSOld(), a1LVRefit_Tauminus, a1LV_Tauminus, isMinusReal, TauplusSecondaryVertex-Ntp->PVRefitBSOld(), a1LVRefit_Tauplus, a1LVRefit_Tauplus, isPlusReal,false);
	solutionsNoBSZNominalOld=tauPairMomentumSolutions(TauminusSecondaryVertex-Ntp->PVRefitNoBSZNominalOld(), a1LVRefit_Tauminus, a1LV_Tauminus, isMinusReal, TauplusSecondaryVertex-Ntp->PVRefitNoBSZNominalOld(), a1LVRefit_Tauplus, a1LVRefit_Tauplus, isPlusReal,false);
	solutionsBSZNominalOld=tauPairMomentumSolutions(TauminusSecondaryVertex-Ntp->PVRefitBSZNominalOld(), a1LVRefit_Tauminus, a1LV_Tauminus, isMinusReal, TauplusSecondaryVertex-Ntp->PVRefitBSZNominalOld(), a1LVRefit_Tauplus, a1LVRefit_Tauplus, isPlusReal,false);
	//cout<<"test0"<<endl;
	// cout<<Ntp->NPVRefitNoBSOneTrackRemovedOld()<<endl;
	// cout<<Ntp->LeptonHashSize()<<endl;
	// solutionsNoBSOneTrackRemovedOld=tauPairMomentumSolutions(TauminusSecondaryVertex-Ntp->PVRefitNoBSOneTrackRemovedOld(Tauminus), a1LVRefit_Tauminus, isMinusReal, TauplusSecondaryVertex-Ntp->PVRefitNoBSOneTrackRemovedOld(Tauplus), a1LVRefit_Tauplus, isPlusReal);
	// cout<<"test0.1"<<endl;
	// solutionsBSOneTrackRemovedOld=tauPairMomentumSolutions(TauminusSecondaryVertex-Ntp->PVRefitBSOneTrackRemovedOld(Tauminus), a1LVRefit_Tauminus, isMinusReal, TauplusSecondaryVertex-Ntp->PVRefitBSOneTrackRemovedOld(Tauplus), a1LVRefit_Tauplus, isPlusReal);
	// solutionsNoBSOneTrackRemovedZNominalOld=tauPairMomentumSolutions(TauminusSecondaryVertex-Ntp->PVRefitNoBSOneTrackRemovedZNominalOld(Tauminus), a1LVRefit_Tauminus, isMinusReal, TauplusSecondaryVertex-Ntp->PVRefitNoBSOneTrackRemovedZNominalOld(Tauplus), a1LVRefit_Tauplus, isPlusReal);	cout<<"test3"<<endl;
	// solutionsBSOneTrackRemovedZNominalOld=tauPairMomentumSolutions(TauminusSecondaryVertex-Ntp->PVRefitBSOneTrackRemovedZNominalOld(Tauminus), a1LVRefit_Tauminus, isMinusReal, TauplusSecondaryVertex-Ntp->PVRefitBSOneTrackRemovedZNominalOld(Tauplus), a1LVRefit_Tauplus, isPlusReal);
	//	    	cout<<"test0.5"<<endl;
	solutionsNoBSTracksRemovedOld=tauPairMomentumSolutions(TauminusSecondaryVertex-tauNoBSTracksRemovedOldPrimaryVertex, a1LVRefit_Tauminus, a1LV_Tauminus, isMinusReal, TauplusSecondaryVertex-tauNoBSTracksRemovedOldPrimaryVertex, a1LVRefit_Tauplus, a1LVRefit_Tauplus, isPlusReal,isRefitNoBSTracksRemovedOld);
	solutionsBSTracksRemovedOld=tauPairMomentumSolutions(TauminusSecondaryVertex-tauBSTracksRemovedOldPrimaryVertex, a1LVRefit_Tauminus, a1LV_Tauminus, isMinusReal, TauplusSecondaryVertex-tauBSTracksRemovedOldPrimaryVertex, a1LVRefit_Tauplus, a1LVRefit_Tauplus, isPlusReal,isRefitBSTracksRemovedOld); 
	solutionsNoBSZNominalTracksRemovedOld=tauPairMomentumSolutions(TauminusSecondaryVertex-tauNoBSZNominalTracksRemovedOldPrimaryVertex, a1LVRefit_Tauminus, a1LV_Tauminus, isMinusReal, TauplusSecondaryVertex-tauNoBSZNominalTracksRemovedOldPrimaryVertex, a1LVRefit_Tauplus, a1LVRefit_Tauplus, isPlusReal,isRefitNoBSZNominalTracksRemovedOld); 
	solutionsBSZNominalTracksRemovedOld=tauPairMomentumSolutions(TauminusSecondaryVertex-tauBSZNominalTracksRemovedOldPrimaryVertex, a1LVRefit_Tauminus, a1LV_Tauminus, isMinusReal, TauplusSecondaryVertex-tauBSZNominalTracksRemovedOldPrimaryVertex, a1LVRefit_Tauplus, a1LVRefit_Tauplus, isPlusReal,isRefitBSZNominalTracksRemovedOld);

	//	    cout<<"test4"<<endl;
	solutionsNoBSNew=tauPairMomentumSolutions(TauminusSecondaryVertex-tauNoBSNewPrimaryVertex, a1LVRefit_Tauminus, a1LV_Tauminus, isMinusReal, TauplusSecondaryVertex-tauNoBSNewPrimaryVertex, a1LVRefit_Tauplus, a1LVRefit_Tauplus, isPlusReal,isRefitNoBSNew); 
	solutionsBSNew=tauPairMomentumSolutions(TauminusSecondaryVertex-tauBSNewPrimaryVertex, a1LVRefit_Tauminus, a1LV_Tauminus, isMinusReal, TauplusSecondaryVertex-tauBSNewPrimaryVertex, a1LVRefit_Tauplus, a1LVRefit_Tauplus, isPlusReal,isRefitBSNew); 
	solutionsNoBSZNominalNew=tauPairMomentumSolutions(TauminusSecondaryVertex-tauNoBSZNominalNewPrimaryVertex, a1LVRefit_Tauminus, a1LV_Tauminus, isMinusReal, TauplusSecondaryVertex-tauNoBSZNominalNewPrimaryVertex, a1LVRefit_Tauplus, a1LVRefit_Tauplus, isPlusReal,isRefitNoBSZNominalNew); 
	solutionsBSZNominalNew=tauPairMomentumSolutions(TauminusSecondaryVertex-tauBSZNominalNewPrimaryVertex, a1LVRefit_Tauminus, a1LV_Tauminus, isMinusReal, TauplusSecondaryVertex-tauBSZNominalNewPrimaryVertex, a1LVRefit_Tauplus, a1LVRefit_Tauplus, isPlusReal,isRefitBSZNominalNew); 
	//	    cout<<"test5"<<endl;
	// solutionsNoBSTIP=tauPairMomentumSolutions(TauminusSecondaryVertex-Ntp->PVRefitNoBSTIP(Tauminus), a1LVRefit_Tauminus, a1LVRefit_Tauminus, isMinusReal, TauplusSecondaryVertex-Ntp->PVRefitNoBSTIP(Tauplus), a1LVRefit_Tauplus, a1LVRefit_Tauplus, isPlusReal); 
	// solutionsBSTIP=tauPairMomentumSolutions(TauminusSecondaryVertex-Ntp->PVRefitBSTIP(Tauminus), a1LVRefit_Tauminus, a1LVRefit_Tauminus, isMinusReal, TauplusSecondaryVertex-Ntp->PVRefitBSTIP(Tauplus), a1LVRefit_Tauplus, a1LVRefit_Tauplus, isPlusReal); 
	// solutionsNoBSZNominalTIP=tauPairMomentumSolutions(TauminusSecondaryVertex-Ntp->PVRefitNoBSZNominalTIP(Tauminus), a1LVRefit_Tauminus, a1LVRefit_Tauminus, isMinusReal, TauplusSecondaryVertex-Ntp->PVRefitNoBSZNominalTIP(Tauplus), a1LVRefit_Tauplus, a1LVRefit_Tauplus, isPlusReal); 
	// solutionsBSZNominalTIP=tauPairMomentumSolutions(TauminusSecondaryVertex-Ntp->PVRefitBSZNominalTIP(Tauminus), a1LVRefit_Tauminus, a1LVRefit_Tauminus, isMinusReal, TauplusSecondaryVertex-Ntp->PVRefitBSZNominalTIP(Tauplus), a1LVRefit_Tauplus, a1LVRefit_Tauplus, isPlusReal);
	//cout<<"test6"<<endl;
	// TauminusDirectionNoBSRefit = TauminusSecondaryVertex - tauNoBSRefitPrimaryVertex;
	// TauplusDirectionNoBSRefit = TauplusSecondaryVertex - tauNoBSRefitPrimaryVertex;
	// TauminusDirectionBSRefit = TauminusSecondaryVertex - tauBSRefitPrimaryVertex;
	// TauplusDirectionBSRefit = TauplusSecondaryVertex - tauBSRefitPrimaryVertex;

	//solutionsNoBSRefit=tauPairMomentumSolutions(TauminusDirectionNoBSRefit, a1LVRefit_Tauminus, isMinusReal, TauplusDirectionNoBSRefit, a1LVRefit_Tauplus, isPlusReal);
	//solutionsBSRefit=tauPairMomentumSolutions(TauminusDirectionBSRefit, a1LVRefit_Tauminus, isMinusReal, TauplusDirectionBSRefit, a1LVRefit_Tauplus, isPlusReal);

	//TauminusPairConstraintNoBS = solutionsNoBSRefit.at(3);
	//TauminusPairConstraintBS = solutionsBSRefit.at(3);
	    
	//TauplusPairConstraintNoBS = solutionsNoBSRefit.at(7);
	//TauplusPairConstraintBS = solutionsBSRefit.at(7);
	    

	//isMinusReal = isPhysical(TauminusDirection, a1LVRefit_Tauminus);
	//isPlusReal = isPhysical(TauplusDirection, a1LVRefit_Tauplus);
	    
	//IsPairPhysical.at(t).Fill(isPlusReal && isMinusReal,w);
	    
	solutions=tauPairMomentumSolutions(TauminusDirection, a1LVRefit_Tauminus, a1LV_Tauminus, isMinusReal, TauplusDirection, a1LVRefit_Tauplus, a1LV_Tauplus, isPlusReal,false);
	  
	// TauminusSmall = solutions.at(0); 
	// TauminusLarge = solutions.at(1); 
	// TauminusMean = solutions.at(2); 
	TauminusPairConstraint = solutions.at(3);
	    
	// TauplusSmall = solutions.at(4);
	// TauplusLarge = solutions.at(5); 
	// TauplusMean = solutions.at(6);
	TauplusPairConstraint = solutions.at(7);

	//cout<<"test7"<<endl;
	TauminusPairConstraintNoBSOld=solutionsNoBSOld.at(3);
	TauminusPairConstraintBSOld=solutionsBSOld.at(3);
	TauminusPairConstraintNoBSZNominalOld=solutionsNoBSZNominalOld.at(3);
	TauminusPairConstraintBSZNominalOld=solutionsBSZNominalOld.at(3);
	// TauminusPairConstraintNoBSOneTrackRemovedOld=solutionsNoBSOneTrackRemovedOld.at(3);
	// TauminusPairConstraintBSOneTrackRemovedOld=solutionsBSOneTrackRemovedOld.at(3);
	// TauminusPairConstraintNoBSZNominalOneTrackRemovedOld=solutionsNoBSOneTrackRemovedZNominalOld.at(3);
	// TauminusPairConstraintBSZNominalOneTrackRemovedOld=solutionsBSOneTrackRemovedZNominalOld.at(3);
	TauminusPairConstraintNoBSTracksRemovedOld=solutionsNoBSTracksRemovedOld.at(3);
	TauminusPairConstraintBSTracksRemovedOld=solutionsBSTracksRemovedOld.at(3);
	TauminusPairConstraintNoBSZNominalTracksRemovedOld=solutionsNoBSZNominalTracksRemovedOld.at(3);
	TauminusPairConstraintBSZNominalTracksRemovedOld=solutionsBSZNominalTracksRemovedOld.at(3);
	TauminusPairConstraintNoBSNew=solutionsNoBSNew.at(3);
	TauminusPairConstraintBSNew=solutionsBSNew.at(3);
	TauminusPairConstraintNoBSZNominalNew=solutionsNoBSZNominalNew.at(3);
	TauminusPairConstraintBSZNominalNew=solutionsBSZNominalNew.at(3);
	// TauminusPairConstraintNoBSTIP=solutionsNoBSTIP.at(3);
	// TauminusPairConstraintBSTIP=solutionsBSTIP.at(3);
	// TauminusPairConstraintNoBSZNominalTIP=solutionsNoBSZNominalTIP.at(3);
	// TauminusPairConstraintBSZNominalTIP=solutionsBSZNominalTIP.at(3);
	// cout<<"test8"<<endl;
	TauplusPairConstraintNoBSOld=solutionsNoBSOld.at(7);
	TauplusPairConstraintBSOld=solutionsBSOld.at(7);
	TauplusPairConstraintNoBSZNominalOld=solutionsNoBSZNominalOld.at(7);
	TauplusPairConstraintBSZNominalOld=solutionsBSZNominalOld.at(7);
	// TauplusPairConstraintNoBSOneTrackRemovedOld=solutionsNoBSOneTrackRemovedOld.at(7);
	// TauplusPairConstraintBSOneTrackRemovedOld=solutionsBSOneTrackRemovedOld.at(7);
	// TauplusPairConstraintNoBSZNominalOneTrackRemovedOld=solutionsNoBSOneTrackRemovedZNominalOld.at(7);
	// TauplusPairConstraintBSZNominalOneTrackRemovedOld=solutionsBSOneTrackRemovedZNominalOld.at(7);
	TauplusPairConstraintNoBSTracksRemovedOld=solutionsNoBSTracksRemovedOld.at(7);
	TauplusPairConstraintBSTracksRemovedOld=solutionsBSTracksRemovedOld.at(7);
	TauplusPairConstraintNoBSZNominalTracksRemovedOld=solutionsNoBSZNominalTracksRemovedOld.at(7);
	TauplusPairConstraintBSZNominalTracksRemovedOld=solutionsBSZNominalTracksRemovedOld.at(7);
	TauplusPairConstraintNoBSNew=solutionsNoBSNew.at(7);
	TauplusPairConstraintBSNew=solutionsBSNew.at(7);
	TauplusPairConstraintNoBSZNominalNew=solutionsNoBSZNominalNew.at(7);
	TauplusPairConstraintBSZNominalNew=solutionsBSZNominalNew.at(7);
	// TauplusPairConstraintNoBSTIP=solutionsNoBSTIP.at(7);
	// TauplusPairConstraintBSTIP=solutionsBSTIP.at(7);
	// TauplusPairConstraintNoBSZNominalTIP=solutionsNoBSZNominalTIP.at(7);
	// TauplusPairConstraintBSZNominalTIP=solutionsBSZNominalTIP.at(7);
	//cout<<"test9"<<endl;

	//if(true/*isPlusReal && isMinusReal*/)
	//{
	// tauminusa1a1MomentumVis.at(t).Fill(Tauminusvis.P(),w);                                  
	// tauplusa1a1MomentumVis.at(t).Fill(Tauplusvis.P(),w);
	// InvariantMasstausa1a1Vis.at(t).Fill((Tauminusvis + Tauplusvis).M(),w);

	// tauminusa1a1MomentumMean.at(t).Fill(TauminusMean.P(),w);
	// tauplusa1a1MomentumMean.at(t).Fill(TauplusMean.P(),w);
	// InvariantMasstausa1a1Mean.at(t).Fill((TauminusMean + TauplusMean).M(),w);
	    
	// tauminusa1a1MomentumPairConstraint.at(t).Fill(TauminusPairConstraint.P(),w);           
	// tauplusa1a1MomentumPairConstraint.at(t).Fill(TauplusPairConstraint.P(),w);
	// InvariantMasstausa1a1PairConstraint.at(t).Fill( (TauminusPairConstraint + TauplusPairConstraint).M() ,w);

	// double widthSolutions_Tauminus = fabs(TauminusLarge.P() - TauminusSmall.P());
	// double widthSolutions_Tauplus = fabs(TauplusLarge.P() - TauplusSmall.P());
	    
	//ZMean=TauminusMean+TauplusMean;
	ZPairConstraint= TauplusPairConstraint+TauminusPairConstraint;
	//}
	// if(a1a1MVA)
	//   {
	//     if (Ntp->isPVtxRefit()) {
	//       // find the vertex among the refitted vertices
	//       for (unsigned int ivertex =0;ivertex<Ntp->NRefitVtx();ivertex++){
	// 	size_t selectionHash = 0;
	// 	boost::hash_combine(selectionHash, Ntp->VertexHashNoBS1(ivertex));
	// 	boost::hash_combine(selectionHash, Ntp->VertexHashNoBS2(ivertex));
	// 	if ( std::find(hashes.begin(), hashes.end(), selectionHash) != hashes.end() ){
	// 	  //cout<<"NoBS MVA Hash Matching!!!"<<endl;
	// 	  tauNoBSRefitMVAPrimaryVertex = Ntp->PVRefitNoBS(ivertex);
	// 	  hasNoBSMVA=true;
	// 	  break;
	// 	}
	//       } // loop over refitted vertices collection
	      
	//       // find the vertex among the refitted vertices calculated w/ beamspot constraint
	//       for (unsigned int ivertex =0;ivertex<Ntp->NRefitVtx();ivertex++){
	// 	size_t selectionHash = 0;
	// 	boost::hash_combine(selectionHash, Ntp->VertexHashBS1(ivertex));
	// 	boost::hash_combine(selectionHash, Ntp->VertexHashBS2(ivertex));
	// 	//cout<<"HashVertex: "<<selectionHash<<endl;
	// 	if ( std::find(hashes.begin(), hashes.end(), selectionHash) != hashes.end() ){
	// 	  //cout<<"BS MVA Hash Matching!!!"<<endl;
	// 	  tauBSRefitMVAPrimaryVertex = Ntp->PVRefitBS(ivertex);
	// 	  hasBSMVA=true;
	// 	  break;
	// 	}
	//       }
	//     } // loop over refitted vertices collection
	    
	//     tauMVAPrimaryVertex = Ntp->PVtx();	    
	//     TauminusDirectionMVA = TauminusSecondaryVertex - tauMVAPrimaryVertex;
	//     TauplusDirectionMVA = TauplusSecondaryVertex - tauMVAPrimaryVertex;
	    
	//     if (hasNoBSMVA || hasBSMVA) {
	//       TauminusDirectionNoBSRefitMVA = TauminusSecondaryVertex - tauNoBSRefitMVAPrimaryVertex;
	//       TauplusDirectionNoBSRefitMVA = TauplusSecondaryVertex - tauNoBSRefitMVAPrimaryVertex;
	//       TauminusDirectionBSRefitMVA = TauminusSecondaryVertex - tauBSRefitMVAPrimaryVertex;
	//       TauplusDirectionBSRefitMVA = TauplusSecondaryVertex - tauBSRefitMVAPrimaryVertex;

	//       solutionsNoBSRefitMVA=tauPairMomentumSolutions(TauminusDirectionNoBSRefitMVA, a1LVRefit_Tauminus, isMinusReal, TauplusDirectionNoBSRefitMVA, a1LVRefit_Tauplus, isPlusReal);
	//       solutionsBSRefitMVA=tauPairMomentumSolutions(TauminusDirectionBSRefitMVA, a1LVRefit_Tauminus, isMinusReal, TauplusDirectionBSRefitMVA, a1LVRefit_Tauplus, isPlusReal);

	//       TauminusPairConstraintNoBSMVA = solutionsNoBSRefitMVA.at(3);
	//       TauminusPairConstraintBSMVA = solutionsBSRefitMVA.at(3);
	    
	//       TauplusPairConstraintNoBSMVA = solutionsNoBSRefitMVA.at(7);
	//       TauplusPairConstraintBSMVA = solutionsBSRefitMVA.at(7);
	//     }
	    
	//     solutionsMVA=tauPairMomentumSolutions(TauminusDirectionMVA, a1LVRefit_Tauminus, isMinusReal, TauplusDirectionMVA, a1LVRefit_Tauplus, isPlusReal);
	    
	//     TauminusPairConstraintMVA = solutionsMVA.at(3);
	//     TauplusPairConstraintMVA = solutionsMVA.at(7);
	    
	//   }
      }
    // if(a1X/* && (id==11 || id==12)*/)
    //   {
    // 	isPlusReal=false;
    // 	isMinusReal=false;
    // 	tauPrimaryVertex = Ntp->PVtx();
    // 	if(a1minus)
    // 	  {
    // 	    isMinusReal=true;//isPhysical(TauminusDirection, a1LV_Tauminus);
    // 	    if(isMinusReal)
    // 	      {
    // 		TauminusSecondaryVertex = Ntp->PFTau_secondaryVertex_pos(Tauminus);
    // 		TauminusDirection = TauminusSecondaryVertex - tauPrimaryVertex;
    // 		a1LV_Tauminus = Ntp->PFTau_PionsP4(Tauminus,0) + Ntp->PFTau_PionsP4(Tauminus,1) + Ntp->PFTau_PionsP4(Tauminus,2);
    // 		TauMinusSolutions = tauMomentumSolutions(TauminusDirection,a1LV_Tauminus,isPlusReal);
    // 		TauminusMean=TauMinusSolutions.at(2);
    // 		InvariantMasstausa1XMean.at(t).Fill((TauminusMean + Tauplusvis).M(),w);
    // 		tauminusa1XMomentumMean.at(t).Fill(TauminusMean.P(),w);
    // 		tauminusa1XMomentumVis.at(t).Fill(Tauminusvis.P(),w);
    // 		//TauplusSecondaryVertex = Ntp->PFTau_secondaryVertex_pos(Tauplus);
    // 		//TauplusDirection = Tauplusvis;
    // 		//a1LV_Tauplus = Ntp->PFTau_PionsP4(Tauplus,0) + Ntp->PFTau_PionsP4(Tauplus,1) + Ntp->PFTau_PionsP4(Tauplus,2);
    // 		//TauPlusSolutions = tauMomentumSolutions(Tauplusvis,a1LV_Tauplus,isMinusReal);
    // 		//TauplusMean=TauPlusSolutions.at(2);
    // 		//TauplusMean.Print();
    // 		//InvariantMasstausa1XMean.at(t).Fill((TauplusMean + Tauminusvis).M(),w);
    // 		//tauplusa1XMomentumMean.at(t).Fill(TauplusDIrection.P(),w);
    // 		tauplusa1XMomentumVis.at(t).Fill(Tauplusvis.P(),w);
    // 		ZMeana1X=Tauplusvis+TauminusMean;
    // 	      }
    // 	  }
    // 	else if(a1plus)
    // 	  {
    // 	    isPlusReal=true; //isPhysical(TauplusDirection, a1LV_Tauplus);
    // 	    if(isPlusReal)
    // 	      {
    // 		//TauminusSecondaryVertex = Ntp->PFTau_secondaryVertex_pos(Tauminus);
    // 		//TauminusDirection = TauminusSecondaryVertex - tauPrimaryVertex;
    // 		//a1LV_Tauminus = Ntp->PFTau_PionsP4(Tauminus,0) + Ntp->PFTau_PionsP4(Tauminus,1) + Ntp->PFTau_PionsP4(Tauminus,2);
    // 		//TauMinusSolutions = tauMomentumSolutions(TauminusDirection,a1LV_Tauminus,isPlusReal);
    // 		//TauminusMean=TauMinusSolutions.at(2);
    // 		//InvariantMasstausa1XMean.at(t).Fill((TauminusMean + Tauplusvis).M(),w);
    // 		//tauminusa1XMomentumMean.at(t).Fill(TauminusMean.P(),w);
    // 		tauminusa1XMomentumVis.at(t).Fill(Tauminusvis.P(),w);
    // 		TauplusSecondaryVertex = Ntp->PFTau_secondaryVertex_pos(Tauplus);
    // 		TauplusDirection = TauplusSecondaryVertex - tauPrimaryVertex;
    // 		a1LV_Tauplus = Ntp->PFTau_PionsP4(Tauplus,0) + Ntp->PFTau_PionsP4(Tauplus,1) + Ntp->PFTau_PionsP4(Tauplus,2);
    // 		TauPlusSolutions = tauMomentumSolutions(TauplusDirection,a1LV_Tauplus,isMinusReal);
    // 		TauplusMean=TauPlusSolutions.at(2);
    // 		InvariantMasstausa1XMean.at(t).Fill((TauplusMean + Tauminusvis).M(),w);
    // 		tauplusa1XMomentumMean.at(t).Fill(TauplusMean.P(),w);
    // 		tauplusa1XMomentumVis.at(t).Fill(Tauplusvis.P(),w);
    // 		ZMeana1X=Tauminusvis+TauplusMean;
    // 	      }
    // 	  }
    // 	InvariantMasstausa1XVis.at(t).Fill((Tauminusvis + Tauplusvis).M(),w);
    //   }
    
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
    // bool gammaa10=false;
    // bool gammaa11=false;
    // bool gammarho0=false;
    // bool gammarho1=false;
    bool purityDM=false;
    bool purityNewMVA=false;

   

    for(int i=0;i<Ntp->NMCTauDecayProducts(0);i++)
      {
	decay0.push_back(Ntp->MCTauandProd_pdgid(0,i));
      }
    for(int i=0;i<Ntp->NMCTauDecayProducts(1);i++)
      {
	decay1.push_back(Ntp->MCTauandProd_pdgid(1,i));
      }

    // for (auto i: decay0){
    //   std::cout << i << ' ';}
    // cout<<endl;
    
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
    
    // if(Ntp->NMCTauDecayProducts(0)==5)a10=abs(Ntp->MCTauandProd_pdgid(0,2))==211 && abs(Ntp->MCTauandProd_pdgid(0,3))==211 && abs(Ntp->MCTauandProd_pdgid(0,4))==211;
    // if(Ntp->NMCTauDecayProducts(1)==5)a11=abs(Ntp->MCTauandProd_pdgid(1,2))==211 && abs(Ntp->MCTauandProd_pdgid(1,3))==211 && abs(Ntp->MCTauandProd_pdgid(1,4))==211;
    // if(Ntp->NMCTauDecayProducts(0)==6)rho0=abs(Ntp->MCTauandProd_pdgid(0,2))==111 && abs(Ntp->MCTauandProd_pdgid(0,3))==22 && abs(Ntp->MCTauandProd_pdgid(0,4))==22 &&  abs(Ntp->MCTauandProd_pdgid(0,5))==211;
    // if(Ntp->NMCTauDecayProducts(1)==6)rho1=abs(Ntp->MCTauandProd_pdgid(1,2))==111 && abs(Ntp->MCTauandProd_pdgid(1,3))==22 && abs(Ntp->MCTauandProd_pdgid(1,4))==22 &&  abs(Ntp->MCTauandProd_pdgid(1,5))==211;
    // if(Ntp->NMCTauDecayProducts(0)==3)pi0=abs(Ntp->MCTauandProd_pdgid(0,2))==211;
    // if(Ntp->NMCTauDecayProducts(1)==3)pi1=abs(Ntp->MCTauandProd_pdgid(1,2))==211;

    // if(Ntp->NMCTauDecayProducts(0)==6)gammaa10=abs(Ntp->MCTauandProd_pdgid(0,5))==22;
    // if(Ntp->NMCTauDecayProducts(1)==6)gammaa11=abs(Ntp->MCTauandProd_pdgid(1,5))==22;
    //gammarho0=abs(Ntp->MCTauandProd_pdgid(0,4))==22;
    //gammarho1=abs(Ntp->MCTauandProd_pdgid(1,4))==22;
    
    
    

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
    
    SCalculator ScalcPVRefitNoBSOld("a1");
    SCalculator ScalcPVRefitBSOld("a1");
    SCalculator ScalcPVRefitNoBSZNominalOld("a1");
    SCalculator ScalcPVRefitBSZNominalOld("a1");

    // SCalculator ScalcPVRefitNoBSOneTrackRemovedOld("a1");
    // SCalculator ScalcPVRefitBSOneTrackRemovedOld("a1");
    // SCalculator ScalcPVRefitNoBSZNominalOneTrackRemovedOld("a1");
    // SCalculator ScalcPVRefitBSZNominalOneTrackRemovedOld("a1");

    SCalculator ScalcPVRefitNoBSTracksRemovedOld("a1");
    SCalculator ScalcPVRefitBSTracksRemovedOld("a1");
    SCalculator ScalcPVRefitNoBSZNominalTracksRemovedOld("a1");
    SCalculator ScalcPVRefitBSZNominalTracksRemovedOld("a1");

    SCalculator ScalcPVRefitNoBSNew("a1");
    SCalculator ScalcPVRefitBSNew("a1");
    SCalculator ScalcPVRefitNoBSZNominalNew("a1");
    SCalculator ScalcPVRefitBSZNominalNew("a1");

    // SCalculator ScalcPVRefitNoBSTIP("a1");
    // SCalculator ScalcPVRefitBSTIP("a1");
    // SCalculator ScalcPVRefitNoBSZNominalTIP("a1");
    // SCalculator ScalcPVRefitBSZNominalTIP("a1");
    
    // SCalculator ScalcMVA("a1");
    // SCalculator ScalcPVRefitNoBSMVA("a1");
    // SCalculator ScalcPVRefitBSMVA("a1");
    // SCalculator ScalcPVRefitNoBSZNominalMVA("a1");
    // SCalculator ScalcPVRefitBSZNominalMVA("a1");

    
  
    // SCalculator ScalcPtTruthA1("a1");
    // SCalculator ScalcMVAPtTruthA1("a1");

    // SCalculator ScalcPtTruthRho1("a1");
    // SCalculator ScalcPtTruthRho2("a1");
    // SCalculator ScalcMVAPtTruthRho1("a1");
    // SCalculator ScalcMVAPtTruthRho2("a1");

    // SCalculator ScalcPtTruthPi1("a1");
    // SCalculator ScalcPtTruthPi2("a1");
    // SCalculator ScalcMVAPtTruthPi1("a1");
    // SCalculator ScalcMVAPtTruthPi2("a1");

    // SCalculator ScalcSVFitA1("a1");
    // SCalculator ScalcMVASVFitA1("a1");

    // SCalculator ScalcSVFitRho1("a1");
    // SCalculator ScalcSVFitRho2("a1");
    // SCalculator ScalcMVASVFitRho1("a1");
    // SCalculator ScalcMVASVFitRho2("a1");
  
    // SCalculator ScalcSVFitPi1("a1");
    // SCalculator ScalcSVFitPi2("a1");
    // SCalculator ScalcMVASVFitPi1("a1");
    // SCalculator ScalcMVASVFitPi2("a1");
  
    // SCalculator ScalcBackground("a1");
    // SCalculator ScalcPVRefitNoBSBackground("a1");
    // SCalculator ScalcPVRefitBSBackground("a1");
    // SCalculator ScalcMVABackground("a1");
    // SCalculator ScalcPVRefitNoBSMVABackground("a1");
    // SCalculator ScalcAllPVRefitBSMVABackground("a1");

    // SCalculator ScalcPVRefitNoBSWithoutWSpin("a1");  

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
  
		    //Scalc1test.SortPions(HadPions_minus, HadPionsCharge_minus);
		    //Scalc2test.SortPions(HadPions_plus, HadPionsCharge_plus);

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

		    //cout<<"Tauminus"<<endl;
		    //for(unsigned int i=0; i<tauandprodminustest.size();i++) {tauandprodminustest.at(i).Print();}
		    //for(unsigned int i=1; i<tauandprodminustest.size();i++){}
		    //cout<<"Tauplus"<<endl;
		    //for(unsigned int i=0; i<tauandprodplustest.size();i++) {tauandprodplustest.at(i).Print();}
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

			//TLorentzVector zeroLV(0,0,0,0);
			//TLorentzVector HadLVMinus(0,0,0,0);
			//TLorentzVector HadLVPlus(0,0,0,0);

			//for (unsigned int i=0;i<sumPionsMinus.size();i++)HadLVMinus+=sumPionsMinus.at(i);
			//for (unsigned int i=0;i<sumPionsPlus.size();i++)HadLVPlus+=sumPionsPlus.at(i);
		
			if(((h1test.Cross(h2test))*(tauminustest_HRF.Vect().Unit()))<=0) test.at(t).Fill(TMath::ATan2((k1test.Cross(k2test)).Mag(),k1test*k2test),Wspin);
			else test.at(t).Fill(2.*TMath::Pi()-TMath::ATan2((k1test.Cross(k2test)).Mag(),k1test*k2test),Wspin);
		      
			    // cout<<"----------------------------------"<<endl;
			    // tauminustest_HRF.Print();
			    // tauplustest_HRF.Print();
			
			polarimetricAcopAngle.at(t).Fill(Scalc.AcopAngle("a1", "a1", TauminusPairConstraint, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraint , HadPions_plus, HadPionsCharge_plus),Wspin);
			    //cout<<"Angle GEF: "<<Scalc.AcopAngle("a1", "a1", TauminusPairConstraint, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraint , HadPions_plus, HadPionsCharge_plus)<<endl;



			    //TLorentzVector testTruth(0,0,0,0);
			if(a1a1)
			  {
				//if((Ntp->GetTruthTauLV(5,0)!=Ntp->GetTruthTauLV(5,1)) && (Ntp->GetTruthTauLV(5,0)!=testTruth) && (Ntp->GetTruthTauLV(5,1)!=testTruth) )
				//{
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
			    
			    if(a1a1 && svfitAlgo1.isValidSolution())
			      {
				if(Tauminussvfit.DeltaR(Tauminustruth)<0.4)
				  {
				    TauSVFitPxResPull.at(t).Fill((Tauminussvfit.X()-Tauminustruth.X())/Tauminustruth.X(),w);
				    TauSVFitPyResPull.at(t).Fill((Tauminussvfit.Y()-Tauminustruth.Y())/Tauminustruth.Y(),w);
				    TauSVFitPzResPull.at(t).Fill((Tauminussvfit.Z()-Tauminustruth.Z())/Tauminustruth.Z(),w);

				    TauSVFitPxResPull.at(t).Fill((Tauplussvfit.X()-Tauplustruth.X())/Tauplustruth.X(),w);
				    TauSVFitPyResPull.at(t).Fill((Tauplussvfit.Y()-Tauplustruth.Y())/Tauplustruth.Y(),w);
				    TauSVFitPzResPull.at(t).Fill((Tauplussvfit.Z()-Tauplustruth.Z())/Tauplustruth.Z(),w);
				  }
				else if(Tauminussvfit.DeltaR(Tauplustruth)<0.4)
				  {
				    TauSVFitPxResPull.at(t).Fill((Tauminussvfit.X()-Tauplustruth.X())/Tauplustruth.X(),w);
				    TauSVFitPyResPull.at(t).Fill((Tauminussvfit.Y()-Tauplustruth.Y())/Tauplustruth.Y(),w);
				    TauSVFitPzResPull.at(t).Fill((Tauminussvfit.Z()-Tauplustruth.Z())/Tauplustruth.Z(),w);
		  
				    TauSVFitPxResPull.at(t).Fill((Tauplussvfit.X()-Tauminustruth.X())/Tauminustruth.X(),w);
				    TauSVFitPyResPull.at(t).Fill((Tauplussvfit.Y()-Tauminustruth.Y())/Tauminustruth.Y(),w);
				    TauSVFitPzResPull.at(t).Fill((Tauplussvfit.Z()-Tauminustruth.Z())/Tauminustruth.Z(),w);
		  
				  }
				if(TauminusPairConstraint.DeltaR(Tauminustruth)<0.4)
				  {
				    TauPxResPull.at(t).Fill((TauminusPairConstraint.X()-Tauminustruth.X())/Tauminustruth.X(),w);
				    TauPyResPull.at(t).Fill((TauminusPairConstraint.Y()-Tauminustruth.Y())/Tauminustruth.Y(),w);
				    TauPzResPull.at(t).Fill((TauminusPairConstraint.Z()-Tauminustruth.Z())/Tauminustruth.Z(),w);

	    
				    TauPxResPull.at(t).Fill((TauplusPairConstraint.X()-Tauplustruth.X())/Tauplustruth.X(),w);
				    TauPyResPull.at(t).Fill((TauplusPairConstraint.Y()-Tauplustruth.Y())/Tauplustruth.Y(),w);
				    TauPzResPull.at(t).Fill((TauplusPairConstraint.Z()-Tauplustruth.Z())/Tauplustruth.Z(),w);

				  }
				else if (TauminusPairConstraint.DeltaR(Tauplustruth)<0.4)
				  {
				    TauPxResPull.at(t).Fill((TauminusPairConstraint.X()-Tauplustruth.X())/Tauplustruth.X(),w);
				    TauPyResPull.at(t).Fill((TauminusPairConstraint.Y()-Tauplustruth.Y())/Tauplustruth.Y(),w);
				    TauPzResPull.at(t).Fill((TauminusPairConstraint.Z()-Tauplustruth.Z())/Tauplustruth.Z(),w);

	    
				    TauPxResPull.at(t).Fill((TauplusPairConstraint.X()-Tauminustruth.X())/Tauminustruth.X(),w);
				    TauPyResPull.at(t).Fill((TauplusPairConstraint.Y()-Tauminustruth.Y())/Tauminustruth.Y(),w);
				    TauPzResPull.at(t).Fill((TauplusPairConstraint.Z()-Tauminustruth.Z())/Tauminustruth.Z(),w);
				  }
			      }
		      }
		  }
		// if(TauminusPairConstraintNoBS!=zeroLV && TauplusPairConstraintNoBS!=zeroLV && ScalcPVRefitNoBS.isOk("a1", "a1", TauminusPairConstraintNoBS, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraintNoBS, HadPions_plus, HadPionsCharge_plus)==true) polarimetricAcopAnglePVRefitNoBS.at(t).Fill(ScalcPVRefitNoBS.AcopAngle("a1", "a1", TauminusPairConstraintNoBS, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraintNoBS, HadPions_plus, HadPionsCharge_plus),Wspin);

		// if(TauminusPairConstraintBS!=zeroLV && TauplusPairConstraintBS!=zeroLV && ScalcPVRefitBS.isOk("a1", "a1", TauminusPairConstraintBS, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraintBS , HadPions_plus, HadPionsCharge_plus)==true) polarimetricAcopAnglePVRefitBS.at(t).Fill(ScalcPVRefitBS.AcopAngle("a1", "a1", TauminusPairConstraintBS, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraintBS , HadPions_plus, HadPionsCharge_plus),Wspin);

		// if(TauminusPairConstraintNoBS!=zeroLV && TauplusPairConstraintNoBS!=zeroLV && ScalcPVRefitNoBSWithoutWSpin.isOk("a1", "a1", TauminusPairConstraintNoBS, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraintNoBS , HadPions_plus, HadPionsCharge_plus)==true) polarimetricAcopAnglePVRefitNoBSWithoutWSpin.at(t).Fill(ScalcPVRefitNoBSWithoutWSpin.AcopAngle("a1", "a1", TauminusPairConstraintNoBS, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraintNoBS , HadPions_plus, HadPionsCharge_plus),w); 

		if(TauminusPairConstraintNoBSOld!=zeroLV && TauplusPairConstraintNoBSOld!=zeroLV && ScalcPVRefitNoBSOld.isOk("a1", "a1", TauminusPairConstraintNoBSOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitNoBSOld.at(t).Fill(ScalcPVRefitNoBSOld.AcopAngle("a1", "a1", TauminusPairConstraintNoBSOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(TauminusPairConstraintBSOld!=zeroLV && TauplusPairConstraintBSOld!=zeroLV && ScalcPVRefitBSOld.isOk("a1", "a1", TauminusPairConstraintBSOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitBSOld.at(t).Fill(ScalcPVRefitBSOld.AcopAngle("a1", "a1", TauminusPairConstraintBSOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(TauminusPairConstraintNoBSZNominalOld!=zeroLV && TauplusPairConstraintNoBSZNominalOld!=zeroLV && ScalcPVRefitNoBSZNominalOld.isOk("a1", "a1", TauminusPairConstraintNoBSZNominalOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominalOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitNoBSZNominalOld.at(t).Fill(ScalcPVRefitNoBSZNominalOld.AcopAngle("a1", "a1", TauminusPairConstraintNoBSZNominalOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominalOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(TauminusPairConstraintBSZNominalOld!=zeroLV && TauplusPairConstraintBSZNominalOld!=zeroLV && ScalcPVRefitBSZNominalOld.isOk("a1", "a1", TauminusPairConstraintBSZNominalOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominalOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitBSZNominalOld.at(t).Fill(ScalcPVRefitBSZNominalOld.AcopAngle("a1", "a1", TauminusPairConstraintBSZNominalOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominalOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);


		
		// if(TauminusPairConstraintNoBSOneTrackRemovedOld!=zeroLV && TauplusPairConstraintNoBSOneTrackRemovedOld!=zeroLV && ScalcPVRefitNoBSOneTrackRemovedOld.isOk("a1", "a1", TauminusPairConstraintNoBSOneTrackRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSOneTrackRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitNoBSOneTrackRemovedOld.at(t).Fill(ScalcPVRefitNoBSOneTrackRemovedOld.AcopAngle("a1", "a1", TauminusPairConstraintNoBSOneTrackRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSOneTrackRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		// if(TauminusPairConstraintBSOneTrackRemovedOld!=zeroLV && TauplusPairConstraintBSOneTrackRemovedOld!=zeroLV && ScalcPVRefitBSOneTrackRemovedOld.isOk("a1", "a1", TauminusPairConstraintBSOneTrackRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSOneTrackRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitBSOneTrackRemovedOld.at(t).Fill(ScalcPVRefitBSOneTrackRemovedOld.AcopAngle("a1", "a1", TauminusPairConstraintBSOneTrackRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSOneTrackRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		// if(TauminusPairConstraintNoBSZNominalOneTrackRemovedOld!=zeroLV && TauplusPairConstraintNoBSZNominalOneTrackRemovedOld!=zeroLV && ScalcPVRefitNoBSZNominalOneTrackRemovedOld.isOk("a1", "a1", TauminusPairConstraintNoBSZNominalOneTrackRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominalOneTrackRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitNoBSOneTrackRemovedZNominalOld.at(t).Fill(ScalcPVRefitNoBSZNominalOneTrackRemovedOld.AcopAngle("a1", "a1", TauminusPairConstraintNoBSZNominalOneTrackRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominalOneTrackRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		// if(TauminusPairConstraintBSZNominalOneTrackRemovedOld!=zeroLV && TauplusPairConstraintBSZNominalOneTrackRemovedOld!=zeroLV && ScalcPVRefitBSZNominalOneTrackRemovedOld.isOk("a1", "a1", TauminusPairConstraintBSZNominalOneTrackRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominalOneTrackRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitBSOneTrackRemovedZNominalOld.at(t).Fill(ScalcPVRefitBSZNominalOneTrackRemovedOld.AcopAngle("a1", "a1", TauminusPairConstraintBSZNominalOneTrackRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominalOneTrackRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);
		// 	cout<<"test12"<<endl;


		if(TauminusPairConstraintNoBSTracksRemovedOld!=zeroLV && TauplusPairConstraintNoBSTracksRemovedOld!=zeroLV && ScalcPVRefitNoBSTracksRemovedOld.isOk("a1", "a1", TauminusPairConstraintNoBSTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitNoBSTracksRemovedOld.at(t).Fill(ScalcPVRefitNoBSTracksRemovedOld.AcopAngle("a1", "a1", TauminusPairConstraintNoBSTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(TauminusPairConstraintBSTracksRemovedOld!=zeroLV && TauplusPairConstraintBSTracksRemovedOld!=zeroLV && ScalcPVRefitBSTracksRemovedOld.isOk("a1", "a1", TauminusPairConstraintBSTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitBSTracksRemovedOld.at(t).Fill(ScalcPVRefitBSTracksRemovedOld.AcopAngle("a1", "a1", TauminusPairConstraintBSTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(TauminusPairConstraintNoBSZNominalTracksRemovedOld!=zeroLV && TauplusPairConstraintNoBSZNominalTracksRemovedOld!=zeroLV && ScalcPVRefitNoBSZNominalTracksRemovedOld.isOk("a1", "a1", TauminusPairConstraintNoBSZNominalTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominalTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitNoBSTracksRemovedZNominalOld.at(t).Fill(ScalcPVRefitNoBSZNominalTracksRemovedOld.AcopAngle("a1", "a1", TauminusPairConstraintNoBSZNominalTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominalTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(TauminusPairConstraintBSZNominalTracksRemovedOld!=zeroLV && TauplusPairConstraintBSZNominalTracksRemovedOld!=zeroLV && ScalcPVRefitBSZNominalTracksRemovedOld.isOk("a1", "a1", TauminusPairConstraintBSZNominalTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominalTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitBSTracksRemovedZNominalOld.at(t).Fill(ScalcPVRefitBSZNominalTracksRemovedOld.AcopAngle("a1", "a1", TauminusPairConstraintBSZNominalTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominalTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		
		if(TauminusPairConstraintNoBSNew!=zeroLV && TauplusPairConstraintNoBSNew!=zeroLV && ScalcPVRefitNoBSNew.isOk("a1", "a1", TauminusPairConstraintNoBSNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSNew, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitNoBSNew.at(t).Fill(ScalcPVRefitNoBSNew.AcopAngle("a1", "a1", TauminusPairConstraintNoBSNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSNew, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(TauminusPairConstraintBSNew!=zeroLV && TauplusPairConstraintBSNew!=zeroLV && ScalcPVRefitBSNew.isOk("a1", "a1", TauminusPairConstraintBSNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSNew, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitBSNew.at(t).Fill(ScalcPVRefitBSNew.AcopAngle("a1", "a1", TauminusPairConstraintBSNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSNew, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(TauminusPairConstraintNoBSZNominalNew!=zeroLV && TauplusPairConstraintNoBSZNominalNew!=zeroLV && ScalcPVRefitNoBSZNominalNew.isOk("a1", "a1", TauminusPairConstraintNoBSZNominalNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominalNew, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitNoBSZNominalNew.at(t).Fill(ScalcPVRefitNoBSZNominalNew.AcopAngle("a1", "a1", TauminusPairConstraintNoBSZNominalNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominalNew, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(TauminusPairConstraintBSZNominalNew!=zeroLV && TauplusPairConstraintBSZNominalNew!=zeroLV && ScalcPVRefitBSZNominalNew.isOk("a1", "a1", TauminusPairConstraintBSZNominalNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominalNew, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitBSZNominalNew.at(t).Fill(ScalcPVRefitBSZNominalNew.AcopAngle("a1", "a1", TauminusPairConstraintBSZNominalNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominalNew, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);
		


		if(isRefitNoBSTracksRemovedOld && TauminusPairConstraintNoBSTracksRemovedOld!=zeroLV && TauplusPairConstraintNoBSTracksRemovedOld!=zeroLV && ScalcPVRefitNoBSTracksRemovedOld.isOk("a1", "a1", TauminusPairConstraintNoBSTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitOnlyNoBSTracksRemovedOld.at(t).Fill(ScalcPVRefitNoBSTracksRemovedOld.AcopAngle("a1", "a1", TauminusPairConstraintNoBSTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(isRefitBSTracksRemovedOld && TauminusPairConstraintBSTracksRemovedOld!=zeroLV && TauplusPairConstraintBSTracksRemovedOld!=zeroLV && ScalcPVRefitBSTracksRemovedOld.isOk("a1", "a1", TauminusPairConstraintBSTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitOnlyBSTracksRemovedOld.at(t).Fill(ScalcPVRefitBSTracksRemovedOld.AcopAngle("a1", "a1", TauminusPairConstraintBSTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(isRefitNoBSZNominalTracksRemovedOld && TauminusPairConstraintNoBSZNominalTracksRemovedOld!=zeroLV && TauplusPairConstraintNoBSZNominalTracksRemovedOld!=zeroLV && ScalcPVRefitNoBSZNominalTracksRemovedOld.isOk("a1", "a1", TauminusPairConstraintNoBSZNominalTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominalTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitOnlyNoBSTracksRemovedZNominalOld.at(t).Fill(ScalcPVRefitNoBSZNominalTracksRemovedOld.AcopAngle("a1", "a1", TauminusPairConstraintNoBSZNominalTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominalTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);


		if(isRefitBSZNominalTracksRemovedOld && TauminusPairConstraintBSZNominalTracksRemovedOld!=zeroLV && TauplusPairConstraintBSZNominalTracksRemovedOld!=zeroLV && ScalcPVRefitBSZNominalTracksRemovedOld.isOk("a1", "a1", TauminusPairConstraintBSZNominalTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominalTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitOnlyBSTracksRemovedZNominalOld.at(t).Fill(ScalcPVRefitBSZNominalTracksRemovedOld.AcopAngle("a1", "a1", TauminusPairConstraintBSZNominalTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominalTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		
		if(isRefitNoBSNew && TauminusPairConstraintNoBSNew!=zeroLV && TauplusPairConstraintNoBSNew!=zeroLV && ScalcPVRefitNoBSNew.isOk("a1", "a1", TauminusPairConstraintNoBSNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSNew, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitOnlyNoBSNew.at(t).Fill(ScalcPVRefitNoBSNew.AcopAngle("a1", "a1", TauminusPairConstraintNoBSNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSNew, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(isRefitBSNew && TauminusPairConstraintBSNew!=zeroLV && TauplusPairConstraintBSNew!=zeroLV && ScalcPVRefitBSNew.isOk("a1", "a1", TauminusPairConstraintBSNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSNew, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitOnlyBSNew.at(t).Fill(ScalcPVRefitBSNew.AcopAngle("a1", "a1", TauminusPairConstraintBSNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSNew, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(isRefitNoBSZNominalNew && TauminusPairConstraintNoBSZNominalNew!=zeroLV && TauplusPairConstraintNoBSZNominalNew!=zeroLV && ScalcPVRefitNoBSZNominalNew.isOk("a1", "a1", TauminusPairConstraintNoBSZNominalNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominalNew, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitOnlyNoBSZNominalNew.at(t).Fill(ScalcPVRefitNoBSZNominalNew.AcopAngle("a1", "a1", TauminusPairConstraintNoBSZNominalNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominalNew, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(isRefitBSZNominalNew && TauminusPairConstraintBSZNominalNew!=zeroLV && TauplusPairConstraintBSZNominalNew!=zeroLV && ScalcPVRefitBSZNominalNew.isOk("a1", "a1", TauminusPairConstraintBSZNominalNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominalNew, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitOnlyBSZNominalNew.at(t).Fill(ScalcPVRefitBSZNominalNew.AcopAngle("a1", "a1", TauminusPairConstraintBSZNominalNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominalNew, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);
		
		// if(TauminusPairConstraintNoBSTIP!=zeroLV && TauplusPairConstraintNoBSTIP!=zeroLV && ScalcPVRefitNoBSTIP.isOk("a1", "a1", TauminusPairConstraintNoBSTIP, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSTIP, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitNoBSTIP.at(t).Fill(ScalcPVRefitNoBSTIP.AcopAngle("a1", "a1", TauminusPairConstraintNoBSTIP, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSTIP, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		// if(TauminusPairConstraintBSTIP!=zeroLV && TauplusPairConstraintBSTIP!=zeroLV && ScalcPVRefitBSTIP.isOk("a1", "a1", TauminusPairConstraintBSTIP, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSTIP, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitBSTIP.at(t).Fill(ScalcPVRefitBSTIP.AcopAngle("a1", "a1", TauminusPairConstraintBSTIP, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSTIP, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		// if(TauminusPairConstraintNoBSZNominalTIP!=zeroLV && TauplusPairConstraintNoBSZNominalTIP!=zeroLV && ScalcPVRefitNoBSZNominalTIP.isOk("a1", "a1", TauminusPairConstraintNoBSZNominalTIP, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominalTIP, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitNoBSZNominalTIP.at(t).Fill(ScalcPVRefitNoBSZNominalTIP.AcopAngle("a1", "a1", TauminusPairConstraintNoBSZNominalTIP, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominalTIP, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		// if(TauminusPairConstraintBSZNominalTIP!=zeroLV && TauplusPairConstraintBSZNominalTIP!=zeroLV && ScalcPVRefitBSZNominalTIP.isOk("a1", "a1", TauminusPairConstraintBSZNominalTIP, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominalTIP, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitBSZNominalTIP.at(t).Fill(ScalcPVRefitBSZNominalTIP.AcopAngle("a1", "a1", TauminusPairConstraintBSZNominalTIP, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominalTIP, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);
		// 	cout<<"test15"<<endl;


		if(a1a1TruthSVFit && svfitAlgo1.isValidSolution() && Scalc.isOk("a1", "a1", Tauminussvfit, HadPions_minus, HadPionsCharge_minus, Tauplussvfit, HadPions_plus, HadPionsCharge_plus)==true){polarimetricAcopAngleSVFit.at(t).Fill(Scalc.AcopAngle("a1", "a1", Tauminussvfit, HadPions_minus, HadPionsCharge_minus, Tauplussvfit, HadPions_plus, HadPionsCharge_plus),Wspin);/*cout<<HadPionsCharge_minus.at(0)<<"   "<<HadPionsCharge_minus.at(1)<<"   "<<HadPionsCharge_minus.at(2)<<endl;cout<<HadPionsCharge_plus.at(0)<<"   "<<HadPionsCharge_plus.at(1)<<"   "<<HadPionsCharge_plus.at(2)<<endl;*/}
	
	       }
	      
	    if(a1a1MVA)
	      {	      
	    	if(TauminusPairConstraint!=zeroLV && TauplusPairConstraint!=zeroLV && Scalc.isOk("a1", "a1", TauminusPairConstraint, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraint , HadPions_plus, HadPionsCharge_plus)==true) polarimetricAcopAngleMVA.at(t).Fill(Scalc.AcopAngle("a1", "a1", TauminusPairConstraint, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraint, HadPions_plus, HadPionsCharge_plus),Wspin);
	    	// if(TauminusPairConstraintNoBSMVA!=zeroLV && TauplusPairConstraintNoBSMVA!=zeroLV && ScalcPVRefitNoBSMVA.isOk("a1", "a1", TauminusPairConstraintNoBSMVA, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraintNoBSMVA , HadPions_plus, HadPionsCharge_plus)==true) polarimetricAcopAnglePVRefitNoBSMVA.at(t).Fill(ScalcPVRefitNoBSMVA.AcopAngle("a1", "a1", TauminusPairConstraintNoBSMVA, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraintNoBSMVA , HadPions_plus, HadPionsCharge_plus),Wspin);
	    	// if(TauminusPairConstraintBSMVA!=zeroLV && TauplusPairConstraintBSMVA!=zeroLV && ScalcPVRefitBSMVA.isOk("a1", "a1", TauminusPairConstraintBSMVA, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraintBSMVA, HadPions_plus, HadPionsCharge_plus)==true) polarimetricAcopAnglePVRefitBSMVA.at(t).Fill(ScalcPVRefitBSMVA.AcopAngle("a1", "a1", TauminusPairConstraintBSMVA, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraintBSMVA, HadPions_plus, HadPionsCharge_plus),Wspin);


		
		if(TauminusPairConstraintNoBSOld!=zeroLV && TauplusPairConstraintNoBSOld!=zeroLV && ScalcPVRefitNoBSOld.isOk("a1", "a1", TauminusPairConstraintNoBSOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitNoBSOldMVA.at(t).Fill(ScalcPVRefitNoBSOld.AcopAngle("a1", "a1", TauminusPairConstraintNoBSOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(TauminusPairConstraintBSOld!=zeroLV && TauplusPairConstraintBSOld!=zeroLV && ScalcPVRefitBSOld.isOk("a1", "a1", TauminusPairConstraintBSOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitBSOldMVA.at(t).Fill(ScalcPVRefitBSOld.AcopAngle("a1", "a1", TauminusPairConstraintBSOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(TauminusPairConstraintNoBSZNominalOld!=zeroLV && TauplusPairConstraintNoBSZNominalOld!=zeroLV && ScalcPVRefitNoBSZNominalOld.isOk("a1", "a1", TauminusPairConstraintNoBSZNominalOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominalOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitNoBSZNominalOldMVA.at(t).Fill(ScalcPVRefitNoBSZNominalOld.AcopAngle("a1", "a1", TauminusPairConstraintNoBSZNominalOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominalOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(TauminusPairConstraintBSZNominalOld!=zeroLV && TauplusPairConstraintBSZNominalOld!=zeroLV && ScalcPVRefitBSZNominalOld.isOk("a1", "a1", TauminusPairConstraintBSZNominalOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominalOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitBSZNominalOldMVA.at(t).Fill(ScalcPVRefitBSZNominalOld.AcopAngle("a1", "a1", TauminusPairConstraintBSZNominalOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominalOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		
		
		// if(TauminusPairConstraintNoBSOneTrackRemovedOld!=zeroLV && TauplusPairConstraintNoBSOneTrackRemovedOld!=zeroLV && ScalcPVRefitNoBSOneTrackRemovedOld.isOk("a1", "a1", TauminusPairConstraintNoBSOneTrackRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSOneTrackRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitNoBSOneTrackRemovedOldMVA.at(t).Fill(ScalcPVRefitNoBSOneTrackRemovedOld.AcopAngle("a1", "a1", TauminusPairConstraintNoBSOneTrackRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSOneTrackRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		// if(TauminusPairConstraintBSOneTrackRemovedOld!=zeroLV && TauplusPairConstraintBSOneTrackRemovedOld!=zeroLV && ScalcPVRefitBSOneTrackRemovedOld.isOk("a1", "a1", TauminusPairConstraintBSOneTrackRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSOneTrackRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitBSOneTrackRemovedOldMVA.at(t).Fill(ScalcPVRefitBSOneTrackRemovedOld.AcopAngle("a1", "a1", TauminusPairConstraintBSOneTrackRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSOneTrackRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		// if(TauminusPairConstraintNoBSZNominalOneTrackRemovedOld!=zeroLV && TauplusPairConstraintNoBSZNominalOneTrackRemovedOld!=zeroLV && ScalcPVRefitNoBSZNominalOneTrackRemovedOld.isOk("a1", "a1", TauminusPairConstraintNoBSZNominalOneTrackRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominalOneTrackRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitNoBSOneTrackRemovedZNominalOldMVA.at(t).Fill(ScalcPVRefitNoBSZNominalOneTrackRemovedOld.AcopAngle("a1", "a1", TauminusPairConstraintNoBSZNominalOneTrackRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominalOneTrackRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		// if(TauminusPairConstraintBSZNominalOneTrackRemovedOld!=zeroLV && TauplusPairConstraintBSZNominalOneTrackRemovedOld!=zeroLV && ScalcPVRefitBSZNominalOneTrackRemovedOld.isOk("a1", "a1", TauminusPairConstraintBSZNominalOneTrackRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominalOneTrackRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitBSOneTrackRemovedZNominalOldMVA.at(t).Fill(ScalcPVRefitBSZNominalOneTrackRemovedOld.AcopAngle("a1", "a1", TauminusPairConstraintBSZNominalOneTrackRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominalOneTrackRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);



		if(TauminusPairConstraintNoBSTracksRemovedOld!=zeroLV && TauplusPairConstraintNoBSTracksRemovedOld!=zeroLV && ScalcPVRefitNoBSTracksRemovedOld.isOk("a1", "a1", TauminusPairConstraintNoBSTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitNoBSTracksRemovedOldMVA.at(t).Fill(ScalcPVRefitNoBSTracksRemovedOld.AcopAngle("a1", "a1", TauminusPairConstraintNoBSTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(TauminusPairConstraintBSTracksRemovedOld!=zeroLV && TauplusPairConstraintBSTracksRemovedOld!=zeroLV && ScalcPVRefitBSTracksRemovedOld.isOk("a1", "a1", TauminusPairConstraintBSTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitBSTracksRemovedOldMVA.at(t).Fill(ScalcPVRefitBSTracksRemovedOld.AcopAngle("a1", "a1", TauminusPairConstraintBSTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(TauminusPairConstraintNoBSZNominalTracksRemovedOld!=zeroLV && TauplusPairConstraintNoBSZNominalTracksRemovedOld!=zeroLV && ScalcPVRefitNoBSZNominalTracksRemovedOld.isOk("a1", "a1", TauminusPairConstraintNoBSZNominalTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominalTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitNoBSTracksRemovedZNominalOldMVA.at(t).Fill(ScalcPVRefitNoBSZNominalTracksRemovedOld.AcopAngle("a1", "a1", TauminusPairConstraintNoBSZNominalTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominalTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(TauminusPairConstraintBSZNominalTracksRemovedOld!=zeroLV && TauplusPairConstraintBSZNominalTracksRemovedOld!=zeroLV && ScalcPVRefitBSZNominalTracksRemovedOld.isOk("a1", "a1", TauminusPairConstraintBSZNominalTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominalTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitBSTracksRemovedZNominalOldMVA.at(t).Fill(ScalcPVRefitBSZNominalTracksRemovedOld.AcopAngle("a1", "a1", TauminusPairConstraintBSZNominalTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominalTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);
		
		
		if(TauminusPairConstraintNoBSNew!=zeroLV && TauplusPairConstraintNoBSNew!=zeroLV && ScalcPVRefitNoBSNew.isOk("a1", "a1", TauminusPairConstraintNoBSNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSNew, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitNoBSNewMVA.at(t).Fill(ScalcPVRefitNoBSNew.AcopAngle("a1", "a1", TauminusPairConstraintNoBSNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSNew, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(TauminusPairConstraintBSNew!=zeroLV && TauplusPairConstraintBSNew!=zeroLV && ScalcPVRefitBSNew.isOk("a1", "a1", TauminusPairConstraintBSNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSNew, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitBSNewMVA.at(t).Fill(ScalcPVRefitBSNew.AcopAngle("a1", "a1", TauminusPairConstraintBSNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSNew, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(TauminusPairConstraintNoBSZNominalNew!=zeroLV && TauplusPairConstraintNoBSZNominalNew!=zeroLV && ScalcPVRefitNoBSZNominalNew.isOk("a1", "a1", TauminusPairConstraintNoBSZNominalNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominalNew, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitNoBSZNominalNewMVA.at(t).Fill(ScalcPVRefitNoBSZNominalNew.AcopAngle("a1", "a1", TauminusPairConstraintNoBSZNominalNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominalNew, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(TauminusPairConstraintBSZNominalNew!=zeroLV && TauplusPairConstraintBSZNominalNew!=zeroLV && ScalcPVRefitBSZNominalNew.isOk("a1", "a1", TauminusPairConstraintBSZNominalNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominalNew, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitBSZNominalNewMVA.at(t).Fill(ScalcPVRefitBSZNominalNew.AcopAngle("a1", "a1", TauminusPairConstraintBSZNominalNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominalNew, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);



		if(isRefitNoBSTracksRemovedOld && TauminusPairConstraintNoBSTracksRemovedOld!=zeroLV && TauplusPairConstraintNoBSTracksRemovedOld!=zeroLV && ScalcPVRefitNoBSTracksRemovedOld.isOk("a1", "a1", TauminusPairConstraintNoBSTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitOnlyNoBSTracksRemovedOldMVA.at(t).Fill(ScalcPVRefitNoBSTracksRemovedOld.AcopAngle("a1", "a1", TauminusPairConstraintNoBSTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(isRefitBSTracksRemovedOld && TauminusPairConstraintBSTracksRemovedOld!=zeroLV && TauplusPairConstraintBSTracksRemovedOld!=zeroLV && ScalcPVRefitBSTracksRemovedOld.isOk("a1", "a1", TauminusPairConstraintBSTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitOnlyBSTracksRemovedOldMVA.at(t).Fill(ScalcPVRefitBSTracksRemovedOld.AcopAngle("a1", "a1", TauminusPairConstraintBSTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(isRefitNoBSZNominalTracksRemovedOld && TauminusPairConstraintNoBSZNominalTracksRemovedOld!=zeroLV && TauplusPairConstraintNoBSZNominalTracksRemovedOld!=zeroLV && ScalcPVRefitNoBSZNominalTracksRemovedOld.isOk("a1", "a1", TauminusPairConstraintNoBSZNominalTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominalTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitOnlyNoBSTracksRemovedZNominalOldMVA.at(t).Fill(ScalcPVRefitNoBSZNominalTracksRemovedOld.AcopAngle("a1", "a1", TauminusPairConstraintNoBSZNominalTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominalTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);


		if(isRefitBSZNominalTracksRemovedOld && TauminusPairConstraintBSZNominalTracksRemovedOld!=zeroLV && TauplusPairConstraintBSZNominalTracksRemovedOld!=zeroLV && ScalcPVRefitBSZNominalTracksRemovedOld.isOk("a1", "a1", TauminusPairConstraintBSZNominalTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominalTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitOnlyBSTracksRemovedZNominalOldMVA.at(t).Fill(ScalcPVRefitBSZNominalTracksRemovedOld.AcopAngle("a1", "a1", TauminusPairConstraintBSZNominalTracksRemovedOld, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominalTracksRemovedOld, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);
		
		
		if(isRefitNoBSNew && TauminusPairConstraintNoBSNew!=zeroLV && TauplusPairConstraintNoBSNew!=zeroLV && ScalcPVRefitNoBSNew.isOk("a1", "a1", TauminusPairConstraintNoBSNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSNew, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitOnlyNoBSNewMVA.at(t).Fill(ScalcPVRefitNoBSNew.AcopAngle("a1", "a1", TauminusPairConstraintNoBSNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSNew, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(isRefitBSNew && TauminusPairConstraintBSNew!=zeroLV && TauplusPairConstraintBSNew!=zeroLV && ScalcPVRefitBSNew.isOk("a1", "a1", TauminusPairConstraintBSNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSNew, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitOnlyBSNewMVA.at(t).Fill(ScalcPVRefitBSNew.AcopAngle("a1", "a1", TauminusPairConstraintBSNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSNew, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(isRefitNoBSZNominalNew && TauminusPairConstraintNoBSZNominalNew!=zeroLV && TauplusPairConstraintNoBSZNominalNew!=zeroLV && ScalcPVRefitNoBSZNominalNew.isOk("a1", "a1", TauminusPairConstraintNoBSZNominalNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominalNew, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitOnlyNoBSZNominalNewMVA.at(t).Fill(ScalcPVRefitNoBSZNominalNew.AcopAngle("a1", "a1", TauminusPairConstraintNoBSZNominalNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominalNew, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		if(isRefitBSZNominalNew && TauminusPairConstraintBSZNominalNew!=zeroLV && TauplusPairConstraintBSZNominalNew!=zeroLV && ScalcPVRefitBSZNominalNew.isOk("a1", "a1", TauminusPairConstraintBSZNominalNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominalNew, HadRefitPions_plus, HadRefitPionsCharge_plus)==true){ polarimetricAcopAnglePVRefitOnlyBSZNominalNewMVA.at(t).Fill(ScalcPVRefitBSZNominalNew.AcopAngle("a1", "a1", TauminusPairConstraintBSZNominalNew, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominalNew, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		  
		  if(a1a1TruthSVFitMVA && svfitAlgo1.isValidSolution() && Scalc.isOk("a1", "a1", Tauminussvfit, HadPions_minus, HadPionsCharge_minus,  Tauplussvfit, HadPions_plus, HadPionsCharge_plus)==true) polarimetricAcopAngleMVASVFit.at(t).Fill(Scalc.AcopAngle("a1", "a1", Tauminussvfit, HadPions_minus, HadPionsCharge_minus,  Tauplussvfit, HadPions_plus, HadPionsCharge_plus),Wspin);

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
		  //}

		  if(a1a1MVA && svfitAlgo1.isValidSolution())
			      {
				if(Tauminussvfit.DeltaR(Tauminustruth)<0.4)
				  {
				    TauSVFitPxResPullMVA.at(t).Fill((Tauminussvfit.X()-Tauminustruth.X())/Tauminustruth.X(),w);
				    TauSVFitPyResPullMVA.at(t).Fill((Tauminussvfit.Y()-Tauminustruth.Y())/Tauminustruth.Y(),w);
				    TauSVFitPzResPullMVA.at(t).Fill((Tauminussvfit.Z()-Tauminustruth.Z())/Tauminustruth.Z(),w);

				    TauSVFitPxResPullMVA.at(t).Fill((Tauplussvfit.X()-Tauplustruth.X())/Tauplustruth.X(),w);
				    TauSVFitPyResPullMVA.at(t).Fill((Tauplussvfit.Y()-Tauplustruth.Y())/Tauplustruth.Y(),w);
				    TauSVFitPzResPullMVA.at(t).Fill((Tauplussvfit.Z()-Tauplustruth.Z())/Tauplustruth.Z(),w);
				  }
				else if(Tauminussvfit.DeltaR(Tauplustruth)<0.4)
				  {
				    TauSVFitPxResPullMVA.at(t).Fill((Tauminussvfit.X()-Tauplustruth.X())/Tauplustruth.X(),w);
				    TauSVFitPyResPullMVA.at(t).Fill((Tauminussvfit.Y()-Tauplustruth.Y())/Tauplustruth.Y(),w);
				    TauSVFitPzResPullMVA.at(t).Fill((Tauminussvfit.Z()-Tauplustruth.Z())/Tauplustruth.Z(),w);
		  
				    TauSVFitPxResPullMVA.at(t).Fill((Tauplussvfit.X()-Tauminustruth.X())/Tauminustruth.X(),w);
				    TauSVFitPyResPullMVA.at(t).Fill((Tauplussvfit.Y()-Tauminustruth.Y())/Tauminustruth.Y(),w);
				    TauSVFitPzResPullMVA.at(t).Fill((Tauplussvfit.Z()-Tauminustruth.Z())/Tauminustruth.Z(),w);
		  
				  }
				if(TauminusPairConstraint.DeltaR(Tauminustruth)<0.4)
				  {
				    TauPxResPullMVA.at(t).Fill((TauminusPairConstraint.X()-Tauminustruth.X())/Tauminustruth.X(),w);
				    TauPyResPullMVA.at(t).Fill((TauminusPairConstraint.Y()-Tauminustruth.Y())/Tauminustruth.Y(),w);
				    TauPzResPullMVA.at(t).Fill((TauminusPairConstraint.Z()-Tauminustruth.Z())/Tauminustruth.Z(),w);

	    
				    TauPxResPullMVA.at(t).Fill((TauplusPairConstraint.X()-Tauplustruth.X())/Tauplustruth.X(),w);
				    TauPyResPullMVA.at(t).Fill((TauplusPairConstraint.Y()-Tauplustruth.Y())/Tauplustruth.Y(),w);
				    TauPzResPullMVA.at(t).Fill((TauplusPairConstraint.Z()-Tauplustruth.Z())/Tauplustruth.Z(),w);

				  }
				else if (TauminusPairConstraint.DeltaR(Tauplustruth)<0.4)
				  {
				    TauPxResPullMVA.at(t).Fill((TauminusPairConstraint.X()-Tauplustruth.X())/Tauplustruth.X(),w);
				    TauPyResPullMVA.at(t).Fill((TauminusPairConstraint.Y()-Tauplustruth.Y())/Tauplustruth.Y(),w);
				    TauPzResPullMVA.at(t).Fill((TauminusPairConstraint.Z()-Tauplustruth.Z())/Tauplustruth.Z(),w);

	    
				    TauPxResPullMVA.at(t).Fill((TauplusPairConstraint.X()-Tauminustruth.X())/Tauminustruth.X(),w);
				    TauPyResPullMVA.at(t).Fill((TauplusPairConstraint.Y()-Tauminustruth.Y())/Tauminustruth.Y(),w);
				    TauPzResPullMVA.at(t).Fill((TauplusPairConstraint.Z()-Tauminustruth.Z())/Tauminustruth.Z(),w);
				  }
			      }
		  
		}

		// if(TauminusPairConstraintNoBSTIP!=zeroLV && TauplusPairConstraintNoBSTIP!=zeroLV && ScalcPVRefitNoBSTIP.isOk("a1", "a1", TauminusPairConstraintNoBSTIP, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSTIP, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitNoBSTIPMVA.at(t).Fill(ScalcPVRefitNoBSTIP.AcopAngle("a1", "a1", TauminusPairConstraintNoBSTIP, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSTIP, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		// if(TauminusPairConstraintBSTIP!=zeroLV && TauplusPairConstraintBSTIP!=zeroLV && ScalcPVRefitBSTIP.isOk("a1", "a1", TauminusPairConstraintBSTIP, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSTIP, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitBSTIPMVA.at(t).Fill(ScalcPVRefitBSTIP.AcopAngle("a1", "a1", TauminusPairConstraintBSTIP, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSTIP, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		// if(TauminusPairConstraintNoBSZNominalTIP!=zeroLV && TauplusPairConstraintNoBSZNominalTIP!=zeroLV && ScalcPVRefitNoBSZNominalTIP.isOk("a1", "a1", TauminusPairConstraintNoBSZNominalTIP, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominalTIP, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitNoBSZNominalTIPMVA.at(t).Fill(ScalcPVRefitNoBSZNominalTIP.AcopAngle("a1", "a1", TauminusPairConstraintNoBSZNominalTIP, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintNoBSZNominalTIP, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);

		// if(TauminusPairConstraintBSZNominalTIP!=zeroLV && TauplusPairConstraintBSZNominalTIP!=zeroLV && ScalcPVRefitBSZNominalTIP.isOk("a1", "a1", TauminusPairConstraintBSZNominalTIP, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominalTIP, HadRefitPions_plus, HadRefitPionsCharge_plus)==true) polarimetricAcopAnglePVRefitBSZNominalTIPMVA.at(t).Fill(ScalcPVRefitBSZNominalTIP.AcopAngle("a1", "a1", TauminusPairConstraintBSZNominalTIP, HadRefitPions_minus, HadRefitPionsCharge_minus, TauplusPairConstraintBSZNominalTIP, HadRefitPions_plus, HadRefitPionsCharge_plus),Wspin);
	      }

	    // if(a1a1TruthSVFit && Ntp->GetTruthTauLV(5,0)!=zeroLV && Ntp->GetTruthTauLV(5,1)!=zeroLV && Ntp->GetTruthTauLV(5,0)!=Ntp->GetTruthTauLV(5,1) && ScalcPtTruthA1.isOk("a1", "a1",Ntp->GetTruthTauLV(5,0) , HadPions_minus, HadPionsCharge_minus,Ntp->GetTruthTauLV(5,1), HadPions_plus, HadPionsCharge_plus)==true) polarimetricAcopAnglePtTruthA1.at(t).Fill(ScalcPtTruthA1.AcopAngle("a1", "a1",Ntp->GetTruthTauLV(5,0) , HadPions_minus, HadPionsCharge_minus,Ntp->GetTruthTauLV(5,1), HadPions_plus, HadPionsCharge_plus),Wspin);
	
	    // if(a1a1TruthSVFitMVA && Ntp->GetTruthTauLV(5,0)!=zeroLV && Ntp->GetTruthTauLV(5,1)!=zeroLV && Ntp->GetTruthTauLV(5,0)!=Ntp->GetTruthTauLV(5,1) && ScalcMVAPtTruthA1.isOk("a1", "a1", Ntp->GetTruthTauLV(5,0), HadPions_minus, HadPionsCharge_minus, Ntp->GetTruthTauLV(5,1) , HadPions_plus, HadPionsCharge_plus)==true) polarimetricAcopAngleMVAPtTruthA1.at(t).Fill(ScalcMVAPtTruthA1.AcopAngle("a1", "a1", Ntp->GetTruthTauLV(5,0), HadPions_minus, HadPionsCharge_minus, Ntp->GetTruthTauLV(5,1) , HadPions_plus, HadPionsCharge_plus),Wspin);
	    
	    //cout<<a1minusTruthSVFit<<"  "<<rhoplusTruthSVFit<<endl;
	    // if(a1minusTruthSVFit && rhoplusTruthSVFit && Ntp->GetTruthTauLV(5,0)!=zeroLV && Ntp->GetTruthTauLV(4,1)!=zeroLV && ScalcPtTruthRho1.isOk("a1", "rho", Ntp->GetTruthTauLV(5,0), HadPions_minus, HadPionsCharge_minus, Ntp->GetTruthTauLV(4,1) , HadPions_plus, HadPionsCharge_plus)==true) /*cout<<"Rentre"<<endl;*/ polarimetricAcopAnglePtTruthRho.at(t).Fill(ScalcPtTruthRho1.AcopAngle("a1", "rho", Ntp->GetTruthTauLV(5,0), HadPions_minus, HadPionsCharge_minus, Ntp->GetTruthTauLV(4,1) , HadPions_plus, HadPionsCharge_plus),Wspin);
	    
	    // if(a1plusTruthSVFit && rhominusTruthSVFit && Ntp->GetTruthTauLV(4,0)!=zeroLV && Ntp->GetTruthTauLV(5,1)!=zeroLV && ScalcPtTruthRho2.isOk("rho", "a1", Ntp->GetTruthTauLV(4,0), HadPions_minus, HadPionsCharge_minus, Ntp->GetTruthTauLV(5,1) , HadPions_plus, HadPionsCharge_plus)==true) polarimetricAcopAnglePtTruthRho.at(t).Fill(ScalcPtTruthRho2.AcopAngle("rho", "a1", Ntp->GetTruthTauLV(4,0), HadPions_minus, HadPionsCharge_minus, Ntp->GetTruthTauLV(5,1) , HadPions_plus, HadPionsCharge_plus),Wspin);
	    
	    // if(a1minusTruthSVFitMVA && rhoplusTruthSVFitMVA && Ntp->GetTruthTauLV(5,0)!=zeroLV && Ntp->GetTruthTauLV(4,1)!=zeroLV && ScalcMVAPtTruthRho1.isOk("a1", "rho", Ntp->GetTruthTauLV(5,0), HadPions_minus, HadPionsCharge_minus, Ntp->GetTruthTauLV(4,1) , HadPions_plus, HadPionsCharge_plus)==true) polarimetricAcopAngleMVAPtTruthRho.at(t).Fill(ScalcMVAPtTruthRho1.AcopAngle("a1", "rho", Ntp->GetTruthTauLV(5,0), HadPions_minus, HadPionsCharge_minus, Ntp->GetTruthTauLV(4,1) , HadPions_plus, HadPionsCharge_plus),Wspin);
	    
	    // if(a1plusTruthSVFitMVA && rhominusTruthSVFitMVA && Ntp->GetTruthTauLV(4,0)!=zeroLV && Ntp->GetTruthTauLV(5,1)!=zeroLV && ScalcMVAPtTruthRho2.isOk("rho", "a1", Ntp->GetTruthTauLV(4,0), HadPions_minus, HadPionsCharge_minus, Ntp->GetTruthTauLV(5,1) , HadPions_plus, HadPionsCharge_plus)==true) polarimetricAcopAngleMVAPtTruthRho.at(t).Fill(ScalcMVAPtTruthRho2.AcopAngle("rho", "a1", Ntp->GetTruthTauLV(4,0), HadPions_minus, HadPionsCharge_minus, Ntp->GetTruthTauLV(5,1) , HadPions_plus, HadPionsCharge_plus),Wspin);


	    // if(a1minusTruthSVFit && piplusTruthSVFit && Ntp->GetTruthTauLV(5,0)!=zeroLV && Ntp->GetTruthTauLV(3,1)!=zeroLV && ScalcPtTruthPi1.isOk("a1", "pion", Ntp->GetTruthTauLV(5,0), HadPions_minus, HadPionsCharge_minus, Ntp->GetTruthTauLV(3,1) , HadPions_plus, HadPionsCharge_plus)==true) polarimetricAcopAnglePtTruthPi.at(t).Fill(ScalcPtTruthPi1.AcopAngle("a1", "pion", Ntp->GetTruthTauLV(5,0), HadPions_minus, HadPionsCharge_minus, Ntp->GetTruthTauLV(3,1) , HadPions_plus, HadPionsCharge_plus),Wspin);
	    // if(a1plusTruthSVFit && piminusTruthSVFit && Ntp->GetTruthTauLV(3,0)!=zeroLV && Ntp->GetTruthTauLV(5,1)!=zeroLV && ScalcPtTruthPi2.isOk("pion", "a1", Ntp->GetTruthTauLV(3,0), HadPions_minus, HadPionsCharge_minus, Ntp->GetTruthTauLV(5,1) , HadPions_plus, HadPionsCharge_plus)==true) polarimetricAcopAnglePtTruthPi.at(t).Fill(ScalcPtTruthPi2.AcopAngle("pion", "a1", Ntp->GetTruthTauLV(3,0), HadPions_minus, HadPionsCharge_minus, Ntp->GetTruthTauLV(5,1) , HadPions_plus, HadPionsCharge_plus),Wspin);
	    // if(a1minusTruthSVFitMVA && piplusTruthSVFitMVA && Ntp->GetTruthTauLV(5,0)!=zeroLV && Ntp->GetTruthTauLV(3,1)!=zeroLV && ScalcMVAPtTruthPi1.isOk("a1", "pion", Ntp->GetTruthTauLV(5,0), HadPions_minus, HadPionsCharge_minus, Ntp->GetTruthTauLV(3,1) , HadPions_plus, HadPionsCharge_plus)==true) polarimetricAcopAngleMVAPtTruthPi.at(t).Fill(ScalcMVAPtTruthPi1.AcopAngle("a1", "pion", Ntp->GetTruthTauLV(5,0), HadPions_minus, HadPionsCharge_minus, Ntp->GetTruthTauLV(3,1) , HadPions_plus, HadPionsCharge_plus),Wspin);
	    // if(a1plusTruthSVFitMVA && piminusTruthSVFitMVA && Ntp->GetTruthTauLV(3,0)!=zeroLV && Ntp->GetTruthTauLV(5,1)!=zeroLV && ScalcMVAPtTruthPi2.isOk("pion", "a1", Ntp->GetTruthTauLV(3,0), HadPions_minus, HadPionsCharge_minus, Ntp->GetTruthTauLV(5,1) , HadPions_plus, HadPionsCharge_plus)==true) polarimetricAcopAngleMVAPtTruthPi.at(t).Fill(ScalcMVAPtTruthPi2.AcopAngle("pion", "a1", Ntp->GetTruthTauLV(3,0), HadPions_minus, HadPionsCharge_minus, Ntp->GetTruthTauLV(5,1) , HadPions_plus, HadPionsCharge_plus),Wspin);
	   
	    // 	//	cout<<"-------------"<<endl;
	     
	     

	    // if(a1minusTruthSVFit && svfitAlgo1.isValidSolution() && rhoplusTruthSVFit && ScalcSVFitRho1.isOk("a1", "rho", Tauminussvfit, HadPions_minus, HadPionsCharge_minus, Tauplussvfit , HadPions_plus, HadPionsCharge_plus)==true) polarimetricAcopAngleSVFitRho.at(t).Fill(ScalcSVFitRho1.AcopAngle("a1", "rho", Tauminussvfit, HadPions_minus, HadPionsCharge_minus, Tauplussvfit , HadPions_plus, HadPionsCharge_plus),Wspin);
	    // if(a1plusTruthSVFit && svfitAlgo1.isValidSolution() && rhominusTruthSVFit && ScalcSVFitRho2.isOk("rho", "a1", Tauminussvfit, HadPions_minus, HadPionsCharge_minus, Tauplussvfit , HadPions_plus, HadPionsCharge_plus)==true) polarimetricAcopAngleSVFitRho.at(t).Fill(ScalcSVFitRho2.AcopAngle("rho", "a1", Tauminussvfit, HadPions_minus, HadPionsCharge_minus, Tauplussvfit , HadPions_plus, HadPionsCharge_plus),Wspin);
	    // if(a1minusTruthSVFitMVA && svfitAlgo1.isValidSolution() && rhoplusTruthSVFitMVA && ScalcMVASVFitRho1.isOk("a1", "rho", Tauminussvfit, HadPions_minus, HadPionsCharge_minus, Tauplussvfit, HadPions_plus, HadPionsCharge_plus)==true) polarimetricAcopAngleMVASVFitRho.at(t).Fill(ScalcMVASVFitRho1.AcopAngle("a1", "rho", Tauminussvfit, HadPions_minus, HadPionsCharge_minus, Tauplussvfit, HadPions_plus, HadPionsCharge_plus),Wspin);
	    // if(a1plusTruthSVFitMVA && svfitAlgo1.isValidSolution() && rhominusTruthSVFitMVA && ScalcMVASVFitRho2.isOk("rho", "a1", Tauminussvfit, HadPions_minus, HadPionsCharge_minus, Tauplussvfit, HadPions_plus, HadPionsCharge_plus)==true) polarimetricAcopAngleMVASVFitRho.at(t).Fill(ScalcMVASVFitRho2.AcopAngle("rho", "a1", Tauminussvfit, HadPions_minus, HadPionsCharge_minus, Tauplussvfit, HadPions_plus, HadPionsCharge_plus),Wspin);

	    // if(a1minusTruthSVFit && svfitAlgo1.isValidSolution() && piplusTruthSVFit && ScalcSVFitPi1.isOk("a1", "pion", Tauminussvfit, HadPions_minus, HadPionsCharge_minus,  Tauplussvfit, HadPions_plus, HadPionsCharge_plus)==true) polarimetricAcopAngleSVFitPi.at(t).Fill(ScalcSVFitPi1.AcopAngle("a1", "pion", Tauminussvfit, HadPions_minus, HadPionsCharge_minus,  Tauplussvfit, HadPions_plus, HadPionsCharge_plus),Wspin);
	    // if(a1plusTruthSVFit && svfitAlgo1.isValidSolution() && piminusTruthSVFit && ScalcSVFitPi2.isOk("pion", "a1", Tauminussvfit, HadPions_minus, HadPionsCharge_minus,  Tauplussvfit, HadPions_plus, HadPionsCharge_plus)==true) polarimetricAcopAngleSVFitPi.at(t).Fill(ScalcSVFitPi2.AcopAngle("pion", "a1", Tauminussvfit, HadPions_minus, HadPionsCharge_minus,  Tauplussvfit, HadPions_plus, HadPionsCharge_plus),Wspin);
	    // if(a1minusTruthSVFitMVA && svfitAlgo1.isValidSolution() && piplusTruthSVFitMVA && ScalcMVASVFitPi1.isOk("a1", "pion", Tauminussvfit, HadPions_minus, HadPionsCharge_minus,  Tauplussvfit, HadPions_plus, HadPionsCharge_plus)==true) polarimetricAcopAngleMVASVFitPi.at(t).Fill(ScalcMVASVFitPi1.AcopAngle("a1", "pion", Tauminussvfit, HadPions_minus, HadPionsCharge_minus,  Tauplussvfit, HadPions_plus, HadPionsCharge_plus),Wspin);
	    // if(a1plusTruthSVFitMVA && svfitAlgo1.isValidSolution() && piminusTruthSVFitMVA && ScalcMVASVFitPi2.isOk("pion", "a1", Tauminussvfit, HadPions_minus, HadPionsCharge_minus,  Tauplussvfit, HadPions_plus, HadPionsCharge_plus)==true) polarimetricAcopAngleMVASVFitPi.at(t).Fill(ScalcMVASVFitPi2.AcopAngle("pion", "a1", Tauminussvfit, HadPions_minus, HadPionsCharge_minus,  Tauplussvfit, HadPions_plus, HadPionsCharge_plus),Wspin);
	     
	    // if(a1a1 && !a1a1_)
	    //   {
	    // 	if(TauminusPairConstraint!=TauplusPairConstraint && TauminusPairConstraint!=zeroLV && TauplusPairConstraint!=zeroLV && ScalcBackground.isOk("a1", "a1", TauminusPairConstraint, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraint , HadPions_plus, HadPionsCharge_plus)==true) polarimetricAcopAngleBackground.at(t).Fill(ScalcBackground.AcopAngle("a1", "a1", TauminusPairConstraint, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraint , HadPions_plus, HadPionsCharge_plus),Wspin);
	    // 	if(TauminusPairConstraintNoBS!=TauplusPairConstraintNoBS && TauminusPairConstraintNoBS!=zeroLV &&  TauplusPairConstraintNoBS!=zeroLV && hasNoBS && ScalcPVRefitNoBSBackground.isOk("a1", "a1", TauminusPairConstraintNoBS, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraintNoBS , HadPions_plus, HadPionsCharge_plus)==true)polarimetricAcopAnglePVRefitNoBSBackground.at(t).Fill(ScalcPVRefitNoBSBackground.AcopAngle("a1", "a1", TauminusPairConstraintNoBS, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraintNoBS , HadPions_plus, HadPionsCharge_plus),Wspin);
	    // 	if(TauminusPairConstraintBS!=TauplusPairConstraintBS && TauminusPairConstraintBS!=zeroLV &&  TauplusPairConstraintBS!=zeroLV && hasNoBS && ScalcPVRefitBSBackground.isOk("a1", "a1", TauminusPairConstraintBS, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraintBS , HadPions_plus, HadPionsCharge_plus)==true)polarimetricAcopAnglePVRefitBSBackground.at(t).Fill(ScalcPVRefitBSBackground.AcopAngle("a1", "a1", TauminusPairConstraintBS, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraintBS , HadPions_plus, HadPionsCharge_plus),Wspin);
	    //   }
	      
	    // if(a1a1MVA && !a1a1_)
	    //   {
	    // 	if(TauminusPairConstraintMVA!=TauplusPairConstraintMVA && TauminusPairConstraintMVA!=zeroLV && TauplusPairConstraintMVA!=zeroLV && ScalcMVABackground.isOk("a1", "a1", TauminusPairConstraintMVA, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraintMVA , HadPions_plus, HadPionsCharge_plus)==true) polarimetricAcopAngleMVABackground.at(t).Fill(ScalcMVABackground.AcopAngle("a1", "a1", TauminusPairConstraintMVA, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraintMVA , HadPions_plus, HadPionsCharge_plus),Wspin);
	    // 	if(TauminusPairConstraintNoBSMVA!=TauplusPairConstraintNoBSMVA && TauminusPairConstraintNoBSMVA!=zeroLV && TauplusPairConstraintNoBSMVA!=zeroLV && hasNoBSMVA && ScalcPVRefitNoBSMVABackground.isOk("a1", "a1", TauminusPairConstraintNoBSMVA, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraintNoBSMVA, HadPions_plus, HadPionsCharge_plus)==true)polarimetricAcopAnglePVRefitNoBSMVABackground.at(t).Fill(ScalcPVRefitNoBSMVABackground.AcopAngle("a1", "a1", TauminusPairConstraintNoBSMVA, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraintNoBSMVA, HadPions_plus, HadPionsCharge_plus),Wspin);
	    // 	if(TauminusPairConstraintBSMVA!=TauplusPairConstraintBSMVA && TauminusPairConstraintBSMVA!=zeroLV && TauplusPairConstraintBSMVA!=zeroLV && hasBSMVA && ScalcAllPVRefitBSMVABackground.isOk("a1", "a1", TauminusPairConstraintBSMVA, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraintBSMVA , HadPions_plus, HadPionsCharge_plus)==true)polarimetricAcopAnglePVRefitBSMVABackground.at(t).Fill(ScalcAllPVRefitBSMVABackground.AcopAngle("a1", "a1", TauminusPairConstraintBSMVA, HadPions_minus, HadPionsCharge_minus, TauplusPairConstraintBSMVA , HadPions_plus, HadPionsCharge_plus),Wspin);
	    // }
	     
	    }

	if(a1a1_ /*|| a1rho || a1pi*/)
	      {
	   
		if(a1a1_)
		  {
		    Pions1=Ntp->GetTruthPionsFromA1(0);
		    Pions1Charge.push_back(1);
		    Pions1Charge.push_back(-1);
		    Pions1Charge.push_back(-1);
		    Pions2=Ntp->GetTruthPionsFromA1(1);
		    Pions2Charge.push_back(-1);
		    Pions2Charge.push_back(1);
		    Pions2Charge.push_back(1);
		  }
		// if(a1rho)
		//   {
		//     if(rho0 && a11)
		//       {
		// 	Pions1=Ntp->GetTruthPionsFromRho(0);
		// 	Pions1Charge.push_back(-1);
		// 	Pions1Charge.push_back(0);
		// 	Pions2=Ntp->GetTruthPionsFromA1(1);
		// 	Pions2Charge.push_back(-1);
		// 	Pions2Charge.push_back(+1);
		// 	Pions2Charge.push_back(+1);
		//       }
		//     else if (rho1 && a10)
		//       {
		// 	Pions1=Ntp->GetTruthPionsFromA1(0);
		// 	Pions1Charge.push_back(1);
		// 	Pions1Charge.push_back(-1);
		// 	Pions1Charge.push_back(-1);
		// 	Pions2=Ntp->GetTruthPionsFromRho(1);
		// 	Pions2Charge.push_back(1);
		// 	Pions2Charge.push_back(0);
		//       }

		//   }

		// if(a1pi)
		//   {

		//     if(pi1 && a10)
		//       {
		// 	Pions1=Ntp->GetTruthPionsFromA1(0);
		// 	Pions1Charge.push_back(1);
		// 	Pions1Charge.push_back(-1);
		// 	Pions1Charge.push_back(-1);
		// 	Pions2.push_back(Ntp->GetTruthTauProductLV(3,211,1));
		//       }
		//     else if (pi0 && a11)
		//       {
		// 	Pions2=Ntp->GetTruthPionsFromA1(1);
		// 	Pions2Charge.push_back(-1);
		// 	Pions2Charge.push_back(1);
		// 	Pions2Charge.push_back(1);
		// 	Pions1.push_back(Ntp->GetTruthTauProductLV(3,211,0));
		//       }

		//     // cout<<endl;
		//     // for(int j=0;j<Ntp->NMCTauDecayProducts(Tauminus);j++)cout<<Ntp->MCTauandProd_pdgid(Tauminus, j)<<"   ";
		//     // cout<<endl;
		//     // cout<<endl;
		//     // for(int j=0;j<Ntp->NMCTauDecayProducts(Tauplus);j++)cout<<Ntp->MCTauandProd_pdgid(Tauplus, j)<<"   ";
		//     // cout<<endl;
		//   }

		    
		// HadPionsTruth_minus.push_back(Pions1.at(0));
		// if(Ntp->CheckDecayID(4,5) || Ntp->CheckDecayID(5,5)) HadPionsTruth_minus.push_back(Pions1.at(1)); 
		// if(Ntp->CheckDecayID(5,5)) HadPionsTruth_minus.push_back(Pions1.at(2));
		// HadPionsTruth_plus.push_back(Pions2.at(0));
		// if(Ntp->CheckDecayID(4,5) || Ntp->CheckDecayID(5,5)) HadPionsTruth_plus.push_back(Pions2.at(1)); 
		// if(Ntp->CheckDecayID(5,5)) HadPionsTruth_plus.push_back(Pions2.at(2));
		    
		// HadPionsChargeTruth_minus.push_back(Ntp->MCTauandProd_charge(0, 2));
		// if(Ntp->CheckDecayID(4,5) || Ntp->CheckDecayID(5,5)) HadPionsChargeTruth_minus.push_back(Ntp->MCTauandProd_charge(0, 3));
		// if(Ntp->CheckDecayID(5,5)) HadPionsChargeTruth_minus.push_back(Ntp->MCTauandProd_charge(0, 4));
		// HadPionsChargeTruth_plus.push_back(Ntp->MCTauandProd_charge(1, 2));
		// if(Ntp->CheckDecayID(4,5) || Ntp->CheckDecayID(5,5)) HadPionsChargeTruth_plus.push_back(Ntp->MCTauandProd_charge(1, 3));
		// if(Ntp->CheckDecayID(5,5)) HadPionsChargeTruth_plus.push_back(Ntp->MCTauandProd_charge(1, 4));
		   
		if(a1a1_)
		  {
		    SCalculator Scalc1Truth("a1");
		    SCalculator Scalc2Truth("a1");

		    Scalc1Truth.SortPions(Pions1, Pions1Charge);
		    Scalc2Truth.SortPions(Pions2, Pions2Charge);
		    
		    if(a1a1)
		      {
			//TauminusPairConstraint.Print();
			//Ntp->GetTruthTauLV(5,0).Print();
			//TauplusPairConstraint.Print();
			//Ntp->GetTruthTauLV(5,1).Print();
			//for(int i=0;i<3;i++){Pions1.at(i).Print();}
			//for(int i=0;i<3;i++){HadPions_minus.at(i).Print();}
		      }
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
		// if(a1rho)
		//   {
		//     if(rho0 && a11)
		//       {
		// 	SCalculator Scalc1Truth("rho");
		// 	SCalculator Scalc2Truth("a1");

		// 	Scalc1Truth.SortPions(Pions1, Pions1Charge);
		// 	Scalc2Truth.SortPions(Pions2, Pions2Charge);
		    
		// 	if(Ntp->MCTau_pdgid(0)==15)
		// 	  {
		// 	    tauandprodTruthminus.push_back(Ntp->GetTruthTauLV(4,0));
		// 	    tauandprodTruthminus.push_back(Pions1.at(0));
		// 	    tauandprodTruthminus.push_back(Pions1.at(1));
		// 	    tauandprodTruthplus.push_back(Ntp->GetTruthTauLV(5,1));   
		// 	    tauandprodTruthplus.push_back(Pions2.at(0));   
		// 	    tauandprodTruthplus.push_back(Pions2.at(1));   
		// 	    tauandprodTruthplus.push_back(Pions2.at(2));   
		// 	  }
		// 	else if (Ntp->MCTau_pdgid(0)==-15)
		// 	  {
		// 	    tauandprodTruthminus.push_back(Ntp->GetTruthTauLV(5,1));
		// 	    tauandprodTruthminus.push_back(Pions2.at(0));
		// 	    tauandprodTruthminus.push_back(Pions2.at(1));
		// 	    tauandprodTruthminus.push_back(Pions2.at(2));
		// 	    tauandprodTruthplus.push_back(Ntp->GetTruthTauLV(4,0));   
		// 	    tauandprodTruthplus.push_back(Pions1.at(0));   
		// 	    tauandprodTruthplus.push_back(Pions1.at(1));    
		// 	  }

		//       }
		//     else if(rho1 && a10)
		//       {
		// 	SCalculator Scalc1Truth("a1");
		// 	SCalculator Scalc2Truth("rho");

		// 	Scalc1Truth.SortPions(Pions1, Pions1Charge);
		// 	Scalc2Truth.SortPions(Pions2, Pions2Charge);
		// 	if(Ntp->MCTau_pdgid(0)==15)
		// 	  {
		// 	    tauandprodTruthminus.push_back(Ntp->GetTruthTauLV(5,0));
		// 	    tauandprodTruthminus.push_back(Pions1.at(0));
		// 	    tauandprodTruthminus.push_back(Pions1.at(1));
		// 	    tauandprodTruthminus.push_back(Pions1.at(2));
		// 	    tauandprodTruthplus.push_back(Ntp->GetTruthTauLV(4,1));   
		// 	    tauandprodTruthplus.push_back(Pions2.at(0));   
		// 	    tauandprodTruthplus.push_back(Pions2.at(1));   
		// 	  }
		// 	else if (Ntp->MCTau_pdgid(0)==-15)
		// 	  {
		// 	    tauandprodTruthminus.push_back(Ntp->GetTruthTauLV(4,1));
		// 	    tauandprodTruthminus.push_back(Pions2.at(0));
		// 	    tauandprodTruthminus.push_back(Pions2.at(1));
		// 	    tauandprodTruthplus.push_back(Ntp->GetTruthTauLV(5,0));   
		// 	    tauandprodTruthplus.push_back(Pions1.at(0));   
		// 	    tauandprodTruthplus.push_back(Pions1.at(1));    
		// 	    tauandprodTruthplus.push_back(Pions1.at(2));   
		// 	  }
		//       }

		   
		//   }

		// if(a1pi)
		//   {
		//     if(pi0 && a11)
		//       {
		// 	SCalculator Scalc1Truth("pion");
		// 	SCalculator Scalc2Truth("a1");

		// 	Scalc2Truth.SortPions(Pions2, Pions2Charge);
		// 	if(Ntp->MCTau_pdgid(0)==15)
		// 	  {
		// 	    tauandprodTruthminus.push_back(Ntp->GetTruthTauLV(3,0));
		// 	    tauandprodTruthminus.push_back(Pions1.at(0));
		// 	    tauandprodTruthplus.push_back(Ntp->GetTruthTauLV(5,1));   
		// 	    tauandprodTruthplus.push_back(Pions2.at(0));   
		// 	    tauandprodTruthplus.push_back(Pions2.at(1));   
		// 	    tauandprodTruthplus.push_back(Pions2.at(2));   
		// 	  }
		// 	else if (Ntp->MCTau_pdgid(0)==-15)
		// 	  {
		// 	    tauandprodTruthminus.push_back(Ntp->GetTruthTauLV(5,1));
		// 	    tauandprodTruthminus.push_back(Pions2.at(0));
		// 	    tauandprodTruthminus.push_back(Pions2.at(1));
		// 	    tauandprodTruthminus.push_back(Pions2.at(2));
		// 	    tauandprodTruthplus.push_back(Ntp->GetTruthTauLV(3,0));   
		// 	    tauandprodTruthplus.push_back(Pions1.at(0));     
		// 	  }
		    
		//       }
		//     else if (pi1 && a10)
		//       {
		// 	SCalculator Scalc1Truth("a1");
		// 	SCalculator Scalc2Truth("pion");

		// 	Scalc1Truth.SortPions(Pions1, Pions1Charge);
		// 	if(Ntp->MCTau_pdgid(0)==15)
		// 	  {
		// 	    tauandprodTruthminus.push_back(Ntp->GetTruthTauLV(5,0));
		// 	    tauandprodTruthminus.push_back(Pions1.at(0));
		// 	    tauandprodTruthminus.push_back(Pions1.at(1));
		// 	    tauandprodTruthminus.push_back(Pions1.at(2));
		// 	    tauandprodTruthplus.push_back(Ntp->GetTruthTauLV(3,1));   
		// 	    tauandprodTruthplus.push_back(Pions2.at(0));   
		// 	  }
		// 	else if (Ntp->MCTau_pdgid(0)==-15)
		// 	  {
		// 	    tauandprodTruthminus.push_back(Ntp->GetTruthTauLV(3,1));
		// 	    tauandprodTruthminus.push_back(Pions2.at(0));
		// 	    tauandprodTruthplus.push_back(Ntp->GetTruthTauLV(5,0));   
		// 	    tauandprodTruthplus.push_back(Pions1.at(0));   
		// 	    tauandprodTruthplus.push_back(Pions1.at(1));    
		// 	    tauandprodTruthplus.push_back(Pions1.at(2));   
		// 	  }

		//       }
		    
		//  }
		if((Pions1!=Pions2) && (Pions1!=VectZeroLV) && (Pions2!=VectZeroLV) && tauandprodTruthminus.at(0)!=zeroLV &&tauandprodTruthplus.at(0)!=zeroLV &&tauandprodTruthminus.at(0)!=tauandprodTruthplus.at(0))
		{
			if(a1a1_)
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
			// if(rho0 && a11)
			//   {
		   
			//     SCalculator Scalc1Truth("rho");
			//     SCalculator Scalc2Truth("a1");

			//     Scalc1Truth.Configure(tauandprodTruthminus,tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0), -1);
			//     TVector3 h1Truth=Scalc1Truth.pv();
			
			//     Scalc2Truth.Configure(tauandprodTruthplus,tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0), +1);
			//     TVector3 h2Truth=Scalc2Truth.pv();
			
			//     double h1TruthNorm=1./h1Truth.Mag();
			//     double h2TruthNorm=1./h2Truth.Mag();
			
			//     if(std::isnan(h1TruthNorm)!=true && std::isnan(h2TruthNorm)!=true)
			//       { 
			// 	TLorentzVector tauminusTruth_HRF = Scalc1Truth.Boost(tauandprodTruthminus.at(0),tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0));
			// 	TLorentzVector tauplusTruth_HRF  = Scalc2Truth.Boost(tauandprodTruthplus.at(0),tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0));
			
			// 	double norm1Truth=1./(((h1Truth*h1TruthNorm).Cross(tauminusTruth_HRF.Vect().Unit())).Mag());
			// 	double norm2Truth=1./(((h2Truth*h2TruthNorm).Cross(tauplusTruth_HRF.Vect().Unit())).Mag());
			// 	TVector3 k1Truth = ((h1Truth*h1TruthNorm).Cross(tauminusTruth_HRF.Vect().Unit()))*norm1Truth;
			// 	TVector3 k2Truth = ((h2Truth*h2TruthNorm).Cross(tauplusTruth_HRF.Vect().Unit()))*norm2Truth;
			
			// 	if(((h1Truth*h1TruthNorm).Cross(h2Truth*h2TruthNorm))*tauminusTruth_HRF.Vect().Unit()<=0)polarimetricAcopAngleTruthRho.at(t).Fill(TMath::ATan2((k1Truth.Cross(k2Truth)).Mag(),k1Truth*k2Truth),Wspin);
			// 	else polarimetricAcopAngleTruthRho.at(t).Fill(2*TMath::Pi()-TMath::ATan2((k1Truth.Cross(k2Truth)).Mag(),k1Truth*k2Truth),Wspin);
			//       }
			//   }
			// if(rho1 && a10)
			//   {
			//     SCalculator Scalc1Truth("a1");
			//     SCalculator Scalc2Truth("rho");

			//     Scalc1Truth.Configure(tauandprodTruthminus,tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0), -1);
			//     TVector3 h1Truth=Scalc1Truth.pv();
			
			//     Scalc2Truth.Configure(tauandprodTruthplus,tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0), +1);
			//     TVector3 h2Truth=Scalc2Truth.pv();
			
			//     double h1TruthNorm=1./h1Truth.Mag();
			//     double h2TruthNorm=1./h2Truth.Mag();
		
			//     if(std::isnan(h1TruthNorm)!=true && std::isnan(h2TruthNorm)!=true)
			//       { 
			// 	TLorentzVector tauminusTruth_HRF = Scalc1Truth.Boost(tauandprodTruthminus.at(0),tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0));
			// 	TLorentzVector tauplusTruth_HRF  = Scalc2Truth.Boost(tauandprodTruthplus.at(0),tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0));
			
			// 	double norm1Truth=1./(((h1Truth*h1TruthNorm).Cross(tauminusTruth_HRF.Vect().Unit())).Mag());
			// 	double norm2Truth=1./(((h2Truth*h2TruthNorm).Cross(tauplusTruth_HRF.Vect().Unit())).Mag());
			// 	TVector3 k1Truth = ((h1Truth*h1TruthNorm).Cross(tauminusTruth_HRF.Vect().Unit()))*norm1Truth;
			// 	TVector3 k2Truth = ((h2Truth*h2TruthNorm).Cross(tauplusTruth_HRF.Vect().Unit()))*norm2Truth;
			
			// 	if(((h1Truth*h1TruthNorm).Cross(h2Truth*h2TruthNorm))*tauminusTruth_HRF.Vect().Unit()<=0)polarimetricAcopAngleTruthRho.at(t).Fill(TMath::ATan2((k1Truth.Cross(k2Truth)).Mag(),k1Truth*k2Truth),Wspin);
			// 	else polarimetricAcopAngleTruthRho.at(t).Fill(2*TMath::Pi()-TMath::ATan2((k1Truth.Cross(k2Truth)).Mag(),k1Truth*k2Truth),Wspin);
			//       }
			//   }
			// if(pi0 && a11)
			//   {
			//     //cout<<"0000000000000000"<<endl;
			//     // SCalculator testPi1;
			//     // SCalculator testPi2;
		    
			//     // test.at(t).Fill(testPi1.AcopAngle("pion", "a1", tauandprodTruthminus.at(0), Pions1, HadPionsCharge_minus, tauandprodTruthplus.at(0) , Pions2, HadPionsCharge_plus),Wspin);
		    
			//     SCalculator Scalc1Truth("pion");
			//     SCalculator Scalc2Truth("a1");

		    
			//     //cout<<"1111111111111111111"<<endl;
			//     Scalc1Truth.Configure(tauandprodTruthminus,tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0), -1);
			//     TVector3 h1Truth=Scalc1Truth.pv();
		    
			//     Scalc2Truth.Configure(tauandprodTruthplus,tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0), +1);
			//     TVector3 h2Truth=Scalc2Truth.pv();

		    
		    
			//     double h1TruthNorm=1./h1Truth.Mag();
			//     double h2TruthNorm=1./h2Truth.Mag();
			//     if(std::isnan(h1TruthNorm)!=true && std::isnan(h2TruthNorm)!=true)
			//       { 
			// 	TLorentzVector tauminusTruth_HRF = Scalc1Truth.Boost(tauandprodTruthminus.at(0),tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0));
			// 	TLorentzVector tauplusTruth_HRF  = Scalc2Truth.Boost(tauandprodTruthplus.at(0),tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0));
			
			// 	double norm1Truth=1./(((h1Truth*h1TruthNorm).Cross(tauminusTruth_HRF.Vect().Unit())).Mag());
			// 	double norm2Truth=1./(((h2Truth*h2TruthNorm).Cross(tauplusTruth_HRF.Vect().Unit())).Mag());
			// 	TVector3 k1Truth = ((h1Truth*h1TruthNorm).Cross(tauminusTruth_HRF.Vect().Unit()))*norm1Truth;
			// 	TVector3 k2Truth = ((h2Truth*h2TruthNorm).Cross(tauplusTruth_HRF.Vect().Unit()))*norm2Truth;
		    
			// 	if(((h1Truth*h1TruthNorm).Cross(h2Truth*h2TruthNorm))*tauminusTruth_HRF.Vect().Unit()<=0) {
			// 	  polarimetricAcopAngleTruthPi.at(t).Fill(TMath::ATan2((k1Truth.Cross(k2Truth)).Mag(),k1Truth*k2Truth),Wspin);
			// 	  // if(TMath::ATan2((k1Truth.Cross(k2Truth)).Mag(),k1Truth*k2Truth)>2.94)
			// 	  // 	{
			// 	  // 	  cout<<endl;
			// 	  // 	  tauminusTruth_HRF.Print();
			// 	  // 	  tauplusTruth_HRF.Print();
			// 	  // 	  k1Truth.Print();
			// 	  // 	  k2Truth.Print();
			// 	  // 	  cout<<endl;
			// 	  // 	}
			// 	}
			// 	else{ polarimetricAcopAngleTruthPi.at(t).Fill(2*TMath::Pi()-TMath::ATan2((k1Truth.Cross(k2Truth)).Mag(),k1Truth*k2Truth),Wspin);
			// 	  // if(2*TMath::Pi()-TMath::ATan2((k1Truth.Cross(k2Truth)).Mag(),k1Truth*k2Truth)<3.34){
			// 	  // 	cout<<endl;
			// 	  // 	tauminusTruth_HRF.Print();
			// 	  // 	tauplusTruth_HRF.Print();
			// 	  // 	k1Truth.Print();
			// 	  // 	k2Truth.Print();
			// 	  // 	cout<<endl;
			// 	  // }
			// 	}
			//       }
			//   }
			// if(pi1 && a10)
			//   {
			//     //cout<<"22222222222222222"<<endl;
			//     SCalculator Scalc1Truth("a1");
			//     SCalculator Scalc2Truth("pion");
		    
		    
			//     //cout<<"3333333333333333"<<endl;
			//     Scalc1Truth.Configure(tauandprodTruthminus,tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0), -1);
			//     TVector3 h1Truth=Scalc1Truth.pv();
			
			//     Scalc2Truth.Configure(tauandprodTruthplus,tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0), +1);
			//     TVector3 h2Truth=Scalc2Truth.pv();
			
			//     double h1TruthNorm=1./h1Truth.Mag();
			//     double h2TruthNorm=1./h2Truth.Mag();
			//     if(std::isnan(h1TruthNorm)!=true && std::isnan(h2TruthNorm)!=true)
			//       { 
			// 	TLorentzVector tauminusTruth_HRF = Scalc1Truth.Boost(tauandprodTruthminus.at(0),tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0));
			// 	TLorentzVector tauplusTruth_HRF  = Scalc2Truth.Boost(tauandprodTruthplus.at(0),tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0));
			
			// 	double norm1Truth=1./(((h1Truth*h1TruthNorm).Cross(tauminusTruth_HRF.Vect().Unit())).Mag());
			// 	double norm2Truth=1./(((h2Truth*h2TruthNorm).Cross(tauplusTruth_HRF.Vect().Unit())).Mag());
			// 	TVector3 k1Truth = ((h1Truth*h1TruthNorm).Cross(tauminusTruth_HRF.Vect().Unit()))*norm1Truth;
			// 	TVector3 k2Truth = ((h2Truth*h2TruthNorm).Cross(tauplusTruth_HRF.Vect().Unit()))*norm2Truth;
			
			// 	if(((h1Truth*h1TruthNorm).Cross(h2Truth*h2TruthNorm))*tauminusTruth_HRF.Vect().Unit()<=0){polarimetricAcopAngleTruthPi.at(t).Fill(TMath::ATan2((k1Truth.Cross(k2Truth)).Mag(),k1Truth*k2Truth),Wspin);
			// 	  // if(TMath::ATan2((k1Truth.Cross(k2Truth)).Mag(),k1Truth*k2Truth)>2.94)
			// 	  // 	{
			// 	  // 	  cout<<endl;
			// 	  // 	  tauminusTruth_HRF.Print();
			// 	  // 	  tauplusTruth_HRF.Print();
			// 	  // 	  k1Truth.Print();
			// 	  // 	  k2Truth.Print();
			// 	  // 	  cout<<endl;
			// 	  // 	}
			// 	}
			// 	else{ polarimetricAcopAngleTruthPi.at(t).Fill(2*TMath::Pi()-TMath::ATan2((k1Truth.Cross(k2Truth)).Mag(),k1Truth*k2Truth),Wspin);
			// 	  // if(2*TMath::Pi()-TMath::ATan2((k1Truth.Cross(k2Truth)).Mag(),k1Truth*k2Truth)<3.34){
			// 	  // 	cout<<endl;
			// 	  // 	tauminusTruth_HRF.Print();
			// 	  // 	tauplusTruth_HRF.Print();
			// 	  // 	k1Truth.Print();
			// 	  // 	k2Truth.Print();
			// 	  // 	cout<<endl;
			// 	  // }
			// 	}
			//       }
			//  }
		}
	      }
      }
    
   
	
    
	//cout<<"minus: "<<A1PionsCharge_minus.at(0)<<" "<<A1PionsCharge_minus.at(1)<<" "<<A1PionsCharge_minus.at(2)<<endl;
	//cout<<"plus: "<<A1PionsCharge_plus.at(0)<<" "<<A1PionsCharge_plus.at(1)<<" "<<A1PionsCharge_plus.at(2)<<endl;
	// vector<TLorentzVector> tauandprodSVFitminus;
	// vector<TLorentzVector> tauandprodSVFitplus;
	// TLorentzVector TauminussvfitLV;
	// TLorentzVector TauplussvfitLV;
	// TauminussvfitLV.SetPxPyPzE(Tauminussvfit.Px(),Tauminussvfit.Py(),Tauminussvfit.Pz(),Tauminussvfit.E());
	// TauplussvfitLV.SetPxPyPzE(Tauplussvfit.Px(),Tauplussvfit.Py(),Tauplussvfit.Pz(),Tauplussvfit.E());
	    
	// tauandprodSVFitminus.push_back(TauminussvfitLV);
	// tauandprodSVFitminus.push_back(A1Pions_minus.at(0));
	// tauandprodSVFitminus.push_back(A1Pions_minus.at(1));
	// tauandprodSVFitminus.push_back(A1Pions_minus.at(2));
	// tauandprodSVFitplus.push_back(TauplussvfitLV);   
	// tauandprodSVFitplus.push_back(A1Pions_plus.at(0));   
	// tauandprodSVFitplus.push_back(A1Pions_plus.at(1));   
	// tauandprodSVFitplus.push_back(A1Pions_plus.at(2));   
	// //
	// ScalcSVFit1.Configure(tauandprodSVFitminus,tauandprodSVFitminus.at(0)+tauandprodSVFitplus.at(0), -1);
	// TVector3 hSVFit1=-ScalcSVFit1.pv();
	// ScalcSVFit2.Configure(tauandprodSVFitplus,tauandprodSVFitminus.at(0)+tauandprodSVFitplus.at(0), +1);
	// TVector3 hSVFit2=-ScalcSVFit2.pv();
            
	// TLorentzVector tauSVFitminus_HRF = ScalcSVFit1.Boost(tauandprodSVFitminus.at(0),tauandprodSVFitminus.at(0)+tauandprodSVFitplus.at(0));
	// TLorentzVector tauSVFitplus_HRF  = ScalcSVFit2.Boost(tauandprodSVFitplus.at(0),tauandprodSVFitminus.at(0)+tauandprodSVFitplus.at(0));
	// //(tauminus_HRF+tauplus_HRF).Print();
	// double hSVFit1Norm=1./hSVFit1.Mag();
	// double hSVFit2Norm=1./hSVFit2.Mag();
	// hSVFit1=hSVFit1*hSVFit1Norm;
	// hSVFit2=hSVFit2*hSVFit2Norm;
	// double kSVFit1Norm=1./((hSVFit1.Cross(tauSVFitminus_HRF.Vect().Unit())).Mag());
	// double kSVFit2Norm=1./((hSVFit2.Cross(tauSVFitplus_HRF.Vect().Unit())).Mag());
	// TVector3 kSVFit1 = (hSVFit1.Cross(tauSVFitminus_HRF.Vect().Unit()))*kSVFit1Norm;
	// TVector3 kSVFit2 = (hSVFit2.Cross(tauSVFitplus_HRF.Vect().Unit()))*kSVFit2Norm;



	// TLorentzVector testTruth(0,0,0,0);
	// // if((TruthDecayFromTau1!=TruthDecayFromTau2) && (TruthDecayFromTau1!=testTruth) && (TruthDecayFromTau2!=testTruth) && Ntp->CheckDecayID(5,5))
	// // 	{
	// SCalculator Scalc1Truth("a1");
	// SCalculator Scalc2Truth("a1");
		 
    
		  
	// vector<double> A1PionsChargeTruth_minus;
	// vector<double> A1PionsChargeTruth_plus;
		  
	// Pions1=Ntp->GetTruthPionsFromA1(0);
	// TruthDecayFromTau1=Pions1.at(0)+Pions1.at(1)+Pions1.at(2);
	// Pions2=Ntp->GetTruthPionsFromA1(1);
	// TruthDecayFromTau2=Pions2.at(0)+Pions2.at(1)+Pions2.at(2);
	      
	// A1PionsTruth_minus.push_back(Pions1.at(0));
	// A1PionsTruth_minus.push_back(Pions1.at(1)); 
	// A1PionsTruth_minus.push_back(Pions1.at(2));

	// A1PionsTruth_plus.push_back(Pions2.at(0));
	// A1PionsTruth_plus.push_back(Pions2.at(1)); 
	// A1PionsTruth_plus.push_back(Pions2.at(2));

		 
		  
	// A1PionsChargeTruth_minus.push_back(Ntp->MCTauandProd_charge(0, 2));
	// A1PionsChargeTruth_minus.push_back(Ntp->MCTauandProd_charge(0, 3));
	// A1PionsChargeTruth_minus.push_back(Ntp->MCTauandProd_charge(0, 4));
	// A1PionsChargeTruth_plus.push_back(Ntp->MCTauandProd_charge(1, 2));
	// A1PionsChargeTruth_plus.push_back(Ntp->MCTauandProd_charge(1, 3));
	// A1PionsChargeTruth_plus.push_back(Ntp->MCTauandProd_charge(1, 4));
		  
		 
	// Scalc1Truth.SortPions(A1PionsTruth_minus, A1PionsChargeTruth_minus);
	// Scalc2Truth.SortPions(A1PionsTruth_plus, A1PionsChargeTruth_plus);
		  		  

	// if(Ntp->MCTau_pdgid(0)==15)
	//   {
	//     tauandprodTruthminus.push_back(Ntp->GetTruthTauLV(5,0));
	//     tauandprodTruthminus.push_back(A1PionsTruth_minus.at(0));
	//     tauandprodTruthminus.push_back(A1PionsTruth_minus.at(1));
	//     tauandprodTruthminus.push_back(A1PionsTruth_minus.at(2));
	//     tauandprodTruthplus.push_back(Ntp->GetTruthTauLV(5,1));   
	//     tauandprodTruthplus.push_back(A1PionsTruth_plus.at(0));   
	//     tauandprodTruthplus.push_back(A1PionsTruth_plus.at(1));   
	//     tauandprodTruthplus.push_back(A1PionsTruth_plus.at(2));   
	//   }
	// else if (Ntp->MCTau_pdgid(0)==-15)
	//   {
	//     tauandprodTruthminus.push_back(Ntp->GetTruthTauLV(5,1));
	//     tauandprodTruthminus.push_back(A1PionsTruth_minus.at(0));
	//     tauandprodTruthminus.push_back(A1PionsTruth_minus.at(1));
	//     tauandprodTruthminus.push_back(A1PionsTruth_minus.at(2));
	//     tauandprodTruthplus.push_back(Ntp->GetTruthTauLV(5,0));   
	//     tauandprodTruthplus.push_back(A1PionsTruth_plus.at(0));   
	//     tauandprodTruthplus.push_back(A1PionsTruth_plus.at(1));   
	//     tauandprodTruthplus.push_back(A1PionsTruth_plus.at(2));   
	//   }
   
	if(a1a1)
	  {
	 
	    // ResolPullTauFroma1a1MZMomentum.at(t).Fill((TauplusPairConstraint.P()-Tauplustruth.P())/Tauplustruth.P(),w);
	    // ResolPullTauFroma1a1MZMomentum.at(t).Fill((TauminusPairConstraint.P()-Tauminustruth.P())/Tauminustruth.P(),w);
	    // ResolPullTauTauFroma1a1MZMomentum.at(t).Fill((ZPairConstraint.P()-TruthZ.P())/TruthZ.P(),w);
	    // ResolPullTauminusFroma1a1MZMomentum.at(t).Fill((TauminusPairConstraint.P()-Tauminustruth.P())/Tauminustruth.P(),w);
	    // ResolPullTauplusFroma1a1MZMomentum.at(t).Fill((TauplusPairConstraint.P()-Tauplustruth.P())/Tauplustruth.P(),w);
			      
			  
	    // AcolAngleTruth.at(t).Fill(acos((h1Truth*h1TruthNorm)*(h2Truth*h2TruthNorm)),Wspin);
			 
	    // if(((h1Truth*h1TruthNorm).Cross(h2Truth*h2TruthNorm))*tauminusTruth_HRF.Vect().Unit()<=0)polarimetricAcopAngleTruth.at(t).Fill(TMath::ATan2((k1Truth.Cross(k2Truth)).Mag(),k1Truth*k2Truth),Wspin);
	    // else polarimetricAcopAngleTruth.at(t).Fill(2*TMath::Pi()-TMath::ATan2((k1Truth.Cross(k2Truth)).Mag(),k1Truth*k2Truth),Wspin);
		      
	    //AcolAngle.at(t).Fill(acos(h1*h2),Wspin);
	      
			  
	    // AcolAngleSVFit.at(t).Fill(acos(hSVFit1*hSVFit2),Wspin);
	    // if(((hSVFit1.Cross(hSVFit2))*(tauSVFitminus_HRF.Vect().Unit()))<=0)polarimetricAcopAngleSVFit.at(t).Fill(acos(kSVFit1*kSVFit2),Wspin);
	    // else polarimetricAcopAngleSVFit.at(t).Fill(2.*TMath::Pi()-acos(kSVFit1*kSVFit2),Wspin);

	    // if(abs((tauandprodminus.at(0).P()-Tauminustruth.P())/Tauminustruth.P())<0.3 && abs((tauandprodplus.at(0).P()-Tauplustruth.P())/Tauplustruth.P())<0.3)
	    // 	{
	    // 	  if(((h1.Cross(h2))*(tauminus_HRF.Vect().Unit()))<=0)polarimetricAcopAngle30.at(t).Fill(acos(k1*k2),Wspin);
	    // 	  else polarimetricAcopAngle30.at(t).Fill(2.*TMath::Pi()-acos(k1*k2),Wspin);
	    // 	}
	    // if(abs((tauandprodminus.at(0).P()-Tauminustruth.P())/Tauminustruth.P())<0.25 && abs((tauandprodplus.at(0).P()-Tauplustruth.P())/Tauplustruth.P())<0.25)
	    // 	{
	    // 	  if(((h1.Cross(h2))*(tauminus_HRF.Vect().Unit()))<=0)polarimetricAcopAngle25.at(t).Fill(acos(k1*k2),Wspin);
	    // 	  else polarimetricAcopAngle25.at(t).Fill(2.*TMath::Pi()-acos(k1*k2),Wspin);
	    // 	}
	    // if(abs((tauandprodminus.at(0).P()-Tauminustruth.P())/Tauminustruth.P())<0.2 && abs((tauandprodplus.at(0).P()-Tauplustruth.P())/Tauplustruth.P())<0.2)
	    // 	{
	    // 	  if(((h1.Cross(h2))*(tauminus_HRF.Vect().Unit()))<=0)polarimetricAcopAngle20.at(t).Fill(acos(k1*k2),Wspin);
	    // 	  else polarimetricAcopAngle20.at(t).Fill(2.*TMath::Pi()-acos(k1*k2),Wspin);
	    // 	}
	    // if(abs((tauandprodminus.at(0).P()-Tauminustruth.P())/Tauminustruth.P())<0.15 && abs((tauandprodplus.at(0).P()-Tauplustruth.P())/Tauplustruth.P())<0.15)
	    // 	{
	    // 	  if(((h1.Cross(h2))*(tauminus_HRF.Vect().Unit()))<=0)polarimetricAcopAngle15.at(t).Fill(acos(k1*k2),Wspin);
	    // 	  else polarimetricAcopAngle15.at(t).Fill(2.*TMath::Pi()-acos(k1*k2),Wspin);
	    // 	}
	    // if(abs((tauandprodminus.at(0).P()-Tauminustruth.P())/Tauminustruth.P())<0.1 && abs((tauandprodplus.at(0).P()-Tauplustruth.P())/Tauplustruth.P())<0.1)
	    // 	{
	    // 	  if(((h1.Cross(h2))*(tauminus_HRF.Vect().Unit()))<=0)polarimetricAcopAngle10.at(t).Fill(acos(k1*k2),Wspin);
	    // 	  else polarimetricAcopAngle10.at(t).Fill(2.*TMath::Pi()-acos(k1*k2),Wspin);
	    // 	}
	    // if(abs((tauandprodminus.at(0).P()-Tauminustruth.P())/Tauminustruth.P())<0.05 && abs((tauandprodplus.at(0).P()-Tauplustruth.P())/Tauplustruth.P())<0.05)
	    // 	{
	    // 	  if(((h1.Cross(h2))*(tauminus_HRF.Vect().Unit()))<=0)polarimetricAcopAngle5.at(t).Fill(acos(k1*k2),Wspin);
	    // 	  else polarimetricAcopAngle5.at(t).Fill(2.*TMath::Pi()-acos(k1*k2),Wspin);
	    // 	}
	   
	    
	    PVXResol.at(t).Fill((Ntp->PVtx().X()-Ntp->PVtx_Gen().X())/Ntp->PVtx_Gen().X(),w);
	    PVXNoBSOldResol.at(t).Fill((Ntp->PVRefitNoBSOld().X()-Ntp->PVtx_Gen().X())/Ntp->PVtx_Gen().X(),w);
	    PVXNoBSTracksRemovedOldResol.at(t).Fill((tauNoBSTracksRemovedOldPrimaryVertex.X()-Ntp->PVtx_Gen().X())/Ntp->PVtx_Gen().X(),w);
	    PVXNoBSNewResol.at(t).Fill((tauNoBSNewPrimaryVertex.X()-Ntp->PVtx_Gen().X())/Ntp->PVtx_Gen().X(),w);
	    PVXBSOldResol.at(t).Fill((Ntp->PVRefitBSOld().X()-Ntp->PVtx_Gen().X())/Ntp->PVtx_Gen().X(),w);
	    PVXBSTracksRemovedOldResol.at(t).Fill((tauBSTracksRemovedOldPrimaryVertex.X()-Ntp->PVtx_Gen().X())/Ntp->PVtx_Gen().X(),w);
	    PVXBSNewResol.at(t).Fill((tauBSNewPrimaryVertex.X()-Ntp->PVtx_Gen().X())/Ntp->PVtx_Gen().X(),w);
	    
	    PVYResol.at(t).Fill((Ntp->PVtx().Y()-Ntp->PVtx_Gen().Y())/Ntp->PVtx_Gen().Y(),w);
	    PVYNoBSOldResol.at(t).Fill((Ntp->PVRefitNoBSOld().Y()-Ntp->PVtx_Gen().Y())/Ntp->PVtx_Gen().Y(),w);
	    PVYNoBSTracksRemovedOldResol.at(t).Fill((tauNoBSTracksRemovedOldPrimaryVertex.Y()-Ntp->PVtx_Gen().Y())/Ntp->PVtx_Gen().Y(),w);
	    PVYNoBSNewResol.at(t).Fill((tauNoBSNewPrimaryVertex.Y()-Ntp->PVtx_Gen().Y())/Ntp->PVtx_Gen().Y(),w);
	    PVYBSOldResol.at(t).Fill((Ntp->PVRefitBSOld().Y()-Ntp->PVtx_Gen().Y())/Ntp->PVtx_Gen().Y(),w);
	    PVYBSTracksRemovedOldResol.at(t).Fill((tauBSTracksRemovedOldPrimaryVertex.Y()-Ntp->PVtx_Gen().Y())/Ntp->PVtx_Gen().Y(),w);
	    PVYBSNewResol.at(t).Fill((tauBSNewPrimaryVertex.Y()-Ntp->PVtx_Gen().Y())/Ntp->PVtx_Gen().Y(),w);
	    
	    PVZResol.at(t).Fill((Ntp->PVtx().Z()-Ntp->PVtx_Gen().Z())/Ntp->PVtx_Gen().Z(),w);
	    PVZNoBSOldResol.at(t).Fill((Ntp->PVRefitNoBSOld().Z()-Ntp->PVtx_Gen().Z())/Ntp->PVtx_Gen().Z(),w);
	    PVZNoBSTracksRemovedOldResol.at(t).Fill((tauNoBSTracksRemovedOldPrimaryVertex.Z()-Ntp->PVtx_Gen().Z())/Ntp->PVtx_Gen().Z(),w);
	    PVZNoBSNewResol.at(t).Fill((tauNoBSNewPrimaryVertex.Z()-Ntp->PVtx_Gen().Z())/Ntp->PVtx_Gen().Z(),w);
	    PVZBSOldResol.at(t).Fill((Ntp->PVRefitBSOld().Z()-Ntp->PVtx_Gen().Z())/Ntp->PVtx_Gen().Z(),w);
	    PVZBSTracksRemovedOldResol.at(t).Fill((tauBSTracksRemovedOldPrimaryVertex.Z()-Ntp->PVtx_Gen().Z())/Ntp->PVtx_Gen().Z(),w);
	    PVZBSNewResol.at(t).Fill((tauBSNewPrimaryVertex.Z()-Ntp->PVtx_Gen().Z())/Ntp->PVtx_Gen().Z(),w);


	    
	    if(isRefitNoBSTracksRemovedOld)PVXNoBSTracksRemovedOldOnlyResol.at(t).Fill((tauNoBSTracksRemovedOldPrimaryVertex.X()-Ntp->PVtx_Gen().X())/Ntp->PVtx_Gen().X(),w);
	    if(isRefitNoBSNew)PVXNoBSNewOnlyResol.at(t).Fill((tauNoBSNewPrimaryVertex.X()-Ntp->PVtx_Gen().X())/Ntp->PVtx_Gen().X(),w);
	    if(isRefitBSTracksRemovedOld)PVXBSTracksRemovedOldOnlyResol.at(t).Fill((tauBSTracksRemovedOldPrimaryVertex.X()-Ntp->PVtx_Gen().X())/Ntp->PVtx_Gen().X(),w);
	    if(isRefitBSNew)PVXBSNewOnlyResol.at(t).Fill((tauBSNewPrimaryVertex.X()-Ntp->PVtx_Gen().X())/Ntp->PVtx_Gen().X(),w);

	    if(isRefitNoBSTracksRemovedOld)PVYNoBSTracksRemovedOldOnlyResol.at(t).Fill((tauNoBSTracksRemovedOldPrimaryVertex.Y()-Ntp->PVtx_Gen().Y())/Ntp->PVtx_Gen().Y(),w);
	    if(isRefitNoBSNew)PVYNoBSNewOnlyResol.at(t).Fill((tauNoBSNewPrimaryVertex.Y()-Ntp->PVtx_Gen().Y())/Ntp->PVtx_Gen().Y(),w);
	    if(isRefitBSTracksRemovedOld)PVYBSTracksRemovedOldOnlyResol.at(t).Fill((tauBSTracksRemovedOldPrimaryVertex.Y()-Ntp->PVtx_Gen().Y())/Ntp->PVtx_Gen().Y(),w);
	    if(isRefitBSNew)PVYBSNewOnlyResol.at(t).Fill((tauBSNewPrimaryVertex.Y()-Ntp->PVtx_Gen().Y())/Ntp->PVtx_Gen().Y(),w);


	    if(isRefitNoBSTracksRemovedOld)PVZNoBSTracksRemovedOldOnlyResol.at(t).Fill((tauNoBSTracksRemovedOldPrimaryVertex.Z()-Ntp->PVtx_Gen().Z())/Ntp->PVtx_Gen().Z(),w);
	    if(isRefitNoBSNew)PVZNoBSNewOnlyResol.at(t).Fill((tauNoBSNewPrimaryVertex.Z()-Ntp->PVtx_Gen().Z())/Ntp->PVtx_Gen().Z(),w);
	    if(isRefitBSTracksRemovedOld)PVZBSTracksRemovedOldOnlyResol.at(t).Fill((tauBSTracksRemovedOldPrimaryVertex.Z()-Ntp->PVtx_Gen().Z())/Ntp->PVtx_Gen().Z(),w);
	    if(isRefitBSNew)PVZBSNewOnlyResol.at(t).Fill((tauBSNewPrimaryVertex.Z()-Ntp->PVtx_Gen().Z())/Ntp->PVtx_Gen().Z(),w);
	    
	    
	    if(Ntp->MCTau_pdgid(0)==15)//nplus==1 && nminus==2)
	      {
	  
		ResolPullXVtxIna1a1.at(t).Fill((Ntp->PFTau_secondaryVertex_pos(Tauminus).X()-Ntp->MCTauandProd_Vertex(0,1).X())/Ntp->MCTauandProd_Vertex(0,1).X());
	  
		ResolPullXVtxIna1a1.at(t).Fill((Ntp->PFTau_secondaryVertex_pos(Tauplus).X()-Ntp->MCTauandProd_Vertex(1,1).X())/Ntp->MCTauandProd_Vertex(1,1).X());
	  
		ResolPullYVtxIna1a1.at(t).Fill((Ntp->PFTau_secondaryVertex_pos(Tauminus).Y()-Ntp->MCTauandProd_Vertex(0,1).Y())/Ntp->MCTauandProd_Vertex(0,1).Y());
		ResolPullYVtxIna1a1.at(t).Fill((Ntp->PFTau_secondaryVertex_pos(Tauplus).Y()-Ntp->MCTauandProd_Vertex(1,1).Y())/Ntp->MCTauandProd_Vertex(1,1).Y());
		ResolPullZVtxIna1a1.at(t).Fill((Ntp->PFTau_secondaryVertex_pos(Tauminus).Z()-Ntp->MCTauandProd_Vertex(0,1).Z())/Ntp->MCTauandProd_Vertex(0,1).Z());
		ResolPullZVtxIna1a1.at(t).Fill((Ntp->PFTau_secondaryVertex_pos(Tauplus).Z()-Ntp->MCTauandProd_Vertex(1,1).Z())/Ntp->MCTauandProd_Vertex(1,1).Z());
	      }
	    else if(Ntp->MCTau_pdgid(0)==-15)
	      {
		ResolPullXVtxIna1a1.at(t).Fill((Ntp->PFTau_secondaryVertex_pos(Tauminus).X()-Ntp->MCTauandProd_Vertex(1,1).X())/Ntp->MCTauandProd_Vertex(1,1).X());
		ResolPullXVtxIna1a1.at(t).Fill((Ntp->PFTau_secondaryVertex_pos(Tauplus).X()-Ntp->MCTauandProd_Vertex(0,1).X())/Ntp->MCTauandProd_Vertex(0,1).X());
		ResolPullYVtxIna1a1.at(t).Fill((Ntp->PFTau_secondaryVertex_pos(Tauminus).Y()-Ntp->MCTauandProd_Vertex(1,1).Y())/Ntp->MCTauandProd_Vertex(1,1).Y());
		ResolPullYVtxIna1a1.at(t).Fill((Ntp->PFTau_secondaryVertex_pos(Tauplus).Y()-Ntp->MCTauandProd_Vertex(0,1).Y())/Ntp->MCTauandProd_Vertex(0,1).Y());
		ResolPullZVtxIna1a1.at(t).Fill((Ntp->PFTau_secondaryVertex_pos(Tauminus).Z()-Ntp->MCTauandProd_Vertex(1,1).Z())/Ntp->MCTauandProd_Vertex(1,1).Z());
		ResolPullZVtxIna1a1.at(t).Fill((Ntp->PFTau_secondaryVertex_pos(Tauplus).Z()-Ntp->MCTauandProd_Vertex(0,1).Z())/Ntp->MCTauandProd_Vertex(0,1).Z());
	      }
	    // //ResolPullTauplusFroma1a1MeanEnergy.at(t).Fill((TauplusMean.E()-Tauplustruth.E())/Tauplustruth.E(),w);
	    // //ResolPullTauplusFroma1a1MZEnergy.at(t).Fill((TauplusPairConstraint.E()-Tauplustruth.E())/Tauplustruth.E(),w);
	    // ResolPullTauplusFroma1a1MeanMomentum.at(t).Fill((TauplusMean.P()-Tauplustruth.P())/Tauplustruth.P(),w);
	      

	    // //ResolPullTauminusFroma1a1MeanEnergy.at(t).Fill((TauminusMean.E()-Tauminustruth.E())/Tauminustruth.E(),w);
	    // //ResolPullTauminusFroma1a1MZEnergy.at(t).Fill((TauminusPairConstraint.E()-Tauminustruth.E())/Tauminustruth.E(),w);
	    // ResolPullTauminusFroma1a1MeanMomentum.at(t).Fill((TauminusMean.P()-Tauminustruth.P())/Tauminustruth.P(),w);
	      

	    // //ResolPullTauTauFroma1a1MeanEnergy.at(t).Fill((ZMean.E()-TruthZ.E())/TruthZ.E(),w);
	    // //ResolPullTauTauFroma1a1MZEnergy.at(t).Fill((ZPairConstraint.E()-TruthZ.E())/TruthZ.E(),w);
	    // ResolPullTauTauFroma1a1MeanMomentum.at(t).Fill((ZMean.P()-TruthZ.P())/TruthZ.P(),w);
	      

	    // //ResolPullTauFroma1a1MeanEnergy.at(t).Fill((TauminusMean.E()-Tauminustruth.E())/Tauminustruth.E(),w);
	    // //ResolPullTauFroma1a1MZEnergy.at(t).Fill((TauminusPairConstraint.E()-Tauminustruth.E())/Tauminustruth.E(),w);
	    // ResolPullTauFroma1a1MeanMomentum.at(t).Fill((TauminusMean.P()-Tauminustruth.P())/Tauminustruth.P(),w);
	    // //ResolPullTauFroma1a1MZMomentum.at(t).Fill((TauminusPairConstraint.P()-Tauminustruth.P())/Tauminustruth.P(),w);

	    // //ResolPullTauFroma1a1MeanEnergy.at(t).Fill((TauplusMean.E()-Tauplustruth.E())/Tauplustruth.E(),w);
	    // //ResolPullTauFroma1a1MZEnergy.at(t).Fill((TauplusPairConstraint.E()-Tauplustruth.E())/Tauplustruth.E(),w);
	    // ResolPullTauFroma1a1MeanMomentum.at(t).Fill((TauplusMean.P()-Tauplustruth.P())/Tauplustruth.P(),w);
	    // //ResolPullTauFroma1a1MZMomentum.at(t).Fill((TauplusPairConstraint.P()-Tauplustruth.P())/Tauplustruth.P(),w);
	      
	    // Phitrutha1a1.at(t).Fill(Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	    // // PhiSvFitResa1a1.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	    // // PhiVisResa1a1.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	  }
	// if(a1X)
	//   {
	//     if(a1minus)
	//       {
	// 	//ResolPullTauminusFroma1XMeanEnergy.at(t).Fill((TauminusMean.E()-Tauminustruth.E())/Tauminustruth.E(),w);
	// 	ResolPullTauminusFroma1XMeanMomentum.at(t).Fill((TauminusMean.P()-Tauminustruth.P())/Tauminustruth.P(),w);
	// 	//ResolPullXplusFroma1XMeanEnergy.at(t).Fill((TauplusMean.E()-Tauplustruth.E())/Tauplustruth.E(),w);
	// 	ResolPullXplusFroma1XMeanMomentum.at(t).Fill((TauplusMean.P()-Tauplustruth.P())/Tauplustruth.P(),w);
		    
	//       }
	//     else if (a1plus)
	//       {
	// 	//ResolPullTauplusFroma1XMeanEnergy.at(t).Fill((TauplusMean.E()-Tauplustruth.E())/Tauplustruth.E(),w);
	// 	ResolPullTauplusFroma1XMeanMomentum.at(t).Fill((TauplusMean.P()-Tauplustruth.P())/Tauplustruth.P(),w);
	// 	//ResolPullXminusFroma1XMeanEnergy.at(t).Fill((TauminusMean.E()-Tauminustruth.E())/Tauminustruth.E(),w);
	// 	ResolPullXminusFroma1XMeanMomentum.at(t).Fill((TauminusMean.P()-Tauminustruth.P())/Tauminustruth.P(),w);
		   
	//       }
	//   }

	if(id==10110333 || id==10210333) {	    
	  //Phitruthlpi.at(t).Fill(Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	  //  PhiSvFitReslpi.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	  // PhiVisReslpi.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	}
	if(id==10110433 || id==10210433) {	   
	  // Phitruthlrho.at(t).Fill(Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	  //  PhiSvFitReslrho.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	  //  PhiVisReslrho.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	}
	if(id==10130533 || id==10230533) {	   
	  //Phitruthla1.at(t).Fill(Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	  //  PhiSvFitResla1.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	  // PhiVisResla1.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	}
	Thetatruth.at(t).Fill(Tauminustruth.Theta(),w);
	    
	//DRTruth.at(t).Fill(Tauplustruth.DeltaR(Tauminustruth));
	//DRFull.at(t).Fill(Tauplustruth.DeltaR(Tauminustruth));
	// DRFullTruth.at(t).Fill(Tauplussvfit.DeltaR(Tauplustruth));
	//DRVisTruth.at(t).Fill(Tauplusvis.DeltaR(VisDecayfromTauplus));
	      
	// TauTauFullPtRes.at(t).Fill((Tauplussvfit+Tauminussvfit).Pt()-(Tauplustruth+Tauminustruth).Pt(),w);   
	// TauTauFullEtaRes.at(t).Fill((Tauplussvfit+Tauminussvfit).Eta()-(Tauplustruth+Tauminustruth).Eta(),w);
	// TauTauFullPhiRes.at(t).Fill(Ntp->DeltaPhi((Tauplussvfit+Tauminussvfit).Phi(),(Tauplustruth+Tauminustruth).Phi()),w);

	// TauTauVisPtRes.at(t).Fill((Tauplusvis+Tauminusvis).Pt()-(TruthDecayFromTau1 + TruthDecayFromTau2).Pt(),w);   
	//  TauTauVisEtaRes.at(t).Fill((Tauplusvis+Tauminusvis).Eta()-(TruthDecayFromTau1 + TruthDecayFromTau2).Eta(),w);
	//  TauTauVisPhiRes.at(t).Fill(Ntp->DeltaPhi((Tauplusvis+Tauminusvis).Phi(),(TruthDecayFromTau1 + TruthDecayFromTau2).Phi()),w);


	//TauTauFullPtResPull.at(t).Fill(((Tauplussvfit+Tauminussvfit).Pt()-(Tauplustruth+Tauminustruth).Pt())/(Tauplustruth+Tauminustruth).Pt(),w);   
	//TauTauFullEtaResPull.at(t).Fill(((Tauplussvfit+Tauminussvfit).Eta()-(Tauplustruth+Tauminustruth).Eta())/(Tauplustruth+Tauminustruth).Eta(),w);
	//TauTauFullPhiResPull.at(t).Fill(Ntp->DeltaPhi((Tauplussvfit+Tauminussvfit).Phi(),(Tauplustruth+Tauminustruth).Phi())/(Tauplustruth+Tauminustruth).Phi(),w);
	    
	// TauTauVisPtResPull.at(t).Fill(((Tauplusvis+Tauminusvis).Pt()-(TruthDecayFromTau1+TruthDecayFromTau2).Pt())/(TruthDecayFromTau1+TruthDecayFromTau2).Pt(),w);   
	//  TauTauVisEtaResPull.at(t).Fill(((Tauplusvis+Tauminusvis).Eta()-(TruthDecayFromTau1+TruthDecayFromTau2).Eta())/(TruthDecayFromTau1+TruthDecayFromTau2).Eta(),w);
	//  TauTauVisPhiResPull.at(t).Fill(Ntp->DeltaPhi((Tauplusvis+Tauminusvis).Phi(),(TruthDecayFromTau1+TruthDecayFromTau2).Phi())/(TruthDecayFromTau1+TruthDecayFromTau2).Phi(),w);


	//TauplusFullPtRes.at(t).Fill(Tauplussvfit.Pt()-Tauplustruth.Pt(),w);
	// TauplusFullEtaRes.at(t).Fill(Tauplussvfit.Eta()-Tauplustruth.Eta(),w);
	// TauplusFullPhiRes.at(t).Fill(Ntp->DeltaPhi(Tauplussvfit.Phi(),Tauplustruth.Phi()),w);
	    
	// TauminusFullPtRes.at(t).Fill(Tauminussvfit.Pt()-Tauminustruth.Pt(),w);
	// TauminusFullEtaRes.at(t).Fill(Tauminussvfit.Eta()-Tauminustruth.Eta(),w);
	// TauminusFullPhiRes.at(t).Fill(Ntp->DeltaPhi(Tauminussvfit.Phi(),Tauminustruth.Phi()),w);
	    
	//TauplusFullPtResPull.at(t).Fill((Tauplussvfit.Pt()-Tauplustruth.Pt())/Tauplustruth.Pt(),w);
	// TauplusFullEtaResPull.at(t).Fill((Tauplussvfit.Eta()-Tauplustruth.Eta())/Tauplustruth.Eta(),w);
	//  TauplusFullPhiResPull.at(t).Fill(Ntp->DeltaPhi(Tauplussvfit.Phi(),Tauplustruth.Phi())/Tauplustruth.Phi(),w);

	//TauminusFullPtResPull.at(t).Fill((Tauminussvfit.Pt()-Tauminustruth.Pt())/Tauminustruth.Pt(),w);
	//TauminusFullEtaResPull.at(t).Fill((Tauminussvfit.Eta()-Tauminustruth.Eta())/Tauminustruth.Eta(),w);
	//TauminusFullPhiResPull.at(t).Fill(Ntp->DeltaPhi(Tauminussvfit.Phi(),Tauminustruth.Phi())/Tauminustruth.Phi(),w);
	   

	// TauplusVisPtRes.at(t).Fill(Tauplusvis.Pt()-VisDecayfromTauplus.Pt(),w);
	// TauplusVisEtaRes.at(t).Fill(Tauplusvis.Eta()-VisDecayfromTauplus.Eta(),w);
	//TauplusVisPhiRes.at(t).Fill(Ntp->DeltaPhi(Tauplusvis.Phi(),VisDecayfromTauplus.Phi()),w);

	//TauminusVisPtRes.at(t).Fill(Tauminusvis.Pt()-VisDecayfromTauminus.Pt(),w);
	// TauminusVisEtaRes.at(t).Fill(Tauminusvis.Eta()-VisDecayfromTauminus.Eta(),w);
	// TauminusVisPhiRes.at(t).Fill(Ntp->DeltaPhi(Tauminusvis.Phi(),VisDecayfromTauminus.Phi()),w);

	//TauplusVisPtResPull.at(t).Fill((Tauplusvis.Pt()-VisDecayfromTauplus.Pt())/VisDecayfromTauplus.Pt(),w);
	//TauplusVisEtaResPull.at(t).Fill((Tauplusvis.Eta()-VisDecayfromTauplus.Eta())/VisDecayfromTauplus.Eta(),w);
	//TauplusVisPhiResPull.at(t).Fill(Ntp->DeltaPhi(Tauplusvis.Phi(),VisDecayfromTauplus.Phi())/VisDecayfromTauplus.Phi(),w);

	//TauminusVisPtResPull.at(t).Fill((Tauminusvis.Pt()-VisDecayfromTauminus.Pt())/VisDecayfromTauminus.Pt(),w);
	//TauminusVisEtaResPull.at(t).Fill((Tauminusvis.Eta()-VisDecayfromTauminus.Eta())/VisDecayfromTauminus.Eta(),w);
	// TauminusVisPhiResPull.at(t).Fill(Ntp->DeltaPhi(Tauminusvis.Phi(),VisDecayfromTauminus.Phi())/VisDecayfromTauminus.Phi(),w);
      
      }
  }
  //}
  //}
  //  This is a function if you want to do something after the event loop
  void HTauTau::Finish() {
  
    // if(mode == RECONSTRUCT) {
    //   std::cout<<" Starting Finish!  " <<std::endl;
    //   std::cout<<"A  Data  "<< NQCD.at(0).GetBinContent(1) << std::endl;
    //   std::cout<<"B  Data  "<< NQCD.at(0).GetBinContent(2) << std::endl;
    //   std::cout<<"C  Data  "<< NQCD.at(0).GetBinContent(3) << std::endl;
    //   std::cout<<"D  Data  "<< NQCD.at(0).GetBinContent(4) << std::endl;
    //   SkimConfig SC;
    //   SC.ApplySkimEfficiency(types,Npassed, Npassed_noweight);
    //   std::vector<double> QCD_Integral_B;
    //   double QCD_IntegralMC_B;
    //   double QCD_Integral_B_Data_minus_MC = 0;
    
    //   std::vector<double> QCD_Integral_C;
    //   double QCD_IntegralMC_C;
    //   double QCD_Integral_C_Data_minus_MC = 0;
    
    //   std::vector<double> QCD_Integral_D;
    //   double QCD_IntegralMC_D;
    //   double QCD_Integral_D_Data_minus_MC = 0;
    //   //Get Yields in ABCD for QCD Scalefactor                                                                                                                                                                  
    //   for(unsigned i=0;i<CrossSectionandAcceptance.size();i++){
    //     QCD_Integral_B.push_back(NQCD.at(i).GetBinContent(2));
    //     QCD_Integral_C.push_back(NQCD.at(i).GetBinContent(3));
    //     QCD_Integral_D.push_back(NQCD.at(i).GetBinContent(4));
    //     if(CrossSectionandAcceptance.at(i)>0){
    // 	QCD_Integral_B.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
    // 	QCD_Integral_C.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
    // 	QCD_Integral_D.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
    //     }
    //   }
    //   for(unsigned i=0;i<CrossSectionandAcceptance.size();i++){
    //     if(HConfig.GetID(i) == DataMCType::Data){
    // 	QCD_Integral_B_Data_minus_MC  += QCD_Integral_B.at(i);
    // 	QCD_Integral_C_Data_minus_MC += QCD_Integral_C.at(i);
    // 	QCD_Integral_D_Data_minus_MC += QCD_Integral_D.at(i);
    //     }
    //     if(CrossSectionandAcceptance.at(i)>0){
    // 	QCD_IntegralMC_B  += QCD_Integral_B.at(i);
    // 	QCD_IntegralMC_C  += QCD_Integral_C.at(i);
    // 	QCD_IntegralMC_D  += QCD_Integral_D.at(i);
    //     }
    //   }
    //   double CDFactor = (QCD_Integral_C_Data_minus_MC  - QCD_IntegralMC_C )/ (QCD_Integral_D_Data_minus_MC - QCD_IntegralMC_D);
    //   double QCD_Signal = QCD_Integral_B_Data_minus_MC *CDFactor;

    //   std::cout << "Factor: " << CDFactor << std::endl;
    //   std::cout << "QCD_Signal: " << QCD_Signal << std::endl;
    //   std::cout << "QCD in B region "<<  QCD_Integral_B_Data_minus_MC <<std::endl;
    //   std::cout << "QCD_Integral_B_Data_minus_MC is: " << QCD_Integral_B_Data_minus_MC << std::endl;
    //   std::cout << "QCD_Integral_C_Data_minus_MC is: " << QCD_Integral_C_Data_minus_MC << std::endl;
    //   std::cout << "QCD_Integral_D_Data_minus_MC is: " << QCD_Integral_D_Data_minus_MC << std::endl;
    //   std::cout << "QCD_IntegralMC_B is: " << QCD_IntegralMC_B << std::endl;
    //   std::cout << "QCD_IntegralMC_C is: " << QCD_IntegralMC_C << std::endl;
    //   std::cout << "QCD_IntegralMC_D is: " << QCD_IntegralMC_D << std::endl;
    //   ScaleAllHistOfType(HConfig.GetType(DataMCType::QCD),QCD_Signal/Nminus0.at(0).at(HConfig.GetType(DataMCType::QCD)).Integral());
    // }
    Selection::Finish();
  }
