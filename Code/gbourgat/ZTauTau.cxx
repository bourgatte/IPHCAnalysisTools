#include "ZTauTau.h"
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
#include "TauPolSoftware/TauDecaysInterface/interface/fonction_a1.h"



ZTauTau::ZTauTau(TString Name_, TString id_):
  Selection(Name_,id_),
  DataMC_Corr(true,true,false),
  tauTrgSF("tight")
{
  ChargeSumDummy = -999;
  selMuon_IsoDummy = 999.;
}

ZTauTau::~ZTauTau(){
  for(unsigned int j=0; j<Npassed.size(); j++){
    Logger(Logger::Info) << "Selection Summary before: "
			 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
			 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  ZTauTau::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==Trigger)             cut.at(Trigger)=1;
    if(i==Id_and_Kin)          cut.at(Id_and_Kin)=1;
    if(i==NPairsFound)         cut.at(NPairsFound)=1;
    if(i==Tau1Isolation)       cut.at(Tau1Isolation)=1.;
    if(i==Tau2Isolation)       cut.at(Tau2Isolation)=1.;
    if(i==LeptonVeto)          cut.at(LeptonVeto)=0.;
    if(i==PairCharge)          cut.at(PairCharge)=1.;
    if(i==PairMass)            cut.at(PairMass)=100.;
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
      title.at(i)="At least 1 good pair with Trig+Matching";
      hlabel="At least 1 good pair with Trig+Matching";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Trigger_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Trigger_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==Id_and_Kin){
      title.at(i)="Id and Kinematic";
      hlabel="Number of Event with good particles";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Id_and_Kin_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Id_and_Kin_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==NPairsFound){
      title.at(i)="Pairs with good DeltaR";
      hlabel="Pairs with good DeltaR";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NPairsFound_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NPairsFound_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==Tau1Isolation){
      title.at(i)="Tau1 Isolation";
      hlabel="Isolation of Tau1";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Tau1Isolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Tau1Isolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==Tau2Isolation){
      title.at(i)="Tau2 Isolation";
      hlabel="Isolation of Tau2";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Tau2Isolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Tau2Isolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==LeptonVeto){
      title.at(i)="Lepton Veto";
      hlabel="Third Lepton Veto  ";
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

  TauTauVisMass=HConfig.GetTH1D(Name+"_TauTauVisMass","Visible invariant mass of a tau pair",40,0,150," M(#tau#tau)_{vis}, GeV","Events");
  TauTauTruthMass=HConfig.GetTH1D(Name+"_TauTauTruthMass","Truth invariant mass of a tau pair",40,0,150," M(#tau#tau)_{truth}, GeV","Events");
  //TauTauFullMass=HConfig.GetTH1D(Name+"_TauTauFullMass","Full invariant mass of a tau pair",40,0,150," M(#tau#tau)_{full}, GeV","Events");

  NQCD=HConfig.GetTH1D(Name+"_NQCD","NQCD",4,0.5,4.5,"NQCD in ABCD","Events");
  // TauTauFullMass_B=HConfig.GetTH1D(Name+"_TauTauFullMass_B","TauTauFullMass_B",40,0,150,"#tau_h#tau_h SVFit Mass in B","Events");
  // TauTauFullMass_C=HConfig.GetTH1D(Name+"_TauTauFullMass_C","TauTauFullMass_C",40,0,150,"#tau_h#tau_h SVFit Mass in C","Events");
  // TauTauFullMass_D=HConfig.GetTH1D(Name+"_TauTauFullMass_D","TauTauFullMass_D",40,0,150,"#tau_h#tau_h SVFit Mass in D","Events");

  TauTauVisMass_B=HConfig.GetTH1D(Name+"_TauTauVisMass_B","TauTauVisMass_B",40,0,150,"#tau_h#tau_h Visible Mass in B","Events");
  TauTauVisMass_C=HConfig.GetTH1D(Name+"_TauTauVisMass_C","TauTauVisMass_C",40,0,150,"#tau_h#tau_h Visible Mass in C","Events");
  TauTauVisMass_D=HConfig.GetTH1D(Name+"_TauTauVisMass_D","TauTauVisMass_D",40,0,150,"#tau_h#tau_h Visible Mass in D","Events");
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
  
  ResolPullTauTauFroma1a1MeanEnergy=HConfig.GetTH1D(Name+"_ResolPullTauTauFroma1a1MeanEnergy","ResolPullTauTauFroma1a1MeanEnergy",30,-1,1,"","Events");
  ResolPullTauTauFroma1a1MZEnergy=HConfig.GetTH1D(Name+"_ResolPullTauTauFroma1a1MZEnergy","ResolPullTauTauFroma1a1MZEnergy",30,-1,1,"","Events");
  ResolPullTauTauFroma1a1MeanMomentum=HConfig.GetTH1D(Name+"_ResolPullTauTauFroma1a1MeanMomentum","ResolPullTauTauFroma1a1MeanMomentum",30,-1,1,"","Events");
  ResolPullTauTauFroma1a1MZMomentum=HConfig.GetTH1D(Name+"_ResolPullTauTauFroma1a1MZMomentum","ResolPullTauTauFroma1a1MZMomentum",30,-1,1,"","Events");

  ResolPullTauminusFroma1a1MeanEnergy=HConfig.GetTH1D(Name+"_ResolPullTauminusFroma1a1MeanEnergy","ResolPullTauminusFroma1a1MeanEnergy",30,-1,1,"","Events");
  ResolPullTauminusFroma1a1MZEnergy=HConfig.GetTH1D(Name+"_ResolPullTauminusFroma1a1MZEnergy","ResolPullTauminusFroma1a1MZEnergy",30,-1,1,"","Events");
  ResolPullTauminusFroma1a1MeanMomentum=HConfig.GetTH1D(Name+"_ResolPullTauminusFroma1a1MeanMomentum","ResolPullTauminusFroma1a1MeanMomentum",30,-1,1,"","Events");
  ResolPullTauminusFroma1a1MZMomentum=HConfig.GetTH1D(Name+"_ResolPullTauminusFroma1a1MZMomentum","ResolPullTauminusFroma1a1MZMomentum",30,-1,1,"","Events");


  ResolPullTauplusFroma1a1MeanEnergy=HConfig.GetTH1D(Name+"_ResolPullTauplusFroma1a1MeanEnergy","ResolPullTauplusFroma1a1MeanEnergy",30,-1,1,"","Events");
  ResolPullTauplusFroma1a1MZEnergy=HConfig.GetTH1D(Name+"_ResolPullTauplusFroma1a1MZEnergy","ResolPullTauplusFroma1a1MZEnergy",30,-1,1,"","Events");
  ResolPullTauplusFroma1a1MeanMomentum=HConfig.GetTH1D(Name+"_ResolPullTauplusFroma1a1MeanMomentum","ResolPullTauplusFroma1a1MeanMomentum",30,-1,1,"","Events");
  ResolPullTauplusFroma1a1MZMomentum=HConfig.GetTH1D(Name+"_ResolPullTauplusFroma1a1MZMomentum","ResolPullTauplusFroma1a1MZMomentum",30,-1,1,"","Events");



  ResolPullTauminusFroma1XMeanEnergy=HConfig.GetTH1D(Name+"_ResolPullTauminusFroma1XMeanEnergy","ResolPullTauminusFroma1XMeanEnergy",30,-1,1,"","Events");
  ResolPullTauminusFroma1XMeanMomentum=HConfig.GetTH1D(Name+"_ResolPullTauminusFroma1XMeanMomentum","ResolPullTauminusFroma1XMeanMomentum",30,-1,1,"","Events");
  ResolPullXminusFroma1XMeanEnergy=HConfig.GetTH1D(Name+"_ResolPullXminusFroma1XMeanEnergy","ResolPullXminusFroma1XMeanEnergy",30,-1,1,"","Events");
  ResolPullXminusFroma1XMeanMomentum=HConfig.GetTH1D(Name+"_ResolPullXminusFroma1XMeanMomentum","ResolPullXminusFroma1XMeanMomentum",30,-1,1,"","Events");


  ResolPullTauplusFroma1XMeanEnergy=HConfig.GetTH1D(Name+"_ResolPullTauplusFroma1XMeanEnergy","ResolPullTauplusFroma1XMeanEnergy",30,-1,1,"","Events");
  ResolPullTauplusFroma1XMeanMomentum=HConfig.GetTH1D(Name+"_ResolPullTauplusFroma1XMeanMomentum","ResolPullTauplusFroma1XMeanMomentum",30,-1,1,"","Events");
  ResolPullXplusFroma1XMeanEnergy=HConfig.GetTH1D(Name+"_ResolPullXplusFroma1XMeanEnergy","ResolPullXplusFroma1XMeanEnergy",30,-1,1,"","Events");
  ResolPullXplusFroma1XMeanMomentum=HConfig.GetTH1D(Name+"_ResolPullXplusFroma1XMeanMomentum","ResolPullXplusFroma1XMeanMomentum",30,-1,1,"","Events");



  tauminusa1a1MomentumVis=HConfig.GetTH1D(Name+"_tauminusa1a1MomentumVis","tauminusa1a1MomentumVis",15,20,200,"","Events");                             
  tauplusa1a1MomentumVis=HConfig.GetTH1D(Name+"_tauplusa1a1MomentumVis","tauplusa1a1MomentumVis",15,20,200,"","Events");
  InvariantMasstausa1a1Vis=HConfig.GetTH1D(Name+"_InvariantMasstausa1a1Vis","InvariantMasstausa1a1Vis",15,20,110,"","Events");
		    
  tauminusa1a1MomentumMean=HConfig.GetTH1D(Name+"_tauminusa1a1MomentumMean","tauminusa1a1MomentumMean",15,20,200,"","Events");                             
  tauplusa1a1MomentumMean=HConfig.GetTH1D(Name+"_tauplusa1a1MomentumMean","tauplusa1a1MomentumMean",15,20,200,"","Events");
  InvariantMasstausa1a1Mean=HConfig.GetTH1D(Name+"_InvariantMasstausa1a1Mean","InvariantMasstausa1a1Mean",15,0,200,"","Events");

  tauminusa1a1MomentumPairConstraint=HConfig.GetTH1D(Name+"_tauminusa1a1MomentumPairConstraint","tauminusa1a1MomentumPairConstraint",15,0,200,"","Events");       
  tauplusa1a1MomentumPairConstraint=HConfig.GetTH1D(Name+"_tauplusa1a1MomentumPairConstraint","tauplusa1a1MomentumPairConstraint",15,0,200,"","Events");
  InvariantMasstausa1a1PairConstraint=HConfig.GetTH1D(Name+"_InvariantMasstausa1a1PairConstraint","InvariantMasstausa1a1PairConstraint",15,20,160,"","Events");
	
  tauminusa1XMomentumVis=HConfig.GetTH1D(Name+"_tauminusa1XMomentumVis","tauminusa1XMomentumVis",15,0,200,"","Events");                             
  tauplusa1XMomentumVis=HConfig.GetTH1D(Name+"_tauplusa1XMomentumVis","tauplusa1XMomentumVis",15,0,200,"","Events");
  InvariantMasstausa1XVis=HConfig.GetTH1D(Name+"_InvariantMasstausa1XVis","InvariantMasstausa1XVis",15,20,110,"","Events");
		    
  tauminusa1XMomentumMean=HConfig.GetTH1D(Name+"_tauminusa1XMomentumMean","tauminusa1XMomentumMean",15,0,200,"","Events");                             
  tauplusa1XMomentumMean=HConfig.GetTH1D(Name+"_tauplusa1XMomentumMean","tauplusa1XMomentumMean",15,0,200,"","Events");
  InvariantMasstausa1XMean=HConfig.GetTH1D(Name+"_InvariantMasstausa1XMean","InvariantMasstausa1XMean",15,20,160,"","Events");

  
  Selection::ConfigureHistograms();   //   do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);  // do not remove
}

 

void  ZTauTau::Store_ExtraDist(){

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
  
  Extradist1d.push_back(&TauTauVisMass_B);
  Extradist1d.push_back(&TauTauVisMass_C);
  Extradist1d.push_back(&TauTauVisMass_D);
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
  // Extradist1d.push_back(&DRFullTruth);
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

  Extradist1d.push_back(&ResolPullTauTauFroma1a1MeanEnergy);
  Extradist1d.push_back(&ResolPullTauTauFroma1a1MZEnergy);
  Extradist1d.push_back(&ResolPullTauTauFroma1a1MeanMomentum);
  Extradist1d.push_back(&ResolPullTauTauFroma1a1MZMomentum);

  Extradist1d.push_back(&ResolPullTauminusFroma1a1MeanEnergy);
  Extradist1d.push_back(&ResolPullTauminusFroma1a1MZEnergy);
  Extradist1d.push_back(&ResolPullTauminusFroma1a1MeanMomentum);
  Extradist1d.push_back(&ResolPullTauminusFroma1a1MZMomentum);

  Extradist1d.push_back(&ResolPullTauplusFroma1a1MeanEnergy);
  Extradist1d.push_back(&ResolPullTauplusFroma1a1MZEnergy);
  Extradist1d.push_back(&ResolPullTauplusFroma1a1MeanMomentum);
  Extradist1d.push_back(&ResolPullTauplusFroma1a1MZMomentum);


  Extradist1d.push_back(&ResolPullTauminusFroma1XMeanEnergy);
  Extradist1d.push_back(&ResolPullTauminusFroma1XMeanMomentum);
  Extradist1d.push_back(&ResolPullXminusFroma1XMeanEnergy);
  Extradist1d.push_back(&ResolPullXminusFroma1XMeanMomentum);

  Extradist1d.push_back(&ResolPullTauplusFroma1XMeanEnergy);
  Extradist1d.push_back(&ResolPullTauplusFroma1XMeanMomentum); 
  Extradist1d.push_back(&ResolPullXplusFroma1XMeanEnergy);
  Extradist1d.push_back(&ResolPullXplusFroma1XMeanMomentum);


  Extradist1d.push_back(&tauminusa1a1MomentumVis);                        
  Extradist1d.push_back(&tauplusa1a1MomentumVis);
  Extradist1d.push_back(&InvariantMasstausa1a1Vis);
		    
  Extradist1d.push_back(&tauminusa1a1MomentumMean);
  Extradist1d.push_back(&tauplusa1a1MomentumMean);
  Extradist1d.push_back(&InvariantMasstausa1a1Mean);

  Extradist1d.push_back(&tauminusa1a1MomentumPairConstraint);
  Extradist1d.push_back(&tauplusa1a1MomentumPairConstraint);
  Extradist1d.push_back(&InvariantMasstausa1a1PairConstraint);
	
  Extradist1d.push_back(&tauminusa1XMomentumVis);
  Extradist1d.push_back(&tauplusa1XMomentumVis);
  Extradist1d.push_back(&InvariantMasstausa1XVis);
		    
  Extradist1d.push_back(&tauminusa1XMomentumMean);
  Extradist1d.push_back(&tauplusa1XMomentumMean);
  Extradist1d.push_back(&InvariantMasstausa1XMean);


}

void  ZTauTau::doEvent()  { //  Method called on every event

  unsigned int t;                // sample type, you may manage in your further analysis, if needed
  int id(Ntp->GetMCID());  //read event ID of a sample
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}  //  gives a warning if list of samples in Histo.txt  and SkimSummary.log do not coincide 
  //  std::cout<<"------------------ New Event -----------------------"<<std::endl;
  Charge = ChargeSumDummy;

  bool trig=0;
  std::vector<int> TauIndex ;
  std::vector<int> TriggerIndexVector ;
  std::vector<TString>  MatchedTriggerNames;
  value.at(Trigger)=0;
  MatchedTriggerNames.push_back("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v");
  MatchedTriggerNames.push_back("HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v");
  TriggerIndexVector=Ntp->GetVectorTriggers(MatchedTriggerNames);
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
  for(unsigned int itrig = 0; itrig < TriggerIndexVector.size(); itrig++){
    if(Ntp->TriggerAccept(TriggerIndexVector.at(itrig))){
      trig=1;
    }
  }
  for(unsigned int iDaughter=0;   iDaughter  <  Ntp->NDaughters() ;iDaughter++ ) {
    if(Ntp->CHECK_BIT(Ntp->Daughters_trgMatched(iDaughter),29) || Ntp->CHECK_BIT(Ntp->Daughters_trgMatched(iDaughter),27))
      {
	TauIndex.push_back(iDaughter);
      }
  }
  if(trig && TauIndex.size()>1)value.at(Trigger)=1;
  pass.at(Trigger)=(value.at(Trigger)==cut.at(Trigger));
  value.at(Id_and_Kin)=0;
  int goodTau_counter=0;
  std::vector<int> thirdLeptonCounter;
  std::vector<int> goodTausIndex;
  for(unsigned int iDaughter=0;   iDaughter  < TauIndex.size()  ;iDaughter++ ) {
    if(Ntp->tauBaselineSelection(iDaughter,36., 2.1, 0,0)){
      goodTausIndex.push_back(TauIndex.at(iDaughter));
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
  value.at(NPairsFound)=0;

  for(int unsigned ipair =0; ipair <goodTausIndex.size(); ipair++)
    {
      for(int unsigned jpair =1; jpair <goodTausIndex.size(); jpair++)
	{
	  if(jpair>ipair)
	    {
	      if((Ntp->Daughters_P4(goodTausIndex.at(ipair)).DeltaR(Ntp->Daughters_P4(goodTausIndex.at(jpair))))>0.5){
		PairsIndexTemp.push_back(j);
		PairsIndexTau1Temp.push_back(goodTausIndex.at(ipair));
		PairsIndexTau2Temp.push_back(goodTausIndex.at(jpair));
		j++;
	      }
	    }
	}
    }
  if(PairsIndexTemp.size()>0)value.at(NPairsFound)=1;
  pass.at(NPairsFound)=(value.at(NPairsFound)==cut.at(NPairsFound));
  if(pass.at(NPairsFound))
    {
      Sorted = Ntp->SortPair(PairsIndexTemp,PairsIndexTau1Temp,PairsIndexTau2Temp);
      Tau1=PairsIndexTau1Temp.at(Sorted.back());
      Tau2=PairsIndexTau2Temp.at(Sorted.back());
      value.at(Tau1Isolation)=0;
      value.at(Tau1Isolation) = (Ntp->isIsolatedTau(Tau1,"Tight"));
      pass.at(Tau1Isolation) = value.at(Tau1Isolation);
      value.at(Tau2Isolation)=0;
      value.at(Tau2Isolation) = (Ntp->isIsolatedTau(Tau2,"Tight"));
      pass.at(Tau2Isolation) = value.at(Tau2Isolation);

      
      value.at(LeptonVeto)=0;
      for(unsigned int iDaughter=0;   iDaughter  <  Ntp->NDaughters() ;iDaughter++ ) {  // loop over all daughters in the event
	if((iDaughter!=Tau1)&&(iDaughter!=Tau2)){
	  if(Ntp->ElectronVeto(iDaughter) || Ntp->MuonVeto(iDaughter))thirdLeptonCounter.push_back(iDaughter);
	}
      }
      value.at(LeptonVeto) = thirdLeptonCounter.size()>0;
      pass.at(LeptonVeto) = (value.at(LeptonVeto)==cut.at(LeptonVeto));

      value.at(PairCharge)=0;
      bool isOS=false;
      isOS=((Ntp->Daughters_charge(Tau1)/abs(Ntp->Daughters_charge(Tau1))) != (Ntp->Daughters_charge(Tau2)/abs(Ntp->Daughters_charge(Tau2))));
      if(isOS)value.at(PairCharge) = 1;
      pass.at(PairCharge) = value.at(PairCharge);
      value.at(PairMass) = 999.;
      //value.at(MTM) = 999.;
      //value.at(MTM) = .;
      value.at(PairMass)=(Ntp->Daughters_P4(Tau1)+Ntp->Daughters_P4(Tau2)).M();
      pass.at(PairMass) = (value.at(PairMass) < cut.at(PairMass));
      //pass.at(MTM) = (value.at(MTM) <= cut.at(MTM));
    }
  // Here you can defined different type of weights you want to apply to events.
  double wobs=1;
  double w=1;
  if(!Ntp->isData() && id!=DataMCType::QCD) {
    //    w *= reweight.weight(2016,26,Ntp->PUNumInteractions());
    w *= reweight.PUweightHTT(Ntp->npu());
    //std::cout<<" pu weigh HTT  "<< reweight.PUweightHTT(Ntp->npu())<<std::endl;
    if(!Ntp->isData() && pass.at(NPairsFound) ){
      double w1 = tauTrgSF.getSF(Ntp->TauP4_Corrected(Tau1).Pt(),  Ntp->decayMode(Tau1)) ;  //from Luca
      double w2 = tauTrgSF.getSF(Ntp->TauP4_Corrected(Tau2).Pt(),  Ntp->decayMode(Tau2)) ;
      w*=w1;
      w*=w2;
    }
    if(!Ntp->isData() && pass.at(NPairsFound) && (id==33 || id == 10110333 || id == 10110433|| id == 10130533|| id ==10210333|| id == 10210433|| id == 10230533|| id ==10310333 || id ==10330533 || id ==10410433 || id == 10410333|| id == 10430533|| id == 30530533)){
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
  if(id == 33){
    if(pass.at(NPairsFound)){
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


  w*=Ntp->stitch_weight();


  std::vector<unsigned int> exclude_cuts;
  exclude_cuts.push_back(Tau1Isolation);
  exclude_cuts.push_back(Tau2Isolation);
  exclude_cuts.push_back(PairCharge);
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
    if(pass.at(PairCharge)){
      if(pass.at(Tau1Isolation) && pass.at(Tau2Isolation) ){
	NQCD.at(t).Fill(1.,w); //A
      }
      if((Ntp->isIsolatedTau(Tau1,"Medium") && !Ntp->isIsolatedTau(Tau2,"Tight") && Ntp->isIsolatedTau(Tau2,"Loose")) || (Ntp->isIsolatedTau(Tau2,"Medium") && !Ntp->isIsolatedTau(Tau1,"Tight") && Ntp->isIsolatedTau(Tau1,"Loose"))){
	NQCD.at(t).Fill(2.,w); //B
	TauTauVisMass_B.at(t).Fill((Tau1P4+Tau2P4).M(),w);
	//TauTauFullMass_B.at(t).Fill((tau1P4+tau2P4).M(),w);
      }
    }
    if(!pass.at(PairCharge)){
      if(pass.at(Tau1Isolation) && pass.at(Tau2Isolation)){
	NQCD.at(t).Fill(3.,w); //C
	TauTauVisMass_C.at(t).Fill((Tau1P4+Tau2P4).M(),w);
	//TauTauFullMass_C.at(t).Fill((tau1P4+tau2P4).M(),w);
      }
      if((Ntp->isIsolatedTau(Tau1,"Medium") && !Ntp->isIsolatedTau(Tau2,"Tight") && Ntp->isIsolatedTau(Tau2,"Loose")) || (Ntp->isIsolatedTau(Tau2,"Medium") && !Ntp->isIsolatedTau(Tau1,"Tight") && Ntp->isIsolatedTau(Tau1,"Loose"))){
	NQCD.at(t).Fill(4.,w); //D
	TauTauVisMass_D.at(t).Fill((Tau1P4+Tau2P4).M(),w);
	//TauTauFullMass_D.at(t).Fill((tau1P4+tau2P4).M(),w);
      }
    }
  }
  bool IsQCDEvent = false;
  if(passAllBut(exclude_cuts)){
    if(pass.at(PairCharge)){
      if((Ntp->isIsolatedTau(Tau1,"Medium") && !Ntp->isIsolatedTau(Tau2,"Tight") && Ntp->isIsolatedTau(Tau2,"Loose")) || (Ntp->isIsolatedTau(Tau2,"Medium") && !Ntp->isIsolatedTau(Tau1,"Tight") && Ntp->isIsolatedTau(Tau1,"Loose"))){
	if(id == DataMCType::Data){
	  QCDShape.at(t).Fill(1,w);
	  t=HConfig.GetType(DataMCType::QCD);
	  IsQCDEvent = true;
	}
      }
    }
  }

  if(IsQCDEvent){ pass.at(PairCharge)= true;pass.at(Tau2Isolation)= true;pass.at(Tau1Isolation)=true;}
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
    double pvx(0);
    pvx =  Ntp->npv();
    // if(id == DataMCType::Data) pvx =  Ntp->npv();
    if(id !=DataMCType::Data && id !=DataMCType::QCD)	  pvx = Ntp->PUNumInteractions();
    NPrimeVtx.at(t).Fill(pvx,w);
    NPU.at(t).Fill(Ntp->npu(),w);
    RHO.at(t).Fill(Ntp->rho(),w);
    
    
    std::vector<int> thirdLepton;

    //---------  svfit ---------------------
    std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;
    classic_svFit::MeasuredTauLepton lep1(1, Tau1P4.Pt(), Tau1P4.Eta(),  Tau1P4.Phi(), Tau1P4.M());
    classic_svFit::MeasuredTauLepton lep2(1, Tau2P4.Pt(), Tau2P4.Eta(),  Tau2P4.Phi(), Tau2P4.M());
 
    measuredTauLeptons.push_back(lep1);
    measuredTauLeptons.push_back(lep2);
    TMatrixD metcov(2,2);
    double metx = Ntp->MET()*cos(Ntp->METphi());
    double mety = Ntp->MET()*sin(Ntp->METphi());

    metcov[0][0] = Ntp->PFMETCov00();
    metcov[1][0] = Ntp->PFMETCov01();
    metcov[0][1] = Ntp->PFMETCov10();
    metcov[1][1] = Ntp->PFMETCov11();

    svfitAlgo1.addLogM_fixed(true,5.0);
    svfitAlgo1.setDiTauMassConstraint(-1.0);
    svfitAlgo1.integrate(measuredTauLeptons,metx,mety, metcov );
    if(svfitAlgo1.isValidSolution()){
      double higgsmass  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svfitAlgo1.getHistogramAdapter())->getMass();
      h_SVFitMass.at(t).Fill(higgsmass,w); 
    }
    ClassicSVfit svfitAlgo2;
    svfitAlgo2.setHistogramAdapter(new classic_svFit::TauTauHistogramAdapter());
    svfitAlgo2.addLogM_fixed(true, 5.);
    svfitAlgo2.integrate(measuredTauLeptons,metx,mety, metcov );
    tau1P4 = static_cast<classic_svFit::TauTauHistogramAdapter*>(svfitAlgo2.getHistogramAdapter())->GetFittedTau1LV();
    tau2P4 = static_cast<classic_svFit::TauTauHistogramAdapter*>(svfitAlgo2.getHistogramAdapter())->GetFittedTau2LV();

    //---------  svfit ---------------------


    //    TLorentzVector taunew(tau1P4.Px(), );

    // svfTau1E.at(t).Fill(tau1P4.E(),w);
    // svfTau2E.at(t).Fill(tau2P4.E(),w);
    
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
    
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> Tauplussvfit;
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> Tauminussvfit;
    TLorentzVector Tauplusvis;
    TLorentzVector Tauminusvis;
    TLorentzVector Pi0RECO;
    TLorentzVector Tauplustruth;
    TLorentzVector Tauminustruth;
    unsigned int Tauplus= 0;
    unsigned int Tauminus=0;

    if(Ntp->Daughters_charge(Tau1)>0)
      {
	Tauplussvfit=tau1P4;
	Tauminussvfit=tau2P4;
	Tauplusvis=Tau1P4;
	Tauminusvis=Tau2P4;
	Tauplus=Tau1;
	Tauminus=Tau2;
      }
    else
      {
	Tauplussvfit=tau2P4;
	Tauminussvfit=tau1P4;
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

    /*
      if(Ntp->decayMode(Tau1)==0 && Ntp->decayMode(Tau2)==0 ){PhiDatavispipi.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),w);}
      if((Ntp->decayMode(Tau1)==0 && Ntp->decayMode(Tau2)==1) || (Ntp->decayMode(Tau2)==0 && Ntp->decayMode(Tau1)==1)){PhiDatavispirho.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),w);}
      if((Ntp->decayMode(Tau1)==0 && Ntp->decayMode(Tau2)==10) || (Ntp->decayMode(Tau2)==0 && Ntp->decayMode(Tau1)==10)){PhiDatavispia1.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),w);}
      if(Ntp->decayMode(Tau1)==1 && Ntp->decayMode(Tau2)==1 ){PhiDatavisrhorho.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),w);}
      if((Ntp->decayMode(Tau1)==1 && Ntp->decayMode(Tau2)==10) || (Ntp->decayMode(Tau2)==1 && Ntp->decayMode(Tau1)==10)){PhiDatavisrhoa1.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),w);}
      if(Ntp->decayMode(Tau1)==10 && Ntp->decayMode(Tau2)==10 ){PhiDatavisa1a1.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),w);}*/

    
    TVector3 tauPrimaryVertex , TauminusSecondaryVertex , TauplusSecondaryVertex;
    TVector3 TauminusDirection , TauplusDirection;
    double thetaGJ_Tauminus , thetaGJ_Tauplus;
    TLorentzVector a1LV_Tauminus , a1LV_Tauplus, piLV_Tauminus , piLV_Tauplus, rhoLV_Tauminus , rhoLV_Tauplus;
    TLorentzVector TauminusSmall, TauminusLarge, TauminusMean, TauminusPairConstraint;
    TLorentzVector TauplusSmall, TauplusLarge, TauplusMean, TauplusPairConstraint;
    bool isPlusReal, isMinusReal, a1X=false, a1a1=false, piminus=false, piplus=false, rhominus=false, rhoplus=false, a1minus=false, a1plus=false;
    std::vector<TLorentzVector> solutions, TauMinusSolutions, TauPlusSolutions ;
    TLorentzVector ZMean,ZMeana1X;
    TLorentzVector ZPairConstraint;

    if(Ntp->decayMode(Tauminus) == 0)piminus=true;
    if(Ntp->decayMode(Tauplus) == 0)piplus=true;
    if(Ntp->decayMode(Tauminus) == 1)rhominus=true;
    if(Ntp->decayMode(Tauplus) == 1)rhoplus=true;
    if(Ntp->decayMode(Tauminus) == 10 &&  Ntp->PFTau_hassecondaryVertex(Tauminus) && Ntp->PFtauHasThreePions(Tauminus))a1minus=true;
    if(Ntp->decayMode(Tauplus) == 10 && Ntp->PFTau_hassecondaryVertex(Tauplus) && Ntp->PFtauHasThreePions(Tauplus))a1plus=true;

    if((a1minus && !a1plus ) ||( !a1minus && a1plus)) a1X=true;  //a1-X
    if(a1minus && a1plus ) a1a1=true;  //a1-a1

    if(a1a1 && id==30530533)
      {

	tauPrimaryVertex = Ntp->PVtx();
	TauminusSecondaryVertex = Ntp->PFTau_secondaryVertex_pos(Tauminus);
	TauplusSecondaryVertex = Ntp->PFTau_secondaryVertex_pos(Tauplus);
	TauminusDirection = TauminusSecondaryVertex - tauPrimaryVertex;
	TauplusDirection = TauplusSecondaryVertex - tauPrimaryVertex;

	a1LV_Tauminus = Ntp->PFTau_PionsP4(Tauminus,0) + Ntp->PFTau_PionsP4(Tauminus,1) + Ntp->PFTau_PionsP4(Tauminus,2);
	a1LV_Tauplus = Ntp->PFTau_PionsP4(Tauplus,0) + Ntp->PFTau_PionsP4(Tauplus,1) + Ntp->PFTau_PionsP4(Tauplus,2); 

	isMinusReal = isPhysical(TauminusDirection, a1LV_Tauminus);
	isPlusReal = isPhysical(TauplusDirection, a1LV_Tauplus);
	
	IsPairPhysical.at(t).Fill(isPlusReal && isMinusReal,w);

	solutions=tauPairMomentumSolutions(TauminusDirection, a1LV_Tauminus, isPlusReal, TauplusDirection, a1LV_Tauplus, isMinusReal);
	TauminusSmall = solutions.at(0); 
	TauminusLarge = solutions.at(1); 
	TauminusMean = solutions.at(2); 
	TauminusPairConstraint = solutions.at(3);

	TauplusSmall = solutions.at(4);
	TauplusLarge = solutions.at(5); 
	TauplusMean = solutions.at(6);
	TauplusPairConstraint = solutions.at(7);
	if(true/*isPlusReal && isMinusReal*/)
	  {
	    tauminusa1a1MomentumVis.at(t).Fill(Tauminusvis.P(),w);                                  
	    tauplusa1a1MomentumVis.at(t).Fill(Tauplusvis.P(),w);
	    InvariantMasstausa1a1Vis.at(t).Fill((Tauminusvis + Tauplusvis).M(),w);

	    tauminusa1a1MomentumMean.at(t).Fill(TauminusMean.P(),w);
	    tauplusa1a1MomentumMean.at(t).Fill(TauplusMean.P(),w);
	    InvariantMasstausa1a1Mean.at(t).Fill((TauminusMean + TauplusMean).M(),w);

	    tauminusa1a1MomentumPairConstraint.at(t).Fill(TauminusPairConstraint.P(),w);           
	    tauplusa1a1MomentumPairConstraint.at(t).Fill(TauplusPairConstraint.P(),w);
	    InvariantMasstausa1a1PairConstraint.at(t).Fill( (TauminusPairConstraint + TauplusPairConstraint).M() ,w);
	  
	    double widthSolutions_Tauminus = fabs(TauminusLarge.P() - TauminusSmall.P());
	    double widthSolutions_Tauplus = fabs(TauplusLarge.P() - TauplusSmall.P());

	    ZMean=TauminusMean+TauplusMean;
	    ZPairConstraint= TauplusPairConstraint+TauminusPairConstraint;
	  }
      }
    if(a1X && (id==10130533 || id==10230533 || id==10330533|| id==10430533))
      {

	isPlusReal=false;
	isMinusReal=false;

	
	tauPrimaryVertex = Ntp->PVtx();
	if(a1minus)
	  {
	    isMinusReal=true;//isPhysical(TauminusDirection, a1LV_Tauminus);
	    if(isMinusReal)
	      {
		TauminusSecondaryVertex = Ntp->PFTau_secondaryVertex_pos(Tauminus);
		TauminusDirection = TauminusSecondaryVertex - tauPrimaryVertex;
		a1LV_Tauminus = Ntp->PFTau_PionsP4(Tauminus,0) + Ntp->PFTau_PionsP4(Tauminus,1) + Ntp->PFTau_PionsP4(Tauminus,2);
		TauMinusSolutions = tauMomentumSolutions(TauminusDirection,a1LV_Tauminus,isPlusReal);
		TauminusMean=TauMinusSolutions.at(2);
		InvariantMasstausa1XMean.at(t).Fill((TauminusMean + Tauplusvis).M(),w);
		tauminusa1XMomentumMean.at(t).Fill(TauminusMean.P(),w);
		tauminusa1XMomentumVis.at(t).Fill(Tauminusvis.P(),w);
		//TauplusSecondaryVertex = Ntp->PFTau_secondaryVertex_pos(Tauplus);
		//TauplusDirection = Tauplusvis;
		//a1LV_Tauplus = Ntp->PFTau_PionsP4(Tauplus,0) + Ntp->PFTau_PionsP4(Tauplus,1) + Ntp->PFTau_PionsP4(Tauplus,2);
		//TauPlusSolutions = tauMomentumSolutions(Tauplusvis,a1LV_Tauplus,isMinusReal);
		//TauplusMean=TauPlusSolutions.at(2);
		//TauplusMean.Print();
		//InvariantMasstausa1XMean.at(t).Fill((TauplusMean + Tauminusvis).M(),w);
		//tauplusa1XMomentumMean.at(t).Fill(TauplusDIrection.P(),w);
		tauplusa1XMomentumVis.at(t).Fill(Tauplusvis.P(),w);
		ZMeana1X=Tauplusvis+TauminusMean;
	      }
	  }
	else if(a1plus)
	  {
	    isPlusReal=true; //isPhysical(TauplusDirection, a1LV_Tauplus);
	    if(isPlusReal)
	      {
		
		//TauminusSecondaryVertex = Ntp->PFTau_secondaryVertex_pos(Tauminus);
		//TauminusDirection = TauminusSecondaryVertex - tauPrimaryVertex;
		//a1LV_Tauminus = Ntp->PFTau_PionsP4(Tauminus,0) + Ntp->PFTau_PionsP4(Tauminus,1) + Ntp->PFTau_PionsP4(Tauminus,2);
		//TauMinusSolutions = tauMomentumSolutions(TauminusDirection,a1LV_Tauminus,isPlusReal);
		//TauminusMean=TauMinusSolutions.at(2);
		//InvariantMasstausa1XMean.at(t).Fill((TauminusMean + Tauplusvis).M(),w);
		//tauminusa1XMomentumMean.at(t).Fill(TauminusMean.P(),w);
		tauminusa1XMomentumVis.at(t).Fill(Tauminusvis.P(),w);
		TauplusSecondaryVertex = Ntp->PFTau_secondaryVertex_pos(Tauplus);
		TauplusDirection = TauplusSecondaryVertex - tauPrimaryVertex;
		a1LV_Tauplus = Ntp->PFTau_PionsP4(Tauplus,0) + Ntp->PFTau_PionsP4(Tauplus,1) + Ntp->PFTau_PionsP4(Tauplus,2);
		TauPlusSolutions = tauMomentumSolutions(TauplusDirection,a1LV_Tauplus,isMinusReal);
		TauplusMean=TauPlusSolutions.at(2);
		InvariantMasstausa1XMean.at(t).Fill((TauplusMean + Tauminusvis).M(),w);
		tauplusa1XMomentumMean.at(t).Fill(TauplusMean.P(),w);
		tauplusa1XMomentumVis.at(t).Fill(Tauplusvis.P(),w);
		ZMeana1X=Tauminusvis+TauplusMean;
	      }
	  }
	InvariantMasstausa1XVis.at(t).Fill((Tauminusvis + Tauplusvis).M(),w);
	


	// if(piminus)
	//   {
	//     piLV_Tauminus = Ntp->PFTau_PionsP4(Tauminus,0);
	//   }
	// if(piplus)
	//   {
	//     piLV_Tauplus = Ntp->PFTau_PionsP4(Tauplus,0);
	//   }
	// if(rhominus)
	//   {
	//     rhoLV_Tauminus = Ntp->PFTau_PionsP4(Tauminus,0)+ Ntp->PFTau_PionsP4(Tauminus,1);
	//   }
	// if(rhoplus)
	//   {
	//     rhoLV_Tauplus = Ntp->PFTau_PionsP4(Tauplus,0)+ Ntp->PFTau_PionsP4(Tauplus,1);
	//   }
	
		    

      }
    if(id == 10110333 || id == 10110433|| id == 10130533|| id ==10210333|| id == 10210433|| id == 10230533|| id ==10310333 || id ==10330533 || id ==10410433 || id == 10410333|| id == 10430533|| id == 30530533)
      {
	TLorentzVector Tau1Truth; 
	TLorentzVector Tau2Truth;
	TLorentzVector TruthDecayFromTau1;
	TLorentzVector TruthDecayFromTau2; 
	std::vector<TLorentzVector> Pions1;
	std::vector<TLorentzVector> Pions2;
	bool decay=0;
	if(Ntp->CheckDecayID(1,3)){
	  Tau1Truth=Ntp->GetTruthTauLV(1,0);
	  Tau2Truth=Ntp->GetTruthTauLV(3,1);
	  TruthDecayFromTau1=Ntp->GetTruthTauProductLV(1,11,0);
	  TruthDecayFromTau2=Ntp->GetTruthTauProductLV(3,211,1);decay=1;if(TruthDecayFromTau1==TruthDecayFromTau2)cout<<"1 3"<<endl;
	  
	}
	if(Ntp->CheckDecayID(1,4)){
	  Tau1Truth=Ntp->GetTruthTauLV(1,0);
	  Tau2Truth=Ntp->GetTruthTauLV(4,1);
	  TruthDecayFromTau1=Ntp->GetTruthTauProductLV(1,11,0);
	  Pions2=Ntp->GetTruthPionsFromRho(1);
	  TruthDecayFromTau2=Pions2.at(0)+Pions2.at(1);decay=1;if(TruthDecayFromTau1==TruthDecayFromTau2)cout<<"1 4"<<endl;
	  if((sqrt((Tau1P4.Eta()-Tau1Truth.Eta())*(Tau1P4.Eta()-Tau1Truth.Eta())+(Tau1P4.Phi()-Tau1Truth.Phi())*(Tau1P4.Phi()-Tau1Truth.Phi())))<0.5 && (sqrt((Tau2P4.Eta()-Tau2Truth.Eta())*(Tau2P4.Eta()-Tau2Truth.Eta())+(Tau2P4.Phi()-Tau2Truth.Phi())*(Tau2P4.Phi()-Tau2Truth.Phi())))<0.5)
	    {
	      Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau2).E()-Pions2.at(1).E(),w);
	      Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau2).E()-Pions2.at(1).E())/Pions2.at(1).E(),w);
	    }
	  else{ Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau1).E()-Pions2.at(1).E(),w);Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau1).E()-Pions2.at(1).E())/Pions2.at(1).E(),w);}
	}
	if(Ntp->CheckDecayID(1,5)){
	  Tau1Truth=Ntp->GetTruthTauLV(1,0);
	  Tau2Truth=Ntp->GetTruthTauLV(5,1);
	  TruthDecayFromTau1=Ntp->GetTruthTauProductLV(1,11,0);
	  Pions2=Ntp->GetTruthPionsFromA1(1);
	  TruthDecayFromTau2=Pions2.at(0)+Pions2.at(1)+Pions2.at(2);decay=1;if(TruthDecayFromTau1==TruthDecayFromTau2)cout<<"1 5"<<endl;

	}


	if(Ntp->CheckDecayID(2,3)){
	  Tau1Truth=Ntp->GetTruthTauLV(2,0);
	  Tau2Truth=Ntp->GetTruthTauLV(3,1);
	  TruthDecayFromTau1=Ntp->GetTruthTauProductLV(2,13,0);
	  TruthDecayFromTau2=Ntp->GetTruthTauProductLV(3,211,1);decay=1;if(TruthDecayFromTau1==TruthDecayFromTau2)cout<<"2 3"<<endl;

	}
	if(Ntp->CheckDecayID(2,4)){
	  Tau1Truth=Ntp->GetTruthTauLV(2,0);
	  Tau2Truth=Ntp->GetTruthTauLV(4,1);
	  TruthDecayFromTau1=Ntp->GetTruthTauProductLV(2,13,0);
	  Pions2=Ntp->GetTruthPionsFromRho(1);
	  TruthDecayFromTau2=Pions2.at(0)+Pions2.at(1);decay=1;if(TruthDecayFromTau1==TruthDecayFromTau2)cout<<"2 4"<<endl;
	  if((sqrt((Tau1P4.Eta()-Tau1Truth.Eta())*(Tau1P4.Eta()-Tau1Truth.Eta())+(Tau1P4.Phi()-Tau1Truth.Phi())*(Tau1P4.Phi()-Tau1Truth.Phi())))<0.5 && (sqrt((Tau2P4.Eta()-Tau2Truth.Eta())*(Tau2P4.Eta()-Tau2Truth.Eta())+(Tau2P4.Phi()-Tau2Truth.Phi())*(Tau2P4.Phi()-Tau2Truth.Phi())))<0.5)
	    {
	      Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau2).E()-Pions2.at(1).E(),w);
	      Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau2).E()-Pions2.at(1).E())/Pions2.at(1).E(),w);
	    }
	  else{ Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau1).E()-Pions2.at(1).E(),w);Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau1).E()-Pions2.at(1).E())/Pions2.at(1).E(),w);}

	}
	if(Ntp->CheckDecayID(2,5)){
	  Tau1Truth=Ntp->GetTruthTauLV(2,0);
	  Tau2Truth=Ntp->GetTruthTauLV(5,1);
	  TruthDecayFromTau1=Ntp->GetTruthTauProductLV(2,13,0);
	  Pions2=Ntp->GetTruthPionsFromA1(1);
	  TruthDecayFromTau2=Pions2.at(0)+Pions2.at(1)+Pions2.at(2);decay=1;if(TruthDecayFromTau1==TruthDecayFromTau2)cout<<"2 5"<<endl;

	}


	if(Ntp->CheckDecayID(3,3)){
	  Tau1Truth=Ntp->GetTruthTauLV(3,0);
	  Tau2Truth=Ntp->GetTruthTauLV(3,1);
	  TruthDecayFromTau1=Ntp->GetTruthTauProductLV(3,211,0);
	  TruthDecayFromTau2=Ntp->GetTruthTauProductLV(3,211,1);decay=1;if(TruthDecayFromTau1==TruthDecayFromTau2)cout<<"3 3"<<endl;

	}
	if(Ntp->CheckDecayID(3,5)){
	  Tau1Truth=Ntp->GetTruthTauLV(3,0);
	  Tau2Truth=Ntp->GetTruthTauLV(5,1);
	  TruthDecayFromTau1=Ntp->GetTruthTauProductLV(3,211,0);
	  Pions2=Ntp->GetTruthPionsFromA1(1);
	  TruthDecayFromTau2=Pions2.at(0)+Pions2.at(1)+Pions2.at(2);decay=1;if(TruthDecayFromTau1==TruthDecayFromTau2)cout<<"3 5"<<endl;

	}


	if(Ntp->CheckDecayID(4,4)){
	  Tau1Truth=Ntp->GetTruthTauLV(4,0);
	  Tau2Truth=Ntp->GetTruthTauLV(4,1);
	  Pions1=Ntp->GetTruthPionsFromRho(0);
	  TruthDecayFromTau1=Pions1.at(0)+Pions1.at(1);
	  Pions2=Ntp->GetTruthPionsFromRho(1);
	  TruthDecayFromTau2=Pions2.at(0)+Pions2.at(1);decay=1;if(TruthDecayFromTau1==TruthDecayFromTau2)cout<<"4 4"<<endl;
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

	}
	if(Ntp->CheckDecayID(4,3)){
	  Tau1Truth=Ntp->GetTruthTauLV(4,0);
	  Tau2Truth=Ntp->GetTruthTauLV(3,1);
	  Pions1=Ntp->GetTruthPionsFromRho(0);
	  TruthDecayFromTau1=Pions1.at(0)+Pions1.at(1);
	  TruthDecayFromTau2=Ntp->GetTruthTauProductLV(3,211,1);decay=1;if(TruthDecayFromTau1==TruthDecayFromTau2)cout<<"4 3"<<endl;
	  if((sqrt((Tau1P4.Eta()-Tau1Truth.Eta())*(Tau1P4.Eta()-Tau1Truth.Eta())+(Tau1P4.Phi()-Tau1Truth.Phi())*(Tau1P4.Phi()-Tau1Truth.Phi())))<0.5 && (sqrt((Tau2P4.Eta()-Tau2Truth.Eta())*(Tau2P4.Eta()-Tau2Truth.Eta())+(Tau2P4.Phi()-Tau2Truth.Phi())*(Tau2P4.Phi()-Tau2Truth.Phi())))<0.5)
	    {
	      Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau1).E()-Pions1.at(1).E(),w);
	      Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau1).E()-Pions1.at(1).E())/Pions1.at(1).E(),w);
	    }
	  else{ Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau2).E()-Pions1.at(1).E(),w);Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau2).E()-Pions1.at(1).E())/Pions1.at(1).E(),w);}

	}
	if(Ntp->CheckDecayID(4,5)){
	  Tau1Truth=Ntp->GetTruthTauLV(4,0);
	  Tau2Truth=Ntp->GetTruthTauLV(5,1);
	  Pions1=Ntp->GetTruthPionsFromRho(0);
	  TruthDecayFromTau1=Pions1.at(0)+Pions1.at(1);
	  Pions2=Ntp->GetTruthPionsFromA1(1);
	  TruthDecayFromTau2=Pions2.at(0)+Pions2.at(1)+Pions2.at(2);decay=1;if(TruthDecayFromTau1==TruthDecayFromTau2){cout<<"4 5"<<endl;TruthDecayFromTau1.Print();TruthDecayFromTau2.Print();}
	  if((sqrt((Tau1P4.Eta()-Tau1Truth.Eta())*(Tau1P4.Eta()-Tau1Truth.Eta())+(Tau1P4.Phi()-Tau1Truth.Phi())*(Tau1P4.Phi()-Tau1Truth.Phi())))<0.5 && (sqrt((Tau2P4.Eta()-Tau2Truth.Eta())*(Tau2P4.Eta()-Tau2Truth.Eta())+(Tau2P4.Phi()-Tau2Truth.Phi())*(Tau2P4.Phi()-Tau2Truth.Phi())))<0.5)
	    {
	      Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau1).E()-Pions1.at(1).E(),w);
	      Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau1).E()-Pions1.at(1).E())/Pions1.at(1).E(),w);
	    }
	  else{ Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau2).E()-Pions1.at(1).E(),w);Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau2).E()-Pions1.at(1).E())/Pions1.at(1).E(),w);}

	}


	if(Ntp->CheckDecayID(5,5)){
	  Tau1Truth=Ntp->GetTruthTauLV(5,0);
	  Tau2Truth=Ntp->GetTruthTauLV(5,1);
	  Pions1=Ntp->GetTruthPionsFromA1(0);
	  TruthDecayFromTau1=Pions1.at(0)+Pions1.at(1)+Pions1.at(2);
	  Pions2=Ntp->GetTruthPionsFromA1(1);
	  TruthDecayFromTau2=Pions2.at(0)+Pions2.at(1)+Pions2.at(2);decay=1;if(TruthDecayFromTau1==TruthDecayFromTau2){cout<<"5 5"<<endl;TruthDecayFromTau1.Print();TruthDecayFromTau2.Print();}

	}
	if (decay==1)
	  {
	    if(Tau1Truth==Tau2Truth)cout<<"Same Taus"<<endl;

	    TLorentzVector TruthZ    = Tau1Truth+Tau2Truth;
	    // TLorentzVector TruthrhoFromTau2 = (Ntp->GetTruthTauProductLV(3,211) + Ntp->GetTruthTauProductLV(3,111));
	    if(TruthDecayFromTau1==TruthDecayFromTau2)cout<<"Same decay particle"<<endl;

	    //double visiblePtTruth = (TruthDecayFromTau1 + TruthDecayFromTau2).Pt();

	    TLorentzVector VisDecayfromTauplus;
	    TLorentzVector VisDecayfromTauminus;
	  
	    TauTauTruthMass.at(t).Fill((Tau1Truth+Tau2Truth).M(),w);
	  
	    if((sqrt((Tauplusvis.Eta()-Tau1Truth.Eta())*(Tauplusvis.Eta()-Tau1Truth.Eta())+(Tauplusvis.Phi()-Tau1Truth.Phi())*(Tauplusvis.Phi()-Tau1Truth.Phi())))<0.5 && (sqrt((Tauminusvis.Eta()-Tau2Truth.Eta())*(Tauminusvis.Eta()-Tau2Truth.Eta())+(Tauminusvis.Phi()-Tau2Truth.Phi())*(Tauminusvis.Phi()-Tau2Truth.Phi())))<0.5)
	      {
		Tauplustruth=Tau1Truth;
		Tauminustruth=Tau2Truth;
		VisDecayfromTauplus=TruthDecayFromTau1;
		VisDecayfromTauminus=TruthDecayFromTau2;
	      }
	    else
	      {
		Tauplustruth=Tau2Truth;
		Tauminustruth=Tau1Truth;
		VisDecayfromTauplus=TruthDecayFromTau2;
		VisDecayfromTauminus=TruthDecayFromTau1;
	      }
	    
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
	    if(a1a1 && id==30530533) {
	      ResolPullTauplusFroma1a1MeanEnergy.at(t).Fill((TauplusMean.E()-Tauplustruth.E())/Tauplustruth.E(),w);
	      ResolPullTauplusFroma1a1MZEnergy.at(t).Fill((TauplusPairConstraint.E()-Tauplustruth.E())/Tauplustruth.E(),w);
	      ResolPullTauplusFroma1a1MeanMomentum.at(t).Fill((TauplusMean.P()-Tauplustruth.P())/Tauplustruth.P(),w);
	      ResolPullTauplusFroma1a1MZMomentum.at(t).Fill((TauplusPairConstraint.P()-Tauplustruth.P())/Tauplustruth.P(),w);

	      ResolPullTauminusFroma1a1MeanEnergy.at(t).Fill((TauminusMean.E()-Tauminustruth.E())/Tauminustruth.E(),w);
	      ResolPullTauminusFroma1a1MZEnergy.at(t).Fill((TauminusPairConstraint.E()-Tauminustruth.E())/Tauminustruth.E(),w);
	      ResolPullTauminusFroma1a1MeanMomentum.at(t).Fill((TauminusMean.P()-Tauminustruth.P())/Tauminustruth.P(),w);
	      ResolPullTauminusFroma1a1MZMomentum.at(t).Fill((TauminusPairConstraint.P()-Tauminustruth.P())/Tauminustruth.P(),w);

	      ResolPullTauTauFroma1a1MeanEnergy.at(t).Fill((ZMean.E()-TruthZ.E())/TruthZ.E(),w);
	      ResolPullTauTauFroma1a1MZEnergy.at(t).Fill((ZPairConstraint.E()-TruthZ.E())/TruthZ.E(),w);
	      ResolPullTauTauFroma1a1MeanMomentum.at(t).Fill((ZMean.P()-TruthZ.P())/TruthZ.P(),w);
	      ResolPullTauTauFroma1a1MZMomentum.at(t).Fill((ZPairConstraint.P()-TruthZ.P())/TruthZ.P(),w);

	      Phitrutha1a1.at(t).Fill(Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	      // PhiSvFitResa1a1.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	      // PhiVisResa1a1.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	    }
	    if(a1X && (id==10130533 || id==10230533 || id==10330533|| id==10430533))
	      {
		if(a1minus)
		  {
		    ResolPullTauminusFroma1XMeanEnergy.at(t).Fill((TauminusMean.E()-Tauminustruth.E())/Tauminustruth.E(),w);
		    ResolPullTauminusFroma1XMeanMomentum.at(t).Fill((TauminusMean.P()-Tauminustruth.P())/Tauminustruth.P(),w);
		    ResolPullXplusFroma1XMeanEnergy.at(t).Fill((TauplusMean.E()-Tauplustruth.E())/Tauplustruth.E(),w);
		    ResolPullXplusFroma1XMeanMomentum.at(t).Fill((TauplusMean.P()-Tauplustruth.P())/Tauplustruth.P(),w);
		    
		  }
		else if (a1plus)
		  {
		    ResolPullTauplusFroma1XMeanEnergy.at(t).Fill((TauplusMean.E()-Tauplustruth.E())/Tauplustruth.E(),w);
		    ResolPullTauplusFroma1XMeanMomentum.at(t).Fill((TauplusMean.P()-Tauplustruth.P())/Tauplustruth.P(),w);
		    ResolPullXminusFroma1XMeanEnergy.at(t).Fill((TauminusMean.E()-Tauminustruth.E())/Tauminustruth.E(),w);
		    ResolPullXminusFroma1XMeanMomentum.at(t).Fill((TauminusMean.P()-Tauminustruth.P())/Tauminustruth.P(),w);
		  }
	      }
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
  }
}
//  This is a function if you want to do something after the event loop
void ZTauTau::Finish() {
  
  if(mode == RECONSTRUCT) {
    std::cout<<" Starting Finish!  " <<std::endl;
    std::cout<<"A  Data  "<< NQCD.at(0).GetBinContent(1) << std::endl;
    std::cout<<"B  Data  "<< NQCD.at(0).GetBinContent(2) << std::endl;
    std::cout<<"C  Data  "<< NQCD.at(0).GetBinContent(3) << std::endl;
    std::cout<<"D  Data  "<< NQCD.at(0).GetBinContent(4) << std::endl;
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
    for(unsigned i=0;i<CrossSectionandAcceptance.size();i++){
      QCD_Integral_B.push_back(NQCD.at(i).GetBinContent(2));
      QCD_Integral_C.push_back(NQCD.at(i).GetBinContent(3));
      QCD_Integral_D.push_back(NQCD.at(i).GetBinContent(4));
      if(CrossSectionandAcceptance.at(i)>0){
  	QCD_Integral_B.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
  	QCD_Integral_C.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
  	QCD_Integral_D.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
      }
    }
    for(unsigned i=0;i<CrossSectionandAcceptance.size();i++){
      if(HConfig.GetID(i) == DataMCType::Data){
  	QCD_Integral_B_Data_minus_MC  += QCD_Integral_B.at(i);
  	QCD_Integral_C_Data_minus_MC += QCD_Integral_C.at(i);
  	QCD_Integral_D_Data_minus_MC += QCD_Integral_D.at(i);
      }
      if(CrossSectionandAcceptance.at(i)>0){
  	QCD_IntegralMC_B  += QCD_Integral_B.at(i);
  	QCD_IntegralMC_C  += QCD_Integral_C.at(i);
  	QCD_IntegralMC_D  += QCD_Integral_D.at(i);
      }
    }
    double CDFactor = (QCD_Integral_C_Data_minus_MC  - QCD_IntegralMC_C )/ (QCD_Integral_D_Data_minus_MC - QCD_IntegralMC_D);
    double QCD_Signal = QCD_Integral_B_Data_minus_MC *CDFactor;

    std::cout << "Factor: " << CDFactor << std::endl;
    std::cout << "QCD_Signal: " << QCD_Signal << std::endl;
    std::cout << "QCD in B region "<<  QCD_Integral_B_Data_minus_MC <<std::endl;
    std::cout << "QCD_Integral_B_Data_minus_MC is: " << QCD_Integral_B_Data_minus_MC << std::endl;
    std::cout << "QCD_Integral_C_Data_minus_MC is: " << QCD_Integral_C_Data_minus_MC << std::endl;
    std::cout << "QCD_Integral_D_Data_minus_MC is: " << QCD_Integral_D_Data_minus_MC << std::endl;
    std::cout << "QCD_IntegralMC_B is: " << QCD_IntegralMC_B << std::endl;
    std::cout << "QCD_IntegralMC_C is: " << QCD_IntegralMC_C << std::endl;
    std::cout << "QCD_IntegralMC_D is: " << QCD_IntegralMC_D << std::endl;
    ScaleAllHistOfType(HConfig.GetType(DataMCType::QCD),QCD_Signal/Nminus0.at(0).at(HConfig.GetType(DataMCType::QCD)).Integral());
  }
  Selection::Finish();
}
