#include "ZMuTau.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "SVFitObject.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"
 
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




ZMuTau::ZMuTau(TString Name_, TString id_):
  Selection(Name_,id_),
  cMu_pt(20),
  cMu_eta(2.1),
  cTau_pt(20),
  cTau_eta(2.1),
  DataMC_Corr(true,true,false),
  DataMC_CorrLeptonIso(true,false,true)
  // tauTrgSF("vtight")
{
  // TString basedir = "";
  // basedir = (TString)std::getenv("workdir")+"/Code/CommonFiles/";
  // myScaleFactor->init_ScaleFactor(basedir+"LeptonEfficiencies/Muon/Run2016BtoH/Muon_IsoMu24_OR_TkIsoMu24_2016BtoH_eff.root");

  ChargeSumDummy = -999;
  selMuon_IsoDummy = 999.;
}

ZMuTau::~ZMuTau(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  ZMuTau::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==Trigger)             cut.at(Trigger)=1.;
    if(i==Id_and_Kin)          cut.at(Id_and_Kin)=1.;
    if(i==NPairsFound)         cut.at(NPairsFound)=1.;
    if(i==MuonIsolation)       cut.at(MuonIsolation)=0.15;  
    if(i==TauIsolation)        cut.at(TauIsolation)=1.;
    if(i==DiMuon_Veto)         cut.at(DiMuon_Veto)=0.;
    if(i==LeptonVeto)          cut.at(LeptonVeto)=0.;
    if(i==PairCharge)          cut.at(PairCharge)=1.;
    //    if(i==PairMass)            cut.at(PairMass)=115.;
    if(i==MTM)                 cut.at(MTM)=40.;
    //if(i==deltaR)              cut.at(deltaR)=3.3;
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
      hlabel="Has good particles";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Id_and_Kin_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Id_and_Kin_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==NPairsFound){
      title.at(i)="At least 1 pair with good Delta R";
      hlabel="At least 1 pair with good Delta R";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NPairsFound_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NPairsFound_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==MuonIsolation){
      title.at(i)="Muon Isolation";
      hlabel="Muon Isolation";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonIsolation_",htitle,14,0.,0.35,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonIsolation_",htitle,14,0.,0.35,hlabel,"Events"));
    }
    else if(i==TauIsolation){
      title.at(i)="Pairs with good tau isolation";
      hlabel="Pairs with good tau isolation";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauIsolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauIsolation_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
     else if(i==DiMuon_Veto){
      title.at(i)="Di-Muon Veto";
      hlabel="Di-Muon";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_DiMuonVeto_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_DiMuonVeto_",htitle,2,-0.5,1.5,hlabel,"Events"));
      }
    else if(i==LeptonVeto){
      title.at(i)="Lepton Veto";
      hlabel="Third Lepton";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_LeptonVeto_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_LeptonVeto_",htitle,2,-0.5,1.5,hlabel,"Events"));
      }
    else if(i==PairCharge){
      title.at(i)="Pair Charge";
      hlabel="is pair OS";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PairCharge_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PairCharge_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    /* else if(i==PairMass){
      title.at(i)="Pair Visible Mass";
      hlabel="M(tau-tau)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PairMass_",htitle,15,0,220,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PairMass_",htitle,15,0,220,hlabel,"Events"));
      }*/

     else if(i==MTM){
      title.at(i)="Missing Transverse Mass";
      hlabel="Missing Transverse Mass";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MTM_",htitle,10,0,100,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MTM_",htitle,10,0,100,hlabel,"Events"));
      }
  } 
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");

  TauPT=HConfig.GetTH1D(Name+"_TauPT","Transverse momentum of selected #tau candidate",25,0,80," P_{T}(#tau), GeV","Events");
  TauE=HConfig.GetTH1D(Name+"_TauE","Energy of selected #tau candidate",20,15,150," E(#tau), GeV","Events");
  TauMass=HConfig.GetTH1D(Name+"_TauMass","Mass of selected #tau candidate",20,0,2.," M(#tau), GeV","Events");
  TauPhi=HConfig.GetTH1D(Name+"_TauPhi","Phi of selected #tau candidate",10,-3.14,3.14," #phi(#tau)","Events");
  TauEta=HConfig.GetTH1D(Name+"_TauEta","Pseudorapidity tau",20,-2.7,2.7," #eta(#tau)","Events");
  Taudz=HConfig.GetTH1D(Name+"_Taudz","Taudz",10,-0.12,0.12,"Taudz","Events");

  MuonPT=HConfig.GetTH1D(Name+"_MuonPT","Transverse momentum of selected #mu candidate",25,0,100," P_{T}(#mu), GeV","Events");
  MuonE=HConfig.GetTH1D(Name+"_MuonE","Energy of selected #mu candidate",20,20,150," E(#mu), GeV","Events");
  MuonMass=HConfig.GetTH1D(Name+"_MuonMass","Mass of selected #mu candidate",10,0.084,0.133," M(#mu), GeV","Events");
  MuonPhi=HConfig.GetTH1D(Name+"_MuonPhi","Phi of selected #mu candidate",10,-3.14,3.14," #phi(#mu)","Events");
  MuonEta=HConfig.GetTH1D(Name+"_MuonEta","Pseudorapidity muon",10,-2.5,2.5," #eta(#mu)","Events");
  Muondz=HConfig.GetTH1D(Name+"_Muondz","Muondz",10,-0.15,0.15,"Muondz","Events");
  Muondxy=HConfig.GetTH1D(Name+"_Muondxy","Muondxy",10,-0.04,0.04,"Muondxy","Events");
  MuonIsol=HConfig.GetTH1D(Name+"_MuonIsol","MuonIsol",20,0,0.15,"iso(#mu)","Events");
  
  /*
  againstElectronVLooseMVA6=HConfig.GetTH1D(Name+"_againstElectronVLooseMVA6","againstElectronVLooseMVA6",2,-0.5,1.5,"againstElectronVLooseMVA6","Events");
  againstElectronLooseMVA6=HConfig.GetTH1D(Name+"_againstElectronLooseMVA6","againstElectronLooseMVA6",2,-0.5,1.5,"againstElectronLooseMVA6","Events");
  againstElectronMediumMVA6=HConfig.GetTH1D(Name+"_againstElectronMediumMVA6","againstElectronMediumMVA6",2,-0.5,1.5,"againstElectronMediumMVA6","Events");
  againstElectronTightMVA6=HConfig.GetTH1D(Name+"_againstElectronTightMVA6","againstElectronTightMVA6",2,-0.5,1.5,"againstElectronTightMVA6","Events");
  againstElectronVTightMVA6=HConfig.GetTH1D(Name+"_againstElectronVTightMVA6","againstElectronVTightMVA6",2,-0.5,1.5,"againstElectronVTightMVA6","Events");
  againstMuonLoose3=HConfig.GetTH1D(Name+"_againstMuonLoose3","againstMuonLoose3",2,-0.5,1.5,"againstMuonLoose3","Events");
  againstMuonTight3=HConfig.GetTH1D(Name+"_againstMuonTight3","againstMuonTight3",2,-0.5,1.5,"againstMuonTight3","Events");
  byCombinedIsolationDeltaBetaCorrRaw3Hits=HConfig.GetTH1D(Name+"_byCombinedIsolationDeltaBetaCorrRaw3Hits","byCombinedIsolationDeltaBetaCorrRaw3Hits",10,0,20,"byCombinedIsolationDeltaBetaCorrRaw3Hits","Events");
  

  DiMuonVeto=HConfig.GetTH1D(Name+"_DiMuonVeto","DiMuonVeto",2,-0.5,1.5,"DiMuonVeto","Events");
  ExtraLeptonVeto=HConfig.GetTH1D(Name+"_ExtraLeptonVeto","ExtraLeptonVeto",2,-0.5,1.5,"ExtraLeptonVeto","Events");
  */
  TauHPSDecayMode=HConfig.GetTH1D(Name+"_TauHPSDecayMode","Decay mode of the selected #tau candidate",11,-0.5,10.5," HPS Mode ","Events");
  QCDShape=HConfig.GetTH1D(Name+"_QCDShape","QCDShape",2,-0.5,1.5,"QCD Shape","");

  TauTauVisMass=HConfig.GetTH1D(Name+"_TauTauVisMass","Visible invariant mass of a tau pair",15,20,100," M(#tau#tau), GeV","Events");
  TauTauTruthMass=HConfig.GetTH1D(Name+"_TauTauTruthMass","Truth invariant mass of a tau pair",40,0,150," M(#tau#tau)_{truth}, GeV","Events");

  SVChi2=HConfig.GetTH1D(Name+"_SVChi2","SV  #chi^{2}",10,0,18,"#chi^{2}","Events");
  SVQuality=HConfig.GetTH1D(Name+"_SVQuality","Track mathicn #DeltaR",10,0,2,"#Sigma#Delta R","Events");
  SVQualityVsSignificance=HConfig.GetTH2D(Name+"_SVQualityVsSignificance","Track mathicn #DeltaR vs significance",25,0,3,31,-0.5,5.5,"","Events");

  NWJets=HConfig.GetTH1D(Name+"_NWJets","NWJets",4,0.5,4.5,"NWJets in ABCD","Events");
  NWJetsRelaxed=HConfig.GetTH1D(Name+"_NWJetsRelaxed","NWJetsRelaxed",2,0.5,2.5,"NWJetsRelaxed in Low and High MT","Events");
  NQCD=HConfig.GetTH1D(Name+"_NQCD","NQCD",4,0.5,4.5,"NQCD in ABCD","Events");
  TauTauVisMassAHigh=HConfig.GetTH1D(Name+"_TauTauVisMassAHigh","TauTauVisMass_AHigh",30,0,200,"#tau_{h}#mu Visible Mass in AHigh","Events");
  TauTauVisMassBHigh=HConfig.GetTH1D(Name+"_TauTauVisMassBHigh","TauTauVisMass_BHigh",30,0,200,"#tau_{h}#mu Visible Mass in BHigh","Events");
  TauTauVisMassBLow=HConfig.GetTH1D(Name+"_TauTauVisMassBLow","TauTauVisMass_BLow",30,0,200,"#tau_{h}#mu Visible Mass in BLow","Events");
  MTAHigh=HConfig.GetTH1D(Name+"_MTAHigh","MTAHigh",20,70,120,"Transverse Mass in AHigh","Events");
  MTBHigh=HConfig.GetTH1D(Name+"_MTBHigh","MTBHigh",20,70,120,"Transverse Mass in BHigh","Events");
  MTBLow=HConfig.GetTH1D(Name+"_MTBLow","MTBLow",20,0,40,"Transverse Mass in BLow","Events");

  METAHigh=HConfig.GetTH1D(Name+"_METAHigh","METAHigh",20,0,200,"Missing Transverse Energy in AHigh","Events");

  PVSVSignificance=HConfig.GetTH1D(Name+"_PVSVSignificance"," PV-SV significance for tau decay mode = 10",10,-0.5,5.5," pType","Events");
  dRTauTau=HConfig.GetTH1D(Name+"_dRTauTau","#Delta R",20,0.,4.," #Delta R","Events");

  MET=HConfig.GetTH1D(Name+"_MET","MET",20,0,200,"MET, GeV","Events");
  METphi=HConfig.GetTH1D(Name+"_METphi","METphi",10,-3.14,3.14,"METphi","Events");
  PUPPImet=HConfig.GetTH1D(Name+"_PUPPImet","PUPPImet",10,0,250,"PUPPImet, GeV","Events");
  PUPPImetphi=HConfig.GetTH1D(Name+"_PUPPImetphi","PUPPImetphi",10,-3.7,3.7,"PUPPImetphi","Events");
  TransverseMass=HConfig.GetTH1D(Name+"_TransverseMass","TransverseMass, GeV",20,0,40,"TransverseMass","Events");
  
  NPrimeVtx=HConfig.GetTH1D(Name+"_NPrimeVtx","NPrimeVtx",10,0,60,"N vtx","Events");
  NPU=HConfig.GetTH1D(Name+"_npu","npu",10,0,50,"N pu","Events");
  RHO=HConfig.GetTH1D(Name+"_rho","rho",10,0,35,"rho","Events");
  
  NbJets=HConfig.GetTH1D(Name+"_NbJets","NbJets",17,0,17,"Number of jets","Events");
  

  Etavis = HConfig.GetTH1D(Name+"_Etavis", "Etavis", 30, -10, 10, "Eta between Tau+ and Tau- in the new xy plan with visible particles");

  //Phivispipi = HConfig.GetTH1D(Name+"_Phivispipi","Phivis #pi#pi",30,-TMath::Pi(),TMath::Pi(),"Visible angle between Tau+ and an initial proton in the new xy plan for #pi#pi channel");
  //Phivispirho = HConfig.GetTH1D(Name+"_Phivispirho","Phivis #pi#rho",30,-TMath::Pi(),TMath::Pi(),"Visible angle between Tau+ and an initial proton in the new xy plan for #pi#rho channel");
   Phivislpi = HConfig.GetTH1D(Name+"_Phivislpi","Phivis l#pi",30,-TMath::Pi(),TMath::Pi(),"Visible angle between Tau+ and an initial proton in the new xy plan for #pi channel");
   Phivislrho = HConfig.GetTH1D(Name+"_Phivislrho","Phivis l#rho",30,-TMath::Pi(),TMath::Pi(),"Visible angle between Tau+ and an initial proton in the new xy plan for l#rho channel");
   //Phivispia1 = HConfig.GetTH1D(Name+"_Phivispia1","Phivis a1#pi",30,-TMath::Pi(),TMath::Pi(),"Visible angle between Tau+ and an initial proton in the new xy plan for a1#pi channel");
   //Phivisrhorho = HConfig.GetTH1D(Name+"_Phivisrhorho","Phivis #rho#rho",30,-TMath::Pi(),TMath::Pi(),"Visible angle between Tau+ and an initial proton in the new xy plan for #rho#rho channel");
   //Phivisrhoa1 = HConfig.GetTH1D(Name+"_Phivisrhoa1","Phivis a1#rho",30,-TMath::Pi(),TMath::Pi(),"Visible angle between Tau+ and an initial proton in the new xy plan for a1#rho channel");
   Phivisla1 = HConfig.GetTH1D(Name+"_Phivisla1","Phivis la1",30,-TMath::Pi(),TMath::Pi(),"Visible angle between Tau+ and an initial proton in the new xy plan for la1 channel");
   //Phivisa1a1 = HConfig.GetTH1D(Name+"_Phivisa1a1","Phivis a1a1",30,-TMath::Pi(),TMath::Pi(),"Visible angle between Tau+ and an initial proton in the new xy plan for a1a1 channel");

  Thetavis = HConfig.GetTH1D(Name+"_Thetavis","Thetavis",30,0.,TMath::Pi(),"Original theta of Tau- with visible particles");
  
  
  //PhiInf8GeVvis=HConfig.GetTH1D(Name+"_PhiInf8GeVvis","PhiInf8GeVvis",30,-TMath::Pi(),TMath::Pi(),"Visible angle between Tau+ and an initial proton in the new xy plan for ZPt<8GeV");
  //PhiSup8GeVvis=HConfig.GetTH1D(Name+"_PhiSup8GeVvis","PhiSup8GeVvis",30,-TMath::Pi(),TMath::Pi(),"Visible angle between Tau+ and an initial proton in the new xy plan for ZPt>8GeV");

  Etatruth = HConfig.GetTH1D(Name+"_Etatruth", "Etatruth", 30, -10, 10, "Real Eta between Tau+ and Tau- in the new xy plan");

  //Phitruthpipi = HConfig.GetTH1D(Name+"_Phitruthpipi","Phitruth #pi#pi",30,-TMath::Pi(),TMath::Pi(),"Real angle between Tau+ and an initial proton in the new xy plan for #pi#pi channel");
  //Phitruthpirho = HConfig.GetTH1D(Name+"_Phitruthpirho","Phitruth #pi#rho",30,-TMath::Pi(),TMath::Pi(),"Real angle between Tau+ and an initial proton in the new xy plan for #pi#rho channel");
  Phitruthlpi = HConfig.GetTH1D(Name+"_Phitruthlpi","Phitruth l#pi",30,-TMath::Pi(),TMath::Pi(),"Real angle between Tau+ and an initial proton in the new xy plan for #pi channel");
  Phitruthlrho = HConfig.GetTH1D(Name+"_Phitruthlrho","Phitruth l#rho",30,-TMath::Pi(),TMath::Pi(),"Real angle between Tau+ and an initial proton in the new xy plan for l#rho channel");
  //Phitruthpia1 = HConfig.GetTH1D(Name+"_Phitruthpia1","Phitruth a1#pi",30,-TMath::Pi(),TMath::Pi(),"Real angle between Tau+ and an initial proton in the new xy plan for a1#pi channel");
  //Phitruthrhorho = HConfig.GetTH1D(Name+"_Phitruthrhorho","Phitruth #rho#rho",30,-TMath::Pi(),TMath::Pi(),"Real angle between Tau+ and an initial proton in the new xy plan for #rho#rho channel");
  //Phitruthrhoa1 = HConfig.GetTH1D(Name+"_Phitruthrhoa1","Phitruth a1#rho",30,-TMath::Pi(),TMath::Pi(),"Real angle between Tau+ and an initial proton in the new xy plan for a1#rho channel");
   Phitruthla1 = HConfig.GetTH1D(Name+"_Phitruthla1","Phitruth la1",30,-TMath::Pi(),TMath::Pi(),"Real angle between Tau+ and an initial proton in the new xy plan for la1 channel");
   //Phitrutha1a1 = HConfig.GetTH1D(Name+"_Phitrutha1a1","Phitruth a1a1",30,-TMath::Pi(),TMath::Pi(),"Real angle between Tau+ and an initial proton in the new xy plan for a1a1 channel");

  Thetatruth = HConfig.GetTH1D(Name+"_Thetatruth","Thetatruth",30,0.,TMath::Pi(),"Real original theta of Tau-");

  Pi0EnergyRes=HConfig.GetTH1D(Name+"_Pi0EnergyRes","Energy resolution of Pi0",100,-50.,50.,"Energy resolution of Pi0, GeV","Events");
  Pi0EnergyResPull=HConfig.GetTH1D(Name+"_Pi0EnergyResPull","Energy Pull Plot of Pi0",100,-50.,50.,"Energy Pull Plot of Pi0, GeV","Events");

  
  NewPhivsDeltaPhi=HConfig.GetTH2D(Name+"_NewPhivsDeltaPhi","New #phi VS old #Delta#phi",30,-TMath::Pi(),TMath::Pi(),30,-TMath::Pi(),TMath::Pi(),"","Events");
  NewPhivsDeltaEta=HConfig.GetTH2D(Name+"_NewPhivsDeltaEta","New #phi VS old #Delta#eta",30,-TMath::Pi(),TMath::Pi(),40,-5.,5.,"","Events");
  //NewPhivsPhiproton=HConfig.GetTH2D(Name+"_NewPhivsPhiproton","New #phi VS old #phi of the proton",30,-TMath::Pi(),TMath::Pi(),30,-TMath::Pi(),TMath::Pi(),"","Events");
  NewPhivsPhiTauplus=HConfig.GetTH2D(Name+"_NewPhivsPhiTauplus","New #phi VS old #phi of the Tau+",30,-TMath::Pi(),TMath::Pi(),30,-TMath::Pi(),TMath::Pi(),"","Events");
  //NewPhivsEtaproton=HConfig.GetTH2D(Name+"_NewPhivsEtaproton","New #phi VS old #eta of the proton",30,-TMath::Pi(),TMath::Pi(),20,-2.7,2.7,"","Events");
  NewPhivsEtaTauplus=HConfig.GetTH2D(Name+"_NewPhivsEtaTauplus","New #phi VS old #eta of the Tau+",30,-TMath::Pi(),TMath::Pi(),20,-2.7,2.7,"","Events");
  NewPhivsZPt=HConfig.GetTH2D(Name+"_NewPhivsZPt","New #phi VS Pt_{Z}",30,-TMath::Pi(),TMath::Pi(),40,0,100,"","Events");

  NewPhiSignal=HConfig.GetTH1D(Name+"_NewPhiSignal","New #phi for all Signal",30,-TMath::Pi(),TMath::Pi(),"","Events");
  NewPhiQCD=HConfig.GetTH1D(Name+"_NewPhiQCD","New #phi for QCD",30,-TMath::Pi(),TMath::Pi(),"","Events");
																			      
  ZPtVis=HConfig.GetTH1D(Name+"_ZPtVis","Visible Pt_{Z}",20,0,40,"","Events");

  NewPhivsPhiTauminus=HConfig.GetTH2D(Name+"_NewPhivsPhiTauminus","New #phi VS old #phi of the Tau-",30,-TMath::Pi(),TMath::Pi(),30,-TMath::Pi(),TMath::Pi(),"","Events");
  NewPhivsEtaTauminus=HConfig.GetTH2D(Name+"_NewPhivsEtaTauminus","New #phi VS old #eta of the Tau-",30,-TMath::Pi(),TMath::Pi(),20,-2.7,2.7,"","Events");

  //NewPhivsTauplusPt=HConfig.GetTH2D(Name+"_NewPhivsTauplusPt","New #phi VS Pt_{#tau+}",30,-TMath::Pi(),TMath::Pi(),30,0,60,"","Events");
  //NewPhivsTauminusPt=HConfig.GetTH2D(Name+"_NewPhivsTauminusPt","New #phi VS Pt_{#tau-}",30,-TMath::Pi(),TMath::Pi(),30,0,60,"","Events");
  
  NewPhivsDeltaEtaplusminus=HConfig.GetTH2D(Name+"_NewPhivsDeltaEtaplusminus","New #phi VS old #Delta#eta for tau+ and tau-",30,-TMath::Pi(),TMath::Pi(),40,-5.,5.,"","Events");

  EtaplusvsEtaminus=HConfig.GetTH2D(Name+"_EtaplusvsEtaminus","#eta_{#tau+} VS #eta_{#tau-}",40,-2.7,2.7,40,-2.7,2.7,"","Events");

  EtaTau1vsEtaTau2=HConfig.GetTH2D(Name+"_EtaTau1vsEtaTau2","#eta_{#tau1} VS #eta_{#tau2}",40,-2.7,2.7,40,-2.7,2.7,"","Events");

  PhiTau1vsPhiTau2=HConfig.GetTH2D(Name+"_PhiTau1vsPhiTau2","#phi_{#tau1} VS #phi_{#tau2}",40,-TMath::Pi(),TMath::Pi(),40,-TMath::Pi(),TMath::Pi(),"","Events");

  DzTau1vsDzTau2=HConfig.GetTH2D(Name+"_DzTau1vsDzTau2","Dz_{#tau1} VS Dz_{#tau2}",20,-0.05,0.05,20,-0.05,0.05,"","Events");

  Selection::ConfigureHistograms();   //   do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);  // do not remove
}

void  ZMuTau::Store_ExtraDist(){

  //every new histo should be addedd to Extradist1d vector, just push it back;
  Extradist1d.push_back(&TauPT);
  Extradist1d.push_back(&TauE);
  Extradist1d.push_back(&TauMass);
  Extradist1d.push_back(&TauPhi);
  Extradist1d.push_back(&TauEta);
  Extradist1d.push_back(&Taudz);

  Extradist1d.push_back(&MuonPT);
  Extradist1d.push_back(&MuonE);
  Extradist1d.push_back(&MuonMass);
  Extradist1d.push_back(&MuonPhi);
  Extradist1d.push_back(&MuonEta);
  Extradist1d.push_back(&Muondz);
  Extradist1d.push_back(&Muondxy);
  Extradist1d.push_back(&MuonIsol);
  /*
  Extradist1d.push_back(&againstElectronVLooseMVA6);
  Extradist1d.push_back(&againstElectronLooseMVA6);
  Extradist1d.push_back(&againstElectronMediumMVA6);
  Extradist1d.push_back(&againstElectronTightMVA6);
  Extradist1d.push_back(&againstElectronVTightMVA6);
  Extradist1d.push_back(&againstMuonLoose3);
  Extradist1d.push_back(&againstMuonTight3);
  Extradist1d.push_back(&byCombinedIsolationDeltaBetaCorrRaw3Hits);
  //Extradist1d.push_back(&byIsolationMVA3oldDMwLTraw);
  //Extradist1d.push_back(&byIsolationMVA3newDMwLTraw);
  
  Extradist1d.push_back(&DiMuonVeto);
  Extradist1d.push_back(&ExtraLeptonVeto);
  */
  Extradist1d.push_back(&TauHPSDecayMode);

  Extradist1d.push_back(&dRTauTau);

  Extradist1d.push_back(&TauTauVisMass);
  Extradist1d.push_back(&TauTauTruthMass);

  Extradist1d.push_back(&SVChi2);
  Extradist1d.push_back(&SVQuality);
  Extradist2d.push_back(&SVQualityVsSignificance);
  Extradist1d.push_back(&QCDShape);
  Extradist1d.push_back(&NQCD);
  Extradist1d.push_back(&NWJets);
  Extradist1d.push_back(&NWJetsRelaxed);
  Extradist1d.push_back(&MTAHigh);
  Extradist1d.push_back(&MTBHigh);
  Extradist1d.push_back(&MTBLow);
  Extradist1d.push_back(&TauTauVisMassAHigh);
  Extradist1d.push_back(&TauTauVisMassBHigh);
  Extradist1d.push_back(&TauTauVisMassBLow);

  Extradist1d.push_back(&METAHigh);

  Extradist1d.push_back(&MET);
  Extradist1d.push_back(&METphi);
  Extradist1d.push_back(&PUPPImet);
  Extradist1d.push_back(&PUPPImetphi);
  Extradist1d.push_back(&TransverseMass);

  Extradist1d.push_back(&NPrimeVtx);
  Extradist1d.push_back(&NPU);
  Extradist1d.push_back(&RHO);

  Extradist1d.push_back(&PVSVSignificance);

  Extradist1d.push_back(&NbJets);

  
  Extradist1d.push_back(&Etavis);
  // Extradist1d.push_back(&Phivispipi);
  // Extradist1d.push_back(&Phivispirho);
  Extradist1d.push_back(&Phivislpi);
   Extradist1d.push_back(&Phivislrho);
   //Extradist1d.push_back(&Phivispia1);
   //Extradist1d.push_back(&Phivisrhorho);
   // Extradist1d.push_back(&Phivisrhoa1);
  Extradist1d.push_back(&Phivisla1);
  //Extradist1d.push_back(&Phivisa1a1);
  Extradist1d.push_back(&Thetavis);
  
  // Extradist1d.push_back(&PhiInf8GeVvis);
  //Extradist1d.push_back(&PhiSup8GeVvis);

  Extradist1d.push_back(&Etatruth);
  //Extradist1d.push_back(&Phitruthpipi);
  //Extradist1d.push_back(&Phitruthpirho);
  Extradist1d.push_back(&Phitruthlpi);
   Extradist1d.push_back(&Phitruthlrho);
   //Extradist1d.push_back(&Phitruthpia1);
   //Extradist1d.push_back(&Phitruthrhorho);
   //Extradist1d.push_back(&Phitruthrhoa1);
  Extradist1d.push_back(&Phitruthla1);
  //Extradist1d.push_back(&Phitrutha1a1);
  Extradist1d.push_back(&Thetatruth);
  
  Extradist1d.push_back(&Pi0EnergyRes);
  Extradist1d.push_back(&Pi0EnergyResPull);

  Extradist1d.push_back(&ZPtVis);

  Extradist2d.push_back(&NewPhivsDeltaPhi);
  Extradist2d.push_back(&NewPhivsDeltaEta);
  //Extradist2d.push_back(&NewPhivsPhiproton);
  Extradist2d.push_back(&NewPhivsPhiTauplus);
  //Extradist2d.push_back(&NewPhivsEtaproton);
  Extradist2d.push_back(&NewPhivsEtaTauplus);
  Extradist2d.push_back(&NewPhivsZPt);
  Extradist1d.push_back(&NewPhiSignal);
  Extradist1d.push_back(&NewPhiQCD);

  Extradist2d.push_back(&NewPhivsPhiTauminus);
  Extradist2d.push_back(&NewPhivsEtaTauminus);
  //Extradist2d.push_back(&NewPhivsTauminusPt);
  //Extradist2d.push_back(&NewPhivsTauplusPt);
  Extradist2d.push_back(&NewPhivsDeltaEtaplusminus);

  Extradist2d.push_back(&EtaplusvsEtaminus);
  Extradist2d.push_back(&EtaTau1vsEtaTau2);

  Extradist2d.push_back(&PhiTau1vsPhiTau2);
  Extradist2d.push_back(&DzTau1vsDzTau2);
}

void  ZMuTau::doEvent(){ //  Method called on every event
  unsigned int t;                // sample type, you may manage in your further analysis, if needed
  int id(Ntp->GetMCID());  //read event ID of a sample
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}  //  gives a warining if list of samples in Histo.txt  and SkimSummary.log do not coincide 
  //  std::cout<<"------------------ New Event -----------------------"<<std::endl;
  Charge = ChargeSumDummy;
  bool trig=0;
  bool match=0;
  value.at(Trigger)=0;
  std::vector<int> MuonIndex;
  std::vector<int> TriggerIndexVector ;
  std::vector<TString>  MatchedTriggerNames;
  MatchedTriggerNames.push_back("HLT_IsoMu24_v");
  TriggerIndexVector=Ntp->GetVectorTriggers(MatchedTriggerNames);
  
  for(unsigned int itrig = 0; itrig < TriggerIndexVector.size(); itrig++){
    if(Ntp->TriggerAccept(TriggerIndexVector.at(itrig))){
      trig=1;
    }
  }
  
  for(unsigned int iDaughter=0;   iDaughter  <  Ntp->NDaughters() ;iDaughter++ ) {
    if(Ntp->CHECK_BIT(Ntp->Daughters_trgMatched(iDaughter),5))
      {
	MuonIndex.push_back(iDaughter);
	match=1;
     }
  }
  if(trig && match)value.at(Trigger)=1;
  pass.at(Trigger)=(value.at(Trigger)==cut.at(Trigger)); 
  bool goodTau=false;
  bool goodMuon=false;
  std::vector<int> goodMuonsIndex;
  std::vector<int> DiMuon;
  std::vector<int> goodTausIndex;
  value.at(Id_and_Kin)=0;
  for(unsigned int iDaughter=0;   iDaughter  <  Ntp->NDaughters() ;iDaughter++ ) {  // loop over all daughters in the event
      if(Ntp->tauBaselineSelection(iDaughter,30., 2.3, 0,1)){
	goodTausIndex.push_back(iDaughter);
	goodTau=true;
      } 
  }
  for(unsigned int index=0;index<MuonIndex.size();index++)
    {
      if(Ntp->isMediumGoodMuon(MuonIndex.at(index))){
	if(Ntp->muonBaselineSelection(MuonIndex.at(index), 25.,2.1,2)){
	  goodMuonsIndex.push_back(MuonIndex.at(index));
	  goodMuon=true;
	}
      }
    }
  if(goodTau && goodMuon)value.at(Id_and_Kin)=1;
  pass.at(Id_and_Kin)=(value.at(Id_and_Kin) == cut.at(Id_and_Kin));
  unsigned int Muon= 999;
  unsigned int Tau= 999;
  std::vector<int>  Sorted;
  std::vector<int>  PairsIndex;
  std::vector<int>  PairsIndexTau;
  std::vector<int>  PairsIndexMuon;
  int i=0;
  value.at(NPairsFound)=0;
  for(unsigned int ipair =0; ipair <goodTausIndex.size(); ipair++)
    {
      for(unsigned int jpair =0; jpair <goodMuonsIndex.size(); jpair++)
	{
	  if((Ntp->Daughters_P4(goodTausIndex.at(ipair)).DeltaR(Ntp->Daughters_P4(goodMuonsIndex.at(jpair))))>0.5){
	    PairsIndex.push_back(i);
	    PairsIndexTau.push_back(goodTausIndex.at(ipair));
	    PairsIndexMuon.push_back(goodMuonsIndex.at(jpair));
	    value.at(NPairsFound)=1;
	    i++;
	  }
	}
    }
  TLorentzVector MuonP4;
  TLorentzVector TauP4;
  pass.at(NPairsFound)=(value.at(NPairsFound)==cut.at(NPairsFound));
  if(pass.at(NPairsFound)) {
    Sorted = Ntp->SortPair(PairsIndex,PairsIndexTau,PairsIndexMuon);
    Muon=PairsIndexMuon.at(Sorted.back());
    Tau=PairsIndexTau.at(Sorted.back());
    MuonP4 = Ntp->Daughters_P4(Muon);
    TauP4 = Ntp->TauP4_Corrected(Tau);
    value.at(TauIsolation) = Ntp->isIsolatedTau(Tau,"Tight");
    value.at(MuonIsolation) = Ntp->combreliso(Muon);
    value.at(DiMuon_Veto)=(Ntp->DiMuonVeto());
    std::vector<int> thirdLeptonCounter;
    for(unsigned int iDaughter=0;   iDaughter  <  Ntp->NDaughters() ;iDaughter++ ) {  // loop over all daughters in the event
      if((iDaughter!=Tau)&&(iDaughter!=Muon)){
	if(Ntp->ElectronVeto(iDaughter) || Ntp->MuonVeto(iDaughter))thirdLeptonCounter.push_back(iDaughter);
      }
    }
    value.at(LeptonVeto) = thirdLeptonCounter.size()>0;
    value.at(PairCharge)=0;
    bool isOS=false;
    isOS=((Ntp->Daughters_charge(Muon)/abs(Ntp->Daughters_charge(Muon))) != (Ntp->Daughters_charge(Tau)/abs(Ntp->Daughters_charge(Tau))));
    if(isOS)value.at(PairCharge) = 1;

    //value.at(PairMass) = 999.;
    value.at(MTM) = 999.;
    value.at(MTM)=Ntp->transverseMass(MuonP4.Pt(), MuonP4.Phi(), Ntp->MET(), Ntp->METphi());
    //value.at(PairMass)=(MuonP4+TauP4).M();
    //pass.at(PairMass) = (value.at(PairMass) <= cut.at(PairMass));
  } 
  pass.at(MTM) = (value.at(MTM) <= cut.at(MTM));
  pass.at(MuonIsolation) = (value.at(MuonIsolation) < cut.at(MuonIsolation));
  //pass.at(againstEleMu)=(value.at(againstEleMu)==cut.at(againstEleMu));
  pass.at(TauIsolation) = (value.at(TauIsolation) == cut.at(TauIsolation));
  pass.at(DiMuon_Veto) = ( value.at(DiMuon_Veto) ==  cut.at(DiMuon_Veto));
  pass.at(LeptonVeto) = ( value.at(LeptonVeto) ==  cut.at(LeptonVeto));
  pass.at(PairCharge) = (value.at(PairCharge) == cut.at(PairCharge));

  // Here you can defined different type of weights you want to apply to events.
  double wobs=1;
  double w=1;
  double MuonSF(1);
  double MuonIso(1);


  if(!Ntp->isData() && id!=DataMCType::QCD && id!=DataMCType::W_lnu){
    //    w *= reweight.weight(2016,26,Ntp->PUNumInteractions());
    w *= reweight.PUweightHTT(Ntp->npu());
      //std::cout<<" pu weigh HTT  "<< reweight.PUweightHTT(Ntp->npu())<<std::endl;
    
    if(!Ntp->isData() && pass.at(NPairsFound) && (id==33|| id == 10110333 || id == 10110433|| id == 10130533|| id ==10210333|| id == 10210433|| id == 10230533|| id ==10310333 || id ==10330533 || id ==10410433 || id == 10410333|| id == 10430533|| id == 30530533)){
      w *= 0.95;  // Tau ID  correction
       MuonSF = DataMC_Corr.get_ScaleFactor(Ntp->Daughters_P4(Muon).Pt(),Ntp->Daughters_P4(Muon).Eta());
       MuonIso =  DataMC_CorrLeptonIso.get_ScaleFactor(Ntp->Daughters_P4(Muon).Pt(),Ntp->Daughters_P4(Muon).Eta());
      //      std::cout<<"MuonSF   "<<MuonSF <<"muon iso   "<<MuonIso <<std::endl;
    }
    w*=MuonSF;
    w*=MuonIso;

    // std::cout<<"weight "<< MuonSF*MuonIso*0.95 <<std::endl;

  TLorentzVector genMomentum(0,0,0,0);
  if( id == 33|| id == 10110333 || id == 10110433|| id == 10130533|| id ==10210333|| id == 10210433|| id == 10230533|| id ==10310333 || id ==10330533 || id ==10410433 || id == 10410333|| id == 10430533|| id == 30530533){
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

  if(!Ntp->isData() && id!=DataMCType::QCD)w*=Ntp->MC_weight(); //generator weight


  w*=Ntp->stitch_weight();

  //  std::cout<<"zpt "<< zptw<<std::endl;

  }


 // QCD ABCD BG Method
  /*******************
   *        |        *
   *    C   |    D   *  MuAntiIso
   *        |        *       
   * ----------------*------ 
   *        |        *       
   *    A   |    B   *  MuIso   
   *        |        *
   *******************
   *   OS   |    SS
   *
   *
   */


// W+Jets ABCD BG Method
  /*******************
   *        |        *
   *    C   |    D   *  High MT
   *        |        *       
   * ----------------*------ 
   *        |        *       
   *    A   |    B   *  Low MT   
   *        |        *
   *******************
   *   OS   |    SS
   *
   *
   */
  std::vector<unsigned int> exclude_cuts;
  std::vector<unsigned int> relaxed_cuts;
  std::vector<unsigned int> WJetsexclude_cuts;
  exclude_cuts.push_back(PairCharge);
  exclude_cuts.push_back(MuonIsolation);
  relaxed_cuts.push_back(TauIsolation);
  relaxed_cuts.push_back(MuonIsolation);
  relaxed_cuts.push_back(MTM);
  WJetsexclude_cuts.push_back(PairCharge);
  WJetsexclude_cuts.push_back(MTM);
  
  // std::cout<<" before  " << pass.at(TriggerOk) << "    " <<   pass.at(PrimeVtx) << "    " <<  pass.at(nGoodPairs)<< "    " <<   pass.at(FirstTauIsolation) << "    " <<  pass.at(SecondTauIsolation) << "    " <<  pass.at(nGoodMuons) << "    " <<  pass.at(PairCharge) << "  passAllBut  " << passAllBut(exclude_cuts) <<std::endl;
  if(passAllBut(relaxed_cuts)){
    MuonP4 = Ntp->Daughters_P4(Muon);
    TauP4 = Ntp->TauP4_Corrected(Tau);
    if(Ntp->combreliso(Muon)<0.3 && Ntp->isMediumGoodTau(Tau) && Ntp->transverseMass(MuonP4.Pt(), MuonP4.Phi(), Ntp->MET(), Ntp->METphi())<40.){
      NWJetsRelaxed.at(t).Fill(1.,w);
    }
    else if(Ntp->combreliso(Muon)<0.3 && Ntp->isMediumGoodTau(Tau) && Ntp->transverseMass(MuonP4.Pt(), MuonP4.Phi(), Ntp->MET(), Ntp->METphi())>70.){
      NWJetsRelaxed.at(t).Fill(2.,w);
    }
  }
  if(passAllBut(WJetsexclude_cuts)) {
    MuonP4 = Ntp->Daughters_P4(Muon);
    TauP4 = Ntp->TauP4_Corrected(Tau);
    if(pass.at(PairCharge)) {
      if(Ntp->transverseMass(MuonP4.Pt(), MuonP4.Phi(), Ntp->MET(), Ntp->METphi())<40.) {
	NWJets.at(t).Fill(1.,w); //A Low
      }
      if(Ntp->transverseMass(MuonP4.Pt(), MuonP4.Phi(), Ntp->MET(), Ntp->METphi())>70.){
	NWJets.at(t).Fill(3.,w); //A High
	//MTAHigh.at(t).Fill(Ntp->transverseMass(MuonP4.Pt(), MuonP4.Phi(), Ntp->MET(), Ntp->METphi()),w);
	//TauTauVisMassAHigh.at(t).Fill((MuonP4+TauP4).M(),w);
      }
    }
    if(!pass.at(PairCharge)){
      if(Ntp->transverseMass(MuonP4.Pt(), MuonP4.Phi(), Ntp->MET(), Ntp->METphi())<40.){
	NWJets.at(t).Fill(2.,w); //B Low
	//MTBLow.at(t).Fill(Ntp->transverseMass(MuonP4.Pt(), MuonP4.Phi(), Ntp->MET(), Ntp->METphi()),w);
	//TauTauVisMassBLow.at(t).Fill((MuonP4+TauP4).M(),w);
      }
      if(Ntp->transverseMass(MuonP4.Pt(), MuonP4.Phi(), Ntp->MET(), Ntp->METphi())>70.){
	NWJets.at(t).Fill(4.,w); //B High
	//MTBHigh.at(t).Fill(Ntp->transverseMass(MuonP4.Pt(), MuonP4.Phi(), Ntp->MET(), Ntp->METphi()),w);
	//TauTauVisMassBHigh.at(t).Fill((MuonP4+TauP4).M(),w);
      }
    }
  }
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
      if(pass.at(MuonIsolation) ){
	NQCD.at(t).Fill(1.,w); //A
      }
      if(!pass.at(MuonIsolation)){
	NQCD.at(t).Fill(3.,w); //C
      }
    }
    
    if(!pass.at(PairCharge)){
      if(pass.at(MuonIsolation)){
	NQCD.at(t).Fill(2.,w); //B
      }
      if(!pass.at(MuonIsolation)){
	NQCD.at(t).Fill(4.,w); //D
      }
    }
  }

  bool IsQCDEvent = false;
  //  if(passAllBut(exclude_cuts)){
  if(!pass.at(PairCharge)){
    if(id == DataMCType::Data){
      QCDShape.at(t).Fill(1,w);
      t=HConfig.GetType(DataMCType::QCD);
      IsQCDEvent = true;
    }
  }
  //  }
  if(IsQCDEvent){    pass.at(PairCharge)= true;}  


  if(passAllBut(WJetsexclude_cuts)) {
    if(pass.at(PairCharge)) {
      if(Ntp->transverseMass(MuonP4.Pt(), MuonP4.Phi(), Ntp->MET(), Ntp->METphi())>70.){ //A High
	MTAHigh.at(t).Fill(Ntp->transverseMass(MuonP4.Pt(), MuonP4.Phi(), Ntp->MET(), Ntp->METphi()),w);
	//cout<<"w:   "<<Ntp->transverseMass(MuonP4.Pt(), MuonP4.Phi(), Ntp->MET(), Ntp->METphi())<<endl;
	TauTauVisMassAHigh.at(t).Fill((MuonP4+TauP4).M(),w);
	  METAHigh.at(t).Fill(Ntp->MET(),w);
    	}
      }
      if(!pass.at(PairCharge)){
    	if(Ntp->transverseMass(MuonP4.Pt(), MuonP4.Phi(), Ntp->MET(), Ntp->METphi())<40.){ //B Low
    	  MTBLow.at(t).Fill(Ntp->transverseMass(MuonP4.Pt(), MuonP4.Phi(), Ntp->MET(), Ntp->METphi()),w);
    	  TauTauVisMassBLow.at(t).Fill((MuonP4+TauP4).M(),w);
    	}
    	if(Ntp->transverseMass(MuonP4.Pt(), MuonP4.Phi(), Ntp->MET(), Ntp->METphi())>70.){ //B High
    	  MTBHigh.at(t).Fill(Ntp->transverseMass(MuonP4.Pt(), MuonP4.Phi(), Ntp->MET(), Ntp->METphi()),w);
    	  TauTauVisMassBHigh.at(t).Fill((MuonP4+TauP4).M(),w);
    	}
      }
    }





  
     

  bool status=AnalysisCuts(t,w,wobs);  // boolean that say whether your event passed critera defined in pass vector. The whole vector must be true for status = true
  ///////////////////////////////////////////////////////////
  // Analyse events  which passed selection
  if(status){
    Sorted = Ntp->SortPair(PairsIndex,PairsIndexTau,PairsIndexMuon);
    Muon=PairsIndexMuon.at(Sorted.back());
    Tau=PairsIndexTau.at(Sorted.back());
    MuonP4 = Ntp->Daughters_P4(Muon);
    TauP4 = Ntp->TauP4_Corrected(Tau);
    double pvx(0);
    pvx =  Ntp->npv();

    cout<<"1   Eta Tau= "<<TauP4.Eta()<<"  Eta Muon= "<<MuonP4.Eta()<<endl;

    // if(id == DataMCType::Data) pvx =  Ntp->npv();
     if(id !=DataMCType::Data && id !=DataMCType::QCD && id !=DataMCType::W_lnu)	  pvx = Ntp->PUNumInteractions();
     if(Ntp->PFTau_hassecondaryVertex(Tau) && Ntp->isPVCovAvailable()){
       if(Ntp->PFTau_secondaryVertex_TracksMatchingQuality(Tau) < 0.01)   PVSVSignificance.at(t).Fill( Ntp->PFTau_FlightLength_significance(Ntp->PVtx(),Ntp->PFTau_TIP_primaryVertex_cov(), Ntp->PFTau_secondaryVertex_pos(Tau), Ntp->PFTau_TIP_secondaryVertex_cov(Tau)),w);
       SVChi2.at(t).Fill(Ntp->PFTau_secondaryVertex_vtxchi2(Tau),w);
       SVQuality.at(t).Fill(Ntp->PFTau_secondaryVertex_TracksMatchingQuality(Tau),w);
SVQualityVsSignificance.at(t).Fill(Ntp->PFTau_secondaryVertex_TracksMatchingQuality(Tau),Ntp->PFTau_FlightLength_significance(Ntp->PVtx(),Ntp->PFTau_TIP_primaryVertex_cov(), Ntp->PFTau_secondaryVertex_pos(Tau), Ntp->PFTau_TIP_secondaryVertex_cov(Tau)));
     }
  NPrimeVtx.at(t).Fill(pvx,w);
  NPU.at(t).Fill(Ntp->npu(),w);
  RHO.at(t).Fill(Ntp->rho(),w);
  
  std::vector<int> thirdLepton;

  TauPT.at(t).Fill(TauP4.Pt(),w);
  TauE.at(t).Fill(TauP4.E(),w);
  TauMass.at(t).Fill(TauP4.M(),w);
  TauPhi.at(t).Fill(TauP4.Phi(),w);
  TauEta.at(t).Fill(TauP4.Eta(),w);
  Taudz.at(t).Fill(Ntp->dz(Tau),w);
  

  MuonPT.at(t).Fill(MuonP4.Pt(),w);
  MuonE.at(t).Fill(MuonP4.E(),w);
  MuonMass.at(t).Fill(MuonP4.M(),w);
  MuonPhi.at(t).Fill(MuonP4.Phi(),w);
  MuonEta.at(t).Fill(MuonP4.Eta(),w);
  Muondz.at(t).Fill(Ntp->dz(Muon),w);
  Muondxy.at(t).Fill(Ntp->dxy(Muon),w);
  MuonIsol.at(t).Fill(Ntp->combreliso(Muon),w);

  /*
  againstElectronVLooseMVA6.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau),Ntp->Bit_againstElectronVLooseMVA6),w);
  againstElectronLooseMVA6.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau),Ntp->Bit_againstElectronLooseMVA6),w);
  againstElectronMediumMVA6.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau),Ntp->Bit_againstElectronMediumMVA6),w);
  againstElectronTightMVA6.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau),Ntp->Bit_againstElectronTightMVA6),w);
  againstElectronVTightMVA6.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau),Ntp->Bit_againstElectronVTightMVA6),w);
  againstMuonLoose3.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau),Ntp->Bit_againstMuonLoose3),w);
  againstMuonTight3.at(t).Fill(Ntp->CHECK_BIT(Ntp->tauID(Tau),Ntp->Bit_againstMuonTight3),w);
  byCombinedIsolationDeltaBetaCorrRaw3Hits.at(t).Fill(Ntp->Daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits(Tau),w);
  
  DiMuonVeto.at(t).Fill(Ntp->DiMuonVeto(),w);
  for(unsigned int iDaughter=0;   iDaughter  <  Ntp->NDaughters() ;iDaughter++ ) {
	if((iDaughter!=Tau)&&(iDaughter!=Muon)){
	  if(Ntp->ElectronVeto(iDaughter) || Ntp->MuonVeto(iDaughter))thirdLepton.push_back(iDaughter);
	}
      }
      if(thirdLepton.size()>0)ExtraLeptonVeto.at(t).Fill(1.,w);
  else ExtraLeptonVeto.at(t).Fill(0.,w);*/
  TauHPSDecayMode.at(t).Fill(Ntp->decayMode(Tau),w);
  
  TauTauVisMass.at(t).Fill((MuonP4+TauP4).M(),w);
  dRTauTau.at(t).Fill(MuonP4.DeltaR(TauP4),w);

  MET.at(t).Fill(Ntp->MET(),w);
  METphi.at(t).Fill(Ntp->METphi(),w);
  PUPPImet.at(t).Fill(Ntp->PUPPImet(),w);
  PUPPImetphi.at(t).Fill(Ntp->PUPPImetphi(),w);
  TransverseMass.at(t).Fill(Ntp->transverseMass(MuonP4.Pt(), MuonP4.Phi(), Ntp->MET(), Ntp->METphi()),w);
  
  int jets_counter=0;
  for(int ijet=0; ijet< Ntp->JetsNumber(); ijet++) {
    if((((Ntp->Jet_P4(ijet)).Pt())>20) && (fabs((Ntp->Jet_P4(ijet).Eta())<4.7))) {
      if((((Ntp->Jet_P4(ijet)).DeltaR(Ntp->Daughters_P4(Muon)))>0.5) && (((Ntp->Jet_P4(ijet)).DeltaR(Ntp->Daughters_P4(Tau)))>0.5))jets_counter++; {
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

  if(Ntp->Daughters_charge(Tau)>0)
    {
      Tauplusvis=TauP4;
      Tauminusvis=MuonP4;
    }
  else
    {
      Tauplusvis=MuonP4;
      Tauminusvis=TauP4;
    }
  TVector3 Tauplus3Dvis(Tauplusvis.X(),Tauplusvis.Y(),Tauplusvis.Z());
  TVector3 Tauminus3Dvis(Tauminusvis.X(),Tauminusvis.Y(),Tauminusvis.Z());
  TVector3 Zvis=(Tauminus3Dvis).Unit();
  TVector3 Yvis=(Tauminus3Dvis.Cross(Tauplus3Dvis)).Unit();
  TVector3 Xvis=(Yvis.Cross(Tauminus3Dvis)).Unit();
  TVector3 protonvis(0,0,1);
  TVector3 xyprotonvis(protonvis.Dot(Xvis),protonvis.Dot(Yvis),0);
  TVector3 xytauplusvis(Tauplus3Dvis.Dot(Xvis),Tauplus3Dvis.Dot(Yvis),0);

  Etavis.at(t).Fill(TMath::ATanH((Tauplus3Dvis*Zvis)/(sqrt((Tauplus3Dvis*Zvis)*(Tauplus3Dvis*Zvis)+(Tauplus3Dvis*Yvis)*(Tauplus3Dvis*Yvis)+(Tauplus3Dvis*Xvis)*(Tauplus3Dvis*Xvis)))),w);

  // if(id==10310333){	    Phivispipi.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),w);}
  // if(id==10410333){	    Phivispirho.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),w);}
  // if(id==10410433){	    Phivisrhorho.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),w);}
  // if(id==10330533){	    Phivispia1.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),w);}
  // if(id==10430533){	    Phivisrhoa1.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),w);}
  // if(id==30530533){	    Phivisa1a1.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),w);}
  if(id==10110333 || id==10210333){	    Phivislpi.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),w);}
  if(id==10110433 || id==10210433){	    Phivislrho.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),w);}
  if(id==10130533 || id==10230533){	    Phivisla1.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),w);}
  
  Thetavis.at(t).Fill(Tauminusvis.Theta(),w);
  
  
  ZPtVis.at(t).Fill((TauP4+MuonP4).Pt(),w);

  if(id==33 || id == 10110333 || id == 10110433|| id == 10130533|| id ==10210333|| id == 10210433|| id == 10230533|| id ==10310333 || id ==10330533 || id ==10410433 || id == 10410333 || id == 10430533 || id == 30530533){
    NewPhivsDeltaPhi.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),Ntp->DeltaPhi(protonvis.Phi(),Tauplus3Dvis.Phi()),w);
    NewPhivsDeltaEta.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),protonvis.Eta()-Tauplus3Dvis.Eta(),w);
    //NewPhivsPhiproton.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),protonvis.Phi(),w);
    NewPhivsPhiTauplus.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),Tauplus3Dvis.Phi(),w);
    //NewPhivsEtaproton.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),protonvis.Eta(),w);
    NewPhivsEtaTauplus.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),Tauplus3Dvis.Eta(),w);
    NewPhivsZPt.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),(TauP4+MuonP4).Pt(),w);
    
    NewPhiSignal.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),w);


    NewPhivsPhiTauminus.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),Tauminus3Dvis.Phi(),w);
    NewPhivsEtaTauminus.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),Tauminus3Dvis.Eta(),w);
    //NewPhivsTauminusPt.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),Tauminusvis.Pt(),w);
    //NewPhivsTauplusPt.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),Tauplusvis.Pt(),w);

    NewPhivsDeltaEtaplusminus.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),Tauminus3Dvis.Eta()-Tauplus3Dvis.Eta(),w);
    
    EtaplusvsEtaminus.at(t).Fill(Tauplus3Dvis.Eta(),Tauminus3Dvis.Eta(),w);
    EtaTau1vsEtaTau2.at(t).Fill(MuonP4.Eta(),TauP4.Eta(),w);
    cout<<"2   Eta Tau= "<<TauP4.Eta()<<"  Eta Muon= "<<MuonP4.Eta()<<endl;
    PhiTau1vsPhiTau2.at(t).Fill(MuonP4.Phi(),TauP4.Phi(),w);

    DzTau1vsDzTau2.at(t).Fill(Ntp->dz(Muon),Ntp->dz(Tau),w);
    
    //if((TauP4+MuonP4).Pt()<=8.)PhiInf8GeVvis.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()));
    //if((TauP4+MuonP4).Pt()>8.)PhiSup8GeVvis.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()));
      
  }
  if(id==DataMCType::QCD || id==DataMCType::Data)NewPhiQCD.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi()),w);							
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
	if((sqrt((TauP4.Eta()-Tau1Truth.Eta())*(TauP4.Eta()-Tau1Truth.Eta())+(TauP4.Phi()-Tau1Truth.Phi())*(TauP4.Phi()-Tau1Truth.Phi())))<0.5 && (sqrt((MuonP4.Eta()-Tau2Truth.Eta())*(MuonP4.Eta()-Tau2Truth.Eta())+(MuonP4.Phi()-Tau2Truth.Phi())*(MuonP4.Phi()-Tau2Truth.Phi())))<0.5)
	  {
	    Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Muon).E()-Pions2.at(1).E(),w);
	    Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Muon).E()-Pions2.at(1).E())/Pions2.at(1).E(),w);
	  }
	else{ Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau).E()-Pions2.at(1).E(),w);Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau).E()-Pions2.at(1).E())/Pions2.at(1).E(),w);}
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
	if((sqrt((TauP4.Eta()-Tau1Truth.Eta())*(TauP4.Eta()-Tau1Truth.Eta())+(TauP4.Phi()-Tau1Truth.Phi())*(TauP4.Phi()-Tau1Truth.Phi())))<0.5 && (sqrt((MuonP4.Eta()-Tau2Truth.Eta())*(MuonP4.Eta()-Tau2Truth.Eta())+(MuonP4.Phi()-Tau2Truth.Phi())*(MuonP4.Phi()-Tau2Truth.Phi())))<0.5)
	  {
	    Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Muon).E()-Pions2.at(1).E(),w);
	    Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Muon).E()-Pions2.at(1).E())/Pions2.at(1).E(),w);
	  }
	else{Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau).E()-Pions2.at(1).E(),w);
	  Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau).E()-Pions2.at(1).E())/Pions2.at(1).E(),w);}
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
	if((sqrt((TauP4.Eta()-Pions1.at(1).Eta())*(TauP4.Eta()-Pions1.at(1).Eta())+(TauP4.Phi()-Pions1.at(1).Phi())*(TauP4.Phi()-Pions1.at(1).Phi())))<0.5 && (sqrt((MuonP4.Eta()-Pions2.at(1).Eta())*(MuonP4.Eta()-Pions2.at(1).Eta())+(MuonP4.Phi()-Pions2.at(1).Phi())*(MuonP4.Phi()-Pions2.at(1).Phi())))<0.5)
	  {
	    Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau).E()-Pions1.at(1).E(),w);
	    Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Muon).E()-Pions2.at(1).E(),w);
	    Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau).E()-Pions1.at(1).E())/Pions1.at(1).E(),w);
	    Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Muon).E()-Pions2.at(1).E())/Pions2.at(1).E(),w);
	  }
	else
	  {
	    Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau).E()-Pions2.at(1).E(),w);
	    Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Muon).E()-Pions1.at(1).E(),w);
	    Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau).E()-Pions2.at(1).E())/Pions2.at(1).E(),w);
	    Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Muon).E()-Pions1.at(1).E())/Pions1.at(1).E(),w);
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
	if((sqrt((TauP4.Eta()-Tau1Truth.Eta())*(TauP4.Eta()-Tau1Truth.Eta())+(TauP4.Phi()-Tau1Truth.Phi())*(TauP4.Phi()-Tau1Truth.Phi())))<0.5 && (sqrt((MuonP4.Eta()-Tau2Truth.Eta())*(MuonP4.Eta()-Tau2Truth.Eta())+(MuonP4.Phi()-Tau2Truth.Phi())*(MuonP4.Phi()-Tau2Truth.Phi())))<0.5)
	  {
	    Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau).E()-Pions1.at(1).E(),w);
	    Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau).E()-Pions1.at(1).E())/Pions1.at(1).E(),w);
	  }
	else{ Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Muon).E()-Pions1.at(1).E(),w);Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Muon).E()-Pions1.at(1).E())/Pions1.at(1).E(),w);}
      }
      if(Ntp->CheckDecayID(4,5)){
	Tau1Truth=Ntp->GetTruthTauLV(4,0);
	Tau2Truth=Ntp->GetTruthTauLV(5,1);
	Pions1=Ntp->GetTruthPionsFromRho(0);
	TruthDecayFromTau1=Pions1.at(0)+Pions1.at(1);
	Pions2=Ntp->GetTruthPionsFromA1(1);
	TruthDecayFromTau2=Pions2.at(0)+Pions2.at(1)+Pions2.at(2);decay=1;if(TruthDecayFromTau1==TruthDecayFromTau2){cout<<"4 5"<<endl;TruthDecayFromTau1.Print();TruthDecayFromTau2.Print();}
	if((sqrt((TauP4.Eta()-Tau1Truth.Eta())*(TauP4.Eta()-Tau1Truth.Eta())+(TauP4.Phi()-Tau1Truth.Phi())*(TauP4.Phi()-Tau1Truth.Phi())))<0.5 && (sqrt((MuonP4.Eta()-Tau2Truth.Eta())*(MuonP4.Eta()-Tau2Truth.Eta())+(MuonP4.Phi()-Tau2Truth.Phi())*(MuonP4.Phi()-Tau2Truth.Phi())))<0.5)
	  {
	    Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Tau).E()-Pions1.at(1).E(),w);
	    Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Tau).E()-Pions1.at(1).E())/Pions1.at(1).E(),w);
	  }
	else{ Pi0EnergyRes.at(t).Fill(Ntp->NeutralDaughters_P4(Muon).E()-Pions1.at(1).E(),w);Pi0EnergyResPull.at(t).Fill((Ntp->NeutralDaughters_P4(Muon).E()-Pions1.at(1).E())/Pions1.at(1).E(),w);}
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
	    //Phitruthpipi.at(t).Fill(Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	    //PhiSvFitRespipi.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	    //PhiVisRespipi.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	  }
	  if(id==10410333) {	  
	    //Phitruthpirho.at(t).Fill(Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	    // PhiSvFitRespirho.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	    //PhiVisRespirho.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	  }
	  if(id==10410433) {	  
	    // Phitruthrhorho.at(t).Fill(Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	    //PhiSvFitResrhorho.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	    // PhiVisResrhorho.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	  }
	  if(id==10330533) {	  
	    //Phitruthpia1.at(t).Fill(Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	    // PhiSvFitRespia1.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	    // PhiVisRespia1.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	  }
	  if(id==10430533) {	  
	    // Phitruthrhoa1.at(t).Fill(Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	    // PhiSvFitResrhoa1.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	    // PhiVisResrhoa1.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	  }
	  if(id==30530533) {	  
	    // Phitrutha1a1.at(t).Fill(Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	    // PhiSvFitResa1a1.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	    // PhiVisResa1a1.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	  }
	  if(id==10110333 || id==10210333) {	    
	    Phitruthlpi.at(t).Fill(Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	    //  PhiSvFitReslpi.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	    // PhiVisReslpi.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	  }
	  if(id==10110433 || id==10210433) {	   
	    Phitruthlrho.at(t).Fill(Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	    //  PhiSvFitReslrho.at(t).Fill(Ntp->DeltaPhi(xyprotonsvfit.Phi(),xytauplussvfit.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	    //  PhiVisReslrho.at(t).Fill(Ntp->DeltaPhi(xyprotonvis.Phi(),xytauplusvis.Phi())-Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
	  }
	  if(id==10130533 || id==10230533) {	   
	    Phitruthla1.at(t).Fill(Ntp->DeltaPhi(xyprotontruth.Phi(),xytauplustruth.Phi()),w);
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
void  ZMuTau::Finish(){
  if(mode == RECONSTRUCT){
    std::cout<<" Starting Finish!  " <<std::endl;
    
    std::cout<<"A  Data  "<< NQCD.at(0).GetBinContent(1) << std::endl;
    std::cout<<"B  Data  "<< NQCD.at(0).GetBinContent(2) << std::endl;
    std::cout<<"C  Data  "<< NQCD.at(0).GetBinContent(3) << std::endl;
    std::cout<<"D  Data  "<< NQCD.at(0).GetBinContent(4) << std::endl;
    SkimConfig SC;
    SC.ApplySkimEfficiency(types,Npassed, Npassed_noweight);
    
    double W_ALow=0;
    std::vector<double> ALow;
    double QCD_ALow=0;
    
    std::vector<double> QCD_Integral_B;
    double QCD_IntegralMC_B=0;
    double QCD_Integral_B_Data = 0;
    std::vector<double> BLow;
    double MC_BLow=0;
    double BLow_Data = 0;
    double W_BLow=0;
    double QCD_BLow=0;
    
    std::vector<double> QCD_Integral_C;
    double QCD_IntegralMC_C=0;
    double QCD_Integral_C_Data = 0;
    std::vector<double> AHigh;
    double MC_AHigh=0;
    double AHigh_Data = 0;
    
    std::vector<double> QCD_Integral_D;
    double QCD_IntegralMC_D=0;
    double QCD_Integral_D_Data = 0;
    std::vector<double> BHigh;
    double MC_BHigh=0;
    double BHigh_Data = 0;
    
    double W_MC_ALow=0;
    double W_MC_BLow=0;
    double W_MC_AHigh=0;
    double W_MC_BHigh=0;
    double fBLW=0;
    double fALW=0;
    double W_Factor_HighLow=0;
    double QCD_Integral_C_W;
    double QCD_Integral_D_W;
    //Get Yields in ABCD for QCD Scalefactor
    for(unsigned i=0;i<CrossSectionandAcceptance.size();i++){
      QCD_Integral_B.push_back(NQCD.at(i).GetBinContent(2));
      QCD_Integral_C.push_back(NQCD.at(i).GetBinContent(3));
      QCD_Integral_D.push_back(NQCD.at(i).GetBinContent(4));
      BLow.push_back(NWJets.at(i).GetBinContent(2));
      AHigh.push_back(NWJets.at(i).GetBinContent(3));
      BHigh.push_back(NWJets.at(i).GetBinContent(4));
      if(CrossSectionandAcceptance.at(i)>0 ){
	QCD_Integral_B.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
	QCD_Integral_C.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
	QCD_Integral_D.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
	BLow.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
	AHigh.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
	BHigh.at(i) *= CrossSectionandAcceptance.at(i)*Lumi/Npassed.at(i).GetBinContent(0);
      }
      if(HConfig.GetID(i) == DataMCType::W_lnu)
	{
	  QCD_Integral_C_W=QCD_Integral_C.at(i);
	  QCD_Integral_D_W=QCD_Integral_D.at(i);
	  W_MC_ALow+=NWJets.at(i).GetBinContent(1);
	  W_MC_BLow+=BLow.at(i);
	  W_MC_AHigh+=AHigh.at(i);
	  W_MC_BHigh+=BHigh.at(i);
	  W_Factor_HighLow+=(NWJetsRelaxed.at(i).GetBinContent(1))/(NWJetsRelaxed.at(i).GetBinContent(2));
	  fBLW+=W_MC_BLow/W_MC_BHigh;
	  fALW+=W_MC_ALow/W_MC_AHigh;
	}
    }
    for(unsigned i=0;i<CrossSectionandAcceptance.size();i++) {
      if(HConfig.GetID(i) == DataMCType::Data){
	QCD_Integral_B_Data += QCD_Integral_B.at(i);
	QCD_Integral_C_Data += QCD_Integral_C.at(i);
	QCD_Integral_D_Data += QCD_Integral_D.at(i);
	BLow_Data += BLow.at(i);
	AHigh_Data += AHigh.at(i);
	BHigh_Data += BHigh.at(i);
      }
      if(CrossSectionandAcceptance.at(i)>0){
	QCD_IntegralMC_B  += QCD_Integral_B.at(i);
	QCD_IntegralMC_C  += QCD_Integral_C.at(i);
	QCD_IntegralMC_D  += QCD_Integral_D.at(i);
	MC_BLow  += BLow.at(i);
	MC_AHigh  += AHigh.at(i);
	MC_BHigh  += BHigh.at(i);
      }
    }
    //W_ALow=(AHigh_Data-MC_AHigh)*fALW;
    W_BLow=W_MC_BLow;//(BHigh_Data-MC_BHigh)*fBLW;
    double OS2SS=(QCD_Integral_C_Data - QCD_IntegralMC_C-QCD_Integral_C_W )/(QCD_Integral_D_Data - QCD_IntegralMC_D-QCD_Integral_D_W);
    double QCD_AHigh=OS2SS*(BHigh_Data-MC_BHigh);
    double W_AHigh=AHigh_Data-MC_AHigh-QCD_AHigh;
    double W_ScaleFactor=W_AHigh/W_MC_AHigh;
    W_ALow=W_Factor_HighLow*W_AHigh;
    QCD_BLow=BLow_Data-MC_BLow-W_BLow*W_ScaleFactor;
    QCD_ALow=QCD_BLow*OS2SS;
    //double QCD_A=(QCD_Integral_B_Data-W_IntegralMC_B+W_MC_B-W_ScaleFactor*W_MC_B)*OS2SS;

    std::cout << "QCD_Integral_B_Data is: " << QCD_Integral_B_Data << std::endl;
    std::cout << "QCD_Integral_C_Data is: " << QCD_Integral_C_Data << std::endl;
    std::cout << "QCD_Integral_D_Data is: " << QCD_Integral_D_Data << std::endl;
    std::cout << "QCD_IntegralMC_B is: " << QCD_IntegralMC_B << std::endl;
    std::cout << "QCD_IntegralMC_C is: " << QCD_IntegralMC_C << std::endl;
    std::cout << "QCD_IntegralMC_D is: " << QCD_IntegralMC_D << std::endl;
    std::cout << "QCD_AHigh: " << QCD_AHigh << std::endl;
    std::cout << "QCD_ALow: " << QCD_ALow << std::endl;
    std::cout << "QCD_BLow: " << QCD_BLow << std::endl;
    std::cout << "OS/SS QCD Sample: " << OS2SS << std::endl;
    std::cout << "BLow_Data: " <<BLow_Data << std::endl;
    std::cout << "AHigh_Data: " << AHigh_Data << std::endl;
    std::cout << "BHigh_Data: " << BHigh_Data << std::endl;
    std::cout << "MC_BLow: " << MC_BLow << std::endl;
    std::cout << "MC_AHigh: " << MC_AHigh << std::endl;
    std::cout << "MC_BHigh: " << MC_BHigh << std::endl;
    std::cout << "W_ALow: " << W_ALow << std::endl;
    std::cout << "W_BLow: " << W_BLow << std::endl;
    std::cout << "W_AHigh: " << W_AHigh << std::endl;
    std::cout << "W_MC_ALow: " <<W_MC_ALow << std::endl;
    std::cout << "W_MC_AHigh: " <<W_MC_AHigh << std::endl;
    std::cout << "W_MC_BLow: " <<W_MC_BLow << std::endl;
    std::cout << "W_MC_BHigh: " <<W_MC_BHigh << std::endl;
    std::cout << "W_ScaleFactor: " << W_ScaleFactor << std::endl;
    
    
    ScaleAllHistOfType(HConfig.GetType(DataMCType::QCD),QCD_ALow/Nminus0.at(0).at(HConfig.GetType(DataMCType::QCD)).Integral());
    ScaleAllHistOfType(HConfig.GetType(DataMCType::W_lnu),W_ALow/Nminus0.at(0).at(HConfig.GetType(DataMCType::W_lnu)).Integral());
  }

  Selection::Finish();
}
