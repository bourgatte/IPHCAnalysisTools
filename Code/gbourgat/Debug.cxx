#include "Debug.h"
#include "TLorentzVector.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
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



Debug::Debug(TString Name_, TString id_):
  Selection(Name_,id_),
  DataMC_Corr(true,true,false),
  tauTrgSF("tight")
{
  ChargeSumDummy = -999;
  selMuon_IsoDummy = 999.;
}

Debug::~Debug(){
  for(unsigned int j=0; j<Npassed.size(); j++){
    Logger(Logger::Info) << "Selection Summary before: "
			 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
			 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  Debug::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==Id_and_Kin)          cut.at(Id_and_Kin)=1;
    if(i==PairCharge)          cut.at(PairCharge)=1.;  
  }
  // Setup cut plots
  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
    if(i==Id_and_Kin){
      title.at(i)="Id and Kinematic";
      hlabel="Number of Event with good particles";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Id_and_Kin_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Id_and_Kin_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==PairCharge){
      title.at(i)="Pair Charge";
      hlabel="is pair OS";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PairCharge_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PairCharge_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  
  h_SVFitMass = HConfig.GetTH1D(Name+"_SVFitMass","SVFitMass",100,70.,160.,"m_{SVfit}(#tau_{h},#tau_{h})/GeV");

  polarimetricAcopAngle=HConfig.GetTH1D(Name+"_polarimetricAcopAngle","acoplanar angle GEF",60,0.,2*TMath::Pi(),"acoplanar angle GEF","Events");  
  polarimetricAcopAngleNoBS=HConfig.GetTH1D(Name+"_polarimetricAcopAngleNoBS","acoplanar angle refitted PV without BS",60,0.,2*TMath::Pi(),"acoplanar angle refitted PV without BS","Events");
  polarimetricAcopAngleNoBSPionsRefit=HConfig.GetTH1D(Name+"_polarimetricAcopAngleNoBSPionsRefit","acoplanar angle refitted PV without BS and refitted pions",60,0.,2*TMath::Pi(),"acoplanar angle refitted PV without BS and refitted pions","Events");
  polarimetricAcopAngleTruthA1=HConfig.GetTH1D(Name+"_polarimetricAcopAngleTruthA1","acoplanar angle Truth",60,0,2*TMath::Pi(),"acoplanar angle Truth","Events");
  polarimetricAcopAngleSVFitA1=HConfig.GetTH1D(Name+"_polarimetricAcopAngleSVFitA1","acoplanar angle FastMTT",60,0.,2*TMath::Pi(),"acoplanar angle FastMTT","Events");
  polarimetricAcopAngleOld=HConfig.GetTH1D(Name+"_polarimetricAcopAngleOld","acoplanar angle old refitted PV with BS",60,0.,2*TMath::Pi(),"acoplanar angle old refitted PV with BS","Events");  
  polarimetricAcopAngleOldPionsRefit=HConfig.GetTH1D(Name+"_polarimetricAcopAngleOldPionsRefit","acoplanar angle old refitted PV with BS and refitted pions",60,0.,2*TMath::Pi(),"acoplanar angle old refitted PV with BS and refitted pions","Events");
  polarimetricAcopAngleTIP=HConfig.GetTH1D(Name+"_polarimetricAcopAngleTIP","acoplanar angle TIP refitted PV with BS",60,0.,2*TMath::Pi(),"acoplanar angle TIP refitted PV with BS","Events");  
  polarimetricAcopAngleTIPPionsRefit=HConfig.GetTH1D(Name+"_polarimetricAcopAngleTIPPionsRefit","acoplanar angle TIP refitted PV with BS and refitted pions",60,0.,2*TMath::Pi(),"acoplanar angle TIP refitted PV with BS and refitted pions","Events");

  PVXResol=HConfig.GetTH1D(Name+"_PVXResol","PV_{X} pull",50,-0.1,0.1,"PV_{X} pull","Events");
  PVXNoBSResol=HConfig.GetTH1D(Name+"_PVXNoBSResol","no BS refitted PV_{X} pull",50,-0.1,0.1,"no BS refitted PV_{X} pull","Events");
  PVXOldRefitResol=HConfig.GetTH1D(Name+"_PVXOldRefitResol","old refitted PV_{X} pull",50,-0.1,0.1,"old refitted PV_{X} pull","Events");
  PVXTIPRefitResol=HConfig.GetTH1D(Name+"_PVXTIPRefitResol","TIP refitted PV_{X} pull",50,-0.1,0.1,"TIP refitted PV_{X} pull","Events");

  PVYResol=HConfig.GetTH1D(Name+"_PVYResol","PV_{Y} pull",50,-0.1,0.1,"PV_{Y} pull","Events");
  PVYNoBSResol=HConfig.GetTH1D(Name+"_PVYNoBSResol","no BS refitted PV_{Y} pull",50,-0.1,0.1,"no BS refitted PV_{Y} pull","Events");
  PVYOldRefitResol=HConfig.GetTH1D(Name+"_PVYOldRefitResol","old refitted PV_{Y} pull",50,-0.1,0.1,"old refitted PV_{Y} pull","Events");
  PVYTIPRefitResol=HConfig.GetTH1D(Name+"_PVYTIPRefitResol","TIP refitted PV_{Y} pull",50,-0.1,0.1,"TIP refitted PV_{Y} pull","Events");

  PVZResol=HConfig.GetTH1D(Name+"_PVZResol","PV_{Z} pull",50,-0.01,0.01,"PV_{Z} pull","Events");
  PVZNoBSResol=HConfig.GetTH1D(Name+"_PVZNoBSResol","no BS refitted PV_{Z} pull",50,-0.01,0.01,"no BS refitted PV_{Z} pull","Events");
  PVZOldRefitResol=HConfig.GetTH1D(Name+"_PVZOldRefitResol","old refitted PV_{Z} pull",50,-0.01,0.01,"old refitted PV_{Z} pull","Events");
  PVZTIPRefitResol=HConfig.GetTH1D(Name+"_PVZTIPRefitResol","TIP refitted PV_{Z} pull",50,-0.01,0.01,"TIP refitted PV_{Z} pull","Events");
  
  TauSVFitPxResPull=HConfig.GetTH1D(Name+"_TauSVFitPxResPull","p_{X} pull of #tau",50,-2.,2.,"p_{X} pull of #tau","Events");
  TauSVFitPyResPull=HConfig.GetTH1D(Name+"_TauSVFitPyResPull","p_{Y} pull of #tau",50,-2.,2.,"p_{Y} pull of #tau","Events");
  TauSVFitPzResPull=HConfig.GetTH1D(Name+"_TauSVFitPzResPull","p_{Z} pull of #tau",50,-2.,2.,"p_{Z} pull of #tau","Events");
  TauSVFitEResPull=HConfig.GetTH1D(Name+"_TauSVFitEResPull","E pull of #tau",50,-2.,2.,"E pull of #tau","Events");

  TauPxResPull=HConfig.GetTH1D(Name+"_TauPxResPull","p_{X} pull of #tau",50,-2.,2.,"p_{X} pull of #tau","Events");
  TauPyResPull=HConfig.GetTH1D(Name+"_TauPyResPull","p_{Y} pull of #tau",50,-2.,2.,"p_{Y} pull of #tau","Events");
  TauPzResPull=HConfig.GetTH1D(Name+"_TauPzResPull","p_{Z} pull of #tau",50,-2.,2.,"p_{Z} pull of #tau","Events");
  TauEResPull=HConfig.GetTH1D(Name+"_TauEResPull","E pull of #tau",50,-2.,2.,"E pull of #tau","Events");

  TauNoBSPxResPull=HConfig.GetTH1D(Name+"_TauNoBSPxResPull","p_{X} pull of #tau NoBS",50,-2.,2.,"p_{X} pull of #tau NoBS","Events");
  TauNoBSPyResPull=HConfig.GetTH1D(Name+"_TauNoBSPyResPull","p_{Y} pull of #tau NoBS",50,-2.,2.,"p_{Y} pull of #tau NoBS","Events");
  TauNoBSPzResPull=HConfig.GetTH1D(Name+"_TauNoBSPzResPull","p_{Z} pull of #tau NoBS",50,-2.,2.,"p_{Z} pull of #tau NoBS","Events");
  TauNoBSEResPull=HConfig.GetTH1D(Name+"_TauNoBSEResPull","E pull of #tau NoBS",50,-2.,2.,"E pull of #tau NoBS","Events");

  TauOldPxResPull=HConfig.GetTH1D(Name+"_TauOldPxResPull","p_{X} pull of #tau Old",50,-2.,2.,"p_{X} pull of #tau Old","Events");
  TauOldPyResPull=HConfig.GetTH1D(Name+"_TauOldPyResPull","p_{Y} pull of #tau Old",50,-2.,2.,"p_{Y} pull of #tau Old","Events");
  TauOldPzResPull=HConfig.GetTH1D(Name+"_TauOldPzResPull","p_{Z} pull of #tau Old",50,-2.,2.,"p_{Z} pull of #tau Old","Events");
  TauOldEResPull=HConfig.GetTH1D(Name+"_TauOldEResPull","E pull of #tau Old",50,-2.,2.,"E pull of #tau Old","Events");
  
  TauTIPPxResPull=HConfig.GetTH1D(Name+"_TauTIPPxResPull","p_{X} pull of #tau TIP",50,-2.,2.,"p_{X} pull of #tau TIP","Events");
  TauTIPPyResPull=HConfig.GetTH1D(Name+"_TauTIPPyResPull","p_{Y} pull of #tau TIP",50,-2.,2.,"p_{Y} pull of #tau TIP","Events");
  TauTIPPzResPull=HConfig.GetTH1D(Name+"_TauTIPPzResPull","p_{Z} pull of #tau TIP",50,-2.,2.,"p_{Z} pull of #tau TIP","Events");
  TauTIPEResPull=HConfig.GetTH1D(Name+"_TauTIPEResPull","E pull of #tau TIP",50,-2.,2.,"E pull of #tau TIP","Events");

  TauNoBSRefittedPionsPxResPull=HConfig.GetTH1D(Name+"_TauNoBSRefittedPionsPxResPull","p_{X} pull of #tau NoBS with refitted pions",50,-2.,2.,"p_{X} pull of #tau NoBS with refitted pions","Events");
  TauNoBSRefittedPionsPyResPull=HConfig.GetTH1D(Name+"_TauNoBSRefittedPionsPyResPull","p_{Y} pull of #tau NoBS with refitted pions",50,-2.,2.,"p_{Y} pull of #tau NoBS with refitted pions","Events");
  TauNoBSRefittedPionsPzResPull=HConfig.GetTH1D(Name+"_TauNoBSRefittedPionsPzResPull","p_{Z} pull of #tau NoBS with refitted pions",50,-2.,2.,"p_{Z} pull of #tau NoBS with refitted pions","Events");
  TauNoBSRefittedPionsEResPull=HConfig.GetTH1D(Name+"_TauNoBSRefittedPionsEResPull","E pull of #tau NoBS with refitted pions",50,-2.,2.,"E pull of #tau NoBS with refitted pions","Events");

  TauOldRefittedPionsPxResPull=HConfig.GetTH1D(Name+"_TauOldRefittedPionsPxResPull","p_{X} pull of #tau Old with refitted pions",50,-2.,2.,"p_{X} pull of #tau Old with refitted pions","Events");
  TauOldRefittedPionsPyResPull=HConfig.GetTH1D(Name+"_TauOldRefittedPionsPyResPull","p_{Y} pull of #tau Old with refitted pions",50,-2.,2.,"p_{Y} pull of #tau Old with refitted pions","Events");
  TauOldRefittedPionsPzResPull=HConfig.GetTH1D(Name+"_TauOldRefittedPionsPzResPull","p_{Z} pull of #tau Old with refitted pions",50,-2.,2.,"p_{Z} pull of #tau Old with refitted pions","Events");
  TauOldRefittedPionsEResPull=HConfig.GetTH1D(Name+"_TauOldRefittedPionsEResPull","E pull of #tau Old with refitted pions",50,-2.,2.,"E pull of #tau Old with refitted pions","Events");
  
  TauTIPRefittedPionsPxResPull=HConfig.GetTH1D(Name+"_TauTIPRefittedPionsPxResPull","p_{X} pull of #tau TIP with refitted pions",50,-2.,2.,"p_{X} pull of #tau TIP with refitted pions","Events");
  TauTIPRefittedPionsPyResPull=HConfig.GetTH1D(Name+"_TauTIPRefittedPionsPyResPull","p_{Y} pull of #tau TIP with refitted pions",50,-2.,2.,"p_{Y} pull of #tau TIP with refitted pions","Events");
  TauTIPRefittedPionsPzResPull=HConfig.GetTH1D(Name+"_TauTIPRefittedPionsPzResPull","p_{Z} pull of #tau TIP with refitted pions",50,-2.,2.,"p_{Z} pull of #tau TIP with refitted pions","Events");
  TauTIPRefittedPionsEResPull=HConfig.GetTH1D(Name+"_TauTIPRefittedPionsEResPull","E pull of #tau TIP with refitted pions",50,-2.,2.,"E pull of #tau TIP with refitted pions","Events");
  
  LeadingPionPxResPull=HConfig.GetTH1D(Name+"_LeadingPionPxResPull","p_{X} pull of leading pion",50,-0.1,0.1,"p_{X} pull of leading pion","Events");
  LeadingPionPyResPull=HConfig.GetTH1D(Name+"_LeadingPionPyResPull","p_{Y} pull of leading pion",50,-0.1,0.1,"p_{Y} pull of leading pion","Events");
  LeadingPionPzResPull=HConfig.GetTH1D(Name+"_LeadingPionPzResPull","p_{Z} pull of leading pion",50,-0.1,0.1,"p_{Z} pull of leading pion","Events");
  LeadingPionEResPull=HConfig.GetTH1D(Name+"_LeadingPionEResPull","E pull of leading pion",50,-0.1,0.1,"E pull of leading pion","Events");

  OppositePionPxResPull=HConfig.GetTH1D(Name+"_OppositePionPxResPull","p_{X} pull of opposite pion",50,-0.1,0.1,"p_{X} pull of opposite pion","Events");
  OppositePionPyResPull=HConfig.GetTH1D(Name+"_OppositePionPyResPull","p_{Y} pull of opposite pion",50,-0.1,0.1,"p_{Y} pull of opposite pion","Events");
  OppositePionPzResPull=HConfig.GetTH1D(Name+"_OppositePionPzResPull","p_{Z} pull of opposite pion",50,-0.1,0.1,"p_{Z} pull of opposite pion","Events");
  OppositePionEResPull=HConfig.GetTH1D(Name+"_OppositePionEResPull","E pull of opposite pion",50,-0.1,0.1,"E pull of opposite pion","Events");

  LeadingRefittedPionPxResPull=HConfig.GetTH1D(Name+"_LeadingRefittedPionPxResPull","p_{X} pull of refitted leading pion",50,-0.1,0.1,"p_{X} pull of refitted leading pion","Events");
  LeadingRefittedPionPyResPull=HConfig.GetTH1D(Name+"_LeadingRefittedPionPyResPull","p_{Y} pull of refitted leading pion",50,-0.1,0.1,"p_{Y} pull of refitted leading pion","Events");
  LeadingRefittedPionPzResPull=HConfig.GetTH1D(Name+"_LeadingRefittedPionPzResPull","p_{Z} pull of refitted leading pion",50,-0.1,0.1,"p_{Z} pull of refitted leading pion","Events");
  LeadingRefittedPionEResPull=HConfig.GetTH1D(Name+"_LeadingRefittedPionEResPull","E pull of refitted leading pion",50,-0.1,0.1,"E pull of refitted leading pion","Events");

  OppositeRefittedPionPxResPull=HConfig.GetTH1D(Name+"_OppositeRefittedPionPxResPull","p_{X} pull of refitted opposite pion",50,-0.1,0.1,"p_{X} pull of refitted opposite pion","Events");
  OppositeRefittedPionPyResPull=HConfig.GetTH1D(Name+"_OppositeRefittedPionPyResPull","p_{Y} pull of refitted opposite pion",50,-0.1,0.1,"p_{Y} pull of refitted opposite pion","Events");
  OppositeRefittedPionPzResPull=HConfig.GetTH1D(Name+"_OppositeRefittedPionPzResPull","p_{Z} pull of refitted opposite pion",50,-0.1,0.1,"p_{Z} pull of refitted opposite pion","Events");
  OppositeRefittedPionEResPull=HConfig.GetTH1D(Name+"_OppositeRefittedPionEResPull","E pull of refitted opposite pion",50,-0.1,0.1,"E pull of refitted opposite pion","Events");
  

  VectPolaSVFitPxResPull=HConfig.GetTH1D(Name+"_VectPolaSVFitPxResPull","h_{X} pull",100,-3.,1.,"h_{X} pull FastMTT","Events");
  VectPolaSVFitPyResPull=HConfig.GetTH1D(Name+"_VectPolaSVFitPyResPull","h_{Y} pull",100,-3.,1.,"h_{Y} pull FastMTT","Events");
  VectPolaSVFitPzResPull=HConfig.GetTH1D(Name+"_VectPolaSVFitPzResPull","h_{Z} pull",100,-3.,1.,"h_{Z} pull FastMTT","Events");
  
  VectPolaPxResPull=HConfig.GetTH1D(Name+"_VectPolaPxResPull","h_{X} pull",100,-3.,1.,"h_{X} pull","Events");
  VectPolaPyResPull=HConfig.GetTH1D(Name+"_VectPolaPyResPull","h_{Y} pull",100,-3.,1.,"h_{Y} pull","Events");
  VectPolaPzResPull=HConfig.GetTH1D(Name+"_VectPolaPzResPull","h_{Z} pull",100,-3.,1.,"h_{Z} pull","Events");

  Unphysical=HConfig.GetTH1D(Name+"_Unphysical","Unphysical events",2,0,2,"Unphysical event","Events");
  
  // DiffTauXLeadingPionXSVFit=HConfig.GetTH1D(Name+"_DiffTauXLeadingPionXSVFit","#tau_{X} - leading pion_{X}",80,0.,80.,"#tau_{X} - leading pion_{X}","Events");
  // DiffTauYLeadingPionYSVFit=HConfig.GetTH1D(Name+"_DiffTauYLeadingPionYSVFit","#tau_{Y} - leading pion_{Y}",80,0.,80.,"#tau_{Y} - leading pion_{Y}","Events");
  // DiffTauZLeadingPionZSVFit=HConfig.GetTH1D(Name+"_DiffTauZLeadingPionZSVFit","#tau_{Z} - leading pion_{Z}",80,0.,80.,"#tau_{Z} - leading pion_{Z}","Events");
  // DiffTauELeadingPionESVFit=HConfig.GetTH1D(Name+"_DiffTauELeadingPionESVFit","#tau_{E} - leading pion_{E}",100,0.,100.,"#tau_{E} - leading pion_{E}","Events");

  // DiffTauXLeadingPionXReco=HConfig.GetTH1D(Name+"_DiffTauXLeadingPionXReco","#tau_{X} - leading pion_{X}",80,0.,80.,"#tau_{X} - leading pion_{X}","Events");
  // DiffTauYLeadingPionYReco=HConfig.GetTH1D(Name+"_DiffTauYLeadingPionYReco","#tau_{Y} - leading pion_{Y}",80,0.,80.,"#tau_{Y} - leading pion_{Y}","Events");
  // DiffTauZLeadingPionZReco=HConfig.GetTH1D(Name+"_DiffTauZLeadingPionZReco","#tau_{Z} - leading pion_{Z}",80,0.,80.,"#tau_{Z} - leading pion_{Z}","Events");
  // DiffTauELeadingPionEReco=HConfig.GetTH1D(Name+"_DiffTauELeadingPionEReco","#tau_{E} - leading pion_{E}",100,0.,100.,"#tau_{E} - leading pion_{E}","Events");
  
  //AcopAngleVSTauPt=HConfig.GetTH2D(Name+"_AcopAngleVSTauPt","acoplanar angle VS Tau Pt",50,0,2.*TMath::Pi(),80,20.,100.,"acoplanar angle","tau Pt");
  //AcopAngleVSTauEta=HConfig.GetTH2D(Name+"_AcopAngleVSTauEta","acoplanar angle VS Tau Eta",50,0,2.*TMath::Pi(),48,-2.4,2.4,"acoplanar angle","tau eta");

  // PxtauSVFitVSPxLeadingPion=HConfig.GetTH2D(Name+"_PxtauSVFitVSPxLeadingPion","tau P_{X} VS leading pion P_{X}",150,-150.,150.,100,-100.,100.,"tau P_{X}","leading pion P_{X}");
  // PytauSVFitVSPyLeadingPion=HConfig.GetTH2D(Name+"_PytauSVFitVSPyLeadingPion","tau P_{Y} VS leading pion P_{Y}",150,-150.,150.,100,-100.,100.,"tau P_{Y}","leading pion P_{Y}");
  // PztauSVFitVSPzLeadingPion=HConfig.GetTH2D(Name+"_PztauSVFitVSPzLeadingPion","tau P_{Z} VS leading pion P_{Z}",400,-400.,400.,300,-300.,300.,"tau P_{Z}","leading pion P_{Z}");
  // EtauSVFitVSELeadingPion=HConfig.GetTH2D(Name+"_EtauSVFitVSELeadingPion","tau energy VS leading pion energy",200,0,400.,150,0.,300.,"tau energy","leading pion energy");
  
  // HxResolSVFitVSDiffPxtauSVFitPxLeadingPion=HConfig.GetTH2D(Name+"_HxResolSVFitVSDiffPxtauSVFitPxLeadingPion"," h_{X} pull VS #tau_{X} - leading pion_{X}",150,-3.,1.,80,0.,80.,"h_{X} pull","#tau_{X} - leading pion_{X}");
  // HxResolSVFitVSDiffPytauSVFitPyLeadingPion=HConfig.GetTH2D(Name+"_HxResolSVFitVSDiffPytauSVFitPyLeadingPion","h_{X} pull VS #tau_{Y} - leading pion_{Y}",150,-3,1,80,0.,80.,"h_{X} pull","#tau_{Y} - leading pion_{Y}");
  // HxResolSVFitVSDiffPztauSVFitPzLeadingPion=HConfig.GetTH2D(Name+"_HxResolSVFitVSDiffPztauSVFitPzLeadingPion","h_{X} pull VS #tau_{Z} - leading pion_{Z}",150,-3,1,80,0.,80.,"h_{X} pull","#tau_{Z} - leading pion_{Z}");
  // HxResolSVFitVSDiffEtauSVFitELeadingPion=HConfig.GetTH2D(Name+"_HxResolSVFitVSDiffEtauSVFitELeadingPion","h_{X} pull VS #tau energy - leading pion energy",150,-3,1,100,0.,100.,"h_{X} pull","#tau energy - leading pion energy");
  
  // HyResolSVFitVSDiffPxtauSVFitPxLeadingPion=HConfig.GetTH2D(Name+"_HyResolSVFitVSDiffPxtauSVFitPxLeadingPion","h_{Y} pull VS #tau_{X} - leading pion_{X}",150,-3.,1.,80,0.,80.,"h_{Y} pull","#tau_{X} - leading pion_{X}");
  // HyResolSVFitVSDiffPytauSVFitPyLeadingPion=HConfig.GetTH2D(Name+"_HyResolSVFitVSDiffPytauSVFitPyLeadingPion","h_{Y} pull VS #tau_{Y} - leading pion_{Y}",150,-3,1,80,0.,80.,"h_{Y} pull","#tau_{Y} - leading pion_{Y}");
  // HyResolSVFitVSDiffPztauSVFitPzLeadingPion=HConfig.GetTH2D(Name+"_HyResolSVFitVSDiffPztauSVFitPzLeadingPion","h_{Y} pull VS #tau_{Z} - leading pion_{Z}",150,-3,1,80,0.,80.,"h_{Y} pull","#tau_{Z} - leading pion_{Z}");
  // HyResolSVFitVSDiffEtauSVFitELeadingPion=HConfig.GetTH2D(Name+"_HyResolSVFitVSDiffEtauSVFitELeadingPion","h_{Y} pull VS #tau energy - leading pion energy",150,-3,1,100,0.,100.,"h_{Y} pull","#tau energy - leading pion energy");

  // HzResolSVFitVSDiffPxtauSVFitPxLeadingPion=HConfig.GetTH2D(Name+"_HzResolSVFitVSDiffPxtauSVFitPxLeadingPion","h_{Z} pull VS #tau_{X} - leading pion_{X}",150,-3.,1.,80,0.,80.,"h_{Z} pull","#tau_{X} - leading pion_{X}");
  // HzResolSVFitVSDiffPytauSVFitPyLeadingPion=HConfig.GetTH2D(Name+"_HzResolSVFitVSDiffPytauSVFitPyLeadingPion","h_{Z} pull VS #tau_{Y} - leading pion_{Y}",150,-3,1,80,0.,80.,"h_{Z} pull","#tau_{Y} - leading pion_{Y}");
  // HzResolSVFitVSDiffPztauSVFitPzLeadingPion=HConfig.GetTH2D(Name+"_HzResolSVFitVSDiffPztauSVFitPzLeadingPion","h_{Z} pull VS #tau_{Z} - leading pion_{Z}",150,-3,1,80,0.,80.,"h_{Z} pull","#tau_{Z} - leading pion_{Z}");
  // HzResolSVFitVSDiffEtauSVFitELeadingPion=HConfig.GetTH2D(Name+"_HzResolSVFitVSDiffEtauSVFitELeadingPion","h_{Z} pull VS #tau energy - leading pion energy",150,-3,1,100,0.,100.,"h_{Z} pull","#tau energy - leading pion energy");


  // PxtauRecoVSPxLeadingPion=HConfig.GetTH2D(Name+"_PxtauRecoVSPxLeadingPion","tau P_{X} VS leading pion P_{X}",150,-150.,150.,100,-100.,100.,"tau P_{X}","leading pion P_{X}");
  // PytauRecoVSPyLeadingPion=HConfig.GetTH2D(Name+"_PytauRecoVSPyLeadingPion","tau P_{Y} VS leading pion P_{Y}",150,-150.,150.,100,-100.,100.,"tau P_{Y}","leading pion P_{Y}");
  // PztauRecoVSPzLeadingPion=HConfig.GetTH2D(Name+"_PztauRecoVSPzLeadingPion","tau P_{Z} VS leading pion P_{Z}",400,-400.,400.,300,-300.,300.,"tau P_{Z}","leading pion P_{Z}");
  // EtauRecoVSELeadingPion=HConfig.GetTH2D(Name+"_EtauRecoVSELeadingPion","tau energy VS leading pion energy",200,0,400.,150,0.,300.,"tau energy","leading pion energy");
  

  // HxResolRecoVSDiffPxtauRecoPxLeadingPion=HConfig.GetTH2D(Name+"_HxResolRecoVSDiffPxtauRecoPxLeadingPion","h_{X} pull VS #tau_{X} - leading pion_{X}",150,-3.,1.,80,0.,80.,"h_{X} pull","#tau_{X} - leading pion_{X}");
  // HxResolRecoVSDiffPytauRecoPyLeadingPion=HConfig.GetTH2D(Name+"_HxResolRecoVSDiffPytauRecoPyLeadingPion","h_{X} pull VS #tau_{Y} - leading pion_{Y}",150,-3,1,80,0.,80.,"h_{X} pull","#tau_{Y} - leading pion_{Y}");
  // HxResolRecoVSDiffPztauRecoPzLeadingPion=HConfig.GetTH2D(Name+"_HxResolRecoVSDiffPztauRecoPzLeadingPion","h_{X} pull VS #tau_{Z} - leading pion_{Z}",150,-3,1,80,0.,80.,"h_{X} pull","#tau_{Z} - leading pion_{Z}");
  // HxResolRecoVSDiffEtauRecoELeadingPion=HConfig.GetTH2D(Name+"_HxResolRecoVSDiffEtauRecoELeadingPion","h_{X} pull VS #tau energy - leading pion energy",150,-3,1,100,0.,100.,"h_{X} pull","#tau energy - leading pion energy");
  
  // HyResolRecoVSDiffPxtauRecoPxLeadingPion=HConfig.GetTH2D(Name+"_HyResolRecoVSDiffPxtauRecoPxLeadingPion","h_{Y} pull VS #tau_{X} - leading pion_{X}",150,-3.,1.,80,0.,80.,"h_{Y} pull","#tau_{X} - leading pion_{X}");
  // HyResolRecoVSDiffPytauRecoPyLeadingPion=HConfig.GetTH2D(Name+"_HyResolRecoVSDiffPytauRecoPyLeadingPion","h_{Y} pull VS #tau_{Y} - leading pion_{Y}",150,-3,1,80,0.,80.,"h_{Y} pull","#tau_{Y} - leading pion_{Y}");
  // HyResolRecoVSDiffPztauRecoPzLeadingPion=HConfig.GetTH2D(Name+"_HyResolRecoVSDiffPztauRecoPzLeadingPion","h_{Y} pull VS #tau_{Z} - leading pion_{Z}",150,-3,1,80,0.,80.,"h_{Y} pull","#tau_{Z} - leading pion_{Z}");
  // HyResolRecoVSDiffEtauRecoELeadingPion=HConfig.GetTH2D(Name+"_HyResolRecoVSDiffEtauRecoELeadingPion","h_{Y} pull VS #tau energy - leading pion energy",150,-3,1,100,0.,100.,"h_{Y} pull","#tau energy - leading pion energy");

  // HzResolRecoVSDiffPxtauRecoPxLeadingPion=HConfig.GetTH2D(Name+"_HzResolRecoVSDiffPxtauRecoPxLeadingPion","h_{Z} pull VS #tau_{X} - leading pion_{X}",150,-3.,1.,80,0.,80.,"h_{Z} pull","#tau_{X} - leading pion_{X}");
  // HzResolRecoVSDiffPytauRecoPyLeadingPion=HConfig.GetTH2D(Name+"_HzResolRecoVSDiffPytauRecoPyLeadingPion","h_{Z} pull VS #tau_{Y} - leading pion_{Y}",150,-3,1,80,0.,80.,"h_{Z} pull","#tau_{Y} - leading pion_{Y}");
  // HzResolRecoVSDiffPztauRecoPzLeadingPion=HConfig.GetTH2D(Name+"_HzResolRecoVSDiffPztauRecoPzLeadingPion","h_{Z} pull VS #tau_{Z} - leading pion_{Z}",150,-3,1,80,0.,80.,"h_{Z} pull","#tau_{Z} - leading pion_{Z}");
  // HzResolRecoVSDiffEtauRecoELeadingPion=HConfig.GetTH2D(Name+"_HzResolRecoVSDiffEtauRecoELeadingPion","h_{Z} pull VS #tau energy - leading pion energy",150,-3,1,100,0.,100.,"h_{Z} pull","#tau energy - leading pion energy");
  
  
  
  // PxtauSVFitVSPxOppositePion=HConfig.GetTH2D(Name+"_PxtauSVFitVSPxOppositePion","PxtauSVFitVSPxOppositePion",150,-150.,150.,100,-100.,100.,"PxtauSVFitVSPxOppositePion","Events");
  // PytauSVFitVSPyOppositePion=HConfig.GetTH2D(Name+"_PytauSVFitVSPyOppositePion","PytauSVFitVSPyOppositePion",150,-150.,150.,100,-100.,100.,"PytauSVFitVSPyOppositePion","Events");
  // PztauSVFitVSPzOppositePion=HConfig.GetTH2D(Name+"_PztauSVFitVSPzOppositePion","PztauSVFitVSPzOppositePion",400,-400.,400.,300,-300.,300.,"PztauSVFitVSPzOppositePion","Events");
  // EtauSVFitVSEOppositePion=HConfig.GetTH2D(Name+"_EtauSVFitVSEOppositePion","EtauSVFitVSEOppositePion",200,0,400.,150,0.,300.,"EtauSVFitVSEOppositePion","Events");
  
  // HxResolSVFitVSDiffPxtauSVFitPxOppositePion=HConfig.GetTH2D(Name+"_HxResolSVFitVSDiffPxtauSVFitPxOppositePion","HxResolSVFitVSDiffPxtauSVFitPxOppositePion",150,-3.,1.,80,0.,80.,"HxResolSVFitVSDiffPxtauSVFitPxOppositePion","Events");
  // HxResolSVFitVSDiffPytauSVFitPyOppositePion=HConfig.GetTH2D(Name+"_HxResolSVFitVSDiffPytauSVFitPyOppositePion","HxResolSVFitVSDiffPytauSVFitPyOppositePion",150,-3.,1.,80,0.,80.,"HxResolSVFitVSDiffPytauSVFitPyOppositePion","Events");
  // HxResolSVFitVSDiffPztauSVFitPzOppositePion=HConfig.GetTH2D(Name+"_HxResolSVFitVSDiffPztauSVFitPzOppositePion","HxResolSVFitVSDiffPztauSVFitPzOppositePion",150,-3.,1.,80,0.,80.,"HxResolSVFitVSDiffPztauSVFitPzOppositePion","Events");
  // HxResolSVFitVSDiffEtauSVFitEOppositePion=HConfig.GetTH2D(Name+"_HxResolSVFitVSDiffEtauSVFitEOppositePion","HxResolSVFitVSDiffEtauSVFitEOppositePion",150,-3,1,100,0.,100.,"HxResolSVFitVSDiffEtauSVFitEOppositePion","Events");

   
  // PxtauRecoVSPxOppositePion=HConfig.GetTH2D(Name+"_PxtauRecoVSPxOppositePion","PxtauRecoVSPxOppositePion",150,-150.,150.,100,-100.,100.,"PxtauRecoVSPxOppositePion","Events");
  // PytauRecoVSPyOppositePion=HConfig.GetTH2D(Name+"_PytauRecoVSPyOppositePion","PytauRecoVSPyOppositePion",150,-150.,150.,100,-100.,100.,"PytauRecoVSPyOppositePion","Events");
  // PztauRecoVSPzOppositePion=HConfig.GetTH2D(Name+"_PztauRecoVSPzOppositePion","PztauRecoVSPzOppositePion",400,-400.,400.,300,-300.,300.,"PztauRecoVSPzOppositePion","Events");
  // EtauRecoVSEOppositePion=HConfig.GetTH2D(Name+"_EtauRecoVSEOppositePion","EtauRecoVSEOppositePion",200,0,400.,150,0.,300.,"EtauRecoVSEOppositePion","Events");
  
  // HxResolRecoVSDiffPxtauRecoPxOppositePion=HConfig.GetTH2D(Name+"_HxResolRecoVSDiffPxtauRecoPxOppositePion","HxResolRecoVSDiffPxtauRecoPxOppositePion",150,-3.,1.,80,0.,80.,"HxResolRecoVSDiffPxtauRecoPxOppositePion","Events");
  // HxResolRecoVSDiffPytauRecoPyOppositePion=HConfig.GetTH2D(Name+"_HxResolRecoVSDiffPytauRecoPyOppositePion","HxResolRecoVSDiffPytauRecoPyOppositePion",150,-3.,1.,80,0.,80.,"HxResolRecoVSDiffPytauRecoPyOppositePion","Events");
  // HxResolRecoVSDiffPztauRecoPzOppositePion=HConfig.GetTH2D(Name+"_HxResolRecoVSDiffPztauRecoPzOppositePion","HxResolRecoVSDiffPztauRecoPzOppositePion",150,-3.,1.,80,0.,80.,"HxResolRecoVSDiffPztauRecoPzOppositePion","Events");
  // HxResolRecoVSDiffEtauRecoEOppositePion=HConfig.GetTH2D(Name+"_HxResolRecoVSDiffEtauRecoEOppositePion","HxResolRecoVSDiffEtauRecoEOppositePion",150,-3,1,100,0.,100.,"HxResolRecoVSDiffEtauRecoEOppositePion","Events");
   
  
  // HxResolSVFitVSTauPt=HConfig.GetTH2D(Name+"_HxResolSVFitVSTauPt","h_{X} pull VS tau Pt",150,-3,1,100,0.,100.,"h_{X} pull","tau Pt");
  // HxResolSVFitVSTauEta=HConfig.GetTH2D(Name+"_HxResolSVFitVSTauEta","h_{X} pull VS tau eta",150,-3,1,48,-2.4,2.4,"h_{X} pull","tau eta");
  // HxResolSVFitVSPxTau=HConfig.GetTH2D(Name+"_HxResolSVFitVSPxTau","h_{X} pull VS tau P_[X}",150,-3,1,150,-150,150,"h_{X} pull","tau P_[X}");
  // HxResolSVFitVSPyTau=HConfig.GetTH2D(Name+"_HxResolSVFitVSPyTau","h_{X} pull VS tau P_{Y}",150,-3,1,150,-150,150,"h_{X} pull","tau P_{Y}");
  // HxResolSVFitVSPzTau=HConfig.GetTH2D(Name+"_HxResolSVFitVSPzTau","h_{X} pull VS tau P_{Z}",150,-3,1,400,-400.,400.,"h_{X} pull","tau P_{Z}");
  // HxResolSVFitVSETau=HConfig.GetTH2D(Name+"_HxResolSVFitVSETau","h_{X} pull VS tau energy",150,-3,1,200,0,400.,"h_{X} pull","tau energy");
  // HxResolSVFitVSGJAngle=HConfig.GetTH2D(Name+"_HxResolSVFitVSGJAngle","h_{X} pull VS ",150,-3,1,100,0,0.03,"h_{X} pull","");
  // HxResolSVFitVSMH=HConfig.GetTH2D(Name+"_HxResolSVFitVSMH","h_{X} pull VS M_{H}",150,-3,1,100,120,130,"h_{X} pull","M_{H}");
  // HxResolSVFitVSSumPtPions=HConfig.GetTH2D(Name+"_HxResolSVFitVSSumPtPions","h_{X} pull VS pions Pt sum",150,-3,1,100,0,100,"h_{X} pull","pions Pt sum");
  // HxResolSVFitVSSumPxPions=HConfig.GetTH2D(Name+"_HxResolSVFitVSSumPxPions","h_{X} pull VS pions_{X} sum",150,-3,1,80,0.,80.,"h_{X} pull","pions_{X} sum");
  // HxResolSVFitVSSumPyPions=HConfig.GetTH2D(Name+"_HxResolSVFitVSSumPyPions","h_{X} pull VS pions_{Y} sum",150,-3,1,80,0.,80.,"h_{X} pull","pions_{Y} sum");
  // HxResolSVFitVSSumPzPions=HConfig.GetTH2D(Name+"_HxResolSVFitVSSumPzPions","h_{X} pull VS pions_{Z} sum",150,-3,1,80,0.,80.,"h_{X} pull","pions_{Z} sum");
  // HxResolSVFitVSSumEPions=HConfig.GetTH2D(Name+"_HxResolSVFitVSSumEPions","h_{X} pull VS pions energy sum",150,-3,1,150,0,300,"h_{X} pull","pions energy sum");

  // HyResolSVFitVSTauPt=HConfig.GetTH2D(Name+"_HyResolSVFitVSTauPt","h_{Y} pull VS tau Pt",150,-3,1,100,0.,100.,"h_{Y} pull","tau Pt");
  // HyResolSVFitVSTauEta=HConfig.GetTH2D(Name+"_HyResolSVFitVSTauEta","h_{Y} pull VS tau eta",150,-3,1,48,-2.4,2.4,"h_{Y} pull","tau eta");
  // HyResolSVFitVSPxTau=HConfig.GetTH2D(Name+"_HyResolSVFitVSPxTau","h_{Y} pull VS tau P_[X}",150,-3,1,150,-150,150,"h_{Y} pull","tau P_[X}");
  // HyResolSVFitVSPyTau=HConfig.GetTH2D(Name+"_HyResolSVFitVSPyTau","h_{Y} pull VS tau P_{Y}",150,-3,1,150,-150,150,"h_{Y} pull","tau P_{Y}");
  // HyResolSVFitVSPzTau=HConfig.GetTH2D(Name+"_HyResolSVFitVSPzTau","h_{Y} pull VS tau P_{Z}",150,-3,1,400,-400.,400.,"h_{Y} pull","tau P_{Z}");
  // HyResolSVFitVSETau=HConfig.GetTH2D(Name+"_HyResolSVFitVSETau","h_{Y} pull VS tau energy",150,-3,1,200,0,400.,"h_{Y} pull","tau energy");
  // HyResolSVFitVSGJAngle=HConfig.GetTH2D(Name+"_HyResolSVFitVSGJAngle","h_{Y} pull VS ",150,-3,1,100,0,0.03,"h_{Y} pull","");
  // HyResolSVFitVSMH=HConfig.GetTH2D(Name+"_HyResolSVFitVSMH","h_{Y} pull VS M_{H}",150,-3,1,100,120,130,"h_{Y} pull","M_{H}");
  // HyResolSVFitVSSumPtPions=HConfig.GetTH2D(Name+"_HyResolSVFitVSSumPtPions","h_{Y} pull VS pions Pt sum",150,-3,1,100,0,100,"h_{Y} pull","pions Pt sum");
  // HyResolSVFitVSSumPxPions=HConfig.GetTH2D(Name+"_HyResolSVFitVSSumPxPions","h_{Y} pull VS pions_{X} sum",150,-3,1,80,0.,80.,"h_{Y} pull","pions_{X} sum");
  // HyResolSVFitVSSumPyPions=HConfig.GetTH2D(Name+"_HyResolSVFitVSSumPyPions","h_{Y} pull VS pions_{Y} sum",150,-3,1,80,0.,80.,"h_{Y} pull","pions_{Y} sum");
  // HyResolSVFitVSSumPzPions=HConfig.GetTH2D(Name+"_HyResolSVFitVSSumPzPions","h_{Y} pull VS pions_{Z} sum",150,-3,1,80,0.,80.,"h_{Y} pull","pions_{Z} sum");
  // HyResolSVFitVSSumEPions=HConfig.GetTH2D(Name+"_HyResolSVFitVSSumEPions","h_{Y} pull VS pions energy sum",150,-3,1,150,0,300,"h_{Y} pull","pions energy sum");


  // HzResolSVFitVSTauPt=HConfig.GetTH2D(Name+"_HzResolSVFitVSTauPt","h_{Z} pull VS tau Pt",150,-3,1,100,0.,100.,"h_{Z} pull","tau Pt");
  // HzResolSVFitVSTauEta=HConfig.GetTH2D(Name+"_HzResolSVFitVSTauEta","h_{Z} pull VS tau eta",150,-3,1,48,-2.4,2.4,"h_{Z} pull","tau eta");
  // HzResolSVFitVSPxTau=HConfig.GetTH2D(Name+"_HzResolSVFitVSPxTau","h_{Z} pull VS tau P_[X}",150,-3,1,150,-150,150,"h_{Z} pull","tau P_[X}");
  // HzResolSVFitVSPyTau=HConfig.GetTH2D(Name+"_HzResolSVFitVSPyTau","h_{Z} pull VS tau P_{Y}",150,-3,1,150,-150,150,"h_{Z} pull","tau P_{Y}");
  // HzResolSVFitVSPzTau=HConfig.GetTH2D(Name+"_HzResolSVFitVSPzTau","h_{Z} pull VS tau P_{Z}",150,-3,1,400,-400.,400.,"h_{Z} pull","tau P_{Z}");
  // HzResolSVFitVSETau=HConfig.GetTH2D(Name+"_HzResolSVFitVSETau","h_{Z} pull VS tau energy",150,-3,1,200,0,400.,"h_{Z} pull","tau energy");
  // HzResolSVFitVSGJAngle=HConfig.GetTH2D(Name+"_HzResolSVFitVSGJAngle","h_{Z} pull VS ",150,-3,1,100,0,0.03,"h_{Z} pull","");
  // HzResolSVFitVSMH=HConfig.GetTH2D(Name+"_HzResolSVFitVSMH","h_{Z} pull VS M_{H}",150,-3,1,100,120,130,"h_{Z} pull","M_{H}");
  // HzResolSVFitVSSumPtPions=HConfig.GetTH2D(Name+"_HzResolSVFitVSSumPtPions","h_{Z} pull VS pions Pt sum",150,-3,1,100,0,100,"h_{Z} pull","pions Pt sum");
  // HzResolSVFitVSSumPxPions=HConfig.GetTH2D(Name+"_HzResolSVFitVSSumPxPions","h_{Z} pull VS pions_{X} sum",150,-3,1,80,0.,80.,"h_{Z} pull","pions_{X} sum");
  // HzResolSVFitVSSumPyPions=HConfig.GetTH2D(Name+"_HzResolSVFitVSSumPyPions","h_{Z} pull VS pions_{Y} sum",150,-3,1,80,0.,80.,"h_{Z} pull","pions_{Y} sum");
  // HzResolSVFitVSSumPzPions=HConfig.GetTH2D(Name+"_HzResolSVFitVSSumPzPions","h_{Z} pull VS pions_{Z} sum",150,-3,1,80,0.,80.,"h_{Z} pull","pions_{Z} sum");
  // HzResolSVFitVSSumEPions=HConfig.GetTH2D(Name+"_HzResolSVFitVSSumEPions","h_{Z} pull VS pions energy sum",150,-3,1,150,0,300,"h_{Z} pull","pions energy sum");

  // HxResolRecoVSTauPt=HConfig.GetTH2D(Name+"_HxResolRecoVSTauPt","h_{X} pull VS tau Pt",150,-3,1,100,0.,100.,"h_{X} pull","tau Pt");
  // HxResolRecoVSTauEta=HConfig.GetTH2D(Name+"_HxResolRecoVSTauEta","h_{X} pull VS tau eta",150,-3,1,48,-2.4,2.4,"h_{X} pull","tau eta");
  // HxResolRecoVSPxTau=HConfig.GetTH2D(Name+"_HxResolRecoVSPxTau","h_{X} pull VS tau P_[X}",150,-3,1,150,-150,150,"h_{X} pull","tau P_[X}");
  // HxResolRecoVSPyTau=HConfig.GetTH2D(Name+"_HxResolRecoVSPyTau","h_{X} pull VS tau P_{Y}",150,-3,1,150,-150,150,"h_{X} pull","tau P_{Y}");
  // HxResolRecoVSPzTau=HConfig.GetTH2D(Name+"_HxResolRecoVSPzTau","h_{X} pull VS tau P_{Z}",150,-3,1,400,-400.,400.,"h_{X} pull","tau P_{Z}");
  // HxResolRecoVSETau=HConfig.GetTH2D(Name+"_HxResolRecoVSETau","h_{X} pull VS tau energy",150,-3,1,200,0,400.,"h_{X} pull","tau energy");
  // HxResolRecoVSGJAngle=HConfig.GetTH2D(Name+"_HxResolRecoVSGJAngle","h_{X} pull VS ",150,-3,1,100,0,0.03,"h_{X} pull","");
  // HxResolRecoVSMH=HConfig.GetTH2D(Name+"_HxResolRecoVSMH","h_{X} pull VS M_{H}",150,-3,1,100,120,130,"h_{X} pull","M_{H}");
  // HxResolRecoVSSumPtPions=HConfig.GetTH2D(Name+"_HxResolRecoVSSumPtPions","h_{X} pull VS pions Pt sum",150,-3,1,100,0,100,"h_{X} pull","pions Pt sum");
  // HxResolRecoVSSumPxPions=HConfig.GetTH2D(Name+"_HxResolRecoVSSumPxPions","h_{X} pull VS pions_{X} sum",150,-3,1,80,0.,80.,"h_{X} pull","pions_{X} sum");
  // HxResolRecoVSSumPyPions=HConfig.GetTH2D(Name+"_HxResolRecoVSSumPyPions","h_{X} pull VS pions_{Y} sum",150,-3,1,80,0.,80.,"h_{X} pull","pions_{Y} sum");
  // HxResolRecoVSSumPzPions=HConfig.GetTH2D(Name+"_HxResolRecoVSSumPzPions","h_{X} pull VS pions_{Z} sum",150,-3,1,80,0.,80.,"h_{X} pull","pions_{Z} sum");
  // HxResolRecoVSSumEPions=HConfig.GetTH2D(Name+"_HxResolRecoVSSumEPions","h_{X} pull VS pions energy sum",150,-3,1,150,0,300,"h_{X} pull","pions energy sum");

  // HyResolRecoVSTauPt=HConfig.GetTH2D(Name+"_HyResolRecoVSTauPt","h_{Y} pull VS tau Pt",150,-3,1,100,0.,100.,"h_{Y} pull","tau Pt");
  // HyResolRecoVSTauEta=HConfig.GetTH2D(Name+"_HyResolRecoVSTauEta","h_{Y} pull VS tau eta",150,-3,1,48,-2.4,2.4,"h_{Y} pull","tau eta");
  // HyResolRecoVSPxTau=HConfig.GetTH2D(Name+"_HyResolRecoVSPxTau","h_{Y} pull VS tau P_[X}",150,-3,1,150,-150,150,"h_{Y} pull","tau P_[X}");
  // HyResolRecoVSPyTau=HConfig.GetTH2D(Name+"_HyResolRecoVSPyTau","h_{Y} pull VS tau P_{Y}",150,-3,1,150,-150,150,"h_{Y} pull","tau P_{Y}");
  // HyResolRecoVSPzTau=HConfig.GetTH2D(Name+"_HyResolRecoVSPzTau","h_{Y} pull VS tau P_{Z}",150,-3,1,400,-400.,400.,"h_{Y} pull","tau P_{Z}");
  // HyResolRecoVSETau=HConfig.GetTH2D(Name+"_HyResolRecoVSETau","h_{Y} pull VS tau energy",150,-3,1,200,0,400.,"h_{Y} pull","tau energy");
  // HyResolRecoVSGJAngle=HConfig.GetTH2D(Name+"_HyResolRecoVSGJAngle","h_{Y} pull VS ",150,-3,1,100,0,0.03,"h_{Y} pull","");
  // HyResolRecoVSMH=HConfig.GetTH2D(Name+"_HyResolRecoVSMH","h_{Y} pull VS M_{H}",150,-3,1,100,120,130,"h_{Y} pull","M_{H}");
  // HyResolRecoVSSumPtPions=HConfig.GetTH2D(Name+"_HyResolRecoVSSumPtPions","h_{Y} pull VS pions Pt sum",150,-3,1,100,0,100,"h_{Y} pull","pions Pt sum");
  // HyResolRecoVSSumPxPions=HConfig.GetTH2D(Name+"_HyResolRecoVSSumPxPions","h_{Y} pull VS pions_{X} sum",150,-3,1,80,0.,80.,"h_{Y} pull","pions_{X} sum");
  // HyResolRecoVSSumPyPions=HConfig.GetTH2D(Name+"_HyResolRecoVSSumPyPions","h_{Y} pull VS pions_{Y} sum",150,-3,1,80,0.,80.,"h_{Y} pull","pions_{Y} sum");
  // HyResolRecoVSSumPzPions=HConfig.GetTH2D(Name+"_HyResolRecoVSSumPzPions","h_{Y} pull VS pions_{Z} sum",150,-3,1,80,0.,80.,"h_{Y} pull","pions_{Z} sum");
  // HyResolRecoVSSumEPions=HConfig.GetTH2D(Name+"_HyResolRecoVSSumEPions","h_{Y} pull VS pions energy sum",150,-3,1,150,0,300,"h_{Y} pull","pions energy sum");


  // HzResolRecoVSTauPt=HConfig.GetTH2D(Name+"_HzResolRecoVSTauPt","h_{Z} pull VS tau Pt",150,-3,1,100,0.,100.,"h_{Z} pull","tau Pt");
  // HzResolRecoVSTauEta=HConfig.GetTH2D(Name+"_HzResolRecoVSTauEta","h_{Z} pull VS tau eta",150,-3,1,48,-2.4,2.4,"h_{Z} pull","tau eta");
  // HzResolRecoVSPxTau=HConfig.GetTH2D(Name+"_HzResolRecoVSPxTau","h_{Z} pull VS tau P_[X}",150,-3,1,150,-150,150,"h_{Z} pull","tau P_[X}");
  // HzResolRecoVSPyTau=HConfig.GetTH2D(Name+"_HzResolRecoVSPyTau","h_{Z} pull VS tau P_{Y}",150,-3,1,150,-150,150,"h_{Z} pull","tau P_{Y}");
  // HzResolRecoVSPzTau=HConfig.GetTH2D(Name+"_HzResolRecoVSPzTau","h_{Z} pull VS tau P_{Z}",150,-3,1,400,-400.,400.,"h_{Z} pull","tau P_{Z}");
  // HzResolRecoVSETau=HConfig.GetTH2D(Name+"_HzResolRecoVSETau","h_{Z} pull VS tau energy",150,-3,1,200,0,400.,"h_{Z} pull","tau energy");
  // HzResolRecoVSGJAngle=HConfig.GetTH2D(Name+"_HzResolRecoVSGJAngle","h_{Z} pull VS ",150,-3,1,100,0,0.03,"h_{Z} pull","");
  // HzResolRecoVSMH=HConfig.GetTH2D(Name+"_HzResolRecoVSMH","h_{Z} pull VS M_{H}",150,-3,1,100,120,130,"h_{Z} pull","M_{H}");
  // HzResolRecoVSSumPtPions=HConfig.GetTH2D(Name+"_HzResolRecoVSSumPtPions","h_{Z} pull VS pions Pt sum",150,-3,1,100,0,100,"h_{Z} pull","pions Pt sum");
  // HzResolRecoVSSumPxPions=HConfig.GetTH2D(Name+"_HzResolRecoVSSumPxPions","h_{Z} pull VS pions_{X} sum",150,-3,1,80,0.,80.,"h_{Z} pull","pions_{X} sum");
  // HzResolRecoVSSumPyPions=HConfig.GetTH2D(Name+"_HzResolRecoVSSumPyPions","h_{Z} pull VS pions_{Y} sum",150,-3,1,80,0.,80.,"h_{Z} pull","pions_{Y} sum");
  // HzResolRecoVSSumPzPions=HConfig.GetTH2D(Name+"_HzResolRecoVSSumPzPions","h_{Z} pull VS pions_{Z} sum",150,-3,1,80,0.,80.,"h_{Z} pull","pions_{Z} sum");
  // HzResolRecoVSSumEPions=HConfig.GetTH2D(Name+"_HzResolRecoVSSumEPions","h_{Z} pull VS pions energy sum",150,-3,1,150,0,300,"h_{Z} pull","pions energy sum");

  
  Selection::ConfigureHistograms();   //   do not remove
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);  // do not remove
}

 

void  Debug::Store_ExtraDist(){

  //every new histo should be addedd to Extradist1d vector, just push it back;
  
  Extradist1d.push_back(&h_SVFitMass);
  
  Extradist1d.push_back(&polarimetricAcopAngle);
  Extradist1d.push_back(&polarimetricAcopAngleNoBS);
  Extradist1d.push_back(&polarimetricAcopAngleNoBSPionsRefit);
  Extradist1d.push_back(&polarimetricAcopAngleTruthA1);
  Extradist1d.push_back(&polarimetricAcopAngleSVFitA1);  
  Extradist1d.push_back(&polarimetricAcopAngleOld);
  Extradist1d.push_back(&polarimetricAcopAngleOldPionsRefit);
  Extradist1d.push_back(&polarimetricAcopAngleTIP);
  Extradist1d.push_back(&polarimetricAcopAngleTIPPionsRefit);
  
  

  Extradist1d.push_back(&PVXResol);
  Extradist1d.push_back(&PVXNoBSResol);
  Extradist1d.push_back(&PVXOldRefitResol);
  Extradist1d.push_back(&PVXTIPRefitResol);
 
  Extradist1d.push_back(&PVYResol);
  Extradist1d.push_back(&PVYNoBSResol);
  Extradist1d.push_back(&PVYOldRefitResol);
  Extradist1d.push_back(&PVYTIPRefitResol);
 
  Extradist1d.push_back(&PVZResol);
  Extradist1d.push_back(&PVZNoBSResol);
  Extradist1d.push_back(&PVZOldRefitResol);
  Extradist1d.push_back(&PVZTIPRefitResol);
   
  Extradist1d.push_back(&TauSVFitPxResPull);
  Extradist1d.push_back(&TauSVFitPyResPull);
  Extradist1d.push_back(&TauSVFitPzResPull);
  Extradist1d.push_back(&TauSVFitEResPull);
 
  Extradist1d.push_back(&TauPxResPull);
  Extradist1d.push_back(&TauPyResPull);
  Extradist1d.push_back(&TauPzResPull);
  Extradist1d.push_back(&TauEResPull);

  Extradist1d.push_back(&TauNoBSPxResPull);
  Extradist1d.push_back(&TauNoBSPyResPull);
  Extradist1d.push_back(&TauNoBSPzResPull);
  Extradist1d.push_back(&TauNoBSEResPull);
  
  Extradist1d.push_back(&TauOldPxResPull);
  Extradist1d.push_back(&TauOldPyResPull);
  Extradist1d.push_back(&TauOldPzResPull);
  Extradist1d.push_back(&TauOldEResPull);

  Extradist1d.push_back(&TauTIPPxResPull);
  Extradist1d.push_back(&TauTIPPyResPull);
  Extradist1d.push_back(&TauTIPPzResPull);
  Extradist1d.push_back(&TauTIPEResPull);

  Extradist1d.push_back(&TauNoBSRefittedPionsPxResPull);
  Extradist1d.push_back(&TauNoBSRefittedPionsPyResPull);
  Extradist1d.push_back(&TauNoBSRefittedPionsPzResPull);
  Extradist1d.push_back(&TauNoBSRefittedPionsEResPull);
  
  Extradist1d.push_back(&TauOldRefittedPionsPxResPull);
  Extradist1d.push_back(&TauOldRefittedPionsPyResPull);
  Extradist1d.push_back(&TauOldRefittedPionsPzResPull);
  Extradist1d.push_back(&TauOldRefittedPionsEResPull);

  Extradist1d.push_back(&TauTIPRefittedPionsPxResPull);
  Extradist1d.push_back(&TauTIPRefittedPionsPyResPull);
  Extradist1d.push_back(&TauTIPRefittedPionsPzResPull);
  Extradist1d.push_back(&TauTIPRefittedPionsEResPull);
  
  Extradist1d.push_back(&LeadingPionPxResPull);
  Extradist1d.push_back(&LeadingPionPyResPull);
  Extradist1d.push_back(&LeadingPionPzResPull);
  Extradist1d.push_back(&LeadingPionEResPull);

  Extradist1d.push_back(&OppositePionPxResPull);
  Extradist1d.push_back(&OppositePionPyResPull);
  Extradist1d.push_back(&OppositePionPzResPull);
  Extradist1d.push_back(&OppositePionEResPull);

  Extradist1d.push_back(&LeadingRefittedPionPxResPull);
  Extradist1d.push_back(&LeadingRefittedPionPyResPull);
  Extradist1d.push_back(&LeadingRefittedPionPzResPull);
  Extradist1d.push_back(&LeadingRefittedPionEResPull);

  Extradist1d.push_back(&OppositeRefittedPionPxResPull);
  Extradist1d.push_back(&OppositeRefittedPionPyResPull);
  Extradist1d.push_back(&OppositeRefittedPionPzResPull);
  Extradist1d.push_back(&OppositeRefittedPionEResPull);

  Extradist1d.push_back(&VectPolaSVFitPxResPull);
  Extradist1d.push_back(&VectPolaSVFitPyResPull);
  Extradist1d.push_back(&VectPolaSVFitPzResPull);

  Extradist1d.push_back(&VectPolaPxResPull);
  Extradist1d.push_back(&VectPolaPyResPull);
  Extradist1d.push_back(&VectPolaPzResPull);

  // Extradist1d.push_back(&DiffTauXLeadingPionXSVFit);
  // Extradist1d.push_back(&DiffTauYLeadingPionYSVFit);
  // Extradist1d.push_back(&DiffTauZLeadingPionZSVFit);
  // Extradist1d.push_back(&DiffTauELeadingPionESVFit);

  // Extradist1d.push_back(&DiffTauXLeadingPionXReco);
  // Extradist1d.push_back(&DiffTauYLeadingPionYReco);
  // Extradist1d.push_back(&DiffTauZLeadingPionZReco);
  // Extradist1d.push_back(&DiffTauELeadingPionEReco);

  //Extradist2d.push_back(&AcopAngleVSTauPt);
  //Extradist2d.push_back(&AcopAngleVSTauEta);

  // Extradist2d.push_back(&PxtauSVFitVSPxLeadingPion);
  // Extradist2d.push_back(&PytauSVFitVSPyLeadingPion);
  // Extradist2d.push_back(&PztauSVFitVSPzLeadingPion);
  // Extradist2d.push_back(&EtauSVFitVSELeadingPion);

  // Extradist2d.push_back(&HxResolSVFitVSDiffPxtauSVFitPxLeadingPion);
  // Extradist2d.push_back(&HxResolSVFitVSDiffPytauSVFitPyLeadingPion);
  // Extradist2d.push_back(&HxResolSVFitVSDiffPztauSVFitPzLeadingPion);
  // Extradist2d.push_back(&HxResolSVFitVSDiffEtauSVFitELeadingPion);
  
  // Extradist2d.push_back(&HyResolSVFitVSDiffPxtauSVFitPxLeadingPion);
  // Extradist2d.push_back(&HyResolSVFitVSDiffPytauSVFitPyLeadingPion);
  // Extradist2d.push_back(&HyResolSVFitVSDiffPztauSVFitPzLeadingPion);
  // Extradist2d.push_back(&HyResolSVFitVSDiffEtauSVFitELeadingPion);

  // Extradist2d.push_back(&HzResolSVFitVSDiffPxtauSVFitPxLeadingPion);
  // Extradist2d.push_back(&HzResolSVFitVSDiffPytauSVFitPyLeadingPion);
  // Extradist2d.push_back(&HzResolSVFitVSDiffPztauSVFitPzLeadingPion);
  // Extradist2d.push_back(&HzResolSVFitVSDiffEtauSVFitELeadingPion);

  
  // Extradist2d.push_back(&PxtauSVFitVSPxOppositePion);
  // Extradist2d.push_back(&PytauSVFitVSPyOppositePion);
  // Extradist2d.push_back(&PztauSVFitVSPzOppositePion);
  // Extradist2d.push_back(&EtauSVFitVSEOppositePion);

  // Extradist2d.push_back(&HxResolSVFitVSDiffPxtauSVFitPxOppositePion);
  // Extradist2d.push_back(&HxResolSVFitVSDiffPytauSVFitPyOppositePion);
  // Extradist2d.push_back(&HxResolSVFitVSDiffPztauSVFitPzOppositePion);
  // Extradist2d.push_back(&HxResolSVFitVSDiffEtauSVFitEOppositePion);

  
  // Extradist2d.push_back(&PxtauRecoVSPxLeadingPion);
  // Extradist2d.push_back(&PytauRecoVSPyLeadingPion);
  // Extradist2d.push_back(&PztauRecoVSPzLeadingPion);
  // Extradist2d.push_back(&EtauRecoVSELeadingPion);

  // Extradist2d.push_back(&HxResolRecoVSDiffPxtauRecoPxLeadingPion);
  // Extradist2d.push_back(&HxResolRecoVSDiffPytauRecoPyLeadingPion);
  // Extradist2d.push_back(&HxResolRecoVSDiffPztauRecoPzLeadingPion);
  // Extradist2d.push_back(&HxResolRecoVSDiffEtauRecoELeadingPion);
  
  // Extradist2d.push_back(&HyResolRecoVSDiffPxtauRecoPxLeadingPion);
  // Extradist2d.push_back(&HyResolRecoVSDiffPytauRecoPyLeadingPion);
  // Extradist2d.push_back(&HyResolRecoVSDiffPztauRecoPzLeadingPion);
  // Extradist2d.push_back(&HyResolRecoVSDiffEtauRecoELeadingPion);

  // Extradist2d.push_back(&HzResolRecoVSDiffPxtauRecoPxLeadingPion);
  // Extradist2d.push_back(&HzResolRecoVSDiffPytauRecoPyLeadingPion);
  // Extradist2d.push_back(&HzResolRecoVSDiffPztauRecoPzLeadingPion);
  // Extradist2d.push_back(&HzResolRecoVSDiffEtauRecoELeadingPion);
  
  // Extradist2d.push_back(&PxtauRecoVSPxOppositePion);
  // Extradist2d.push_back(&PytauRecoVSPyOppositePion);
  // Extradist2d.push_back(&PztauRecoVSPzOppositePion);
  // Extradist2d.push_back(&EtauRecoVSEOppositePion);

  // Extradist2d.push_back(&HxResolRecoVSDiffPxtauRecoPxOppositePion);
  // Extradist2d.push_back(&HxResolRecoVSDiffPytauRecoPyOppositePion);
  // Extradist2d.push_back(&HxResolRecoVSDiffPztauRecoPzOppositePion);
  // Extradist2d.push_back(&HxResolRecoVSDiffEtauRecoEOppositePion);


  // Extradist2d.push_back(&HxResolSVFitVSTauPt);
  // Extradist2d.push_back(&HxResolSVFitVSPxTau);
  // Extradist2d.push_back(&HxResolSVFitVSPyTau);
  // Extradist2d.push_back(&HxResolSVFitVSPzTau);
  // Extradist2d.push_back(&HxResolSVFitVSETau);
  // Extradist2d.push_back(&HxResolSVFitVSTauEta);
  // Extradist2d.push_back(&HxResolSVFitVSMH);
  // Extradist2d.push_back(&HxResolSVFitVSSumPtPions);
  // Extradist2d.push_back(&HxResolSVFitVSSumPxPions);
  // Extradist2d.push_back(&HxResolSVFitVSSumPyPions);
  // Extradist2d.push_back(&HxResolSVFitVSSumPzPions);
  // Extradist2d.push_back(&HxResolSVFitVSSumEPions);

  // Extradist2d.push_back(&HxResolRecoVSTauPt);
  // Extradist2d.push_back(&HxResolRecoVSPxTau);
  // Extradist2d.push_back(&HxResolRecoVSPyTau);
  // Extradist2d.push_back(&HxResolRecoVSPzTau);
  // Extradist2d.push_back(&HxResolRecoVSETau);
  // Extradist2d.push_back(&HxResolRecoVSSumEPions);
  // Extradist2d.push_back(&HxResolRecoVSTauEta);
  // Extradist2d.push_back(&HxResolRecoVSMH);
  // Extradist2d.push_back(&HxResolRecoVSSumPtPions);
  // Extradist2d.push_back(&HxResolRecoVSSumPxPions);
  // Extradist2d.push_back(&HxResolRecoVSSumPyPions);
  // Extradist2d.push_back(&HxResolRecoVSSumPzPions);
  // Extradist2d.push_back(&HxResolRecoVSSumEPions);
  
  // Extradist2d.push_back(&HyResolSVFitVSTauPt);
  // Extradist2d.push_back(&HyResolSVFitVSPxTau);
  // Extradist2d.push_back(&HyResolSVFitVSPyTau);
  // Extradist2d.push_back(&HyResolSVFitVSPzTau);
  // Extradist2d.push_back(&HyResolSVFitVSETau);
  // Extradist2d.push_back(&HyResolSVFitVSSumEPions);
  // Extradist2d.push_back(&HyResolSVFitVSTauEta);
  // Extradist2d.push_back(&HyResolSVFitVSMH);
  // Extradist2d.push_back(&HyResolSVFitVSSumPtPions);
  // Extradist2d.push_back(&HyResolSVFitVSSumPxPions);
  // Extradist2d.push_back(&HyResolSVFitVSSumPyPions);
  // Extradist2d.push_back(&HyResolSVFitVSSumPzPions);
  // Extradist2d.push_back(&HyResolSVFitVSSumEPions);

  // Extradist2d.push_back(&HyResolRecoVSTauPt);
  // Extradist2d.push_back(&HyResolRecoVSPxTau);
  // Extradist2d.push_back(&HyResolRecoVSPyTau);
  // Extradist2d.push_back(&HyResolRecoVSPzTau);
  // Extradist2d.push_back(&HyResolRecoVSETau);
  // Extradist2d.push_back(&HyResolRecoVSSumEPions);
  // Extradist2d.push_back(&HyResolRecoVSTauEta);
  // Extradist2d.push_back(&HyResolRecoVSMH);
  // Extradist2d.push_back(&HyResolRecoVSSumPtPions);
  // Extradist2d.push_back(&HyResolRecoVSSumPxPions);
  // Extradist2d.push_back(&HyResolRecoVSSumPyPions);
  // Extradist2d.push_back(&HyResolRecoVSSumPzPions);
  // Extradist2d.push_back(&HyResolRecoVSSumEPions);
  
  // Extradist2d.push_back(&HzResolSVFitVSTauPt);
  // Extradist2d.push_back(&HzResolSVFitVSPxTau);
  // Extradist2d.push_back(&HzResolSVFitVSPyTau);
  // Extradist2d.push_back(&HzResolSVFitVSPzTau);
  // Extradist2d.push_back(&HzResolSVFitVSETau);
  // Extradist2d.push_back(&HzResolSVFitVSSumEPions);
  // Extradist2d.push_back(&HzResolSVFitVSTauEta);
  // Extradist2d.push_back(&HzResolSVFitVSMH);
  // Extradist2d.push_back(&HzResolSVFitVSSumPtPions);
  // Extradist2d.push_back(&HzResolSVFitVSSumPxPions);
  // Extradist2d.push_back(&HzResolSVFitVSSumPyPions);
  // Extradist2d.push_back(&HzResolSVFitVSSumPzPions);
  // Extradist2d.push_back(&HzResolSVFitVSSumEPions);

  // Extradist2d.push_back(&HzResolRecoVSTauPt);
  // Extradist2d.push_back(&HzResolRecoVSPxTau);
  // Extradist2d.push_back(&HzResolRecoVSPyTau);
  // Extradist2d.push_back(&HzResolRecoVSPzTau);
  // Extradist2d.push_back(&HzResolRecoVSETau);
  // Extradist2d.push_back(&HzResolRecoVSSumEPions);
  // Extradist2d.push_back(&HzResolRecoVSTauEta);
  // Extradist2d.push_back(&HzResolRecoVSMH);
  // Extradist2d.push_back(&HzResolRecoVSSumPtPions);
  // Extradist2d.push_back(&HzResolRecoVSSumPxPions);
  // Extradist2d.push_back(&HzResolRecoVSSumPyPions);
  // Extradist2d.push_back(&HzResolRecoVSSumPzPions);
  // Extradist2d.push_back(&HzResolRecoVSSumEPions);

  //Extradist2d.push_back(&HxResolSVFitVSGJAngle);
  //Extradist2d.push_back(&HyResolSVFitVSGJAngle);
  //Extradist2d.push_back(&HzResolSVFitVSGJAngle);
  //Extradist2d.push_back(&HxResolRecoVSGJAngle);
  //Extradist2d.push_back(&HyResolRecoVSGJAngle);
  //Extradist2d.push_back(&HzResolRecoVSGJAngle);

  Extradist1d.push_back(&Unphysical);
}

void  Debug::doEvent()  { //  Method called on every event

  unsigned int t;                // sample type, you may manage in your further analysis, if needed
  int id(Ntp->GetMCID());  //read event ID of a sample
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}  //  gives a warning if list of samples in Histo.txt  and SkimSummary.log do not coincide 

  value.at(Id_and_Kin)=0;
  int goodTau_counter=0;
  std::vector<int> thirdLeptonCounter;
  std::vector<int> goodTausIndex;
  for(unsigned int iDaughter=0;   iDaughter  <Ntp->NDaughters();iDaughter++ ) {
    if(Ntp->tauBaselineSelection(iDaughter,20., 2.4, 0,0)){
      goodTausIndex.push_back(iDaughter);
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

  for(int unsigned ipair =0; ipair <goodTausIndex.size(); ipair++)
    {
      for(int unsigned jpair =1; jpair <goodTausIndex.size(); jpair++)
	{
	  if(jpair>ipair)
	    {
	      PairsIndexTemp.push_back(j);
	      PairsIndexTau1Temp.push_back(goodTausIndex.at(ipair));
	      PairsIndexTau2Temp.push_back(goodTausIndex.at(jpair));
	      j++;
	    }
	}
    }
  if(PairsIndexTemp.size()>0)
    {
      Sorted = Ntp->SortPair(PairsIndexTemp,PairsIndexTau1Temp,PairsIndexTau2Temp);
      Tau1=PairsIndexTau1Temp.at(Sorted.back());
      Tau2=PairsIndexTau2Temp.at(Sorted.back());

      value.at(PairCharge)=0;
      bool isOS=false;
      isOS=((Ntp->Daughters_charge(Tau1)/abs(Ntp->Daughters_charge(Tau1))) != (Ntp->Daughters_charge(Tau2)/abs(Ntp->Daughters_charge(Tau2))));
      if(isOS)value.at(PairCharge) = 1;
      pass.at(PairCharge) = value.at(PairCharge);
    }

  // Here you can defined different type of weights you want to apply to events.
  double wobs=1;
  double w=1;
  if(!Ntp->isData() && id!=DataMCType::QCD) {
    w *= reweight.PUweightHTT(Ntp->npu());
    if(!Ntp->isData() &&PairsIndexTemp.size()>0){
      double w1 = tauTrgSF.getSF(Ntp->TauP4_Corrected(Tau1).Pt(),  Ntp->decayMode(Tau1)) ;  //from Luca
      double w2 = tauTrgSF.getSF(Ntp->TauP4_Corrected(Tau2).Pt(),  Ntp->decayMode(Tau2)) ;
      w*=w1;
      w*=w2;
    }
    if(!Ntp->isData() && PairsIndexTemp.size()>0 && (id==33 || id == 10110333 || id == 10110433|| id == 10130533|| id ==10210333|| id == 10210433|| id == 10230533|| id ==10310333 || id ==10330533 || id ==10410433 || id == 10410333|| id == 10430533|| id == 30530533 || id ==11 || id ==12)){
      w *= 0.95*0.95;
    }
  }

  TLorentzVector genMomentum(0,0,0,0);
  if(id==33 || id == 10110333 || id == 10110433|| id == 10130533|| id ==10210333|| id == 10210433|| id == 10230533|| id ==10310333 || id ==10330533 || id ==10410433 || id == 10410333|| id == 10430533|| id == 30530533){
    for(unsigned int imc=0; imc < Ntp->NGenParts(); imc++){
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
    if(PairsIndexTemp.size()>0){
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
  w*=wAgainstMuon1;
  w*=wAgainstElectron1;
  w*=wAgainstMuon2;
  w*=wAgainstElectron2;
  if(!Ntp->isData() && id!=DataMCType::QCD)w*=Ntp->MC_weight(); //generator weight because negative weights for this samples
  
  if(id==11)w*=0.2514;
  if(id==12)w*=0.2825;
  w*=Ntp->stitch_weight();
  
  bool status=AnalysisCuts(t,w,wobs);  // boolean that say whether your event passed critera defined in pass vector. The whole vector must be true for status = true
  ///////////////////////////////////////////////////////////
  // Analyse events which passed selection
  if(status) {

    TLorentzVector Tau1P4;
    TLorentzVector Tau2P4;
    
    Tau1P4 = Ntp->TauP4_Corrected(Tau1);
    Tau2P4 = Ntp->TauP4_Corrected(Tau2);
    
    unsigned int Tauplus=0;
    unsigned int Tauminus=0;
  
    if(Ntp->Daughters_charge(Tau1)>0)
      {
    	Tauplus=Tau1;
    	Tauminus=Tau2;
      }
    else
      {
    	Tauplus=Tau2;
    	Tauminus=Tau1;
      }
    
    TVector3 tauPrimaryVertex, TauminusSecondaryVertex, TauplusSecondaryVertex, tauPrimaryVertexNoBS, tauPrimaryVertexOld, tauPrimaryVertexTIP;
    TVector3 TauminusDirection , TauplusDirection, TauminusDirectionNoBS, TauplusDirectionNoBS, TauminusDirectionOld, TauplusDirectionOld, TauminusDirectionTIP, TauplusDirectionTIP;
    TLorentzVector a1LV_Tauminus , a1LV_Tauplus,a1LVNoBS_Tauminus , a1LVNoBS_Tauplus, a1LVRefittedPions_Tauminus , a1LVRefittedPions_Tauplus;
    TLorentzVector TauminusPairConstraint, TauplusPairConstraint, TauminusPairConstraintNoBS,TauplusPairConstraintNoBS, TauminusPairConstraintOld,TauplusPairConstraintOld, TauminusPairConstraintTIP,TauplusPairConstraintTIP, TauminusPairConstraintNoBSRefittedPions,TauplusPairConstraintNoBSRefittedPions, TauminusPairConstraintOldRefittedPions,TauplusPairConstraintOldRefittedPions, TauminusPairConstraintTIPRefittedPions,TauplusPairConstraintTIPRefittedPions;
    bool isPlusReal=true, isMinusReal=true, a1a1=false,  a1minus=false, a1plus=false, piminus=false, piplus=false, a1pi=false;
    std::vector<TLorentzVector> solutions, solutionsNoBS, solutionsNoBSRefittedPions, solutionsOld, solutionsOldRefittedPions, solutionsTIP, solutionsTIPRefittedPions;

    std::vector<size_t> hashes;
    size_t hash = 0;

    bool hasNoBS=false;

    if(Ntp->decayMode(Tauminus) == 10 &&  Ntp->PFTau_hassecondaryVertex(Tauminus) && Ntp->PFtauHasThreePions(Tauminus))a1minus=true;
    if(Ntp->decayMode(Tauplus) == 10 && Ntp->PFTau_hassecondaryVertex(Tauplus) && Ntp->PFtauHasThreePions(Tauplus))a1plus=true;
    if(Ntp->decayMode(Tauminus) == 0 )piminus=true;
    if(Ntp->decayMode(Tauplus) == 0 )piplus=true;


    if(a1plus && a1minus) a1a1=true;
    if((a1plus && piminus) ||(piplus && a1minus)) a1pi=true;

    double Spin_WT=Ntp->TauSpinerGet(TauSpinerInterface::Spin);	 
    double Wspin=w*Spin_WT; 
    
    classic_svFit::LorentzVector tau1P4;
    classic_svFit::LorentzVector tau2P4;
  
    TLorentzVector Tauplussvfit;
    TLorentzVector Tauminussvfit;
    
    ClassicSVfit svfitAlgo1;
    FastMTT FastMTTAlgo;
    double higgsmass;

    if(a1a1 && std::isnan(Wspin)!=true)
      {
    	// // //---------  svfit ---------------------
    	std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;
    	classic_svFit::MeasuredTauLepton lep1(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Tau1P4.Pt(), Tau1P4.Eta(),  Tau1P4.Phi(), Tau1P4.M(),10);
    	classic_svFit::MeasuredTauLepton lep2(classic_svFit::MeasuredTauLepton::kTauToHadDecay, Tau2P4.Pt(), Tau2P4.Eta(),  Tau2P4.Phi(), Tau2P4.M(),10);
	
    	measuredTauLeptons.push_back(lep1);
    	measuredTauLeptons.push_back(lep2);
    	TMatrixD metcov(2,2);
    	//double metx = Ntp->PUPPImet()*cos(Ntp->PUPPImetphi());
    	//double mety = Ntp->PUPPImet()*sin(Ntp->PUPPImetphi());
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
	  
	  //h_SVFitMass.at(t).Fill(FastMTTAlgo.getBestP4().M(),w);
    	higgsmass  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svfitAlgo1.getHistogramAdapter())->getMass();
    	h_SVFitMass.at(t).Fill(higgsmass,w); 
	  
    	//tau1P4 = FastMTTAlgo.getTau1P4();
    	//tau2P4 = FastMTTAlgo.getTau2P4();

    	tau1P4 = static_cast<classic_svFit::TauTauHistogramAdapter*>(svfitAlgo1.getHistogramAdapter())->GetFittedTau1LV();
    	tau2P4 = static_cast<classic_svFit::TauTauHistogramAdapter*>(svfitAlgo1.getHistogramAdapter())->GetFittedTau2LV();
	
    	// // //---------  svfit ---------------------
	  
    	//classic_svFit::TauTauHistogramAdapter* HistogramAdapter = static_cast<classic_svFit::TauTauHistogramAdapter*>(svfitAlgo1.getHistogramAdapter());
	  
    	if(Ntp->Daughters_charge(Tau1)>0)
    	  {
    	    //Tauplussvfit.SetPtEtaPhiM(HistogramAdapter->getPt(),HistogramAdapter->getEta(),HistogramAdapter->getPhi(), HistogramAdapter->getMass());
    	    //Tauminussvfit.SetPtEtaPhiM(HistogramAdapter->getPt(),HistogramAdapter->getEta(),HistogramAdapter->getPhi(), HistogramAdapter->getMass());
	      
	      
    	    Tauplussvfit.SetPxPyPzE(tau1P4.Px(),tau1P4.Py(),tau1P4.Pz(),sqrt( pow(1.77685, 2.0) + pow(tau1P4.P(), 2.0) )/*tau1P4.E()*/);
    	    Tauminussvfit.SetPxPyPzE(tau2P4.Px(),tau2P4.Py(),tau2P4.Pz(),sqrt( pow(1.77685, 2.0) + pow(tau2P4.P(), 2.0) )/*tau2P4.E()*/);
	      
    	    // cout<<"Vector SVFit"<<endl;
    	    // Tau1P4.Print();
    	    // Tauplussvfit.Print();
    	    // Tau2P4.Print();
    	    // Tauminussvfit.Print();
    	  }
    	else
    	  {
	      
    	    Tauplussvfit.SetPxPyPzE(tau2P4.Px(),tau2P4.Py(),tau2P4.Pz(),sqrt( pow(1.77685, 2.0) + pow(tau2P4.P(), 2.0) )/*tau2P4.E()*/);
    	    Tauminussvfit.SetPxPyPzE(tau1P4.Px(),tau1P4.Py(),tau1P4.Pz(),sqrt( pow(1.77685, 2.0) + pow(tau1P4.P(), 2.0) )/*tau1P4.E()*/);
    	    // cout<<"Vector SVFit"<<endl;
    	    // Tau1P4.Print();
    	    // Tauminussvfit.Print();
    	    // Tau2P4.Print();
    	    // Tauplussvfit.Print();
    	  }
    	//Tauplussvfit.Print();
    	//cout<<"M: "<<tau1P4.E()*tau1P4.E()-tau1P4.Px()*tau1P4.Px()-tau1P4.Py()*tau1P4.Py()-tau1P4.Pz()*tau1P4.Pz()<<endl;
    	}
	
	
    	TLorentzVector zeroLV(0,0,0,0);
    	std::vector<TLorentzVector> VectZeroLV;
    	VectZeroLV.push_back(zeroLV);
    	VectZeroLV.push_back(zeroLV);
    	VectZeroLV.push_back(zeroLV);
	  
    	TLorentzVector Tau1Truth; 
    	TLorentzVector Tau2Truth;
    	TLorentzVector TruthDecayFromTau1;
    	TLorentzVector TruthDecayFromTau2; 
	
    	std::vector<TLorentzVector> Pions1Truth;
    	std::vector<TLorentzVector> Pions2Truth;
    	std::vector<double> Pions1TruthCharge;
    	std::vector<double> Pions2TruthCharge;
    	vector<TLorentzVector> HadPions_plus;    
    	vector<double> HadPionsCharge_plus;
    	vector<TLorentzVector> HadPions_minus; 
    	vector<double> HadPionsCharge_minus;
    	vector<TLorentzVector> HadRefittedPions_plus;    
    	vector<double> HadRefittedPionsCharge_plus;
    	vector<TLorentzVector> HadRefittedPions_minus; 
    	vector<double> HadRefittedPionsCharge_minus;
    	vector<double> HadPionsNoBSCharge_minus;
    	vector<TLorentzVector> HadPionsNoBS_plus;    
    	vector<double> HadPionsNoBSCharge_plus;
    	vector<TLorentzVector> HadPionsNoBS_minus;
    	vector<double> HadRefittedPionsNoBSCharge_minus;
    	vector<TLorentzVector> HadRefittedPionsNoBS_plus;    
    	vector<double> HadRefittedPionsNoBSCharge_plus;
    	vector<TLorentzVector> HadRefittedPionsNoBS_minus;
	
    	std::vector<TLorentzVector> PionsSVFit1;
    	std::vector<TLorentzVector> PionsSVFit2;
    	std::vector<double> PionsSVFit1Charge;
    	std::vector<double> PionsSVFit2Charge;
    
    	vector<TLorentzVector> tauandprodTruthminus;
    	vector<TLorentzVector> tauandprodTruthplus;
    	vector<TLorentzVector> tauandprodminus;
    	vector<TLorentzVector> tauandprodplus;
    	vector<TLorentzVector> tauandprodNoBSminus;
    	vector<TLorentzVector> tauandprodNoBSplus;
    	vector<TLorentzVector> tauandprodOldminus;
    	vector<TLorentzVector> tauandprodOldplus;
    	vector<TLorentzVector> tauandprodTIPminus;
    	vector<TLorentzVector> tauandprodTIPplus;
    	vector<TLorentzVector> tauandprodNoBSRefittedPionsminus;
    	vector<TLorentzVector> tauandprodNoBSRefittedPionsplus;
    	vector<TLorentzVector> tauandprodOldRefittedPionsminus;
    	vector<TLorentzVector> tauandprodOldRefittedPionsplus;
    	vector<TLorentzVector> tauandprodTIPRefittedPionsminus;
    	vector<TLorentzVector> tauandprodTIPRefittedPionsplus;
    	vector<TLorentzVector> tauandprodSVFitminus;
    	vector<TLorentzVector> tauandprodSVFitplus;

    	bool passParticlesTruth=false;
    	bool passParticlesReco=false;
    	bool passParticlesOld=false;
    	bool passParticlesTIP=false;
    	bool passParticlesNoBS=false;
    	bool passParticlesOldRefittedPions=false;
    	bool passParticlesTIPRefittedPions=false;
    	bool passParticlesNoBSRefittedPions=false;
    	bool passParticlesSVFit=false;
	bool unphysicalBool=false;

    	bool passVariablesTruth=false;
    	bool passVariablesReco=false;
    	bool passVariablesOld=false;
    	bool passVariablesTIP=false;
    	bool passVariablesNoBS=false;
    	bool passVariablesOldRefittedPions=false;
    	bool passVariablesTIPRefittedPions=false;
    	bool passVariablesNoBSRefittedPions=false;
    	bool passVariablesSVFit=false;

    	// Truth
	  
    	// for(int i=0;i<Ntp->NMCTauDecayProducts(0);i++)
    	//   {
    	//     cout<<"pdgid0: "<<Ntp->MCTauandProd_pdgid(0,i);
    	//   }
    	// cout<<endl;
    	// for(int i=0;i<Ntp->NMCTauDecayProducts(1);i++)
    	//   {
    	//     cout<<"pdgid1: "<<Ntp->MCTauandProd_pdgid(1,i);
    	//   }
    	// cout<<endl;
	
    	Tau1Truth=Ntp->GetTruthTauLV(5,0);
    	Tau2Truth=Ntp->GetTruthTauLV(5,1);
      
    	Pions1Truth=Ntp->GetTruthPionsFromA1(0);
    	TruthDecayFromTau1=Pions1Truth.at(0)+Pions1Truth.at(1)+Pions1Truth.at(2);
    	Pions1TruthCharge.push_back(1);
    	Pions1TruthCharge.push_back(-1);
    	Pions1TruthCharge.push_back(-1);
      
    	Pions2Truth=Ntp->GetTruthPionsFromA1(1);
    	TruthDecayFromTau2=Pions2Truth.at(0)+Pions2Truth.at(1)+Pions2Truth.at(2);
    	Pions2TruthCharge.push_back(-1);
    	Pions2TruthCharge.push_back(1);
    	Pions2TruthCharge.push_back(1);
	
    	SCalculator Scalc1Truth("a1");
    	SCalculator Scalc2Truth("a1");
	
    	Scalc1Truth.SortPions(Pions1Truth, Pions1TruthCharge);
    	Scalc2Truth.SortPions(Pions2Truth, Pions2TruthCharge);
			  
    	tauandprodTruthminus.push_back(Ntp->GetTruthTauLV(5,0));
    	tauandprodTruthminus.push_back(Pions1Truth.at(0));
    	tauandprodTruthminus.push_back(Pions1Truth.at(1));
    	tauandprodTruthminus.push_back(Pions1Truth.at(2));
    	tauandprodTruthplus.push_back(Ntp->GetTruthTauLV(5,1));   
    	tauandprodTruthplus.push_back(Pions2Truth.at(0));   
    	tauandprodTruthplus.push_back(Pions2Truth.at(1));   
    	tauandprodTruthplus.push_back(Pions2Truth.at(2));   
	    
    	if((Pions1Truth!=Pions2Truth) && (Pions1Truth!=VectZeroLV) && (Pions2Truth!=VectZeroLV) && tauandprodTruthminus.at(0)!=zeroLV && tauandprodTruthplus.at(0)!=zeroLV && tauandprodTruthminus.at(0)!=tauandprodTruthplus.at(0)) passParticlesTruth=true;
	 
    
    	// Reco NoBS
    	// create hashes from lepton selection
    	if (Ntp->isPVtxRefit()) {
    	  boost::hash_combine(hash, Ntp->LeptonHash(Tau1));
     	  boost::hash_combine(hash, Ntp->LeptonHash(Tau2));

    	  hashes.push_back(hash);
    	  hash = 0;

    	  boost::hash_combine(hash, Ntp->LeptonHash(Tau2));
    	  boost::hash_combine(hash, Ntp->LeptonHash(Tau1));
    	  hashes.push_back(hash);
    	}

    	if (Ntp->isPVtxRefit()) {
    	  // find the vertex among the refitted vertices
    	  for (unsigned int ivertex =0;ivertex<Ntp->NRefitVtx();ivertex++){
    	    size_t selectionHash = 0;
    	    boost::hash_combine(selectionHash, Ntp->VertexHashNoBS1(ivertex));
    	    boost::hash_combine(selectionHash, Ntp->VertexHashNoBS2(ivertex));
    	    if ( std::find(hashes.begin(), hashes.end(), selectionHash) != hashes.end() ){
    	      tauPrimaryVertexNoBS = Ntp->PVRefitNoBS(ivertex);
    	      hasNoBS=true;
    	      break;
    	    }
    	  } // loop over refitted vertices collection
	  
    	} // loop over refitted vertices collection

    	SCalculator Scalc1NoBS("a1");
    	SCalculator Scalc2NoBS("a1");
    	TauminusSecondaryVertex = Ntp->PFTau_secondaryVertex_pos(Tauminus);
    	TauplusSecondaryVertex = Ntp->PFTau_secondaryVertex_pos(Tauplus);
    	if(hasNoBS)
    	  {
    	    HadPionsNoBS_minus.push_back(Ntp->PFTau_PionsP4(Tauminus,0));
    	    HadPionsNoBS_minus.push_back(Ntp->PFTau_PionsP4(Tauminus,1));
    	    HadPionsNoBS_minus.push_back(Ntp->PFTau_PionsP4(Tauminus,2));

    	    HadPionsNoBSCharge_minus.push_back(Ntp->PFTau_PionsCharge(Tauminus, 0));
    	    HadPionsNoBSCharge_minus.push_back(Ntp->PFTau_PionsCharge(Tauminus, 1));
    	    HadPionsNoBSCharge_minus.push_back(Ntp->PFTau_PionsCharge(Tauminus, 2));

	
    	    HadPionsNoBS_plus.push_back(Ntp->PFTau_PionsP4(Tauplus,0));
    	    HadPionsNoBS_plus.push_back(Ntp->PFTau_PionsP4(Tauplus,1));
    	    HadPionsNoBS_plus.push_back(Ntp->PFTau_PionsP4(Tauplus,2));
	
    	    HadPionsNoBSCharge_plus.push_back(Ntp->PFTau_PionsCharge(Tauplus, 0));
    	    HadPionsNoBSCharge_plus.push_back(Ntp->PFTau_PionsCharge(Tauplus, 1));
    	    HadPionsNoBSCharge_plus.push_back(Ntp->PFTau_PionsCharge(Tauplus, 2));
	  
    	    TauminusDirectionNoBS = TauminusSecondaryVertex - tauPrimaryVertexNoBS;
    	    TauplusDirectionNoBS = TauplusSecondaryVertex - tauPrimaryVertexNoBS;
		  
    	    a1LVNoBS_Tauminus = Ntp->PFTau_PionsP4(Tauminus,0) + Ntp->PFTau_PionsP4(Tauminus,1) + Ntp->PFTau_PionsP4(Tauminus,2);	
    	    a1LVNoBS_Tauplus = Ntp->PFTau_PionsP4(Tauplus,0) + Ntp->PFTau_PionsP4(Tauplus,1) + Ntp->PFTau_PionsP4(Tauplus,2);
	    cout<<"NoBS:"<<endl;
    	    solutionsNoBS=tauPairMomentumSolutions(TauminusDirectionNoBS, a1LVNoBS_Tauminus, isMinusReal, TauplusDirectionNoBS, a1LVNoBS_Tauplus, isPlusReal);
	
    	    TauminusPairConstraintNoBS = solutionsNoBS.at(3);
	  
    	    TauplusPairConstraintNoBS = solutionsNoBS.at(7);
	  

    	    Scalc1NoBS.SortPions(HadPionsNoBS_minus, HadPionsNoBSCharge_minus);
    	    Scalc2NoBS.SortPions(HadPionsNoBS_plus, HadPionsNoBSCharge_plus);

    	    tauandprodNoBSminus.push_back(TauminusPairConstraintNoBS);
    	    tauandprodNoBSminus.push_back(HadPionsNoBS_minus.at(0));
    	    tauandprodNoBSminus.push_back(HadPionsNoBS_minus.at(1));
    	    tauandprodNoBSminus.push_back(HadPionsNoBS_minus.at(2));

    	    tauandprodNoBSplus.push_back(TauplusPairConstraintNoBS); 
    	    tauandprodNoBSplus.push_back(HadPionsNoBS_plus.at(0));  
    	    tauandprodNoBSplus.push_back(HadPionsNoBS_plus.at(1)); 
    	    tauandprodNoBSplus.push_back(HadPionsNoBS_plus.at(2)); 

    	    if(HadPionsNoBS_minus!=HadPionsNoBS_plus && HadPionsNoBS_minus!=VectZeroLV && HadPionsNoBS_plus!=VectZeroLV && TauminusPairConstraintNoBS!=TauplusPairConstraintNoBS && TauminusPairConstraintNoBS!=zeroLV && TauplusPairConstraintNoBS!=zeroLV && hasNoBS) passParticlesNoBS=true;
    	  }

    	SCalculator Scalc1NoBSRefittedPions("a1");
    	SCalculator Scalc2NoBSRefittedPions("a1");
	
    	HadRefittedPions_minus.push_back(Ntp->PFTauRefit_PionsP4(Tauminus,0));
    	HadRefittedPions_minus.push_back(Ntp->PFTauRefit_PionsP4(Tauminus,1));
    	HadRefittedPions_minus.push_back(Ntp->PFTauRefit_PionsP4(Tauminus,2));

    	HadRefittedPionsCharge_minus.push_back(Ntp->PFTau_PionsCharge(Tauminus, 0));
    	HadRefittedPionsCharge_minus.push_back(Ntp->PFTau_PionsCharge(Tauminus, 1));
    	HadRefittedPionsCharge_minus.push_back(Ntp->PFTau_PionsCharge(Tauminus, 2));
	
    	HadRefittedPions_plus.push_back(Ntp->PFTauRefit_PionsP4(Tauplus,0));
    	HadRefittedPions_plus.push_back(Ntp->PFTauRefit_PionsP4(Tauplus,1));
    	HadRefittedPions_plus.push_back(Ntp->PFTauRefit_PionsP4(Tauplus,2));
	
    	HadRefittedPionsCharge_plus.push_back(Ntp->PFTau_PionsCharge(Tauplus, 0));
    	HadRefittedPionsCharge_plus.push_back(Ntp->PFTau_PionsCharge(Tauplus, 1));
    	HadRefittedPionsCharge_plus.push_back(Ntp->PFTau_PionsCharge(Tauplus, 2));
	
    	if(hasNoBS)
    	  {
		  
    	    a1LVRefittedPions_Tauminus = Ntp->PFTauRefit_PionsP4(Tauminus,0) + Ntp->PFTauRefit_PionsP4(Tauminus,1) + Ntp->PFTauRefit_PionsP4(Tauminus,2);	
    	    a1LVRefittedPions_Tauplus = Ntp->PFTauRefit_PionsP4(Tauplus,0) + Ntp->PFTauRefit_PionsP4(Tauplus,1) + Ntp->PFTauRefit_PionsP4(Tauplus,2);
	    cout<<"NoBS Refit:"<<endl;
    	    solutionsNoBSRefittedPions=tauPairMomentumSolutions(TauminusDirectionNoBS, a1LVRefittedPions_Tauminus, isMinusReal, TauplusDirectionNoBS, a1LVRefittedPions_Tauplus, isPlusReal);
	
    	    TauminusPairConstraintNoBSRefittedPions = solutionsNoBSRefittedPions.at(3);
	  
    	    TauplusPairConstraintNoBSRefittedPions = solutionsNoBSRefittedPions.at(7);
	  

    	    Scalc1NoBSRefittedPions.SortPions(HadRefittedPions_minus, HadRefittedPionsCharge_minus);
    	    Scalc2NoBSRefittedPions.SortPions(HadRefittedPions_plus, HadRefittedPionsCharge_plus);

    	    tauandprodNoBSRefittedPionsminus.push_back(TauminusPairConstraintNoBSRefittedPions);
    	    tauandprodNoBSRefittedPionsminus.push_back(HadRefittedPions_minus.at(0));
    	    tauandprodNoBSRefittedPionsminus.push_back(HadRefittedPions_minus.at(1));
    	    tauandprodNoBSRefittedPionsminus.push_back(HadRefittedPions_minus.at(2));

    	    tauandprodNoBSRefittedPionsplus.push_back(TauplusPairConstraintNoBSRefittedPions); 
    	    tauandprodNoBSRefittedPionsplus.push_back(HadRefittedPions_plus.at(0));  
    	    tauandprodNoBSRefittedPionsplus.push_back(HadRefittedPions_plus.at(1)); 
    	    tauandprodNoBSRefittedPionsplus.push_back(HadRefittedPions_plus.at(2)); 

    	    if(HadRefittedPions_minus!=HadRefittedPions_plus && HadRefittedPions_minus!=VectZeroLV && HadRefittedPions_plus!=VectZeroLV && TauminusPairConstraintNoBSRefittedPions!=TauplusPairConstraintNoBSRefittedPions && TauminusPairConstraintNoBSRefittedPions!=zeroLV && TauplusPairConstraintNoBSRefittedPions!=zeroLV && hasNoBS) passParticlesNoBSRefittedPions=true;
    	  }
	
    	//Reco

    	HadPions_minus.push_back(Ntp->PFTau_PionsP4(Tauminus,0));
    	HadPions_minus.push_back(Ntp->PFTau_PionsP4(Tauminus,1));
    	HadPions_minus.push_back(Ntp->PFTau_PionsP4(Tauminus,2));
	
    	HadPionsCharge_minus.push_back(Ntp->PFTau_PionsCharge(Tauminus, 0));
    	HadPionsCharge_minus.push_back(Ntp->PFTau_PionsCharge(Tauminus, 1));
    	HadPionsCharge_minus.push_back(Ntp->PFTau_PionsCharge(Tauminus, 2));
	
	
    	HadPions_plus.push_back(Ntp->PFTau_PionsP4(Tauplus,0));
    	HadPions_plus.push_back(Ntp->PFTau_PionsP4(Tauplus,1));
    	HadPions_plus.push_back(Ntp->PFTau_PionsP4(Tauplus,2));
	
    	HadPionsCharge_plus.push_back(Ntp->PFTau_PionsCharge(Tauplus, 0));
    	HadPionsCharge_plus.push_back(Ntp->PFTau_PionsCharge(Tauplus, 1));
    	HadPionsCharge_plus.push_back(Ntp->PFTau_PionsCharge(Tauplus, 2));
	  
    	tauPrimaryVertex = Ntp->PVtx();
    	TauminusDirection = TauminusSecondaryVertex - tauPrimaryVertex;
    	TauplusDirection = TauplusSecondaryVertex - tauPrimaryVertex;
		  
    	a1LV_Tauminus = Ntp->PFTau_PionsP4(Tauminus,0) + Ntp->PFTau_PionsP4(Tauminus,1) + Ntp->PFTau_PionsP4(Tauminus,2);	
    	a1LV_Tauplus = Ntp->PFTau_PionsP4(Tauplus,0) + Ntp->PFTau_PionsP4(Tauplus,1) + Ntp->PFTau_PionsP4(Tauplus,2);
	// cout<<"Dir: ";
	// TauminusDirection.Unit().Print();
	// Ntp->GetTruthTauLV(5,0).Vect().Unit().Print();
	// cout<<"Truth a1-: ";
	// (Pions1Truth.at(0)+Pions1Truth.at(1)+Pions1Truth.at(2)).Vect().Unit().Print();
	// cout<<"Truth a1+: ";
	//(Pions2Truth.at(0)+Pions2Truth.at(1)+Pions2Truth.at(2)).Vect().Unit().Print();
    	solutions=tauPairMomentumSolutions(TauminusDirection, a1LV_Tauminus, isMinusReal, TauplusDirection, a1LV_Tauplus, isPlusReal);
	

	TauminusPairConstraint = solutions.at(3);
	TauminusPairConstraint.Vect().Unit().Print();
    	TauplusPairConstraint = solutions.at(7);

	
	
    	SCalculator Scalc1("a1");
    	SCalculator Scalc2("a1");

    	Scalc1.SortPions(HadPions_minus, HadPionsCharge_minus);
    	Scalc2.SortPions(HadPions_plus, HadPionsCharge_plus);

    	tauandprodminus.push_back(TauminusPairConstraint);
    	tauandprodminus.push_back(HadPions_minus.at(0));
    	tauandprodminus.push_back(HadPions_minus.at(1));
    	tauandprodminus.push_back(HadPions_minus.at(2));

    	tauandprodplus.push_back(TauplusPairConstraint); 
    	tauandprodplus.push_back(HadPions_plus.at(0));  
    	tauandprodplus.push_back(HadPions_plus.at(1)); 
    	tauandprodplus.push_back(HadPions_plus.at(2)); 

    	if (HadPions_minus!=HadPions_plus && HadPions_minus!=VectZeroLV && HadPions_plus!=VectZeroLV && TauminusPairConstraint!=TauplusPairConstraint && TauminusPairConstraint!=zeroLV && TauplusPairConstraint!=zeroLV) passParticlesReco=true;

    
    	//Old
	  
    	tauPrimaryVertexOld = Ntp->PVtx_Refit();
    	TauminusDirectionOld = TauminusSecondaryVertex - tauPrimaryVertexOld;
    	TauplusDirectionOld = TauplusSecondaryVertex - tauPrimaryVertexOld;
	
    	solutionsOld=tauPairMomentumSolutions(TauminusDirectionOld, a1LV_Tauminus, isMinusReal, TauplusDirectionOld, a1LV_Tauplus, isPlusReal);

    	TauminusPairConstraintOld = solutionsOld.at(3);
	
    	TauplusPairConstraintOld = solutionsOld.at(7);
	  
	
    	SCalculator Scalc1Old("a1");
    	SCalculator Scalc2Old("a1");

    	tauandprodOldminus.push_back(TauminusPairConstraintOld);
    	tauandprodOldminus.push_back(HadPions_minus.at(0));
    	tauandprodOldminus.push_back(HadPions_minus.at(1));
    	tauandprodOldminus.push_back(HadPions_minus.at(2));

    	tauandprodOldplus.push_back(TauplusPairConstraintOld); 
    	tauandprodOldplus.push_back(HadPions_plus.at(0));  
    	tauandprodOldplus.push_back(HadPions_plus.at(1)); 
    	tauandprodOldplus.push_back(HadPions_plus.at(2)); 

    	if (HadPions_minus!=HadPions_plus && HadPions_minus!=VectZeroLV && HadPions_plus!=VectZeroLV && TauminusPairConstraintOld!=TauplusPairConstraintOld && TauminusPairConstraintOld!=zeroLV && TauplusPairConstraintOld!=zeroLV) passParticlesOld=true;

    	//Old Refitted Pions
	  

    	solutionsOldRefittedPions=tauPairMomentumSolutions(TauminusDirectionOld, a1LVRefittedPions_Tauminus, isMinusReal, TauplusDirectionOld, a1LVRefittedPions_Tauplus, isPlusReal);

    	TauminusPairConstraintOldRefittedPions = solutionsOldRefittedPions.at(3);
	
    	TauplusPairConstraintOldRefittedPions = solutionsOldRefittedPions.at(7);
	  
	
    	SCalculator Scalc1OldRefittedPions("a1");
    	SCalculator Scalc2OldRefittedPions("a1");

    	tauandprodOldRefittedPionsminus.push_back(TauminusPairConstraintOldRefittedPions);
    	tauandprodOldRefittedPionsminus.push_back(HadRefittedPions_minus.at(0));
    	tauandprodOldRefittedPionsminus.push_back(HadRefittedPions_minus.at(1));
    	tauandprodOldRefittedPionsminus.push_back(HadRefittedPions_minus.at(2));

    	tauandprodOldRefittedPionsplus.push_back(TauplusPairConstraintOldRefittedPions); 
    	tauandprodOldRefittedPionsplus.push_back(HadRefittedPions_plus.at(0));  
    	tauandprodOldRefittedPionsplus.push_back(HadRefittedPions_plus.at(1)); 
    	tauandprodOldRefittedPionsplus.push_back(HadRefittedPions_plus.at(2)); 

    	if (HadRefittedPions_minus!=HadRefittedPions_plus && HadRefittedPions_minus!=VectZeroLV && HadRefittedPions_plus!=VectZeroLV && TauminusPairConstraintOldRefittedPions!=TauplusPairConstraintOldRefittedPions && TauminusPairConstraintOldRefittedPions!=zeroLV && TauplusPairConstraintOldRefittedPions!=zeroLV) passParticlesOldRefittedPions=true;

    	//TIP
	  	
    	SCalculator Scalc1TIP("a1");
    	SCalculator Scalc2TIP("a1");
    	SCalculator Scalc1TIPRefittedPions("a1");
    	SCalculator Scalc2TIPRefittedPions("a1");
    	if(Ntp->PFTau_TIP_primaryVertex_pos_Size(Tauminus)==3 && Ntp->PFTau_TIP_primaryVertex_pos_Size(Tauplus)==3)
    	  {
	    
    	    TauminusDirectionTIP = TauminusSecondaryVertex - Ntp->PFTau_TIP_primaryVertex_pos(Tauminus);
	    
    	    TauplusDirectionTIP = TauplusSecondaryVertex - Ntp->PFTau_TIP_primaryVertex_pos(Tauplus);
	   
    	    solutionsTIP=tauPairMomentumSolutions(TauminusDirectionTIP, a1LV_Tauminus, isMinusReal, TauplusDirectionTIP, a1LV_Tauplus, isPlusReal);
	    
    	    TauminusPairConstraintTIP = solutionsTIP.at(3);
	
    	    TauplusPairConstraintTIP = solutionsTIP.at(7);
	    
	    
    	    tauandprodTIPminus.push_back(TauminusPairConstraintTIP);
    	    tauandprodTIPminus.push_back(HadPions_minus.at(0));
    	    tauandprodTIPminus.push_back(HadPions_minus.at(1));
    	    tauandprodTIPminus.push_back(HadPions_minus.at(2));

    	    tauandprodTIPplus.push_back(TauplusPairConstraintTIP); 
    	    tauandprodTIPplus.push_back(HadPions_plus.at(0));  
    	    tauandprodTIPplus.push_back(HadPions_plus.at(1)); 
    	    tauandprodTIPplus.push_back(HadPions_plus.at(2)); 

    	    if (HadPions_minus!=HadPions_plus && HadPions_minus!=VectZeroLV && HadPions_plus!=VectZeroLV && TauminusPairConstraintTIP!=TauplusPairConstraintTIP && TauminusPairConstraintTIP!=zeroLV && TauplusPairConstraintTIP!=zeroLV) passParticlesTIP=true;
	    

    	    //TIP Refitted Pions
	  

    	    solutionsTIPRefittedPions=tauPairMomentumSolutions(TauminusDirectionTIP, a1LVRefittedPions_Tauminus, isMinusReal, TauplusDirectionTIP, a1LVRefittedPions_Tauplus, isPlusReal);

    	    TauminusPairConstraintTIPRefittedPions = solutionsTIPRefittedPions.at(3);
	
    	    TauplusPairConstraintTIPRefittedPions = solutionsTIPRefittedPions.at(7);
	  

    	    tauandprodTIPRefittedPionsminus.push_back(TauminusPairConstraintTIPRefittedPions);
    	    tauandprodTIPRefittedPionsminus.push_back(HadRefittedPions_minus.at(0));
    	    tauandprodTIPRefittedPionsminus.push_back(HadRefittedPions_minus.at(1));
    	    tauandprodTIPRefittedPionsminus.push_back(HadRefittedPions_minus.at(2));

    	    tauandprodTIPRefittedPionsplus.push_back(TauplusPairConstraintTIPRefittedPions); 
    	    tauandprodTIPRefittedPionsplus.push_back(HadRefittedPions_plus.at(0));  
    	    tauandprodTIPRefittedPionsplus.push_back(HadRefittedPions_plus.at(1)); 
    	    tauandprodTIPRefittedPionsplus.push_back(HadRefittedPions_plus.at(2)); 

    	    if (HadRefittedPions_minus!=HadRefittedPions_plus && HadRefittedPions_minus!=VectZeroLV && HadRefittedPions_plus!=VectZeroLV && TauminusPairConstraintTIPRefittedPions!=TauplusPairConstraintTIPRefittedPions && TauminusPairConstraintTIPRefittedPions!=zeroLV && TauplusPairConstraintTIPRefittedPions!=zeroLV) passParticlesTIPRefittedPions=true;
    	  }

    	// SVFit	

    	PionsSVFit1.push_back(Ntp->PFTau_PionsP4(Tauminus,0));
    	PionsSVFit1.push_back(Ntp->PFTau_PionsP4(Tauminus,1));
    	PionsSVFit1.push_back(Ntp->PFTau_PionsP4(Tauminus,2));
	
    	PionsSVFit1Charge.push_back(Ntp->PFTau_PionsCharge(Tauminus, 0));
    	PionsSVFit1Charge.push_back(Ntp->PFTau_PionsCharge(Tauminus, 1));
    	PionsSVFit1Charge.push_back(Ntp->PFTau_PionsCharge(Tauminus, 2));

    	PionsSVFit2.push_back(Ntp->PFTau_PionsP4(Tauplus,0));
    	PionsSVFit2.push_back(Ntp->PFTau_PionsP4(Tauplus,1));
    	PionsSVFit2.push_back(Ntp->PFTau_PionsP4(Tauplus,2));
	
    	PionsSVFit2Charge.push_back(Ntp->PFTau_PionsCharge(Tauplus, 0));
    	PionsSVFit2Charge.push_back(Ntp->PFTau_PionsCharge(Tauplus, 1));
    	PionsSVFit2Charge.push_back(Ntp->PFTau_PionsCharge(Tauplus, 2));
	
	
    	SCalculator ScalcSVFit1("a1");
    	SCalculator ScalcSVFit2("a1");

    	// cout<<"Nu: "<<endl;
    	// ScalcSVFit1.Boost((Tauminussvfit-(Ntp->PFTau_PionsP4(Tauminus,0)+Ntp->PFTau_PionsP4(Tauminus,1)+Ntp->PFTau_PionsP4(Tauminus,2))),Tauminussvfit+Tauplussvfit).Print();
    	// cout<<"NuRefit: "<<endl;
    	// ScalcSVFit1.Boost((Tauminussvfit-(Ntp->PFTauRefit_PionsP4(Tauminus,0)+Ntp->PFTauRefit_PionsP4(Tauminus,1)+Ntp->PFTauRefit_PionsP4(Tauminus,2))),Tauminussvfit+Tauplussvfit).Print();

    	ScalcSVFit1.SortPions(PionsSVFit1, PionsSVFit1Charge);
    	ScalcSVFit2.SortPions(PionsSVFit2, PionsSVFit2Charge);
	       

    	tauandprodSVFitminus.push_back(Tauminussvfit);
    	cout.precision(10);

    	tauandprodSVFitminus.push_back(PionsSVFit1.at(0));
    	tauandprodSVFitminus.push_back(PionsSVFit1.at(1));
    	tauandprodSVFitminus.push_back(PionsSVFit1.at(2));
    	tauandprodSVFitplus.push_back(Tauplussvfit);
    	tauandprodSVFitplus.push_back(PionsSVFit2.at(0));   
    	tauandprodSVFitplus.push_back(PionsSVFit2.at(1));   
    	tauandprodSVFitplus.push_back(PionsSVFit2.at(2));   

    	if((PionsSVFit1!=PionsSVFit2) && (PionsSVFit1!=VectZeroLV) && (PionsSVFit2!=VectZeroLV) && tauandprodSVFitminus.at(0)!=zeroLV && tauandprodSVFitplus.at(0)!=zeroLV && tauandprodSVFitminus.at(0)!=tauandprodSVFitplus.at(0)  && svfitAlgo1.isValidSolution()) passParticlesSVFit=true;


    	if(passParticlesTruth && passParticlesReco && passParticlesNoBS && passParticlesNoBSRefittedPions && passParticlesOld && passParticlesOldRefittedPions && passParticlesTIP && passParticlesTIPRefittedPions &&passParticlesSVFit)
    	  {
	
    	    // Truth
    	    cout<<"Truth"<<endl;
    	    Scalc1Truth.Configure(tauandprodTruthminus,tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0), -1);
    	    TVector3 h1Truth=Scalc1Truth.pv();
	  
    	    Scalc2Truth.Configure(tauandprodTruthplus,tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0), +1);
    	    TVector3 h2Truth=Scalc2Truth.pv();
	      
    	    TLorentzVector tauminusTruth_HRF = Scalc1Truth.Boost(tauandprodTruthminus.at(0),tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0));
    	    TLorentzVector tauplusTruth_HRF  = Scalc2Truth.Boost(tauandprodTruthplus.at(0),tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0));
		
    	    double h1TruthNorm=1./h1Truth.Mag();
    	    double h2TruthNorm=1./h2Truth.Mag();
		
    	    h1Truth=h1Truth*h1TruthNorm;
    	    h2Truth=h2Truth*h2TruthNorm;
		
    	    double normk1Truth=1./((h1Truth.Cross(tauminusTruth_HRF.Vect().Unit())).Mag());
    	    double normk2Truth=1./((h2Truth.Cross(tauplusTruth_HRF.Vect().Unit())).Mag());
		
    	    TVector3 k1Truth = (h1Truth.Cross(tauminusTruth_HRF.Vect().Unit()))*normk1Truth;
    	    TVector3 k2Truth = (h2Truth.Cross(tauplusTruth_HRF.Vect().Unit()))*normk2Truth;

    	    if(std::isnan(h1TruthNorm)!=true && std::isnan(h2TruthNorm)!=true && std::isnan(normk1Truth)!=true && std::isnan(normk2Truth)!=true) passVariablesTruth=true;
	    
    	    // Reco
    	    cout<<"Reco"<<endl;
    	    Scalc1.Configure(tauandprodminus,tauandprodminus.at(0)+tauandprodplus.at(0), -1);
    	    TVector3 h1=Scalc1.pv();

    	    Scalc2.Configure(tauandprodplus,tauandprodminus.at(0)+tauandprodplus.at(0), +1);
    	    TVector3 h2=Scalc2.pv();

    	    TLorentzVector tauminus_HRF = Scalc1.Boost(tauandprodminus.at(0),tauandprodminus.at(0)+tauandprodplus.at(0));
    	    TLorentzVector tauplus_HRF  = Scalc2.Boost(tauandprodplus.at(0),tauandprodminus.at(0)+tauandprodplus.at(0));

    	    double h1Norm=1./h1.Mag();
    	    double h2Norm=1./h2.Mag();
		
    	    h1=h1*h1Norm;
    	    h2=h2*h2Norm;

    	    double k1Norm=1./((h1.Cross(tauminus_HRF.Vect().Unit())).Mag());
    	    double k2Norm=1./((h2.Cross(tauplus_HRF.Vect().Unit())).Mag());
    	    TVector3 k1 = (h1.Cross(tauminus_HRF.Vect().Unit()))*k1Norm;
    	    TVector3 k2 = (h2.Cross(tauplus_HRF.Vect().Unit()))*k2Norm;

    	    if(std::isnan(h1Norm)!=true && std::isnan(h2Norm)!=true && std::isnan(k1Norm)!=true && std::isnan(k2Norm)!=true) passVariablesReco=true;
	    
    	    // Reco NoBS
    	    cout<<"Reco NoBS"<<endl;
    	    Scalc1NoBS.Configure(tauandprodNoBSminus,tauandprodNoBSminus.at(0)+tauandprodNoBSplus.at(0), -1);
    	    TVector3 h1NoBS=Scalc1NoBS.pv();

    	    Scalc2NoBS.Configure(tauandprodNoBSplus,tauandprodNoBSminus.at(0)+tauandprodNoBSplus.at(0), +1);
    	    TVector3 h2NoBS=Scalc2NoBS.pv();

    	    TLorentzVector tauminusNoBS_HRF = Scalc1NoBS.Boost(tauandprodNoBSminus.at(0),tauandprodNoBSminus.at(0)+tauandprodNoBSplus.at(0));
    	    TLorentzVector tauplusNoBS_HRF  = Scalc2NoBS.Boost(tauandprodNoBSplus.at(0),tauandprodNoBSminus.at(0)+tauandprodNoBSplus.at(0));

    	    double h1NormNoBS=1./h1NoBS.Mag();
    	    double h2NormNoBS=1./h2NoBS.Mag();
		
    	    h1NoBS=h1NoBS*h1NormNoBS;
    	    h2NoBS=h2NoBS*h2NormNoBS;

    	    double k1NormNoBS=1./((h1NoBS.Cross(tauminusNoBS_HRF.Vect().Unit())).Mag());
    	    double k2NormNoBS=1./((h2NoBS.Cross(tauplusNoBS_HRF.Vect().Unit())).Mag());
    	    TVector3 k1NoBS = (h1NoBS.Cross(tauminusNoBS_HRF.Vect().Unit()))*k1NormNoBS;
    	    TVector3 k2NoBS = (h2NoBS.Cross(tauplusNoBS_HRF.Vect().Unit()))*k2NormNoBS;

    	    if(std::isnan(h1NormNoBS)!=true && std::isnan(h2NormNoBS)!=true && std::isnan(k1NormNoBS)!=true && std::isnan(k2NormNoBS)!=true) passVariablesNoBS=true;
    	    
    	    // Reco NoBS Refitted Pions
    	    cout<<"Reco NoBS Refitted Pions"<<endl;
    	    Scalc1NoBSRefittedPions.Configure(tauandprodNoBSRefittedPionsminus,tauandprodNoBSRefittedPionsminus.at(0)+tauandprodNoBSRefittedPionsplus.at(0), -1);
    	    TVector3 h1NoBSRefittedPions=Scalc1NoBSRefittedPions.pv();

    	    Scalc2NoBSRefittedPions.Configure(tauandprodNoBSRefittedPionsplus,tauandprodNoBSRefittedPionsminus.at(0)+tauandprodNoBSRefittedPionsplus.at(0), +1);
    	    TVector3 h2NoBSRefittedPions=Scalc2NoBSRefittedPions.pv();

    	    TLorentzVector tauminusNoBSRefittedPions_HRF = Scalc1NoBSRefittedPions.Boost(tauandprodNoBSRefittedPionsminus.at(0),tauandprodNoBSRefittedPionsminus.at(0)+tauandprodNoBSRefittedPionsplus.at(0));
    	    TLorentzVector tauplusNoBSRefittedPions_HRF  = Scalc2NoBSRefittedPions.Boost(tauandprodNoBSRefittedPionsplus.at(0),tauandprodNoBSRefittedPionsminus.at(0)+tauandprodNoBSRefittedPionsplus.at(0));

    	    double h1NormNoBSRefittedPions=1./h1NoBSRefittedPions.Mag();
    	    double h2NormNoBSRefittedPions=1./h2NoBSRefittedPions.Mag();
		
    	    h1NoBSRefittedPions=h1NoBSRefittedPions*h1NormNoBSRefittedPions;
    	    h2NoBSRefittedPions=h2NoBSRefittedPions*h2NormNoBSRefittedPions;

    	    double k1NormNoBSRefittedPions=1./((h1NoBSRefittedPions.Cross(tauminusNoBSRefittedPions_HRF.Vect().Unit())).Mag());
    	    double k2NormNoBSRefittedPions=1./((h2NoBSRefittedPions.Cross(tauplusNoBSRefittedPions_HRF.Vect().Unit())).Mag());
    	    TVector3 k1NoBSRefittedPions = (h1NoBSRefittedPions.Cross(tauminusNoBSRefittedPions_HRF.Vect().Unit()))*k1NormNoBSRefittedPions;
    	    TVector3 k2NoBSRefittedPions = (h2NoBSRefittedPions.Cross(tauplusNoBSRefittedPions_HRF.Vect().Unit()))*k2NormNoBSRefittedPions;

    	    if(std::isnan(h1NormNoBSRefittedPions)!=true && std::isnan(h2NormNoBSRefittedPions)!=true && std::isnan(k1NormNoBSRefittedPions)!=true && std::isnan(k2NormNoBSRefittedPions)!=true) passVariablesNoBSRefittedPions=true;
	    

    	    // Reco Old
    	    cout<<"Reco Old"<<endl;
    	    Scalc1Old.Configure(tauandprodOldminus,tauandprodOldminus.at(0)+tauandprodOldplus.at(0), -1);
    	    TVector3 h1Old=Scalc1Old.pv();

    	    Scalc2Old.Configure(tauandprodOldplus,tauandprodOldminus.at(0)+tauandprodOldplus.at(0), +1);
    	    TVector3 h2Old=Scalc2Old.pv();

    	    TLorentzVector tauminusOld_HRF = Scalc1Old.Boost(tauandprodOldminus.at(0),tauandprodOldminus.at(0)+tauandprodOldplus.at(0));
    	    TLorentzVector tauplusOld_HRF  = Scalc2Old.Boost(tauandprodOldplus.at(0),tauandprodOldminus.at(0)+tauandprodOldplus.at(0));

    	    double h1NormOld=1./h1Old.Mag();
    	    double h2NormOld=1./h2Old.Mag();
		
    	    h1Old=h1Old*h1NormOld;
    	    h2Old=h2Old*h2NormOld;

    	    double k1NormOld=1./((h1Old.Cross(tauminusOld_HRF.Vect().Unit())).Mag());
    	    double k2NormOld=1./((h2Old.Cross(tauplusOld_HRF.Vect().Unit())).Mag());
    	    TVector3 k1Old = (h1Old.Cross(tauminusOld_HRF.Vect().Unit()))*k1NormOld;
    	    TVector3 k2Old = (h2Old.Cross(tauplusOld_HRF.Vect().Unit()))*k2NormOld;

    	    if(std::isnan(h1NormOld)!=true && std::isnan(h2NormOld)!=true && std::isnan(k1NormOld)!=true && std::isnan(k2NormOld)!=true) passVariablesOld=true;
    	    
    	    // Reco Old Refitted Pions
    	    cout<<"Reco Old Refitted Pions"<<endl;
    	    Scalc1OldRefittedPions.Configure(tauandprodOldRefittedPionsminus,tauandprodOldRefittedPionsminus.at(0)+tauandprodOldRefittedPionsplus.at(0), -1);
    	    TVector3 h1OldRefittedPions=Scalc1OldRefittedPions.pv();

    	    Scalc2OldRefittedPions.Configure(tauandprodOldRefittedPionsplus,tauandprodOldRefittedPionsminus.at(0)+tauandprodOldRefittedPionsplus.at(0), +1);
    	    TVector3 h2OldRefittedPions=Scalc2OldRefittedPions.pv();

    	    TLorentzVector tauminusOldRefittedPions_HRF = Scalc1OldRefittedPions.Boost(tauandprodOldRefittedPionsminus.at(0),tauandprodOldRefittedPionsminus.at(0)+tauandprodOldRefittedPionsplus.at(0));
    	    TLorentzVector tauplusOldRefittedPions_HRF  = Scalc2OldRefittedPions.Boost(tauandprodOldRefittedPionsplus.at(0),tauandprodOldRefittedPionsminus.at(0)+tauandprodOldRefittedPionsplus.at(0));

    	    double h1NormOldRefittedPions=1./h1OldRefittedPions.Mag();
    	    double h2NormOldRefittedPions=1./h2OldRefittedPions.Mag();
		
    	    h1OldRefittedPions=h1OldRefittedPions*h1NormOldRefittedPions;
    	    h2OldRefittedPions=h2OldRefittedPions*h2NormOldRefittedPions;

    	    double k1NormOldRefittedPions=1./((h1OldRefittedPions.Cross(tauminusOldRefittedPions_HRF.Vect().Unit())).Mag());
    	    double k2NormOldRefittedPions=1./((h2OldRefittedPions.Cross(tauplusOldRefittedPions_HRF.Vect().Unit())).Mag());
    	    TVector3 k1OldRefittedPions = (h1OldRefittedPions.Cross(tauminusOldRefittedPions_HRF.Vect().Unit()))*k1NormOldRefittedPions;
    	    TVector3 k2OldRefittedPions = (h2OldRefittedPions.Cross(tauplusOldRefittedPions_HRF.Vect().Unit()))*k2NormOldRefittedPions;

    	    if(std::isnan(h1NormOldRefittedPions)!=true && std::isnan(h2NormOldRefittedPions)!=true && std::isnan(k1NormOldRefittedPions)!=true && std::isnan(k2NormOldRefittedPions)!=true) passVariablesOldRefittedPions=true;

	    
    	    // Reco TIP
    	    cout<<"Reco TIP"<<endl;
    	    Scalc1TIP.Configure(tauandprodTIPminus,tauandprodTIPminus.at(0)+tauandprodTIPplus.at(0), -1);
    	    TVector3 h1TIP=Scalc1TIP.pv();

    	    Scalc2TIP.Configure(tauandprodTIPplus,tauandprodTIPminus.at(0)+tauandprodTIPplus.at(0), +1);
    	    TVector3 h2TIP=Scalc2TIP.pv();

    	    TLorentzVector tauminusTIP_HRF = Scalc1TIP.Boost(tauandprodTIPminus.at(0),tauandprodTIPminus.at(0)+tauandprodTIPplus.at(0));
    	    TLorentzVector tauplusTIP_HRF  = Scalc2TIP.Boost(tauandprodTIPplus.at(0),tauandprodTIPminus.at(0)+tauandprodTIPplus.at(0));

    	    double h1NormTIP=1./h1TIP.Mag();
    	    double h2NormTIP=1./h2TIP.Mag();
		
    	    h1TIP=h1TIP*h1NormTIP;
    	    h2TIP=h2TIP*h2NormTIP;

    	    double k1NormTIP=1./((h1TIP.Cross(tauminusTIP_HRF.Vect().Unit())).Mag());
    	    double k2NormTIP=1./((h2TIP.Cross(tauplusTIP_HRF.Vect().Unit())).Mag());
    	    TVector3 k1TIP = (h1TIP.Cross(tauminusTIP_HRF.Vect().Unit()))*k1NormTIP;
    	    TVector3 k2TIP = (h2TIP.Cross(tauplusTIP_HRF.Vect().Unit()))*k2NormTIP;

    	    if(std::isnan(h1NormTIP)!=true && std::isnan(h2NormTIP)!=true && std::isnan(k1NormTIP)!=true && std::isnan(k2NormTIP)!=true) passVariablesTIP=true;
    	    
    	    // Reco TIP Refitted Pions
    	    cout<<"Reco TIP Refitted Pions"<<endl;
    	    Scalc1TIPRefittedPions.Configure(tauandprodTIPRefittedPionsminus,tauandprodTIPRefittedPionsminus.at(0)+tauandprodTIPRefittedPionsplus.at(0), -1);
    	    TVector3 h1TIPRefittedPions=Scalc1TIPRefittedPions.pv();

    	    Scalc2TIPRefittedPions.Configure(tauandprodTIPRefittedPionsplus,tauandprodTIPRefittedPionsminus.at(0)+tauandprodTIPRefittedPionsplus.at(0), +1);
    	    TVector3 h2TIPRefittedPions=Scalc2TIPRefittedPions.pv();

    	    TLorentzVector tauminusTIPRefittedPions_HRF = Scalc1TIPRefittedPions.Boost(tauandprodTIPRefittedPionsminus.at(0),tauandprodTIPRefittedPionsminus.at(0)+tauandprodTIPRefittedPionsplus.at(0));
    	    TLorentzVector tauplusTIPRefittedPions_HRF  = Scalc2TIPRefittedPions.Boost(tauandprodTIPRefittedPionsplus.at(0),tauandprodTIPRefittedPionsminus.at(0)+tauandprodTIPRefittedPionsplus.at(0));

    	    double h1NormTIPRefittedPions=1./h1TIPRefittedPions.Mag();
    	    double h2NormTIPRefittedPions=1./h2TIPRefittedPions.Mag();
		
    	    h1TIPRefittedPions=h1TIPRefittedPions*h1NormTIPRefittedPions;
    	    h2TIPRefittedPions=h2TIPRefittedPions*h2NormTIPRefittedPions;

    	    double k1NormTIPRefittedPions=1./((h1TIPRefittedPions.Cross(tauminusTIPRefittedPions_HRF.Vect().Unit())).Mag());
    	    double k2NormTIPRefittedPions=1./((h2TIPRefittedPions.Cross(tauplusTIPRefittedPions_HRF.Vect().Unit())).Mag());
    	    TVector3 k1TIPRefittedPions = (h1TIPRefittedPions.Cross(tauminusTIPRefittedPions_HRF.Vect().Unit()))*k1NormTIPRefittedPions;
    	    TVector3 k2TIPRefittedPions = (h2TIPRefittedPions.Cross(tauplusTIPRefittedPions_HRF.Vect().Unit()))*k2NormTIPRefittedPions;

    	    if(std::isnan(h1NormTIPRefittedPions)!=true && std::isnan(h2NormTIPRefittedPions)!=true && std::isnan(k1NormTIPRefittedPions)!=true && std::isnan(k2NormTIPRefittedPions)!=true) passVariablesTIPRefittedPions=true;
	    
    	    // SVFit
    	    cout<<"SVFit"<<endl;
    	    ScalcSVFit1.Configure(tauandprodSVFitminus,tauandprodSVFitminus.at(0)+tauandprodSVFitplus.at(0), -1);
    	    TVector3 h1SVFit=ScalcSVFit1.pv();

    	    ScalcSVFit2.Configure(tauandprodSVFitplus,tauandprodSVFitminus.at(0)+tauandprodSVFitplus.at(0), +1);
    	    TVector3 h2SVFit=ScalcSVFit2.pv();
	      
    	    TLorentzVector tauSVFitminus_HRF = ScalcSVFit1.Boost(tauandprodSVFitminus.at(0),tauandprodSVFitminus.at(0)+tauandprodSVFitplus.at(0));
    	    TLorentzVector tauSVFitplus_HRF  = ScalcSVFit2.Boost(tauandprodSVFitplus.at(0),tauandprodSVFitminus.at(0)+tauandprodSVFitplus.at(0));

    	    double h1SVFitNorm=1./h1SVFit.Mag();
    	    double h2SVFitNorm=1./h2SVFit.Mag();
	      
    	    h1SVFit=h1SVFit*h1SVFitNorm;
    	    h2SVFit=h2SVFit*h2SVFitNorm;

    	    double k1SVFitNorm=1./((h1SVFit.Cross(tauSVFitminus_HRF.Vect().Unit())).Mag());
    	    double k2SVFitNorm=1./((h2SVFit.Cross(tauSVFitplus_HRF.Vect().Unit())).Mag());
    	    TVector3 k1SVFit = (h1SVFit.Cross(tauSVFitminus_HRF.Vect().Unit()))*k1SVFitNorm;
    	    TVector3 k2SVFit = (h2SVFit.Cross(tauSVFitplus_HRF.Vect().Unit()))*k2SVFitNorm;
		    
    	    if(std::isnan(h1SVFitNorm)!=true && std::isnan(h2SVFitNorm)!=true && std::isnan(k1SVFitNorm)!=true && std::isnan(k2SVFitNorm)!=true) passVariablesSVFit=true;
		  
    	    if(passVariablesTruth && passVariablesReco && passVariablesNoBS && passVariablesNoBSRefittedPions &&  passVariablesOld && passVariablesOldRefittedPions &&  passVariablesTIP && passVariablesTIPRefittedPions && passVariablesSVFit)
    	      {

		//cout<<GetThetaGJ(TauminusDirection,a1LV_Tauminus)<<"  "<<GetThetaGJMax(a1LV_Tauminus)<<endl;
		//cout<<(!isPhysical(TauminusPairConstraint.Vect(),a1LV_Tauminus))<<endl;
		if(!isPhysical(TauminusPairConstraint.Vect(),a1LV_Tauminus) || !isPhysical(TauplusPairConstraint.Vect(),a1LV_Tauplus))
		  {
		    //Unphysical.at(t).Fill(1.,w);
		    unphysicalBool=true;
		  }
		//else {unphysicalBool=false;Unphysical.at(t).Fill(0.,w);}
	

		if(!isPhysical(Tauminussvfit.Vect(),PionsSVFit1.at(0)+PionsSVFit1.at(1)+PionsSVFit1.at(2)) || !isPhysical(Tauplussvfit.Vect(),PionsSVFit2.at(0)+PionsSVFit2.at(1)+PionsSVFit2.at(2)))
		  {
	   
		    //Unphysical.at(t).Fill(1.,w);
		    unphysicalBool=true;
	      
		  }

		if (unphysicalBool==false)Unphysical.at(t).Fill(0.,w);
		if (unphysicalBool==true) Unphysical.at(t).Fill(1.,w);
		if(!unphysicalBool)
		  {
    		// Truth
		    
    		if(((h1Truth).Cross(h2Truth))*tauminusTruth_HRF.Vect().Unit()<=0){polarimetricAcopAngleTruthA1.at(t).Fill(TMath::ATan2((k1Truth.Cross(k2Truth)).Mag(),k1Truth*k2Truth),Wspin);}
    		else {polarimetricAcopAngleTruthA1.at(t).Fill(2*TMath::Pi()-TMath::ATan2((k1Truth.Cross(k2Truth)).Mag(),k1Truth*k2Truth),Wspin);}
		    
		
		  
    		// Reco

    		if(((h1.Cross(h2))*(tauminus_HRF.Vect().Unit()))<=0){ polarimetricAcopAngle.at(t).Fill(TMath::ATan2((k1.Cross(k2)).Mag(),k1*k2),Wspin);}
    		else {polarimetricAcopAngle.at(t).Fill(2.*TMath::Pi()-TMath::ATan2((k1.Cross(k2)).Mag(),k1*k2),Wspin);}

    		// NoBS
    		if(((h1NoBS.Cross(h2NoBS))*(tauminusNoBS_HRF.Vect().Unit()))<=0){ polarimetricAcopAngleNoBS.at(t).Fill(TMath::ATan2((k1NoBS.Cross(k2NoBS)).Mag(),k1NoBS*k2NoBS),Wspin);}
    		else {polarimetricAcopAngleNoBS.at(t).Fill(2.*TMath::Pi()-TMath::ATan2((k1NoBS.Cross(k2NoBS)).Mag(),k1NoBS*k2NoBS),Wspin);}

    		// NoBS Refitted Pions

    		if(((h1NoBSRefittedPions.Cross(h2NoBSRefittedPions))*(tauminusNoBSRefittedPions_HRF.Vect().Unit()))<=0){ polarimetricAcopAngleNoBSPionsRefit.at(t).Fill(TMath::ATan2((k1NoBSRefittedPions.Cross(k2NoBSRefittedPions)).Mag(),k1NoBSRefittedPions*k2NoBSRefittedPions),Wspin);}
    		else {polarimetricAcopAngleNoBSPionsRefit.at(t).Fill(2.*TMath::Pi()-TMath::ATan2((k1NoBSRefittedPions.Cross(k2NoBSRefittedPions)).Mag(),k1NoBSRefittedPions*k2NoBSRefittedPions),Wspin);}

    		// Old
    		if(((h1Old.Cross(h2Old))*(tauminusOld_HRF.Vect().Unit()))<=0){ polarimetricAcopAngleOld.at(t).Fill(TMath::ATan2((k1Old.Cross(k2Old)).Mag(),k1Old*k2Old),Wspin);}
    		else {polarimetricAcopAngleOld.at(t).Fill(2.*TMath::Pi()-TMath::ATan2((k1Old.Cross(k2Old)).Mag(),k1Old*k2Old),Wspin);}

    		// Old Refitted Pions

    		if(((h1OldRefittedPions.Cross(h2OldRefittedPions))*(tauminusOldRefittedPions_HRF.Vect().Unit()))<=0){ polarimetricAcopAngleOldPionsRefit.at(t).Fill(TMath::ATan2((k1OldRefittedPions.Cross(k2OldRefittedPions)).Mag(),k1OldRefittedPions*k2OldRefittedPions),Wspin);}
    		else {polarimetricAcopAngleOldPionsRefit.at(t).Fill(2.*TMath::Pi()-TMath::ATan2((k1OldRefittedPions.Cross(k2OldRefittedPions)).Mag(),k1OldRefittedPions*k2OldRefittedPions),Wspin);}


    		// TIP
    		if(((h1TIP.Cross(h2TIP))*(tauminusTIP_HRF.Vect().Unit()))<=0){ polarimetricAcopAngleTIP.at(t).Fill(TMath::ATan2((k1TIP.Cross(k2TIP)).Mag(),k1TIP*k2TIP),Wspin);}
    		else {polarimetricAcopAngleTIP.at(t).Fill(2.*TMath::Pi()-TMath::ATan2((k1TIP.Cross(k2TIP)).Mag(),k1TIP*k2TIP),Wspin);}
		    
    		// TIP Refitted Pions

    		if(((h1TIPRefittedPions.Cross(h2TIPRefittedPions))*(tauminusTIPRefittedPions_HRF.Vect().Unit()))<=0){ polarimetricAcopAngleTIPPionsRefit.at(t).Fill(TMath::ATan2((k1TIPRefittedPions.Cross(k2TIPRefittedPions)).Mag(),k1TIPRefittedPions*k2TIPRefittedPions),Wspin);}
    		else {polarimetricAcopAngleTIPPionsRefit.at(t).Fill(2.*TMath::Pi()-TMath::ATan2((k1TIPRefittedPions.Cross(k2TIPRefittedPions)).Mag(),k1TIPRefittedPions*k2TIPRefittedPions),Wspin);}
		
    		// SVFit

    		if(((h1SVFit).Cross(h2SVFit))*tauSVFitminus_HRF.Vect().Unit()<=0){polarimetricAcopAngleSVFitA1.at(t).Fill(TMath::ATan2((k1SVFit.Cross(k2SVFit)).Mag(),k1SVFit*k2SVFit),Wspin);}
    		else {polarimetricAcopAngleSVFitA1.at(t).Fill(2.*TMath::Pi()-TMath::ATan2((k1SVFit.Cross(k2SVFit)).Mag(),k1SVFit*k2SVFit),Wspin);}		    


    		PVXResol.at(t).Fill((Ntp->PVtx().X()-Ntp->PVtx_Gen().X())/Ntp->PVtx_Gen().X(),w);
    		PVXNoBSResol.at(t).Fill((tauPrimaryVertexNoBS.X()-Ntp->PVtx_Gen().X())/Ntp->PVtx_Gen().X(),w);
    		PVXOldRefitResol.at(t).Fill((Ntp->PVtx_Refit().X()-Ntp->PVtx_Gen().X())/Ntp->PVtx_Gen().X(),w);
    		PVXTIPRefitResol.at(t).Fill((Ntp->PFTau_TIP_primaryVertex_pos(Tauminus).X()-Ntp->PVtx_Gen().X())/Ntp->PVtx_Gen().X(),w);
    		PVXTIPRefitResol.at(t).Fill((Ntp->PFTau_TIP_primaryVertex_pos(Tauplus).X()-Ntp->PVtx_Gen().X())/Ntp->PVtx_Gen().X(),w);

    		PVYResol.at(t).Fill((Ntp->PVtx().Y()-Ntp->PVtx_Gen().Y())/Ntp->PVtx_Gen().Y(),w);
    		PVYNoBSResol.at(t).Fill((tauPrimaryVertexNoBS.Y()-Ntp->PVtx_Gen().Y())/Ntp->PVtx_Gen().Y(),w);
    		PVYOldRefitResol.at(t).Fill((Ntp->PVtx_Refit().Y()-Ntp->PVtx_Gen().Y())/Ntp->PVtx_Gen().Y(),w);
    		PVYTIPRefitResol.at(t).Fill((Ntp->PFTau_TIP_primaryVertex_pos(Tauminus).Y()-Ntp->PVtx_Gen().Y())/Ntp->PVtx_Gen().Y(),w);
    		PVYTIPRefitResol.at(t).Fill((Ntp->PFTau_TIP_primaryVertex_pos(Tauplus).Y()-Ntp->PVtx_Gen().Y())/Ntp->PVtx_Gen().Y(),w);

    		PVZResol.at(t).Fill((Ntp->PVtx().Z()-Ntp->PVtx_Gen().Z())/Ntp->PVtx_Gen().Z(),w);
    		PVZNoBSResol.at(t).Fill((tauPrimaryVertexNoBS.Z()-Ntp->PVtx_Gen().Z())/Ntp->PVtx_Gen().Z(),w);
    		PVZOldRefitResol.at(t).Fill((Ntp->PVtx_Refit().Z()-Ntp->PVtx_Gen().Z())/Ntp->PVtx_Gen().Z(),w);
    		PVZTIPRefitResol.at(t).Fill((Ntp->PFTau_TIP_primaryVertex_pos(Tauminus).Z()-Ntp->PVtx_Gen().Z())/Ntp->PVtx_Gen().Z(),w);
    		PVZTIPRefitResol.at(t).Fill((Ntp->PFTau_TIP_primaryVertex_pos(Tauplus).Z()-Ntp->PVtx_Gen().Z())/Ntp->PVtx_Gen().Z(),w);

    		TauSVFitPxResPull.at(t).Fill((tauandprodSVFitplus.at(0).X()-tauandprodTruthplus.at(0).X())/tauandprodTruthplus.at(0).X(),w);
    		TauSVFitPyResPull.at(t).Fill((tauandprodSVFitplus.at(0).Y()-tauandprodTruthplus.at(0).Y())/tauandprodTruthplus.at(0).Y(),w);
    		TauSVFitPzResPull.at(t).Fill((tauandprodSVFitplus.at(0).Z()-tauandprodTruthplus.at(0).Z())/tauandprodTruthplus.at(0).Z(),w);
    		TauSVFitEResPull.at(t).Fill((tauandprodSVFitplus.at(0).E()-tauandprodTruthplus.at(0).E())/tauandprodTruthplus.at(0).E(),w);
    		TauSVFitPxResPull.at(t).Fill((tauandprodSVFitminus.at(0).X()-tauandprodTruthminus.at(0).X())/tauandprodTruthminus.at(0).X(),w);
    		TauSVFitPyResPull.at(t).Fill((tauandprodSVFitminus.at(0).Y()-tauandprodTruthminus.at(0).Y())/tauandprodTruthminus.at(0).Y(),w);
    		TauSVFitPzResPull.at(t).Fill((tauandprodSVFitminus.at(0).Z()-tauandprodTruthminus.at(0).Z())/tauandprodTruthminus.at(0).Z(),w);
    		TauSVFitEResPull.at(t).Fill((tauandprodSVFitminus.at(0).E()-tauandprodTruthminus.at(0).E())/tauandprodTruthminus.at(0).E(),w);

    		TauPxResPull.at(t).Fill((tauandprodplus.at(0).X()-tauandprodTruthplus.at(0).X())/tauandprodTruthplus.at(0).X(),w);
    		TauPyResPull.at(t).Fill((tauandprodplus.at(0).Y()-tauandprodTruthplus.at(0).Y())/tauandprodTruthplus.at(0).Y(),w);
    		TauPzResPull.at(t).Fill((tauandprodplus.at(0).Z()-tauandprodTruthplus.at(0).Z())/tauandprodTruthplus.at(0).Z(),w);
    		TauEResPull.at(t).Fill((tauandprodplus.at(0).E()-tauandprodTruthplus.at(0).E())/tauandprodTruthplus.at(0).E(),w);
    		TauPxResPull.at(t).Fill((tauandprodminus.at(0).X()-tauandprodTruthminus.at(0).X())/tauandprodTruthminus.at(0).X(),w);
    		TauPyResPull.at(t).Fill((tauandprodminus.at(0).Y()-tauandprodTruthminus.at(0).Y())/tauandprodTruthminus.at(0).Y(),w);
    		TauPzResPull.at(t).Fill((tauandprodminus.at(0).Z()-tauandprodTruthminus.at(0).Z())/tauandprodTruthminus.at(0).Z(),w);
    		TauEResPull.at(t).Fill((tauandprodminus.at(0).E()-tauandprodTruthminus.at(0).E())/tauandprodTruthminus.at(0).E(),w);

    		TauNoBSPxResPull.at(t).Fill((tauandprodNoBSplus.at(0).X()-tauandprodTruthplus.at(0).X())/tauandprodTruthplus.at(0).X(),w);
    		TauNoBSPyResPull.at(t).Fill((tauandprodNoBSplus.at(0).Y()-tauandprodTruthplus.at(0).Y())/tauandprodTruthplus.at(0).Y(),w);
    		TauNoBSPzResPull.at(t).Fill((tauandprodNoBSplus.at(0).Z()-tauandprodTruthplus.at(0).Z())/tauandprodTruthplus.at(0).Z(),w);
    		TauNoBSEResPull.at(t).Fill((tauandprodNoBSplus.at(0).E()-tauandprodTruthplus.at(0).E())/tauandprodTruthplus.at(0).E(),w);
    		TauNoBSPxResPull.at(t).Fill((tauandprodNoBSminus.at(0).X()-tauandprodTruthminus.at(0).X())/tauandprodTruthminus.at(0).X(),w);
    		TauNoBSPyResPull.at(t).Fill((tauandprodNoBSminus.at(0).Y()-tauandprodTruthminus.at(0).Y())/tauandprodTruthminus.at(0).Y(),w);
    		TauNoBSPzResPull.at(t).Fill((tauandprodNoBSminus.at(0).Z()-tauandprodTruthminus.at(0).Z())/tauandprodTruthminus.at(0).Z(),w);
    		TauNoBSEResPull.at(t).Fill((tauandprodNoBSminus.at(0).E()-tauandprodTruthminus.at(0).E())/tauandprodTruthminus.at(0).E(),w);

    		TauNoBSRefittedPionsPxResPull.at(t).Fill((tauandprodNoBSRefittedPionsplus.at(0).X()-tauandprodTruthplus.at(0).X())/tauandprodTruthplus.at(0).X(),w);
    		TauNoBSRefittedPionsPyResPull.at(t).Fill((tauandprodNoBSRefittedPionsplus.at(0).Y()-tauandprodTruthplus.at(0).Y())/tauandprodTruthplus.at(0).Y(),w);
    		TauNoBSRefittedPionsPzResPull.at(t).Fill((tauandprodNoBSRefittedPionsplus.at(0).Z()-tauandprodTruthplus.at(0).Z())/tauandprodTruthplus.at(0).Z(),w);
    		TauNoBSRefittedPionsEResPull.at(t).Fill((tauandprodNoBSRefittedPionsplus.at(0).E()-tauandprodTruthplus.at(0).E())/tauandprodTruthplus.at(0).E(),w);
    		TauNoBSRefittedPionsPxResPull.at(t).Fill((tauandprodNoBSRefittedPionsminus.at(0).X()-tauandprodTruthminus.at(0).X())/tauandprodTruthminus.at(0).X(),w);
    		TauNoBSRefittedPionsPyResPull.at(t).Fill((tauandprodNoBSRefittedPionsminus.at(0).Y()-tauandprodTruthminus.at(0).Y())/tauandprodTruthminus.at(0).Y(),w);
    		TauNoBSRefittedPionsPzResPull.at(t).Fill((tauandprodNoBSRefittedPionsminus.at(0).Z()-tauandprodTruthminus.at(0).Z())/tauandprodTruthminus.at(0).Z(),w);
    		TauNoBSRefittedPionsEResPull.at(t).Fill((tauandprodNoBSRefittedPionsminus.at(0).E()-tauandprodTruthminus.at(0).E())/tauandprodTruthminus.at(0).E(),w);

    		TauOldPxResPull.at(t).Fill((tauandprodOldplus.at(0).X()-tauandprodTruthplus.at(0).X())/tauandprodTruthplus.at(0).X(),w);
    		TauOldPyResPull.at(t).Fill((tauandprodOldplus.at(0).Y()-tauandprodTruthplus.at(0).Y())/tauandprodTruthplus.at(0).Y(),w);
    		TauOldPzResPull.at(t).Fill((tauandprodOldplus.at(0).Z()-tauandprodTruthplus.at(0).Z())/tauandprodTruthplus.at(0).Z(),w);
    		TauOldEResPull.at(t).Fill((tauandprodOldplus.at(0).E()-tauandprodTruthplus.at(0).E())/tauandprodTruthplus.at(0).E(),w);
    		TauOldPxResPull.at(t).Fill((tauandprodOldminus.at(0).X()-tauandprodTruthminus.at(0).X())/tauandprodTruthminus.at(0).X(),w);
    		TauOldPyResPull.at(t).Fill((tauandprodOldminus.at(0).Y()-tauandprodTruthminus.at(0).Y())/tauandprodTruthminus.at(0).Y(),w);
    		TauOldPzResPull.at(t).Fill((tauandprodOldminus.at(0).Z()-tauandprodTruthminus.at(0).Z())/tauandprodTruthminus.at(0).Z(),w);
    		TauOldEResPull.at(t).Fill((tauandprodOldminus.at(0).E()-tauandprodTruthminus.at(0).E())/tauandprodTruthminus.at(0).E(),w);
		
    		TauOldRefittedPionsPxResPull.at(t).Fill((tauandprodOldRefittedPionsplus.at(0).X()-tauandprodTruthplus.at(0).X())/tauandprodTruthplus.at(0).X(),w);
    		TauOldRefittedPionsPyResPull.at(t).Fill((tauandprodOldRefittedPionsplus.at(0).Y()-tauandprodTruthplus.at(0).Y())/tauandprodTruthplus.at(0).Y(),w);
    		TauOldRefittedPionsPzResPull.at(t).Fill((tauandprodOldRefittedPionsplus.at(0).Z()-tauandprodTruthplus.at(0).Z())/tauandprodTruthplus.at(0).Z(),w);
    		TauOldRefittedPionsEResPull.at(t).Fill((tauandprodOldRefittedPionsplus.at(0).E()-tauandprodTruthplus.at(0).E())/tauandprodTruthplus.at(0).E(),w);
    		TauOldRefittedPionsPxResPull.at(t).Fill((tauandprodOldRefittedPionsminus.at(0).X()-tauandprodTruthminus.at(0).X())/tauandprodTruthminus.at(0).X(),w);
    		TauOldRefittedPionsPyResPull.at(t).Fill((tauandprodOldRefittedPionsminus.at(0).Y()-tauandprodTruthminus.at(0).Y())/tauandprodTruthminus.at(0).Y(),w);
    		TauOldRefittedPionsPzResPull.at(t).Fill((tauandprodOldRefittedPionsminus.at(0).Z()-tauandprodTruthminus.at(0).Z())/tauandprodTruthminus.at(0).Z(),w);
    		TauOldRefittedPionsEResPull.at(t).Fill((tauandprodOldRefittedPionsminus.at(0).E()-tauandprodTruthminus.at(0).E())/tauandprodTruthminus.at(0).E(),w);


    		TauTIPPxResPull.at(t).Fill((tauandprodTIPplus.at(0).X()-tauandprodTruthplus.at(0).X())/tauandprodTruthplus.at(0).X(),w);
    		TauTIPPyResPull.at(t).Fill((tauandprodTIPplus.at(0).Y()-tauandprodTruthplus.at(0).Y())/tauandprodTruthplus.at(0).Y(),w);
    		TauTIPPzResPull.at(t).Fill((tauandprodTIPplus.at(0).Z()-tauandprodTruthplus.at(0).Z())/tauandprodTruthplus.at(0).Z(),w);
    		TauTIPEResPull.at(t).Fill((tauandprodTIPplus.at(0).E()-tauandprodTruthplus.at(0).E())/tauandprodTruthplus.at(0).E(),w);
    		TauTIPPxResPull.at(t).Fill((tauandprodTIPminus.at(0).X()-tauandprodTruthminus.at(0).X())/tauandprodTruthminus.at(0).X(),w);
    		TauTIPPyResPull.at(t).Fill((tauandprodTIPminus.at(0).Y()-tauandprodTruthminus.at(0).Y())/tauandprodTruthminus.at(0).Y(),w);
    		TauTIPPzResPull.at(t).Fill((tauandprodTIPminus.at(0).Z()-tauandprodTruthminus.at(0).Z())/tauandprodTruthminus.at(0).Z(),w);
    		TauTIPEResPull.at(t).Fill((tauandprodTIPminus.at(0).E()-tauandprodTruthminus.at(0).E())/tauandprodTruthminus.at(0).E(),w);
		
    		TauTIPRefittedPionsPxResPull.at(t).Fill((tauandprodTIPRefittedPionsplus.at(0).X()-tauandprodTruthplus.at(0).X())/tauandprodTruthplus.at(0).X(),w);
    		TauTIPRefittedPionsPyResPull.at(t).Fill((tauandprodTIPRefittedPionsplus.at(0).Y()-tauandprodTruthplus.at(0).Y())/tauandprodTruthplus.at(0).Y(),w);
    		TauTIPRefittedPionsPzResPull.at(t).Fill((tauandprodTIPRefittedPionsplus.at(0).Z()-tauandprodTruthplus.at(0).Z())/tauandprodTruthplus.at(0).Z(),w);
    		TauTIPRefittedPionsEResPull.at(t).Fill((tauandprodTIPRefittedPionsplus.at(0).E()-tauandprodTruthplus.at(0).E())/tauandprodTruthplus.at(0).E(),w);
    		TauTIPRefittedPionsPxResPull.at(t).Fill((tauandprodTIPRefittedPionsminus.at(0).X()-tauandprodTruthminus.at(0).X())/tauandprodTruthminus.at(0).X(),w);
    		TauTIPRefittedPionsPyResPull.at(t).Fill((tauandprodTIPRefittedPionsminus.at(0).Y()-tauandprodTruthminus.at(0).Y())/tauandprodTruthminus.at(0).Y(),w);
    		TauTIPRefittedPionsPzResPull.at(t).Fill((tauandprodTIPRefittedPionsminus.at(0).Z()-tauandprodTruthminus.at(0).Z())/tauandprodTruthminus.at(0).Z(),w);
    		TauTIPRefittedPionsEResPull.at(t).Fill((tauandprodTIPRefittedPionsminus.at(0).E()-tauandprodTruthminus.at(0).E())/tauandprodTruthminus.at(0).E(),w);

		
    		VectPolaSVFitPxResPull.at(t).Fill((h1SVFit.X()-h1Truth.X())/h1Truth.X(),w);
    		VectPolaSVFitPyResPull.at(t).Fill((h1SVFit.Y()-h1Truth.Y())/h1Truth.Y(),w);
    		VectPolaSVFitPzResPull.at(t).Fill((h1SVFit.Z()-h1Truth.Z())/h1Truth.Z(),w);
    		VectPolaSVFitPxResPull.at(t).Fill((h2SVFit.X()-h2Truth.X())/h2Truth.X(),w);
    		VectPolaSVFitPyResPull.at(t).Fill((h2SVFit.Y()-h2Truth.Y())/h2Truth.Y(),w);
    		VectPolaSVFitPzResPull.at(t).Fill((h2SVFit.Z()-h2Truth.Z())/h2Truth.Z(),w);

    		VectPolaPxResPull.at(t).Fill((h1.X()-h1Truth.X())/h1Truth.X(),w);
    		VectPolaPyResPull.at(t).Fill((h1.Y()-h1Truth.Y())/h1Truth.Y(),w);
    		VectPolaPzResPull.at(t).Fill((h1.Z()-h1Truth.Z())/h1Truth.Z(),w);
    		VectPolaPxResPull.at(t).Fill((h2.X()-h2Truth.X())/h2Truth.X(),w);
    		VectPolaPyResPull.at(t).Fill((h2.Y()-h2Truth.Y())/h2Truth.Y(),w);
    		VectPolaPzResPull.at(t).Fill((h2.Z()-h2Truth.Z())/h2Truth.Z(),w);
    
		
    		TLorentzVector LeadingPionMinusTruth=Pions1Truth.at(0), LeadingPionPlusTruth=Pions2Truth.at(0),LeadingPionMinusSVFit=PionsSVFit1.at(0), LeadingPionPlusSVFit=PionsSVFit2.at(0),LeadingPionMinusReco=HadPions_minus.at(0), LeadingPionPlusReco=HadPions_plus.at(0), LeadingRefittedPionMinus=HadRefittedPions_minus.at(0), LeadingRefittedPionPlus=HadRefittedPions_plus.at(0);
		
    		for(int i=1;i<3;i++)
    		  {
    		    if(LeadingPionMinusTruth.Pt()<Pions1Truth.at(i).Pt())LeadingPionMinusTruth=Pions1Truth.at(i);
    		  }
    		for(int i=1;i<3;i++)
    		  {
    		    if(LeadingPionPlusTruth.Pt()<Pions2Truth.at(i).Pt())LeadingPionPlusTruth=Pions2Truth.at(i);
    		  }
    		for(int i=1;i<3;i++)
    		  {
    		    if(LeadingPionMinusSVFit.Pt()<PionsSVFit1.at(i).Pt())LeadingPionMinusSVFit=PionsSVFit1.at(i);
    		  }

    		for(int i=1;i<3;i++)
    		  {
    		    if(LeadingPionPlusSVFit.Pt()<PionsSVFit2.at(i).Pt())LeadingPionPlusSVFit=PionsSVFit2.at(i);
    		  }

    		for(int i=1;i<3;i++)
    		  {
    		    if(LeadingPionMinusReco.Pt()<HadPions_minus.at(i).Pt())LeadingPionMinusReco=HadPions_minus.at(i);
    		  }

    		for(int i=1;i<3;i++)
    		  {
    		    if(LeadingPionPlusReco.Pt()<HadPions_plus.at(i).Pt())LeadingPionPlusReco=HadPions_plus.at(i);
    		  }
    		for(int i=1;i<3;i++)
    		  {
    		    if(LeadingRefittedPionMinus.Pt()<HadRefittedPions_minus.at(i).Pt())LeadingRefittedPionMinus=HadRefittedPions_minus.at(i);
    		  }

    		for(int i=1;i<3;i++)
    		  {
    		    if(LeadingRefittedPionPlus.Pt()<HadRefittedPions_plus.at(i).Pt())LeadingRefittedPionPlus=HadRefittedPions_plus.at(i);
    		  }

    		LeadingPionPxResPull.at(t).Fill((LeadingPionPlusReco.X()-LeadingPionPlusTruth.X())/LeadingPionPlusTruth.X(),w);
    		LeadingPionPyResPull.at(t).Fill((LeadingPionPlusReco.Y()-LeadingPionPlusTruth.Y())/LeadingPionPlusTruth.Y(),w);
    		LeadingPionPzResPull.at(t).Fill((LeadingPionPlusReco.Z()-LeadingPionPlusTruth.Z())/LeadingPionPlusTruth.Z(),w);
    		LeadingPionEResPull.at(t).Fill((LeadingPionPlusReco.E()-LeadingPionPlusTruth.E())/LeadingPionPlusTruth.E(),w);
    		LeadingPionPxResPull.at(t).Fill((LeadingPionMinusReco.X()-LeadingPionMinusTruth.X())/LeadingPionMinusTruth.X(),w);
    		LeadingPionPyResPull.at(t).Fill((LeadingPionMinusReco.Y()-LeadingPionMinusTruth.Y())/LeadingPionMinusTruth.Y(),w);
    		LeadingPionPzResPull.at(t).Fill((LeadingPionMinusReco.Z()-LeadingPionMinusTruth.Z())/LeadingPionMinusTruth.Z(),w);
    		LeadingPionEResPull.at(t).Fill((LeadingPionMinusReco.E()-LeadingPionMinusTruth.E())/LeadingPionMinusTruth.E(),w);
    		LeadingRefittedPionPxResPull.at(t).Fill((LeadingRefittedPionPlus.X()-LeadingPionPlusTruth.X())/LeadingPionPlusTruth.X(),w);
    		LeadingRefittedPionPyResPull.at(t).Fill((LeadingRefittedPionPlus.Y()-LeadingPionPlusTruth.Y())/LeadingPionPlusTruth.Y(),w);
    		LeadingRefittedPionPzResPull.at(t).Fill((LeadingRefittedPionPlus.Z()-LeadingPionPlusTruth.Z())/LeadingPionPlusTruth.Z(),w);
    		LeadingRefittedPionEResPull.at(t).Fill((LeadingRefittedPionPlus.E()-LeadingPionPlusTruth.E())/LeadingPionPlusTruth.E(),w);
    		LeadingRefittedPionPxResPull.at(t).Fill((LeadingRefittedPionMinus.X()-LeadingPionMinusTruth.X())/LeadingPionMinusTruth.X(),w);
    		LeadingRefittedPionPyResPull.at(t).Fill((LeadingRefittedPionMinus.Y()-LeadingPionMinusTruth.Y())/LeadingPionMinusTruth.Y(),w);
    		LeadingRefittedPionPzResPull.at(t).Fill((LeadingRefittedPionMinus.Z()-LeadingPionMinusTruth.Z())/LeadingPionMinusTruth.Z(),w);
    		LeadingRefittedPionEResPull.at(t).Fill((LeadingRefittedPionMinus.E()-LeadingPionMinusTruth.E())/LeadingPionMinusTruth.E(),w);
		

    		OppositePionPxResPull.at(t).Fill((HadPions_plus.at(0).X()-Pions2Truth.at(0).X())/Pions2Truth.at(0).X(),w);
    		OppositePionPyResPull.at(t).Fill((HadPions_plus.at(0).Y()-Pions2Truth.at(0).Y())/Pions2Truth.at(0).Y(),w);
    		OppositePionPzResPull.at(t).Fill((HadPions_plus.at(0).Z()-Pions2Truth.at(0).Z())/Pions2Truth.at(0).Z(),w);
    		OppositePionEResPull.at(t).Fill((HadPions_plus.at(0).E()-Pions2Truth.at(0).E())/Pions2Truth.at(0).E(),w);
    		OppositePionPxResPull.at(t).Fill((HadPions_minus.at(0).X()-Pions1Truth.at(0).X())/Pions1Truth.at(0).X(),w);
    		OppositePionPyResPull.at(t).Fill((HadPions_minus.at(0).Y()-Pions1Truth.at(0).Y())/Pions1Truth.at(0).Y(),w);
    		OppositePionPzResPull.at(t).Fill((HadPions_minus.at(0).Z()-Pions1Truth.at(0).Z())/Pions1Truth.at(0).Z(),w);
    		OppositePionEResPull.at(t).Fill((HadPions_minus.at(0).E()-Pions1Truth.at(0).E())/Pions1Truth.at(0).E(),w);

    		OppositeRefittedPionPxResPull.at(t).Fill((HadRefittedPions_plus.at(0).X()-Pions2Truth.at(0).X())/Pions2Truth.at(0).X(),w);
    		OppositeRefittedPionPyResPull.at(t).Fill((HadRefittedPions_plus.at(0).Y()-Pions2Truth.at(0).Y())/Pions2Truth.at(0).Y(),w);
    		OppositeRefittedPionPzResPull.at(t).Fill((HadRefittedPions_plus.at(0).Z()-Pions2Truth.at(0).Z())/Pions2Truth.at(0).Z(),w);
    		OppositeRefittedPionEResPull.at(t).Fill((HadRefittedPions_plus.at(0).E()-Pions2Truth.at(0).E())/Pions2Truth.at(0).E(),w);

    		OppositeRefittedPionPxResPull.at(t).Fill((HadRefittedPions_minus.at(0).X()-Pions1Truth.at(0).X())/Pions1Truth.at(0).X(),w);
    		OppositeRefittedPionPyResPull.at(t).Fill((HadRefittedPions_minus.at(0).Y()-Pions1Truth.at(0).Y())/Pions1Truth.at(0).Y(),w);
    		OppositeRefittedPionPzResPull.at(t).Fill((HadRefittedPions_minus.at(0).Z()-Pions1Truth.at(0).Z())/Pions1Truth.at(0).Z(),w);
    		OppositeRefittedPionEResPull.at(t).Fill((HadRefittedPions_minus.at(0).E()-Pions1Truth.at(0).E())/Pions1Truth.at(0).E(),w);
		

    		// DiffTauXLeadingPionXSVFit.at(t).Fill(Tauminussvfit.X()-LeadingPionMinusSVFit.X(),w);
    		// DiffTauXLeadingPionXSVFit.at(t).Fill(Tauplussvfit.X()-LeadingPionPlusSVFit.X(),w);
    		// DiffTauYLeadingPionYSVFit.at(t).Fill(Tauminussvfit.Y()-LeadingPionMinusSVFit.Y(),w);
    		// DiffTauYLeadingPionYSVFit.at(t).Fill(Tauplussvfit.Y()-LeadingPionPlusSVFit.Y(),w);
    		// DiffTauZLeadingPionZSVFit.at(t).Fill(Tauminussvfit.Z()-LeadingPionMinusSVFit.Z(),w);
    		// DiffTauZLeadingPionZSVFit.at(t).Fill(Tauplussvfit.Z()-LeadingPionPlusSVFit.Z(),w);
    		// DiffTauELeadingPionESVFit.at(t).Fill(Tauminussvfit.E()-LeadingPionMinusSVFit.E(),w);
    		// DiffTauELeadingPionESVFit.at(t).Fill(Tauplussvfit.E()-LeadingPionPlusSVFit.E(),w);
		

    		// DiffTauXLeadingPionXReco.at(t).Fill(TauminusPairConstraint.X()-LeadingPionMinusReco.X(),w);
    		// DiffTauXLeadingPionXReco.at(t).Fill(TauplusPairConstraint.X()-LeadingPionPlusReco.X(),w);
    		// DiffTauYLeadingPionYReco.at(t).Fill(TauminusPairConstraint.Y()-LeadingPionMinusReco.Y(),w);
    		// DiffTauYLeadingPionYReco.at(t).Fill(TauplusPairConstraint.Y()-LeadingPionPlusReco.Y(),w);
    		// DiffTauZLeadingPionZReco.at(t).Fill(TauminusPairConstraint.Z()-LeadingPionMinusReco.Z(),w);
    		// DiffTauZLeadingPionZReco.at(t).Fill(TauplusPairConstraint.Z()-LeadingPionPlusReco.Z(),w);
    		// DiffTauELeadingPionEReco.at(t).Fill(TauminusPairConstraint.E()-LeadingPionMinusReco.E(),w);
    		// DiffTauELeadingPionEReco.at(t).Fill(TauplusPairConstraint.E()-LeadingPionPlusReco.E(),w);
		
    		// PxtauSVFitVSPxLeadingPion.at(t).Fill(Tauminussvfit.X(),LeadingPionMinusSVFit.X(),w);
    		// PxtauSVFitVSPxLeadingPion.at(t).Fill(Tauplussvfit.X(),LeadingPionPlusSVFit.X(),w);
    		// PytauSVFitVSPyLeadingPion.at(t).Fill(Tauminussvfit.Y(),LeadingPionMinusSVFit.Y(),w);
    		// PytauSVFitVSPyLeadingPion.at(t).Fill(Tauplussvfit.Y(),LeadingPionPlusSVFit.Y(),w);
    		// PztauSVFitVSPzLeadingPion.at(t).Fill(Tauminussvfit.Z(),LeadingPionMinusSVFit.Z(),w);
    		// PztauSVFitVSPzLeadingPion.at(t).Fill(Tauplussvfit.Z(),LeadingPionPlusSVFit.Z(),w);
    		// EtauSVFitVSELeadingPion.at(t).Fill(Tauminussvfit.E(),LeadingPionMinusSVFit.E(),w);
    		// EtauSVFitVSELeadingPion.at(t).Fill(Tauplussvfit.E(),LeadingPionPlusSVFit.E(),w);
		
    		double h1xresolSVFit=(h1SVFit.X()-h1Truth.X())/h1Truth.X();
    		double h2xresolSVFit=(h2SVFit.X()-h2Truth.X())/h2Truth.X();
    		double h1yresolSVFit=(h1SVFit.Y()-h1Truth.Y())/h1Truth.Y();
    		double h2yresolSVFit=(h2SVFit.Y()-h2Truth.Y())/h2Truth.Y();
    		double h1zresolSVFit=(h1SVFit.Z()-h1Truth.Z())/h1Truth.Z();
    		double h2zresolSVFit=(h2SVFit.Z()-h2Truth.Z())/h2Truth.Z();
    		double h1xresolReco=(h1.X()-h1Truth.X())/h1Truth.X();
    		double h2xresolReco=(h2.X()-h2Truth.X())/h2Truth.X();
    		double h1yresolReco=(h1.Y()-h1Truth.Y())/h1Truth.Y();
    		double h2yresolReco=(h2.Y()-h2Truth.Y())/h2Truth.Y();
    		double h1zresolReco=(h1.Z()-h1Truth.Z())/h1Truth.Z();
    		double h2zresolReco=(h2.Z()-h2Truth.Z())/h2Truth.Z();
    		// cout<<"Tauminussvfit.Z()-a1LV_Tauminus.Z(): "<<Tauminussvfit.Z()-a1LV_Tauminus.Z()<<endl;
    		// cout<<"Tauplussvfit.Z()-a1LV_Tauplus.Z(): "<<Tauplussvfit.Z()-a1LV_Tauplus.Z()<<endl;
    		// cout<<"Tauminussvfit.E()-a1LV_Tauminus.E(): "<<Tauminussvfit.E()-a1LV_Tauminus.E()<<endl;
    		// cout<<"Tauplussvfit.E()-a1LV_Tauplus.E(): "<<Tauplussvfit.E()-a1LV_Tauplus.E()<<endl;
		
    		// HxResolSVFitVSDiffPxtauSVFitPxLeadingPion.at(t).Fill(h1xresolSVFit,Tauminussvfit.X()-LeadingPionMinusSVFit.X(),w);
    		// HxResolSVFitVSDiffPxtauSVFitPxLeadingPion.at(t).Fill(h2xresolSVFit,Tauplussvfit.X()-LeadingPionPlusSVFit.X(),w);
    		// HxResolSVFitVSDiffPytauSVFitPyLeadingPion.at(t).Fill(h1xresolSVFit,Tauminussvfit.Y()-LeadingPionMinusSVFit.Y(),w);
    		// HxResolSVFitVSDiffPytauSVFitPyLeadingPion.at(t).Fill(h2xresolSVFit,Tauplussvfit.Y()-LeadingPionPlusSVFit.Y(),w);
    		// HxResolSVFitVSDiffPztauSVFitPzLeadingPion.at(t).Fill(h1xresolSVFit,Tauminussvfit.Z()-LeadingPionMinusSVFit.Z(),w);
    		// HxResolSVFitVSDiffPztauSVFitPzLeadingPion.at(t).Fill(h2xresolSVFit,Tauplussvfit.Z()-LeadingPionPlusSVFit.Z(),w);
    		// HxResolSVFitVSDiffEtauSVFitELeadingPion.at(t).Fill(h1xresolSVFit,Tauminussvfit.E()-LeadingPionMinusSVFit.E(),w);
    		// HxResolSVFitVSDiffEtauSVFitELeadingPion.at(t).Fill(h2xresolSVFit,Tauplussvfit.E()-LeadingPionPlusSVFit.E(),w);

    		// HyResolSVFitVSDiffPxtauSVFitPxLeadingPion.at(t).Fill(h1yresolSVFit,Tauminussvfit.X()-LeadingPionMinusSVFit.X(),w);
    		// HyResolSVFitVSDiffPxtauSVFitPxLeadingPion.at(t).Fill(h2yresolSVFit,Tauplussvfit.X()-LeadingPionPlusSVFit.X(),w);
    		// HyResolSVFitVSDiffPytauSVFitPyLeadingPion.at(t).Fill(h1yresolSVFit,Tauminussvfit.Y()-LeadingPionMinusSVFit.Y(),w);
    		// HyResolSVFitVSDiffPytauSVFitPyLeadingPion.at(t).Fill(h2yresolSVFit,Tauplussvfit.Y()-LeadingPionPlusSVFit.Y(),w);
    		// HyResolSVFitVSDiffPztauSVFitPzLeadingPion.at(t).Fill(h1yresolSVFit,Tauminussvfit.Z()-LeadingPionMinusSVFit.Z(),w);
    		// HyResolSVFitVSDiffPztauSVFitPzLeadingPion.at(t).Fill(h2yresolSVFit,Tauplussvfit.Z()-LeadingPionPlusSVFit.Z(),w);
    		// HyResolSVFitVSDiffEtauSVFitELeadingPion.at(t).Fill(h1yresolSVFit,Tauminussvfit.E()-LeadingPionMinusSVFit.E(),w);
    		// HyResolSVFitVSDiffEtauSVFitELeadingPion.at(t).Fill(h2yresolSVFit,Tauplussvfit.E()-LeadingPionPlusSVFit.E(),w);

    		// HzResolSVFitVSDiffPxtauSVFitPxLeadingPion.at(t).Fill(h1zresolSVFit,Tauminussvfit.X()-LeadingPionMinusSVFit.X(),w);
    		// HzResolSVFitVSDiffPxtauSVFitPxLeadingPion.at(t).Fill(h2zresolSVFit,Tauplussvfit.X()-LeadingPionPlusSVFit.X(),w);
    		// HzResolSVFitVSDiffPytauSVFitPyLeadingPion.at(t).Fill(h1zresolSVFit,Tauminussvfit.Y()-LeadingPionMinusSVFit.Y(),w);
    		// HzResolSVFitVSDiffPytauSVFitPyLeadingPion.at(t).Fill(h2zresolSVFit,Tauplussvfit.Y()-LeadingPionPlusSVFit.Y(),w);
    		// HzResolSVFitVSDiffPztauSVFitPzLeadingPion.at(t).Fill(h1zresolSVFit,Tauminussvfit.Z()-LeadingPionMinusSVFit.Z(),w);
    		// HzResolSVFitVSDiffPztauSVFitPzLeadingPion.at(t).Fill(h2zresolSVFit,Tauplussvfit.Z()-LeadingPionPlusSVFit.Z(),w);
    		// HzResolSVFitVSDiffEtauSVFitELeadingPion.at(t).Fill(h1zresolSVFit,Tauminussvfit.E()-LeadingPionMinusSVFit.E(),w);
    		// HzResolSVFitVSDiffEtauSVFitELeadingPion.at(t).Fill(h2zresolSVFit,Tauplussvfit.E()-LeadingPionPlusSVFit.E(),w);
		
    		// PxtauSVFitVSPxOppositePion.at(t).Fill(Tauminussvfit.X(),PionsSVFit1.at(0).X(),w);
    		// PxtauSVFitVSPxOppositePion.at(t).Fill(Tauplussvfit.X(),PionsSVFit2.at(0).X(),w);
    		// PytauSVFitVSPyOppositePion.at(t).Fill(Tauminussvfit.Y(),PionsSVFit1.at(0).Y(),w);
    		// PytauSVFitVSPyOppositePion.at(t).Fill(Tauplussvfit.Y(),PionsSVFit2.at(0).Y(),w);
    		// PztauSVFitVSPzOppositePion.at(t).Fill(Tauminussvfit.Z(),PionsSVFit1.at(0).Z(),w);
    		// PztauSVFitVSPzOppositePion.at(t).Fill(Tauplussvfit.Z(),PionsSVFit2.at(0).Z(),w);
    		// EtauSVFitVSEOppositePion.at(t).Fill(Tauminussvfit.E(),PionsSVFit1.at(0).E(),w);
    		// EtauSVFitVSEOppositePion.at(t).Fill(Tauplussvfit.E(),PionsSVFit2.at(0).E(),w);
		
		
    		// HxResolSVFitVSDiffPxtauSVFitPxOppositePion.at(t).Fill(h1xresolSVFit,Tauminussvfit.X()-PionsSVFit1.at(0).X(),w);
    		// HxResolSVFitVSDiffPxtauSVFitPxOppositePion.at(t).Fill(h2xresolSVFit,Tauplussvfit.X()-PionsSVFit2.at(0).X(),w);
    		// HxResolSVFitVSDiffPytauSVFitPyOppositePion.at(t).Fill(h1xresolSVFit,Tauminussvfit.Y()-PionsSVFit1.at(0).Y(),w);
    		// HxResolSVFitVSDiffPytauSVFitPyOppositePion.at(t).Fill(h2xresolSVFit,Tauplussvfit.Y()-PionsSVFit2.at(0).Y(),w);
    		// HxResolSVFitVSDiffPztauSVFitPzOppositePion.at(t).Fill(h1xresolSVFit,Tauminussvfit.Z()-PionsSVFit1.at(0).Z(),w);
    		// HxResolSVFitVSDiffPztauSVFitPzOppositePion.at(t).Fill(h2xresolSVFit,Tauplussvfit.Z()-PionsSVFit2.at(0).Z(),w);
    		// HxResolSVFitVSDiffEtauSVFitEOppositePion.at(t).Fill(h1xresolSVFit,Tauminussvfit.E()-PionsSVFit1.at(0).E(),w);
    		// HxResolSVFitVSDiffEtauSVFitEOppositePion.at(t).Fill(h2xresolSVFit,Tauplussvfit.E()-PionsSVFit2.at(0).E(),w);
		
		
    		// PxtauRecoVSPxLeadingPion.at(t).Fill(TauminusPairConstraint.X(),LeadingPionMinusReco.X(),w);
    		// PxtauRecoVSPxLeadingPion.at(t).Fill(TauplusPairConstraint.X(),LeadingPionPlusReco.X(),w);
    		// PytauRecoVSPyLeadingPion.at(t).Fill(TauminusPairConstraint.Y(),LeadingPionMinusReco.Y(),w);
    		// PytauRecoVSPyLeadingPion.at(t).Fill(TauplusPairConstraint.Y(),LeadingPionPlusReco.Y(),w);
    		// PztauRecoVSPzLeadingPion.at(t).Fill(TauminusPairConstraint.Z(),LeadingPionMinusReco.Z(),w);
    		// PztauRecoVSPzLeadingPion.at(t).Fill(TauplusPairConstraint.Z(),LeadingPionPlusReco.Z(),w);
    		// EtauRecoVSELeadingPion.at(t).Fill(TauminusPairConstraint.E(),LeadingPionMinusReco.E(),w);
    		// EtauRecoVSELeadingPion.at(t).Fill(TauplusPairConstraint.E(),LeadingPionPlusReco.E(),w);
		
	
		
    		// HxResolRecoVSDiffPxtauRecoPxLeadingPion.at(t).Fill(h1xresolReco,TauminusPairConstraint.X()-LeadingPionMinusReco.X(),w);
    		// HxResolRecoVSDiffPxtauRecoPxLeadingPion.at(t).Fill(h2xresolReco,TauplusPairConstraint.X()-LeadingPionPlusReco.X(),w);
    		// HxResolRecoVSDiffPytauRecoPyLeadingPion.at(t).Fill(h1xresolReco,TauminusPairConstraint.Y()-LeadingPionMinusReco.Y(),w);
    		// HxResolRecoVSDiffPytauRecoPyLeadingPion.at(t).Fill(h2xresolReco,TauplusPairConstraint.Y()-LeadingPionPlusReco.Y(),w);
    		// HxResolRecoVSDiffPztauRecoPzLeadingPion.at(t).Fill(h1xresolReco,TauminusPairConstraint.Z()-LeadingPionMinusReco.Z(),w);
    		// HxResolRecoVSDiffPztauRecoPzLeadingPion.at(t).Fill(h2xresolReco,TauplusPairConstraint.Z()-LeadingPionPlusReco.Z(),w);
    		// HxResolRecoVSDiffEtauRecoELeadingPion.at(t).Fill(h1xresolReco,TauminusPairConstraint.E()-LeadingPionMinusReco.E(),w);
    		// HxResolRecoVSDiffEtauRecoELeadingPion.at(t).Fill(h2xresolReco,TauplusPairConstraint.E()-LeadingPionPlusReco.E(),w);

    		// HyResolRecoVSDiffPxtauRecoPxLeadingPion.at(t).Fill(h1yresolReco,TauminusPairConstraint.X()-LeadingPionMinusReco.X(),w);
    		// HyResolRecoVSDiffPxtauRecoPxLeadingPion.at(t).Fill(h2yresolReco,TauplusPairConstraint.X()-LeadingPionPlusReco.X(),w);
    		// HyResolRecoVSDiffPytauRecoPyLeadingPion.at(t).Fill(h1yresolReco,TauminusPairConstraint.Y()-LeadingPionMinusReco.Y(),w);
    		// HyResolRecoVSDiffPytauRecoPyLeadingPion.at(t).Fill(h2yresolReco,TauplusPairConstraint.Y()-LeadingPionPlusReco.Y(),w);
    		// HyResolRecoVSDiffPztauRecoPzLeadingPion.at(t).Fill(h1yresolReco,TauminusPairConstraint.Z()-LeadingPionMinusReco.Z(),w);
    		// HyResolRecoVSDiffPztauRecoPzLeadingPion.at(t).Fill(h2yresolReco,TauplusPairConstraint.Z()-LeadingPionPlusReco.Z(),w);
    		// HyResolRecoVSDiffEtauRecoELeadingPion.at(t).Fill(h1yresolReco,TauminusPairConstraint.E()-LeadingPionMinusReco.E(),w);
    		// HyResolRecoVSDiffEtauRecoELeadingPion.at(t).Fill(h2yresolReco,TauplusPairConstraint.E()-LeadingPionPlusReco.E(),w);

    		// HzResolRecoVSDiffPxtauRecoPxLeadingPion.at(t).Fill(h1zresolReco,TauminusPairConstraint.X()-LeadingPionMinusReco.X(),w);
    		// HzResolRecoVSDiffPxtauRecoPxLeadingPion.at(t).Fill(h2zresolReco,TauplusPairConstraint.X()-LeadingPionPlusReco.X(),w);
    		// HzResolRecoVSDiffPytauRecoPyLeadingPion.at(t).Fill(h1zresolReco,TauminusPairConstraint.Y()-LeadingPionMinusReco.Y(),w);
    		// HzResolRecoVSDiffPytauRecoPyLeadingPion.at(t).Fill(h2zresolReco,TauplusPairConstraint.Y()-LeadingPionPlusReco.Y(),w);
    		// HzResolRecoVSDiffPztauRecoPzLeadingPion.at(t).Fill(h1zresolReco,TauminusPairConstraint.Z()-LeadingPionMinusReco.Z(),w);
    		// HzResolRecoVSDiffPztauRecoPzLeadingPion.at(t).Fill(h2zresolReco,TauplusPairConstraint.Z()-LeadingPionPlusReco.Z(),w);
    		// HzResolRecoVSDiffEtauRecoELeadingPion.at(t).Fill(h1zresolReco,TauminusPairConstraint.E()-LeadingPionMinusReco.E(),w);
    		// HzResolRecoVSDiffEtauRecoELeadingPion.at(t).Fill(h2zresolReco,TauplusPairConstraint.E()-LeadingPionPlusReco.E(),w);
		
    		// PxtauRecoVSPxOppositePion.at(t).Fill(TauminusPairConstraint.X(),HadPions_minus.at(0).X(),w);
    		// PxtauRecoVSPxOppositePion.at(t).Fill(TauplusPairConstraint.X(),HadPions_plus.at(0).X(),w);
    		// PytauRecoVSPyOppositePion.at(t).Fill(TauminusPairConstraint.Y(),HadPions_minus.at(0).Y(),w);
    		// PytauRecoVSPyOppositePion.at(t).Fill(TauplusPairConstraint.Y(),HadPions_plus.at(0).Y(),w);
    		// PztauRecoVSPzOppositePion.at(t).Fill(TauminusPairConstraint.Z(),HadPions_minus.at(0).Z(),w);
    		// PztauRecoVSPzOppositePion.at(t).Fill(TauplusPairConstraint.Z(),HadPions_plus.at(0).Z(),w);
    		// EtauRecoVSEOppositePion.at(t).Fill(TauminusPairConstraint.E(),HadPions_minus.at(0).E(),w);
    		// EtauRecoVSEOppositePion.at(t).Fill(TauplusPairConstraint.E(),HadPions_plus.at(0).E(),w);
		
		
    		// HxResolRecoVSDiffPxtauRecoPxOppositePion.at(t).Fill(h1xresolReco,TauminusPairConstraint.X()-HadPions_minus.at(0).X(),w);
    		// HxResolRecoVSDiffPxtauRecoPxOppositePion.at(t).Fill(h2xresolReco,TauplusPairConstraint.X()-HadPions_plus.at(0).X(),w);
    		// HxResolRecoVSDiffPytauRecoPyOppositePion.at(t).Fill(h1xresolReco,TauminusPairConstraint.Y()-HadPions_minus.at(0).Y(),w);
    		// HxResolRecoVSDiffPytauRecoPyOppositePion.at(t).Fill(h2xresolReco,TauplusPairConstraint.Y()-HadPions_plus.at(0).Y(),w);
    		// HxResolRecoVSDiffPztauRecoPzOppositePion.at(t).Fill(h1xresolReco,TauminusPairConstraint.Z()-HadPions_minus.at(0).Z(),w);
    		// HxResolRecoVSDiffPztauRecoPzOppositePion.at(t).Fill(h2xresolReco,TauplusPairConstraint.Z()-HadPions_plus.at(0).Z(),w);
    		// HxResolRecoVSDiffEtauRecoEOppositePion.at(t).Fill(h1xresolReco,TauminusPairConstraint.E()-HadPions_minus.at(0).E(),w);
    		// HxResolRecoVSDiffEtauRecoEOppositePion.at(t).Fill(h2xresolReco,TauplusPairConstraint.E()-HadPions_plus.at(0).E(),w);


    		// HxResolRecoVSTauPt.at(t).Fill(h1xresolReco,TauminusPairConstraint.Pt(),w);
    		// HxResolRecoVSTauPt.at(t).Fill(h2xresolReco,TauplusPairConstraint.Pt(),w);
    		// HxResolRecoVSPxTau.at(t).Fill(h1xresolReco,TauminusPairConstraint.X(),w);
    		// HxResolRecoVSPxTau.at(t).Fill(h2xresolReco,TauplusPairConstraint.X(),w);
    		// HxResolRecoVSPyTau.at(t).Fill(h1xresolReco,TauminusPairConstraint.Y(),w);
    		// HxResolRecoVSPyTau.at(t).Fill(h2xresolReco,TauplusPairConstraint.Y(),w);
    		// HxResolRecoVSPzTau.at(t).Fill(h1xresolReco,TauminusPairConstraint.Z(),w);
    		// HxResolRecoVSPzTau.at(t).Fill(h2xresolReco,TauplusPairConstraint.Z(),w);
    		// HxResolRecoVSETau.at(t).Fill(h1xresolReco,TauminusPairConstraint.E(),w);
    		// HxResolRecoVSETau.at(t).Fill(h2xresolReco,TauplusPairConstraint.E(),w);
    		// HxResolRecoVSTauEta.at(t).Fill(h1xresolReco,TauminusPairConstraint.Eta(),w);
    		// HxResolRecoVSTauEta.at(t).Fill(h2xresolReco,TauplusPairConstraint.Eta(),w);
    		// HxResolRecoVSGJAngle.at(t).Fill(h1xresolReco,GetThetaGJ(TauminusPairConstraint.Vect(),a1LV_Tauminus),w);
    		// HxResolRecoVSGJAngle.at(t).Fill(h2xresolReco,GetThetaGJ(TauplusPairConstraint.Vect(),a1LV_Tauplus),w);
    		// HxResolRecoVSMH.at(t).Fill(h1xresolReco,higgsmass,w);
    		// HxResolRecoVSMH.at(t).Fill(h2xresolReco,higgsmass,w);
    		// HxResolRecoVSSumPtPions.at(t).Fill(h1xresolReco,a1LV_Tauminus.Pt(),w);
    		// HxResolRecoVSSumPtPions.at(t).Fill(h2xresolReco,a1LV_Tauplus.Pt(),w);
    		// HxResolRecoVSSumPxPions.at(t).Fill(h1xresolReco,a1LV_Tauminus.X(),w);
    		// HxResolRecoVSSumPxPions.at(t).Fill(h2xresolReco,a1LV_Tauplus.X(),w);
    		// HxResolRecoVSSumPyPions.at(t).Fill(h1xresolReco,a1LV_Tauminus.Y(),w);
    		// HxResolRecoVSSumPyPions.at(t).Fill(h2xresolReco,a1LV_Tauplus.Y(),w);
    		// HxResolRecoVSSumPzPions.at(t).Fill(h1xresolReco,a1LV_Tauminus.Z(),w);
    		// HxResolRecoVSSumPzPions.at(t).Fill(h2xresolReco,a1LV_Tauplus.Z(),w);
    		// HxResolRecoVSSumEPions.at(t).Fill(h1xresolReco,a1LV_Tauminus.E(),w);
    		// HxResolRecoVSSumEPions.at(t).Fill(h2xresolReco,a1LV_Tauplus.E(),w);

    		// HyResolRecoVSTauPt.at(t).Fill(h1yresolReco,TauminusPairConstraint.Pt(),w);
    		// HyResolRecoVSTauPt.at(t).Fill(h2yresolReco,TauplusPairConstraint.Pt(),w);
    		// HyResolRecoVSPxTau.at(t).Fill(h1yresolReco,TauminusPairConstraint.X(),w);
    		// HyResolRecoVSPxTau.at(t).Fill(h2yresolReco,TauplusPairConstraint.X(),w);
    		// HyResolRecoVSPyTau.at(t).Fill(h1yresolReco,TauminusPairConstraint.Y(),w);
    		// HyResolRecoVSPyTau.at(t).Fill(h2yresolReco,TauplusPairConstraint.Y(),w);
    		// HyResolRecoVSPzTau.at(t).Fill(h1yresolReco,TauminusPairConstraint.Z(),w);
    		// HyResolRecoVSPzTau.at(t).Fill(h2yresolReco,TauplusPairConstraint.Z(),w);
    		// HyResolRecoVSETau.at(t).Fill(h1yresolReco,TauminusPairConstraint.E(),w);
    		// HyResolRecoVSETau.at(t).Fill(h2yresolReco,TauplusPairConstraint.E(),w);
    		// HyResolRecoVSTauEta.at(t).Fill(h1yresolReco,TauminusPairConstraint.Eta(),w);
    		// HyResolRecoVSTauEta.at(t).Fill(h2yresolReco,TauplusPairConstraint.Eta(),w);
    		// HyResolRecoVSGJAngle.at(t).Fill(h1yresolReco,GetThetaGJ(TauminusPairConstraint.Vect(),a1LV_Tauminus),w);
    		// HyResolRecoVSGJAngle.at(t).Fill(h2yresolReco,GetThetaGJ(TauplusPairConstraint.Vect(),a1LV_Tauplus),w);
    		// HyResolRecoVSMH.at(t).Fill(h1yresolReco,higgsmass,w);
    		// HyResolRecoVSMH.at(t).Fill(h2yresolReco,higgsmass,w);
    		// HyResolRecoVSSumPtPions.at(t).Fill(h1yresolReco,a1LV_Tauminus.Pt(),w);
    		// HyResolRecoVSSumPtPions.at(t).Fill(h2yresolReco,a1LV_Tauplus.Pt(),w);
    		// HyResolRecoVSSumPxPions.at(t).Fill(h1yresolReco,a1LV_Tauminus.X(),w);
    		// HyResolRecoVSSumPxPions.at(t).Fill(h2yresolReco,a1LV_Tauplus.X(),w);
    		// HyResolRecoVSSumPyPions.at(t).Fill(h1yresolReco,a1LV_Tauminus.Y(),w);
    		// HyResolRecoVSSumPyPions.at(t).Fill(h2yresolReco,a1LV_Tauplus.Y(),w);
    		// HyResolRecoVSSumPzPions.at(t).Fill(h1yresolReco,a1LV_Tauminus.Z(),w);
    		// HyResolRecoVSSumPzPions.at(t).Fill(h2yresolReco,a1LV_Tauplus.Z(),w);
    		// HyResolRecoVSSumEPions.at(t).Fill(h1yresolReco,a1LV_Tauminus.E(),w);
    		// HyResolRecoVSSumEPions.at(t).Fill(h2yresolReco,a1LV_Tauplus.E(),w);

    		// HzResolRecoVSTauPt.at(t).Fill(h1zresolReco,TauminusPairConstraint.Pt(),w);
    		// HzResolRecoVSTauPt.at(t).Fill(h2zresolReco,TauplusPairConstraint.Pt(),w);
    		// HzResolRecoVSPxTau.at(t).Fill(h1zresolReco,TauminusPairConstraint.X(),w);
    		// HzResolRecoVSPxTau.at(t).Fill(h2zresolReco,TauplusPairConstraint.X(),w);
    		// HzResolRecoVSPyTau.at(t).Fill(h1zresolReco,TauminusPairConstraint.Y(),w);
    		// HzResolRecoVSPyTau.at(t).Fill(h2zresolReco,TauplusPairConstraint.Y(),w);
    		// HzResolRecoVSPzTau.at(t).Fill(h1zresolReco,TauminusPairConstraint.Z(),w);
    		// HzResolRecoVSPzTau.at(t).Fill(h2zresolReco,TauplusPairConstraint.Z(),w);
    		// HzResolRecoVSETau.at(t).Fill(h1zresolReco,TauminusPairConstraint.E(),w);
    		// HzResolRecoVSETau.at(t).Fill(h2zresolReco,TauplusPairConstraint.E(),w);
    		// HzResolRecoVSTauEta.at(t).Fill(h1zresolReco,TauminusPairConstraint.Eta(),w);
    		// HzResolRecoVSTauEta.at(t).Fill(h2zresolReco,TauplusPairConstraint.Eta(),w);
    		// HzResolRecoVSGJAngle.at(t).Fill(h1zresolReco,GetThetaGJ(TauminusPairConstraint.Vect(),a1LV_Tauminus),w);
    		// HzResolRecoVSGJAngle.at(t).Fill(h2zresolReco,GetThetaGJ(TauplusPairConstraint.Vect(),a1LV_Tauplus),w);
    		// HzResolRecoVSMH.at(t).Fill(h1zresolReco,higgsmass,w);
    		// HzResolRecoVSMH.at(t).Fill(h2zresolReco,higgsmass,w);
    		// HzResolRecoVSSumPtPions.at(t).Fill(h1zresolReco,a1LV_Tauminus.Pt(),w);
    		// HzResolRecoVSSumPtPions.at(t).Fill(h2zresolReco,a1LV_Tauplus.Pt(),w);
    		// HzResolRecoVSSumPxPions.at(t).Fill(h1zresolReco,a1LV_Tauminus.X(),w);
    		// HzResolRecoVSSumPxPions.at(t).Fill(h2zresolReco,a1LV_Tauplus.X(),w);
    		// HzResolRecoVSSumPyPions.at(t).Fill(h1zresolReco,a1LV_Tauminus.Y(),w);
    		// HzResolRecoVSSumPyPions.at(t).Fill(h2zresolReco,a1LV_Tauplus.Y(),w);
    		// HzResolRecoVSSumPzPions.at(t).Fill(h1zresolReco,a1LV_Tauminus.Z(),w);
    		// HzResolRecoVSSumPzPions.at(t).Fill(h2zresolReco,a1LV_Tauplus.Z(),w);
    		// HzResolRecoVSSumEPions.at(t).Fill(h1zresolReco,a1LV_Tauminus.E(),w);
    		// HzResolRecoVSSumEPions.at(t).Fill(h2zresolReco,a1LV_Tauplus.E(),w);
		

    		// HxResolSVFitVSTauPt.at(t).Fill(h1xresolSVFit,Tauminussvfit.Pt(),w);
    		// HxResolSVFitVSTauPt.at(t).Fill(h2xresolSVFit,Tauplussvfit.Pt(),w);
    		// HxResolSVFitVSPxTau.at(t).Fill(h1xresolSVFit,Tauminussvfit.X(),w);
    		// HxResolSVFitVSPxTau.at(t).Fill(h2xresolSVFit,Tauplussvfit.X(),w);
    		// HxResolSVFitVSPyTau.at(t).Fill(h1xresolSVFit,Tauminussvfit.Y(),w);
    		// HxResolSVFitVSPyTau.at(t).Fill(h2xresolSVFit,Tauplussvfit.Y(),w);
    		// HxResolSVFitVSPzTau.at(t).Fill(h1xresolSVFit,Tauminussvfit.Z(),w);
    		// HxResolSVFitVSPzTau.at(t).Fill(h2xresolSVFit,Tauplussvfit.Z(),w);
    		// HxResolSVFitVSETau.at(t).Fill(h1xresolSVFit,Tauminussvfit.E(),w);
    		// HxResolSVFitVSETau.at(t).Fill(h2xresolSVFit,Tauplussvfit.E(),w);
    		// HxResolSVFitVSTauEta.at(t).Fill(h1xresolSVFit,Tauminussvfit.Eta(),w);
    		// HxResolSVFitVSTauEta.at(t).Fill(h2xresolSVFit,Tauplussvfit.Eta(),w);
    		// HxResolSVFitVSGJAngle.at(t).Fill(h1xresolSVFit,GetThetaGJ(Tauminussvfit.Vect(),a1LV_Tauminus),w);
    		// HxResolSVFitVSGJAngle.at(t).Fill(h2xresolSVFit,GetThetaGJ(Tauplussvfit.Vect(),a1LV_Tauplus),w);
    		// HxResolSVFitVSMH.at(t).Fill(h1xresolSVFit,higgsmass,w);
    		// HxResolSVFitVSMH.at(t).Fill(h2xresolSVFit,higgsmass,w);
    		// HxResolSVFitVSSumPtPions.at(t).Fill(h1xresolSVFit,a1LV_Tauminus.Pt(),w);
    		// HxResolSVFitVSSumPtPions.at(t).Fill(h2xresolSVFit,a1LV_Tauplus.Pt(),w);
    		// HxResolSVFitVSSumPxPions.at(t).Fill(h1xresolSVFit,a1LV_Tauminus.X(),w);
    		// HxResolSVFitVSSumPxPions.at(t).Fill(h2xresolSVFit,a1LV_Tauplus.X(),w);
    		// HxResolSVFitVSSumPyPions.at(t).Fill(h1xresolSVFit,a1LV_Tauminus.Y(),w);
    		// HxResolSVFitVSSumPyPions.at(t).Fill(h2xresolSVFit,a1LV_Tauplus.Y(),w);
    		// HxResolSVFitVSSumPzPions.at(t).Fill(h1xresolSVFit,a1LV_Tauminus.Z(),w);
    		// HxResolSVFitVSSumPzPions.at(t).Fill(h2xresolSVFit,a1LV_Tauplus.Z(),w);
    		// HxResolSVFitVSSumEPions.at(t).Fill(h1xresolSVFit,a1LV_Tauminus.E(),w);
    		// HxResolSVFitVSSumEPions.at(t).Fill(h2xresolSVFit,a1LV_Tauplus.E(),w);

    		// HyResolSVFitVSTauPt.at(t).Fill(h1yresolSVFit,Tauminussvfit.Pt(),w);
    		// HyResolSVFitVSTauPt.at(t).Fill(h2yresolSVFit,Tauplussvfit.Pt(),w);
    		// HyResolSVFitVSPxTau.at(t).Fill(h1yresolSVFit,Tauminussvfit.X(),w);
    		// HyResolSVFitVSPxTau.at(t).Fill(h2yresolSVFit,Tauplussvfit.X(),w);
    		// HyResolSVFitVSPyTau.at(t).Fill(h1yresolSVFit,Tauminussvfit.Y(),w);
    		// HyResolSVFitVSPyTau.at(t).Fill(h2yresolSVFit,Tauplussvfit.Y(),w);
    		// HyResolSVFitVSPzTau.at(t).Fill(h1yresolSVFit,Tauminussvfit.Z(),w);
    		// HyResolSVFitVSPzTau.at(t).Fill(h2yresolSVFit,Tauplussvfit.Z(),w);
    		// HyResolSVFitVSETau.at(t).Fill(h1yresolSVFit,Tauminussvfit.E(),w);
    		// HyResolSVFitVSETau.at(t).Fill(h2yresolSVFit,Tauplussvfit.E(),w);
    		// HyResolSVFitVSTauEta.at(t).Fill(h1yresolSVFit,Tauminussvfit.Eta(),w);
    		// HyResolSVFitVSTauEta.at(t).Fill(h2yresolSVFit,Tauplussvfit.Eta(),w);
    		// HyResolSVFitVSGJAngle.at(t).Fill(h1yresolSVFit,GetThetaGJ(Tauminussvfit.Vect(),a1LV_Tauminus),w);
    		// HyResolSVFitVSGJAngle.at(t).Fill(h2yresolSVFit,GetThetaGJ(Tauplussvfit.Vect(),a1LV_Tauplus),w);
    		// HyResolSVFitVSMH.at(t).Fill(h1yresolSVFit,higgsmass,w);
    		// HyResolSVFitVSMH.at(t).Fill(h2yresolSVFit,higgsmass,w);
    		// HyResolSVFitVSSumPtPions.at(t).Fill(h1yresolSVFit,a1LV_Tauminus.Pt(),w);
    		// HyResolSVFitVSSumPtPions.at(t).Fill(h2yresolSVFit,a1LV_Tauplus.Pt(),w);
    		// HyResolSVFitVSSumPxPions.at(t).Fill(h1yresolSVFit,a1LV_Tauminus.X(),w);
    		// HyResolSVFitVSSumPxPions.at(t).Fill(h2yresolSVFit,a1LV_Tauplus.X(),w);
    		// HyResolSVFitVSSumPyPions.at(t).Fill(h1yresolSVFit,a1LV_Tauminus.Y(),w);
    		// HyResolSVFitVSSumPyPions.at(t).Fill(h2yresolSVFit,a1LV_Tauplus.Y(),w);
    		// HyResolSVFitVSSumPzPions.at(t).Fill(h1yresolSVFit,a1LV_Tauminus.Z(),w);
    		// HyResolSVFitVSSumPzPions.at(t).Fill(h2yresolSVFit,a1LV_Tauplus.Z(),w);
    		// HyResolSVFitVSSumEPions.at(t).Fill(h1yresolSVFit,a1LV_Tauminus.E(),w);
    		// HyResolSVFitVSSumEPions.at(t).Fill(h2yresolSVFit,a1LV_Tauplus.E(),w);

    		// HzResolSVFitVSTauPt.at(t).Fill(h1zresolSVFit,Tauminussvfit.Pt(),w);
    		// HzResolSVFitVSTauPt.at(t).Fill(h2zresolSVFit,Tauplussvfit.Pt(),w);
    		// HzResolSVFitVSPxTau.at(t).Fill(h1zresolSVFit,Tauminussvfit.X(),w);
    		// HzResolSVFitVSPxTau.at(t).Fill(h2zresolSVFit,Tauplussvfit.X(),w);
    		// HzResolSVFitVSPyTau.at(t).Fill(h1zresolSVFit,Tauminussvfit.Y(),w);
    		// HzResolSVFitVSPyTau.at(t).Fill(h2zresolSVFit,Tauplussvfit.Y(),w);
    		// HzResolSVFitVSPzTau.at(t).Fill(h1zresolSVFit,Tauminussvfit.Z(),w);
    		// HzResolSVFitVSPzTau.at(t).Fill(h2zresolSVFit,Tauplussvfit.Z(),w);
    		// HzResolSVFitVSETau.at(t).Fill(h1zresolSVFit,Tauminussvfit.E(),w);
    		// HzResolSVFitVSETau.at(t).Fill(h2zresolSVFit,Tauplussvfit.E(),w);
    		// HzResolSVFitVSTauEta.at(t).Fill(h1zresolSVFit,Tauminussvfit.Eta(),w);
    		// HzResolSVFitVSTauEta.at(t).Fill(h2zresolSVFit,Tauplussvfit.Eta(),w);
    		// HzResolSVFitVSGJAngle.at(t).Fill(h1zresolSVFit,GetThetaGJ(Tauminussvfit.Vect(),a1LV_Tauminus),w);
    		// HzResolSVFitVSGJAngle.at(t).Fill(h2zresolSVFit,GetThetaGJ(Tauplussvfit.Vect(),a1LV_Tauplus),w);
    		// HzResolSVFitVSMH.at(t).Fill(h1zresolSVFit,higgsmass,w);
    		// HzResolSVFitVSMH.at(t).Fill(h2zresolSVFit,higgsmass,w);
    		// HzResolSVFitVSSumPtPions.at(t).Fill(h1zresolSVFit,a1LV_Tauminus.Pt(),w);
    		// HzResolSVFitVSSumPtPions.at(t).Fill(h2zresolSVFit,a1LV_Tauplus.Pt(),w);
    		// HzResolSVFitVSSumPxPions.at(t).Fill(h1zresolSVFit,a1LV_Tauminus.X(),w);
    		// HzResolSVFitVSSumPxPions.at(t).Fill(h2zresolSVFit,a1LV_Tauplus.X(),w);
    		// HzResolSVFitVSSumPyPions.at(t).Fill(h1zresolSVFit,a1LV_Tauminus.Y(),w);
    		// HzResolSVFitVSSumPyPions.at(t).Fill(h2zresolSVFit,a1LV_Tauplus.Y(),w);
    		// HzResolSVFitVSSumPzPions.at(t).Fill(h1zresolSVFit,a1LV_Tauminus.Z(),w);
    		// HzResolSVFitVSSumPzPions.at(t).Fill(h2zresolSVFit,a1LV_Tauplus.Z(),w);
    		// HzResolSVFitVSSumEPions.at(t).Fill(h1zresolSVFit,a1LV_Tauminus.E(),w);
    		// HzResolSVFitVSSumEPions.at(t).Fill(h2zresolSVFit,a1LV_Tauplus.E(),w);

		
    		// if(((h1SVFit).Cross(h2SVFit))*tauSVFitminus_HRF.Vect().Unit()<=0)
    		//   {
    		// 	AcopAngleVSTauPt.at(t).Fill(TMath::ATan2((k1SVFit.Cross(k2SVFit)).Mag(),k1SVFit*k2SVFit),tauandprodSVFitminus.at(0).Pt(),w);
    		// 	AcopAngleVSTauPt.at(t).Fill(TMath::ATan2((k1SVFit.Cross(k2SVFit)).Mag(),k1SVFit*k2SVFit),tauandprodSVFitplus.at(0).Pt(),w);
    		// 	AcopAngleVSTauEta.at(t).Fill(TMath::ATan2((k1SVFit.Cross(k2SVFit)).Mag(),k1SVFit*k2SVFit),tauandprodSVFitminus.at(0).Eta(),w);
    		// 	AcopAngleVSTauEta.at(t).Fill(TMath::ATan2((k1SVFit.Cross(k2SVFit)).Mag(),k1SVFit*k2SVFit),tauandprodSVFitplus.at(0).Eta(),w);
    		//   }
    		// else 
    		//   {
    		//     AcopAngleVSTauPt.at(t).Fill(2.*TMath::Pi()-TMath::ATan2((k1SVFit.Cross(k2SVFit)).Mag(),k1SVFit*k2SVFit),tauandprodSVFitminus.at(0).Pt(),w);
    		//     AcopAngleVSTauPt.at(t).Fill(2.*TMath::Pi()-TMath::ATan2((k1SVFit.Cross(k2SVFit)).Mag(),k1SVFit*k2SVFit),tauandprodSVFitplus.at(0).Pt(),w);
    		//     AcopAngleVSTauEta.at(t).Fill(2.*TMath::Pi()-TMath::ATan2((k1SVFit.Cross(k2SVFit)).Mag(),k1SVFit*k2SVFit),tauandprodSVFitminus.at(0).Eta(),w);
    		//     AcopAngleVSTauEta.at(t).Fill(2.*TMath::Pi()-TMath::ATan2((k1SVFit.Cross(k2SVFit)).Mag(),k1SVFit*k2SVFit),tauandprodSVFitplus.at(0).Eta(),w);
    		//   }
		 

		
		
    		// TVector3 hmoins(-0.889958,-0.447263,-0.089056);
    		// TVector3 taumoins(55.375824,26.615817,4.943766);
    		// //hmoins.Unit().Cross(taumoins.Unit()).Print();
    		// //hmoins=-hmoins;
    		// //hmoins.Unit().Cross(taumoins.Unit()).Print();
    		// TVector3 hplus(0.910199,0.408791,0.066536);
    		// TVector3 tauplus(-55.375824,-26.615817,-4.943766);
    		// TVector3 kmoins=hmoins.Unit().Cross(taumoins.Unit());
    		// TVector3 kplus=hplus.Unit().Cross(tauplus.Unit());
    		// if(((hmoins.Unit()).Cross(hplus.Unit()))*taumoins.Unit()<=0){cout<<"Result1: "<<TMath::ATan2((kmoins.Unit().Cross(kplus)).Mag(),kmoins*kplus)<<endl;}
    		// else {cout<<"Result1: "<<2.*TMath::Pi()-TMath::ATan2((kmoins.Unit().Cross(kplus)).Mag(),kmoins*kplus)<<endl;}
    		// hmoins=-hmoins;
    		// if(((hmoins.Unit()).Cross(hplus.Unit()))*taumoins.Unit()<=0){cout<<"Result2: "<<TMath::ATan2((kmoins.Unit().Cross(kplus)).Mag(),kmoins*kplus)<<endl;}
		
    		// if((((h1SVFit.X()-h1Truth.X())/h1Truth.X())>-2.5) && (((h1SVFit.X()-h1Truth.X())/h1Truth.X())<-1.5))
    		//   {
		    
    		// cout<<"Truth Reco SVFit"<<endl;
    		// cout<<"P-: "<<endl;
    		// tauandprodTruthminus.at(0).Print();
    		// tauandprodminus.at(0).Print();
    		// tauandprodSVFitminus.at(0).Print();
    		// cout<<"P+: "<<endl;
    		// tauandprodTruthplus.at(0).Print();
    		// tauandprodplus.at(0).Print();
    		// tauandprodSVFitplus.at(0).Print();
    		// cout<<"P- boosted: "<<endl;
    		// tauminusTruth_HRF.Print();
    		// tauminus_HRF.Print();
    		// tauSVFitminus_HRF.Print();
    		// cout<<"P+ boosted: "<<endl;
    		// tauplusTruth_HRF.Print();
    		// tauplus_HRF.Print();
    		// tauSVFitplus_HRF.Print();
    		// cout<<"Pions-: "<<endl;
    		// tauandprodTruthminus.at(1).Print(); tauandprodTruthminus.at(2).Print(); tauandprodTruthminus.at(3).Print();
    		// tauandprodminus.at(1).Print(); tauandprodminus.at(2).Print(); tauandprodminus.at(3).Print();
    		// tauandprodSVFitminus.at(1).Print(); tauandprodSVFitminus.at(2).Print(); tauandprodSVFitminus.at(3).Print();
    		// if(hasNoBS){HadPionsNoBS_minus.at(0).Print();HadPionsNoBS_minus.at(1).Print();HadPionsNoBS_minus.at(2).Print();}
    		// cout<<"Pions+: "<<endl;
    		// tauandprodTruthplus.at(1).Print(); tauandprodTruthplus.at(2).Print(); tauandprodTruthplus.at(3).Print();
    		// tauandprodplus.at(1).Print(); tauandprodplus.at(2).Print(); tauandprodplus.at(3).Print();
    		// tauandprodSVFitplus.at(1).Print(); tauandprodSVFitplus.at(2).Print(); tauandprodSVFitplus.at(3).Print();
    		// if(hasNoBS){HadPionsNoBS_plus.at(0).Print();HadPionsNoBS_plus.at(1).Print();HadPionsNoBS_plus.at(2).Print();}
    		// cout<<"Pions- charges"<<endl;
    		// cout<<Pions1TruthCharge.at(0)<<" "<<Pions1TruthCharge.at(1)<<" "<<Pions1TruthCharge.at(2)<<endl;
    		// cout<<HadPionsCharge_minus.at(0)<<" "<<HadPionsCharge_minus.at(1)<<" "<<HadPionsCharge_minus.at(2)<<endl;
    		// cout<<PionsSVFit1Charge.at(0)<<" "<<PionsSVFit1Charge.at(1)<<" "<<PionsSVFit1Charge.at(2)<<endl;
    		// //cout<<Ntp->PFTauRefit_PionsCharge(Tauminus,0)<<" "<<Ntp->PFTauRefit_PionsCharge(Tauminus,1)<<" "<<Ntp->PFTauRefit_PionsCharge(Tauminus,2)<<endl;
    		// cout<<"Pions+ charges"<<endl;
    		// cout<<Pions2TruthCharge.at(0)<<" "<<Pions2TruthCharge.at(1)<<" "<<Pions2TruthCharge.at(2)<<endl;
    		// cout<<HadPionsCharge_plus.at(0)<<" "<<HadPionsCharge_plus.at(1)<<" "<<HadPionsCharge_plus.at(2)<<endl;
    		// cout<<PionsSVFit2Charge.at(0)<<" "<<PionsSVFit2Charge.at(1)<<" "<<PionsSVFit2Charge.at(2)<<endl;
    		// //cout<<Ntp->PFTauRefit_PionsCharge(Tauplus,0)<<" "<<Ntp->PFTauRefit_PionsCharge(Tauplus,1)<<" "<<Ntp->PFTauRefit_PionsCharge(Tauplus,2)<<endl;
    		// cout<<"h-"<<endl;
    		// h1Truth.Print();
    		// h1.Print();
    		// h1SVFit.Print();
    		// cout<<"h+"<<endl;
    		// h2Truth.Print();
    		// h2.Print();
    		// h2SVFit.Print();
    		// cout<<"k-"<<endl;
    		// k1Truth.Print();
    		// k1.Print();
    		// k1SVFit.Print();
    		// cout<<"k+"<<endl;
    		// k2Truth.Print();
    		// k2.Print();
    		// k2SVFit.Print();
    		// cout<<"Angle"<<endl;
    		// if(((h1Truth*h1TruthNorm).Cross(h2Truth*h2TruthNorm))*tauminusTruth_HRF.Vect().Unit()<=0)cout<<TMath::ATan2((k1Truth.Cross(k2Truth)).Mag(),k1Truth*k2Truth)<<endl;
    		// else cout<<2*TMath::Pi()-TMath::ATan2((k1Truth.Cross(k2Truth)).Mag(),k1Truth*k2Truth)<<endl;
    		// if(((h1*h1Norm).Cross(h2*h2Norm))*tauminus_HRF.Vect().Unit()<=0)cout<<TMath::ATan2((k1.Cross(k2)).Mag(),k1*k2)<<endl;
    		// else cout<<2*TMath::Pi()-TMath::ATan2((k1.Cross(k2)).Mag(),k1*k2)<<endl;
    		// if(((h1SVFit*h1SVFitNorm).Cross(h2SVFit*h2SVFitNorm))*tauSVFitminus_HRF.Vect().Unit()<=0)cout<<TMath::ATan2((k1SVFit.Cross(k2SVFit)).Mag(),k1SVFit*k2SVFit)<<endl;
    		// else cout<<2*TMath::Pi()-TMath::ATan2((k1SVFit.Cross(k2SVFit)).Mag(),k1SVFit*k2SVFit)<<endl;
    		// cout<<"TEESSTTT"<<endl;
    		//}
		  }
    	      }
    	  }
      }





    // ClassicSVfit svfitAlgo1;
    // if(a1pi && std::isnan(Wspin)!=true)
    //   {
    // 	// // //---------  svfit ---------------------
    // 	std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;
    // 	classic_svFit::MeasuredTauLepton lep1(1, Tau1P4.Pt(), Tau1P4.Eta(),  Tau1P4.Phi(), Tau1P4.M());
    // 	classic_svFit::MeasuredTauLepton lep2(1, Tau2P4.Pt(), Tau2P4.Eta(),  Tau2P4.Phi(), Tau2P4.M());
	
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
    // 	svfitAlgo1.integrate(measuredTauLeptons,metx,mety, metcov );
    // 	if(svfitAlgo1.isValidSolution()){
	  
    // 	  double higgsmass  = static_cast<classic_svFit::DiTauSystemHistogramAdapter*>(svfitAlgo1.getHistogramAdapter())->getMass();
    // 	  h_SVFitMass.at(t).Fill(higgsmass,w); 
	  
    // 	  tau1P4 = static_cast<classic_svFit::TauTauHistogramAdapter*>(svfitAlgo1.getHistogramAdapter())->GetFittedTau1LV();
    // 	  tau2P4 = static_cast<classic_svFit::TauTauHistogramAdapter*>(svfitAlgo1.getHistogramAdapter())->GetFittedTau2LV();
	
    // 	  // // //---------  svfit ---------------------
	  
    // 	  if(Ntp->Daughters_charge(Tau1)>0)
    // 	    {
    // 	      Tauplussvfit.SetPxPyPzE(tau1P4.Px(),tau1P4.Py(),tau1P4.Pz(),sqrt( pow(1.77685, 2.0) + pow(tau1P4.P(), 2.0) )/*tau1P4.E()*/);
    // 	      Tauminussvfit.SetPxPyPzE(tau2P4.Px(),tau2P4.Py(),tau2P4.Pz(),sqrt( pow(1.77685, 2.0) + pow(tau2P4.P(), 2.0) )/*tau2P4.E()*/);
    // 	    }
    // 	  else
    // 	    {
    // 	      Tauplussvfit.SetPxPyPzE(tau2P4.Px(),tau2P4.Py(),tau2P4.Pz(),sqrt( pow(1.77685, 2.0) + pow(tau2P4.P(), 2.0) )/*tau2P4.E()*/);
    // 	      Tauminussvfit.SetPxPyPzE(tau1P4.Px(),tau1P4.Py(),tau1P4.Pz(),sqrt( pow(1.77685, 2.0) + pow(tau1P4.P(), 2.0) )/*tau1P4.E()*/);
    // 	    }
    // 	  //Tauplussvfit.Print();
    // 	  //cout<<"M: "<<tau1P4.E()*tau1P4.E()-tau1P4.Px()*tau1P4.Px()-tau1P4.Py()*tau1P4.Py()-tau1P4.Pz()*tau1P4.Pz()<<endl;
    // 	}
	
	
    // 	TLorentzVector zeroLV(0,0,0,0);
	  
    // 	TLorentzVector Tau1Truth; 
    // 	TLorentzVector Tau2Truth;
    // 	TLorentzVector TruthDecayFromTau1;
    // 	TLorentzVector TruthDecayFromTau2; 
	
    // 	std::vector<TLorentzVector> Pions1Truth;
    // 	std::vector<TLorentzVector> Pions2Truth;
    // 	std::vector<double> Pions1TruthCharge;
    // 	std::vector<double> Pions2TruthCharge;
    // 	vector<TLorentzVector> HadPions_plus;    
    // 	vector<double> HadPionsCharge_plus;
    // 	vector<TLorentzVector> HadPions_minus;   
    // 	vector<double> HadPionsNoBSCharge_minus;
    // 	vector<TLorentzVector> HadPionsNoBS_plus;    
    // 	vector<double> HadPionsNoBSCharge_plus;
    // 	vector<TLorentzVector> HadPionsNoBS_minus;   
    // 	vector<double> HadPionsCharge_minus;
    // 	std::vector<TLorentzVector> PionsSVFit1;
    // 	std::vector<TLorentzVector> PionsSVFit2;
    // 	std::vector<double> PionsSVFit1Charge;
    // 	std::vector<double> PionsSVFit2Charge;
    
    // 	vector<TLorentzVector> tauandprodTruthminus;
    // 	vector<TLorentzVector> tauandprodTruthplus;
    // 	vector<TLorentzVector> tauandprodminus;
    // 	vector<TLorentzVector> tauandprodplus;
    // 	vector<TLorentzVector> tauandprodSVFitminus;
    // 	vector<TLorentzVector> tauandprodSVFitplus;

    // 	bool passParticlesTruth=false;
    // 	bool passParticlesReco=false;
    // 	bool passParticlesSVFit=false;

    // 	bool passVariablesTruth=false;
    // 	bool passVariablesReco=false;
    // 	bool passVariablesSVFit=false;

	  
    // 	// Truth
	  
    // 	// for(int i=0;i<Ntp->NMCTauDecayProducts(0);i++)
    // 	//   {
    // 	//     cout<<"pdgid0: "<<Ntp->MCTauandProd_pdgid(0,i);
    // 	//   }
    // 	// cout<<endl;
    // 	// for(int i=0;i<Ntp->NMCTauDecayProducts(1);i++)
    // 	//   {
    // 	//     cout<<"pdgid1: "<<Ntp->MCTauandProd_pdgid(1,i);
    // 	//   }
    // 	// cout<<endl;
    // 	SCalculator ScalcTrutha11("a1");
    // 	SCalculator ScalcTruthpion2("pion");
    // 	SCalculator ScalcTruthpion1("pion");
    // 	SCalculator ScalcTrutha12("a1");
    // 	if(a1minus && piplus)
    // 	  {

    // 	    Tau1Truth=Ntp->GetTruthTauLV(5,0);
    // 	    Tau2Truth=Ntp->GetTruthTauLV(3,1);

    // 	    Pions1Truth=Ntp->GetTruthPionsFromA1(0);
    // 	    TruthDecayFromTau1=Pions1Truth.at(0)+Pions1Truth.at(1)+Pions1Truth.at(2);
    // 	    Pions1TruthCharge.push_back(1);
    // 	    Pions1TruthCharge.push_back(-1);
    // 	    Pions1TruthCharge.push_back(-1);

    // 	    Pions2Truth.push_back(Ntp->GetTruthTauProductLV(3,211,0));
    // 	    Pions2TruthCharge.push_back(1);

    // 	    ScalcTrutha11.SortPions(Pions1Truth, Pions1TruthCharge);
    // 	    ScalcTruthpion2.SortPions(Pions2Truth, Pions2TruthCharge);
			  
    // 	    tauandprodTruthminus.push_back(Tau1Truth);
    // 	    tauandprodTruthminus.push_back(Pions1Truth.at(0));
    // 	    tauandprodTruthminus.push_back(Pions1Truth.at(1));
    // 	    tauandprodTruthminus.push_back(Pions1Truth.at(2));
    // 	    tauandprodTruthplus.push_back(Tau2Truth);   
    // 	    tauandprodTruthplus.push_back(Pions2Truth.at(0));     

    // 	    if((Pions1Truth!=VectZeroLV) && (Pions2Truth.at(0)!=zeroLV) && tauandprodTruthminus.at(0)!=zeroLV && tauandprodTruthplus.at(0)!=zeroLV && tauandprodTruthminus.at(0)!=tauandprodTruthplus.at(0)) passParticlesTruth=true;
    // 	  }
    // 	else if(piminus && a1plus)
    // 	  {
    // 	    Tau1Truth=Ntp->GetTruthTauLV(3,0);
    // 	    Tau2Truth=Ntp->GetTruthTauLV(5,1);
      
    // 	    Pions2Truth=Ntp->GetTruthPionsFromA1(1);
    // 	    TruthDecayFromTau2=Pions2Truth.at(0)+Pions2Truth.at(1)+Pions2Truth.at(2);
    // 	    Pions2TruthCharge.push_back(-1);
    // 	    Pions2TruthCharge.push_back(1);
    // 	    Pions2TruthCharge.push_back(1);
      
    // 	    Pions1Truth.push_back(Ntp->GetTruthTauProductLV(3,211,0));
    // 	    Pions1TruthCharge.push_back(-1);
	    
    // 	    ScalcTruthpion1.SortPions(Pions1Truth, Pions1TruthCharge);
    // 	    ScalcTrutha12.SortPions(Pions2Truth, Pions2TruthCharge);
			  
    // 	    tauandprodTruthplus.push_back(Tau2Truth);
    // 	    tauandprodTruthplus.push_back(Pions2Truth.at(0));
    // 	    tauandprodTruthplus.push_back(Pions2Truth.at(1));
    // 	    tauandprodTruthplus.push_back(Pions2Truth.at(2));
    // 	    tauandprodTruthminus.push_back(Tau1Truth);   
    // 	    tauandprodTruthminus.push_back(Pions1Truth.at(0));     

    // 	    if((Pions1Truth.at(0)!=zeroLV) && (Pions2Truth!=VectZeroLV) && tauandprodTruthminus.at(0)!=zeroLV && tauandprodTruthplus.at(0)!=zeroLV && tauandprodTruthminus.at(0)!=tauandprodTruthplus.at(0)) passParticlesTruth=true;
    // 	  }

	
    // 	passParticlesReco=true;

    
    // 	// SVFit
    // 	SCalculator ScalcSVFita11("a1");
    // 	SCalculator ScalcSVFitpion2("pion");
    // 	SCalculator ScalcSVFitpion1("pion");
    // 	SCalculator ScalcSVFita12("a1");
    // 	if (a1minus && piplus)
    // 	  {

    // 	    PionsSVFit1.push_back(Ntp->PFTau_PionsP4(Tauminus,0));
    // 	    PionsSVFit1.push_back(Ntp->PFTau_PionsP4(Tauminus,1));
    // 	    PionsSVFit1.push_back(Ntp->PFTau_PionsP4(Tauminus,2));
	    
    // 	    PionsSVFit1Charge.push_back(Ntp->PFTau_PionsCharge(Tauminus, 0));
    // 	    PionsSVFit1Charge.push_back(Ntp->PFTau_PionsCharge(Tauminus, 1));
    // 	    PionsSVFit1Charge.push_back(Ntp->PFTau_PionsCharge(Tauminus, 2));

    // 	    PionsSVFit2.push_back(Ntp->Daughters_P4(Tauplus));

    // 	    PionsSVFit2Charge.push_back(Ntp->Daughters_charge(Tauplus));
	  

    // 	    ScalcSVFita11.SortPions(PionsSVFit1, PionsSVFit1Charge);
    // 	    ScalcSVFitpion2.SortPions(PionsSVFit2, PionsSVFit2Charge);
	    
	
    // 	    tauandprodSVFitminus.push_back(Tauminussvfit);
    // 	    cout.precision(10);

    // 	    tauandprodSVFitminus.push_back(PionsSVFit1.at(0));
    // 	    tauandprodSVFitminus.push_back(PionsSVFit1.at(1));
    // 	    tauandprodSVFitminus.push_back(PionsSVFit1.at(2));
    // 	    tauandprodSVFitplus.push_back(Tauplussvfit);
    // 	    tauandprodSVFitplus.push_back(PionsSVFit2.at(0));   

    // 	    if((PionsSVFit1!=VectZeroLV) && (PionsSVFit2.at(0)!=zeroLV) && tauandprodSVFitminus.at(0)!=zeroLV && tauandprodSVFitplus.at(0)!=zeroLV && tauandprodSVFitminus.at(0)!=tauandprodSVFitplus.at(0) && svfitAlgo1.isValidSolution()) passParticlesSVFit=true;
		
    // 	  }
    // 	else if (a1plus && piminus)
    // 	  {
    // 	    PionsSVFit2.push_back(Ntp->PFTau_PionsP4(Tauplus,0));
    // 	    PionsSVFit2.push_back(Ntp->PFTau_PionsP4(Tauplus,1));
    // 	    PionsSVFit2.push_back(Ntp->PFTau_PionsP4(Tauplus,2));
	    
    // 	    PionsSVFit2Charge.push_back(Ntp->PFTau_PionsCharge(Tauplus, 0));
    // 	    PionsSVFit2Charge.push_back(Ntp->PFTau_PionsCharge(Tauplus, 1));
    // 	    PionsSVFit2Charge.push_back(Ntp->PFTau_PionsCharge(Tauplus, 2));
	    
    // 	    PionsSVFit1.push_back(Ntp->Daughters_P4(Tauminus));

    // 	    PionsSVFit1Charge.push_back(Ntp->Daughters_charge(Tauminus));
	  

    // 	    ScalcSVFitpion1.SortPions(PionsSVFit1, PionsSVFit1Charge);
    // 	    ScalcSVFita12.SortPions(PionsSVFit2, PionsSVFit2Charge);
	
	    
    // 	    tauandprodSVFitplus.push_back(Tauplussvfit);
    // 	    cout.precision(10);
	    
    // 	    tauandprodSVFitplus.push_back(PionsSVFit2.at(0));
    // 	    tauandprodSVFitplus.push_back(PionsSVFit2.at(1));
    // 	    tauandprodSVFitplus.push_back(PionsSVFit2.at(2));
    // 	    tauandprodSVFitminus.push_back(Tauminussvfit);
    // 	    tauandprodSVFitminus.push_back(PionsSVFit1.at(0));   

	    
    // 	    if((PionsSVFit2!=VectZeroLV) && (PionsSVFit1.at(0)!=zeroLV) && tauandprodSVFitminus.at(0)!=zeroLV && tauandprodSVFitplus.at(0)!=zeroLV && tauandprodSVFitminus.at(0)!=tauandprodSVFitplus.at(0) && svfitAlgo1.isValidSolution()) passParticlesSVFit=true;
    // 	  }
    // 	// cout<<"Nu: "<<endl;
    // 	// ScalcSVFit1.Boost((Tauminussvfit-(Ntp->PFTau_PionsP4(Tauminus,0)+Ntp->PFTau_PionsP4(Tauminus,1)+Ntp->PFTau_PionsP4(Tauminus,2))),Tauminussvfit+Tauplussvfit).Print();
    // 	// cout<<"NuRefit: "<<endl;
    // 	// ScalcSVFit1.Boost((Tauminussvfit-(Ntp->PFTauRefit_PionsP4(Tauminus,0)+Ntp->PFTauRefit_PionsP4(Tauminus,1)+Ntp->PFTauRefit_PionsP4(Tauminus,2))),Tauminussvfit+Tauplussvfit).Print();

	  

    // 	if(passParticlesTruth && passParticlesReco && passParticlesSVFit)
    // 	  {
    // 	    TVector3 h1Truth, h2Truth, k1Truth, k2Truth;
    // 	    TLorentzVector tauminusTruth_HRF,tauplusTruth_HRF;
    // 	    double h1TruthNorm, h2TruthNorm, normk1Truth, normk2Truth;
    // 	    if(a1minus && piplus)
    // 	      {
    // 		// Truth
    // 		cout<<"Truth"<<endl;
    // 		ScalcTrutha11.Configure(tauandprodTruthminus,tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0), -1);
    // 		h1Truth=ScalcTrutha11.pv();
		
    // 		ScalcTruthpion2.Configure(tauandprodTruthplus,tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0), +1);
    // 		h2Truth=ScalcTruthpion2.pv();
		
    // 		tauminusTruth_HRF = ScalcTrutha11.Boost(tauandprodTruthminus.at(0),tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0));
    // 		tauplusTruth_HRF  = ScalcTruthpion2.Boost(tauandprodTruthplus.at(0),tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0));
		
    // 		h1TruthNorm=1./h1Truth.Mag();
    // 		h2TruthNorm=1./h2Truth.Mag();
		
    // 		h1Truth=h1Truth*h1TruthNorm;
    // 		h2Truth=h2Truth*h2TruthNorm;
		
    // 		normk1Truth=1./((h1Truth.Cross(tauminusTruth_HRF.Vect().Unit())).Mag());
    // 		normk2Truth=1./((h2Truth.Cross(tauplusTruth_HRF.Vect().Unit())).Mag());
		
    // 		k1Truth = (h1Truth.Cross(tauminusTruth_HRF.Vect().Unit()))*normk1Truth;
    // 		k2Truth = (h2Truth.Cross(tauplusTruth_HRF.Vect().Unit()))*normk2Truth;
		
    // 		if(std::isnan(h1TruthNorm)!=true && std::isnan(h2TruthNorm)!=true && std::isnan(normk1Truth)!=true && std::isnan(normk2Truth)!=true) passVariablesTruth=true;
		
    // 		passVariablesReco=true;
    // 	      }
    // 	    else if(a1plus && piminus)
    // 	      {
    // 		// Truth
    // 		cout<<"Truth"<<endl;
    // 		ScalcTruthpion1.Configure(tauandprodTruthminus,tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0), -1);
    // 		h1Truth=ScalcTruthpion1.pv();

    // 		ScalcTrutha12.Configure(tauandprodTruthplus,tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0), +1);
    // 		h2Truth=ScalcTrutha12.pv();
		
    // 		tauminusTruth_HRF = ScalcTruthpion1.Boost(tauandprodTruthminus.at(0),tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0));
    // 		tauplusTruth_HRF  = ScalcTrutha12.Boost(tauandprodTruthplus.at(0),tauandprodTruthminus.at(0)+tauandprodTruthplus.at(0));
		
    // 		h1TruthNorm=1./h1Truth.Mag();
    // 		h2TruthNorm=1./h2Truth.Mag();
		
    // 		h1Truth=h1Truth*h1TruthNorm;
    // 		h2Truth=h2Truth*h2TruthNorm;
		
    // 		normk1Truth=1./((h1Truth.Cross(tauminusTruth_HRF.Vect().Unit())).Mag());
    // 		normk2Truth=1./((h2Truth.Cross(tauplusTruth_HRF.Vect().Unit())).Mag());
		
    // 		k1Truth = (h1Truth.Cross(tauminusTruth_HRF.Vect().Unit()))*normk1Truth;
    // 		k2Truth = (h2Truth.Cross(tauplusTruth_HRF.Vect().Unit()))*normk2Truth;
		
    // 		if(std::isnan(h1TruthNorm)!=true && std::isnan(h2TruthNorm)!=true && std::isnan(normk1Truth)!=true && std::isnan(normk2Truth)!=true) passVariablesTruth=true;
		
    // 		passVariablesReco=true;
    // 	      }
    // 	    // SVFit
    // 	    cout<<"SVFit"<<endl;
    // 	    TVector3 h1SVFit, h2SVFit, k1SVFit, k2SVFit;
    // 	    TLorentzVector tauSVFitminus_HRF,tauSVFitplus_HRF;
    // 	    double h1SVFitNorm, h2SVFitNorm, k1SVFitNorm, k2SVFitNorm;
	    
    // 	    if (a1minus && piplus)
    // 	      {
    // 		ScalcSVFita11.Configure(tauandprodSVFitminus,tauandprodSVFitminus.at(0)+tauandprodSVFitplus.at(0), -1);
    // 		h1SVFit=ScalcSVFita11.pv();

    // 		ScalcSVFitpion2.Configure(tauandprodSVFitplus,tauandprodSVFitminus.at(0)+tauandprodSVFitplus.at(0), +1);
    // 		h2SVFit=ScalcSVFitpion2.pv();
	      
    // 		tauSVFitminus_HRF = ScalcSVFita11.Boost(tauandprodSVFitminus.at(0),tauandprodSVFitminus.at(0)+tauandprodSVFitplus.at(0));
    // 		tauSVFitplus_HRF  = ScalcSVFitpion2.Boost(tauandprodSVFitplus.at(0),tauandprodSVFitminus.at(0)+tauandprodSVFitplus.at(0));

    // 		h1SVFitNorm=1./h1SVFit.Mag();
    // 		h2SVFitNorm=1./h2SVFit.Mag();
	      
    // 		h1SVFit=h1SVFit*h1SVFitNorm;
    // 		h2SVFit=h2SVFit*h2SVFitNorm;

    // 		k1SVFitNorm=1./((h1SVFit.Cross(tauSVFitminus_HRF.Vect().Unit())).Mag());
    // 		k2SVFitNorm=1./((h2SVFit.Cross(tauSVFitplus_HRF.Vect().Unit())).Mag());
    // 		k1SVFit = (h1SVFit.Cross(tauSVFitminus_HRF.Vect().Unit()))*k1SVFitNorm;
    // 		k2SVFit = (h2SVFit.Cross(tauSVFitplus_HRF.Vect().Unit()))*k2SVFitNorm;
		    
    // 		if(std::isnan(h1SVFitNorm)!=true && std::isnan(h2SVFitNorm)!=true && std::isnan(k1SVFitNorm)!=true && std::isnan(k2SVFitNorm)!=true) passVariablesSVFit=true;
    // 	      }
    // 	    else if(a1plus && piminus)
    // 	      {
    // 		ScalcSVFitpion1.Configure(tauandprodSVFitminus,tauandprodSVFitminus.at(0)+tauandprodSVFitplus.at(0), -1);
    // 		h1SVFit=ScalcSVFitpion1.pv();

    // 		ScalcSVFita12.Configure(tauandprodSVFitplus,tauandprodSVFitminus.at(0)+tauandprodSVFitplus.at(0), +1);
    // 		h2SVFit=ScalcSVFita12.pv();
	      
    // 		tauSVFitminus_HRF = ScalcSVFitpion1.Boost(tauandprodSVFitminus.at(0),tauandprodSVFitminus.at(0)+tauandprodSVFitplus.at(0));
    // 		tauSVFitplus_HRF  = ScalcSVFita12.Boost(tauandprodSVFitplus.at(0),tauandprodSVFitminus.at(0)+tauandprodSVFitplus.at(0));

    // 		h1SVFitNorm=1./h1SVFit.Mag();
    // 		h2SVFitNorm=1./h2SVFit.Mag();

    // 		h1SVFit=h1SVFit*h1SVFitNorm;
    // 		h2SVFit=h2SVFit*h2SVFitNorm;

    // 		k1SVFitNorm=1./((h1SVFit.Cross(tauSVFitminus_HRF.Vect().Unit())).Mag());
    // 		k2SVFitNorm=1./((h2SVFit.Cross(tauSVFitplus_HRF.Vect().Unit())).Mag());
    // 		k1SVFit = (h1SVFit.Cross(tauSVFitminus_HRF.Vect().Unit()))*k1SVFitNorm;
    // 		k2SVFit = (h2SVFit.Cross(tauSVFitplus_HRF.Vect().Unit()))*k2SVFitNorm;

    // 		if(std::isnan(h1SVFitNorm)!=true && std::isnan(h2SVFitNorm)!=true && std::isnan(k1SVFitNorm)!=true && std::isnan(k2SVFitNorm)!=true) passVariablesSVFit=true;
    // 	      }
    // 	    if(passVariablesTruth && passVariablesReco && passVariablesSVFit)
    // 	      {

    // 		// Truth
		    
    // 		if(((h1Truth).Cross(h2Truth))*tauminusTruth_HRF.Vect().Unit()<=0){polarimetricAcopAngleTruthA1.at(t).Fill(TMath::ATan2((k1Truth.Cross(k2Truth)).Mag(),k1Truth*k2Truth),Wspin);}
    // 		else {polarimetricAcopAngleTruthA1.at(t).Fill(2*TMath::Pi()-TMath::ATan2((k1Truth.Cross(k2Truth)).Mag(),k1Truth*k2Truth),Wspin);}
		    

    // 		// SVFit

    // 		if(((h1SVFit).Cross(h2SVFit))*tauSVFitminus_HRF.Vect().Unit()<=0){polarimetricAcopAngleSVFitA1.at(t).Fill(TMath::ATan2((k1SVFit.Cross(k2SVFit)).Mag(),k1SVFit*k2SVFit),Wspin);}
    // 		else {polarimetricAcopAngleSVFitA1.at(t).Fill(2.*TMath::Pi()-TMath::ATan2((k1SVFit.Cross(k2SVFit)).Mag(),k1SVFit*k2SVFit),Wspin);}		    

    // 		if(a1plus && piminus)
    // 		  {

    // 		    TauminusSVFitPxResPull.at(t).Fill((tauandprodSVFitminus.at(0).X()-tauandprodTruthminus.at(0).X())/tauandprodTruthminus.at(0).X(),w);
    // 		    TauminusSVFitPyResPull.at(t).Fill((tauandprodSVFitminus.at(0).Y()-tauandprodTruthminus.at(0).Y())/tauandprodTruthminus.at(0).Y(),w);
    // 		    TauminusSVFitPzResPull.at(t).Fill((tauandprodSVFitminus.at(0).Z()-tauandprodTruthminus.at(0).Z())/tauandprodTruthminus.at(0).Z(),w);

    // 		    // TauplusSVFitPxResPull.at(t).Fill((tauandprodSVFitplus.at(0).X()-tauandprodTruthplus.at(0).X())/tauandprodTruthplus.at(0).X(),w);
    // 		    // TauplusSVFitPyResPull.at(t).Fill((tauandprodSVFitplus.at(0).Y()-tauandprodTruthplus.at(0).Y())/tauandprodTruthplus.at(0).Y(),w);
    // 		    // TauplusSVFitPzResPull.at(t).Fill((tauandprodSVFitplus.at(0).Z()-tauandprodTruthplus.at(0).Z())/tauandprodTruthplus.at(0).Z(),w);

    // 		    // TauminusPxResPull.at(t).Fill((tauandprodminus.at(0).X()-tauandprodTruthminus.at(0).X())/tauandprodTruthminus.at(0).X(),w);
    // 		    // TauminusPyResPull.at(t).Fill((tauandprodminus.at(0).Y()-tauandprodTruthminus.at(0).Y())/tauandprodTruthminus.at(0).Y(),w);
    // 		    // TauminusPzResPull.at(t).Fill((tauandprodminus.at(0).Z()-tauandprodTruthminus.at(0).Z())/tauandprodTruthminus.at(0).Z(),w);

    // 		    // TauplusPxResPull.at(t).Fill((tauandprodplus.at(0).X()-tauandprodTruthplus.at(0).X())/tauandprodTruthplus.at(0).X(),w);
    // 		    // TauplusPyResPull.at(t).Fill((tauandprodplus.at(0).Y()-tauandprodTruthplus.at(0).Y())/tauandprodTruthplus.at(0).Y(),w);
    // 		    // TauplusPzResPull.at(t).Fill((tauandprodplus.at(0).Z()-tauandprodTruthplus.at(0).Z())/tauandprodTruthplus.at(0).Z(),w);
		

    // 		     VectPolaminusSVFitPxResPull.at(t).Fill((h1SVFit.X()-h1Truth.X())/h1Truth.X(),w);
    // 		     VectPolaminusSVFitPyResPull.at(t).Fill((h1SVFit.Y()-h1Truth.Y())/h1Truth.Y(),w);
    // 		     VectPolaminusSVFitPzResPull.at(t).Fill((h1SVFit.Z()-h1Truth.Z())/h1Truth.Z(),w);
		
    // 		    // VectPolaplusSVFitPxResPull.at(t).Fill((h2SVFit.X()-h2Truth.X())/h2Truth.X(),w);
    // 		    // VectPolaplusSVFitPyResPull.at(t).Fill((h2SVFit.Y()-h2Truth.Y())/h2Truth.Y(),w);
    // 		    // VectPolaplusSVFitPzResPull.at(t).Fill((h2SVFit.Z()-h2Truth.Z())/h2Truth.Z(),w);

    // 		    // VectPolaminusPxResPull.at(t).Fill((h1.X()-h1Truth.X())/h1Truth.X(),w);
    // 		    // VectPolaminusPyResPull.at(t).Fill((h1.Y()-h1Truth.Y())/h1Truth.Y(),w);
    // 		    // VectPolaminusPzResPull.at(t).Fill((h1.Z()-h1Truth.Z())/h1Truth.Z(),w);

    // 		    // VectPolaplusPxResPull.at(t).Fill((h2.X()-h2Truth.X())/h2Truth.X(),w);
    // 		    // VectPolaplusPyResPull.at(t).Fill((h2.Y()-h2Truth.Y())/h2Truth.Y(),w);
    // 		    // VectPolaplusPzResPull.at(t).Fill((h2.Z()-h2Truth.Z())/h2Truth.Z(),w);
    // 		  }


    // 		// if(((h1SVFit).Cross(h2SVFit))*tauSVFitminus_HRF.Vect().Unit()<=0)
    // 		//   {
    // 		//     if(TMath::ATan2((k1SVFit.Cross(k2SVFit)).Mag(),k1SVFit*k2SVFit)>2.1)
    // 		//       {
    // 		// 	BumpPtminus.at(t).Fill(TMath::ATan2((k1SVFit.Cross(k2SVFit)).Mag(),k1SVFit*k2SVFit),tauandprodSVFitminus.at(0).Pt(),w);
    // 		// 	BumpPtplus.at(t).Fill(TMath::ATan2((k1SVFit.Cross(k2SVFit)).Mag(),k1SVFit*k2SVFit),tauandprodSVFitplus.at(0).Pt(),w);
    // 		// 	BumpEtaminus.at(t).Fill(TMath::ATan2((k1SVFit.Cross(k2SVFit)).Mag(),k1SVFit*k2SVFit),tauandprodSVFitminus.at(0).Eta(),w);
    // 		// 	BumpEtaplus.at(t).Fill(TMath::ATan2((k1SVFit.Cross(k2SVFit)).Mag(),k1SVFit*k2SVFit),tauandprodSVFitplus.at(0).Eta(),w);

    // 		//       }
    // 		//   }
    // 		// else 
    // 		//   {
    // 		//     if((2.*TMath::Pi()-TMath::ATan2((k1SVFit.Cross(k2SVFit)).Mag(),k1SVFit*k2SVFit))<4.1)
    // 		//       {
    // 		// 	BumpPtminus.at(t).Fill(2.*TMath::Pi()-TMath::ATan2((k1SVFit.Cross(k2SVFit)).Mag(),k1SVFit*k2SVFit),tauandprodSVFitminus.at(0).Pt(),w);
    // 		// 	BumpPtplus.at(t).Fill(2.*TMath::Pi()-TMath::ATan2((k1SVFit.Cross(k2SVFit)).Mag(),k1SVFit*k2SVFit),tauandprodSVFitplus.at(0).Pt(),w);
    // 		// 	BumpEtaminus.at(t).Fill(2.*TMath::Pi()-TMath::ATan2((k1SVFit.Cross(k2SVFit)).Mag(),k1SVFit*k2SVFit),tauandprodSVFitminus.at(0).Eta(),w);
    // 		// 	BumpEtaplus.at(t).Fill(2.*TMath::Pi()-TMath::ATan2((k1SVFit.Cross(k2SVFit)).Mag(),k1SVFit*k2SVFit),tauandprodSVFitplus.at(0).Eta(),w);

    // 		//       }
    // 		//}
		 

		
		
    // 		// TVector3 hmoins(-0.889958,-0.447263,-0.089056);
    // 		// TVector3 taumoins(55.375824,26.615817,4.943766);
    // 		// //hmoins.Unit().Cross(taumoins.Unit()).Print();
    // 		// //hmoins=-hmoins;
    // 		// //hmoins.Unit().Cross(taumoins.Unit()).Print();
    // 		// TVector3 hplus(0.910199,0.408791,0.066536);
    // 		// TVector3 tauplus(-55.375824,-26.615817,-4.943766);
    // 		// TVector3 kmoins=hmoins.Unit().Cross(taumoins.Unit());
    // 		// TVector3 kplus=hplus.Unit().Cross(tauplus.Unit());
    // 		// if(((hmoins.Unit()).Cross(hplus.Unit()))*taumoins.Unit()<=0){cout<<"Result1: "<<TMath::ATan2((kmoins.Unit().Cross(kplus)).Mag(),kmoins*kplus)<<endl;}
    // 		// else {cout<<"Result1: "<<2.*TMath::Pi()-TMath::ATan2((kmoins.Unit().Cross(kplus)).Mag(),kmoins*kplus)<<endl;}
    // 		// hmoins=-hmoins;
    // 		// if(((hmoins.Unit()).Cross(hplus.Unit()))*taumoins.Unit()<=0){cout<<"Result2: "<<TMath::ATan2((kmoins.Unit().Cross(kplus)).Mag(),kmoins*kplus)<<endl;}
		
    // 		if(a1plus && piminus)
    // 		  {
    // 		    //		    if((((h1SVFit.X()-h1Truth.X())/h1Truth.X())>-2.5) && (((h1SVFit.X()-h1Truth.X())/h1Truth.X())<-1.5))
    // 		    //{

    // 			cout<<"Truth Reco SVFit"<<endl;
    // 			cout<<"P-: "<<endl;
    // 			tauandprodTruthminus.at(0).Print();
		    
    // 			tauandprodSVFitminus.at(0).Print();
    // 			cout<<"P+: "<<endl;
    // 			tauandprodTruthplus.at(0).Print();
		    
    // 			tauandprodSVFitplus.at(0).Print();
    // 			cout<<"P- boosted: "<<endl;
    // 			tauminusTruth_HRF.Print();
		    
    // 			tauSVFitminus_HRF.Print();
    // 			cout<<"P+ boosted: "<<endl;
    // 			tauplusTruth_HRF.Print();
		    
    // 			tauSVFitplus_HRF.Print();

    // 			cout<<"Pions-: "<<endl;
    // 			tauandprodTruthminus.at(1).Print();
		    
    // 			tauandprodSVFitminus.at(1).Print();
		    
    // 			cout<<"Pions+: "<<endl;
    // 			tauandprodTruthplus.at(1).Print(); tauandprodTruthplus.at(2).Print(); tauandprodTruthplus.at(3).Print();
		    
    // 			tauandprodSVFitplus.at(1).Print(); tauandprodSVFitplus.at(2).Print(); tauandprodSVFitplus.at(3).Print();

    // 			cout<<"Pions- charges"<<endl;
    // 			cout<<Pions1TruthCharge.at(0)<<endl;
    // 			cout<<PionsSVFit1Charge.at(0)<<endl;
		  
    // 			cout<<"Pions+ charges"<<endl;
    // 			cout<<Pions2TruthCharge.at(0)<<" "<<Pions2TruthCharge.at(1)<<" "<<Pions2TruthCharge.at(2)<<endl;
    // 			cout<<PionsSVFit2Charge.at(0)<<" "<<PionsSVFit2Charge.at(1)<<" "<<PionsSVFit2Charge.at(2)<<endl;
		  
    // 			cout<<"h-"<<endl;
    // 			h1Truth.Print();
		    
    // 			h1SVFit.Print();
    // 			cout<<"h+"<<endl;
    // 			h2Truth.Print();
		    
    // 			h2SVFit.Print();
    // 			cout<<"k-"<<endl;
    // 			k1Truth.Print();
		    
    // 			k1SVFit.Print();
    // 			cout<<"k+"<<endl;
    // 			k2Truth.Print();

    // 			k2SVFit.Print();
    // 			cout<<"Angle"<<endl;
    // 			if(((h1Truth*h1TruthNorm).Cross(h2Truth*h2TruthNorm))*tauminusTruth_HRF.Vect().Unit()<=0)cout<<TMath::ATan2((k1Truth.Cross(k2Truth)).Mag(),k1Truth*k2Truth)<<endl;
    // 			else cout<<2*TMath::Pi()-TMath::ATan2((k1Truth.Cross(k2Truth)).Mag(),k1Truth*k2Truth)<<endl;
    // 			if(((h1SVFit*h1SVFitNorm).Cross(h2SVFit*h2SVFitNorm))*tauSVFitminus_HRF.Vect().Unit()<=0)cout<<TMath::ATan2((k1SVFit.Cross(k2SVFit)).Mag(),k1SVFit*k2SVFit)<<endl;
    // 			else cout<<2*TMath::Pi()-TMath::ATan2((k1SVFit.Cross(k2SVFit)).Mag(),k1SVFit*k2SVFit)<<endl;
    // 			cout<<"TEESSTTT"<<endl;
    // 			//		      }
    // 		  }
		
    // 	      }
    // 	  }
    //   }
     
  }
}
//  This is a function if you want to do something after the event loop
void Debug::Finish() {
  Selection::Finish();
}
