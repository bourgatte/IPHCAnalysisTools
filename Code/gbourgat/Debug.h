#ifndef Debug_h
#define Debug_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "boost/functional/hash.hpp"
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

class Debug : public Selection {

 public:
  Debug(TString Name_, TString id_);
  virtual ~Debug();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {Id_and_Kin=0, 
	     PairCharge, 
	     NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();
  ReferenceScaleFactors *RSF;
  int TriggerOkDummy, selVertexDummy, selMuon_IsoDummy, selMuon_AntiIsoDummy, selTauDummy, ChargeSumDummy;
  double MTDummy, MvisDummy, TauFLSigmaDummy;

  int Charge;

  PUReweight reweight;
  DataMCCorrections DataMC_Corr;
  tauTrigSFreader tauTrgSF;

 private:
  // Selection Variables and Histos

  std::vector<TH1D> h_SVFitMass; 

  std::vector<TH1D> polarimetricAcopAngle;
  std::vector<TH1D> polarimetricAcopAngleNoBS;
  std::vector<TH1D> polarimetricAcopAngleNoBSPionsRefit;
  std::vector<TH1D> polarimetricAcopAngleTruthA1;
  std::vector<TH1D> polarimetricAcopAngleSVFitA1;  
  std::vector<TH1D> polarimetricAcopAngleOld;
  std::vector<TH1D> polarimetricAcopAngleOldPionsRefit;
  std::vector<TH1D> polarimetricAcopAngleTIP;
  std::vector<TH1D> polarimetricAcopAngleTIPPionsRefit;

  std::vector<TH1D> Unphysical;
  
  std::vector<TH1D> PVXResol;
  std::vector<TH1D> PVXNoBSResol;
  std::vector<TH1D> PVXOldRefitResol;
  std::vector<TH1D> PVXTIPRefitResol;
 
  std::vector<TH1D> PVYResol;
  std::vector<TH1D> PVYNoBSResol;
  std::vector<TH1D> PVYOldRefitResol;
  std::vector<TH1D> PVYTIPRefitResol;
 
  std::vector<TH1D> PVZResol;
  std::vector<TH1D> PVZNoBSResol;
  std::vector<TH1D> PVZOldRefitResol;
  std::vector<TH1D> PVZTIPRefitResol;
  
  std::vector<TH1D> TauSVFitPxResPull;
  std::vector<TH1D> TauSVFitPyResPull;
  std::vector<TH1D> TauSVFitPzResPull;
  std::vector<TH1D> TauSVFitEResPull;

  std::vector<TH1D> TauPxResPull;
  std::vector<TH1D> TauPyResPull;
  std::vector<TH1D> TauPzResPull;
  std::vector<TH1D> TauEResPull;

  std::vector<TH1D> TauNoBSPxResPull;
  std::vector<TH1D> TauNoBSPyResPull;
  std::vector<TH1D> TauNoBSPzResPull;
  std::vector<TH1D> TauNoBSEResPull;

  std::vector<TH1D> TauNoBSRefittedPionsPxResPull;
  std::vector<TH1D> TauNoBSRefittedPionsPyResPull;
  std::vector<TH1D> TauNoBSRefittedPionsPzResPull;
  std::vector<TH1D> TauNoBSRefittedPionsEResPull;
  
  std::vector<TH1D> TauOldPxResPull;
  std::vector<TH1D> TauOldPyResPull;
  std::vector<TH1D> TauOldPzResPull;
  std::vector<TH1D> TauOldEResPull;

  std::vector<TH1D> TauOldRefittedPionsPxResPull;
  std::vector<TH1D> TauOldRefittedPionsPyResPull;
  std::vector<TH1D> TauOldRefittedPionsPzResPull;
  std::vector<TH1D> TauOldRefittedPionsEResPull;

  std::vector<TH1D> TauTIPPxResPull;
  std::vector<TH1D> TauTIPPyResPull;
  std::vector<TH1D> TauTIPPzResPull;
  std::vector<TH1D> TauTIPEResPull;

  std::vector<TH1D> TauTIPRefittedPionsPxResPull;
  std::vector<TH1D> TauTIPRefittedPionsPyResPull;
  std::vector<TH1D> TauTIPRefittedPionsPzResPull;
  std::vector<TH1D> TauTIPRefittedPionsEResPull;

  std::vector<TH1D> LeadingPionPxResPull;
  std::vector<TH1D> LeadingPionPyResPull;
  std::vector<TH1D> LeadingPionPzResPull;
  std::vector<TH1D> LeadingPionEResPull;

  std::vector<TH1D> OppositePionPxResPull;
  std::vector<TH1D> OppositePionPyResPull;
  std::vector<TH1D> OppositePionPzResPull;
  std::vector<TH1D> OppositePionEResPull;

  std::vector<TH1D> LeadingRefittedPionPxResPull;
  std::vector<TH1D> LeadingRefittedPionPyResPull;
  std::vector<TH1D> LeadingRefittedPionPzResPull;
  std::vector<TH1D> LeadingRefittedPionEResPull;

  std::vector<TH1D> OppositeRefittedPionPxResPull;
  std::vector<TH1D> OppositeRefittedPionPyResPull;
  std::vector<TH1D> OppositeRefittedPionPzResPull;
  std::vector<TH1D> OppositeRefittedPionEResPull;
  
  std::vector<TH1D> VectPolaSVFitPxResPull;
  std::vector<TH1D> VectPolaSVFitPyResPull;
  std::vector<TH1D> VectPolaSVFitPzResPull;

  
  std::vector<TH1D> VectPolaPxResPull;
  std::vector<TH1D> VectPolaPyResPull;
  std::vector<TH1D> VectPolaPzResPull;


  std::vector<TH1D> DiffTauXLeadingPionXSVFit;
  std::vector<TH1D> DiffTauYLeadingPionYSVFit;
  std::vector<TH1D> DiffTauZLeadingPionZSVFit;
  std::vector<TH1D> DiffTauELeadingPionESVFit;

  std::vector<TH1D> DiffTauXLeadingPionXReco;
  std::vector<TH1D> DiffTauYLeadingPionYReco;
  std::vector<TH1D> DiffTauZLeadingPionZReco;
  std::vector<TH1D> DiffTauELeadingPionEReco;

  
  std::vector<TH2D> AcopAngleVSTauPt;
  std::vector<TH2D> AcopAngleVSTauEta;
  
  std::vector<TH2D> PxtauSVFitVSPxLeadingPion;
  std::vector<TH2D> PytauSVFitVSPyLeadingPion;
  std::vector<TH2D> PztauSVFitVSPzLeadingPion;
  std::vector<TH2D> EtauSVFitVSELeadingPion;

  std::vector<TH2D> HxResolSVFitVSDiffPxtauSVFitPxLeadingPion;
  std::vector<TH2D> HxResolSVFitVSDiffPytauSVFitPyLeadingPion;
  std::vector<TH2D> HxResolSVFitVSDiffPztauSVFitPzLeadingPion;
  std::vector<TH2D> HxResolSVFitVSDiffEtauSVFitELeadingPion;

  std::vector<TH2D> HyResolSVFitVSDiffPxtauSVFitPxLeadingPion;
  std::vector<TH2D> HyResolSVFitVSDiffPytauSVFitPyLeadingPion;
  std::vector<TH2D> HyResolSVFitVSDiffPztauSVFitPzLeadingPion;
  std::vector<TH2D> HyResolSVFitVSDiffEtauSVFitELeadingPion;

  std::vector<TH2D> HzResolSVFitVSDiffPxtauSVFitPxLeadingPion;
  std::vector<TH2D> HzResolSVFitVSDiffPytauSVFitPyLeadingPion;
  std::vector<TH2D> HzResolSVFitVSDiffPztauSVFitPzLeadingPion;
  std::vector<TH2D> HzResolSVFitVSDiffEtauSVFitELeadingPion;


  /* std::vector<TH2D> PxtauSVFitVSPxOppositePion; */
  /* std::vector<TH2D> PytauSVFitVSPyOppositePion; */
  /* std::vector<TH2D> PztauSVFitVSPzOppositePion; */
  /* std::vector<TH2D> EtauSVFitVSEOppositePion; */

  /* std::vector<TH2D> HxResolSVFitVSDiffPxtauSVFitPxOppositePion; */
  /* std::vector<TH2D> HxResolSVFitVSDiffPytauSVFitPyOppositePion; */
  /* std::vector<TH2D> HxResolSVFitVSDiffPztauSVFitPzOppositePion; */
  /* std::vector<TH2D> HxResolSVFitVSDiffEtauSVFitEOppositePion; */


  std::vector<TH2D> PxtauRecoVSPxLeadingPion;
  std::vector<TH2D> PytauRecoVSPyLeadingPion;
  std::vector<TH2D> PztauRecoVSPzLeadingPion;
  std::vector<TH2D> EtauRecoVSELeadingPion;

  std::vector<TH2D> HxResolRecoVSDiffPxtauRecoPxLeadingPion;
  std::vector<TH2D> HxResolRecoVSDiffPytauRecoPyLeadingPion;
  std::vector<TH2D> HxResolRecoVSDiffPztauRecoPzLeadingPion;
  std::vector<TH2D> HxResolRecoVSDiffEtauRecoELeadingPion;

  std::vector<TH2D> HyResolRecoVSDiffPxtauRecoPxLeadingPion;
  std::vector<TH2D> HyResolRecoVSDiffPytauRecoPyLeadingPion;
  std::vector<TH2D> HyResolRecoVSDiffPztauRecoPzLeadingPion;
  std::vector<TH2D> HyResolRecoVSDiffEtauRecoELeadingPion;

  std::vector<TH2D> HzResolRecoVSDiffPxtauRecoPxLeadingPion;
  std::vector<TH2D> HzResolRecoVSDiffPytauRecoPyLeadingPion;
  std::vector<TH2D> HzResolRecoVSDiffPztauRecoPzLeadingPion;
  std::vector<TH2D> HzResolRecoVSDiffEtauRecoELeadingPion;



  std::vector<TH2D> PxtauRecoVSPxOppositePion;
  std::vector<TH2D> PytauRecoVSPyOppositePion;
  std::vector<TH2D> PztauRecoVSPzOppositePion;
  std::vector<TH2D> EtauRecoVSEOppositePion;

  std::vector<TH2D> HxResolRecoVSDiffPxtauRecoPxOppositePion;
  std::vector<TH2D> HxResolRecoVSDiffPytauRecoPyOppositePion;
  std::vector<TH2D> HxResolRecoVSDiffPztauRecoPzOppositePion;
  std::vector<TH2D> HxResolRecoVSDiffEtauRecoEOppositePion;
  
  std::vector<TH2D> HyResolRecoVSDiffPxtauRecoPxOppositePion;
  std::vector<TH2D> HyResolRecoVSDiffPytauRecoPyOppositePion;
  std::vector<TH2D> HyResolRecoVSDiffPztauRecoPzOppositePion;
  std::vector<TH2D> HyResolRecoVSDiffEtauRecoEOppositePion;
  
  std::vector<TH2D> HzResolRecoVSDiffPxtauRecoPxOppositePion;
  std::vector<TH2D> HzResolRecoVSDiffPytauRecoPyOppositePion;
  std::vector<TH2D> HzResolRecoVSDiffPztauRecoPzOppositePion;
  std::vector<TH2D> HzResolRecoVSDiffEtauRecoEOppositePion;


  std::vector<TH2D> HxResolSVFitVSTauPt;
  std::vector<TH2D> HxResolSVFitVSPxTau;
  std::vector<TH2D> HxResolSVFitVSPyTau;
  std::vector<TH2D> HxResolSVFitVSPzTau;
  std::vector<TH2D> HxResolSVFitVSETau;
  std::vector<TH2D> HxResolSVFitVSTauEta;
  std::vector<TH2D> HxResolSVFitVSGJAngle;
  std::vector<TH2D> HxResolSVFitVSMH;
  std::vector<TH2D> HxResolSVFitVSSumPtPions;
  std::vector<TH2D> HxResolSVFitVSSumPxPions;
  std::vector<TH2D> HxResolSVFitVSSumPyPions;
  std::vector<TH2D> HxResolSVFitVSSumPzPions;
  std::vector<TH2D> HxResolSVFitVSSumEPions;

  std::vector<TH2D> HxResolRecoVSTauPt;
  std::vector<TH2D> HxResolRecoVSPxTau;
  std::vector<TH2D> HxResolRecoVSPyTau;
  std::vector<TH2D> HxResolRecoVSPzTau;
  std::vector<TH2D> HxResolRecoVSETau;
  std::vector<TH2D> HxResolRecoVSTauEta;
  std::vector<TH2D> HxResolRecoVSGJAngle;
  std::vector<TH2D> HxResolRecoVSMH;
  std::vector<TH2D> HxResolRecoVSSumPtPions;
  std::vector<TH2D> HxResolRecoVSSumPxPions;
  std::vector<TH2D> HxResolRecoVSSumPyPions;
  std::vector<TH2D> HxResolRecoVSSumPzPions;
  std::vector<TH2D> HxResolRecoVSSumEPions;
  
  std::vector<TH2D> HyResolSVFitVSTauPt;
  std::vector<TH2D> HyResolSVFitVSPxTau;
  std::vector<TH2D> HyResolSVFitVSPyTau;
  std::vector<TH2D> HyResolSVFitVSPzTau;
  std::vector<TH2D> HyResolSVFitVSETau;
  std::vector<TH2D> HyResolSVFitVSTauEta;
  std::vector<TH2D> HyResolSVFitVSGJAngle;
  std::vector<TH2D> HyResolSVFitVSMH;
  std::vector<TH2D> HyResolSVFitVSSumPtPions;
  std::vector<TH2D> HyResolSVFitVSSumPxPions;
  std::vector<TH2D> HyResolSVFitVSSumPyPions;
  std::vector<TH2D> HyResolSVFitVSSumPzPions;
  std::vector<TH2D> HyResolSVFitVSSumEPions;

  std::vector<TH2D> HyResolRecoVSTauPt;
  std::vector<TH2D> HyResolRecoVSPxTau;
  std::vector<TH2D> HyResolRecoVSPyTau;
  std::vector<TH2D> HyResolRecoVSPzTau;
  std::vector<TH2D> HyResolRecoVSETau;
  std::vector<TH2D> HyResolRecoVSTauEta;
  std::vector<TH2D> HyResolRecoVSGJAngle;
  std::vector<TH2D> HyResolRecoVSMH;
  std::vector<TH2D> HyResolRecoVSSumPtPions;
  std::vector<TH2D> HyResolRecoVSSumPxPions;
  std::vector<TH2D> HyResolRecoVSSumPyPions;
  std::vector<TH2D> HyResolRecoVSSumPzPions;
  std::vector<TH2D> HyResolRecoVSSumEPions;
  
  std::vector<TH2D> HzResolSVFitVSTauPt;
  std::vector<TH2D> HzResolSVFitVSPxTau;
  std::vector<TH2D> HzResolSVFitVSPyTau;
  std::vector<TH2D> HzResolSVFitVSPzTau;
  std::vector<TH2D> HzResolSVFitVSETau;
  std::vector<TH2D> HzResolSVFitVSTauEta;
  std::vector<TH2D> HzResolSVFitVSGJAngle;
  std::vector<TH2D> HzResolSVFitVSMH;
  std::vector<TH2D> HzResolSVFitVSSumPtPions;
  std::vector<TH2D> HzResolSVFitVSSumPxPions;
  std::vector<TH2D> HzResolSVFitVSSumPyPions;
  std::vector<TH2D> HzResolSVFitVSSumPzPions;
  std::vector<TH2D> HzResolSVFitVSSumEPions;

  std::vector<TH2D> HzResolRecoVSTauPt;
  std::vector<TH2D> HzResolRecoVSPxTau;
  std::vector<TH2D> HzResolRecoVSPyTau;
  std::vector<TH2D> HzResolRecoVSPzTau;
  std::vector<TH2D> HzResolRecoVSETau;
  std::vector<TH2D> HzResolRecoVSTauEta;
  std::vector<TH2D> HzResolRecoVSGJAngle;
  std::vector<TH2D> HzResolRecoVSMH;
  std::vector<TH2D> HzResolRecoVSSumPtPions;
  std::vector<TH2D> HzResolRecoVSSumPxPions;
  std::vector<TH2D> HzResolRecoVSSumPyPions;
  std::vector<TH2D> HzResolRecoVSSumPzPions;
  std::vector<TH2D> HzResolRecoVSSumEPions;
  
};
#endif
