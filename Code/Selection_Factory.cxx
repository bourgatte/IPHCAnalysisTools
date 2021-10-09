#include "Selection_Factory.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"

#include "Example.h"
#include "TauSpinExample.h"
#ifdef USE_cherepanov
#include "cherepanov/MyTest.h"
#include "cherepanov/NtupleValidation.h"
#include "cherepanov/SkimmingNtuples.h"
#include "cherepanov/ControlSample.h"
#include "cherepanov/SkimNtupleDiTauHTrigger.h"
#include "cherepanov/ZTauHTauH.h"
#include "cherepanov/ZTauTau.h"
#include "cherepanov/HTauTau.h"
#include "cherepanov/ZMuTau.h"
#include "cherepanov/SingleMuSkim.h"
#include "cherepanov/ZTauMuTauH.h"
#include "cherepanov/TTBar.h"
#include "cherepanov/TemplateStudyMT.h"
#endif

#ifdef USE_goe

#endif

#ifdef USE_lebihan

#endif

#ifdef USE_gbourgat

//#include "gbourgat/HTauTau.h"
#include "gbourgat/HCPTauTau.h"
#include "gbourgat/HCPMuTau.h"
#include "gbourgat/HCPPiTau.h"

#endif

#ifdef USE_cgrimault

#endif

#ifdef USE_msessini

#include "msessini/HCPTauTau.h"

#endif

Selection_Factory::Selection_Factory(){
}

Selection_Factory::~Selection_Factory(){
}

Selection_Base* Selection_Factory::Factory(TString Analysis, TString UncertType,int mode, int runtype, double lumi){
  Selection_Base* s;
  Analysis.ToLower();

  // ensuring code will compile independently of user code
  // WARNING: be aware of the consequences of "Contains". Make sure that Class "foo" is put after "foobar".
  if(Analysis.Contains("example"))s=new Example(Analysis,UncertType);
  else if(Analysis.Contains("tauspin"))s=new TauSpinExample(Analysis,UncertType);
#ifdef USE_cherepanov
  else if(Analysis.Contains("mytest"))s=new MyTest(Analysis,UncertType);
  else if(Analysis.Contains("ztauhtauh"))s=new ZTauHTauH(Analysis,UncertType);
  else if(Analysis.Contains("ztautau"))s=new ZTauTau(Analysis,UncertType);
  else if(Analysis.Contains("htautau"))s=new HTauTau(Analysis,UncertType);
  else if(Analysis.Contains("zmutau"))s=new ZMuTau(Analysis,UncertType);
  else if(Analysis.Contains("singlemuskim"))s=new SingleMuSkim(Analysis,UncertType);
  else if(Analysis.Contains("ztaumutauh"))s=new ZTauMuTauH(Analysis,UncertType);
  else if(Analysis.Contains("ttbar"))s=new TTBar(Analysis,UncertType);
  else if(Analysis.Contains("templatestudymt"))s=new TemplateStudyMT(Analysis,UncertType);



#endif
// #ifdef USE_goe
//   else if(Analysis.Contains("bla"))s=new Bla(Analysis,UncertType);
// #endif

// #ifdef USE_lebihan
//   else if(Analysis.Contains("bla"))s=new Bla(Analysis,UncertType);
// #endif

// #ifdef USE_cgrimault
//   else if(Analysis.Contains("bla"))s=new Bla(Analysis,UncertType);
// #endif



#ifdef USE_gbourgat

  //else if(Analysis.Contains("htautau"))s=new HTauTau(Analysis,UncertType);
  else if(Analysis.Contains("hcptautau"))s=new HCPTauTau(Analysis,UncertType);
  else if(Analysis.Contains("hcpmutau"))s=new HCPMuTau(Analysis,UncertType);
  else if(Analysis.Contains("hcppitau"))s=new HCPPiTau(Analysis,UncertType);

#endif

#ifdef USE_msessini

  else if(Analysis.Contains("hcptautau"))s=new HCPTauTau(Analysis,UncertType);

#endif

  else{
	Logger(Logger::Error)<< "Invalid Analysis type \"" << Analysis << "\". Using default <Example.h> " << std::endl;
    s=new Example(Analysis,UncertType);
  }
  s->SetMode(mode);
  s->SetRunType(runtype);
  s->SetLumi(lumi);
  s->Configure();
  return s;
}
