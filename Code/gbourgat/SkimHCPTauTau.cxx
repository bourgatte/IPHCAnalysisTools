#include "SkimHCPTauTau.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "SVFitObject.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"
 
SkimHCPTauTau::SkimHCPTauTau(TString Name_, TString id_):
  Selection(Name_,id_)
{
}

SkimHCPTauTau::~SkimHCPTauTau(){
  for(unsigned int j=0; j<Npassed.size(); j++){
	 Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "complete." << std::endl;
}

void  SkimHCPTauTau::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)      cut.at(TriggerOk)=1;
    if(i==PrimeVtx)       cut.at(PrimeVtx)=1;
    if(i==ntaus)          cut.at(ntaus)=1;

  }
  // Setup cut plots
  TString hlabel;
  TString htitle;
  for(int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i;
    if(i==PrimeVtx){
      title.at(i)="Number of Prime Vertices $(N>$";
      title.at(i)+=cut.at(PrimeVtx);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of Prime Vertices";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PrimeVtx_",htitle,31,-0.5,30.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PrimeVtx_",htitle,31,-0.5,30.5,hlabel,"Events"));
    }
    else if(i==TriggerOk){
      title.at(i)="Trigger ";
      hlabel="Trigger ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    else if(i==ntaus){
      title.at(i)="Number of tau leptons ";
      hlabel="Number of tau leptons  ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ntaus_",htitle,10,-0.5,9.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ntaus_",htitle,10,-0.5,9.5,hlabel,"Events"));
    }
  } 
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  TauDecayMode=HConfig.GetTH1D(Name+"_TauDecayMode","HPS decay type",11,-0.5,10.5," pType","Events");

  TauPT=HConfig.GetTH1D(Name+"_TauPT","Transverse momentum #tau",50,15,75," P_{T}, GeV","Events");
  JetPT=HConfig.GetTH1D(Name+"_JetPT","Transverse momentum jet",50,15,125," P_{T}, GeV","Events");
  MissingTEnergy=HConfig.GetTH1D(Name+"_MissingTEnergy","Missing Transverse Energy",50,-0.05,99.5," M_{T}, GeV","Events");
  OSPairMass=HConfig.GetTH1D(Name+"_OSPairMass"," Mass of two OS candidates",50,50,100," M_{pair}, GeV","Events");
  SSPairMass=HConfig.GetTH1D(Name+"_SSPairMass"," Mass of two SS candidates",50,50,100," M_{pair}, GeV","Events");



  isPairCandOS=HConfig.GetTH1D(Name+"_isPairCandOS","is pair OS ",2,-0.5,1.5,"OS/SS","Events");
    Selection::ConfigureHistograms();   //   do not remove
    HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);  // do not remove
} 

 

void  SkimHCPTauTau::Store_ExtraDist(){

  //every new histo should be addedd to Extradist1d vector, just push it back;
  Extradist1d.push_back(&TauDecayMode);
  Extradist1d.push_back(&MissingTEnergy);
  Extradist1d.push_back(&OSPairMass);
  Extradist1d.push_back(&SSPairMass);
  Extradist1d.push_back(&isPairCandOS);
  Extradist1d.push_back(&TauPT);
  Extradist1d.push_back(&JetPT);

}

void  SkimHCPTauTau::doEvent(){ //  Method called on every event
  

  unsigned int t;                // sample type, you may manage in your further analysis, if needed
  int id(Ntp->GetMCID());  //read event ID of a sample
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ Logger(Logger::Error) << "failed to find id" <<std::endl; return;}  //  gives a warining if list of samples in Histo.txt  and SkimSummary.log do not coincide 
  //  std::cout<<"------------------ New Event -----------------------"<<std::endl;
  // bool PassedTrigger(false);
  // int triggerindex;
  // std::vector<int> TriggerIndex; 
  // std::vector<int>TriggerIndexVector ;
  // std::vector<TString>  MatchedTriggerNames;

  //MatchedTriggerNames.push_back("HLT_DoubleMediumIsoPFTau");
  //MatchedTriggerNames.push_back("HLT_DoubleMediumCombinedIsoPFTau35");
  //TriggerIndexVector=Ntp->GetVectorTriggers(MatchedTriggerNames);
   
 
   // for(unsigned int itr=0; itr<Ntp->NTriggers(); itr++){

   //   if(Ntp->TriggerAccept(itr)) std::cout<<"   itr   "<< itr <<"   trig name  "<< Ntp->TriggerName(itr) << "  accept?    "<<Ntp->TriggerAccept(itr) << std::endl;
   // }
  //int triggerMask(0);
  //triggerMask |= (1<< Ntp->getBitOfGivenTrigger("HLT_IsoMu22_eta2p1_v"));
  // PassedTrigger = (( Ntp->triggerbit() & triggerMask ) == triggerMask);
  
  // for(unsigned int itrig = 0; itrig < TriggerIndexVector.size(); itrig++){
  //   if(Ntp->TriggerAccept(TriggerIndexVector.at(itrig))){  
  //     //      std::cout<<"  Name  "<< Ntp->TriggerName(TriggerIndexVector.at(itrig)) << "   status   "<< Ntp->TriggerAccept(TriggerIndexVector.at(itrig)) <<std::endl;
  //     PassedTrigger =Ntp->TriggerAccept(TriggerIndexVector.at(itrig)); }
  // }


  std::vector<int> Tau1IndexVect;
  std::vector<int> Tau2IndexVect;
  std::vector<int> IndexSelected;
  for(int iDaughter=0;   iDaughter  <  Ntp->SelectedPairs() ;iDaughter++ ) {
    if(Ntp->trg_doubletau(iDaughter))
      {
  	Tau1IndexVect.push_back(Ntp->tau1IndexVect(iDaughter));
	Tau2IndexVect.push_back(Ntp->tau2IndexVect(iDaughter));
	IndexSelected.push_back(iDaughter);
      }
  }
  value.at(TriggerOk)=0;
  if(Tau1IndexVect.size()>0 && Tau2IndexVect.size()>0)value.at(TriggerOk)=1;
  pass.at(TriggerOk)=(value.at(TriggerOk)==cut.at(TriggerOk));


  std::vector<int> GoodTausIndex;

  for(unsigned int iDaughter=0;   iDaughter  <  IndexSelected.size() ;iDaughter++ ){  // loop over all daughters in the event
      if(Ntp->tauBaselineSelection(Tau1IndexVect.at(iDaughter),38., 2.1, 0,0,0) && Ntp->tauBaselineSelection(Tau2IndexVect.at(iDaughter),38., 2.1, 0,0,0))
	{
	  GoodTausIndex.push_back(iDaughter);
	}
  }
  


  // for(int itrig = 0; itrig < TriggerIndex.size(); itrig++){
  //   if(Ntp->GetTriggerIndex((TString)"HLT_IsoMu", TriggerIndex.at(itrig))) {
  //     if(Ntp->TriggerAccept(TriggerIndex.at(itrig))){   PassedTrigger =Ntp->TriggerAccept(TriggerIndex.at(itrig)); }
  //   }
  // }

 

  // Apply Selection
  // two vectors value and pass  are used to apply selection, store value.at(<A cut variable defined in .h as enumerator>) a quantitity you want to use in your selection.
  // Define in the pass.at(<A cut variable defined in .h as enumerator>) a passing condidtion, true/false
  value.at(PrimeVtx)=Ntp->NVtx();
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx));
  
      //value.at(TriggerOk)=PassedTrigger;
      //pass.at(TriggerOk)=PassedTrigger;
  
  value.at(ntaus)=GoodTausIndex.size();
  pass.at(ntaus) = (value.at(ntaus) >= cut.at(ntaus));


  // Here you can defined different type of weights you want to apply to events. At the moment only PU weight is considered if event is not data
  double wobs=1;
  double w=1;

  bool status=AnalysisCuts(t,w,wobs);  // boolean that say whether your event passed critera defined in pass vector. The whole vector must be true for status = true
  ///////////////////////////////////////////////////////////
  // Analyse events  which passed selection
  if(status){

    for(unsigned int ipair=0; ipair < Ntp->NPairCandidates(); ipair++){
      isPairCandOS.at(t).Fill(Ntp->isOSCand(ipair),w);
      if(Ntp->isOSCand(ipair)==1)      OSPairMass.at(t).Fill((  (Ntp->Daughters_P4(Ntp->indexDau1(ipair)) +  Ntp->Daughters_P4(Ntp->indexDau2(ipair) ) ).M()  ),w);
      if(Ntp->isOSCand(ipair)==0)      SSPairMass.at(t).Fill((  (Ntp->Daughters_P4(Ntp->indexDau1(ipair)) +  Ntp->Daughters_P4( Ntp->indexDau2(ipair) ) ).M()  ),w);
    }
    

    MissingTEnergy.at(t).Fill(Ntp->MET(),w);
    //if(GoodTausIndex.size()!=0)   TauDecayMode.at(t).Fill(Ntp->decayMode(GoodTausIndex.at(0)),w);

    //if(GoodTausIndex.size()!=0)   TauPT.at(t).Fill(Ntp->Daughters_P4(GoodTausIndex.at(0)).Pt(),w);

    //    if(GoodTausIndex.size()!=0)      if(Ntp->Daughters_charge(GoodTausIndex.at(0))*Ntp->Daughters_charge(GoodMuonIndex.at(0))==-1) LeadPairMass.at(t).Fill(  (Ntp->Daughters_P4(GoodTausIndex.at(0)) + Ntp->Daughters_P4(GoodMuonIndex.at(0))).M(),w   );
    
         
    
    if(Ntp->NJets()!=0)JetPT.at(t).Fill(Ntp->Jet_P4(0).Pt(),w);
  }
}




//  This is a function if you want to do something after the event loop
void  SkimHCPTauTau::Finish(){

  Selection::Finish();

	for (unsigned int i = 0; i < CrossSectionandAcceptance.size(); i++) {
	  std::cout<< "  i  "<< i<< "  HConfig.GetID(i) "<< HConfig.GetID(i) << " Npassed.GetBinContent(0) "<< Npassed.at(i).GetBinContent(0) << " Npassed.GetBinContent(1) "<< Npassed.at(i).GetBinContent(1) << " Npassed.GetBinContent(NCuts) "<< Npassed.at(i).GetBinContent(NCuts) <<  "   no weight 0  "<< Npassed_noweight.at(i).GetBinContent(0)<<  "   no weight 1  "<< Npassed_noweight.at(i).GetBinContent(1)<<  "   no weight Ncuts   "<< Npassed_noweight.at(i).GetBinContent(NCuts)<< " scale   "<< Selection::scaleFactorToLumi(HConfig.GetID(i))  << std::endl;

	  int nbins = Npassed.at(0).GetNbinsX();

	  for(unsigned int ibn=0; ibn < nbins+1; ibn ++){
	    std::cout<<"   "<<  HConfig.GetID(i) << "  Npassed  "<< Npassed.at(i).GetBinContent(ibn)<< "  nowe  " << Npassed_noweight.at(i).GetBinContent(ibn)<<std::endl;

	  }

	}



	


}





