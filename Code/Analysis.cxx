#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <time.h>
#include <cstdlib>
#include <stdlib.h> 

#include "SimpleFits/FitSoftware/interface/Logger.h"
#include <string.h>
#include <istream>
#include <strstream>

#include "TApplication.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TStyle.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"

#include "HistoConfig.h"
#include "SkimConfig.h"
#include "Ntuple_Controller.h"
#include "Selection_Base.h"
#include "Selection_Factory.h"
#include "DoubleEventRemoval.h"
#include "Parameters.h"
#include "Plots.h"

int main() {
  Logger::Instance()->SetLevel(Logger::Info);
  Logger(Logger::Info) << "\n\n\n================" << std::endl;
  Logger(Logger::Info) << "Program is Starting" << endl;
  system("date");
  float neventsproc(0);
  time_t startTime, endTime, beforeLoop, afterLoop;
  time(&startTime);

  //TApplication App(argv[0],&argc,argv);
  gROOT->ProcessLine("#include <vector>");
  gSystem->Load("$ROOTSYS/lib/libPhysics.so");
  gROOT->SetBatch(kTRUE);

  ////////////////////////////////////////////////
  // Get Input
  TString base = std::getenv("PWD");
  string basetemp=std::getenv("PWD");
  base += "/";
  Logger(Logger::Info) << "Working Dir: " << base << std::endl;
  Parameters Par(base + "Input.txt");
  std::vector<TString> Files, UncertType, UncertList, Analysis;
  std::vector<double> UncertW;
  bool thin, skim;
  int mode, runtype;
  TString mode_str, runType_str, histofile, skimfile, PlotStyle, PlotLabel;
  double Lumi;
  Par.GetVectorString("File:", Files);
  Par.GetBool("Thin:", thin, "False");
  Par.GetBool("Skim:", skim, "False");
  Par.GetString("Mode:", mode_str, "RECONSTRUCT");  // RECONSTRUCT/ANALYSIS
  Par.GetString("RunType:", runType_str, "LOCAL");  // GRID/LOCAL
  Par.GetString("SkimInfo:", skimfile, "");
  Par.GetString("HistoFile:", histofile, "InputData/Histo.txt");
  Par.GetVectorString("Analysis:", Analysis);
  Par.GetVectorString("UncertType:", UncertType, "<default>");
  Par.GetVectorStringDouble("UncertList:", UncertList, UncertW);
  Par.GetDouble("Lumi:", Lumi, 1);
  Par.GetString("PlotStyle:", PlotStyle, "style1");
  Par.GetString("PlotLabel:", PlotLabel, "none");
  /////////////////////////////////////////////////
  // Check Input
  HistoConfig H;
  SkimConfig SC;
  if (!H.Load(base + histofile)) {
    Logger(Logger::Fatal) << "No Histogram File or Histograms!!!" << endl;
    exit(6);
  }
  if (!SC.Load(base + skimfile)) {
    Logger(Logger::Fatal) << "Invalid Skim File!!!" << endl;
    exit(6);
  }
  if (Analysis.size() == 0) {
    Logger(Logger::Fatal) << "No Analysis Selected!!!!" << endl;
    exit(6);
  }
  if (UncertType.size() == 0) {
    Logger(Logger::Fatal) << "Uncertainty Configuration!!!!" << endl;
    exit(6);
  }
  if (Files.size() == 0) {
    Logger(Logger::Fatal) << "No Input Files!!!!" << endl;
    exit(6);
  }
  mode_str.ToUpper();
  if (mode_str == "ANALYSIS") {
    Logger(Logger::Info) << "Using Mode: ANALYSIS" << std::endl;
    mode = Selection_Base::ANALYSIS;
  } else if (mode_str == "RECONSTRUCT") {
    Logger(Logger::Info) << "Using Mode: RECONSTRUCT" << std::endl;
    mode = Selection_Base::RECONSTRUCT;
  } else {
    Logger(Logger::Fatal) << "No Mode Type!!!!" << endl;
    exit(6);
  }
  runType_str.ToUpper();
  if (runType_str == "GRID") {
    Logger(Logger::Info) << "Using RunType: GRID" << std::endl;
    runtype = Selection_Base::GRID;
  } else if (runType_str == "LOCAL") {
    Logger(Logger::Info) << "Using RunType: Local" << std::endl;
    runtype = Selection_Base::Local;
  } else {
    Logger(Logger::Fatal) << "No RunType!!!!" << endl;
    exit(6);
  }
  Plots P;
  P.Set_Plot_Type(PlotStyle, PlotLabel);
  //////////////////////////////////////////////////
  // Configure Analysis
  Selection_Factory SF;
  vector<Selection_Base*> selections;
  for (unsigned int i = 0; i < UncertType.size(); i++) {
    Logger(Logger::Info) << "Configuring systematic " << UncertType.at(i) << endl;
    for (unsigned int j = 0; j < Analysis.size(); j++) {
      Logger(Logger::Info) << "Configuring Selection " << Analysis.at(j) << endl;
      selections.push_back(SF.Factory(Analysis.at(j), UncertType.at(i), mode, runtype, Lumi));
    }
  }
  Logger(Logger::Info) << "Selection modules Setup NSelection= " << selections.size() << endl;

  //////////////////////////////////////////////////
  // Run Analysis on Ntuples
  string DPM_PATH="";
  string SkimDirectory="";
  string InputName="";
  string OutputPath="";
  if (mode == Selection_Base::ANALYSIS) {
    Logger(Logger::Info) << "Setuping up Ntuple control" << endl;
    Ntuple_Controller Ntp(Files);
    if(skim) {
      DPM_PATH = "/dpm/in2p3.fr/home/cms/phedex/store/user/gbourgat/";
      SkimDirectory="Skim2016_180220/";
      InputName=Files.at(0);
      stringstream ss(InputName);
      string SubString;
      vector<string> elements;
      while (std::getline(ss, SubString, '/'))
	{
	  elements.push_back(SubString);
	}

      system(("rfmkdir "+DPM_PATH + SkimDirectory + elements.at(12)).c_str());
      system(("rfmkdir "+DPM_PATH + SkimDirectory + elements.at(12) + "/" + elements.at(13)).c_str());
      system(("rfmkdir "+DPM_PATH + SkimDirectory + elements.at(12) + "/" + elements.at(13) + "/" + elements.at(14)).c_str());
      
      std::vector<string>::iterator it;

      it = std::find (elements.begin(), elements.end(), "alebihan");
      if (it != elements.end()){
	system(("rfmkdir "+DPM_PATH + SkimDirectory + elements.at(12) + "/" + elements.at(13) + "/" + elements.at(14)+ "/" + elements.at(15) ).c_str());
	system(("rfmkdir "+DPM_PATH + SkimDirectory + elements.at(12) + "/" + elements.at(13) + "/" + elements.at(14)+ "/" + elements.at(15) + "/0000").c_str());
	OutputPath="root://sbgse1.in2p3.fr/" + DPM_PATH + SkimDirectory + elements.at(12) + "/" + elements.at(13) + "/" + elements.at(14)+ "/" + elements.at(15)+"/0000/";
      }
      else{
	system(("rfmkdir "+DPM_PATH + SkimDirectory + elements.at(12) + "/" + elements.at(13) + "/" + elements.at(14)+ "/0000").c_str());
	OutputPath="root://sbgse1.in2p3.fr/" + DPM_PATH + SkimDirectory + elements.at(12) + "/" + elements.at(13) + "/" + elements.at(14)+ "/0000/";
      }
    }
    if (thin)
      Ntp.ThinTree();
    size_t last_index=0;
    string Set="";
    if (skim){
      last_index = basetemp.find_last_not_of("0123456789");
      Set = basetemp.substr(last_index + 1);
      Ntp.CloneTree(OutputPath+"SKIMMED_NTUP_"+Set);}
    for (unsigned int j = 0; j < selections.size(); j++) {
      selections[j]->Set_Ntuple(&Ntp);
    }

    //////////////////////////////////////
    // Event Loop
    Logger(Logger::Info) << "Starting Event Loop" << endl;
    Int_t nentries = Ntp.Get_Entries();
    std::vector<TString> ListOfFilesRead;
    std::vector<int> EventsReadFromFile;
    int RunNumber(-999);
    int num = 1000;
    Logger(Logger::Info) << "Will run over " << nentries << " events" << endl;
    int i = 0;
    unsigned int k = 0;
    int p = 0;
    DoubleEventRemoval DER;
    time(&beforeLoop);
    neventsproc = (float) nentries;
    for (i = 0; i < nentries; i++) {
      TString fileName = Ntp.Get_File_Name();
      if (i == 0 || fileName != ListOfFilesRead[k]) {
	ListOfFilesRead.push_back(fileName);
	EventsReadFromFile.push_back(0);
	k = ListOfFilesRead.size() - 1;
      }
      EventsReadFromFile[k]++;
      num++;
      if (num >= 1000) {
	Logger(Logger::Info) << "Starting event:" << i << " out of " << nentries << endl;
	num = 0;
      }
      Ntp.Get_Event(i);
      if (Ntp.isData()) {
	if (!DER.CheckDoubleEvents(Ntp.RunNumber(), Ntp.EventNumber()))
	  continue;
	if (RunNumber != Ntp.RunNumber()) {
	  RunNumber = Ntp.RunNumber();
	  Logger(Logger::Debug) << "RunNumber: " << RunNumber << endl;
	}
      }
      bool passed = false;

      for (unsigned int j = 0; j < selections.size(); j++) {
	selections.at(j)->Event();
	if (selections.at(j)->Passed())
	  passed = true;
      }
      if (skim && passed) {

	Ntp.AddEventToCloneTree();

	p++;
      }
    }

    time(&afterLoop);
    if (skim){

      Ntp.SaveCloneTree();

    }
    if (ListOfFilesRead.at(0) == ListOfFilesRead.at(ListOfFilesRead.size() - 1))
      ListOfFilesRead.pop_back();
    if (i == nentries && Files.size() == ListOfFilesRead.size()) {
      Logger(Logger::Info) << "Number of events read: " << i << " out of " << nentries << " SUCCESSFULL" << endl;
    } else {
      Logger(Logger::Info) << "Number of events read: " << i << " out of " << nentries << " ; Number of files read: " << ListOfFilesRead.size() << " out of " << Files.size() << " FAILED"
			   << endl;
    }
    Logger(Logger::Info) << "Number of skimmed events in SkimmedD3PD: " << p << endl;
    Logger(Logger::Info) << "Number of files read: " << ListOfFilesRead.size() << endl;
    for (k = 0; k < ListOfFilesRead.size(); k++) {
      Logger(Logger::Verbose) << "Number of events: " << EventsReadFromFile[k] << " from file: " << ListOfFilesRead[k] << endl;
    }
    Logger(Logger::Info) << "Event Loop done" << endl;
  }
  ///////////////////////////////////////////
  // Reconstruct Analysis from Stored Histograms
  else if (mode == Selection_Base::RECONSTRUCT) {
    Logger(Logger::Info) << "Reconstructing histograms: Loading files" << endl;
    for (unsigned int j = 0; j < selections.size(); j++) {
      selections[j]->LoadResults(Files);
    }
    Logger(Logger::Info) << "Loading Files Complete" << std::endl;
    if (runtype == Selection_Base::Local) {
      Selection_Factory SF;
      for (unsigned int j = 0; j < selections.size(); j++) {
	if (selections.at(j)->Get_SysType() == "default") {
	  for (unsigned int i = 0; i < UncertList.size(); i++) {
	    TString n = selections.at(j)->Get_Name();
	    Logger(Logger::Info) << "Adding Systematic Uncertainty " << UncertType.at(i) << endl;
	    Selection_Base *s = SF.Factory(selections.at(j)->Get_Analysis(), UncertList.at(i), mode, runtype, Lumi);
	    s->LoadResults(Files);
	    selections[j]->EvaluateSystematics(s, 1.0);
	    delete s;
	  }
	}
      }
    }
  }
  Logger(Logger::Verbose) << "Finishing" << endl;
  for (unsigned int j = 0; j < selections.size(); j++) {
    Logger(Logger::Info) << "Finishing: " << selections[j]->Get_Name() << endl;
    selections[j]->Finish();
  }
  Logger(Logger::Verbose) << "Finished" << endl;

  Logger(Logger::Info) << "Cleaning up" << endl;
  for (unsigned int j = 0; j < selections.size(); j++) {
    delete selections[j];
    selections[j] = NULL;
  }
  selections.clear();

  if (mode != Selection_Base::RECONSTRUCT) {
    double dt = difftime(afterLoop, beforeLoop);
    int hour = (int) (dt / (60 * 60));
    int min = (int) (dt / (60)) - 60 * hour;
    int sec = (int) (dt - min * 60 - hour * 60 * 60);
    Logger(Logger::Info) << "Time in Ntuple Loop: " << hour << "hr " << min << "min " << sec << "s" << endl;
    float msec = 1000 * dt / neventsproc;
    Logger(Logger::Info) << "Time in per event in Ntuple Loop: " << msec << "ms" << endl;
  }

  time(&endTime);
  double dt = difftime(endTime, startTime);
  int hour = (int) (dt / (60 * 60));
  int min = (int) (dt / (60)) - 60 * hour;
  int sec = (int) (dt - min * 60 - hour * 60 * 60);

  Logger(Logger::Info) << "Run Time: " << hour << "hr " << min << "min " << sec << "s" << endl;
  Logger(Logger::Info) << "Program is Finished" << endl;
  system("date");
}

