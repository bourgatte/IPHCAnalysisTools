#include "HistoConfig.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"

#include <cstdlib>
#include <algorithm>                                     
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include "Riostream.h"
#include <fstream>
#include <istream>
#include <strstream>
#include <cstdlib>
#include "TH1D.h"
#include <math.h>

// Static var
std::vector<Long64_t>      HistoConfig::ID;
std::vector<double>       HistoConfig::CS;
std::vector<TString>      HistoConfig::HistoName;
std::vector<TString>      HistoConfig::HistoLegend;
std::vector<int>          HistoConfig::HistoColour;
bool                      HistoConfig::loaded=false;

using namespace std;

HistoConfig::HistoConfig(){
}


bool HistoConfig::Load(){
  return HistoConfig::Load("InputData/Par.dat");
}

bool HistoConfig::Load(TString Name_)  
{
  if(loaded) return true;
  Logger(Logger::Verbose) << "HistoConfig::Load("<< Name_ <<")" << std::endl;
  ID.clear();
  HistoName.clear();
  HistoLegend.clear();

  // Open File
  ifstream input_file;
  char *file_=(char*)Name_.Data();
  input_file.open(file_, std::ios::in);
  if (!(input_file)){
	  Logger(Logger::Error) << "Opening xml file "<< Name_ <<" for HistoConfig has failed." << std::endl;
    return false;
  }
  loaded=true;
  Logger(Logger::Info) << "Opened HistoConfig xml file: "<< Name_ << std::endl;

  std::string s;
  unsigned int a=0;
  while(getline(input_file, s)){
    a++;
    if(a>250) break;
    std::stringstream line(s); 
    TString type;
    Long64_t id;
    double cs;
    TString name;
    TString leg;
    int colour;
    line >> type >> id >> cs >> name >> leg >> colour;
    type.ToLower();
    if(!type.Contains("histo:")) continue;
    bool isnew=true;
    for(unsigned int i=0; i<ID.size();i++){
      if(ID.at(i)==id) isnew=false;
    }
    if(isnew){
      ID.push_back(id);
      CS.push_back(cs);
      HistoName.push_back(name);
      HistoLegend.push_back(leg);
      HistoColour.push_back(colour);
    }
  }
  input_file.close();
  for(unsigned int i=0; i<ID.size();i++){
	  Logger(Logger::Verbose) << "Histogram Data/MC ID: " << ID.at(i) << " CS: " << CS.at(i) << " Name: " <<  HistoName.at(i) << " Legend: " <<  HistoLegend.at(i) << " Colour: " << HistoColour.at(i) << std::endl;
  }

  Logger(Logger::Verbose) << "HistoConfig::Load("<< Name_ <<") complete" << std::endl;
  if(HistoName.size()>=1) return true;
  return false;
}

HistoConfig::~HistoConfig(){

}

//utility Functions
bool HistoConfig::GetHisto(bool isdata,Long64_t id,unsigned int &histo){
  if(isdata && id!=36){
    id=1;
  }
  for(unsigned int i=0; i<ID.size(); i++){
    if(ID.at(i)==id){
      histo=i;
      return true;
    } 
  }
  return false;
}

double HistoConfig::GetCrossSection(Long64_t id){
  for(unsigned int i=0; i<ID.size(); i++){
    if(ID.at(i)==id){
      return CS.at(i);
    }
  }
  return 0;
}

bool HistoConfig::SetCrossSection(Long64_t id, double xsec){
	for(unsigned int i=0; i<ID.size(); i++){
		if(ID.at(i)==id){
			CS.at(i) = xsec;
			return true;
	    }
	}
	return false;
}

void HistoConfig::GetHistoInfo(std::vector<Long64_t> &types,std::vector<float> &CrossSectionandAcceptance,std::vector<TString> &legend,std::vector<int> &colour){
  types=ID;
  legend=HistoLegend;
  colour=HistoColour;
  CrossSectionandAcceptance.clear();
  for(unsigned int i=0; i<HistoName.size();i++){
    CrossSectionandAcceptance.push_back(CS.at(i));
  }
}

unsigned int HistoConfig::GetNHisto(){
  return HistoName.size();
}

TString HistoConfig::GetName(unsigned int i){
  if(i>=0 && i<HistoName.size()) return HistoName.at(i);
  return "unknown";
}

TString HistoConfig::GetLeg(unsigned int i){
  if(i>=0 && i<HistoLegend.size()) return HistoLegend.at(i);
  return "unknown";
}

std::vector<TH1D> HistoConfig::GetTH1D(TString name,TString title, int nbins, double min, double max, TString xaxis, TString yaxis){
  std::vector<TH1D> histos;
  Logger(Logger::Verbose) << "Adding TH1D " << name << " " << title << std::endl;
  for(unsigned int i=0;i<HistoName.size();i++){
    histos.push_back(TH1D(name+HistoName.at(i),HistoLegend.at(i),nbins,min,max));
    histos.at(i).Sumw2();
    histos.at(i).SetXTitle(xaxis);
    histos.at(i).SetYTitle(yaxis);
  }
  return histos;
}

std::vector<TH1D> HistoConfig::GetTH1D(TString name,TString title, int nbins, double* xbins, TString xaxis,TString yaxis){
  std::vector<TH1D> histos;
  Logger(Logger::Verbose) << "Adding TH1D " << name << " " << title << std::endl;
  for(unsigned int i=0;i<HistoName.size();i++){
    histos.push_back(TH1D(name+HistoName.at(i),HistoLegend.at(i),nbins,xbins));
    histos.at(i).Sumw2();
    histos.at(i).SetXTitle(xaxis);
    histos.at(i).SetYTitle(yaxis);
  }
  return histos;
}

std::vector<TH2D> HistoConfig::GetTH2D(TString name,TString title,int nbinsx, double minx, double maxx, 
				       int nbinsy, double miny, double maxy, TString xaxis, TString yaxis){
  std::vector<TH2D> histos;
  Logger(Logger::Verbose) << "Adding TH2D " << name << " " << title << std::endl;
  for(unsigned int i=0;i<HistoName.size();i++){
    histos.push_back(TH2D(name+HistoName.at(i),HistoLegend.at(i),nbinsx,minx,maxx, nbinsy,miny,maxy));
    histos.at(i).Sumw2();
    histos.at(i).SetXTitle(xaxis);
    histos.at(i).SetYTitle(yaxis);
  }
  return histos;
}

std::vector<TH3F> HistoConfig::GetTH3F(TString name,TString title, int nbinsx, double minx, double maxx,
			  int nbinsy,double miny, double maxy,int nbinsz,double minz,double maxz,
			  TString xaxis, TString yaxis,TString zaxis){

  std::vector<TH3F> histos;
  Logger(Logger::Verbose) << "Adding TH2D " << name << " " << title << std::endl;
  for(unsigned int i=0;i<HistoName.size();i++){
    histos.push_back(TH3F(name+HistoName.at(i),HistoLegend.at(i),nbinsx,minx,maxx,nbinsy,miny,maxy,nbinsz,minz,maxz));
    histos.at(i).Sumw2();
    histos.at(i).SetXTitle(xaxis);
    histos.at(i).SetYTitle(yaxis);
    histos.at(i).SetZTitle(zaxis);
  }
  return histos;
}



bool HistoConfig::hasID(Long64_t id_){
  for(unsigned int i=0; i<ID.size();i++){
    if(ID.at(i)==id_) return true;
  }
  return false;
}


int HistoConfig::GetID(unsigned int i){
  if(ID.size()>i) return ID.at(i);
  return -999;
}

int HistoConfig::GetType(Long64_t id){
	for(unsigned int i=0; i<ID.size(); i++){
	    if(ID.at(i)==id){
	    	return i;
	    }
	}
	return -1;
}
