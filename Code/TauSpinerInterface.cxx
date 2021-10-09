#include "TauSpinerInterface.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"

#ifdef USE_TauSpinner
#include "Tauola.h"
#include "LHAPDF/LHAPDF.h"
#include "tau_reweight_lib.h"
#include "read_particles_from_TAUOLA.h"
#include<iostream>
#include "TLorentzVector.h"

int TauSpinerInterface::signalcharge=-1;
bool TauSpinerInterface::initialized=false;
#endif

TauSpinerInterface::TauSpinerInterface(){
}

TauSpinerInterface::~TauSpinerInterface(){

}
#ifdef USE_TauSpinner
void TauSpinerInterface::Initialize(){
  Tauolapp::Tauola::initialize();
  string name="NNPDF30_nlo_as_0118";
  LHAPDF::initPDFSetByName(name);
  double CMSENE = 13000.0; // center of mass system energy.
                          // used in PDF calculation. For pp collisions only
  bool Ipp = true;  // for pp collisions 
  // Initialize TauSpinner
  //Ipol - polarization of input sample
  //nonSM2 - nonstandard model calculations
  //nonSMN
  int Ipol=0,nonSM2=0,nonSMN=0;

  initialize_spinner(Ipp,Ipol,nonSM2,nonSMN,CMSENE);
  setHiggsParametersTR(-1.0, 1.0, 0.0, 0.0); // Scalar Higgs (Even)
  //setHiggsParametersTR(1.0, -1.0, 0.0, 0.0); // Pseudo-Scalar Higgs (Odd)
  //double theta=TMath::Pi()/4.;
  //setHiggsParametersTR(-cos(2*theta),cos(2*theta) ,-sin(2*theta),-sin(2*theta));
}



double TauSpinerInterface::Get(int type, SimpleParticle X, SimpleParticle tau, std::vector<SimpleParticle> tau_daughters,SimpleParticle tau2, std::vector<SimpleParticle> tau_daughters2){
  if(!initialized){
    Initialize();
    initialized=true;
  }


  double WT    = 1.0; // assume that there is 1 bosons decaying into                                                                                                                                                                        
  // Calculate weight for first boson
  if( abs(X.pdgid())==24 ||  abs(X.pdgid())==37 ){  
    WT = calculateWeightFromParticlesWorHpn(X, tau, tau2, tau_daughters); // note that tau2 is tau neutrino
  }
  else if( X.pdgid()==25 || X.pdgid()==36 || X.pdgid()==22 || X.pdgid()==23 ){
    /*
    std::cout << "simpleParticle - 1 " << std::endl;
    std::cout << "simpleParticle " << X.pdgid() << " " << X.px() << " " <<X.py() << " " <<X.pz() << " " << X.e() << std::endl; 
    std::cout << "simpleParticle " <<  tau.pdgid() << " " << tau.px() << " " <<tau.py() << " " <<tau.pz() << " " << tau.e() 
	      << " size " << tau_daughters.size() << std::endl;
    TLorentzVector Tau(0,0,0,0);
    for(unsigned int i=0;i<tau_daughters.size();i++){
      std::cout << "simpleParticle" << tau_daughters.at(i).pdgid() << " " << tau_daughters.at(i).px() << " " 
		<<tau_daughters.at(i).py() << " " <<tau_daughters.at(i).pz() << " " << tau_daughters.at(i).e() << std::endl;
      Tau+=TLorentzVector(tau_daughters.at(i).px(),tau_daughters.at(i).py(),tau_daughters.at(i).pz(),tau_daughters.at(i).e());
    }
    std::cout << "Taureco " << Tau.Px() << " " << Tau.Py() << " " << Tau.Pz() << " " << Tau.E() << " " << Tau.M() <<  std::endl;
    std::cout << "simpleParticle - 2 " << std::endl;
    std::cout << "simpleParticle " << tau2.pdgid() << " " << tau2.px() << " " <<tau2.py() << " " <<tau2.pz() << " " << tau2.e() 
	      << " size " << tau_daughters2.size() <<std::endl;
    TLorentzVector Tau2(0,0,0,0);
    for(unsigned int i=0;i<tau_daughters2.size();i++){
      std::cout << "simpleParticle" << tau_daughters2.at(i).pdgid() << " " << tau_daughters2.at(i).px() << " "
		<<tau_daughters2.at(i).py() << " " <<tau_daughters2.at(i).pz() << " " << tau_daughters2.at(i).e() << std::endl;
      Tau2+=TLorentzVector(tau_daughters.at(i).px(),tau_daughters.at(i).py(),tau_daughters.at(i).pz(),tau_daughters.at(i).e());
    }
    std::cout << "Tau2reco " << Tau2.Px() << " " << Tau2.Py() << " " << Tau2.Pz() << " " << Tau2.E() << " " << Tau2.M() << std::endl;
    */
    WT = calculateWeightFromParticlesH(X, tau, tau2, tau_daughters,tau_daughters2);
    double polSM=getTauSpin();
    //std::cout << "polSM=getTauSpin() " <<  polSM << " " << getTauSpin() << " WT " << WT << std::endl;
  if(type==hminus || type==hplus) WT=1.0;
    if(type==hminus && polSM>0.0) WT=0; // sign definition flipped in TauSpiner
    if(type==hplus && polSM<0.0) WT=0;
  }

  if(Spin==type || type==hplus ||  type==hminus) return WT;
  if(UnSpin==type) return 1.0/WT;
  if(FlipSpin==type) return (2.0-WT)/(WT);
  Logger(Logger::Warning) << "TauSpinerWeight TauSpinerType " << type << " is INVALID. Returning Spin WT=1.0." << std::endl;
  return 1.0;
}

#endif
