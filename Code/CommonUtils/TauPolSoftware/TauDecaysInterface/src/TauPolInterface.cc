#include "TauPolSoftware/TauDecaysInterface/interface/TauPolInterface.h"

#include <iostream>

#include <TVector3.h>

TauPolInterface::TauPolInterface():
  debug(false)
{ }  
TauPolInterface::TauPolInterface(vector<TLorentzVector> TauAndProd, string type, int TauCharge):
  debug(false)
 {
  if(type=="lepton" || type =="pion"){
    if(TauAndProd.size()!=2){ 
      std::cout<<" Warning!! Size of input vector  !=2  !! "<< " type:  "<<type<<std::endl;
    }
  }
  if(type=="rho"){
    if(TauAndProd.size()!=3){
      std::cout<<" Warning!! Size of input vector  !=3  !! "<< " type:  "<<type<<std::endl;
    }
  }
  if(type=="a1"){
    if(TauAndProd.size()!=4){
      std::cout<<" Warning!! Size of 2nd  input vector  !=4 !! "<< " type:  "<<type<<std::endl;
    }
  }
  
  Configure(TauAndProd,  type,TauCharge);
 }
TauPolInterface::TauPolInterface(vector<TLorentzVector> TauAndProd1, string type1,vector<TLorentzVector> TauAndProd2, string type2,int TauCharge1,int TauCharge2):
  debug(false)
 {
  if(type1=="lepton" || type1 =="pion"){
    if(TauAndProd1.size()!=2){ 
      std::cout<<" Warning!! Size of 1st  input vector  !=2  !! "<< " type:  "<<type1<<std::endl;
    }
  }
  if(type1=="rho"){
    if(TauAndProd1.size()!=3){
      std::cout<<" Warning!! Size of 1st  input vector  !=3  !! "<< " type:  "<<type1<<std::endl;
    }
  }
  if(type1=="a1"){
    if(TauAndProd1.size()!=4){
      std::cout<<" Warning!! Size of 2nd  input vector  !=4 !! "<< " type:  "<<type1<<std::endl;
    }
  }
  if(type2=="lepton" || type2 =="pion"){
    if(TauAndProd2.size()!=2){ 
      std::cout<<" Warning!! Size of 2nd  input vector  !=2  !! "<< " type:  "<<type2<<std::endl;
    }
  }
  if(type2=="rho"){
    if(TauAndProd2.size()!=3){
      std::cout<<" Warning!! Size of 2nd  input vector  !=3  !! "<< " type:  "<<type2<<std::endl;
    }
  }
  if(type2=="a1"){
    if(TauAndProd2.size()!=4){
      std::cout<<" Warning!! Size of 2nd  input vector  !=4 !! "<< " type:  "<<type2<<std::endl;
    }
  }
  
  ConfigurePair(TauAndProd1,type1,TauAndProd2,type2,TauCharge1,TauCharge2);
 }

void 
TauPolInterface::Configure(vector<TLorentzVector> TauAndProd, string type,int TauCharge){
  TauLV = TauAndProd.at(0);
  type_=type;
  if(type_ == "lepton" || type_=="pion"){
    ProductLV= TauAndProd.at(1);
  }
  if(type_ == "rho" ){
    TauRhoPi  = TauAndProd.at(1);
    TauRhoPi0= TauAndProd.at(2);
    ProductLV = TauRhoPi+TauRhoPi0;
    DPF_TauRhoPi = Boost(TauRhoPi,ProductLV);
    DPF_TauRhoPi0 =  Boost(TauRhoPi0,ProductLV);
  }
  if(type_=="a1"){
    TauA1OSPi  = TauAndProd.at(1);
    TauA1SSPi1 = TauAndProd.at(2);
    TauA1SSPi2 = TauAndProd.at(3);
    ProductLV = TauA1OSPi + TauA1SSPi1 + TauA1SSPi2;  
    taucharge_ = TauCharge;
  }
  InvisibleLV = TauLV - ProductLV;
  DPF_TauLV=  Boost(TauLV,ProductLV);
  DPF_InvisibleLV=  Boost(InvisibleLV,ProductLV);
  //  std::cout<<"Configure type 1 "<< type_; ProductLV.Print();

}

void 
TauPolInterface::ConfigurePair(vector<TLorentzVector> TauAndProd1, string type1,  vector<TLorentzVector> TauAndProd2, string type2,int TauCharge1,int TauCharge2){
  TauLV1 = TauAndProd1.at(0);
  type1_=type1;
  if(type1_ == "lepton" || type1_=="pion"){
    ProductLV1= TauAndProd1.at(1);
  }
  if(type1_ == "rho" ){
    TauRhoPi1      = TauAndProd1.at(1);
    TauRhoPi01     = TauAndProd1.at(2);
    ProductLV1     = TauRhoPi1+TauRhoPi01;
    DPF_TauRhoPi1  = Boost(TauRhoPi1,ProductLV1);
    DPF_TauRhoPi01 = Boost(TauRhoPi01,ProductLV1);
  }

  if(type1_=="a1"){
    TauA1OSPi1  = TauAndProd1.at(1);
    TauA1SSPi11  = TauAndProd1.at(2);
    TauA1SSPi21  = TauAndProd1.at(3);
    ProductLV1 = TauA1OSPi1 + TauA1SSPi11 + TauA1SSPi21;  
    taucharge_=TauCharge1;
  }
  InvisibleLV1 = TauLV1 - ProductLV1;
  DPF_TauLV1=  Boost(TauLV1,ProductLV1);
  DPF_InvisibleLV1=  Boost(InvisibleLV1,ProductLV1);


  TauLV2 = TauAndProd2.at(0);
  type2_=type2;
  if(type2_ == "lepton" || type2_=="pion"){
    ProductLV2= TauAndProd2.at(1);
  }
  if(type2_ == "rho" ){
    TauRhoPi2      = TauAndProd2.at(1);
    TauRhoPi02     = TauAndProd2.at(2);
    ProductLV2     = TauRhoPi2+TauRhoPi02;
    DPF_TauRhoPi2  = Boost(TauRhoPi2,ProductLV2);
    DPF_TauRhoPi02 = Boost(TauRhoPi02,ProductLV2);
  }
  if(type2_=="a1"){
    TauA1OSPi2  = TauAndProd2.at(1);
    TauA1SSPi12  = TauAndProd2.at(2);
    TauA1SSPi22  = TauAndProd2.at(3);
    ProductLV2 = TauA1OSPi2 + TauA1SSPi12 + TauA1SSPi22;  
    taucharge_=TauCharge2;
  }
  InvisibleLV2 = TauLV2 - ProductLV2;
  DPF_TauLV2=  Boost(TauLV2,ProductLV2);
  DPF_InvisibleLV2=  Boost(InvisibleLV2,ProductLV2);

  TauLV2 = TauAndProd2.at(0);
  type2_=type2;

}


void 
TauPolInterface::SetupLeg(string which){
    if(which=="first")
    {
      TauLV = TauLV1;
      type_=type1_;
      if(type_ == "lepton" || type_=="pion" ){ProductLV= ProductLV1;  }


      if(type_ == "rho" )
	{
	  TauRhoPi      = TauRhoPi1;
	  TauRhoPi0     = TauRhoPi01;
	  ProductLV     = TauRhoPi+TauRhoPi0;
	  DPF_TauRhoPi  = Boost(TauRhoPi,ProductLV);
	  DPF_TauRhoPi0 = Boost(TauRhoPi0,ProductLV);
	}
      if(type_=="a1")
	{
	  TauA1OSPi  =   TauA1OSPi1;
	  TauA1SSPi1  =  TauA1SSPi11;
	  TauA1SSPi2  =  TauA1SSPi21;
	  ProductLV = TauA1OSPi + TauA1SSPi1 + TauA1SSPi2;  
	}
      InvisibleLV = TauLV - ProductLV;
      DPF_TauLV=  Boost(TauLV,ProductLV);
      DPF_InvisibleLV=  Boost(InvisibleLV,ProductLV);
      //      std::cout<<"type 1 "<< type_ <<std::endl; ProductLV.Print();
    }
  if(which=="second")
    {
      TauLV = TauLV2;
      type_=type2_;
      if(type_ == "lepton" || type_=="pion" ){ProductLV= ProductLV2;  }
      if(type_ == "rho" )
	{
	  TauRhoPi      = TauRhoPi2;
	  TauRhoPi0     = TauRhoPi02;
	  ProductLV     = TauRhoPi+TauRhoPi0;
	  DPF_TauRhoPi  = Boost(TauRhoPi,ProductLV);
	  DPF_TauRhoPi0 = Boost(TauRhoPi0,ProductLV);
	}

      if(type_=="a1")
	{
	  TauA1OSPi  =   TauA1OSPi2;
	  TauA1SSPi1  =  TauA1SSPi12;
	  TauA1SSPi2  =  TauA1SSPi22;
	  ProductLV = TauA1OSPi + TauA1SSPi1 + TauA1SSPi2;  
	}

      InvisibleLV = TauLV - ProductLV;
      DPF_TauLV=  Boost(TauLV,ProductLV);
      DPF_InvisibleLV=  Boost(InvisibleLV,ProductLV);
      //      std::cout<<"type 2 "<< type_ <<std::endl; ProductLV.Print();

    }
}

bool  
TauPolInterface::isConfigured(){if(TauLV.E()!=0 && ProductLV.E()!=0) return true; return false;}

bool  
TauPolInterface::isPairConfigured(){if(TauLV1.E()!=0 && ProductLV1.E()!=0 && TauLV2.E()!=0 && ProductLV2.E()!=0) return true; return false;}


TauPolInterface::~TauPolInterface(){}



TLorentzVector 
TauPolInterface::Boost(TLorentzVector pB, TLorentzVector frame)
   {
     TVector3 boostVector = frame.BoostVector();
     TLorentzVector result(pB);
     result.Boost(boostVector.X(), boostVector.Y(), boostVector.Z());
     return result;
   }

double
TauPolInterface::getOmega(string which)
   { 
     double omega=-999;
     SetupLeg(which);
     if(type_=="pion" || type_=="lepton") // theta
       {
	 omega = 2.0 * std::min(ProductLV.E()/TauLV.E(), 1.0) - 1.0;
       }
     if(type_=="rho") // beta + theta
       {
	 std::vector<TLorentzVector> particles; // tau, pi, pi0
	 particles.push_back(TauLV);
	 particles.push_back(TauRhoPi);
	 particles.push_back(TauRhoPi0);
	 rhoHelper rho;
	 rho.Configure(particles);
	 omega=rho.getOmegaRho();
       }
     if(type_=="a1") // beta + theta + gamma (Kuehn model)
       {
	 std::vector<TLorentzVector> particles;// tau, os,ss1,ss2
	 particles.push_back(TauLV);
	 particles.push_back(TauA1OSPi);
	 particles.push_back(TauA1SSPi1);
	 particles.push_back(TauA1SSPi2);
	 a1Helper a1;
	 a1.Configure(particles,taucharge_);
 	 omega = a1.getOmegaA1();
       }
     return omega;
   }

double
TauPolInterface::getVisibleOmega(string which)
   { 
     double omega(-999.);
     SetupLeg(which);
     if(type_!="rho") {
       //std::cout<<"This observable  is available for rho decay only (beta angle). a1 will be implemented later" << std::endl;
     }
     else if(type_=="rho") // charged-neutral energy asymmetry
       {
	 std::vector<TLorentzVector> particles; // tau, pi, pi0
	 particles.push_back(TauLV);
	 particles.push_back(TauRhoPi);
	 particles.push_back(TauRhoPi0);
	 rhoHelper rho;
	 rho.Configure(particles);
	 omega=rho.getCosbetaRho();
       }
     return omega;
   }


double
TauPolInterface::getOmegabar(string which){
  double omega=-999;
  SetupLeg(which);

  if(type_=="pion" || type_=="lepton"){ // theta
        omega = 2.0 * std::min(ProductLV.E()/TauLV.E(), 1.0) - 1.0;
  }
  if(type_=="rho") // beta + theta + alpha
      {
	std::vector<TLorentzVector> particles; // tau, pi, pi0
	particles.push_back(TauLV);
	particles.push_back(TauRhoPi);
	particles.push_back(TauRhoPi0);
	rhoHelper rho;
	rho.Configure(particles);
	omega=rho.getOmegaRhoBar();
      }

     if(type_=="a1") // beta + theta + gamma + alpha (CLEO model)
       {
	 std::vector<TLorentzVector> particles;// tau, os,ss1,ss2
	 particles.push_back(TauLV);
	 particles.push_back(TauA1OSPi);
	 particles.push_back(TauA1SSPi1);
	 particles.push_back(TauA1SSPi2);
	 a1Helper a1;
	 a1.Configure(particles,taucharge_);
	 omega = a1.getOmegaA1Bar();
       }
  return omega;
}


double
TauPolInterface::getCombOmega(){
  double omega1 = getOmega("first");
  double omega2 = getOmega("second");
  double Omega=999.;
  Omega = (omega1 + omega2)/(1 + omega1*omega2);
  if(  std::isinf(std::fabs(Omega)) ||  std::isnan(std::fabs(Omega))) Omega  = -999.;
  return Omega;
}



double
TauPolInterface::getCombOmegaBar(){
  double omega1 = getOmegabar("first");
  double omega2 = getOmegabar("second");
  double Omega=999.;
  Omega = (omega1 + omega2)/(1 + omega1*omega2);
  if(  std::isinf(std::fabs(Omega)) ||  std::isnan(std::fabs(Omega))) Omega  = -999.;
  return Omega;
}

double
TauPolInterface::getCombVisibleOmega(){
  double omega1 = getVisibleOmega("first");
  double omega2 = getVisibleOmega("second");
  double Omega=999.;
  if ((omega1 > -990.) && (omega2 > -990.))
  {
    Omega = (omega1 + omega2)/(1 + omega1*omega2);
  }
  if(  std::isinf(std::fabs(Omega)) ||  std::isnan(std::fabs(Omega))) Omega  = -999.;
  return Omega;
}


TLorentzVector
TauPolInterface::getVisiblePairLV()
   { 
     if(!isPairConfigured()) return TLorentzVector(0,0,0,0);
     return ProductLV1+ProductLV2;
   }

TMatrixT<double> TauPolInterface::convertToMatrix(TVectorT<double> V){
  TMatrixT<double> M(V.GetNrows(),1);
  for(int i=0; i < M.GetNrows(); i++){
    M(i,0)=V(i);
  } return M;
}
