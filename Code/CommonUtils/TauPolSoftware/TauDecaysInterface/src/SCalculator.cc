#include "TauPolSoftware/TauDecaysInterface/interface/SCalculator.h"
#include <iostream>

SCalculator::SCalculator(){
}

SCalculator::SCalculator(string type):
  type_(type)
{
}

void 
SCalculator::Configure(vector<TLorentzVector> TauAndProd, TLorentzVector Frame, int charge){
   for(unsigned int i=0; i<TauAndProd.size(); i++){
     TauAndProd_HRF.push_back(Boost(TauAndProd.at(i), Frame));
   }
   charge_=charge;
}

bool
SCalculator::isConfigured(){
  if(TauAndProd_LF.size()!=2){ std::cout<<"Error:   SCalculator is not Configured! Check  the size of input vector!  Size =  "<< TauAndProd_LF.size() <<std::endl; return false;} return true;
}

SCalculator::~SCalculator(){
}

TLorentzVector 
SCalculator::Boost(TLorentzVector pB, TLorentzVector frame){
   TMatrixT<double> transform(4,4);
   TMatrixT<double> result(4,1);
   TVectorT<double> vec(4); 
   TVector3 b;
   if(frame.Vect().Mag()==0){ std::cout<<"RH Boost is not set, perfrom calculation in the Lab Frame   "<<std::endl; return pB;}
    if(frame.E()==0){ std::cout<<" Caution: Please check that you perform boost correctly!  " <<std::endl; return pB;} 
   else   b=frame.Vect()*(1/frame.E());
   vec(0)  = pB.E();    vec(1)  = pB.Px();
   vec(2)  = pB.Py();  vec(3)  = pB.Pz();
   double gamma  = 1/sqrt( 1 - b.Mag2());
   transform(0,0)=gamma; transform(0,1) =- gamma*b.X() ;  transform(0,2) =  - gamma*b.Y();  transform(0,3) = - gamma*b.Z(); 
   transform(1,0)=-gamma*b.X(); transform(1,1) =(1+ (gamma-1)*b.X()*b.X()/b.Mag2()) ;  transform(1,2) = ((gamma-1)*b.X()*b.Y()/b.Mag2());  transform(1,3) = ((gamma-1)*b.X()*b.Z()/b.Mag2());
   transform(2,0)=-gamma*b.Y(); transform(2,1) = ((gamma-1)*b.Y()*b.X()/b.Mag2());  transform(2,2) = (1 + (gamma-1)*b.Y()*b.Y()/b.Mag2());  transform(2,3) =  ((gamma-1)*b.Y()*b.Z()/b.Mag2()); 
   transform(3,0)=-gamma*b.Z(); transform(3,1) =((gamma-1)*b.Z()*b.X()/b.Mag2()) ;  transform(3,2) = ((gamma-1)*b.Z()*b.Y()/b.Mag2());  transform(3,3) = (1 + (gamma-1)*b.Z()*b.Z()/b.Mag2()); 
   result=transform*convertToMatrix(vec);
   return TLorentzVector(result(1,0), result(2,0) ,result(3,0), result(0,0));
}
TMatrixT<double> SCalculator::convertToMatrix(TVectorT<double> V){
  TMatrixT<double> M(V.GetNrows(),1);
  for(int i=0; i < M.GetNrows(); i++){
    M(i,0)=V(i);
  } return M;
}

TVector3
SCalculator::Rotate(TVector3 LVec, TVector3 Rot){
  TVector3 vec = LVec;
  vec.RotateZ(0.5*TMath::Pi() - Rot.Phi());  // not 0.5, to avoid warnings about 0 pT
  vec.RotateX(Rot.Theta());
  return vec;
}

TVector3
SCalculator::pv(){
  TVector3 out(0,0,0);
  if(type_=="pion") out = TauAndProd_HRF.at(1).Vect();
  if(type_=="rho"){
    TLorentzVector pi  = TauAndProd_HRF.at(1);
    TLorentzVector pi0 = TauAndProd_HRF.at(2);
    TLorentzVector Tau = TauAndProd_HRF.at(0);
    TLorentzVector q= pi  - pi0;
    TLorentzVector P= Tau;
    TLorentzVector N= Tau - pi - pi0;
    out = P.M()*(2*(q*N)*q.Vect() - q.Mag2()*N.Vect()) * (1/ (2*(q*N)*(q*P) - q.Mag2()*(N*P)));
  }
  if(type_=="a1"){
    PolarimetricA1  a1pol;
    a1pol.Configure(TauAndProd_HRF,charge_);
    out = a1pol.PVC().Vect();
  }
  return out;
}

void SCalculator::SortPions(std::vector<TLorentzVector>& pionsvec, std::vector<double>& charges)
{
  
  int npim(0),npip(0);
  int OSMCPionIndex(0);
  int SSMCPion1Index(0);
  int SSMCPion2Index(0);
  int OSCharge(0);
  int SS1Charge(0);
  int SS2Charge(0);
 
  TLorentzVector os;
  TLorentzVector ss1;
  TLorentzVector ss2;
    for(unsigned int i=0; i<charges.size(); i++){
      if( charges.at(i)== 1) npip++;
      if( charges.at(i)==-1) npim++;
    }
    if(npip == 1 && npim == 2){
      int nss=0;
      for(unsigned int i=0; i<charges.size(); i++){
	if(charges.at(i)== 1){
	  OSCharge=1;
	  OSMCPionIndex=i;
	}
	if(charges.at(i)== -1 && nss == 0){
	  nss++;
	  SS1Charge=-1;
	  SSMCPion1Index=i;
	}
	if(charges.at(i)== -1 && nss == 1){
	  SS2Charge=-1;
	  SSMCPion2Index=i;
	}
      }
    }
    if( npip== 2 && npim==1){
      int nss=0;
      for(unsigned int i=0; i<charges.size(); i++){
	if(charges.at(i) == -1){
	  OSCharge=-1;
	  OSMCPionIndex=i;
	}
	if(charges.at(i) == 1 && nss ==0){
	  nss++;
	  SS1Charge=1;
	  SSMCPion1Index=i;
	}
	if(charges.at(i) == 1 && nss == 1){
	  SS2Charge=1;
	  SSMCPion2Index=i;
	}
      }
    }
    
    os=pionsvec.at(OSMCPionIndex);
    ss1=pionsvec.at(SSMCPion1Index);
    ss2=pionsvec.at(SSMCPion2Index);
    
    charges.clear();
    charges.push_back(OSCharge);
    charges.push_back(SS1Charge);
    charges.push_back(SS2Charge);
    
    pionsvec.clear();
    pionsvec.push_back(os);
    pionsvec.push_back(ss1);
    pionsvec.push_back(ss2);
    //return pionsvec;
}
