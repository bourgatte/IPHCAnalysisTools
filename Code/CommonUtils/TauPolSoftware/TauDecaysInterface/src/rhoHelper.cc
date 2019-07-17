#include "TauPolSoftware/TauDecaysInterface/interface/rhoHelper.h"
#include <iostream>

rhoHelper::rhoHelper(){
}

rhoHelper::rhoHelper(vector<TLorentzVector> TauRhoandProd)
{
  if(TauRhoandProd.size()!=3){
    std::cout<<" Warning!! Size of rho input vector != 4 !! "<<std::endl;
  }
  TLorentzVector LabFrame = TauRhoandProd.at(1) + TauRhoandProd.at(2);
  Setup(TauRhoandProd, LabFrame);
}

void 
rhoHelper::Setup(vector<TLorentzVector> TauRhoandProd, TLorentzVector ReferenceFrame){
   mpi  = 0.13957018; // GeV 
   mpi0 = 0.1349766;   // GeV
   mtau = 1.77687; // GeV
   coscab = 0.975; 
   mrho = 0.773; // GeV
   debug=false;
   for(unsigned int i=0; i<TauRhoandProd.size(); i++){
     TauRhoandProd_RF.push_back(Boost(TauRhoandProd.at(i),ReferenceFrame));
   }

   TauLV         = TauRhoandProd.at(0);
   TauRhoPi      = TauRhoandProd.at(1);
   TauRhoPi0     = TauRhoandProd.at(2);
   ProductLV     = TauRhoPi+TauRhoPi0;

   DPF_TauRhoPi  = Boost(TauRhoPi,ProductLV);
   DPF_TauRhoPi0 = Boost(TauRhoPi0,ProductLV);
   InvisibleLV = TauLV - ProductLV;
   DPF_TauLV=  Boost(TauLV,ProductLV);

   DPF_InvisibleLV=  Boost(InvisibleLV,ProductLV);

   TVector3 RotVector = DPF_TauLV.Vect();
   DPF_TauLV.SetVect(Rotate(DPF_TauLV.Vect(),RotVector));
   DPF_TauRhoPi.SetVect(Rotate(DPF_TauRhoPi.Vect(),RotVector));
   DPF_TauRhoPi0.SetVect(Rotate(DPF_TauRhoPi0.Vect(),RotVector));
   ProductLV.SetVect(Rotate(ProductLV.Vect(),RotVector)); //  rotate nL
}

void 
rhoHelper::Configure(vector<TLorentzVector> TauRhoandProd){
  if(TauRhoandProd.size()!=3){
    std::cout<<" Warning!! Size of input vector != 4 !! "<<std::endl;
  }
  TLorentzVector LabFrame = TauRhoandProd.at(1) + TauRhoandProd.at(2);
  Setup(TauRhoandProd,LabFrame);
}

bool
rhoHelper::isConfigured(){
  if(TauRhoandProd_RF.size()!=3){ std::cout<<"Error:   rhoHelper is not Configured! Check  the size of input vector!  Size =  "<< TauRhoandProd_RF.size() <<std::endl; return false;} return true;
}

rhoHelper::~rhoHelper(){
}

TLorentzVector 
rhoHelper::Boost(TLorentzVector pB, TLorentzVector frame){
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
TVector3
rhoHelper::Rotate(TVector3 LVec, TVector3 Rot){
  TVector3 vec = LVec;
  vec.RotateZ(0.5*TMath::Pi() - Rot.Phi());  // not 0.5, to avoid warnings about 0 pT
  vec.RotateX(Rot.Theta());
  return vec;
}

double
rhoHelper::getCosbetaRho(){
  double cb=-999;
  cb = (mrho/sqrt(mrho*mrho - 4*mpi*mpi))  *   (TauRhoPi.E() -   TauRhoPi0.E())/(TauRhoPi.E() +  TauRhoPi0.E()) ;
  if(std::fabs(cb) - 1 >0  && std::fabs(cb) < 1.01 &&  cb> 0  )  cb= 1;
  if(std::fabs(cb) - 1 >0  && std::fabs(cb) < 1.01 &&  cb< 0  )  cb=-1;
  if(std::fabs(cb) > 1 ){if(debug){std::cout<<"Warning! Cos beta > 1:  "<<cb <<std::endl;  }}
  return cb;
}

double
rhoHelper::getSinbetaRho(){
  double sb=-999;
  if(std::fabs(getCosbetaRho()) > 1 ){if(debug){std::cout<<"Warning! Cos beta > 1:  "<< getCosbetaRho()<<std::endl;  }return 0;}
  sb = sqrt(1- getCosbetaRho()*getCosbetaRho());
  return sb;
}

double
rhoHelper::getCosthetaRho(){
  double ct=-999;
  double QQ = ProductLV.M2();
  double x = ProductLV.E()/TauLV.E();
  double s = 4*TauLV.E()*TauLV.E();
  if(std::fabs(ct) - 1 >0  && std::fabs(ct) < 1.01 &&  ct> 0  )  ct = 1;
  if(std::fabs(ct) - 1 >0  && std::fabs(ct) < 1.01 &&  ct< 0  )  ct=-1;
  if( 1 - 4*mtau*mtau/s  <= 0 ){std::cout<<"Warning! In costheta root square <=0! return -999"<<std::endl; return ct;}
  ct= (2*x*mtau*mtau - mtau*mtau - QQ)/( (mtau*mtau - QQ)*sqrt(1 - 4*mtau*mtau/s) );
  if(std::fabs(ct) > 1 ){if(debug){std::cout<<"Warning! Cos theta > 1:  "<<ct<<std::endl; }}
  return ct;
}

double
rhoHelper::getSinthetaRho(){
  double st=-999;
   st = sqrt(1- getCosthetaRho()*getCosthetaRho());
   if(std::fabs(getCosthetaRho()) > 1 ){if(debug){std::cout<<"Warning! Cos theta > 1"<<std::endl; }}
  return st;
}

double 
rhoHelper::getUltrarel_cospsiLF(){
  double cos=-999;
  double QQ = ProductLV.M2();
  cos = (getCosthetaRho()*(mtau*mtau  + QQ)   + (mtau*mtau  - QQ))/(getCosthetaRho()*(mtau*mtau  - QQ)   + (mtau*mtau  + QQ));
  if(std::fabs(cos) - 1 >0  && std::fabs(cos) < 1.01 && cos > 0  ) cos = 1;
  if(std::fabs(cos) - 1 >0  && std::fabs(cos) < 1.01 && cos < 0  ) cos =-1;
  if(std::fabs(cos) > 1 )if(debug){std::cout<<"Warning! Cos psi > 1:  "<<cos<<std::endl; }
  return cos;
}

double 
rhoHelper::getSinpsiLF(){
  double sin = -999;
  sin = sqrt(1 - getUltrarel_cospsiLF()*getUltrarel_cospsiLF());
  if(getUltrarel_cospsiLF()*getUltrarel_cospsiLF() > 1  ){if(debug){std::cout<<"Warning! getultrarel_cospsiLF > 1"<<std::endl;}}
  return    sin;
}

TMatrixT<double> rhoHelper::convertToMatrix(TVectorT<double> V){
  TMatrixT<double> M(V.GetNrows(),1);
  for(int i=0; i < M.GetNrows(); i++){
    M(i,0)=V(i);
  } return M;
}
 
TLorentzVector
rhoHelper::sLV(){
  double QQ = ProductLV.M2();
  double l0 = 0.5*(mtau*mtau + QQ)/sqrt(QQ);
  double l   = 0.5*(mtau*mtau  - QQ)/sqrt(QQ);
  TLorentzVector tspin(getSinthetaRho(),0,-l0*getCosthetaRho()/mtau,-l*getCosthetaRho()/mtau);
  return tspin;
}

TVector3
rhoHelper::nPerp(){
  return   DPF_TauRhoPi.Vect()*(1/DPF_TauRhoPi.Vect().Mag());
}

TVector3
rhoHelper::ns(){  
  return   sLV().Vect()*(1/sLV().Vect().Mag());
}

TVector3
rhoHelper::nL(){
  return   -ProductLV.Vect()*(1/ProductLV.Vect().Mag());
}
 
TVector3
rhoHelper::nT(){
  return   DPF_TauLV.Vect()*(1/DPF_TauLV.Vect().Mag());
}

TVector3
rhoHelper::nTRotated(){
  TVector3 vec= DPF_TauLV.Vect();
  vec.RotateZ(0.5*TMath::Pi() - vec.Phi());
  vec.RotateX(vec.Theta());
  return   vec;
}
TVector3
rhoHelper::nPerpRotated(){
  TVector3 vec= DPF_TauRhoPi.Vect();
  vec.RotateZ(0.5*TMath::Pi() - DPF_TauLV.Vect().Phi());
  vec.RotateX(DPF_TauLV.Vect().Theta());
  return   vec;
}

double 
rhoHelper::TFK_cosbeta(){
    return nT()*nPerp();
}

double 
rhoHelper::TFK_sinbeta(){
  if(std::fabs(TFK_cosbeta())> 1){std::cout<<" TFK_cosbeta > 1"<< std::endl; return 0;}
  return sqrt(1-TFK_cosbeta()*TFK_cosbeta());
}

double rhoHelper::TFK_costheta(){
  return -ns()*nT();
}

double rhoHelper::TFK_sintheta(){
  if(std::fabs(TFK_costheta())> 1){std::cout<<" TFK_costheta > 1"<< std::endl; return 0;}
  return sqrt(1-TFK_costheta()*TFK_costheta());
}
 
double  
rhoHelper::DPF_cosalpha(){
     TVector3 nTCrossns  = nT().Cross(ns());
     TVector3 nTCrossnPerp  = nT().Cross(nPerp());

    if(nTCrossns.Mag() ==0 || nTCrossnPerp.Mag() ==0){if(debug){std::cout<<" Can not compute cos alpha, one denominator is 0, return DPF cos alpha =0  "<< std::endl; }return 0;}
   return nTCrossns.Dot(nTCrossnPerp)/nTCrossns.Mag()/nTCrossnPerp.Mag();
}

double  
rhoHelper::DPF_sinalpha(){
    TVector3 nTCrossns  = nT().Cross(ns());
    TVector3 nTCrossnPerp  = nT().Cross(nPerp());
    if(nTCrossns.Mag() ==0 || nTCrossnPerp.Mag() ==0){if(debug){std::cout<<" Can not compute sin alpha, one denominator is 0, return DPF sin alpha =0  "<< std::endl; }return 0;}
    return -ns().Dot(nTCrossnPerp)/nTCrossns.Mag()/nTCrossnPerp.Mag();
 }

//----------------------  angles beta  + theta ----------
double 
rhoHelper::getOmegaRho(){
  double omega;
  double QQ =  ProductLV.M2();
  double Be = 0.5*(3*getCosbetaRho()*getCosbetaRho() -1);
  double Ps = 0.5*(3*getUltrarel_cospsiLF() -1);
  double RR = mtau*mtau/QQ;
  double R = sqrt(RR);
  omega = ((-2 + RR + 2*(1+RR)*Ps*Be)*getCosthetaRho() + 3*R*Be*getSinthetaRho()*2*getUltrarel_cospsiLF()*getSinpsiLF() )  /  ( 2 +RR - 2*(1-RR)*Ps*Be);
  if(std::isinf(std::fabs(omega)) || std::isnan(std::fabs(omega))) omega  = -999.;
  return omega;
}

//----------------------  angles beta  + theta  + alpha ----------
double 
rhoHelper::getOmegaRhoBar(){
 TVector3 Rot1 = TauLV.Vect();
  TLorentzVector tauLabR1    = TauLV;
  TLorentzVector piLabR1     = TauRhoPi;
  TLorentzVector pi0LabR1    = TauRhoPi0;
  tauLabR1.SetVect(Rotate(tauLabR1.Vect(),Rot1));
  piLabR1.SetVect(Rotate(piLabR1.Vect(),Rot1));
  pi0LabR1.SetVect(Rotate(pi0LabR1.Vect(),Rot1));
  
  TLorentzVector tauTau= Boost(tauLabR1,tauLabR1);
  TLorentzVector piTau= Boost(piLabR1,tauLabR1);
  TLorentzVector pi0Tau= Boost(pi0LabR1,tauLabR1);
  
  TLorentzVector q= piTau  - pi0Tau;
  TLorentzVector P= tauTau;
  TLorentzVector N= tauTau - piTau - pi0Tau;
  TVector3 h = P.M()*(2*(q*N)*q.Vect() - q.Mag2()*N.Vect()) * (1/ (2*(q*N)*(q*P) - q.Mag2()*(N*P)));
  h *= (1.0 / h.Mag()); // normalisation needed at reconstruction level
  TVector3 TauLabDir =tauLabR1.Vect()*(1/tauLabR1.Vect().Mag());
  return h*TauLabDir;
}
