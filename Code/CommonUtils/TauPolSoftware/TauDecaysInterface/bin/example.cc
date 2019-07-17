#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include "TLorentzVector.h"
#include "TauPolSoftware/TauDecaysInterface/interface/a1Helper.h"
#include "TauPolSoftware/TauDecaysInterface/interface/rhoHelper.h"
#include "TauPolSoftware/TauDecaysInterface/interface/TauPolInterface.h"


int main(int argc, const char* argv[]) {

  std::cout<<"\n\n\n------------------------------ Single Muon --------------------------------------------------" <<std::endl;
  //  Momentum example from pythia 
  // (x,y,z,t)=(4.535169,43.546793,-13.992208,45.998151) (P,eta,phi,E)=(45.963821,-0.314382,1.467026,45.998151)   // tau
  // (x,y,z,t)=(-0.024813,1.699042,-0.626012,1.813950) (P,eta,phi,E)=(1.810870,-0.360548,1.585399,1.813950)  // muon
  
  TauPolInterface  TauSingleMu;
  vector<TLorentzVector> tauandprod_SingleMu;
  tauandprod_SingleMu.push_back(TLorentzVector(4.535169,43.546793,-13.992208,45.998151)); //  order matters! First should be tau
  tauandprod_SingleMu.push_back(TLorentzVector(-0.024813,1.699042,-0.626012,1.813950)); //   muon
  TauSingleMu.Configure(tauandprod_SingleMu,"lepton");
  if(TauSingleMu.isConfigured())std::cout<<" Single Muon: omega =       "<< TauSingleMu.getOmega() <<std::endl;
  
  std::cout<<"\n\n\n------------------------------ Single Pion --------------------------------------------------"<<std::endl;
  //  Momentum example from pythia 
  // (x,y,z,t)=(18.556395,-1.697132,42.019311,45.999996) (P,eta,phi,E)=(45.965667,1.552187,-0.091204,45.999996)// tau
  // (x,y,z,t)=(13.259153,-1.264031,31.915783,34.583806) (P,eta,phi,E)=(34.583524,1.607980,-0.095045,34.583806)// pi
  
  TauPolInterface  TauSinglePi;
  vector<TLorentzVector> tauandprod_SinglePi;
  tauandprod_SinglePi.push_back(TLorentzVector(18.556395,-1.697132,42.019311,45.999996)); //  order matters! First should be tau
  tauandprod_SinglePi.push_back(TLorentzVector(13.259153,-1.264031,31.915783,34.583806)); // pion
  TauSinglePi.Configure(tauandprod_SinglePi,"pion");
  if(TauSinglePi.isConfigured())std::cout<<" Single Pion: omega =   "<< TauSinglePi.getOmega() <<std::endl;
  

  std::cout<<"\n\n\n------------------------------ Single Rho ---------------------------------------------------"<<std::endl;
  //  Momentum example from pythia 
  // (x,y,z,t)=(41.565165,-7.157718,-18.152869,45.951927) (P,eta,phi,E)=(45.917562,-0.418109,-0.170532,45.951927)   -  tau
  // (x,y,z,t)=(13.955560,-2.078894,-5.975415,15.323333) (P,eta,phi,E)=(15.322697,-0.411766,-0.147878,15.323333)    - pi charged  
  // (x,y,z,t)=(5.780525,-0.674187,-2.142884,6.203158) (P,eta,phi,E)=(6.201689,-0.360361,-0.116106,6.203158)        - pi0
  
  TauPolInterface  TauSingleRho;
  vector<TLorentzVector> tauandprod_SingleRho;
  
  tauandprod_SingleRho.push_back(TLorentzVector(41.565165,-7.157718,-18.152869,45.951927)); //  order matters! First should be tau
  tauandprod_SingleRho.push_back(TLorentzVector(13.955560,-2.078894,-5.975415,15.323333)); // pi charged
  tauandprod_SingleRho.push_back(TLorentzVector(5.780525,-0.674187,-2.142884,6.203158)); //  pi0
  TauSingleRho.Configure(tauandprod_SingleRho,"rho");
  if(TauSinglePi.isConfigured())
    {
      std::cout<<" Single Rho: omega    =   "<< TauSingleRho.getOmega() <<std::endl;  // this returns only cosbeta angle (regualr charge-neutral energy assymetry)
      std::cout<<" Single Rho: omegabar =   "<< TauSingleRho.getOmegabar() <<std::endl;  // cosbeta + costheta + cospsi. This is still not the full observable, implementation of the last angle alpha is in progress ...
    }
  
  std::cout<<"\n\n\n------------------------------ Single A1  --------------------------------------------------- "<<std::endl;
  //  Momentum example from pythia 
  // (x,y,z,t)=(18.670063,-37.280665,19.350028,45.999999) (P,eta,phi,E)=(45.965670,0.448867,-1.106511,45.999999)
  // (x,y,z,t)=(1.141691,-2.936690,1.490553,3.488386) (P,eta,phi,E)=(3.485592,0.456996,-1.200010,3.488386)
  // (x,y,z,t)=(6.270693,-11.695517,6.032508,14.577974) (P,eta,phi,E)=(14.577306,0.440222,-1.078639,14.577974)
  // (x,y,z,t)=(7.184014,-14.091224,7.801724,17.636866) (P,eta,phi,E)=(17.636313,0.475170,-1.099322,17.636866)
  
  
  TauPolInterface  TauSingleA1;
  vector<TLorentzVector> tauandprod_SingleA1;
  
  tauandprod_SingleA1.push_back(TLorentzVector(18.670063,-37.280665,19.350028,45.999999));  //  order matters! First should be tau
  tauandprod_SingleA1.push_back(TLorentzVector(1.141691,-2.936690,1.490553,3.488386));    //  OS pi
  tauandprod_SingleA1.push_back(TLorentzVector(6.270693,-11.695517,6.032508,14.577974)); //  SS1 pi 
  tauandprod_SingleA1.push_back(TLorentzVector(7.184014,-14.091224,7.801724,17.636866));   //  SS2 pi  
  
  
  TauSingleA1.Configure(tauandprod_SingleA1,"a1",-1);  // in case of a1 the last argument requires the charge of the tau lepton (sum charge of three pions)
  if(TauSingleA1.isConfigured())
    {
      std::cout<<" Single A1: omega    =   "<< TauSingleA1.getOmega() <<std::endl;     // 
      std::cout<<" Single A1: omegabar =   "<< TauSingleA1.getOmegabar() <<std::endl;  // cosbeta + costheta + cosalpha + cosgamma. 
    }
  
  std::cout<<"\n\n\n------------------------------  Pion - Pion Pair --------------------------------------------"<<std::endl;
  //  Momentum example from pythia 
  //   First _______________
  // (x,y,z,t)=(18.556395,-1.697132,42.019311,45.999996) (P,eta,phi,E)=(45.965667,1.552187,-0.091204,45.999996)  - tau 
  // (x,y,z,t)=(13.259153,-1.264031,31.915783,34.583806) (P,eta,phi,E)=(34.583524,1.607980,-0.095045,34.583806)  - pi
  //   Second  _______________
  // (x,y,z,t)=(-18.556395,1.697132,-42.019309,45.999994) (P,eta,phi,E)=(45.965665,-1.552187,3.050388,45.999994)  - tau 
  // (x,y,z,t)=(-7.890180,-0.016431,-16.596710,18.377309) (P,eta,phi,E)=(18.376779,-1.488969,-3.139510,18.377309)  - pi
  
  
  TauPolInterface  TauPiPiPair;
  vector<TLorentzVector> tauandprod_Pi1;
  vector<TLorentzVector> tauandprod_Pi2;
  
  
  tauandprod_Pi1.push_back(TLorentzVector(18.556395,-1.697132,42.019311,45.999996)); //  order matters! First should be tau
  tauandprod_Pi1.push_back(TLorentzVector(13.259153,-1.264031,31.915783,34.583806)); // pi
  
  tauandprod_Pi2.push_back(TLorentzVector(-18.556395,1.697132,-42.019309,45.999994)); //  order matters! First should be tau
  tauandprod_Pi2.push_back(TLorentzVector(-7.890180,-0.016431,-16.596710,18.377309)); // pi
  
  TauPiPiPair.ConfigurePair(tauandprod_Pi1,"pion",tauandprod_Pi2,"pion");
  if(TauPiPiPair.isPairConfigured())
    {
      std::cout<<" Combined Pi-Pi: omega =   "<< TauPiPiPair.getCombOmega() <<std::endl;
      std::cout<<" Visible Pair Mass     =   "<< TauPiPiPair.getVisiblePairLV().M() <<std::endl;
    }
  
  
  std::cout<<"\n\n\n------------------------------  Pion - Rho Pair ---------------------------------------------"<<std::endl;
  //  Momentum example from pythia 
  //   First _______________
  // (x,y,z,t)=(-41.565165,7.157718,18.221388,45.979038) (P,eta,phi,E)=(45.944693,0.419601,2.971060,45.979038)  - tau
  // (x,y,z,t)=(-10.011427,1.198892,3.770703,10.765858) (P,eta,phi,E)=(10.764953,0.365758,3.022408,10.765858)  - pi 
  //   Second  _______________
  // (x,y,z,t)=(41.565165,-7.157718,-18.152869,45.951927) (P,eta,phi,E)=(45.917562,-0.418109,-0.170532,45.951927)  - tau
  // (x,y,z,t)=(13.955560,-2.078894,-5.975415,15.323333) (P,eta,phi,E)=(15.322697,-0.411766,-0.147878,15.323333)  - pi charged 
  // (x,y,z,t)=(5.780525,-0.674187,-2.142884,6.203158) (P,eta,phi,                                                 - pi0
  
  
  
  TauPolInterface  TauPiRhoPair;
  vector<TLorentzVector> tauandprod_Pi;
  vector<TLorentzVector> tauandprod_Rho;
  
  
  tauandprod_Pi.push_back(TLorentzVector(-41.565165,7.157718,18.221388,45.979038)); //  order matters! First should be tau
  tauandprod_Pi.push_back(TLorentzVector(-10.011427,1.198892,3.770703,10.765858)); // pi 
  
  tauandprod_Rho.push_back(TLorentzVector(41.565165,-7.157718,-18.152869,45.951927)); //  order matters! First should be tau
  tauandprod_Rho.push_back(TLorentzVector(13.955560,-2.078894,-5.975415,15.323333)); //   pi charged 
  tauandprod_Rho.push_back(TLorentzVector(5.780525,-0.674187,-2.142884,6.203158)); //     pi0 
  
  TauPiRhoPair.ConfigurePair(tauandprod_Pi,"pion",tauandprod_Rho,"rho");
  if(TauPiRhoPair.isPairConfigured())
    {
      std::cout<<" Combined Pi-Rho: omega =   "<< TauPiRhoPair.getCombOmegaBar() <<std::endl;
      std::cout<<" Visible Pair Mass      =   "<< TauPiRhoPair.getVisiblePairLV().M() <<std::endl;
    }

  
  std::cout<<"\n\n\n------------------------------  Pion - A1 Pair ----------------------------------------------"<<std::endl;
  //  Momentum example from pythia 
  //   First _______________
  // (x,y,z,t)=(-0.239895,35.808984,28.817160,45.999203) (P,eta,phi,E)=(45.964873,0.736356,1.577496,45.999203)
  // (x,y,z,t)=(0.546745,7.638778,6.673710,10.159121) (P,eta,phi,E)=(10.158163,0.787482,1.499343,10.159121)
  //   Second  _______________
  // (x,y,z,t)=(-0.239895,35.808984,28.817160,45.999203) (P,eta,phi,E)=(45.964873,0.736356,1.577496,45.999203)
  // (x,y,z,t)=(0.045646,-4.407360,-3.941665,5.914652) (P,eta,phi,E)=(5.913005,-0.804616,-1.560440,5.914652)
  // (x,y,z,t)=(0.116904,-18.704540,-14.602361,23.730190) (P,eta,phi,E)=(23.729779,-0.717503,-1.564546,23.730190)
  // (x,y,z,t)=(-0.137735,-7.375450,-5.496942,9.200657) (P,eta,phi,E)=(9.199598,-0.689281,-1.589469,9.200657)
  
  TauPolInterface  TauPiA1Pair;
  vector<TLorentzVector> tauandprod__Pi;
  vector<TLorentzVector> tauandprod_A1;
  
  tauandprod__Pi.push_back(TLorentzVector(-0.239895,35.808984,28.817160,45.999203)); //  order matters! First should be tau
  tauandprod__Pi.push_back(TLorentzVector(0.546745,7.638778,6.673710,10.159121)); //   pi 
  
  tauandprod_A1.push_back(TLorentzVector(18.670063,-37.280665,19.350028,45.999999));  //  order matters! First should be tau
  tauandprod_A1.push_back(TLorentzVector(1.141691,-2.936690,1.490553,3.488386));    //  OS pi
  tauandprod_A1.push_back(TLorentzVector(6.270693,-11.695517,6.032508,14.577974)); //  SS1 pi 
  tauandprod_A1.push_back(TLorentzVector(7.184014,-14.091224,7.801724,17.636866));   //  SS2 pi 


  
  TauPiA1Pair.ConfigurePair(tauandprod__Pi,"pion",tauandprod_A1,"a1",1,-1);
  if(TauPiA1Pair.isPairConfigured())
    {
      std::cout<<" Combined Pi-A1: omega    =   " << TauPiA1Pair.getCombOmega() <<std::endl;
      std::cout<<" Combined Pi-A1: omegaBar =   " << TauPiA1Pair.getCombOmegaBar() <<std::endl;
      std::cout<<" Visible Pair Mass        =    "<< TauPiA1Pair.getVisiblePairLV().M() <<std::endl;
    }
}


