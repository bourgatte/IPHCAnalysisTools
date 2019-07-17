
#include "TH1.h"
#include "TF1.h"

void SensitivityScan(){

  TFile *_file0 = TFile::Open("Combined.root");
  double pol=-0.0;

  TH1F *plus = (TH1F *)_file0->Get("Omegapipi_plus");
  TH1F *minus = (TH1F *)_file0->Get("Omegapipi_minus");
 
  plus->Scale(1./plus->Integral());
  minus->Scale(1./minus->Integral());

  TH1F *hf = (TH1F*)plus->Clone("hf");
  TH1F *hg = (TH1F*)plus->Clone("hg");

  hf->Add(minus,1);
  hg->Add(minus,-1);

  hf->Scale(0.5);  
  hg->Scale(0.5);
 
  hf->Add(hg,pol);
  hg->Multiply(hg);

  hg->Divide(hf);
   
  hf->Draw();
  hg->Draw("same");

  std::cout<<"Sens =  " << sqrt(hg->Integral()) << std::endl;

}
