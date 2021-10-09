#include "BDTClassification.h"

BDTClassification::BDTClassification(){
  ;
}

BDTClassification::~BDTClassification() {
  ;
}

std::vector<float> BDTClassification::read_mva_scores(unsigned isEven, std::vector<float> vars) {
  std::vector<float> score = {};

  var0_=vars[0], var1_=vars[1], var2_=vars[2], var3_=vars[3], var4_=vars[4], var5_=vars[5], var6_=vars[6], var7_=vars[7], var8_=vars[8], var9_=vars[9];

  if(isEven) score = reader_even_->EvaluateMulticlass("Multi"); 
  else       score = reader_odd_->EvaluateMulticlass("Multi");

  return score;
}

// Add this function to take the maximum score and corresponding index
// from the vector of scores
std::pair<float,int> BDTClassification::getMaxScoreWithIndex(std::vector<float> vec) {
  if (vec.empty()) return std::make_pair(0., 0);
  float max_score = vec[0];
  int max_index = 0;
  unsigned ind = 0;
  for (auto s : vec) {
    if (s > max_score) {
      max_score = s;
      max_index = ind;
    }
    ++ind;
  }
  std::pair<float, int> out_pair = std::make_pair(max_score, max_index);
  return out_pair;
}

int BDTClassification::PreAnalysis() {

  reader_even_ = new TMVA::Reader();
  reader_odd_ = new TMVA::Reader();

  // fold0 is trained on even, so apply on odd, and vice versa

  TString filename_even = "";
  TString filename_odd  = "";
  //  if (year_ == 2016) {
    filename_even = (std::string)std::getenv("workdir")+"Code/BDT/multi_fold1_sm_tt_tauspinner_2016_xgb.xml"; // apply to even here
    filename_odd  = (std::string)std::getenv("workdir") + "Code/BDT/multi_fold0_sm_tt_tauspinner_2016_xgb.xml"; // apply to odd
  // } 
  // if (year_ == 2017) {
  //   filename_even = (std::string)std::getenv("workdir")+"Code/BDT/multi_fold1_sm_tt_tauspinner_2017_xgb.xml"; // apply to even here
  //   filename_odd  = (std::string)std::getenv("workdir") + "Code/BDT/multi_fold0_sm_tt_tauspinner_2017_xgb.xml"; // apply to odd
  // } 
  // else if (year_ == 2018) {
  //   filename_even = (std::string)std::getenv("workdir") + "Code/BDT/multi_fold1_sm_tt_tauspinner_2018_xgb.xml"; // apply to even here
  //   filename_odd  = (std::string)std::getenv("workdir") + "Code/BDT/multi_fold0_sm_tt_tauspinner_2018_xgb.xml"; // apply to odd
  // }

     reader_even_->AddVariable( "jdeta",   & var0_ );
     reader_even_->AddVariable( "jpt_1",   & var1_ );
     reader_even_->AddVariable( "m_vis",   & var2_ );
     reader_even_->AddVariable( "met",     & var3_ );
     reader_even_->AddVariable( "mjj",     & var4_ );
     reader_even_->AddVariable( "n_jets",  & var5_ );
     reader_even_->AddVariable( "pt_1",    & var6_ );
     reader_even_->AddVariable( "pt_tt",   & var7_);
     reader_even_->AddVariable( "pt_vis",  & var8_);
     reader_even_->AddVariable( "svfit_mass",    & var9_ );
    
     reader_odd_->AddVariable( "jdeta",    & var0_ );
     reader_odd_->AddVariable( "jpt_1",    & var1_ );
     reader_odd_->AddVariable( "m_vis",    & var2_ );
     reader_odd_->AddVariable( "met",      & var3_ );
     reader_odd_->AddVariable( "mjj",      & var4_ );
     reader_odd_->AddVariable( "n_jets",   & var5_ );
     reader_odd_->AddVariable( "pt_1",     & var6_ );
     reader_odd_->AddVariable( "pt_tt",    & var7_);
     reader_odd_->AddVariable( "pt_vis",   & var8_);
     reader_odd_->AddVariable( "svfit_mass",     & var9_ );
  
     reader_even_->BookMVA( "Multi", filename_even );
     reader_odd_->BookMVA( "Multi", filename_odd );
     
     return 0;
}

int BDTClassification::Execute(double jdeta,double jpt_1,double m_vis,double met,double mjj,unsigned n_jets,double pt_1,double pt_tt,double pt_vis,double svfit_mass, unsigned long long evt_, std::vector<float> &score, std::pair<float, int> &max_pair) {

  //EventInfo const* eventInfo = event->GetPtr<EventInfo>("eventInfo");
  isEven_ = evt_ % 2 == 0; // if even then event_ = 1, odd = 0
  //evt_ = Ntp->EventNumber();
  event_ = (float)isEven_;

  // if (n_lowpt_jets_ >= 1) jpt_1_ = lowpt_jets[0]->pt();
  // if (n_lowpt_jets_ >= 2) {
    // dijetpt_ =  (lowpt_jets[0]->vector() + lowpt_jets[1]->vector()).pt();
    // jdeta_ = fabs(lowpt_jets[0]->eta() - lowpt_jets[1]->eta());
    // mjj_ = (lowpt_jets[0]->vector() + lowpt_jets[1]->vector()).M();
  //}

  //year_=year;
  jdeta_=jdeta;
  jpt_1_=jpt_1;
  mjj_=mjj;
  svfit_mass_=svfit_mass;
  m_vis_=m_vis;
  met_=met;
  n_jets_=n_jets;
  pt_1_=pt_1;
  pt_tt_=pt_tt;
  pt_vis_=pt_vis;
  
  std::vector<float> inputs = {};
  //if (channel_ == channel::tt) {
    inputs.resize(10);
    inputs[0]  = float(jdeta_);
    inputs[1]  = float(jpt_1_);
    inputs[2]  = float(m_vis_);
    inputs[3]  = float(met_);
    inputs[4]  = float(mjj_);
    inputs[5]  = unsigned(n_jets_);
    inputs[6]  = float(pt_1_);
    inputs[7] = float(pt_tt_);
    inputs[8] = float(pt_vis_);
    inputs[9] = float(svfit_mass_);
    // }

  score = read_mva_scores(isEven_,inputs);
  //std::cout<< score[0]<< " "<<score[1]<<"  "<<score[2]<<std::endl;
  // event->Add("higgs_score",    score[0]);
  // event->Add("jetFakes_score", score[1]);
  // event->Add("zttEmbed_score", score[2]);

  max_pair = getMaxScoreWithIndex(score);
  //std::cout<<"Max Score: "<< max_pair.first<<" Index: "<<max_pair.second<<std::endl;
  // event->Add("IC_BDT_max_score", max_pair.first); 
  // event->Add("IC_BDT_max_index", max_pair.second);

  return 0;
}
