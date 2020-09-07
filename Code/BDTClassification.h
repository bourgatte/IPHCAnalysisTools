#ifndef BDTClassification_h
#define BDTClassification_h

#include <string>
#include "boost/format.hpp"
#include "boost/bind.hpp"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

class BDTClassification{
 private:
  
  //int year_;
  //double dijetpt_;
  double jdeta_;
  double jpt_1_;
  double svfit_mass_;
  double m_vis_;
  double met_;
  double mjj_;
  unsigned n_jets_;
  double pt_1_;
  // double pt_2_;
  double pt_tt_;
  double pt_vis_;

  TMVA::Reader *reader_even_;
  TMVA::Reader *reader_odd_;

 public:
  BDTClassification();
  virtual ~BDTClassification();
  virtual std::vector<float> read_mva_scores(unsigned isEven, std::vector<float> vars);
  virtual std::pair<float,int> getMaxScoreWithIndex(std::vector<float> vec);
  virtual int PreAnalysis();
  virtual int Execute(double jdeta,double jpt_1,double m_vis,double met,double mjj,unsigned n_jets,double pt_1, double pt_tt,double pt_vis, double svfit_mass, unsigned long long evt_, std::vector<float> &score, std::pair<float, int> &max_pair);

  //int year_;
  unsigned isEven_;
  float event_;
  unsigned long long evt_;

  /* double dijetpt_; */
  /* double jdeta_; */
  /* double jpt_1_; */
  /* double m_sv_; */
  /* double m_vis_; */
  /* double met_; */
  /* double mjj_; */
  /* unsigned n_jets_; */
  /* double pt_1_; */
  /* double pt_2_; */
  /* double pt_tt_; */
  /* double pt_vis_; */

  float var0_, var1_, var2_, var3_, var4_, var5_, var6_, var7_, var8_, var9_;//, var10_, var11_;

};

#endif
