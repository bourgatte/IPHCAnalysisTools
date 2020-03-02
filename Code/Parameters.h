#ifndef Parameters_h
#define Parameters_h

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
#include <math.h>

#include <vector>
#include "TString.h"

class Parameters {

 public:
  Parameters();
  Parameters(TString f);
  virtual ~Parameters();

  void    SetFile(TString f);
  TString GetFile();

  void GetString(TString p, TString &v, TString dv="");
  void GetBool(TString p, bool &v, bool dv=false);
  void GetInt(TString p, int &v, int dv=0);
  void GetDouble(TString p, double &v, double dv=0.0);
  void GetVectorString(TString p, std::vector<TString> &v, TString dv="");
  void GetVectorStringDouble(TString p, std::vector<TString> &v1, std::vector<double> &v2);

 private:
  static TString file;
  
  template<typename T> void GetParameter(TString p, T &v, T dv);

};
#endif
