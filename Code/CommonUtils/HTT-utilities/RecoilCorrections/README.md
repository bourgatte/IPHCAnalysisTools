## MET recoil corrections and uncertainties

### How to set up the MET Recoil Correction Interface

```bash
cd ${CMSSW_BASE}/src
cmsenv
git clone https://github.com/KIT-CMS/RecoilCorrections.git HTT-utilities/RecoilCorrections 
scram b
```

### How to apply MET recoil corrections (supported only in C++)

```cpp
// add the header file to your source file
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"

...

// create instances of class RecoilCorrection and
// load MET recoil response and resolution histograms before looping over events
// the path to files is defined relative to ${CMSSW_BASE}/src directory

// use these ROOT files when correcting Type I PF and Puppi MET
RecoilCorrector recoilPFMetCorrector2016("HTT-utilities/RecoilCorrections/data/Type1_PFMET_2016.root"); // Type I PF MET 2016
RecoilCorrector recoilPFMetCorrector2017("HTT-utilities/RecoilCorrections/data/Type1_PFMET_2017.root"); // Type I PF MET 2017
RecoilCorrector recoilPFMetCorrector2018("HTT-utilities/RecoilCorrections/data/Type1_PFMET_2018.root"); // Type I PF MET 2018

RecoilCorrector recoilPuppiMetCorrector2016("HTT-utilities/RecoilCorrections/data/Type1_PuppiMET_2016.root"); // Type I Puppi MET 2016
RecoilCorrector recoilPuppiMetCorrector2017("HTT-utilities/RecoilCorrections/data/Type1_PuppiMET_2017.root"); // Type I Puppi MET 2017
RecoilCorrector recoilPuppiMetCorrector2018("HTT-utilities/RecoilCorrections/data/Type1_PuppiMET_2018.root"); // Type I Puppi MET 2018

...

// apply recoil corrections on event-by-event basis (Type I PF MET 2016 example with quantile mapping on histograms)
recoilPFMetCorrector2016.CorrectWithHist(
    pfmet_ex, // uncorrected type I pf met px (float)
    pfmet_ey, // uncorrected type I pf met py (float)
    genPx, // generator Z/W/Higgs px (float)
    genPy, // generator Z/W/Higgs py (float)
    visPx, // generator visible Z/W/Higgs px (float)
    visPy, // generator visible Z/W/Higgs py (float)
    njets,  // number of jets (hadronic jet multiplicity) (int)
    pfmetcorr_ex, // corrected type I pf met px (float)
    pfmetcorr_ey  // corrected type I pf met py (float)
);
```

### How to perform MET recoil uncertainty shifts (supported only in C++)

```cpp
// add the header file to your source file
#include "HTT-utilities/RecoilCorrections/interface/MEtSys.h"

...

// create instances of class MEtSys and
// load recoil response histograms before looping over events
// the path to files is defined relative to ${CMSSW_BASE}/src directory

// use this ROOT file when shifting Type I PF and Puppi MET
RecoilCorrector recoilPFMetShifter2016("HTT-utilities/RecoilCorrections/data/PFMETSys_2016.root"); // Type I PF MET 2016
RecoilCorrector recoilPFMetShifter2017("HTT-utilities/RecoilCorrections/data/PFMETSys_2017.root"); // Type I PF MET 2017
RecoilCorrector recoilPFMetShifter2018("HTT-utilities/RecoilCorrections/data/PFMETSys_2018.root"); // Type I PF MET 2018

RecoilCorrector recoilPuppiMetShifter2016("HTT-utilities/RecoilCorrections/data/PuppiMETSys_2016.root"); // Type I Puppi MET 2016
RecoilCorrector recoilPuppiMetShifter2017("HTT-utilities/RecoilCorrections/data/PuppiMETSys_2017.root"); // Type I Puppi MET 2017
RecoilCorrector recoilPuppiMetShifter2018("HTT-utilities/RecoilCorrections/data/PuppiMETSys_2018.root"); // Type I Puppi MET 2018

...

// apply recoil corrections on event-by-event basis (Type I Puppi MET 2018 example, upward shift for response)
recoilPuppiMetShifter2018.ApplyMEtSys(
    puppimet_ex, // uncorrected type I puppi met px (float)
    puppimet_ey, // uncorrected type I puppi met py (float)
    genPx, // generator Z/W/Higgs px (float)
    genPy, // generator Z/W/Higgs py (float)
    visPx, // generator visible Z/W/Higgs px (float)
    visPy, // generator visible Z/W/Higgs py (float)
    njets,  // number of jets (hadronic jet multiplicity) (int)
    MEtSys::SysType::Response, // shift for hadronic recoil response
    MEtSys::SysShift::Up, // upward shift
    puppimetshifted_ex, // shifted type I puppi met px (float)
    puppimetshifted_ey  // shifted type I puppi met py (float)
);
```
