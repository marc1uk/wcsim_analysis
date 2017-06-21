{
gROOT->ProcessLine(".L /home/marc/anniegpvm/wcsim/root_work/analysiscaller.cxx++g");
//gSystem->Load("/home/marc/anniegpvm/wcsim/root_work/analysiscaller_cxx.s");
//gROOT->ProcessLine("WCSimAnalysis* theana = analysiscaller(\"/home/marc/anniegpvm/wcsim/root_work/in/MRD_muon_sample\")");
gROOT->ProcessLine("WCSimAnalysis* theana = analysiscaller(\"/home/marc/anniegpvm/wcsim/root_work/in/temp\")");
}

