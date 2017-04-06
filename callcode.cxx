{
TString pwd = gSystem->Getenv("PWD");
TString compilestring = ".L " + pwd + "/analysiscaller.cxx++g";
//".L ~/LinuxSystemFiles/WCSim/gitver/root_work/analysiscaller.cxx++g"
gROOT->ProcessLine(compilestring);
gROOT->ProcessLine("WCSimAnalysis* theana = analysiscaller(\"/pnfs/annie/persistent/users/moflaher/wcsim\")");
gROOT->ProcessLine("delete theana");
}
