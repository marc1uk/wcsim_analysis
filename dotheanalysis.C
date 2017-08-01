{
TString compilerflags = "-Wno-sign-compare -Wno-overloaded-virtual -Wno-literal-conversion";
gSystem->SetFlagsOpt(compilerflags.Data());
gSystem->SetFlagsDebug(compilerflags.Data());
TString pwd = gSystem->Getenv("PWD");

// we can either load the pre-compiled library
//TString analysislibrarypath = pwd+"/analysiscaller_cxx.s";
//gSystem->Load(analysislibrarypath.Data());

// or we can compile it again
TString callcommand = ".L "+pwd+"/analysiscaller.cxx++g";
gROOT->ProcessLine(callcommand.Data());

TString analyzecommand = "WCSimAnalysis* theana = analysiscaller";

// various paths to different analysis directories
//TString filestoanalyze = pwd+"/in/MRD_muon_sample";
//TString filestoanalyze = pwd+"/in/temp";
TString filestoanalyze="/home/marc/LinuxSystemFiles/WCSim/gitver/root_work/in/temp";
//TString filestoanalyze = "/pnfs/annie/persistent/users/moflaher/wcsim_tankonly_17-06-17";
TString outputdir = "/home/marc/LinuxSystemFiles/WCSim/gitver/root_work/in/temp";
//TString outputdir = "/pnfs/annie/persistent/users/moflaher/wcsim_tankonly_17-06-17_ana";
//TString outputdir = pwd;
// XXX: to limit files added to the chain, edit the pattern in utilityfuncs.cxx line 64 marked by the XXX

// call analysiscaller() method with the file path. 
// This constructs the WCSimAnalysis and calls DoAnalysis() on it. 
TString fullcall = analyzecommand+"(\""+filestoanalyze+"\", \""+outputdir+"\")";
gROOT->ProcessLine(fullcall.Data());
}

