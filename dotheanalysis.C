{
//TString compilerflags = "-Wno-sign-compare -Wno-overloaded-virtual -Wno-literal-conversion";
//gSystem->SetFlagsOpt(compilerflags.Data());
//gSystem->SetFlagsDebug(compilerflags.Data());
TString pwd = gSystem->Getenv("PWD");

TString wcsimlibraryloc = pwd + "/../wcsim/";
TString wcsimincludeloc = wcsimlibraryloc + "/include/";
gInterpreter->AddIncludePath(wcsimincludeloc);
// n.b. for later versions of root gInterpreter->AddInclude(path); will accept path prefixed with -I for consistency
//gSystem->AddIncludePath(TString::Format("%s%s","-I",wcsimincludeloc.Data()));

//TString includeline = "#include \"" + wcsimincludeloc;
//gInterpreter->Declare(includeline + "WCSimRootEvent.hh\"");
//gInterpreter->Declare(includeline + "WCSimRootGeom.hh\"");
//gInterpreter->Declare(includeline + "WCSimPmtInfo.hh\"");
//gInterpreter->Declare(includeline + "WCSimLAPPDInfo.hh\"");
//gInterpreter->Declare(includeline + "WCSimEnumerations.hh\"");
//gInterpreter->Declare(includeline + "WCSimRootOptions.hh\"");
//gInterpreter->Declare(includeline + "WCSimRootLinkDef.hh\"");

TString wcsimlibrarypath = wcsimlibraryloc + "libWCSimRoot.so";
gSystem->Load(wcsimlibrarypath);
gSystem->Load("libWCSimRoot");

// we can either load the pre-compiled library
//TString analysislibrarypath = pwd+"/analysiscaller_cxx.so";
//gSystem->Load(analysislibrarypath.Data());

// or we can compile it again
TString callcommand = ".L "+pwd+"/analysiscaller.cxx++g";
gROOT->ProcessLine(callcommand.Data());

TString analyzecommand = "WCSimAnalysis* theana = analysiscaller";

// various paths to different analysis directories
//TString filestoanalyze = pwd+"/in/MRD_muon_sample";
//TString filestoanalyze = pwd+"/in/temp";
TString filestoanalyze = "/pnfs/annie/persistent/users/moflaher/wcsim_tankonly_03-05-17_BNB_World_10k_29-06-17/wcsim_0.4000.root";
//TString filestoanalyze="/home/marc/LinuxSystemFiles/WCSim/gitver/root_work/in/temp";
//TString filestoanalyze = "/pnfs/annie/persistent/users/moflaher/wcsim_tankonly_17-06-17";
//TString outputdir = "/home/marc/LinuxSystemFiles/WCSim/gitver/root_work/in/temp";
TString outputdir = pwd+"/out/temp";
//TString outputdir = "/pnfs/annie/persistent/users/moflaher/wcsim_tankonly_17-06-17_ana";
// XXX: to limit files added to the chain, edit the pattern in utilityfuncs.cxx line 64 marked by the XXX

// call analysiscaller() method with the file path. 
// This constructs the WCSimAnalysis and calls DoAnalysis() on it. 
TString fullcall = analyzecommand+"(\""+filestoanalyze+"\", \""+outputdir+"\")";
gROOT->ProcessLine(fullcall.Data());
}

