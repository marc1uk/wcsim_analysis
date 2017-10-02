{
/* vim:set noexpandtab tabstop=4 wrap */
#ifndef PARTICLEGUNEVENTS
//#define PARTICLEGUNEVENTS
//#define NOGENIE
#endif
	TString curr_dir = gSystem->pwd();
#ifndef NOGENIE
	TString script_dir = gSystem->Getenv("GENIE");
	script_dir += "/src/scripts/gcint/";
	gSystem->cd(script_dir.Data());
	gROOT->ProcessLine(".x loadincs.C");
	gROOT->ProcessLine(".x loadlibs.C");
	gSystem->cd(curr_dir.Data());
#endif
	std::string linetoprocess = ".L " + std::string(curr_dir.Data()) + "/plottruthmrdtracks.C++g";
	// seems you can combine strings with + only if one of them is a std::string.
	gROOT->ProcessLine(linetoprocess.c_str());
	//gROOT->ProcessLine(".L /annie/app/users/moflaher/wcsim/root_work/plottruthmrdtracks.C++g");
	// double + is required in case a compilation fails due to bad linking
	truthtracks()
}
