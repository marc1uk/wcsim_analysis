void truthcaller(const char* wcsimpathin="", const char* dirtpathin="", const char* geniepathin="", const char* outpathin=""){
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
	// disable some prolific compilation warnings
	// this _may_ not completely work as make command is a 2-part string: "make objects; make exe" so appending
	// the extra flags to the end of the string will only affect the latter part. Maybe use splits on ';' and append
	// to all...? or sed replace ';' with 'extraflags+;'?
//	std::string compileflags=" -fdiagnostics-color=always -Wno-sign-compare -Wno-overloaded-virtual -Wno-literal-conversion";
//	std::string makesharedlib = gSystem->GetMakeSharedLib();
//	makesharedlib+=compileflags;
//	gSystem->SetMakeSharedLib(makesharedlib.c_str());
//	std::string makeexe = gSystem->GetMakeExe();
//	makeexe+=compileflags;
//	gSystem->SetMakeExe(makeexe.c_str());
	
	std::string linetoprocess = ".L " + std::string(curr_dir.Data()) + "/plottruthmrdtracks.C++g";
	// seems you can combine strings with + only if one of them is a std::string.
	gROOT->ProcessLine(linetoprocess.c_str());
	//gROOT->ProcessLine(".L /annie/app/users/moflaher/wcsim/root_work/plottruthmrdtracks.C++g");
	// double + is required in case a compilation fails due to bad linking
	TString runprogram = TString::Format("truthtracks(\"%s\",\"%s\",\"%s\",\"%s\")", 
								wcsimpathin, dirtpathin, geniepathin, outpathin);
	gROOT->ProcessLine(runprogram);
}
