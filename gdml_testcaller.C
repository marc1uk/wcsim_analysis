{
	/* vim:set noexpandtab tabstop=4 wrap */
	TString curr_dir = gSystem->pwd();
	gSystem->cd(curr_dir.Data());
	TString macropath = curr_dir + "/gdml_test.cxx";
	TString commandstring = ".L "+macropath+"++g";
	gROOT->ProcessLine("#include \"TEveLine.h\"");
	gROOT->ProcessLine(commandstring.Data());
	// double + is required in case a compilation fails due to bad linking
	gdml_test();
}
