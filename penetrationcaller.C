{
/* vim:set noexpandtab tabstop=4 wrap */
	TString script_dir = gSystem->Getenv("GENIE");
	script_dir += "/src/scripts/gcint/";
	TString curr_dir = gSystem->pwd();
	gSystem->cd(script_dir.Data());
	gROOT->ProcessLine(".x loadincs.C");
	gROOT->ProcessLine(".x loadlibs.C");
	gSystem->cd(curr_dir.Data());
	gROOT->ProcessLine(".L /annie/app/users/moflaher/wcsim/root_work/plot_mrd_penetration2.C++g");
	// double + is required in case a compilation fails due to bad linking
	truthtracks()
}
