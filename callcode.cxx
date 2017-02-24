{gROOT->ProcessLine(".L ~/LinuxSystemFiles/WCSim/gitver/root_work/analysiscaller.cxx++g");
gROOT->ProcessLine("WCSimAnalysis* theana = analysiscaller()");
gROOT->ProcessLine("delete theana");
}
