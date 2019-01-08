{
root out/q2comparison/30ps_proper_split/pmts_only/q2comparison.root out/q2comparison/30ps_proper_split/pmts_and_lappds/q2comparison.root
TH1D* hp = (TH1D*)_file0->Get("q2fracerrhist")
hp->SetLineColor(kRed)
hp->SetTitle("PMTs Only")
TH1D* hpl = (TH1D*)_file1->Get("q2fracerrhist")
hpl->SetLineColor(kBlue)
hpl->SetTitle("PMTs+LAPPDs")
hpl->DrawNormalized()
hp->DrawNormalized("same")
c1->BuildLegend()
Preliminary("Simulation");
}
