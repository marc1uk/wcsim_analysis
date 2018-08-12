std::vector<double> timeseps;
int currentlappdid=-1;
int nextlappdid=-1;
std::vector<ROOT::Math::XYZTVector>  filedigitvertices;
std::vector<ROOT::Math::XYZTVector>* filedigitverticesp = &filedigitvertices;
vertextreenocuts->SetBranchAddress("DigitVertices",&filedigitverticesp);
std::vector<std::string> filedigitsensortypes;
std::vector<std::string>* filedigitsensortypesp=&filedigitsensortypes;
vertextreenocuts->SetBranchAddress("DigitWhichDet",&filedigitsensortypesp);
std::vector<int> filedigitPMTIDs;
std::vector<int>* filedigitPMTIDsp=&filedigitPMTIDs;
vertextreenocuts->SetBranchAddress("DigitPmtId",&filedigitPMTIDsp);
double hitone;
double hittwo;

for(int entryi=0; entryi<vertextreenocuts->GetEntries(); entryi++){
  vertextreenocuts->GetEntry(entryi);
  for(int hiti=0; hiti<filedigitvertices.size(); hiti++){
    if(filedigitsensortypes.at(hiti)!="lappd_v0"){
      currentlappdid=-1;
      continue;
    }
    ROOT::Math::XYZTVector adigitvector = filedigitvertices.at(hiti);
    double digitt = adigitvector.T();
    nextlappdid=filedigitPMTIDs.at(hiti);
    if(nextlappdid!=currentlappdid){
      hittwo=digitt;
      currentlappdid=nextlappdid;
    } else {
      hittone=hittwo;
      hittwo=digitt;
      double hittdiff = abs(hittwo-hittone);
      timeseps.push_back(hittdiff);
      //cout<<"hitone="<<hitone<<", hittwo="<<hittwo<<"hittdiff="<<hittdiff<<endl;
      //if(timeseps.size()>10) break;
    }
  }
  currentlappdid=-1;
  //if(timeseps.size()>10) break;
}

double hitsepmax = *std::max(timeseps.begin(),timeseps.end());
TH1D* hittimeseps = new TH1D("hittimeseps","Time Separation Between Hits on LAPPDs",200,0,hitsepmax);
for(auto ats : timeseps) hittimeseps->Fill(ats);
hittimeseps->Draw();

TH1D* hittimeseps2 = new TH1D("hittimeseps2","Time Separation Between Hits on LAPPDs",200,0,20);
for(auto ats : timeseps) hittimeseps2->Fill(ats);
hittimeseps2->Draw();

TH1D* hittimeseps3 = new TH1D("hittimeseps3","Time Separation Between Hits on LAPPDs",200,0,1);
for(auto ats : timeseps) hittimeseps3->Fill(ats);
hittimeseps3->Draw();
