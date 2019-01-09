.L /annie/app/users/moflaher/loadmaps.C+
TFile *_file0 = TFile::Open("/pnfs/annie/persistent/users/moflaher/wcsim/lappd/tankonly/wcsim_lappd_tankonly_24-09-17_BNB_Water_10k_22-05-17_truthana/EventDistributions.0.0.root");
TTree* t = treeout;
std::map<std::string,bool> thetypes;
std::map<string,bool>* thetypesp=&thetypes;
t->SetBranchAddress("TypesMap",&thetypesp);
std::vector<std::map<std::string,bool>> eventmaps; for(int i=0; i<t->GetEntries(); i++){t->GetEntry(i); eventmaps.push_back(thetypes);}
std::vector<std::string> typestrings;
for(auto anen : eventmaps){ std::string thisstring; if(anen.at("IsWeakCC")) thisstring="CC_"; else thisstring="NC_"; if(anen.at("IsQuasiElastic")) thisstring+="QE"; if(anen.at("IsResonant")) thisstring+="RES";  if(anen.at("IsDeepInelastic")) thisstring+="DIS";  if(anen.at("IsCoherent")) thisstring+="COH";  if(anen.at("IsSingleKaon")) thisstring+="KAON";  if(anen.at("IsNuElectronElastic")) thisstring+="e-_EL";  if(anen.at("IsMEC")) thisstring+="MEC"; typestrings.push_back(thisstring); }
std::map<string,int> sti;
std::map<int,string> ist;
for(auto as : typestrings){ if(sti.count(as)) continue; sti.emplace(as,sti.size()); ist.emplace(ist.size(),as); }
for(auto an : sti){ cout<<an.first<<"="<<an.second<<endl; }
TH1I htypes("htypes","Interaction Types",4,0,4);
for(auto anen : typestrings){ if(sti.count(anen)) htypes.Fill(sti.at(anen)); }
for(int i=0; i<sti.size(); i++){ htypes.GetXaxis()->SetBinLabel(i+1,ist.at(i).c_str()); }
std::map<int,int> bincountvsindex;
for(int bini=1; bini<htypes.GetNbinsX()+1; bini++){ bincountvsindex.emplace(htypes.GetBinContent(bini),bini); }
TH1I htypesorted("htypesorted","Interaction Types",4,0,4);
// use the fact that maps are inherently sorted by ascending key, which in this case is # entries
int bini=1;
for(auto abin=bincountvsindex.rbegin(); abin!=bincountvsindex.rend(); abin++){ htypesorted.SetBinContent(bini,abin->first); htypesorted.GetXaxis()->SetBinLabel(bini,htypes.GetXaxis()->GetBinLabel(abin->second)); bini++;}
htypesorted.Draw()
