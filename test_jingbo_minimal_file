TFile* f = TFile::Open("/pnfs/annie/scratch/users/moflaher/trueQEvertexinfo_02_Mar_17.root");
TFile* f = TFile::Open("/annie/app/users/moflaher/wcsim/root_work/trueQEvertexinfo.root");
TTree* t = (TTree*)f->Get("vertextreenocuts");
ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > adigit;
std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >> eventdigits;
std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>* eventdigitsp = &eventdigits;
TBranch* digitsb;
t->SetBranchAddress("DigitVertices", &eventdigitsp, &digitsb);
TLorentzVector mustart;
TLorentzVector* mustartp=&mustart;
TBranch* mustartb;
t->SetBranchAddress("MuonStartVertex",&mustartp, &mustartb);
mustartb->GetEntry(0);
digitsb->GetEntry(0);
eventdigits.size();    // << Number of digits in this event
adigit=eventdigits.at(0);
cout<<"(x, y, z, t) = ("<<adigit.X()<<", "<<adigit.Y()<<", "<<adigit.Z()<<", "<<adigit.T()<<")"<<endl;
// returns (x,y,z,t)=(114.399261,-71.719872,107.832298,942.000000)
TH1D digitimeshist = TH1D("digitimeshist","Digi Times",100,0,10000);
TH1D mutimeshist = TH1D("mutimeshist","Mu Times",100,0,10000);
for(int i=0; i<digitsb->GetEntries(); i++){ digitsb->GetEntry(i); mustartb->GetEntry(i); for(int j=0; j<eventdigits.size(); j++){ digit=eventdigits.at(j); digitimeshist.Fill(digit.T()); } mutimeshist.Fill(mustart.T()); }
TCanvas c1("c1","A Canv", 600,600);
digitimeshist->Draw();
mutimeshist->SetLineColour(kRed);
mutimeshist->Draw("same");

