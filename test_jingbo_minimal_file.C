{
gROOT->ProcessLine("#include <thread>");
gROOT->ProcessLine("#include <chrono>");
gROOT->ProcessLine("#include <time.h>");
//TFile* f = TFile::Open("/pnfs/annie/scratch/users/moflaher/trueQEvertexinfo_02_Mar_17.root");
//TFile* f = TFile::Open("/annie/app/users/moflaher/wcsim/root_work/trueQEvertexinfo.root");
TFile* f = TFile::Open("~/anniegpvm/wcsim/root_work/temp/trueQEvertexinfo.root");
//TFile* f = TFile::Open("root://moflaher@anniegpvm02.fnal.gov//annie/app/users/moflaher/wcsim/root_work/temp/trueQEvertexinfo.root");

//TFile* f2 = TFile::Open("~/anniegpvm/wcsim/build/wcsim_lappd_0.1000.root");
//int lappd_numhitsthisevt;
//std::vector<double> lappdhit_globalcoorx;
//std::vector<double> lappdhit_globalcoory;
//std::vector<double> lappdhit_globalcoorz;
//std::vector<double>* lappdhit_globalcoorxp=&lappdhit_globalcoorx;
//std::vector<double>* lappdhit_globalcooryp=&lappdhit_globalcoory;
//std::vector<double>* lappdhit_globalcoorzp=&lappdhit_globalcoorz;
//
//TTree* lappdtree = (TTree*)f2->Get("LAPPDTree");
//branchesok =lappdtree->SetBranchAddress("lappd_numhits",          &lappd_numhitsthisevt);
//if(branchesok<0) cerr<<"lappd_numhits="<<branchesok<<endl;
//branchesok =lappdtree->SetBranchAddress("lappdhit_globalcoorx",    &lappdhit_globalcoorxp);
//if(branchesok<0) cerr<<"lappd_globalcoorx="<<branchesok<<endl;
//branchesok =lappdtree->SetBranchAddress("lappdhit_globalcoory",    &lappdhit_globalcooryp);
//if(branchesok<0) cerr<<"lappd_globalcoory="<<branchesok<<endl;
//branchesok =lappdtree->SetBranchAddress("lappdhit_globalcoorz",    &lappdhit_globalcoorzp);
//if(branchesok<0) cerr<<"lappd_globalcoorz="<<branchesok<<endl;

cout<<"f="<<f<<endl;
if(!f){ cerr<<"no file"<<endl; exit;}
TTree* t = (TTree*)f->Get("vertextreenocuts");
cout<<"t="<<t<<endl;
if(!t){ cerr<<"no tree"<<endl; exit;}

TLorentzVector mustart;
TLorentzVector* mustartp=&mustart;
TLorentzVector mustop;
TLorentzVector* mustopp=&mustop;
TVector3 mudir;
TVector3* mudirp=&mudir;
ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > adigit;
std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >> eventdigits;
std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>* eventdigitsp = &eventdigits;
std::vector<double> digitqs;
std::vector<double>* digitqsp=&digitqs;
std::vector<double> digittsmears;
std::vector<double>* digittsmearsp=&digittsmears;
std::vector<string> digitpstypes;
std::vector<string>* digitpstypesp=&digitpstypes;

std::vector<double> intileposx;
std::vector<double> intileposy;
std::vector<double> poserrx;
std::vector<double> poserry;
std::vector<double> poserrz;
std::vector<int> tileorient;
std::vector<double>* intileposxp=&intileposx;
std::vector<double>* intileposyp=&intileposy;
std::vector<double>* poserrxp=&poserrx;
std::vector<double>* poserryp=&poserry;
std::vector<double>* poserrzp=&poserrz;
std::vector<int>* tileorientp=&tileorient;

TBranch* digitsb;
TBranch* mustartb;
TBranch* mustopb;
TBranch* mudirb;
TBranch* digitqsb;
TBranch* digittsmearsb;
TBranch* digitpstypesb;

cout<<"setting branch addresses"<<endl;
t->SetBranchAddress("MuonStartVertex",&mustartp, &mustartb);
t->SetBranchAddress("MuonStopVertex",&mustopp, &mustopb);
t->SetBranchAddress("MuonDirection",&mudirp, &mudirb);
t->SetBranchAddress("DigitVertices", &eventdigitsp, &digitsb);
t->SetBranchAddress("DigitCharges",&digitqsp, &digitqsb);
t->SetBranchAddress("DigitTimeSmear",&digittsmearsp, &digittsmearsb);
t->SetBranchAddress("DigitWhichDet",&digitpstypesp, &digitpstypesb);

cout<<" mustartb="<<mustartb<<endl;
cout<<" mustopb="<<mustopb<<endl;
cout<<" mudirb="<<mudirb<<endl;
cout<<" digitsb="<<digitsb<<endl;
cout<<" digitqsb="<<digitqsb<<endl;
cout<<" digittsmearsb="<<digittsmearsb<<endl;
cout<<" digitpstypesb="<<digitpstypesb<<endl;
if(mustartb==0||mustopb==0||mudirb==0||digitsb==0||digitqsb==0||digittsmearsb==0||digitpstypesb==0) exit;

for(int eventi=0; eventi<5; eventi++){
  t->GetEntry(eventi);
  cout<<"event "<<eventi<<" has "<<eventdigits.size()<<" digits"<<endl;
  cout<<"vector of digit charges has "<<digitqs.size()<<" elements"<<endl;
  cout<<"vector of digit time smearings has "<<digittsmears.size()<<" elements"<<endl;
  cout<<"vector of digit pmt types has "<<digitpstypes.size()<<" elements"<<endl;;
  for(int j=0; j<5; j++){ // loop over 5 digits in this event
    adigit=eventdigits.at(j);
    cout<<"(x, y, z, t) = ("<<adigit.X()<<", "<<adigit.Y()<<", "<<adigit.Z()<<", "<<adigit.T()<<")"<<endl;
    cout<<"sensortype is "<<digitpstypes.at(j)<<endl;
    cout<<"hit time smearing is "<<digittsmears.at(j)<<endl;
    cout<<"hit charge is "<<digitqs.at(j)<<endl;
    // returns (x,y,z,t)=(114.399261,-71.719872,107.832298,942.000000)
  }
}
TH1D digiposxpmthist = TH1D("digiposxpmthist","PMT Digi Pos X",100,-200,200);
TH1D digiposxlappdhist = TH1D("digiposxlappdhist","PMT Digi Pos X",100,-200,200);
TH1D digiposypmthist = TH1D("digiposypmthist","PMT Digi Pos Y",100,-220,220);
TH1D digiposylappdhist = TH1D("digiposylappdhist","LAPPD Digi Pos Y",100,-220,220);
TH1D digiposzpmthist = TH1D("digiposzpmthist","LAPPD Digi Pos Z",100,0,320);
TH1D digiposzlappdhist = TH1D("digiposzlappdhist","LAPPD Digi Pos Z",100,0,320);

TH2D digitposxvszpmthist = TH2D("digitposxvszpmthist","PMT Digi Pos X-Z",100,-200,200,100,0,320);
TH2D digitposyvszpmthist = TH2D("digitposyvszpmthist","PMT Digi Pos Y-Z",100,-220,220,100,0,320);
//TH2D digitposxvszlappdhist = TH2D("digitposxvszlappdhist","LAPPD Digi Pos X-Z",100,-200,200,100,0,320);
//TH2D digitposyvszlappdhist = TH2D("digitposyvszlappdhist","LAPPD Digi Pos Y-Z",100,-220,220,100,0,320);

// to get an unbinned histogram we need to get a TPolyMarker3D (the points) and a TH3F (the frame, axes, title)
// we can get these from a tree draw... 
//TH3D digitpositionlappdhist = TH3D("digitpositionlappdhist","LAPPD Digi Pos",50,-200,200,50,0,320,50,-220,220);

TH1D digitimespmthist = TH1D("digitimespmthist","PMT Digi Times",1000,0,10000);
TH1D digitimeslappdhist = TH1D("digitimeslappdhist","LAPPD Digi Times",1000,0,10000);
TH1D mutimeshist = TH1D("mutieshist","Mu Times",100,0,10000);

TH1D pmttsmearingh = TH1D("pmttsmearingh","PMT T Smearing",100,0,3); // 0-3ns smearing
TH1D lappdtsmearingh = TH1D("lappdtsmearingh1","LAPPD T Smearing",100,0,0.1); // 0-100ps smearing
TH1D digitqspmth = TH1D("digitqspmth","Digit Charges PMT",1000,0,30);
TH1D digitqslappdh = TH1D("digitqslappdh","Digit Charges LAPPD",1000,0,30);
int pmtcount=0, lappdcount=0;
cout<<"t has "<<t->GetEntries()<<" entries"<<endl;
for(int i=0; i<std::min((double)t->GetEntries(),100.); i++){
  cout<<"getting entry "<<i<<endl;
  t->GetEntry(i);
  //lappdtree->GetEntry(i);
  
  //cout<<eventdigits.size()<<" digits in entry "<<i<<endl;
  for(int j=0; j<eventdigits.size(); j++){
    //cout<<"getting digit "<<j<<endl;
    adigit=eventdigits.at(j);
    
    //cout<<"getting pmttype"<<endl;
    std::string sensortypestring=digitpstypes.at(j);
    //cout<<"its "<<sensortypestring<<endl;
    if(sensortypestring.find("lappd")==std::string::npos){
      //cout<<"filling pmttsmearingh"<<endl;
      pmttsmearingh.Fill(digittsmears.at(j));
      //cout<<"filling charge"<<endl;
      digitqspmth.Fill(digitqs.at(j));
      digiposxpmthist.Fill(adigit.X());
      digiposypmthist.Fill(adigit.Y());
      digiposzpmthist.Fill(adigit.Z());
      digitposxvszpmthist.Fill(adigit.X(),adigit.Z());
      digitposyvszpmthist.Fill(adigit.Y(),adigit.Z());
      //cout<<"filling time"<<endl;
      digitimespmthist.Fill(adigit.T());
      pmtcount++;
    } else {
      //cout<<"filling lappdtsmearingh"<<endl;
      lappdtsmearingh.Fill(digittsmears.at(j));
      //cout<<"filling charge"<<endl;
      digitqslappdh.Fill(digitqs.at(j));
      digiposxlappdhist.Fill(adigit.X());
      digiposylappdhist.Fill(adigit.Y());
      digiposzlappdhist.Fill(adigit.Z());
      
//      digitposxvszlappdhist.Fill(adigit.X(),adigit.Z());
//      digitposyvszlappdhist.Fill(adigit.Y(),adigit.Z());
      //if(abs((adigit.Y())+14.46)<190)
      //digitpositionlappdhist.Fill(adigit.X(),adigit.Z(),adigit.Y());
      //cout<<"filling time"<<endl;
      digitimeslappdhist.Fill(adigit.T());
      lappdcount++;
    }
  }
  //cout<<"filling mutimeshist"<<endl;
  mutimeshist.Fill(mustart.T());
  //std::this_thread::sleep_for (std::chrono::seconds(15));
}
cout<<"over first 1000 hits saw "<<pmtcount<<" hits on PMTs and "<<lappdcount<<" hits on LAPPDs"<<endl;
TCanvas c1("c1","A Canv", 600,600);
digitimeslappdhist.SetLineColor(kRed);
digitimeslappdhist.Draw();
digitimespmthist.Draw("same");
TCanvas c2;
lappdtsmearingh.SetLineColor(kRed);
lappdtsmearingh.Draw();
TCanvas c3;
pmttsmearingh.Draw();
TCanvas c4;
digitqslappdh.SetLineColor(kRed);
digitqslappdh.Draw();
digitqspmth.Draw("same"); // more entries in lappd histo so draw this after
TCanvas c5;
digiposxlappdhist.Draw();
//digitposxvszpmthist.SetMarkerStyle(2);
//digitposxvszpmthist.Draw();
//digitposxvszlappdhist.SetMarkerColor(kRed);
//digitposxvszlappdhist.SetMarkerStyle(6);
//digitposxvszlappdhist.Draw("same");
TCanvas c6;
digiposylappdhist.Draw();
//digitposyvszpmthist.Draw();
//digitposyvszpmthist.SetMarkerStyle(2);
//digitposyvszlappdhist.SetMarkerColor(kRed);
//digitposyvszlappdhist.SetMarkerStyle(6);
//digitposyvszlappdhist.Draw("same");
TCanvas c7;
digiposzlappdhist.Draw();
TCanvas c8;
// Draw("Zval","Yval","Xval") << this is how the AXES are FILLED
// so hist->GetXaxis() will retrieve the axis on which the LAST value are plotted!
// the Z axis is plotted VERTICALLY. the others are effectively interchangable
t->Draw("DigitVertices.Y()+140.:DigitVertices.Z():(DigitVertices.X())>>htemp","DigitWhichDet==\"lappd_v0\"");
TH3F digitpositionlappdhist(*((TH3F*)c8.GetPrimitive("htemp")));
TPolyMarker3D digitpositionlappdhistpm(*((TPolyMarker3D*)c8.GetListOfPrimitives()->At(2)));
digitpositionlappdhist.Reset(); // clear the histogram binned points
digitpositionlappdhist.GetXaxis()->SetTitle("X");
digitpositionlappdhist.GetYaxis()->SetTitle("Z"); // << Z and Y are swapped because ANNIE uses z as beam
digitpositionlappdhist.GetZaxis()->SetTitle("Y");
digitpositionlappdhist.GetXaxis()->CenterTitle();
digitpositionlappdhist.GetYaxis()->CenterTitle();
digitpositionlappdhist.GetZaxis()->CenterTitle();
digitpositionlappdhist.Draw();
//digitpositionlappdhistpm.SetMarkerStyle(6);
digitpositionlappdhistpm.Draw("same");

}
