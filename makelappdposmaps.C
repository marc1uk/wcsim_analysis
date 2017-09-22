{
cout<<"starting..."<<endl;
std::string wcsimfilepath="../build/wcsim_0.root";
TFile* wcsimfile = TFile::Open(wcsimfilepath.c_str());
TTree* geotree = (TTree*)wcsimfile->Get("wcsimGeoT");
if(geotree==0){ cerr<<"NO GEOMETRY IN FIRST FILE?"<<endl; assert(false); }
WCSimRootGeom* geo = 0;
geotree->SetBranchAddress("wcsimrootgeom", &geo);
if (geotree->GetEntries() == 0) { cerr<<"geotree has no entries!"<<endl; exit(9); }
geotree->GetEntry(0);
int numpmts = geo->GetWCNumPMT();
int numlappds = geo->GetWCNumLAPPD();
std::map<int,TVector3> lappdposmap;
std::map<int, double> lappdanglemap;
std::map<int, int> lappdlocmap;

bool getpmts=true; // get lappds otherwise

cout<<"beginning loop"<<endl;
for(int LAPPDID=1; LAPPDID<numlappds+1; LAPPDID++){
    cout<<"lAPPDID="<<LAPPDID<<endl;
    WCSimRootPMT pmt;
    if(getpmts) pmt = geo->GetLAPPD(LAPPDID-1);
    else        pmt = geo->GetLAPPD(LAPPDID-1);
    double pmtx = pmt.GetPosition(0);
    double pmty = pmt.GetPosition(1);
    double pmtz = pmt.GetPosition(2);
    lappdposmap.emplace(LAPPDID,TVector3(pmtx,pmty,pmtz));
    
    int thepmtsloc = pmt.GetCylLoc();
    lappdlocmap.emplace(LAPPDID,thepmtsloc);
    
    if(getpmts) continue;
    Float_t tank_start = 15.70;          // front face of the tank in cm
    Float_t tank_radius = 152.4;         // tank radius in cm
    double Rinnerstruct=270.9/2.;
    double Rthresh=Rinnerstruct*pow(2.,-0.5);
    double tileangle;
    int theoctagonside;
    switch (thepmtsloc){
        case 0: break; // top cap
        case 2: break; // bottom cap
        case 1: // wall
            // we need to account for the angle of the LAPPD within the tank
            // determine the angle based on it's position
            double octangle1=TMath::Pi()*(3./8.);
            double octangle2=TMath::Pi()*(1./8.);
            pmtz+=-tank_start-tank_radius;
                 if(pmtx<-Rthresh&&pmtz<0)         {tileangle=-octangle1; theoctagonside=0;}
            else if(-Rthresh<pmtx&&pmtx<0&&pmtz<0) {tileangle=-octangle2; theoctagonside=1;}
            else if(0<pmtx&&pmtx<Rthresh&&pmtz<0)  {tileangle=octangle2;  theoctagonside=2;}
            else if(Rthresh<pmtx&&pmtz<0)          {tileangle=octangle1;  theoctagonside=3;}
            else if(pmtx<-Rthresh&&pmtz>0)         {tileangle=octangle1;  theoctagonside=4;}
            else if(-Rthresh<pmtx&&pmtx<0&&pmtz>0) {tileangle=octangle2;  theoctagonside=5;}
            else if(0<pmtx&&pmtx<Rthresh&&pmtz>0)  {tileangle=-octangle2; theoctagonside=6;}
            else if(Rthresh<pmtx&&pmtz>0)          {tileangle=-octangle1; theoctagonside=7;}
            break;
    }
    lappdanglemap.emplace(LAPPDID, tileangle);
}

TFile* outfile;
TTree* outtree;
if(getpmts)  outfile = TFile::Open("pmtposfile.root","RECREATE");
else         outfile = TFile::Open("lappdposfile.root","RECREATE");
if(getpmts)  outtree = new TTree("pmttree","pmt positions");
else         outtree = new TTree("lappdtree","lappd positions");
double posx, posy, posz;
int loc;
double ang;
int id;
outtree->Branch("id",&id);
outtree->Branch("posx",&posx);
outtree->Branch("posy",&posy);
outtree->Branch("posz",&posz);
outtree->Branch("location",&loc);
if(!getpmts) outtree->Branch("angle",&ang);
cout<<"filling tree"<<endl;
cout<<"filling entry ";
for(int i=1; i<lappdposmap.size()+1; i++){
  cout<<i<<", ";
  posx=lappdposmap.at(i).X();
  posy=lappdposmap.at(i).Y();
  posz=lappdposmap.at(i).X();
  loc=lappdlocmap.at(i);
  if(!getpmts) ang=lappdanglemap.at(i);
  id=i;
  if(!getpmts) id+=numpmts;
  outtree->Fill();
}
cout<<endl<<"done"<<endl;
outfile->Write();
outfile->Close();

wcsimfile->Close();
}

