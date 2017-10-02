{
//TFile *_file0 = TFile::Open("/annie/data/users/wjingbo/reconstruction/reco25_5LAPPDs+120PMTs_ext2_0.root")
//TTree* tj=(TTree*)_file0->Get("vertextree")
TChain* tj = new TChain("vertextree");
tj->Add("/annie/data/users/wjingbo/reconstruction/reco25_*");
tj->SetBranchStatus("*",0);
tj->SetBranchStatus("TrueNeutrinoEnergy",1);
tj->SetBranchStatus("event",1);
int treenum=-1;
TFile* f = TFile::Open("/annie/app/users/moflaher/wcsim/root_work/temp/trueQEvertexinfo.root");
TTree* t = (TTree*)f->Get("vertextreenocuts");
int eventnumj; double eventnumm;
double truenuej, truenuem;
double mutracklengthinMRD;
t->SetBranchStatus("NeutrinoEnergy",1);
t->SetBranchStatus("EventNum",1);
t->SetBranchStatus("TrackLengthInMrd",1);
t->SetBranchAddress("NeutrinoEnergy",&truenuem);
t->SetBranchAddress("EventNum",&eventnumm);
t->SetBranchAddress("TrackLengthInMrd",&mutracklengthinMRD);
TFile* fout = new TFile("fixedmrtracklengths.root","RECREATE");
TTree* treeout = new TTree("fixtree","Corrected MRD Track Lengths");
double mutracklengthinmrdfixed;
treeout->Branch("TrackLengthInMrdFixed",&mutracklengthinmrdfixed);
int j=0;
for(int i=0; i<tj->GetEntries(); i++){
        int localentry = tj->LoadTree(i);
        int nextreenum=tj->GetTreeNumber();
        if(nextreenum!=treenum){
                tj->SetBranchAddress("TrueNeutrinoEnergy",&truenuej);
                tj->SetBranchAddress("event",&eventnumj);
                treenum=nextreenum;
        }
        tj->GetEntry(i);
        //cout<<"tj entry "<<i<<" has eventnumj "<<eventnumj<<endl;
        do{
                t->GetEntry(j);
                //cout<<"t entry "<<j<<" has eventnumm "<<eventnumm<<endl;
                j++;
        } while(eventnumm!=eventnumj);
         cout<<"j entry "<<i<<" had event num "<<eventnumj<<" and truenue="
             <<truenuej<<", m entry "<<j<<" had event num "<<eventnumm
             <<", truenuem="<<truenuem<<endl;
        mutracklengthinmrdfixed=mutracklengthinMRD;
        treeout->Fill();
}
fout->Write("",TObject::kOverwrite);
treeout->ResetBranchAddresses();
tj->ResetBranchAddresses();
t->ResetBranchAddresses();
fout->Close();
f->Close();
//_file0->Close();
}

