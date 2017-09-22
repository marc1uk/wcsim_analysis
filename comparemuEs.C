/* vim:set noexpandtab tabstop=4 wrap */
#include "TROOT.h"
#include "TSystem.h"
#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include "TApplication.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TNtuple.h"
#include <string>
#include <iostream>

void comparemuEs(){
    // muon:
    TFile *f1 = TFile::Open("/annie/data/users/edrakopo/OUTrecomuonE.root");
    TTree* t1 = (TTree*)f1->Get("tuple");
    t1->SetBranchStatus("*",0);
    t1->SetBranchStatus("trueKE",1);
    t1->SetBranchStatus("recoKE",1);
    t1->SetBranchStatus("dirX",1);
    t1->SetBranchStatus("dirY",1);
    t1->SetBranchStatus("dirZ",1);
    TVector3 recomuondirection;
    double truemuE, recomuE, recomudirx, recomudiry, recomudirz;
    t1->SetBranchAddress("trueKE",&truemuE);
    t1->SetBranchAddress("recoKE",&recomuE);
    t1->SetBranchAddress("dirX",&recomudirx);
    t1->SetBranchAddress("dirY",&recomudiry);
    t1->SetBranchAddress("dirZ",&recomudirz);
    
    // neutrino:
    TFile *f2 = TFile::Open("/annie/data/users/edrakopo/OUTreconeutrinoE.root");
    TTree* t2 = (TTree*)f2->Get("tuple");
    t2->SetBranchStatus("*",0);
    t2->SetBranchStatus("trueKE",1);
    t2->SetBranchStatus("recoKE",1);
    t2->SetBranchStatus("vtxX",1);
    t2->SetBranchStatus("vtxY",1);
    t2->SetBranchStatus("vtxZ",1);
    double truenuE, reconuE, recovtxx, recovtxy, recovtxz;
    t2->SetBranchAddress("trueKE",&truenuE);
    t2->SetBranchAddress("recoKE",&reconuE);
    t2->SetBranchAddress("vtxX",&recovtxx);
    t2->SetBranchAddress("vtxY",&recovtxy);
    t2->SetBranchAddress("vtxZ",&recovtxz);
    
    // pull the true Q2 information for comparison. 
    // We can only compare distributions for now, not event-by-event, as event # or true Q2 were not saved
    TFile* f3 = TFile::Open("/annie/data/users/moflaher/trueQEvertexinfo_extv2.root");
    TTree* t3 = (TTree*)f3->Get("vertextreenocuts");
    t3->SetBranchStatus("*",0);
    t3->SetBranchStatus("NeutrinoEnergy",1);
    t3->SetBranchStatus("MuonEnergy",1);
    t3->SetBranchAddress("MuonAngle",1);
    t3->SetBranchAddress("MomentumTransfer",1);
    t3->SetBranchStatus("TrackLengthInMrd",1);
    t3->SetBranchStatus("MuonStartVertex",1);
    double truenuEall, truemuEall, trueangleall, trueQ2all, trackleninmrd;
    TVector3 truevertexall;
    t3->SetBranchAddress("NeutrinoEnergy",&truenuEall);
    t3->SetBranchAddress("MuonEnergy",&truemuEall);
    t3->SetBranchAddress("MuonAngle",&trueangleall);
    t3->SetBranchAddress("MomentumTransfer",&trueQ2all);
    t3->SetBranchAddress("TrackLengthInMrd",&trackleninmrd);
    t3->SetBranchAddress("MuonStartVertex",&truevertexall);
    
    // output file
    TFile* fileout = new TFile("q2comparison.root","RECREATE");
    TTree* treeout = new TTree("q2tree","Reco vs True Comparison of Q2");
    double outtruenuE, outtruemuE, outreconuE, outrecomuE, outtrueangle, outrecoangle, outtrueq2, outrecoq2, 
    outvtxerr, outangerr, outmuEerr, outnuEerr, outq2err;
    bool outinfidvol, outhasmrdtrack;
    TVector3 outtruevtx, outrecovtx, outtruedir, outrecodir;
    
    treeout->Branch("TrueNuE",&outtruenuE);
    treeout->Branch("RecoNuE",&outreconuE);
    treeout->Branch("TrueMuE",&outtruemuE);
    treeout->Branch("RecoMuE",&outrecomuE);
    treeout->Branch("TrueAngle",&outtrueangle);
    treeout->Branch("RecoAngle",&outrecoangle);
    treeout->Branch("TrueQ2",&outtrueq2);
    treeout->Branch("RecoQ2",&outrecoq2);
    treeout->Branch("VtxErr",&outvtxerr);
    treeout->Branch("AngErr",&outangerr);
    treeout->Branch("MuEerr",&outmuEerr);
    treeout->Branch("NuEerr",&outnuEerr);
    treeout->Branch("Q2Err",&outq2err);
    treeout->Branch("InFidVol",&outinfidvol);
    treeout->Branch("HasMrdTrack",&outhasmrdtrack);
    treeout->Branch("TrueVtx",&outtruevtx);
    treeout->Branch("RecoVtx",&outrecovtx);
    treeout->Branch("TrueDir",&outtruedir);
    treeout->Branch("RecoDir",&outrecodir);
    
    // histograms to compare distributions
    /* true sample histograms */
    TH1D truenuKEhist("truenuKEhist","True Neutrino KE",100,0,2000);
    t2->Draw("trueKE>>truenuKEhist");
    TH1D truemuKEhist("truemuKEhist","True Muon KE",100,0,2000);
    t1->Draw("trueKE>>truemuKEhist");
    TH1D trueanglehist("trueanglehist","True Scattering Angle",100,0,TMath::Pi());
    // do not yet have
    TH1D trueQ2hist("trueQ2hist","True Momentum Transfer",100,0,2000);
    // do not yet have
    TH3D truevertexhist("truevertexhist","True Vertex",100,-150,150,100,-200,200,100,50,350);
    // do not yet have
    TH3D truedirhist("truedirhist","True Direction",100,-1,1,100,-1,1,100,-1,1);
    // do not yet have
    
    /* reco sample histograms */
    TH1D reconuKEhist("reconuKEhist","Reco Neutrino KE",100,0,2000);
    t2->Draw("recoKE>>reconuKEhist");
    TH1D recomuKEhist("recomuKEhist","Reco Muon KE",100,0,2000);
    t1->Draw("recoKE>>recomuKEhist");
    TH1D recoanglehist("recoanglehist","Reco Scattering Angle",100,0,TMath::Pi());
    // fill in loop
    TH1D recoQ2hist("recoQ2hist","Reco Momentum Transfer",100,0,2000);
    // fill in loop
    TH3D recovertexhist("recovertexhist","Reco Vertex",100,-150,150,100,-200,200,100,50,350);
    t2->Draw("vtxX:vtxY:vtxZ>>recovertexhist");
    TH3D recodirhist("recodirhist","Reco Direction",100,-1,1,100,-1,1,100,-1,1);
    t2->Draw("dirX:dirY:dirZ>>recodirhist");
    
    /* true all histograms */
    TH1D truenuKEallhist("truenuKEallhist","All True Neutrino KE",100,0,2000);
    t3->Draw("NeutrinoEnergy>>truenuKEallhist");
    TH1D truemuKEallhist("truemuKEallhist","All True Muon KE",100,0,2000);
    t3->Draw("MuonEnergy>>truemuKEallhist");
    TH1D trueangleallhist("trueangleallhist","True Scattering Angle",100,0,TMath::Pi());
    t3->Draw("MuonAngle>>trueangleallhist");
    TH1D trueQ2allhist("trueQ2allhist","True Momentum Transfer",100,0,2000);
    t3->Draw("MomentumTransfer>>trueQ2allhist");
    TH3D truevertexallhist("truevertexallhist","True Vertex",100,-150,150,100,-200,200,100,50,350);
    t3->Draw("MuonStartVertex.X():MuonStartVertex.Y():MuonStartVertex.Z()>>truevertexallhist");
    TH3D truedirallhist("truedirallhist","True Direction",100,-1,1,100,-1,1,100,-1,1);
    t3->Draw("MuonDirection.X():MuonDirection.Y():MuonDirection.Z()>>truedirallhist");
    
    
    // loop over reco trees and combine them:
    for(int i=0; i<t2->GetEntries();i++){
        t1->GetEntry(i);
        t2->GetEntry(i);
        outtruenuE=truenuE;
        outreconuE=reconuE;
        outtruemuE=truemuE;
        outrecomuE=recomuE;
        outtruevtx=TVector3(-1,-1,-1); // do not yet have
        outrecovtx=TVector3(recovtxx, recovtxy, recovtxz);
        outtruedir=TVector3(-1,-1,-1); // do not yet have
        outrecodir=TVector3(recomudirx, recomudiry, recomudirz);
        outtrueangle=-1;               // do not yet have
        outrecoangle=outrecodir.Angle(TVector3(0,0,1));
        recoanglehist.Fill(outrecoangle);
        outtrueq2=-1;                  
        outrecoq2=CalculateEventQ2(recomuE, reconuE, outrecoangle);
        recoQ2hist.Fill(outrecoq2);
        outnuEerr=outtruenuE-outreconuE;
        outmuEerr=outtruemuE-outrecomuE;
        outangerr=outtruedir.Angle(outrecodir);
        outvtxerr=(outtruevtx-outrecovtx).Mag();
        outq2err=outtrueq2-outrecoq2;
        outinfidvol=((sqrt(pow(outtruevtx.X(),2.)+pow(outtruevtx.Z()-tank_start-tank_radius,2.))<tank_radius)&&
                     (abs(outtruevtx.Y()-tank_yoffset)<tank_halfheight)&&((outtruevtx.Z()-tank_start-tank_radius)<0));
        outhasmrdtrack=(trackleninmrd>0);
        
        treeout->Fill();
    }
    
    truenuKEhist.Write();
    truemuKEhist.Write();
    trueanglehist.Write();
    trueQ2hist.Write();
    truevertexhist.Write();
    truedirhist.Write();
    
    reconuKEhist.Write();
    recomuKEhist.Write();
    recoanglehist.Write();
    recoQ2hist.Write();
    recovertexhist.Write();
    recodirhist.Write();
    
    truenuKEallhist.Write();
    truemuKEallhist.Write();
    trueangleallhist.Write();
    trueQ2allhist.Write();
    truevertexallhist.Write();
    truedirallhist.Write();
    
//    // loop over trueQevertexinfo file, because draw is acting weird.
//    for(int i=0; i<t3->GetEntries();i++){
//        mrdtlb->GetEntry(i);
//        cout<<"("<<mrdtracklen<<", "
//        cout<<"("<<muE;
//        //if(mrdtracklen>0)
//        ahist.Fill(muE*1000);
//        cout<<", "<<i<<") ";
//    }
    
    t1->ResetBranchAddresses();
    f1->Close();
    t2->ResetBranchAddresses();
    f2->Close();
    t3->ResetBranchAddresses();
    f3->Close();
    treeout->Write();
    treeout->ResetBranchAddresses();
    fileout->Close();
    delete fileout;
}

double CalculateEventQ2(double recoMuonEnergy, double recoNeutrinoEnergy, double recoMuonAngle){
    TDatabasePDG db;
    Double_t neutronmass = (db.GetParticle(2112)->Mass())*1000.;  // converted to MeV
    Double_t protonmass = (db.GetParticle(2212)->Mass())*1000.;   // converted to MeV
    Double_t muonmass = (db.GetParticle(13)->Mass())*1000.;       // converted to MeV
    Double_t O16bindingEnergy = 7.9762086875;                     // MeV (per nucleon), from http://tinyurl.com/y8m9s4z6
    Double_t boundneutronmass = neutronmass-O16bindingEnergy;

    double part1 = recoMuonEnergy - (sqrt(pow(recoMuonEnergy,2.)-pow(muonmass,2.))*TMath::Cos(recoMuonAngle));
    double eventq2 = -pow(muonmass,2.) + 2.*recoNeutrinoEnergy*part1;
    return eventq2;
}
