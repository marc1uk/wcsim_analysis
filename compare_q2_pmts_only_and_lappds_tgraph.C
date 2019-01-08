#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TString.h"
#include "TMath.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TText.h"
#include "TColor.h"
#include "TStyle.h"
#include "TF1.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include <exception>	// for stdexcept
#include <vector>
#include <map>
#include <string>
#include <algorithm>	// remove and remove_if
#include <iostream>
#include <iomanip>
#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds
#include <time.h>         // clock_t, clock, CLOCKS_PER_SEC
//#include "/annie/app/users/moflaher/annie_rootlogon.C"

void TCenterTitles(TGraphErrors* histo);
void TSimulation();

void compare_q2_pmts_only_and_lappds_tgraph(bool relative=false, bool makeintermediates=false){
TFile* _file0 = TFile::Open("out/q2comparison/30ps_proper_split/pmts_only/q2comparison.root");
TFile* _file1 = TFile::Open("out/q2comparison/30ps_proper_split/pmts_and_lappds/q2comparison.root");
TTree* tp = (TTree*)_file0->Get("q2tree");
TTree* tpl = (TTree*)_file1->Get("q2tree");
TString outdirtstring = TString::Format("%s/out/q2comparison/30ps_proper_split/temp",gSystem->pwd());
cout<<"outdirtstring="<<outdirtstring.Data()<<endl;
const char* outdir = outdirtstring.Data();
cout<<"outdir="<<outdir<<endl;

std::vector<TTree*> configtrees{tp,tpl};
std::vector<std::string> confignames{"PMTs","PMTS+LAPPDs"};
std::vector<TGraphErrors*> tgraphs(configtrees.size());

std::vector<double> q2binedges{0.,0.25,0.5,0.75};
double abserrfor68pc;
double onesigmathresh=0.68;
int runningsum;
TH1D* errh=nullptr;
int nbinserrh=50;
double q2errlowlim=0., q2erruplim=1.2; // range of histogram to plot absolute q2 errors
double relq2errlowlim=0., relq2erruplim=3.;
TCanvas* c1 = new TCanvas();
TCanvas* c2 = new TCanvas();
std::vector<TH1D*> errhists(confignames.size(),nullptr);
std::vector<TH1D*> allhists;
double maxyaxisheight=0;
int tallestindex=-1;

for(int tgraphbini=0; tgraphbini<q2binedges.size(); tgraphbini++){
  for(int config=0; config<configtrees.size(); config++){
    TTree* nextree = configtrees.at(config);
    assert(nextree!=nullptr);
    TGraphErrors* nextgraph = tgraphs.at(config);
    assert(tgraphbini==0||nextgraph!=nullptr);
    if(nextgraph==nullptr) nextgraph = new TGraphErrors(q2binedges.size());
    tgraphs.at(config) = nextgraph;
    
    double upperlimit, lowerlimit;
    lowerlimit = q2binedges.at(tgraphbini);
    if((tgraphbini+1)<q2binedges.size()) upperlimit = q2binedges.at(tgraphbini+1);
    else upperlimit = -1.;
    cout<<"bin "<<tgraphbini<<" from "<<lowerlimit<<" to "<<upperlimit<<endl;
    
    if (errh==nullptr){ 
      if(relative==false){
        errh = new TH1D("errh","Error Hist",nbinserrh,q2errlowlim,q2erruplim);
      } else {
        errh = new TH1D("errh","Error Hist",nbinserrh,relq2errlowlim,relq2erruplim);
      }
    }
    
    TString cutstring;
    if(upperlimit>0){
      cutstring = TString::Format("TrueQ2>%f&&TrueQ2<%f",lowerlimit,upperlimit);
    } else{
      cutstring = TString::Format("TrueQ2>%f",lowerlimit);
    }
    cout<<"cutstring is: "<<cutstring.Data()<<endl;
    c1->cd();
    if(relative==false){
      nextree->Draw("(abs(Q2Err))>>errh",cutstring);
    } else {
      nextree->Draw("((abs(Q2Err))/TrueQ2)>>errh",cutstring);
    }
    double entriesthishist68pc = errh->GetEntries()*onesigmathresh;
    cout<<"histogram has "<<errh->GetEntries()<<" entries, setting 68% at "<<entriesthishist68pc<<endl;
    cout<<"underflow contents: "<<errh->GetBinContent(0)<<", overflow contents: "<<errh->GetBinContent(nbinserrh+1)<<endl;
    //c1->Update();
    //std::this_thread::sleep_for (std::chrono::seconds(3));
    //gPad->WaitPrimitive();
    
    runningsum=0;
    int histbini;
    for(histbini=1; histbini<nbinserrh; histbini++){
      runningsum += errh->GetBinContent(histbini);
      if(runningsum>=(entriesthishist68pc)) break;
    }
    abserrfor68pc = errh->GetBinLowEdge(histbini);
    if(upperlimit<0) upperlimit = lowerlimit + (lowerlimit-*(q2binedges.end()-2));
    
    
    
    double xvalerr = (upperlimit-lowerlimit)/2.;
      //sqrt(pq/n)*central value
    double yvalerr = sqrt((onesigmathresh*(1.-onesigmathresh))/errh->GetEntries()) * abserrfor68pc/*(upperlimit+lowerlimit)/2.*/;
    cout<<"broke loop at bin "<<histbini<<"/"<<nbinserrh<<" with "<<runningsum<<" entries accumulated"<<endl;
    cout<<"point errors are "<<xvalerr<<" x and "<<yvalerr<<" y"<<endl;
    cout<<"making point at (" <<(upperlimit+lowerlimit)/2.<<", "<<abserrfor68pc<<")"<<endl;
    nextgraph->SetPoint(tgraphbini, (upperlimit+lowerlimit)/2.,abserrfor68pc);
    nextgraph->SetPointError(tgraphbini, xvalerr,yvalerr);
    
    
    if(makeintermediates){
      errh->SetTitle(confignames.at(config).c_str());
      //errh->SetName(confignames.at(config).c_str()); //  TTREE->DRAW(STUFF>>ERRH) STOPS WORKING
      if(relative==false){
        errh->GetXaxis()->SetTitle("|#Delta Q^{2}| [(GeV/c)^{2}]");
      } else {
        //errh->GetXaxis()->SetTitle("Relative |#Delta Q^{2}| [no units]");
        errh->GetXaxis()->SetTitle("|(Q^{2}_{true}-Q^{2}_{reco})| / Q^{2}_{true}");
      }
      errh->GetYaxis()->SetTitle("Normalized Num Events");
      if(config==0) errh->SetLineColor(kBlue);
      else errh->SetLineColor(kRed);
      
      if(config==0) maxyaxisheight=0;
      double thismax = (errh->GetBinContent(errh->GetMaximumBin()))/(errh->Integral());
      cout<<"maximum bin contents are "<<errh->GetBinContent(errh->GetMaximumBin())
          <<", integral is "<<errh->Integral()<<", thismax is "<<thismax<<", oldmax is "<<maxyaxisheight<<endl;
      if(thismax>maxyaxisheight){ maxyaxisheight=thismax; tallestindex = config; }
      
      if(config==0&&tgraphbini!=0){
        for(auto ahist : errhists) delete ahist;
      }
      errhists.at(config) = new TH1D(*errh);
      
      c2->cd();
//      // can't do this here because we need to call Draw after setting y axis, and we don't know what
//      // range to use until (or order to draw) we've drawn them all.
//      if(config==0) errhists.at(config)->DrawNormalized("",100);
//      else errhists.at(config)->DrawNormalized("same",100);
      if(config==(confignames.size()-1)){
        errhists.at(tallestindex)->DrawNormalized("",100);
        for(int ehisti=0; ehisti<errhists.size(); ehisti++){
          if(ehisti==tallestindex) continue;
          else errhists.at(ehisti)->DrawNormalized("same",100);
        }
        TLegend* leg = c2->BuildLegend();
        leg->SetFillStyle(0);
        leg->SetLineStyle(0);
        c2->Modified();
        c2->Update();
        //TString outloc = TString::Format("%s/Q2_Subplot_%.2fMeVc-%.2fMeVc.png",outdir,lowerlimit,upperlimit);
        TString outloc = TString::Format("%s/Q2_Subplot_%d_%d.png",outdir,tgraphbini,config);
        c2->SaveAs(outloc);
      }
    }
    
  }
  
}

TGraphErrors* pq2e = tgraphs.at(0);
TGraphErrors* plq2e = tgraphs.at(1);
pq2e->SetMarkerColor(kBlue);
pq2e->SetMarkerStyle(25); //open box
pq2e->SetFillStyle(0);
pq2e->SetTitle("128 PMTs");
plq2e->SetMarkerColor(kRed);
plq2e->SetMarkerStyle(20);
plq2e->SetMarkerSize(0.8);
plq2e->SetFillStyle(0);
plq2e->SetTitle("5 LAPPDs + 128 PMTs");
if(relative==false){
  pq2e->GetYaxis()->SetTitle("|#Delta Q^{2}| [(GeV/c)^{2}]");
  pq2e->GetYaxis()->SetRangeUser(0.,0.6);
} else {
  //pq2e->GetYaxis()->SetTitle("Relative |#Delta Q^{2}| [no units]");
  pq2e->GetYaxis()->SetTitle("|(Q^{2}_{true}-Q^{2}_{reco})| / Q^{2}_{true}");
  pq2e->GetYaxis()->SetRangeUser(0.,1.5);
}
pq2e->GetXaxis()->SetTitle("True Q^{2} [(GeV/c)^{2}]");
TCenterTitles(pq2e);
pq2e->Draw("ap");
plq2e->Draw("p  same");
TLegend* leg = c1->BuildLegend();
leg->SetFillStyle(0);
leg->SetLineStyle(0);
for(int legi=0; legi<leg->GetListOfPrimitives()->GetEntries(); legi++){
  TLegendEntry* lege=(TLegendEntry*)leg->GetListOfPrimitives()->At(legi);
  lege->SetLineStyle(0);
  lege->SetLineColor(0);
  lege->SetLineWidth(0);
}
TSimulation();

}

void TCenterTitles(TGraphErrors* histo){
  histo->GetXaxis()->CenterTitle();
  histo->GetYaxis()->CenterTitle();
  //histo->GetZaxis()->CenterTitle();
}

void TSimulation(){
  TLatex* prelim = new TLatex(.9, .95, "ANNIE Simulation");
  prelim->SetTextColor(kGray+1);
  prelim->SetNDC();
  prelim->SetTextSize(2/30.);
  prelim->SetTextAlign(32);
  prelim->Draw();

}

