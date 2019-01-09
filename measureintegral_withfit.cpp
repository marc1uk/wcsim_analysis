/* vim:set noexpandtab tabstop=4 wrap */
// compile with:
// g++ measureintegral_withfit.cpp langaus_TGraph.C `root-config --libs` `root-config --cflags` -g -o measureintegral_withfit
// call with:
// ./measureintegral wave0.txt

#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <thread>
#include <chrono>
#include "TCanvas.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TF1.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"

#ifndef DRAWGRAPHS
//#define DRAWGRAPHS
#endif

int LoadFile(std::string filepath, std::vector<double> &data);
TF1 *langaufit(TGraph *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF);
// function which represents the actual mapping of x->y (the 'fit formula')
Double_t langaufun(Double_t *x, Double_t *par);
Double_t combinedfunc(Double_t* x, Double_t* par);

int main(int argc, char* argv[]){
	
	std::string waveformfile = argv[1];
	
	int myargc=0;
	char anarg[] = "spudleymunchkin";
	char *myargv[1] = {anarg};
	TApplication *app = new TApplication("app", &myargc, myargv);
	
	std::vector<double> data;
	int loadok = LoadFile(waveformfile, data);
	std::cout<<"file loaded"<<std::endl;
	
	double maxvalue = (*(std::max_element(data.begin(),data.end())));
	double minvalue = (*(std::min_element(data.begin(),data.end())));
	TH1D* hwaveform = new TH1D("hwaveform","Sample Distribution",(maxvalue-minvalue),minvalue,maxvalue);
	for(double sampleval : data) hwaveform->Fill(sampleval);
#ifdef DRAWGRAPHS
	TCanvas c2;
	hwaveform->Draw();
#else
	hwaveform->Draw("goff");
#endif
	hwaveform->Fit("gaus");
	TF1* gausfit = hwaveform->GetFunction("gaus");
	double gauscentre = gausfit->GetParameter(1); // 0 is amplitude, 1 is centre, 2 is width
	std::cout<<"pedestal is "<<gauscentre<<std::endl;
	
	std::vector<double> datainverted;
	for(auto datait = data.begin(); datait!=data.end(); datait++){ *datait -= gauscentre; *datait *= -1; }
	std::cout<<"pedestal subtracted"<<std::endl;
	
	uint16_t peaksample = std::distance(data.begin(), std::max_element(data.begin(),data.end()));
	double peakamplitude = data.at(peaksample);
	std::cout<<"peak sample is at "<<peaksample<<std::endl;
	
	// for *SPEED* trim the data window to just around the pulse
	double loosefitsamples = 100;
	std::vector<double> datatrimmed;
	for(int samplenum=std::max(0.,peaksample-loosefitsamples);
			samplenum<std::min(double(data.size()),peaksample+loosefitsamples);
			samplenum++){
		datatrimmed.push_back(data.at(samplenum));
	}
	data = datatrimmed;
	peaksample = std::distance(data.begin(), std::max_element(data.begin(),data.end()));
	std::cout<<"made trimmed version"<<std::endl;
	
	// from the peak sample, scan backwards for the start sample
	int startmode = 2;
	double thresholdfraction = 0.1;  // 20% of max
	// 3 ways to do this:
	// 0. scan away from peak until the first sample that changes direction
	// 1. scan away from peak until the first sample that crosses pedestal
	// 2. scan away from peak until we drop below a given fraction of the amplitude
	uint16_t startsample = peaksample-1;
	while(startsample>0){
		if( ( (startmode==0) && (data.at(startsample) > data.at(startsample-1)) ) ||
			( (startmode==1) && (data.at(startsample) > 0) ) ||
			( (startmode==2) && (data.at(startsample) > (peakamplitude*thresholdfraction)) ) ){
			startsample--;
		} else {
			break;
		}
	}
	std::cout<<"pulse starts at startsample "<<startsample<<std::endl;
	
	// from peak sample, scan forwards for the last sample
	int endmode = 2;
	uint16_t endsample = peaksample+1;
	while(uint16_t(endsample+1)<data.size()){
		if( ( (endmode==0) && (data.at(endsample+1) > data.at(endsample)) ) ||
			( (endmode==1) && (data.at(endsample) > 0) ) ||
			( (endmode==2) && (data.at(endsample) > (peakamplitude*thresholdfraction)) ) ){
			endsample++;
		} else {
			break;
		}
	}
	std::cout<<"pulse ends at endsample "<<endsample<<std::endl;
	
	// Landaugaus fit
	std::cout<<"doing langaus fit"<<std::endl;
	std::vector<double> numberline(data.size());
	std::iota (numberline.begin(),numberline.end(),0);
	TGraph* waveform = new TGraph(data.size(),numberline.data(), data.data());
#ifdef DRAWGRAPHS
	TCanvas c1;
	waveform->Draw("alp");
#else
	waveform->Draw("goff");
#endif
	
	// fit range: lower limit, upper limit
	std::vector<double> fitrange;
	// Fit parameters: landau width, landau centre, integral, width of gaussian component
	std::vector<double> start_vals;
	// lower limits on fit parameter ranges
	std::vector<double> param_lowerlims;
	// upper limits on fit parameter ranges
	std::vector<double> param_upperlims;
	
	fitrange=std::vector<double>{startsample-loosefitsamples,endsample+loosefitsamples};
	start_vals = std::vector<double>{5,static_cast<double>(peaksample),waveform->Integral(),1.};
	param_lowerlims = std::vector<double>{0.,startsample-10.,
										  0.1*waveform->Integral(),0.01};
	param_upperlims = std::vector<double>{10.,endsample+loosefitsamples,
										  100.*waveform->Integral(),20.};
	
	// fit results
	std::vector<double> fit_vals(4), fit_val_errors(4);
	Double_t chisqr;                  // fit chi^2
	Int_t ndf;                        // fit NDF
	
	// do the langaus fit
	TF1 *fitlangaus = langaufit(waveform,fitrange.data(),start_vals.data(),param_lowerlims.data(),
							param_upperlims.data(),fit_vals.data(),fit_val_errors.data(),&chisqr,&ndf);
	//fitlangaus->Draw("same");
	
	// add a linear component, fitting over a shorter range to better fit the local variation of pedestal
	std::cout<<"re-fitting with a linear component"<<std::endl;
	double loosefitsamples2 = 30;
	std::vector<double> fitrange2=std::vector<double>{startsample-loosefitsamples2,endsample+loosefitsamples2};
	double fitm = 0;
	double fitc = -peaksample;
	std::vector<double> allparams(fit_vals);
	allparams.push_back(fitm);
	allparams.push_back(fitc);
	TF1* fit = new TF1("fit",combinedfunc,fitrange2.front(),fitrange2.back(),allparams.size());
	fit->SetParameters(allparams.data());
	fit->SetParName(0,"landau width");
	fit->SetParName(1,"landau centre");
	fit->SetParName(2,"integral of landaugaus");
	fit->SetParName(3,"width of centre gaus component");
	fit->SetParName(4,"linear gradient");
	fit->SetParName(5,"linear offset");
	TFitResultPtr fitresult = waveform->Fit("fit","MRS"); // M: better fit, R: use range, S: get resultsptr
	std::cout<<"got second fit"<<std::endl;
	
	// calculate the integral
	// TODO
	
#ifdef DRAWGRAPHS
	c1.cd();
	fit->Draw("same");
	//gPad->WaitPrimitive();
	while(1){
		std::this_thread::sleep_for (std::chrono::seconds(1));
		gSystem->ProcessEvents();
	}
#endif
	
	// cleanup
	delete fit;
	delete waveform;
	delete hwaveform;
	app->Terminate();
	delete app;
	
	std::cout<<"returning"<<std::endl;
	return 1;
}

int LoadFile(std::string filepath, std::vector<double> &data){
	std::ifstream fin (filepath.c_str());
	std::string Line;
	std::stringstream ssL;
	
	int value;
	while (getline(fin, Line))
	{
		if (Line.empty()) continue;
		//if (Line.find("EndChannelMask") != std::string::npos) maskSettings = false;
		//if (Line[0] == '#') continue;
		else
		{
			ssL.str("");
			ssL.clear();
			ssL << Line;
			if (ssL.str() != "")
			{
				ssL >> value;
				data.push_back(value);
			}
		}
	}
	fin.close();
	return 1;
}

Double_t combinedfunc(Double_t* x, Double_t* par){
	return (langaufun(x, par) + par[4]*((*x)-par[5]));
}
