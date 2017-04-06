/* vim:set noexpandtab tabstop=4 wrap */
#include <string>
#include <regex>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "TSystem.h"
using namespace std;
//.x /annie/app/users/moflaher/wcsim/root_work/RegexTest.C+	<< standalone call

int mrdcluster::StripMrdPositions(std::string fname="/annie/app/users/moflaher/wcsim/root_work/MRD_positions_raw"){
	
	TString pwd = gSystem->Getenv("PWD");
	TString rawfilename = pwd + "/MRD_positions_raw";
	
	// input file reading infrastructure
	//std::ifstream input_file(fname);
	std::ifstream input_file(rawfilename);
		if ( input_file.fail() ) {
		cerr << "Failed to open file " << rawfilename << "!" << endl;
		exit(-1);
	}
	char linebuffer[256];
	std::string fileline;
	std::vector<std::string> filelines;
	
	// read the file and put the lines into a vector of std::strings
	while ( !input_file.eof() ) {
		input_file.getline(linebuffer,256);
		//istringstream streamconverter(linebuffer);
		//streamconverter >> fileline;	// pulls out each element by space delimitation
		fileline = std::string(linebuffer);
		filelines.push_back(fileline);
	}
	
	// lines are of the format:
	// PMT 0 : Orientation H : Layer 0 : Origin (737.5,-1219.5,-586.65) : Extent (1.5→1473.5, -1319.5→-1119.5, 2665.35→2671.35)
	std::string theexpressionstring = "PMT ([0-9]+) : Orientation (.) : Layer ([0-9]+) : Origin \\(([0-9\\.\\+\\-]+),([0-9\\.\\+\\-]+),([0-9\\.\\+\\-]+)\\) : Extent \\(([0-9\\.\\+\\-]+)→([0-9\\.\\+\\-]+), ([0-9\\.\\+\\-]+)→([0-9\\.\\+\\-]+), ([0-9\\.\\+\\-]+)→([0-9\\.\\+\\-]+)\\)";
//	cout<<"pattern to match is "<<theexpressionstring<<endl;
	
	// declare the output vectors to put things in
//	std::vector<Int_t> mrd_pmt_ids;
//	std::vector<Int_t> mrd_orientations;	// H is 0, V is 1
//	std::vector<Int_t> mrd_layers;
//	std::vector<Double_t> mrd_originx;
//	std::vector<Double_t> mrd_originy;
//	std::vector<Double_t> mrd_originz;
//	std::vector<std::pair<Double_t,Double_t> > mrd_extentsx;
//	std::vector<std::pair<Double_t,Double_t> > mrd_extentsy;
//	std::vector<std::pair<Double_t,Double_t> > mrd_extentsz;
	
	// use regex to extract the information
	std::match_results<string::const_iterator> submatches;
	std::regex theexpression (theexpressionstring);
	int numlayers=0;
	for(int linenum=0; linenum<filelines.size(); linenum++){
		std::string stringtomatch = filelines.at(linenum);
//		cout<<"next string is "<<stringtomatch<<endl;
		if(stringtomatch=="") break;
		std::regex_match (stringtomatch, submatches, theexpression);
		if(submatches.size()==0) { 
			cout<<stringtomatch<<" was not matched"<<endl; return 0;
		} else {
//			cout << "whole match was "<<(std::string)submatches[0]<<endl;
//			for(int i=1; i<submatches.size(); i++){
//				std::string tval = (std::string)submatches[i];
//				cout <<"submatch " <<i<<" is "<<tval<<endl;
//			}
			int pmt_id=std::stoi(submatches[1]);
			//pmt_ids.push_back(pmt_id);	// 1:1 mapping, no need
			//cout<<"setting stats for pmt "<<pmt_id<<endl;
			(submatches[2]=="H") ? paddle_orientations.push_back(0) : paddle_orientations.push_back(1);
			paddle_layers.at(pmt_id)=(std::stoi(submatches[3]));
			paddle_originx.at(pmt_id)=(std::stod(submatches[4]));
			paddle_originy.at(pmt_id)=(std::stod(submatches[5]));
			paddle_originz.at(pmt_id)=(std::stod(submatches[6]));
			paddle_extentsx.at(pmt_id)=(std::pair<Double_t, Double_t>(std::stod(submatches[7]),std::stod(submatches[8])));
			paddle_extentsy.at(pmt_id)=(std::pair<Double_t, Double_t>(std::stod(submatches[9]),std::stod(submatches[10])));
			paddle_extentsz.at(pmt_id)=(std::pair<Double_t, Double_t>(std::stod(submatches[11]),std::stod(submatches[12])));
		}
	}
	// success should return for each match:
	//submatch 1 is 0
	//submatch 2 is H
	//submatch 3 is 0
	//submatch 4 is 737.5
	//submatch 5 is -1219.5
	//submatch 6 is -586.65
	//submatch 7 is 1.5
	//submatch 8 is 1473.5
	//submatch 9 is -1319.5
	//submatch 10 is -1119.5
	//submatch 11 is 2665.35
	//submatch 12 is 2671.35
	
//	std::ofstream mrdpositions;
//	mrdpositions.open("mrdpositions.txt", std::ios::out);
//	mrdpositions<<"pmtids is "<<pmt_ids<<endl;
//	mrdpositions<<"orientations is "<<orientations<<endl;
//	mrdpositions<<"layers is "<<layers<<endl;
//	mrdpositions<<"originxs is "<<originx<<endl;
//	mrdpositions<<"originys is "<<originy<<endl;
//	mrdpositions<<"originzs is "<<originz<<endl;
//	mrdpositions<<"extentsx is "<<extentsx<<endl;
//	mrdpositions<<"extentsy is "<<extentsy<<endl;
//	mrdpositions<<"extentsz is "<<extentsz<<endl;
//	mrdpositions.close();
	
//	for(int i=0; i<paddle_orientations.size(); i++){
//		cout<<"orientations is "<<paddle_orientations.at(i)<<endl;
//		cout<<"layers is "<<paddle_layers.at(i)<<endl;
//		cout<<"originxs is "<<paddle_originx.at(i)<<endl;
//		cout<<"originys is "<<paddle_originy.at(i)<<endl;
//		cout<<"originzs is "<<paddle_originz.at(i)<<endl;
//		cout<<"extentsx is "<<paddle_extentsx.at(i).first<<", "<<paddle_extentsx.at(i).second<<endl;
//		cout<<"extentsy is "<<paddle_extentsy.at(i).first<<", "<<paddle_extentsx.at(i).second<<endl;
//		cout<<"extentsz is "<<paddle_extentsz.at(i).first<<", "<<paddle_extentsx.at(i).second<<endl;
//	}
	return 1;
}
