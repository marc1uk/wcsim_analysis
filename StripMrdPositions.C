/* vim:set noexpandtab tabstop=4 wrap */
#include <string>
#include <regex>
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;
//.x /annie/app/users/moflaher/wcsim/root_work/RegexTest.C+	<< standalone call

void cMRDTrack::StripMrdPositions(std::string fname="/annie/app/users/moflaher/wcsim/root_work/MRD_positions_raw"){
	
	// input file reading infrastructure
	std::ifstream input_file(fname);
		if ( input_file.fail() ) {
		cerr << "File open error!" << endl;
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
	cout<<"pattern to match is "<<theexpressionstring<<endl;
	
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
	for(int linenum=0; linenum<filelines.size(); linenum++){
		std::string stringtomatch = filelines.at(linenum);
		std::regex_match (stringtomatch, submatches, theexpression);
		if(submatches.size()==0) { 
			cout<<stringtomatch<<" was not matched"<<endl; return;
		} else { 
			//cout << "whole match was "<<(std::string)submatches[0]<<endl; }
//			for(int i=1; i<submatches.size(); i++){
//				std::string tval = (std::string)submatches[i];
//				cout <<"submatch " <<i<<" is "<<tval<<endl;
//			}
			pmt_its.push_back(std::stoi(submatches[2]));
			(submatches[3]=="H") ? orientations.push_back(0) : orientations.push_back(1);
			layers.push_back(std::stoi(submatches[4]));
			originx.push_back(std::stod(submatches[4]));
			originy.push_back(std::stod(submatches[5]));
			originz.push_back(std::stod(submatches[6]));
			extentsx.push_back(std::pair<Double_t, Double_t>(std::stod(submatches[7]),std::stod(submatches[8])));
			extentsx.push_back(std::pair<Double_t, Double_t>(std::stod(submatches[9]),std::stod(submatches[10])));
			extentsx.push_back(std::pair<Double_t, Double_t>(std::stod(submatches[11]),std::stod(submatches[12])));
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
	
}
