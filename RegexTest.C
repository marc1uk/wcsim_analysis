#include <string>
#include <regex>
#include <iostream>
#include <iomanip>
// call: .x /annie/app/users/moflaher/wcsim/root_work/RegexTest.C+
void RegexTest(std::string fname="/pnfs/annie/persistent/users/moflaher/g4dirt/annie_tank_flux.2000.root"){
  // NOTE: ESCAPE SEQUENCES REQUIRE DOUBLE BACKSLASHES
  //std::string theexpressionstring = ".*/[^0-9]+\\.([0-9]+)\\.root";	// success!
  //std::string theexpressionstring = ".*/?[^\\.]+\\.([0-9]+)\\.root";
  std::string theexpressionstring = ".*/?[^\\.]+\\.([0-9]+)\\.?([0-9]+)?\\.root";
  cout<<"testing "<<fname<<" against "<<theexpressionstring<<endl;
  std::match_results<string::const_iterator> submatches;
  std::regex theexpression (theexpressionstring);
  std::regex_match ((std::string)fname, submatches, theexpression);
  
  if(submatches.size()>0){ cout << "whole match was "<<(std::string)submatches[0]<<endl; }
  else { cout<<"not matched"<<endl; return; }
  for(int i=1; i<submatches.size(); i++){
       std::string tval = (std::string)submatches[i];
       cout <<"submatch " <<i<<" is "<<tval<<endl;
  }
}
