#include <string>
#include <iostream>
#include <sstream>

int dome(){
  std::string stringtomatch = "PMT 0 : Orientation H : Layer 0 : Origin (737.5,-1218,3366.8) : Extent (1.5→1473.5, -1318→-1118, 3360.8→3366.8)";
//  std::string stringtomatch = "PMT 0";
  std::string theexpressionstring = "PMT %d : Orientation %c : Layer %d : Origin (%f,%f,%f) : Extent (%f→%f, %f→%f, %f→%f)";
//  std::string theexpressionstring = "PMT %d";

  int pmt_id=1;
  char next_paddle_orientation='?';
  int next_paddle_layer=-1;
  double next_paddle_originx=-1;
  double next_paddle_originy=-1;
  double next_paddle_originz=-1;
  double next_paddle_extentsx1=-1, next_paddle_extentsx2=-1;
  double next_paddle_extentsy1=-1, next_paddle_extentsy2=-1;
  double next_paddle_extentsz1=-1, next_paddle_extentsz2=-1;

  int nmatched = sscanf(stringtomatch.c_str(), theexpressionstring.c_str(), &pmt_id, &next_paddle_orientation, &next_paddle_layer, &next_paddle_originx, &next_paddle_originy, &next_paddle_originz, &next_paddle_extentsx1, &next_paddle_extentsx2, &next_paddle_extentsy1, &next_paddle_extentsy2, &next_paddle_extentsz1, &next_paddle_extentsz2);
//  int nmatched = sscanf(stringtomatch.c_str(), theexpressionstring.c_str(), &pmt_id);

  std::cout<<"matched "<<nmatched<<" elements"<<std::endl;
  return 0;

}
