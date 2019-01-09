#include <iostream>     // std::cout
#include <algorithm>    // std::equal_range, std::sort
#include <vector>       // std::vector

template <typename T>
std::vector<T> valsinrange(typename std::vector<T>::iterator rangebegin, typename std::vector<T>::iterator rangeend, T lowerlimit, T upperlimit){
  auto lowerit = std::lower_bound(rangebegin, rangeend, lowerlimit);
  auto upperit = std::upper_bound(rangebegin, rangeend, upperlimit)-1;
  return typename std::vector<T> (lowerit, upperit);
}

void mymacro(){
  std::vector<double> myvec{0.,1.,2.,3.2,4.3,5.1,6.7,7.,8.};
  std::vector<double> myvec2 = valsinrange(myvec.begin(), myvec.end(), 2.5,6.8);
  cout<<"myvec2 is :"<<endl;
  for(auto anel : myvec2) cout<<anel<<", ";
  cout<<endl;
}
