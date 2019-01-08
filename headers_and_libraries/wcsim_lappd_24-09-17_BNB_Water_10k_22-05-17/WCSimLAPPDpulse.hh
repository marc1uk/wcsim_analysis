//#ifndef LAPPDPULSE_HH
//#define LAPPDPULSE_HH
#ifndef WCSimLAPPDPULSE_h
#define WCSimLAPPDPULSE_h 1

#include "TObject.h"
#include "TH1.h"
#include "WCSimLAPPDInfo.hh"

class WCSimLAPPDpulse : public TObject {

 public: 

  WCSimLAPPDpulse(double plusetime, double lefttime, double righttime, double peakvalue, int stripnum);

  ~WCSimLAPPDpulse();

  Double_t Getpulsetime() { return _pulsetime; }
  Double_t Getlefttime() { return _lefttime; }
  Double_t Getrighttime() { return _righttime; }
  Double_t Getpeakvalue() { return _peakvalue; }
  Double_t Getstripnum(){ return _stripnum; }


 private: 

  double _pulsetime,_lefttime,_righttime,_peakvalue;
  int _stripnum; 
 
  ClassDef(WCSimLAPPDpulse,0)

};

#endif
