#ifndef WCSimLAPPDRESPONSE_HH
#define WCSimLAPPDRESPONSE_HH

#include "WCSimLAPPDpulse.hh"
#include "WCSimLAPPDpulseCluster.hh"
#include "TObject.h"
#include "TH1.h"
#include "TRandom3.h"
#include "globals.hh"

class WCSimLAPPDresponse : public TObject {

 public: 

  WCSimLAPPDresponse();

  ~WCSimLAPPDresponse();

  void AddSinglePhotonTrace(double trans, double para, double time);
	
  TH1D* GetTrace(int CHnumber, int parity, double starttime, double samplesize, int numsamples, double thenoise);

  int FindStripNumber(double trans);

  double StripCoordinate(int stripnumber);

  WCSimLAPPDpulseCluster* GetPulseCluster() {return _pulseCluster;}

 private: 



  //relevant to a particular event
  double _freezetime;

  //input parameters and distributions
  TH1D* _templatepulse;
  TH1D* _PHD;
  TH1D* _pulsewidth; 

  //output responses
  TH1D** StripResponse_neg;
  TH1D** StripResponse_pos;

  WCSimLAPPDpulseCluster* _pulseCluster;

  //randomizer
  TRandom3* mrand;	

  //useful functions
  int FindNearestStrip(double trans);
  double TransStripCenter(int CHnum);
  
  ClassDef(WCSimLAPPDresponse,0)

};

#endif
