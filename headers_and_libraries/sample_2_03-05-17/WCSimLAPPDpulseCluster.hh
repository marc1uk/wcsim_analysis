#ifndef WCSimLAPPDPULSECLUSTER_HH
#define WCSimLAPPDPULSECLUSTER_HH

#include "WCSimLAPPDpulse.hh"
#include "TObject.h"
#include <map>
#include "WCSimLAPPDInfo.hh"

class WCSimLAPPDpulseCluster : public TObject {

 public: 

  WCSimLAPPDpulseCluster();
  ~WCSimLAPPDpulseCluster();

  void Reset();

  void AddPulse(WCSimLAPPDpulse* pulse);
  int GetNPulsesStrip(int stripnum);
  int GetPulseNum(int stripnum, int pulsenum);
  WCSimLAPPDpulse* GetPulse(Int_t n);
  Int_t GetNPulses();

 private: 

  std::vector<WCSimLAPPDpulse*> fLAPPDpulseList;
  std::vector<double> fLAPPDstripPulseCount;
  std::vector<std::vector<int> > fLAPPDstripPulseCoordinate;
  //  std::map<int,std::vector<double> > fLAPPDstripPulseCoordinate;
  int _nchannels; //number of uniques channels

  ClassDef(WCSimLAPPDpulseCluster,0)

};

#endif
