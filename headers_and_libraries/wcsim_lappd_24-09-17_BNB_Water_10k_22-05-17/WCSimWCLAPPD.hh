#ifndef WCSimWCLAPPD_h
#define WCSimWCLAPPD_h 1

//#include "WCSimDarkRateMessenger.hh"
#include "WCSimDetectorConstruction.hh"
#include "G4VDigitizerModule.hh"
#include "WCSimWCDigi.hh"
#include "WCSimWCHit.hh"
#include "globals.hh"
#include "Randomize.hh"
#include <map>
#include <vector>

#include "WCSimLAPPDpulse.hh"
#include "WCSimLAPPDpulseCluster.hh"
#include "TObject.h"
#include "TH1.h"
#include "TFile.h"
#include "TRandom3.h"
#include "WCSimLAPPDInfo.hh"
#include "WCSimLAPPDObject.hh"

//This class is a copy of  LAPPDresponse.hh 
class WCSimWCLAPPD  : public G4VDigitizerModule{//: public TObject {

public:
  WCSimWCLAPPD(G4String name, WCSimDetectorConstruction*);
  ~WCSimWCLAPPD();

  void ReInitialize() { DigiHitMapLAPPD.clear(); TriggerTimes.clear(); }
  
  void AddSinglePhotonTrace(double trans, double para, double time);
  TH1D GetTrace(int CHnumber, int parity, double starttime, double samplesize, int numsamples, double thenoise);
  int FindStripNumber(double trans);
  double StripCoordinate(int stripnumber);
  WCSimLAPPDpulseCluster* GetPulseCluster() {return _pulseCluster;}
  
public:
  
  void AddLAPPDDarkRate(WCSimWCDigitsCollection*);
  void MakePeCorrection_lappd(WCSimWCHitsCollection*);
  void Digitize();
  G4double GetTriggerTime(int i) { return TriggerTimes[i];}
  // void SetConversion(double iconvrate){ ConvRate = iconvrate; }
  //  static G4double GetLongTime() { return LongTime;}
  
  G4double rn1pe();
  G4double peSmeared;
  // double ConvRate; // kHz
  
  std::vector<G4double> TriggerTimes;
  std::map<int,int> DigiHitMapLAPPD; // need to check if a hit already exists..
	
  WCSimWCDigitsCollection*  DigitsCollection;  
  WCSimDetectorConstruction* myDetector;
  
  static TFile tf;
  static int lappobjectcounter;
 
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
  
  //ClassDef(WCSimWCLAPPD,0)
};

/*class WCSimWCLAPPD : public G4VDigitizerModule
{
public:
  
  WCSimWCLAPPD(G4String name, WCSimDetectorConstruction*);
  ~WCSimWCLAPPD();
  
   void ReInitialize() { DigiHitMapPMT.clear(); TriggerTimes.clear(); }
    
   
public:
  
  void AddPMTDarkRate(WCSimWCDigitsCollection*);
  void MakePeCorrection(WCSimWCHitsCollection*);
  void Digitize();
  G4double GetTriggerTime(int i) { return TriggerTimes[i];}
  // void SetConversion(double iconvrate){ ConvRate = iconvrate; }
  //  static G4double GetLongTime() { return LongTime;}
  
  G4double rn1pe();
  G4double peSmeared;
  // double ConvRate; // kHz
  std::vector<G4double> TriggerTimes;
  std::map<int,int> DigiHitMapPMT; // need to check if a hit already exists..

  WCSimWCDigitsCollection*  DigitsCollection;  
  WCSimDetectorConstruction* myDetector;

};
*/
#endif








