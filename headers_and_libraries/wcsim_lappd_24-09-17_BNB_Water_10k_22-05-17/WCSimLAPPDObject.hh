#ifndef WCSimWCLAPPDObject_h
#define WCSimWCLAPPDObject_h 1

#include "WCSimDetectorConstruction.hh"
#include "globals.hh"
#include "Randomize.hh"
#include <map>
#include <vector>


class WCSimLAPPDObject
{

public:
  virtual G4String GetLAPPDName()=0;
  virtual G4double GetExposeHeight()=0;
  virtual G4double GetRadius()=0;
  virtual G4float* Getqpe()=0;
  virtual G4float* GetQE()=0;
  virtual G4float* GetQEWavelength()=0;
  virtual G4float  GetmaxQE()=0;
  virtual G4float  GetCollectionEfficiency(float);
  virtual float    HitTimeSmearing(float)=0;
  virtual G4double GetLAPPDGlassThickness()=0;
  virtual G4float  GetDarkRate()=0;
  virtual G4float  GetDarkRateConversionFactor()=0;
protected:
  virtual G4float* GetCollectionEfficiencyArray();
  virtual G4float* GetCollectionEfficiencyAngle();
  G4float Interpolate_func(G4float, G4int, G4float*, G4float*);
};

class LAPPD : public WCSimLAPPDObject
{

public:
  
  LAPPD();
  ~LAPPD();
 
public:
  G4String GetLAPPDName() ;
  G4double GetExposeHeight();
  G4double GetRadius();
  G4float* Getqpe();
  G4float* GetQE();
  G4float* GetQEWavelength();
  G4float  GetmaxQE();
  float    HitTimeSmearing(float);
  G4double GetLAPPDGlassThickness();
  G4float  GetDarkRate();
  G4float  GetDarkRateConversionFactor();


};

#endif
