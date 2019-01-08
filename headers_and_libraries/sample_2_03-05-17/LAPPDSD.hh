// ====================================================================
//   LAPPDSD.hh
//
//   
// ====================================================================
#ifndef LAPPD_SD_H
#define LAPPD_SD_H
 
#include "G4VSensitiveDetector.hh"
//#include "G4EmCalculator.hh"
#include "LAPPDHit.hh"
//#include "LAPPDResponse.hh"
#include "WCSimDetectorConstruction.hh"  

class G4HCofThisEvent;
class G4Step;				// ?? why do we need these?
class G4TouchableHistory;		// ??

class LAPPDSD : public G4VSensitiveDetector {

private:

  LAPPDHitsCollection* hitsCollection;
  LAPPDHit* newHit;
  G4int hcid;
  //WCLiteDetectorConstruction* detector;
  
  //LAPPDResponse* mrdresp;
  //SBsimLAPPDDB* mrddb;

  G4int         nhit;

  
  /*
  G4double      dedxtimesstep[NBUF_LAPPD];
  G4double      steplength[NBUF_LAPPD];
  G4double      ADCcount[NBUF_LAPPD];

  G4EmCalculator emcal;
  */
  
public:

  LAPPDSD(const G4String& name);
  virtual ~LAPPDSD();

  // virtual methods
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory*);
  virtual void Initialize(G4HCofThisEvent* HCTE);
  virtual void EndOfEvent(G4HCofThisEvent* HCTE);

  virtual void DrawAll();
  virtual void PrintAll(); 
 
};
#endif
