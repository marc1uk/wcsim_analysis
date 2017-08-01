// ====================================================================
//   SBsimLAPPDSD.hh
//
//   2006/03/03 K. Hiraide
// ====================================================================
#ifndef SBSIM_LAPPD_SD_H
#define SBSIM_LAPPD_SD_H
 
#include "G4VSensitiveDetector.hh"
#include "G4EmCalculator.hh"
#include "SBsimLAPPDHit.hh"
#include "SBsimLAPPDResponse.hh"

#include "SBsimLAPPDDB.hh" //added

#define NBUF_LAPPD 512

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

class SBsimLAPPDSD : public G4VSensitiveDetector {
private:
  SBsimLAPPDHitsCollection* hitsCollection;
  SBsimLAPPDResponse* mrdresp;
  SBsimLAPPDHit* newHit;
  //SBsimDetectorConstruction* detector;
  WCSimDetectorConstruction* detector;
  SBsimLAPPDDB* mrddb;

  G4int         nhit;
  G4int         id;
  G4int         initflag[NBUF_LAPPD];

  //  G4double      edepbuf[NBUF_LAPPD];
  G4double      edepbuf;
  G4double      edeposit[NBUF_LAPPD];
  G4double      edeposit_resp[NBUF_LAPPD];
  G4int         itrack[NBUF_LAPPD];
  G4int         ipart[NBUF_LAPPD];
  G4String      parname[NBUF_LAPPD];
  G4double      hittime[NBUF_LAPPD];
  G4double      detecttime[NBUF_LAPPD];
  G4ThreeVector hitpos[NBUF_LAPPD];
  
  G4double      dedxtimesstep[NBUF_LAPPD];
  G4double      steplength[NBUF_LAPPD];
  G4double      ADCcount[NBUF_LAPPD];

  G4EmCalculator emcal;

public:
  SBsimLAPPDSD(const G4String& name);
  virtual ~SBsimLAPPDSD();

  // virtual methods
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  virtual void Initialize(G4HCofThisEvent* HCTE);
  virtual void EndOfEvent(G4HCofThisEvent* HCTE);

  virtual void DrawAll();
  virtual void PrintAll(); 
 
};

#endif
