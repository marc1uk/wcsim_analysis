#ifndef WCSimSteppingAction_h
#define WCSimSteppingAction_h 1

#include "G4Event.hh"
#include "G4UserSteppingAction.hh"
#include "G4ThreeVector.hh"
#include "G4OpBoundaryProcess.hh"

class G4HCofThisEvent;
class G4Event;

class WCSimSteppingAction : public G4UserSteppingAction
{

public:
  WCSimSteppingAction() : fExpectedNextStatus(Undefined)
  {};

  ~WCSimSteppingAction()
  {};

  void UserSteppingAction(const G4Step*);

  G4int G4ThreeVectorToWireTime(G4ThreeVector *pos3d,
				G4ThreeVector lArPos,
				G4ThreeVector start,
				G4int i);
  
  void Distortion(G4double x,
		  G4double y);

  G4double FieldLines(G4double x,
		      G4double y,
		      G4int xy);

private:

  G4double ret[2];
  G4OpBoundaryProcessStatus fExpectedNextStatus;
  G4String ToName(G4OpBoundaryProcessStatus boundaryStatus);
  G4String ToName2(G4StepStatus stepStatus);

};


#endif
