// ====================================================================
//   SBsimLAPPDHit.hh
//
//   2006/03/03 K. Hiraide
// ====================================================================
#ifndef SBSIM_LAPPD_HIT_H
#define SBSIM_LAPPD_HIT_H

#include "G4ThreeVector.hh"
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class SBsimLAPPDHit : public G4VHit {

public:
  SBsimLAPPDHit();
  ~SBsimLAPPDHit();

  // copy constructor & assignment operator
  SBsimLAPPDHit(const SBsimMRDHit& right);

  const SBsimLAPPDHit& operator=(const SBsimMRDHit& right);
  G4int operator==(const SBsimLAPPDHit& right) const;
  
  // new/delete operators
  void* operator new(size_t);
  void operator delete(void* aHit);
  
  // set/get functions
  void SetEdeposit (G4double aedeposit)    { edeposit = aedeposit;}
  void SetID       (G4int aid)             { id       = aid;      }
  void SetITrack   (G4int aitrack)         { itrack   = aitrack;  }
  void SetIPart    (G4int aipart)          { ipart    = aipart;   }
  void SetParticle (G4String aparname)     { parname  = aparname; }
  void SetHitTime  (G4double ahittime)     { hittime  = ahittime; }
  void SetDetectTime  (G4double adetecttime) { detecttime  = adetecttime; }
  void SetHitPos   (G4ThreeVector ahitpos) { hitpos   = ahitpos;  }
  void SetRot(G4RotationMatrix rmat) { rot = rmat; }

  void SetLogV(G4LogicalVolume* val) { pLogV = val; }
  const G4LogicalVolume* GetLogV() const { return pLogV; }

  void SetADCcount  (G4double aadccount)    { adccount  = aadccount; }
  void SetTotalPath (G4double atotalpath)   { totalpath = atotalpath;}


  G4double      GetEdeposit()  { return edeposit;}
  G4int         GetID()        { return id;      }
  G4int         GetTrackID()   { return itrack;  }
  G4int         GetParticleID(){ return ipart;   }
  G4String      GetParticle()  { return parname; }
  G4double      GetHitTime()   { return hittime; }
  G4double      GetDetectTime()   { return detecttime; }
  G4ThreeVector GetHitPos()    { return hitpos;  }
  G4RotationMatrix GetRot() const { return rot; }

  G4double      GetADCcount()  { return adccount; }
  G4double      GetTotalPath() { return totalpath;}

  // methods
  virtual void Draw();
  virtual void Print();

private:
  G4double edeposit;
  G4int id;
  G4int itrack;
  G4int ipart;
  G4String parname;
  G4double hittime;
  G4double detecttime;
  G4ThreeVector hitpos;

  G4double adccount;
  G4double totalpath;
  
  const G4LogicalVolume* pLogV;
  G4RotationMatrix rot;

};

// ====================================================================
// inline functions
// ====================================================================
inline SBsimLAPPDHit::SBsimMRDHit
       (const SBsimLAPPDHit& right)
  : G4VHit()
{
  edeposit= right.edeposit;
  id      = right.id;
  itrack  = right.itrack;
  ipart   = right.ipart;
  parname = right.parname;
  hittime = right.hittime;
  detecttime = right.detecttime;
  hitpos  = right.hitpos;
  rot     = right.rot;
  pLogV   = right.pLogV;

  adccount = right.adccount;
  totalpath= right.totalpath;
}

inline const SBsimLAPPDHit& SBsimMRDHit::operator=
       (const SBsimLAPPDHit& right)
{
  edeposit= right.edeposit;
  id      = right.id;
  itrack  = right.itrack;
  ipart   = right.ipart;
  parname = right.parname;
  hittime = right.hittime;
  detecttime = right.detecttime;
  hitpos  = right.hitpos;
  rot     = right.rot;
  pLogV   = right.pLogV;

  adccount = right.adccount;
  totalpath= right.totalpath;

  return *this;
}

inline G4int SBsimLAPPDHit::operator==
       (const SBsimLAPPDHit& right) const 
{
   return (this==&right) ? 1 : 0; 
}

// externally instanciated.
typedef G4THitsCollection<SBsimLAPPDHit> SBsimMRDHitsCollection;
extern G4Allocator<SBsimLAPPDHit> SBsimMRDHitAllocator; 

inline void* SBsimLAPPDHit::operator new(size_t)
{
  void* aHit= (void*)SBsimLAPPDHitAllocator.MallocSingle();
  return aHit;
}

inline void SBsimLAPPDHit::operator delete(void* aHit)
{
  SBsimLAPPDHitAllocator.FreeSingle((SBsimMRDHit*) aHit);
}

#endif
