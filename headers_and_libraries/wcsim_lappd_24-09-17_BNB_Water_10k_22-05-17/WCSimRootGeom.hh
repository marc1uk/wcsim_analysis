#ifndef WCSim_RootGeom
#define WCSim_RootGeom

//////////////////////////////////////////////////////////////////////////
//                                                                      
// WCSim_RootGeom                                                      
//                                                                      
// This class contains information needed to be passed to reconstruction
//     routines.  It's just simple right now-- only the bare-bones  
//     WC info
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TClonesArray.h"

class TDirectory;

//////////////////////////////////////////////////////////////////////////

class WCSimRootPMT : public TObject {

private:
  Int_t fTubeNo;
  //Int_t fLAPPDNo;
  Int_t fCylLoc;  // endcap1, wall, endcap2
  Float_t fOrientation[3];
  Float_t fPosition[3];
  std::string fPmtTypeName;

public:
  WCSimRootPMT();
  WCSimRootPMT(Int_t tubeNo, Int_t cylLoc, Float_t orientation[3], Float_t position[3], std::string PmtType);
  virtual ~WCSimRootPMT();

  void  SetTubeNo(Int_t i) {fTubeNo=i;}
  //void  SetLAPPDNo(Int_t i) {fTubeNo=i;}
  void  SetCylLoc(Int_t i) {fCylLoc=i;}
  void  SetOrientation(Int_t i, Float_t f) {fOrientation[i]= ( (i<3) ? f : 0);}
  void  SetPosition(Int_t i, Float_t f) {fPosition[i]= ( (i<3) ? f : 0);}
  void SetType(std::string PmtType) {fPmtTypeName=PmtType;}

  Int_t GetTubeNo() const {return fTubeNo;}
  //Int_t GetLAPPDNo() const {return fLAPPDNo;}
  Int_t GetCylLoc() const {return fCylLoc;}
  Float_t GetOrientation(Int_t i=0) {return (i<3) ? fOrientation[i] : 0;}
  Float_t GetPosition(Int_t i=0) {return (i<3) ? fPosition[i] : 0;}
  std::string GetName(){ return fPmtTypeName;}

  ClassDef(WCSimRootPMT,1)  //WCSimPMT structure
};


//////////////////////////////////////////////////////////////////////////

class WCSimRootGeom : public TObject {

private:

	static const Int_t     maxNumPMT = 40000;
  static const Int_t     maxNumLAPPD = 100;
  Float_t                fWCCylRadius;  // Radius of WC tank
  Float_t                fWCCylLength;  // Length of WC tank
  
  Int_t                  fgeo_type;  // mailbox or cylinder?

  Float_t                fWCPMTRadius; // Radius of PMT
  Int_t                  fWCNumPMT;   // Number of PMTs
  Float_t                fWCLAPPDRadius; // Radius of LAPPD
  Int_t                  fWCNumLAPPD;   // Number of LAPPDs
  Float_t                fMRDPMTRadius; // Radius of MRD PMTs
  Int_t                  fWCNumMrdPMT;  // Number of MRD PMTs
  Float_t                fFACCPMTRadius; //Radius of FACC PMTs
  Int_t                  fWCNumFaccPMT;  // Number of FACC PMTs
  Float_t                fWCOffset[3]; // Offset of barrel center in global coords
  Int_t                  fOrientation; //Orientation o detector, 0 is 2km horizontal, 1 is Upright

  // Could make a TClonesArray of PMTs but let's keep it simple
  //   since the arrays just won't be that large
  //WCSimRootPMT          fPMTArray[maxNumPMT];  // Array of PMTs
  TClonesArray           *fPMTArray;
  TClonesArray           *fLAPPDArray;
  TClonesArray           *fMRDPMTArray;
  TClonesArray           *fFACCPMTArray;

public:

  WCSimRootGeom();
  virtual ~WCSimRootGeom();

  // Sets and gets

  void  SetWCCylRadius(Float_t f) {fWCCylRadius=f;}
  void  SetWCCylLength(Float_t f) {fWCCylLength=f;}

  void SetGeo_Type(Int_t f){fgeo_type = f;}

  void  SetWCNumPMT(Int_t i) {fWCNumPMT= i;}
  void  SetWCPMTRadius(Float_t f) {fWCPMTRadius = f;}
  void  SetWCNumLAPPD(Int_t i) {fWCNumLAPPD= i;}
  void  SetWCLAPPDRadius(Float_t f) {fWCLAPPDRadius = f;}
  void  SetWCNumMrdPMT(Int_t i) {fWCNumMrdPMT = i;}
  void  SetMRDPMTRadius(Float_t f) {fMRDPMTRadius = f;}
  void  SetWCNumFaccPMT(Int_t i) {fWCNumFaccPMT = i;}
  void  SetFACCPMTRadius(Float_t f) {fFACCPMTRadius = f;}
  void  SetWCOffset(Float_t x, Float_t y, Float_t z) 
           {fWCOffset[0]=x; fWCOffset[1]=y; fWCOffset[2] = z;}
  void  SetPMT(Int_t i, Int_t tubeno, Int_t cyl_loc, Float_t rot[3], Float_t pos[3], std::string PmtType, bool expand=true);
  void  SetLAPPD(Int_t i, Int_t lappdno, Int_t cyl_loc, Float_t rot[3], Float_t pos[3], std::string PmtType, bool expand=true);
  void  SetOrientation(Int_t o) {fOrientation = o;}

  Float_t GetWCCylRadius() const {return fWCCylRadius;}
  Float_t GetWCCylLength() const {return fWCCylLength;}

  Int_t GetGeo_Type() const {return fgeo_type;}
  

  Int_t GetWCNumPMT() const {return fWCNumPMT;}
  Float_t GetWCPMTRadius() const {return fWCPMTRadius;}
  Int_t GetWCNumLAPPD() const {return fWCNumLAPPD;}
  Float_t GetWCLAPPDRadius() const {return fWCLAPPDRadius;}
  Int_t GetWCNumMRDPMT() const {return fWCNumMrdPMT;}
  Float_t GetMRDPMTRadius() const {return fMRDPMTRadius;}
  Int_t GetWCNumFACCPMT() const {return fWCNumFaccPMT;}
  Float_t GetFACCPMTRadius() const {return fFACCPMTRadius;}
  Float_t GetWCOffset(Int_t i) const {return (i<3) ? fWCOffset[i] : 0.;}
  Int_t GetOrientation() { return fOrientation; }
  //WCSimRootPMT GetPMT(Int_t i){return *(new WCSimRootPMT());}
  WCSimRootPMT GetPMT(Int_t i){return *(WCSimRootPMT*)(*fPMTArray)[i];}
  WCSimRootPMT GetLAPPD(Int_t i){return *(WCSimRootPMT*)(*fLAPPDArray)[i];}
  WCSimRootPMT GetMRDPMT(Int_t i){return *(WCSimRootPMT*)(*fMRDPMTArray)[i];}
  WCSimRootPMT GetFACCPMT(Int_t i){return *(WCSimRootPMT*)(*fFACCPMTArray)[i];} 

  ClassDef(WCSimRootGeom,1)  //WCSimRootEvent structure
};


#endif
