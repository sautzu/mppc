//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: LXeMainVolume.hh 77486 2013-11-25 10:14:16Z gcosmo $
//
/// \file optical/LXe/include/LXeMainVolume.hh
/// \brief Definition of the LXeMainVolume class
//
#ifndef LXeMainVolume_H
#define LXeMainVolume_H 1

#include "G4PVPlacement.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4SubtractionSolid.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4OpticalSurface.hh"
#include "G4VPhysicalVolume.hh"

#include "LXeDetectorConstruction.hh"

class LXeMainVolume : public G4PVPlacement
{
  public:

    LXeMainVolume(G4RotationMatrix *pRot,
                 const G4ThreeVector &tlate,
                 G4LogicalVolume *pMotherLogical,
                 G4bool pMany,
                 G4int pCopyNo,
                 LXeDetectorConstruction* c);

    //G4LogicalVolume* GetLogPhotoCath() {return fPmt_log;}
    G4LogicalVolume* GetLogScint()     {return fScint_log;}
    G4LogicalVolume* GetLogScintbig()     {return fScintbig_log;}
    G4LogicalVolume* GetLogMPPC()     {return fMPPC_log;}

    std::vector<G4ThreeVector> GetPmtPositions() {return fPmtPositions;}

  private:

    void VisAttributes();
    void SurfaceProperties();

    void PlacePMTs(G4LogicalVolume* pmt_Log,
                   G4RotationMatrix* rot,
                   G4double &a, G4double &b, G4double da,
                   G4double db, G4double amin, G4double bmin,
                   G4int na, G4int nb,
                   G4double &x, G4double &y, G4double &z, G4int &k);

    void CopyValues();

    LXeDetectorConstruction* fConstructor;

    G4double fScint_x;
    G4double fScint_y;
    G4double fScint_z;
    G4double fD_mtl;
		G4int fNx;
		G4int fNy;
		G4int fNz;
    G4double fOuterRadius_pmt;
    G4bool fSphereOn;
    G4double fRefl;

    //Basic Volumes
    G4Box* fScint_box;
    G4Box* fScintbig_box;
    G4Box* fWood_box;
//    G4Box* fWood1_box;
//    G4Box* fWood2_box;
    G4Box* fFe_box;
    G4Box* fABS1_box;
    G4Box* fABS2_box;
    G4Box* fGlice_box;
    G4Box* fGlicebig_box;
    G4Box* fGlice2_box;
    G4Box* fGlice3_box;
    G4Box* fAir_box;
    G4Box* fAir_in_box;
    G4Box* fMPPC_box;
    G4Box* fWin_box;
    G4Box* fEnve_box;
//    G4Box* fTef_box;
    G4Box* fSup_box;
    G4Box* fMPPC_out_pre;
		G4SubtractionSolid* fMPPC_out;
		G4Box* fWindow_box;
		G4Cons* fCollimator;
		G4Tubs* fRadiation;
		G4Tubs* fRadi_case;

    // Logical volumes
    G4LogicalVolume* fScint_log;
    G4LogicalVolume* fScintbig_log;
    G4LogicalVolume* fWood_log;
//    G4LogicalVolume* fWood1_log;
//    G4LogicalVolume* fWood2_log;
    G4LogicalVolume* fFe_log;
    G4LogicalVolume* fABS1_log;
    G4LogicalVolume* fABS2_log;
    G4LogicalVolume* fGlice_log;
    G4LogicalVolume* fGlicebig_log;
    G4LogicalVolume* fGlice2_log;
    G4LogicalVolume* fGlice3_log;
    G4LogicalVolume* fAir_log;
    G4LogicalVolume* fAir_in_log;
    G4LogicalVolume* fMPPC_log;
    G4LogicalVolume* fWin_log;
    G4LogicalVolume* fEnve_log;
//    G4LogicalVolume* fTef_log;
    G4LogicalVolume* fSup_log;
    G4LogicalVolume* fOut_log;
    G4LogicalVolume* fWindow_log;
    G4LogicalVolume* fCollimator_log;
    G4LogicalVolume* fSource_log;
    G4LogicalVolume* fRadiation_log;

    // Physical volumes
		G4VPhysicalVolume* fScint_vol[3];
		G4VPhysicalVolume* fScintbig_vol;
    G4VPhysicalVolume *fGuide_vol[2];
    G4VPhysicalVolume* fGlice_vol;
		G4VPhysicalVolume* fGlicebig_vol;
		G4VPhysicalVolume* fGlice2_vol;
		G4VPhysicalVolume* fGlice3_vol;
		G4VPhysicalVolume* fMPPC_vol;
		G4VPhysicalVolume* fWin_vol[3];
		G4VPhysicalVolume* fEnve_vol[3];
	        G4VPhysicalVolume* fSup_vol;
	        G4VPhysicalVolume* fWindow1_vol;
		G4VPhysicalVolume* fWindow2_vol;
		G4VPhysicalVolume* fWindow3_vol;
		G4VPhysicalVolume* fWindow4_vol;
		G4VPhysicalVolume* fWindow5_vol;
		G4VPhysicalVolume* fWindow6_vol;
		G4VPhysicalVolume* fWindow7_vol;
		G4VPhysicalVolume* fWindow8_vol;
		G4VPhysicalVolume* fOut1_vol;
		G4VPhysicalVolume* fOut2_vol;
		G4VPhysicalVolume* fOut3_vol;
		G4VPhysicalVolume* fOut4_vol;
		G4VPhysicalVolume* fOut5_vol;
		G4VPhysicalVolume* fOut6_vol;
		G4VPhysicalVolume* fOut7_vol;
		G4VPhysicalVolume* fOut8_vol;
		G4VPhysicalVolume* fAir_vol;

    // Sensitive Detectors positions
    std::vector<G4ThreeVector> fPmtPositions;

};

#endif
