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
// $Id: LXeMainVolume.cc 82853 2014-07-14 09:07:11Z gcosmo $
//
/// \file optical/LXe/src/LXeMainVolume.cc
/// \brief Implementation of the LXeMainVolume class
//
//
#include "globals.hh"
#include "LXeMainVolume.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
LXeMainVolume::LXeMainVolume(G4RotationMatrix *pRot,
                             const G4ThreeVector &tlate,
                             G4LogicalVolume *pMotherLogical,
                             G4bool pMany,
                             G4int pCopyNo,
                             LXeDetectorConstruction* c)
  //Pass info to the G4PVPlacement constructor
  :G4PVPlacement(pRot,tlate,
                 //Temp logical volume must be created here
                 new G4LogicalVolume(new G4Box("temp",1,1,1),
                                     G4Material::GetMaterial("G4_Galactic"),
                                     "temp",0,0,0),
                 "housing",pMotherLogical,pMany,pCopyNo),fConstructor(c), fMaterials(NULL)
{
	fMaterials = LXeDetectorConstruction::GetInstance();

  CopyValues();

  G4double mppc1_x,mppc1_y,mppc1_z;
  G4double mppc2_x,mppc2_y,mppc2_z;
  G4RotationMatrix* rm_y2 = new G4RotationMatrix();
 
  width_pmt = 3.0*mm;
  height_pmt = 0.4*mm;
  air_x=fScint_x+fD_mtl;
  air_y=fScint_y+fD_mtl;
  air_z=fScint_z+fD_mtl;
  housing_x=fScint_x+2*height_pmt+fD_mtl;
  housing_y=fScint_y+2*height_pmt+fD_mtl;
  housing_z=fScint_z+2*height_pmt+fD_mtl;

	mppc1_x=0;
	mppc1_y=0;
	mppc1_z=-fScint_z/2.-height_pmt/2.;
	mppc2_x=-fScint_x/2.-height_pmt/2.;
	mppc2_y=0;
	mppc2_z=fScint_z/2.-width_pmt/2.+height_pmt/2.;
  rm_y2->rotateY(-90*deg);

  //*************************** housing and scintillator

  fScint_box 			= new G4Box("scint_box",								fScint_x/2.,	fScint_y/2.,	fScint_z/2.);
  fAir_box 				= new G4Box("air_box",									air_x/2.,			air_y/2.,			air_z/2.);
  fPmt 						= new G4Box("pmt_box",									width_pmt/2.,	width_pmt/2.,	height_pmt/2.);
  fPhotocath 			= new G4Box("photocath_box",						width_pmt/2.,	width_pmt/2.,	height_pmt/4.);
	fPre_solid 			= new G4UnionSolid("Pre_solid",fAir_box,fPmt,0,G4ThreeVector(mppc1_x,mppc1_y,mppc1_z));
 	f_solid 				= new G4UnionSolid("solid",fPre_solid,fPmt,rm_y2,G4ThreeVector(mppc2_x,mppc2_y,mppc2_z));

  fScint_log 			= new G4LogicalVolume(fScint_box,			FindMaterial("NaI"),"scint_log",0,0,0);
	fAir_log 				= new G4LogicalVolume(f_solid,		FindMaterial("G4_AIR"),"air_log",0,0,0);
  fPmt_log 				= new G4LogicalVolume(fPmt,						FindMaterial("Silicone"),"pmt_log");
  fPhotocath_log 	= new G4LogicalVolume(fPhotocath,			FindMaterial("G4_Al"),"photocath_log");
 
  new G4PVPlacement(0,G4ThreeVector(),fScint_log,"scintillator",fAir_log,false,0);

 //****************** Build PMTs
 
  //the "photocathode" is a metal slab at the back of the glass that
  //is only a very rough approximation of the real thing since it only
  //absorbs or detects the photons based on the efficiency set below
 
	new G4PVPlacement(0,G4ThreeVector(0,0,-height_pmt/4.),fPhotocath_log,"photocath",fPmt_log,false,0);

  //***********Arrange pmts around the outside of housing**********
  new G4PVPlacement(0,G4ThreeVector(mppc1_x,mppc1_y,mppc1_z),fPmt_log,"pmt",fAir_log,false,0);//front1//default
      fPmtPositions.push_back(G4ThreeVector(mppc1_x,mppc1_y,mppc1_z));
  new G4PVPlacement(rm_y2,G4ThreeVector(mppc2_x,mppc2_y,mppc2_z),fPmt_log,"pmt",fAir_log,false,1);
      fPmtPositions.push_back(G4ThreeVector(mppc2_x,mppc2_y,mppc2_z));
 
  VisAttributes();
  SurfaceProperties();
	delete rm_y2;

  SetLogicalVolume(fAir_log);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void LXeMainVolume::CopyValues(){
  fScint_x=fConstructor->GetScintX();
  fScint_y=fConstructor->GetScintY();
  fScint_z=fConstructor->GetScintZ();
  fD_mtl=fConstructor->GetHousingThickness();
  fNx=fConstructor->GetNX();
  fNy=fConstructor->GetNY();
  fNz=fConstructor->GetNZ();
  fOuterRadius_pmt=fConstructor->GetPMTRadius();
  fRefl=fConstructor->GetHousingReflectivity();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void LXeMainVolume::VisAttributes(){
  G4VisAttributes* scint_va = new G4VisAttributes(G4Colour(0.8,0.8,0.2));
  fScint_log->SetVisAttributes(scint_va);

  G4VisAttributes* pmt_va = new G4VisAttributes(G4Colour(0.2,0.4,0.2));
  G4VisAttributes* photo_va = new G4VisAttributes(G4Colour(0.4,0.8,0.4));
  fPmt_log->SetVisAttributes(pmt_va);
  fPhotocath_log->SetVisAttributes(photo_va);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  //G4OpticalSurface("name",unified/glisur,polished/ground,dielectric&    metal)
  //glisur : polish = 1.0(default)で鏡面反射
  //glisur : polish = 0.0でランダム反射
  //unified : polish~ +alphaパラメータで散乱具合を変更
  //unified : ground~ でランバート反射
  //polished : 鏡面反射
  //polishedfrontpainted : 前面塗料
  //polishedbackpaint : 背面塗料
  //ground : 荒い
  //groundfrontpainted : 荒い前面
  //groundbackpainted : 荒い後面
  //diele -> metal : 光子は透過せず金属に吸収されるか不導体の方へ反射>    される
  //diele -> diele : 光子は通過するか反射されるか吸収される
void LXeMainVolume::SurfaceProperties(){
  G4double ephoton[] = {2.76*eV, 3.06*eV};
  const G4int num = sizeof(ephoton)/sizeof(G4double);
  //**GSOの反射設定 (GSO周り)
  //反射率
  G4OpticalSurface* GSO_Surface =
    new G4OpticalSurface("GSOSurface",glisur,polished,dielectric_dielectric);
  //**空気層の周り
  G4double Air_REF[] = {0.9, 0.9};
  assert(sizeof(Air_REF) == sizeof(ephoton));
  //光電子効率？
  G4double Air_EFF[] = {0.0, 0.0};
  assert(sizeof(Air_EFF) == sizeof(ephoton));
  G4MaterialPropertiesTable* Air_PT = new G4MaterialPropertiesTable();
  Air_PT->AddProperty("REFLECTIVITY", ephoton, Air_REF, num);
  Air_PT->AddProperty("EFFICIENCY", ephoton, Air_EFF, num);
  G4OpticalSurface* Air_Surface =
    new G4OpticalSurface("AirSurface",glisur,groundfrontpainted,dielectric_dielectric);
  Air_Surface->SetMaterialPropertiesTable(Air_PT);
 
  //**mppcの反射設定 
  G4double pmt_REF[] = {0.8,0.8};
  assert(sizeof(pmt_REF) == sizeof(ephoton));
  G4double pmt_EFF[]={0.0,0.0}; //enables 'detection' of photons
  assert(sizeof(pmt_EFF) == sizeof(ephoton));
  G4MaterialPropertiesTable* pmt_mt = new G4MaterialPropertiesTable();
  pmt_mt->AddProperty("REFLECTIVITY", ephoton, pmt_REF, num);
  pmt_mt->AddProperty("EFFICIENCY",ephoton,pmt_EFF,num);
  G4OpticalSurface* pmt_opsurf=
    new G4OpticalSurface("pmt_opsurf",glisur,polished,
                         dielectric_dielectric);
  pmt_opsurf->SetMaterialPropertiesTable(pmt_mt);
  //**mppcの反射設定(光子感度) 
  G4double photocath_ReR[]={1.92,1.92};
  assert(sizeof(photocath_ReR) == sizeof(ephoton));
  G4double photocath_ImR[]={1.69,1.69};
  assert(sizeof(photocath_ImR) == sizeof(ephoton));
  G4double photocath_EFF[]={0.4,0.4}; //enables 'detection' of photons
  assert(sizeof(photocath_EFF) == sizeof(ephoton));
  G4MaterialPropertiesTable* photocath_mt = new G4MaterialPropertiesTable();
  photocath_mt->AddProperty("REALRINDEX",ephoton,photocath_ReR,num);
  photocath_mt->AddProperty("IMAGINARYRINDEX",ephoton,photocath_ImR,num);
  photocath_mt->AddProperty("EFFICIENCY",ephoton,photocath_EFF,num);
  G4OpticalSurface* photocath_opsurf=
    new G4OpticalSurface("photocath_opsurf",glisur,polished,
                         dielectric_metal);
  photocath_opsurf->SetMaterialPropertiesTable(photocath_mt);

  //**Create logical skin surfaces
  new G4LogicalSkinSurface("GSO_surf",fScint_log,GSO_Surface);
  new G4LogicalSkinSurface("Air_surf",fAir_log,
                           Air_Surface);
  new G4LogicalSkinSurface("photocath_surf",fPhotocath_log,photocath_opsurf);
  new G4LogicalSkinSurface("pmt_surf",fPmt_log,pmt_opsurf);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*void LXeMainVolume::SurfaceProperties(){
  G4double ephoton[] = {7.0*eV, 7.14*eV};
  const G4int num = sizeof(ephoton)/sizeof(G4double);

  //Scintillator housing properties
  G4double reflectivity[] = {fRefl, fRefl};
  assert(sizeof(reflectivity) == sizeof(ephoton));
  G4double efficiency[] = {0.0, 0.0};
  assert(sizeof(efficiency) == sizeof(ephoton));
  G4MaterialPropertiesTable* scintHsngPT = new G4MaterialPropertiesTable();
  scintHsngPT->AddProperty("REFLECTIVITY", ephoton, reflectivity, num);
  scintHsngPT->AddProperty("EFFICIENCY", ephoton, efficiency, num);
  G4OpticalSurface* OpScintHousingSurface =
    new G4OpticalSurface("HousingSurface",unified,polished,dielectric_metal);
  OpScintHousingSurface->SetMaterialPropertiesTable(scintHsngPT);
 
  //Photocathode surface properties
  G4double photocath_EFF[]={1.,1.}; //Enables 'detection' of photons
  assert(sizeof(photocath_EFF) == sizeof(ephoton));
  G4double photocath_ReR[]={1.92,1.92};
  assert(sizeof(photocath_ReR) == sizeof(ephoton));
  G4double photocath_ImR[]={1.69,1.69};
  assert(sizeof(photocath_ImR) == sizeof(ephoton));
  G4MaterialPropertiesTable* photocath_mt = new G4MaterialPropertiesTable();
  photocath_mt->AddProperty("EFFICIENCY",ephoton,photocath_EFF,num);
  photocath_mt->AddProperty("REALRINDEX",ephoton,photocath_ReR,num);
  photocath_mt->AddProperty("IMAGINARYRINDEX",ephoton,photocath_ImR,num);
  G4OpticalSurface* photocath_opsurf=
    new G4OpticalSurface("photocath_opsurf",glisur,polished,
                         dielectric_metal);
  photocath_opsurf->SetMaterialPropertiesTable(photocath_mt);

  //Create logical skin surfaces
  new G4LogicalSkinSurface("photocath_surf",fHousing_log,
                           OpScintHousingSurface);
  new G4LogicalSkinSurface("photocath_surf",fPhotocath_log,photocath_opsurf);
}*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4Material* LXeMainVolume::FindMaterial(G4String name)
{
	G4Material* material = G4Material::GetMaterial(name,true);
	return material;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
