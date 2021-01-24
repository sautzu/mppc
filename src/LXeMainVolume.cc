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

G4double m_posi = 3.0 * mm;
G4double x_posi = 0.0 * mm;
G4double z_posi = 0.0 * mm;
G4double dx = 0.0 * mm; //コリメーターの中心からのずれ

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeMainVolume::LXeMainVolume(G4RotationMatrix *pRot,
							 const G4ThreeVector &tlate,
							 G4LogicalVolume *pMotherLogical,
							 G4bool pMany,
							 G4int pCopyNo,
							 LXeDetectorConstruction *c)
	//Pass info to the G4PVPlacement constructor
	: G4PVPlacement(pRot, tlate,
					//Temp logical volume must be created here
					new G4LogicalVolume(new G4Box("temp", 1, 1, 1),
										G4Material::GetMaterial("Air"),
										"temp", 0, 0, 0),
					"housing", pMotherLogical, pMany, pCopyNo),
	  fConstructor(c)
{
	CopyValues();

	// Logic Size
	// air_housing
	G4double air_x = fScint_x + 2 * m;
	G4double air_y = fScint_y + 2 * m;
	G4double air_z = fScint_z + 2 * m;
	G4double air_in_x = fScint_x + 1 * m;
	G4double air_in_y = fScint_y + 1 * m;
	G4double air_in_z = fScint_z + 1 * m;
	// Collimator
	G4double colli_s_radius = 3.0 * mm;
	G4double colli_l_radius = 1.0 * mm;
	G4double colli_radius = 30.0 * mm;
	G4double colli_width = 80.0 * mm;
	//if(colli)
	//colli_width = 100.0*mm;
	//if(!colli)
	//colli_width = 1.0*mm;
	// Source
	G4double radi_radius = 3.0 * mm;
	G4double radi_width = 1.0 * mm;
	G4double radi_case_radius = 12.5 * mm;
	G4double radi_case_width = 6.0 * mm;
	//Wood
	G4double wood_height = 600.0 * mm;
	G4double wood_width = 15. * mm;
	//Fe
	G4double fe_height = 150.0 * mm;
	G4double fe_width = 2. * mm;
	//Supporter
	G4double sup_width = 2. * mm;
	G4double sup_height = 240. * mm;
	G4double sup_yoko = 150. * mm;

	G4double d = 0.2;
	G4double slide = 70. * mm;
	//	G4double slide = -colli_width+70.*mm;
	// Logic Position

	// Collimator
	G4double colli_posi_x = x_posi;
	G4double colli_posi_y = -colli_width / 2. - d - fScint_y / 2.; //5.5*mm;
	G4double colli_posi_z = z_posi;
	// Wood1
	G4double wood1_posi_x = 0. * mm;
	G4double wood1_posi_y = fScint_y / 2. + wood_width / 2. + 20. * mm;
	G4double wood1_posi_z = 0. * mm;
	// Wood2
	G4double wood2_posi_x = 0. * mm;
	G4double wood2_posi_y = fScint_y / 2. - wood_width / 2. - 255. * mm;
	G4double wood2_posi_z = 0. * mm;
	// Fe1
	G4double Fe1_posi_x = 0. * mm;
	//G4double Fe1_posi_y = -fScint_y / 2. - fe_width / 2. - 0.2 * mm;
	G4double Fe1_posi_y = -fScint_y / 2. - fe_width / 2. - colli_width - d;
	G4double Fe1_posi_z = 0. * mm;
	// Fe2
	G4double Fe2_posi_x = 0. * mm;
	G4double Fe2_posi_y = fScint_y / 2. + fe_width / 2. - 255. * mm;
	G4double Fe2_posi_z = 0. * mm;

	// Source
	G4double source_posi_x = x_posi;
	//G4double source_posi_y = wood1_posi_y + wood_width / 2. + radi_case_width / 2. + slide + 0.01;
	G4double source_posi_y = -colli_width - radi_case_width / 2. - fScint_y / 2. - fe_width - d;
	//	G4double source_posi_y = -colli_width - radi_case_width/2. - fScint_y/2. - fe_width - d/2.-18.*mm;
	//	G4double source_posi_y = wood1_posi_y + wood_width/2. + radi_case_width/2. + 0.01;//-0.02-d -fScint_y/2. -colli_width - radi_case_width/2.;//-colli_width + 5.5*mm + fScint_y/2. + fD_mtl + radi_case_width/2.;
	G4double source_posi_z = z_posi;

	//Supporter
	G4double sup_posi_x = 0. * mm;
	G4double sup_posi_y = Fe2_posi_y + slide + sup_height / 2. * mm;
	G4double sup_posi_z = fe_height / 2. + 50. * mm;

	// center
	G4ThreeVector center_posi(0, 0, 0);

	// Collimator
	G4ThreeVector colli_posi(colli_posi_x, colli_posi_y, colli_posi_z);

	// Wood
	G4ThreeVector wood1_posi(wood1_posi_x, wood1_posi_y + slide, wood1_posi_z);
	G4ThreeVector wood2_posi(wood2_posi_x, wood2_posi_y + slide, wood2_posi_z);

	// Fe
	G4ThreeVector fe1_posi(Fe1_posi_x, Fe1_posi_y, Fe1_posi_z);
	G4ThreeVector fe2_posi(Fe2_posi_x, Fe2_posi_y + slide, Fe2_posi_z);

	// Source
	G4ThreeVector source_posi(source_posi_x, source_posi_y, source_posi_z);
	G4ThreeVector radi_posi(0. * mm, 0. * mm, radi_case_width / 3.);

	// Supporter
	G4ThreeVector sup_posi(sup_posi_x, sup_posi_y, sup_posi_z);

	G4RotationMatrix *rm_x1 = new G4RotationMatrix();
	rm_x1->rotateX(90. * deg);
	G4RotationMatrix *rm_x2 = new G4RotationMatrix();
	rm_x2->rotateX(180. * deg);
	G4RotationMatrix *rm_y1 = new G4RotationMatrix();
	rm_y1->rotateY(90. * deg);
	G4RotationMatrix *rm_y2 = new G4RotationMatrix();
	rm_y2->rotateY(180. * deg);
	G4RotationMatrix *rm_z1 = new G4RotationMatrix();
	rm_z1->rotateZ(90. * deg);
	G4RotationMatrix *rm_z2 = new G4RotationMatrix();
	rm_z2->rotateZ(180. * deg);

	//*************************** housing and scintillator
	fWood_box = new G4Box("wood_box", wood_height / 2., wood_width / 2., wood_height / 2.);
	//  fWood1_box 		= new G4Box("wood_box",wood_width/2.,wood_height/2. - wood_width/2.,wood_height/2. - wood_width/2.); //追加
	//  fWood2_box 		= new G4Box("wood_box",wood_height/2. - wood_width/2.,wood_height/2. - wood_width/2.,wood_width/2.); //追加
	fFe_box = new G4Box("fe_box", fe_height / 2., fe_width / 2. * mm, 120. / 2.);
	fAir_box = new G4Box("air_box", air_x / 2., air_y / 2., air_z / 2.);
	fAir_in_box = new G4Box("air_in_box", air_in_x / 2., air_in_y / 2., air_in_z / 2.);
	fCollimator = new G4Cons("colli_con", colli_s_radius, colli_radius, colli_l_radius, colli_radius, colli_width / 2., 0. * degree, 360. * degree);
	fSup_box = new G4Box("Sup_box", sup_yoko / 2., sup_height / 2., sup_width / 2.);
	fRadiation = new G4Tubs("radi_tub", 0, radi_radius, radi_width / 2., 0. * degree, 360. * degree);
	fRadi_case = new G4Tubs("radi_tub", 0, radi_case_radius, radi_case_width / 2., 0. * degree, 360. * degree);

	fAir_log = new G4LogicalVolume(fAir_box, G4Material::GetMaterial("Air"), "scint_log", 0, 0, 0);
	fAir_in_log = new G4LogicalVolume(fAir_in_box, G4Material::GetMaterial("Air"), "scint_log", 0, 0, 0);
	fWood_log = new G4LogicalVolume(fWood_box, G4Material::GetMaterial("Wood"), "wood_log", 0, 0, 0);
	//	fWood1_log = new G4LogicalVolume(fWood1_box,	G4Material::GetMaterial("Wood"),	"wood_log", 0,0,0);
	//	fWood2_log = new G4LogicalVolume(fWood2_box,	G4Material::GetMaterial("Wood"),	"wood_log", 0,0,0);
	fFe_log = new G4LogicalVolume(fFe_box, G4Material::GetMaterial("Al2"), "fe_log", 0, 0, 0);
	fCollimator_log = new G4LogicalVolume(fCollimator, G4Material::GetMaterial("Pb"), "colli_log", 0, 0, 0);
	fSup_log = new G4LogicalVolume(fSup_box, G4Material::GetMaterial("Al2"), "sup_log", 0, 0, 0);
	fSource_log = new G4LogicalVolume(fRadi_case, G4Material::GetMaterial("PMMA"), "source_log", 0, 0, 0);
	fRadiation_log = new G4LogicalVolume(fRadiation, G4Material::GetMaterial("Air"), "radiation_log", 0, 0, 0);

	new G4PVPlacement(0, radi_posi, fRadiation_log, "radiation", fSource_log, false, 0);

	G4ThreeVector detector_posi(0, 40. * mm - 2.6 * mm, dx);
	detector_vol.reset(new VolumeDetector());
	fDetector_vol = new G4PVPlacement(0, detector_posi, detector_vol->getLogicalVolume(), "Detector", fAir_in_log, false, 0, false);

	//collimator
	//if(colli)
	new G4PVPlacement(rm_x1, colli_posi, fCollimator_log, "collimator", fAir_in_log, false, 0);
	new G4PVPlacement(0, wood1_posi, fWood_log, "Wood", fAir_in_log, false, 0);
	new G4PVPlacement(0, wood2_posi, fWood_log, "wood", fAir_in_log, false, 0);

	new G4PVPlacement(0, fe1_posi, fFe_log, "fe", fAir_in_log, false, 0); //no colli
	new G4PVPlacement(0, fe2_posi, fFe_log, "fe", fAir_in_log, false, 0);

	//Supporter
	////	new G4PVPlacement(0,sup_posi,fSup_log,"supporter",fAir_in_log,false,0);

	//radiation source
	new G4PVPlacement(rm_x1, source_posi, fSource_log, "source", fAir_in_log, false, 0);

	//****************** Build PMTs
	//the "photocathode" is a metal slab at the back of the epoxy that
	//is only a very rough approximation of the real thing since it only
	//absorbs or detects the photons based on the efficiency set below
	fAir_vol = new G4PVPlacement(0, center_posi, fAir_in_log, "airin", fAir_log, false, 8);

	VisAttributes();
	SurfaceProperties();

	SetLogicalVolume(fAir_log);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeMainVolume::CopyValues()
{
	fScint_x = fConstructor->GetScintX();
	fScint_y = fConstructor->GetScintY();
	fScint_z = fConstructor->GetScintZ();
	fD_mtl = fConstructor->GetHousingThickness();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeMainVolume::VisAttributes()
{
	detector_vol->VisAttributes();
	G4VisAttributes *colli_va = new G4VisAttributes(G4Colour(0.4, 0.4, 0.6));
	fCollimator_log->SetVisAttributes(colli_va);
	G4VisAttributes *radi_va = new G4VisAttributes(G4Colour(0.8, 0.8, 0.8));
	fSource_log->SetVisAttributes(radi_va);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//G4OpticalSurface("name",unified/glisur,polished/ground,dielectric&metal)
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
//diele -> metal : 光子は透過せず金属に吸収されるか不導体の方へ反射される
//diele -> diele : 光子は通過するか反射されるか吸収される

void LXeMainVolume::SurfaceProperties()
{
	G4double ephoton[] = {2.76 * eV, 3.06 * eV};
	const G4int num = sizeof(ephoton) / sizeof(G4double);

	//**GSOの反射設定 (GSO周り)
	//	G4double GSO_RIND[] = {1.0,1.0};
	//	G4double GSO_REF[] = {1.0,1.0};
	//	G4double GSO2_RIND[] = {1.0,1.0};
	//	G4double GSO2_REF[] = {1.0,1.0};
	//	G4MaterialPropertiesTable* GSO_PT = new G4MaterialPropertiesTable();
	//	G4MaterialPropertiesTable* GSO2_PT = new G4MaterialPropertiesTable();
	//	GSO_PT->AddProperty("RINDEX",ephoton,GSO_RIND,num);
	//	GSO_PT->AddProperty("REFLECTIVITY",ephoton,GSO_REF,num);
	//	GSO2_PT->AddProperty("RINDEX",ephoton,GSO2_RIND,num);
	//	GSO2_PT->AddProperty("REFLECTIVITY",ephoton,GSO2_REF,num);
	//  G4OpticalSurface* GSO_Surface =
	//    new G4OpticalSurface("GSOSurface",glisur,groundbackpainted,dielectric_dielectric);
	//	GSO_Surface->SetMaterialPropertiesTable(GSO_PT);
	//  G4OpticalSurface* GSO2_Surface = new G4OpticalSurface("GSO2Surface",glisur,polished,dielectric_dielectric);
	//	GSO2_Surface->SetMaterialPropertiesTable(GSO2_PT);
	//
	//  for(int i=0;i<scint_num;i++){
	//  new G4LogicalBorderSurface("GSO_surf",fScint_vol[i],fAir_vol,GSO_Surface);
	//  }
	//  new G4LogicalBorderSurface("GSO_surf",fScintbig_vol,fAir_vol,GSO_Surface);

	//**CsIの反射設定 (CsI周り)
	//	G4double CsI_RIND[] = {1.0,1.0};
	//	G4double CsI_REF[] = {1.0,1.0};
	//	G4double CsI2_RIND[] = {1.0,1.0};
	//	G4double CsI2_REF[] = {1.0,1.0};
	//
	//	G4MaterialPropertiesTable* CsI_PT = new G4MaterialPropertiesTable();
	//	G4MaterialPropertiesTable* CsI2_PT = new G4MaterialPropertiesTable();
	//
	//	CsI_PT->AddProperty("RINDEX",ephoton,CsI_RIND,num);
	//	CsI_PT->AddProperty("REFLECTIVITY",ephoton,CsI_REF,num);
	//	CsI2_PT->AddProperty("RINDEX",ephoton,CsI2_RIND,num);
	//	CsI2_PT->AddProperty("REFLECTIVITY",ephoton,CsI2_REF,num);
	//
	//  G4OpticalSurface* CsI_Surface = new G4OpticalSurface("CsISurface",glisur,groundbackpainted,dielectric_dielectric);
	//	CsI_Surface->SetMaterialPropertiesTable(CsI_PT);
	//  G4OpticalSurface* CsI2_Surface = new G4OpticalSurface("CsI2Surface",glisur,polished,dielectric_dielectric);
	//	CsI2_Surface->SetMaterialPropertiesTable(CsI2_PT);
	//
	//  for(int i=0;i<scint_num;i++){
	//  new G4LogicalBorderSurface("CsI_surf",fScint_vol[i],fAir_vol,CsI_Surface);
	//	}
	//
	//  new G4LogicalBorderSurface("CsI_surf",fScintbig_vol,fAir_vol,CsI_Surface);

	detector_vol->SurfaceProperties(fAir_vol);
}
