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
// $Id: LXeDetectorConstruction.cc 101905 2016-12-07 11:34:39Z gunter $
//
/// \file optical/LXe/src/LXeDetectorConstruction.cc
/// \brief Implementation of the LXeDetectorConstruction class
//
//
#include "LXeDetectorConstruction.hh"
#include "LXePMTSD.hh"
#include "LXeScintSD.hh"
#include "LXeDetectorMessenger.hh"
#include "LXeMainVolume.hh"

#include "G4SDManager.hh"
#include "G4RunManager.hh"

#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"

#include "G4OpticalSurface.hh"
#include "G4MaterialTable.hh"
#include "G4VisAttributes.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4UImanager.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeDetectorConstruction::LXeDetectorConstruction()
: gso_mt(NULL), gagg_mt(NULL), csi_mt(NULL), vacuum_mt(NULL), air_mt(NULL), ptfe_mt(NULL),
pmma_mt(NULL), silicone_mt(NULL), al_mt(NULL)
{
  fExperimentalHall_box = NULL;
  fExperimentalHall_log = NULL;
  fExperimentalHall_phys = NULL;

  fAlumi = fAlumi2 = fNaI = fAir = fVacuum = fEpoxy = fSilicone = fGSO = fPTFE = fABS = fSilicone = fSilica = fPb = fCa = fGAGG = fCsI = NULL;

  fN = fO = fC = fH = fGd = fSi = fF = fNa = fI = fCs = fGa = fAl =  NULL;

  SetDefaults();

  fDetectorMessenger = new LXeDetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeDetectorConstruction::~LXeDetectorConstruction() {
	delete fTetu; 			delete fAlumi;		delete fAlumi2;					delete fNaI;					delete fAir;
	delete fVacuum;			delete fEpoxy;				delete fGSO;		delete fABS;
	delete fPTFE;				delete fSilicone;		delete fSilica;		delete fN;
	delete fGAGG;	delete fCsI;	delete fWood;
	delete fO;					delete fC;						delete fH;
	delete fGd;					delete fSi;						delete fF;
	delete fNa;					delete fI;						delete fCs;
	delete fGa;			delete fAl;
	delete gso_mt;			delete fe_mt;	delete wood_mt; 	delete gagg_mt;		delete csi_mt; 	delete vacuum_mt;			delete air_mt;
	delete ptfe_mt;			delete pmma_mt;				delete silicone_mt;
	delete al_mt;			delete fPb;	delete fCa;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::DefineMaterials(){
	//***definition of G4int***//
  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density;

  //***Elements***//
  fH = new G4Element("H", "H", z=1., a=1.01*g/mole);
  fC = new G4Element("C", "C", z=6., a=12.01*g/mole);
  fN = new G4Element("N", "N", z=7., a= 14.01*g/mole);
  fO = new G4Element("O", "O", z=8., a= 16.00*g/mole);
  fGd = new G4Element("Gd", "Gd", z=64., a= 157.25*g/mole);
  fSi = new G4Element("Si", "Si", z=14., a= 28.09*g/mole);
  fNa = new G4Element("Na", "Na", z=11., a= 22.99*g/mole);
  fI = new G4Element("I", "I", z=53., a= 126.90*g/mole);
  fCs = new G4Element("Cs", "Cs", z=55., a= 132.91*g/mole);
  fF = new G4Element("F", "F", z=9., a= 19.00*g/mole);
  fGa = new G4Element("Ga", "Ga", z=31., a= 69.70*g/mole);
  fAl = new G4Element("Al", "Al", z=13., a= 26.98*g/mole);

  //***Materials***//
  //Vacuum
  fVacuum = new G4Material("Vacuum",z=1.,a=1.01*g/mole,
                          density=universe_mean_density,kStateGas,0.1*kelvin,
                          1.e-19*pascal);
  //Air
  fAir = new G4Material("Air", density= 1.29*mg/cm3, 2);
  fAir->AddElement(fN, 78*perCent);
  fAir->AddElement(fO, 22*perCent);
  //NaI
  fNaI = new G4Material("NaI", density= 3.67*g/cm3, 2);
  fNaI->AddElement(fNa, 1);
  fNaI->AddElement(fI, 1);
  //Almimnium
  fAlumi = new G4Material("Al",z=13.,a=26.98*g/mole,density=2.7*g/cm3);
  fAlumi2 = new G4Material("Al2",z=13.,a=26.98*g/mole,density=2.7*g/cm3);
  fTetu = new G4Material("Fe",z=26.,a=55.845*g/mole,density=7.874*g/cm3);
  //Lead
  fPb = new G4Material("Pb",z=82.,a=207.19*g/mole,density=11.35*g/cm3);
  fCa = new G4Material("Ca",z=20.,a=40.0*g/mole,density=1.550*g/cm3);
  //Epoxy
  fEpoxy = new G4Material("Epoxy", density=1.032*g/cm3,2);
  fEpoxy->AddElement(fC,91.533*perCent);
  fEpoxy->AddElement(fH,8.467*perCent);
	//GSO
	fGSO = new G4Material("GSO", density=6.71*g/cm3,3);
	fGSO->AddElement(fGd, 2);
	fGSO->AddElement(fSi, 1);
	fGSO->AddElement(fO, 5);
	//PTFE
	fPTFE = new G4Material("PTFE", density=2.20*g/cm3, 2);
	fPTFE->AddElement(fC, 2);
	fPTFE->AddElement(fF, 4);
	//PMMA
	fPMMA = new G4Material("PMMA", density=1.18*g/cm3, 3);
	fPMMA->AddElement(fC, 5);
	fPMMA->AddElement(fO, 2);
	fPMMA->AddElement(fH, 8);
	//Silicone
	fSilicone = new G4Material("Silicone", density=1.06*g/cm3, 2);
	fSilicone->AddElement(fC, 2);
	fSilicone->AddElement(fH, 6);
	//Silica
	fSilica = new G4Material("Silica", density=1.06*g/cm3, 2);
	fSilica->AddElement(fSi, 1);
	fSilica->AddElement(fO, 2);
	//GAGG
	fGAGG = new G4Material("GAGG", density=6.63*g/cm3, 4);
	//fGAGG->AddElement(fGd, 1);
	fGAGG->AddElement(fGd, 3);
	fGAGG->AddElement(fAl, 2);
	fGAGG->AddElement(fGa, 3);
	fGAGG->AddElement(fO, 12);
	//CsI
	fCsI = new G4Material("CsI", density=4.51*g/cm3, 2);
	fCsI->AddElement(fCs, 1);
	fCsI->AddElement(fI,1);
	//WOOD
	fWood = new G4Material("Wood", density=1.0*g/cm3, 3);
	fWood->AddElement(fC, 6);
	fWood->AddElement(fO, 5);
	fWood->AddElement(fH, 10);
	//ABS
	fABS = new G4Material("ABS",density=1.06*g/cm3,3);
        fABS->AddElement(fC,15);
        fABS->AddElement(fH,17);
        fABS->AddElement(fN,1);
//        //TEFLON
//        fTef = new G4Material("Teflon",density=2.13*g/cm3,2);
//        fTef->AddElement(fC,1);
//        fTef->AddElement(fF,1);

  //***Material properties tables***//

	G4double photon_energy[] = 
	{ 2.08*eV, 2.18*eV, 2.28*eV, 2.38*eV, 2.51*eV, 2.64*eV, 2.79*eV, 2.95*eV, 3.14*eV, 3.35*eV };
  const G4int num = sizeof(photon_energy)/sizeof(G4double);
	
	G4double epoxy_RIND[10];
	G4double ptfe_RIND[10];
	G4double pmma_RIND[10];
	G4double silicone_RIND[10];
	G4double silica_RIND[10];
	G4double gso_RIND[10];
	G4double gagg_RIND[10];
	G4double csi_RIND[10];
	G4double al_RIND[10];
	G4double vacuum_RIND[10];
	G4double air_RIND[10];
	G4double gso_ABSL[10];
	G4double gagg_ABSL[10];
	G4double csi_ABSL[10];
	G4double ptfe_ABSL[10];
	G4double epoxy_ABSL[10];
	G4double vacuum_ABSL[10];
	G4double air_ABSL[10];

	for(int i=0;i<10;i++){
		epoxy_RIND[i] = 1.453;
		ptfe_RIND[i] = 1.35;
		pmma_RIND[i] = 1.50;
		silicone_RIND[i] = 1.453;
		silica_RIND[i] = 1.453;
		gso_RIND[i] = 1.85;
		gagg_RIND[i] = 1.93;
//		gagg_RIND[i] = 1.0;
		csi_RIND[i] = 1.79;
		al_RIND[i] = 1.48;
		vacuum_RIND[i] = 1.000293;
		air_RIND[i] = 1.0003;
		gso_ABSL[i] = 0.5*m;
		gagg_ABSL[i] = 0.5*m;
		csi_ABSL[i] = 0.5*m;
		ptfe_ABSL[i] = 20.0*cm;
		epoxy_ABSL[i] = 20.0*cm;
		air_ABSL[i] = 1.*m;
		vacuum_ABSL[i] = 1.*m;
	}
	
	G4double gso_SCINT[] = { 0.01,0.01,0.10,0.20,0.45,0.70,0.95,0.90,0.50,0.25 };
	G4double gagg_SCINT[] = { 0.40,0.70,0.90,1.00,0.50,0.25,0.10,0.01,0.01,0.01 };
	G4double csi_SCINT[] = { 0.80,1.00,0.95,0.60,0.30,0.18,0.15,0.10,0.05,0.01 };

	//***EPOXY***
  epoxy_mt = new G4MaterialPropertiesTable();
  epoxy_mt->AddProperty("ABSLENGTH",photon_energy,epoxy_ABSL,num);
  epoxy_mt->AddProperty("RINDEX",photon_energy,epoxy_RIND,num);
  fEpoxy->SetMaterialPropertiesTable(epoxy_mt);

	//***PTFE***
	//屈折率
  ptfe_mt = new G4MaterialPropertiesTable();
  ptfe_mt->AddProperty("ABSLENGTH",photon_energy,ptfe_ABSL,num);
  ptfe_mt->AddProperty("RINDEX",photon_energy,ptfe_RIND,num);
  fPTFE->SetMaterialPropertiesTable(ptfe_mt);

	//***PMMA***
	//屈折率
 // pmma_mt = new G4MaterialPropertiesTable();
 // pmma_mt->AddProperty("ABSLENGTH",photon_energy,ptfe_ABSL,num);
 // pmma_mt->AddProperty("RINDEX",photon_energy,pmma_RIND,num);
 // fPMMA->SetMaterialPropertiesTable(pmma_mt);

	//***SILICONE***
	//屈折率
  silicone_mt = new G4MaterialPropertiesTable();
  silicone_mt->AddProperty("ABSLENGTH",photon_energy,ptfe_ABSL,num);
  silicone_mt->AddProperty("RINDEX",photon_energy,silicone_RIND,num);
  fSilicone->SetMaterialPropertiesTable(silicone_mt);
  
	//***SILICA***
	//屈折率
  silica_mt = new G4MaterialPropertiesTable();
  silica_mt->AddProperty("ABSLENGTH",photon_energy,ptfe_ABSL,num);
  silica_mt->AddProperty("RINDEX",photon_energy,silica_RIND,num);
  fSilica->SetMaterialPropertiesTable(silica_mt);

	 //***GSO***
  gso_mt = new G4MaterialPropertiesTable();
  gso_mt->AddProperty("ABSLENGTH",photon_energy,gso_ABSL,num);
  gso_mt->AddProperty("RINDEX",photon_energy,gso_RIND,num);
	gso_mt->AddProperty("FASTCOMPONENT",photon_energy, gso_SCINT, num);
  gso_mt->AddProperty("SLOWCOMPONENT", photon_energy, gso_SCINT, num);
  //gso_mt->AddProperty("SCINTILLATION", photon_energy, gso_SCINT, num);
	gso_mt->AddConstProperty("SCINTILLATIONYIELD", 10000./MeV );
  gso_mt->AddConstProperty("RESOLUTIONSCALE", 0.01 );//check
  gso_mt->AddConstProperty("FASTTIMECONSTANT", 30.*ns );
  gso_mt->AddConstProperty("SLOWTIMECONSTANT", 600.*ns );
  gso_mt->AddConstProperty("YIELDRATIO",0.9);
  fGSO->SetMaterialPropertiesTable(gso_mt);
  fGSO->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

	//***GAGG***
  gagg_mt = new G4MaterialPropertiesTable();
  gagg_mt->AddProperty("ABSLENGTH",photon_energy,gagg_ABSL,num);
  gagg_mt->AddProperty("RINDEX",photon_energy,gagg_RIND,num);
	gagg_mt->AddProperty("FASTCOMPONENT",photon_energy, gagg_SCINT, num);
  gagg_mt->AddProperty("SLOWCOMPONENT", photon_energy, gagg_SCINT, num);
  //gagg_mt->AddProperty("SCINTILLATION", photon_energy, gagg_SCINT, num);
	gagg_mt->AddConstProperty("SCINTILLATIONYIELD", 65000./MeV);//65000./MeV );
//	gagg_mt->AddConstProperty("SCINTILLATIONYIELD", 0./MeV);//65000./MeV );
  gagg_mt->AddConstProperty("RESOLUTIONSCALE", 8 );//check
  gagg_mt->AddConstProperty("FASTTIMECONSTANT", 88.*ns );
  gagg_mt->AddConstProperty("SLOWTIMECONSTANT", 258.*ns );//258
  gagg_mt->AddConstProperty("YIELDRATIO",0.9);
  fGAGG->SetMaterialPropertiesTable(gagg_mt);
  fGAGG->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

	//***CsI***
  csi_mt = new G4MaterialPropertiesTable();
  csi_mt->AddProperty("ABSLENGTH",photon_energy,csi_ABSL,num);
  csi_mt->AddProperty("RINDEX",photon_energy,csi_RIND,num);
	csi_mt->AddProperty("FASTCOMPONENT",photon_energy, csi_SCINT, num);
  csi_mt->AddProperty("SLOWCOMPONENT", photon_energy, csi_SCINT, num);
  //csi_mt->AddProperty("SCINTILLATION", photon_energy, csi_SCINT, num);
	csi_mt->AddConstProperty("SCINTILLATIONYIELD", 80000./MeV );
  csi_mt->AddConstProperty("RESOLUTIONSCALE", 12 );//check
  csi_mt->AddConstProperty("FASTTIMECONSTANT", 1000.*ns );
  csi_mt->AddConstProperty("SLOWTIMECONSTANT", 3000.*ns );//258  //nazo
  csi_mt->AddConstProperty("YIELDRATIO",0.9);
  fCsI->SetMaterialPropertiesTable(csi_mt);
  fCsI->GetIonisation()->SetBirksConstant(0.126*mm/MeV);



 
	//***vacuum***
  vacuum_mt = new G4MaterialPropertiesTable();
  vacuum_mt->AddProperty("ABSLENGTH", photon_energy, vacuum_ABSL,num);
  vacuum_mt->AddProperty("RINDEX", photon_energy, vacuum_RIND,num);
  fVacuum->SetMaterialPropertiesTable(vacuum_mt);
	//***air***
  air_mt = new G4MaterialPropertiesTable();
  air_mt->AddProperty("ABSLENGTH", photon_energy, air_ABSL,num);
  air_mt->AddProperty("RINDEX", photon_energy, air_RIND,num);
  fAir->SetMaterialPropertiesTable(air_mt);//Give air the same rindex
	//***alminium***
  al_mt = new G4MaterialPropertiesTable();
  al_mt->AddProperty("RINDEX", photon_energy, al_RIND,num);
  fAlumi->SetMaterialPropertiesTable(al_mt);//Give air the same rindex

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* LXeDetectorConstruction::Construct(){

  if (fExperimentalHall_phys) {
     G4GeometryManager::GetInstance()->OpenGeometry();
     G4PhysicalVolumeStore::GetInstance()->Clean();
     G4LogicalVolumeStore::GetInstance()->Clean();
     G4SolidStore::GetInstance()->Clean();
     G4LogicalSkinSurface::CleanSurfaceTable();
     G4LogicalBorderSurface::CleanSurfaceTable();
  }

  DefineMaterials();
  return ConstructDetector();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* LXeDetectorConstruction::ConstructDetector()
{
  //The experimental hall walls are all 1m away from housing walls
  G4double expHall_x = fScint_x+fD_mtl+1.*m;
  G4double expHall_y = fScint_y+fD_mtl+1.*m;
  G4double expHall_z = fScint_z+fD_mtl+1.*m;

  //Create experimental hall
  fExperimentalHall_box
    = new G4Box("expHall_box",expHall_x,expHall_y,expHall_z);
  fExperimentalHall_log = new G4LogicalVolume(fExperimentalHall_box,
                                             fAir,"expHall_log",0,0,0);
  fExperimentalHall_phys = new G4PVPlacement(0,G4ThreeVector(),
                              fExperimentalHall_log,"expHall",0,false,0);

  fExperimentalHall_log->SetVisAttributes(G4VisAttributes::GetInvisible());

  //Place the main volume
  if(fMainVolumeOn){
    fMainVolume
      = new LXeMainVolume(0,G4ThreeVector(),fExperimentalHall_log,false,0,this);
  }

  return fExperimentalHall_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::ConstructSDandField() {

  if (!fMainVolume) return;

  // PMT SD
  if (!fPmt_SD.Get()) {
    //Created here so it exists as pmts are being placed
    G4cout << "Construction /LXeDet/pmtSD" << G4endl;
    LXePMTSD* pmt_SD = new LXePMTSD("/LXeDet/pmtSD");
    fPmt_SD.Put(pmt_SD);

    pmt_SD->InitPMTs(4); //let pmtSD know # of pmts
    pmt_SD->SetPmtPositions(fMainVolume->GetPmtPositions());
  }
  G4SDManager::GetSDMpointer()->AddNewDetector(fPmt_SD.Get());
  //sensitive detector is not actually on the photocathode.
  //processHits gets done manually by the stepping action.
  //It is used to detect when photons hit and get absorbed&detected at the
  //boundary to the photocathode (which doesnt get done by attaching it to a
  //logical volume.
  //It does however need to be attached to something or else it doesnt get
  //reset at the begining of events

  SetSensitiveDetector(fMainVolume->GetLogMPPC(), fPmt_SD.Get());

  // Scint SD

  if (!fScint_SD.Get()) {
    G4cout << "Construction /LXeDet/scintSD" << G4endl;
    LXeScintSD* scint_SD = new LXeScintSD("/LXeDet/scintSD");
    fScint_SD.Put(scint_SD);
  }
  G4SDManager::GetSDMpointer()->AddNewDetector(fScint_SD.Get());
  SetSensitiveDetector(fMainVolume->GetLogScint(), fScint_SD.Get());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetDimensions(G4ThreeVector dims) {
  this->fScint_x=dims[0];
  this->fScint_y=dims[1];
  this->fScint_z=dims[2];
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void LXeDetectorConstruction::SetDefaults() {

  //Resets to default values
	fD_mtl = 0.2*mm;
  fScint_x = 20.0*mm; // Scintillator length x
  fScint_y = 5.0*mm; // Scintillator length y
  fScint_z = 20.0*mm; // Scintillator length z

  fMainVolumeOn=true;
  fMainVolume=NULL;

  G4UImanager::GetUIpointer()
    ->ApplyCommand("/LXe/detector/scintYieldFactor 1.");

  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetMainVolumeOn(G4bool b) {
  fMainVolumeOn=b;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
