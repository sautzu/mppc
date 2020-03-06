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
// $Id: LXeDetectorConstruction.cc 82853 2014-07-14 09:07:11Z gcosmo $
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

LXeDetectorConstruction* LXeDetectorConstruction::fInstance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeDetectorConstruction::LXeDetectorConstruction()
: fGSO_mt(NULL), fNaI_mt(NULL), fvacuum_mt(NULL), fAir_mt(NULL),
	fPTFE_mt(NULL), fPMMA_mt(NULL), fSilicone_mt(NULL), fAl_mt(NULL)
{
	fNistMan = G4NistManager::Instance();
	fNistMan->SetVerbose(2);

  fExperimentalHall_box = NULL;
  fExperimentalHall_log = NULL;
  fExperimentalHall_phys = NULL;

	fAl = fAir = fVacuum = fNaI = fGSO = fPTFE = fPMMA = fSilicone = NULL;

  SetDefaults();
	DefineMaterials();

  fDetectorMessenger = new LXeDetectorMessenger(this);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
LXeDetectorConstruction::~LXeDetectorConstruction() {
	delete fAl;						delete fAir;						delete fVacuum;
	delete fNaI;					delete fGSO;						delete fPTFE;
	delete fPMMA;					delete fSilicone;				delete fAl_mt;
	delete fAir_mt;				delete fvacuum_mt;			delete fNaI_mt;
	delete fGSO_mt;				delete fPTFE_mt;				delete fPMMA_mt;
	delete fSilicone_mt;	
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
LXeDetectorConstruction* LXeDetectorConstruction::GetInstance()
{
	if(fInstance==0)
	{
		fInstance = new LXeDetectorConstruction();
	}
	return fInstance;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4Material* LXeDetectorConstruction::GetMaterial(const G4String material)
{
	G4Material* mat = fNistMan->FindOrBuildMaterial(material);
	if(!mat) mat = G4Material::GetMaterial(material);
	if(!mat){
		std::ostringstream o;
		o << "Material " << material << " not found!";
		G4Exception("LXeDetectorConstruction","",
								FatalException,o.str().c_str());
	}
	return mat;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void LXeDetectorConstruction::DefineMaterials(){
  G4double density;
	G4int ncomponents;
	G4double fractionmass;
	std::vector<G4int> natoms;
	std::vector<G4double>	fractionMass;
	std::vector<G4String>	elements;

  //+++Materials
	//Vacuum
	fVacuum 		= fNistMan->FindOrBuildMaterial("G4_Galactic");
	fAir 				= fNistMan->FindOrBuildMaterial("G4_AIR");
	fAl 				= fNistMan->FindOrBuildMaterial("G4_Al");
	//NaI
	elements.push_back("Na");	natoms.push_back(1);
	elements.push_back("I");	natoms.push_back(1);
	density = 3.67*g/cm3;
	fNaI=fNistMan->ConstructNewMaterial("NaI",elements,natoms,density);
	elements.clear();
	natoms.clear();
	//PTFE
	elements.push_back("C");	natoms.push_back(2);
	elements.push_back("F");	natoms.push_back(4);
	density = 2.2*g/cm3;
	fPTFE=fNistMan->ConstructNewMaterial("PTFE",elements,natoms,density);
	elements.clear();
	natoms.clear();
	//GSO
	elements.push_back("Gd");	natoms.push_back(2);
	elements.push_back("Si");	natoms.push_back(1);
	elements.push_back("O");	natoms.push_back(5);
	density = 6.71*g/cm3;
	fGSO=fNistMan->ConstructNewMaterial("GSO",elements,natoms,density);
	elements.clear();
	natoms.clear();
	//PMMA
	elements.push_back("C");	natoms.push_back(5);
	elements.push_back("O");	natoms.push_back(2);
	elements.push_back("H");	natoms.push_back(8);
	density = 1.18*g/cm3;
	fPMMA=fNistMan->ConstructNewMaterial("PMMA",elements,natoms,density);
	elements.clear();
	natoms.clear();
	//Silicone
	elements.push_back("C");	natoms.push_back(2);
	elements.push_back("H");	natoms.push_back(6);
	density = 1.060*g/cm3;
	fSilicone=fNistMan->ConstructNewMaterial("Silicone",elements,natoms,density);
	elements.clear();
	natoms.clear();

  
	//+++Material properties tables
	G4double photon_Energy[50];

	G4double vacuum_RIND[50];
	G4double air_RIND[50];
	G4double nai_RIND[50];
	G4double gso_RIND[50];
	G4double ptfe_RIND[50];
	G4double pmma_RIND[50];
	G4double silicone_RIND[50];
	G4double al_RIND[50];
	
	G4double nai_ABSL[50];
	G4double gso_ABSL[50];
	G4double ptfe_ABSL[50];

	for(int i=0;i<50;i++){
		photon_Energy[i] = 1240.0/(620.0-5.0*i);
		vacuum_RIND[i] 	 = 1.000293;
		air_RIND[i] 		 = 1.0003;
		nai_RIND[i]			 = 1.85;
		gso_RIND[i]			 = 1.85;
		ptfe_RIND[i]		 = 1.35;
		pmma_RIND[i]		 = 1.50;
		silicone_RIND[i] = 1.41; 
		al_RIND[i]			 = 1.48;
		nai_ABSL[i]			 = 2.59*m;/////
		gso_ABSL[i]			 = 1.38*m;/////
		ptfe_ABSL[i] 		 = 20.0*m;
	}
	const G4int nEntries = sizeof(photon_Energy)/sizeof(G4double);
	G4double nai_SCINT[] = 
	{
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.02, 0.04, 0.06, 0.08, 0.1,
    0.12, 0.14, 0.16, 0.18, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45,
    0.5,  0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95,
    0.98, 1.00, 0.98, 0.95, 0.9,  0.85, 0.8,  0.7,  0.6,  0.5
	};
	G4double gso_SCINT[] =
	{
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.02, 0.04, 0.06, 0.08, 0.1,  0.12, 0.14, 0.16, 0.18, 0.2,
    0.25, 0.3,  0.35, 0.4,  0.45, 0.5,  0.55, 0.6,  0.65, 0.7,
    0.75, 0.8,  0.85, 0.9,  0.95, 0.98, 1.00, 0.98, 0.95, 0.9,
    0.85, 0.8,  0.7,  0.6,  0.5,  0.45, 0.4,  0.35, 0.3,  0.25
	};

	//GSO
 	fGSO_mt = new G4MaterialPropertiesTable();
  fGSO_mt->AddProperty("ABSLENGTH",photon_Energy,gso_ABSL,nEntries);
  fGSO_mt->AddProperty("RINDEX",photon_Energy,gso_RIND,nEntries);
	//fGSO_mt->AddProperty("FASTCOMPONENT", photon_Energy, gso_SCINT, nEntries);
	fGSO_mt->AddProperty("SCINTILLATION", photon_Energy, gso_SCINT, nEntries);
	fGSO_mt->AddConstProperty("SCINTILLATIONYIELD", 10000./MeV );//1MeVあたりの出力光子数
  fGSO_mt->AddConstProperty("RESOLUTIONSCALE", 3.1 );//(R/2.35)*(E*ScintillationYield)^(1/2)//分解能
  fGSO_mt->AddConstProperty("FASTTIMECONSTANT", 30.*ns );//時定数(早)
  fGSO_mt->AddConstProperty("SLOWTIMECONSTANT", 600.*ns );//時定数(遅)
  fGSO_mt->AddConstProperty("YIELDRATIO",0.9);//時定数の比
  fGSO->SetMaterialPropertiesTable(fGSO_mt);
  fGSO->GetIonisation()->SetBirksConstant(0.126*mm/MeV);//Birk定数
	//NaI
  fNaI_mt = new G4MaterialPropertiesTable();
  fNaI_mt->AddProperty("ABSLENGTH",photon_Energy,nai_ABSL,nEntries);
  fNaI_mt->AddProperty("RINDEX",photon_Energy,nai_RIND,nEntries);
	//fNaI_mt->AddProperty("FASTCOMPONENT", photon_Energy, nai_SCINT, nEntries);
	fNaI_mt->AddProperty("SCINTILLATION", photon_Energy, nai_SCINT, nEntries);
	fNaI_mt->AddConstProperty("SCINTILLATIONYIELD", 38000./MeV );
  fNaI_mt->AddConstProperty("RESOLUTIONSCALE", 5.0 );
  fNaI_mt->AddConstProperty("FASTTIMECONSTANT", 230.*ns );
  fNaI_mt->AddConstProperty("YIELDRATIO",1.0);
  fNaI->SetMaterialPropertiesTable(fNaI_mt);
  fNaI->GetIonisation()->SetBirksConstant(0.126*mm/MeV);
	//PTFE
  fPTFE_mt = new G4MaterialPropertiesTable();
  fPTFE_mt->AddProperty("ABSLENGTH",photon_Energy,ptfe_ABSL,nEntries);
  fPTFE_mt->AddProperty("RINDEX",photon_Energy,ptfe_RIND,nEntries);
  fPTFE->SetMaterialPropertiesTable(fPTFE_mt);
	//PMMA
  fPMMA_mt = new G4MaterialPropertiesTable();
  fPMMA_mt->AddProperty("ABSLENGTH",photon_Energy,ptfe_ABSL,nEntries);
  fPMMA_mt->AddProperty("RINDEX",photon_Energy,pmma_RIND,nEntries);
  fPMMA->SetMaterialPropertiesTable(fPMMA_mt);
	//Silicone
  fSilicone_mt = new G4MaterialPropertiesTable();
  fSilicone_mt->AddProperty("ABSLENGTH",photon_Energy,ptfe_ABSL,nEntries);
  fSilicone_mt->AddProperty("RINDEX",photon_Energy,silicone_RIND,nEntries);
  fSilicone->SetMaterialPropertiesTable(fSilicone_mt);
	//Silicone
  fAl_mt = new G4MaterialPropertiesTable();
  fAl_mt->AddProperty("RINDEX",photon_Energy,al_RIND,nEntries);
  fAl->SetMaterialPropertiesTable(fAl_mt);
	//vacuum
  fvacuum_mt = new G4MaterialPropertiesTable();
  fvacuum_mt->AddProperty("RINDEX", photon_Energy, vacuum_RIND, nEntries);
  fVacuum->SetMaterialPropertiesTable(fvacuum_mt);
	//air
  fAir_mt = new G4MaterialPropertiesTable();
  fAir_mt->AddProperty("RINDEX", photon_Energy, air_RIND, nEntries);
  fAir->SetMaterialPropertiesTable(fAir_mt);//Give air the same rindex
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
  //DefineMaterials();
  return ConstructDetector();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* LXeDetectorConstruction::ConstructDetector()
{
  //The experimental hall walls are all 1m away from housing walls
  expHall_x = fScint_x+fD_mtl+1.*m;
  expHall_y = fScint_y+fD_mtl+1.*m;
  expHall_z = fScint_z+fD_mtl+1.*m;

  //Create experimental hall
  fExperimentalHall_box
    = new G4Box("expHall_box",expHall_x,expHall_y,expHall_z);
  fExperimentalHall_log = new G4LogicalVolume(fExperimentalHall_box,
                                             fAir,"expHall_log",0,0,0);
  fExperimentalHall_phys = new G4PVPlacement(0,G4ThreeVector(),
                              fExperimentalHall_log,"expHall",0,false,0);

  fExperimentalHall_log->SetVisAttributes(G4VisAttributes::Invisible);
  //Place the main volume
  if(fMainVolumeOn){
    fMainVolume
      = new LXeMainVolume(0,G4ThreeVector(),fExperimentalHall_log,false,0,this);
  }
	//delete fExperimentalHall_box;
	//delete fExperimentalHall_log;
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
    pmt_SD->InitPMTs((fNx*fNy+fNx*fNz+fNy*fNz)*2); //let pmtSD know # of pmts
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
  SetSensitiveDetector(fMainVolume->GetLogPhotoCath(), fPmt_SD.Get());
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

void LXeDetectorConstruction::SetHousingThickness(G4double d_mtl) {
  this->fD_mtl=d_mtl;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetNX(G4int nx) {
  this->fNx=nx;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetNY(G4int ny) {
  this->fNy=ny;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetNZ(G4int nz) {
  this->fNz=nz;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetPMTRadius(G4double outerRadius_pmt) {
  this->fOuterRadius_pmt=outerRadius_pmt;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void LXeDetectorConstruction::SetDefaults() {

  //Resets to default values
  fD_mtl=0.2*mm;
	expHall_x = 1*m;
	expHall_y = 1*m;
	expHall_z = 1*m;
  fScint_x = 10.*mm;
  fScint_y = 10.*mm;
  fScint_z = 40.*mm;
  fNx = 1;
  fNy = 1;
  fNz = 1;
  fOuterRadius_pmt = 2.3*cm;
  fRefl=1.0;
  fMainVolumeOn=true;
  fMainVolume=NULL;
  G4UImanager::GetUIpointer()
    ->ApplyCommand("/LXe/detector/scintYieldFactor 1.");
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetHousingReflectivity(G4double r) {
  fRefl=r;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

