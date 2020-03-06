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
// $Id: LXeScintSD.cc 68752 2013-04-05 10:23:47Z gcosmo $
//
/// \file optical/LXe/src/LXeScintSD.cc
/// \brief Implementation of the LXeScintSD class
//
//
#include "LXeScintSD.hh"
#include "LXeScintHit.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"


#include "LXeUserEventInformation.hh"
#include "G4EventManager.hh"
#include "G4ProcessManager.hh"
#include "G4Event.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeScintSD::LXeScintSD(G4String name)
  : G4VSensitiveDetector(name)
{
  fScintCollection = NULL;
  collectionName.insert("scintCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeScintSD::~LXeScintSD() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeScintSD::Initialize(G4HCofThisEvent* hitsCE){
  fScintCollection = new LXeScintHitsCollection
                      (SensitiveDetectorName,collectionName[0]);
  //A way to keep all the hits of this event in one place if needed
  static G4int hitsCID = -1;
  if(hitsCID<0){
    hitsCID = GetCollectionID(0);
  }
  hitsCE->AddHitsCollection( hitsCID, fScintCollection );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4bool LXeScintSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ){
  if(aStep->GetTrack()->GetDefinition()->GetParticleName()=="opticalphoton") return false;

	LXeUserEventInformation* eventInformation
  =(LXeUserEventInformation*)G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetUserInformation();
	G4int eventnumber = eventInformation->GetEventNumber(); 

 
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep==0.) return false; //No edep so dont count as hit

  G4StepPoint* thePrePoint = aStep->GetPreStepPoint();
  G4TouchableHistory* theTouchable =
    (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  G4VPhysicalVolume* thePrePV = theTouchable->GetVolume();
//	G4int copynum = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();
 
  G4StepPoint* thePostPoint = aStep->GetPostStepPoint();

  //Get the average position of the hit
  //G4ThreeVector pos = thePrePoint->GetPosition() + thePostPoint->GetPosition();
  //pos/=2.;

  LXeScintHit* scintHit = new LXeScintHit(thePrePV);

  scintHit->SetEdep(edep);
  //scintHit->SetPos(pos);

  fScintCollection->insert(scintHit);



//  std::ofstream ofs("/home/kinoshita/simulation/kekka/Eventtest.txt",std::ios::out | std::ios::app);
////  if(edep!=0){
//  G4String p  = aStep->GetPreStepPoint()->GetMaterial()->GetName();
//  G4String Pv = aStep->GetPostStepPoint()->GetMaterial()->GetName();
//  G4String pn = aStep->GetTrack()->GetDefinition()->GetParticleName();
//  G4String pr = //aStep->GetTrack()->GetCreatorProcess()->GetProcessName();
//		      aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
// 
////	G4ThreeVector pos = thePrePoint->GetPosition() + thePostPoint->GetPosition();
////	pos/=2.;
//
//  G4ThreeVector pos = thePostPoint->GetPosition();
//
////  if(pn!="opticalphoton"){
////  if(pr=="phot"){
//    ofs<<eventnumber<<" "<<pn<<" "<<pr<<" "<<pos.getX()<<" "<<pos.getY()<<" "<<pos.getZ()<<" "<<edep<<std::endl;
////	        ofs<<eventnumber<<" "<<pn<<" "<<pr<<" "<<p<<" "<<Pv<<" "<<edep<<std::endl;
////    }
////  }
//// 	std::cout<<p<<std::endl;
////  }
////		" "<<Pv<<" "<<pn<<" "<<pr<<std::endl;
  return true;
}
 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeScintSD::EndOfEvent(G4HCofThisEvent* ) {
		//if(edep!=0){
			//std::cout<<"finished"<<std::endl;
		//}
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeScintSD::clear() {} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeScintSD::DrawAll() {} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeScintSD::PrintAll() {} 
