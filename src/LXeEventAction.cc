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
// $Id: LXeEventAction.cc 68752 2013-04-05 10:23:47Z gcosmo $
//
/// \file optical/LXe/src/LXeEventAction.cc
/// \brief Implementation of the LXeEventAction class
//
//
#include "LXeEventAction.hh"
#include "LXeScintHit.hh"
#include "LXePMTHit.hh"
#include "LXeUserEventInformation.hh"
#include "LXeTrajectory.hh"
#include "LXeRecorderBase.hh"

#include "G4EventManager.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4UImanager.hh"
#include "G4SystemOfUnits.hh"

#include <iostream>
#include <fstream>

int eventnumber = 0; 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeEventAction::LXeEventAction(LXeRecorderBase* r)
  : fRecorder(r),fSaveThreshold(0),fScintCollID(-1),fPMTCollID(-1),fVerbose(0),
   fPMTThreshold(1),fForcedrawphotons(false),fForcenophotons(false)
{
  fEventMessenger = new LXeEventMessenger(this);
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeEventAction::~LXeEventAction(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeEventAction::BeginOfEventAction(const G4Event* anEvent){
 
  //New event, add the user information object
  G4EventManager::
    GetEventManager()->SetUserInformation(new LXeUserEventInformation);

  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  if(fScintCollID<0)
    fScintCollID=SDman->GetCollectionID("scintCollection");
  if(fPMTCollID<0)
    fPMTCollID=SDman->GetCollectionID("pmtHitCollection");

  if(fRecorder)fRecorder->RecordBeginOfEvent(anEvent);
	
	LXeUserEventInformation* eventInformation
	=(LXeUserEventInformation*)anEvent->GetUserInformation();

	eventnumber++;
//        std::cout<< eventnumber << std::endl;
  eventInformation->SetEventNumber(eventnumber);

	
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeEventAction::EndOfEventAction(const G4Event* anEvent){

	LXeUserEventInformation* eventInformation
	=(LXeUserEventInformation*)anEvent->GetUserInformation();
 	
	G4TrajectoryContainer* trajectoryContainer=anEvent->GetTrajectoryContainer();
 
	G4int n_trajectories = 0;
	if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

	LXeScintHitsCollection* scintHC = 0;
	LXePMTHitsCollection* pmtHC = 0;
	G4HCofThisEvent* hitsCE = anEvent->GetHCofThisEvent();
 
	//Get the hit collections
	if(hitsCE){
		if(fScintCollID>=0)scintHC 
			= (LXeScintHitsCollection*)(hitsCE->GetHC(fScintCollID));
		if(fPMTCollID>=0)pmtHC 
			= (LXePMTHitsCollection*)(hitsCE->GetHC(fPMTCollID));
	}

	//Hits in scintillator
	if(scintHC){
		int n_hit = scintHC->entries();
		G4ThreeVector  eWeightPos(0.);
		G4double edep;
		
		for(int i=0;i<n_hit;i++){ //gather info on hits in scintillator
			edep=(*scintHC)[i]->GetEdep();
			eventInformation->IncEDep(edep); //sum up the edep
		}
	}
	G4double EDEP = eventInformation->GetEDep() / keV ;
	const G4int mppc_num = 3;
	G4int Photon_Sum[mppc_num];
	for(int i=0;i<mppc_num;i++)
		Photon_Sum[i] = 0;
	//Hits in pmt 
	if(pmtHC){
		G4int pmts=pmtHC->entries();
//		std::cout<<pmtHC<<std::endl;
		//Gather info from all PMTs
		for(G4int i=0;i<pmts;i++){
			eventInformation->IncHitCount((*pmtHC)[i]->GetPhotonCount());
			Photon_Sum[(*pmtHC)[i]->GetPMTNumber()-1] = (*pmtHC)[i]->GetPhotonCount();
//			std::cout<<(*pmtHC)[i]->GetPMTNumber()-1<<std::endl;
			
			
			if((*pmtHC)[i]->GetPhotonCount()>=fPMTThreshold){
				eventInformation->IncPMTSAboveThreshold();
			}else{//wasnt above the threshold, turn it back off
				(*pmtHC)[i]->SetDrawit(false);
			}
		}
		pmtHC->DrawAllHits();
	}
	
	G4int TotalCount = eventInformation->GetHitCount();
	G4int Produce = eventInformation->GetPhotonCount_Scint();
	

	//If we have set the flag to save 'special' events, save here
	if(fSaveThreshold&&eventInformation->GetPhotonCount() <= fSaveThreshold)
		G4RunManager::GetRunManager()->rndmSaveThisEvent();
	
	if(fRecorder)fRecorder->RecordEndOfEvent(anEvent);

	if(TotalCount!=0){
	//if(EDEP!=0){
		std::ofstream ofs1("../result/test.txt",std::ios::out | std::ios::app);
		//EDEP = EDEP/2.; //EDEPは2回カウントしてるみたい
		//EDEP = CLHEP::RandGauss::shoot(EDEP,0.1*EDEP);
		ofs1 
//<< TotalCount << " " 
<< eventnumber << " "
<< Produce << " " 
//<< EDEP << " " 
<< Photon_Sum[0] << " " 
<< Photon_Sum[1] << " " 
<< Photon_Sum[2] << " " 
//<< Photon_Sum[3] << " " 
//<< Photon_Sum[4] << " " 
//<< Photon_Sum[5] << " " 
//<< Photon_Sum[6] << " " 
//<< Photon_Sum[7] << " " 
//<< Photon_Sum[8] << " " 
//<< Photon_Sum[9] << " " 
//<< Photon_Sum[10] << " " 
//<< Photon_Sum[11] << " " 
//<< Photon_Sum[12] << " " 
//<< Photon_Sum[13] << " " 
//<< Photon_Sum[14] << " " 
//<< Photon_Sum[15]
	<< std::endl;
	}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeEventAction::SetSaveThreshold(G4int save){
/*Sets the save threshold for the random number seed. If the number of photons
generated in an event is lower than this, then save the seed for this event
in a file called run###evt###.rndm
*/
  fSaveThreshold=save;
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  G4RunManager::GetRunManager()->SetRandomNumberStoreDir("random/");
  //  G4UImanager::GetUIpointer()->ApplyCommand("/random/setSavingFlag 1");
}
