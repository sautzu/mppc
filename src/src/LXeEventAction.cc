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

G4int TotalCount;
G4int Produce;
G4float EDEP;
G4int Photon_Sum0;
G4int Photon_Sum1;

G4int i,j,k;
G4int pmts;
G4int n_trajectories;
G4int n_hit;

G4double edep;
G4double edepMax;

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
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeEventAction::EndOfEventAction(const G4Event* anEvent){

  LXeUserEventInformation* eventInformation
    = (LXeUserEventInformation*)anEvent->GetUserInformation();
 
  G4TrajectoryContainer* trajectoryContainer = anEvent->GetTrajectoryContainer();
 
  n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

  // extract the trajectories and draw them
  if (G4VVisManager::GetConcreteInstance()){
    for (G4int i=0; i<n_trajectories; i++){
      LXeTrajectory* trj = (LXeTrajectory*)
        ((*(anEvent->GetTrajectoryContainer()))[i]);
      if(trj->GetParticleName()=="opticalphoton"){
        trj->SetForceDrawTrajectory(fForcedrawphotons);
        trj->SetForceNoDrawTrajectory(fForcenophotons);
      }
      trj->DrawTrajectory();
    }
  }
  LXeScintHitsCollection* scintHC = 0;
  LXePMTHitsCollection* pmtHC = 0;
  G4HCofThisEvent* hitsCE = anEvent->GetHCofThisEvent();
 
  //Get the hit collections
  if(hitsCE){
    if(fScintCollID>=0)scintHC = (LXeScintHitsCollection*)(hitsCE->GetHC(fScintCollID));
    if(fPMTCollID>=0)pmtHC = (LXePMTHitsCollection*)(hitsCE->GetHC(fPMTCollID));
  }

  //Hits in scintillator
  if(scintHC){
    n_hit = scintHC->entries();
    G4ThreeVector  eWeightPos(0.);
    edepMax=0;

    for(int i=0;i<n_hit;i++){ //gather info on hits in scintillator
      edep=(*scintHC)[i]->GetEdep();
      eventInformation->IncEDep(edep); //sum up the edep
      eWeightPos += (*scintHC)[i]->GetPos()*edep;//calculate energy weighted pos
      if(edep>edepMax){
        edepMax=edep;//store max energy deposit
        G4ThreeVector posMax=(*scintHC)[i]->GetPos();
        eventInformation->SetPosMax(posMax,edep);
      }
    }
    if(eventInformation->GetEDep()==0.){
    }else{
      //Finish calculation of energy weighted position
      eWeightPos/=eventInformation->GetEDep();
      eventInformation->SetEWeightPos(eWeightPos);
    }
  }
 
	EDEP = eventInformation->GetEDep() / keV ;

	//hits in pmt
  if(pmtHC){
    G4ThreeVector reconPos(0.,0.,0.);
    G4int pmts=pmtHC->entries();
    //Gather info from all PMTs
    for(G4int i=0;i<pmts;i++){
      eventInformation->IncHitCount((*pmtHC)[i]->GetPhotonCount());
			if((*pmtHC)[i]->GetPMTNumber()==0){
				Photon_Sum0 = (*pmtHC)[i]->GetPhotonCount();
				//G4cout << "PMT0   " << Photon_Sum0 << G4endl;
			}else if((*pmtHC)[i]->GetPMTNumber()==1){
				Photon_Sum1 = (*pmtHC)[i]->GetPhotonCount();
				//G4cout << "PMT1   " << Photon_Sum1 << G4endl;
			}

      reconPos+=(*pmtHC)[i]->GetPMTPos()*(*pmtHC)[i]->GetPhotonCount();
      if((*pmtHC)[i]->GetPhotonCount()>=fPMTThreshold){
        eventInformation->IncPMTSAboveThreshold();
      }
      else{//wasnt above the threshold, turn it back off
        (*pmtHC)[i]->SetDrawit(false);
      }
    }
 
    if(eventInformation->GetHitCount()>0){//dont bother unless there were hits
      reconPos/=eventInformation->GetHitCount();
      if(fVerbose>0){
        G4cout << "\tReconstructed position of hits in LXe : "
               << reconPos/mm << G4endl;
      }
      eventInformation->SetReconPos(reconPos);
    }
    pmtHC->DrawAllHits();
  }
	
	TotalCount = eventInformation->GetHitCount();
	Produce = eventInformation->GetPhotonCount_Scint();
  //If we have set the flag to save 'special' events, save here
  if(fSaveThreshold&&eventInformation->GetPhotonCount() <= fSaveThreshold)
    G4RunManager::GetRunManager()->rndmSaveThisEvent();

  if(fRecorder)fRecorder->RecordEndOfEvent(anEvent);

	if(TotalCount!=0){
		std::ofstream ofs1("/Users/Alisa/cern/geant4/Simulations/COPYFILE/compile/result/CountValue.txt",std::ios::out | std::ios::app);
		EDEP = EDEP/2.; //EDEPは2回カウントしてるみたい
		EDEP = 0;
//		EDEP = CLHEP::RandGauss::shoot(EDEP,0.1*EDEP);
//		Photon_Sum0 = CLHEP::RandGauss::shoot(Photon_Sum0,0.01*Photon_Sum0);
//		Photon_Sum1 = CLHEP::RandGauss::shoot(Photon_Sum1,0.01*Photon_Sum1);
		ofs1 << TotalCount << " " << Produce << " " << EDEP << " " << 
		Photon_Sum0 << " " << Photon_Sum1 << " " << (Photon_Sum0+Photon_Sum1) 
		<< std::endl; 

		G4cout << "PMT0   " << Photon_Sum0 << G4endl;
		G4cout << "PMT1   " << Photon_Sum1 << G4endl;
	}
	TotalCount = 0;
	Produce = 0;
	EDEP = 0;
	Photon_Sum0 = 0;
	Photon_Sum1 = 0;

	G4cout << "EndOfEvent" << G4endl;
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
