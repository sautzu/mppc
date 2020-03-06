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
// $Id: LXeSteppingAction.cc 73915 2013-09-17 07:32:26Z gcosmo $
//
/// \file optical/LXe/src/LXeSteppingAction.cc
/// \brief Implementation of the LXeSteppingAction class
//
//
#include "LXeSteppingAction.hh"
#include "LXeEventAction.hh"
#include "LXeTrackingAction.hh"
#include "LXeTrajectory.hh"
#include "LXePMTSD.hh"
#include "LXeUserTrackInformation.hh"
#include "LXeUserEventInformation.hh"
#include "LXeSteppingMessenger.hh"
#include "LXeRecorderBase.hh"

#include "G4SteppingManager.hh"
#include "G4SDManager.hh"
#include "G4EventManager.hh"
#include "G4ProcessManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeSteppingAction::LXeSteppingAction(LXeRecorderBase* r)
  : fRecorder(r),fOneStepPrimaries(false)
{
  fSteppingMessenger = new LXeSteppingMessenger(this);

  fExpectedNextStatus = Undefined;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeSteppingAction::~LXeSteppingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeSteppingAction::UserSteppingAction(const G4Step * theStep){

  G4Track* theTrack = theStep->GetTrack();
  if ( theTrack->GetCurrentStepNumber() == 1 ) fExpectedNextStatus = Undefined;
	//現在のStep番号==1:初めのStep
 
  LXeUserTrackInformation* trackInformation
    =(LXeUserTrackInformation*)theTrack->GetUserInformation();
  LXeUserEventInformation* eventInformation
    =(LXeUserEventInformation*)G4EventManager::GetEventManager()
    ->GetConstCurrentEvent()->GetUserInformation();

  G4StepPoint* 				thePrePoint 	= theStep->GetPreStepPoint();
  G4VPhysicalVolume* 	thePrePV 			= thePrePoint->GetPhysicalVolume();
  G4StepPoint* 				thePostPoint 	= theStep->GetPostStepPoint();
  G4VPhysicalVolume* 	thePostPV 		= thePostPoint->GetPhysicalVolume();

  G4OpBoundaryProcessStatus boundaryStatus=Undefined;
  static G4ThreadLocal G4OpBoundaryProcess* boundary=NULL;

  //find the boundary process only once
  if(!boundary){
    G4ProcessManager* pm
      = theStep->GetTrack()->GetDefinition()->GetProcessManager();
    G4int nprocesses = pm->GetProcessListLength();
    G4ProcessVector* pv = pm->GetProcessList();
    G4int i;
    for( i=0;i<nprocesses;i++){
      if((*pv)[i]->GetProcessName()=="OpBoundary"){
        boundary = (G4OpBoundaryProcess*)(*pv)[i];
        break;
      }
    }
  }

  if(theTrack->GetParentID()==0){
    //This is a primary track
    G4TrackVector* fSecondary=fpSteppingManager->GetfSecondary();
    G4int tN2ndariesTot = fpSteppingManager->GetfN2ndariesAtRestDoIt()
      + fpSteppingManager->GetfN2ndariesAlongStepDoIt()
      + fpSteppingManager->GetfN2ndariesPostStepDoIt();

  	//G4String pn = theTrack->GetDefinition()->GetParticleName();
		//G4String processName = thePostPoint->GetProcessDefinedStep()->GetProcessName();
		//if(processName=="phot"||processName=="compt"||processName=="conv"){
		//	if(thePrePV->GetName()=="scintillator")
    //  	theTrack->SetTrackStatus(fKillTrackAndSecondaries);
		//}
 
  	//if(processName=="compt" && thePrePV->GetName()=="scintillator" && pn=="gamma" ){
  	//	theTrack->SetTrackStatus(fStopAndKill);
		//}
	
	}

  if(!thePostPV){//out of world
    fExpectedNextStatus=Undefined;
    return;
  }

  G4ParticleDefinition* particleType = theTrack->GetDefinition();
  if(particleType==G4OpticalPhoton::OpticalPhotonDefinition()){
  	if(thePrePV->GetName()=="Silica" && thePostPV->GetName()=="Air"){	
			theTrack->SetTrackStatus(fStopAndKill);
  	}
	}


	if(particleType==G4OpticalPhoton::OpticalPhotonDefinition()){
    //Optical photon only

    //if(thePostPV->GetName()=="expHall"){
      //Kill photons entering expHall from something other than Slab
      //theTrack->SetTrackStatus(fStopAndKill);
      //theTrack->SetTrackStatus(fStopAndKill);
	//	}

    boundaryStatus=boundary->GetStatus();

    //Check to see if the partcile was actually at a boundary
    //Otherwise the boundary status may not be valid
    //Prior to Geant4.6.0-p1 this would not have been enough to check
    if(thePostPoint->GetStepStatus()==fGeomBoundary){
      if(fExpectedNextStatus==StepTooSmall){
        if(boundaryStatus!=StepTooSmall){
          G4ExceptionDescription ed;
          ed << "LXeSteppingAction::UserSteppingAction(): "
                << "No reallocation step after reflection!"
                << G4endl;
          G4Exception("LXeSteppingAction::UserSteppingAction()", "LXeExpl01",
          FatalException,ed,
          "Something is wrong with the surface normal or geometry");
        }
      }
      fExpectedNextStatus=Undefined;
      switch(boundaryStatus){
      case Detection: //Note, this assumes that the volume causing detection
                      //is the photocathode because it is the only one with
                      //non-zero efficiency
        {
        //Triger sensitive detector manually since photon is
        //absorbed but status was Detection
        G4SDManager* SDman = G4SDManager::GetSDMpointer();
        G4String sdName="/LXeDet/pmtSD";
        LXePMTSD* pmtSD = (LXePMTSD*)SDman->FindSensitiveDetector(sdName);
        if(pmtSD){
					pmtSD->ProcessHits_constStep(theStep,NULL);
				}
				trackInformation->AddTrackStatusFlag(hitPMT);
        break;
        }
      case FresnelReflection:
      case TotalInternalReflection:
      case LambertianReflection:
      case LobeReflection:
      case SpikeReflection:
      case BackScattering:
        trackInformation->IncReflections();
        fExpectedNextStatus=StepTooSmall;
        break;
      default:
        break;
      }
    }
  }
  if(fRecorder)fRecorder->RecordStep(theStep);
}
