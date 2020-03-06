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
// $Id: LXePrimaryGeneratorAction.cc 68752 2013-04-05 10:23:47Z gcosmo $
//
/// \file optical/LXe/src/LXePrimaryGeneratorAction.cc
/// \brief Implementation of the LXePrimaryGeneratorAction class
//
//
#include "LXePrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

#include <iostream>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include "G4RandomDirection.hh"
#include "Randomize.hh"

G4long seed = time(NULL);
G4double theta;
G4double x, y, z, u, f;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXePrimaryGeneratorAction::LXePrimaryGeneratorAction(){
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);
 
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
 
  G4String particleName;
  fParticleGun->SetParticleDefinition(particleTable->
                                     FindParticle(particleName="gamma"));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXePrimaryGeneratorAction::~LXePrimaryGeneratorAction(){
    delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXePrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent){

	seed = seed + 1;
	CLHEP::HepRandom::setTheSeed(seed);

/*	theta = CLHEP::RandFlat::shoot(0.0,2.0*M_PI);
	u = CLHEP::RandFlat::shoot(-1.0,1.0);
	f = CLHEP::RandFlat::shoot(-1.0,1.0);
*/
/*	x = CLHEP::RandFlat::shoot(-3.0,3.0);
	y = CLHEP::RandFlat::shoot(-sqrt(3.0*3.0 - x*x),sqrt(3.0*3.0 - x*x));
	z = CLHEP::RandFlat::shoot(-0.5,0.5);
*/
 
  //Default energy,position,momentum
  fParticleGun->SetParticleEnergy(662.*keV);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,-1.,0.));

	x = CLHEP::RandFlat::shoot(-3.0,3.0);
	y = CLHEP::RandFlat::shoot(-sqrt(3.0*3.0-x*x),sqrt(3.0*3.0-x*x));
	z = CLHEP::RandFlat::shoot(-0.5,0.5);

  fParticleGun->SetParticlePosition(G4ThreeVector(x*cm , 5.*cm, z*cm));
  //fParticleGun->SetParticlePosition(G4ThreeVector(0*cm , 5.*cm, 0*cm));

  fParticleGun->GeneratePrimaryVertex(anEvent);
}
