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
//G4long seed ;
//G4long runCounter;
G4double theta, psi ,psi2;
G4double x,y,z,u,f;
G4double X,Y,Z;
G4double k1,k2,k3;
G4double x_test,z_test;
G4double r_test,theta_test;
G4double scint_y = 5.0/2.*mm;
G4double scint_z = 20.0/2.*mm;

//G4double y_position = -86.7*mm;
G4double y_position = -89.6*mm; //下
//G4double y_position = -scint_y - 86.6-0.25-18.*mm;//下
//G4double y_position = -scint_y - 86.6-0.25*mm;//下
//G4double y_position = scint_y + 36.01 + 70.*mm; //上
G4double x_position = 0.0*mm;
G4double z_position = 0.0*mm;//-scint_z/2.; //3.1*mm;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


LXePrimaryGeneratorAction::LXePrimaryGeneratorAction(){
for(int i=0;i<1;i++){
  G4int n_particle = 1;
  fParticleGun[i] = new G4ParticleGun(n_particle);
 
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
 
  G4String particleName;
  fParticleGun[i]->SetParticleDefinition(particleTable->
                                     FindParticle("gamma"));
                                     //FindParticle("opticalphoton"));
}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
LXePrimaryGeneratorAction::~LXePrimaryGeneratorAction(){
for(int i=0;i<1;i++){
    delete fParticleGun[i];
}
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXePrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent){
  
for(int i=0;i<1;i++){
	
//	std::fstream iF(("/home/kinoshita/simulation/result/runCounter"), std::ios::in);
//	if (iF) {
//	iF >> runCounter;
//	}
//	else{
//	runCounter = 0;
//	}
//	iF.close();
//	
//	seed = runCounter; 
//	runCounter = runCounter + 1;  
//	
//	
//	std::fstream oF(("/home/kinoshita/simulation/result/runCounter"), std::ios::out);
//	oF << runCounter << std::endl;
//	oF.close();
//

  seed = seed + 1;
//        std::cout<<" "<< seed <<std::endl;
	CLHEP::HepRandom::setTheSeed(seed);
//	CLHEP::HepRandom::setTheSeed(1573110511);

//	psi = CLHEP::RandFlat::shoot(0.,2.0*M_PI);
	psi2 = CLHEP::RandFlat::shoot(0.25*M_PI,0.75*M_PI);
//	u = CLHEP::RandFlat::shoot(-1.0,1.0);

//	r_test = CLHEP::RandFlat::shoot(8.0,14.2);
//	theta_test = CLHEP::RandFlat::shoot(0.,2.0*M_PI);

//	x_test = r_test*cos(theta_test);
//	z_test = r_test*sin(theta_test);

	k1 =u;
	k2 =(sqrt(1-(u*u)))*sin(psi2);
	k3 =(sqrt(1-(u*u)))*cos(psi2);

	X = CLHEP::RandFlat::shoot(-3.0,3.0);
//	//X = CLHEP::RandFlat::shoot(-0.5,0.5);
	Y = CLHEP::RandFlat::shoot(-0.5,0.5);
	Z = CLHEP::RandFlat::shoot(-sqrt(3.*3.-X*X),sqrt(3.*3.-X*X));
	
	x = x_position + X;
	y = y_position + Y;
	z = z_position + Z;


	G4ThreeVector position(x*mm, y*mm, z*mm);
//	G4ThreeVector position(0.0*mm, y_position, z_position*mm);
//	G4ThreeVector position(0, -scint_y-86.6*mm, -7.0/3.0);
//	G4ThreeVector position(0, -scint_y-86.6*mm, 0.0);
//	G4ThreeVector position(0, -scint_y-86.6*mm, scint_z-3.);
	G4ThreeVector momentum(k1*mm, k2*mm, k3*mm);
//	G4ThreeVector momentum(x_test*mm,-106.01 *mm, z_test*mm);
//	G4ThreeVector momentum(0, 1, 0);

	fParticleGun[i]->SetParticleEnergy(662*keV);
	fParticleGun[i]->SetParticlePosition(position);
	fParticleGun[i]->SetParticleMomentumDirection(momentum);
	
	fParticleGun[i]->GeneratePrimaryVertex(anEvent);
}
}
