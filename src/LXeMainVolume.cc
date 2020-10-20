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


G4double m_posi = 3.0*mm;
G4double x_posi = 0.0*mm;
G4double z_posi = 0.0*mm;
G4bool colli = false;
const G4int scint_num = 3;
const G4int scintsq_num = 3;
const G4int mppc_num = 3;
      G4double si = 5.;
	
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
                                     G4Material::GetMaterial("Air"),
                                     "temp",0,0,0),
                 "housing",pMotherLogical,pMany,pCopyNo),fConstructor(c)
{
  CopyValues();

	// Logic Size
	// air_housing
	G4double air_x = fScint_x	+ 2*m;
	G4double air_y = fScint_y	+ 2*m;
	G4double air_z = fScint_z	+ 2*m;
	G4double air_in_x = fScint_x	+ 1*m;
	G4double air_in_y = fScint_y	+ 1*m;
	G4double air_in_z = fScint_z	+ 1*m;
	// MPPC
	G4double mppc_height = 3.0*mm;
	G4double mppc_width = 0.05*mm;
	G4double w_posi;
		w_posi = 0.0*mm;
	// window
	G4double win_longside = 6.0*mm;
	G4double win_shortside = 6.0*mm;
	G4double win_width = 0.1*mm;
	//Envelop
	G4double enve_longside = 6.0*mm;
	G4double enve_shortside = 6.0*mm;
	G4double enve_width = 2.0*mm;
	// Collimator
	G4double colli_s_radius = 3.0*mm;
	G4double colli_l_radius = 1.0*mm;
	G4double colli_radius = 30.0*mm;
	G4double colli_width = 80.0*mm;
	//if(colli)
		//colli_width = 100.0*mm;
	//if(!colli)
		//colli_width = 1.0*mm;
	// Source
	G4double radi_radius = 3.0*mm;
	G4double radi_width = 1.0*mm;
	G4double radi_case_radius = 12.5*mm;
	G4double radi_case_width = 6.0*mm;
	//Wood
	G4double wood_height = 600.0*mm;
	G4double wood_width = 15.*mm;
	//Fe
	G4double fe_height = 150.0*mm;
	G4double fe_width = 2.*mm;
	//Supporter
	G4double sup_width = 2.*mm;
	G4double sup_height = 240.*mm;
	G4double sup_yoko = 150.*mm;
	

	// The number of scintirator
	G4double num = scintsq_num;
	G4double d = 0.2;
	G4double X = 0.*mm;
	G4double Z = 0.*mm;

	G4double slide = 70.*mm;
//	G4double slide = -colli_width+70.*mm;
	// Logic Position


 //Scint(複数個)
    G4double scint_posi_x[scint_num];
    G4double scint_posi_y[scint_num];
    G4double scint_posi_z[scint_num];



//		scint_posi_x[0]=0.;
//		scint_posi_y[0]=-d/2.;
//		scint_posi_z[0]=-fScint_z-d;
//
	scint_posi_x[1]=0.;
   	scint_posi_y[0]=-d/2.;
   	scint_posi_z[1]=0.;

	scint_posi_x[2]=0.;
	scint_posi_y[2]=-d/2.;
	scint_posi_z[2]=fScint_z+d;


    scint_posi_x[0]=0.;
    scint_posi_y[0]=-d/2.;
	scint_posi_z[0]=0.;


 //Scint(大きいの１つ)
     G4double scintbig_posi_x;
     G4double scintbig_posi_y;
     G4double scintbig_posi_z;

		scintbig_posi_x = 0.;
		scintbig_posi_y = -d/2.;
		scintbig_posi_z = 0.;


	//Envelope
         G4double enve_posi_x[mppc_num];
         G4double enve_posi_y[mppc_num];
         G4double enve_posi_z[mppc_num];
 

	if(mppc_num == 1){
		for(int i=0;i<mppc_num;i++){
		  enve_posi_x[i] = scint_posi_x[i];
		  enve_posi_y[i] = fScint_y/2. + d/2. + enve_width/2. + win_width;
		  enve_posi_z[i] = scint_posi_z[i];
	      	}
	}
	else{
		for(int i=0;i<mppc_num;i++){
		  enve_posi_x[i] = scint_posi_x[i];
		  enve_posi_y[i] = fScint_y/2. + d/2. + enve_width/2. + win_width;
		  enve_posi_z[i] = (-1+i)*fScint_z/2. - (-1+i)*mppc_height;
		}
	}
	

//	for(int i=0;i<3;i++){
//		enve_posi_x[i] = scint_posi_x[i];
//		enve_posi_y[i] = fScint_y/2. + d/2. + enve_width/2. + win_width;
//		enve_posi_z[i] = scint_posi_z[i]; // - 2*(-1+i)*mppc_height;
//	}
//	for(int i=3;i<5;i++){
//		enve_posi_x[i] = scint_posi_x[0] - 2*(7-2*i)*mppc_height;
//		enve_posi_y[i] = fScint_y/2. + d/2. + enve_width/2. + win_width/2.;
//		enve_posi_z[i] = scint_posi_z[0];
//	}


         //4MPPC

//         for(int i=0;i<2;i++){
//           enve_posi_x[i] = scint_posi_x[0] + si;
//           enve_posi_y[i] = fScint_y/2. + d/2. + enve_width/2. + win_width;
//           enve_posi_z[i] = scint_posi_z[0] - 2*(-1./2.+i)*si;
//         }
//         for(int i=2;i<4;i++){
//           enve_posi_x[i] = scint_posi_x[0] - si;
//           enve_posi_y[i] = fScint_y/2. + d/2. + enve_width/2. + win_width;
//           enve_posi_z[i] = scint_posi_z[0] - 2*(5./2.-i)*si;
//         }




	// MPPC
	G4double mppc_posi_x[mppc_num]; 
	G4double mppc_posi_y[mppc_num]; 
	G4double mppc_posi_z[mppc_num]; 


	for(int i=0;i<mppc_num;i++){
		mppc_posi_x[i] = 0.;
		mppc_posi_y[i] = 0.;
		mppc_posi_z[i] = -win_width/2.+mppc_width/2.;
	}


	//window
	G4double win_posi_x[mppc_num];
	G4double win_posi_y[mppc_num];
	G4double win_posi_z[mppc_num];


	for(int i=0;i<mppc_num;i++){ 
          win_posi_x[i] = scint_posi_x[i];
          win_posi_y[i] = fScint_y/2. + d/2. + win_width/2.;
          win_posi_z[i] = enve_posi_z[i];
	}
	

         //4MPPC
//       for(int i=0;i<2;i++){
//         win_posi_x[i] = scint_posi_x[0] + si;
//         win_posi_y[i] = fScint_y/2. + d/2. + win_width/2.;
//         win_posi_z[i] = scint_posi_z[0] - 2*(-1./2.+i)*si;
//       }
//       for(int i=2;i<4;i++){
//         win_posi_x[i] = scint_posi_x[0] - si;
//         win_posi_y[i] = fScint_y/2. + d/2. + win_width/2.;
//         win_posi_z[i] = scint_posi_z[0] - 2*(5./2.-i)*si;
//       }

	
	
	// Collimator
	G4double colli_posi_x = x_posi;
	G4double colli_posi_y = -0.01-colli_width/2. -d - fScint_y/2.;//5.5*mm;
	G4double colli_posi_z = z_posi;
	// Wood1
	G4double wood1_posi_x = 0.*mm;
	G4double wood1_posi_y = fScint_y/2. + wood_width/2. + 20.*mm;
	G4double wood1_posi_z = 0.*mm;
	// Wood2
	G4double wood2_posi_x = 0.*mm;
	G4double wood2_posi_y = fScint_y/2.-wood_width/2. - 255.*mm;
	G4double wood2_posi_z = 0.*mm;
	// Fe1
	G4double Fe1_posi_x = 0.*mm;
	G4double Fe1_posi_y = -fScint_y/2. - fe_width/2. - 0.2*mm;
//	G4double Fe1_posi_y = -fScint_y/2.-fe_width/2.-colli_width-0.2*mm;
	G4double Fe1_posi_z = 0.*mm;
	// Fe2
	G4double Fe2_posi_x = 0.*mm;
	G4double Fe2_posi_y = fScint_y/2. + fe_width/2. - 255.*mm;
	G4double Fe2_posi_z = 0.*mm;
  

	// Source
	G4double source_posi_x = x_posi;
	G4double source_posi_y = wood1_posi_y + wood_width/2. + radi_case_width/2. + slide + 0.01;
//	G4double source_posi_y = -colli_width - radi_case_width/2. - fScint_y/2. - fe_width - d/2.;
//	G4double source_posi_y = -colli_width - radi_case_width/2. - fScint_y/2. - fe_width - d/2.-18.*mm;
//	G4double source_posi_y = wood1_posi_y + wood_width/2. + radi_case_width/2. + 0.01;//-0.02-d -fScint_y/2. -colli_width - radi_case_width/2.;//-colli_width + 5.5*mm + fScint_y/2. + fD_mtl + radi_case_width/2.; 
	G4double source_posi_z = z_posi;


	//Supporter
	G4double sup_posi_x = 0.*mm;
	G4double sup_posi_y = Fe2_posi_y + slide + sup_height/2.*mm;
	G4double sup_posi_z = fe_height/2. + 50. *mm;



	// center
	G4ThreeVector center_posi(0,0,0);



 	//Glice	
	G4ThreeVector glice_posi(X,0.,Z);



	// Collimator
	G4ThreeVector colli_posi(colli_posi_x, colli_posi_y, colli_posi_z);
	
	// Wood
	G4ThreeVector wood1_posi(wood1_posi_x, wood1_posi_y+slide, wood1_posi_z);
	G4ThreeVector wood2_posi(wood2_posi_x, wood2_posi_y+slide, wood2_posi_z);
	
	// Fe
	G4ThreeVector fe1_posi(Fe1_posi_x, Fe1_posi_y, Fe1_posi_z);
	G4ThreeVector fe2_posi(Fe2_posi_x, Fe2_posi_y+slide, Fe2_posi_z);
	
	// Source
	G4ThreeVector source_posi(source_posi_x, source_posi_y, source_posi_z);
	G4ThreeVector radi_posi(0.*mm, 0.*mm, -radi_case_width/3.);

	// Supporter
	G4ThreeVector sup_posi(sup_posi_x,sup_posi_y,sup_posi_z);


	G4RotationMatrix* rm_x1 = new G4RotationMatrix();
	rm_x1->rotateX(90.*deg);
	G4RotationMatrix* rm_x2 = new G4RotationMatrix();
	rm_x2->rotateX(180.*deg);
	G4RotationMatrix* rm_y1 = new G4RotationMatrix();
	rm_y1->rotateY(90.*deg);
	G4RotationMatrix* rm_y2 = new G4RotationMatrix();
	rm_y2->rotateY(180.*deg);
	G4RotationMatrix* rm_z1 = new G4RotationMatrix();
	rm_z1->rotateZ(90.*deg);
	G4RotationMatrix* rm_z2 = new G4RotationMatrix();
	rm_z2->rotateZ(180.*deg);
 
  //*************************** housing and scintillator
  fScint_box 		= new G4Box("scint_box",fScint_x/2.,fScint_y/2.,fScint_z/2.);
  fScintbig_box 	= new G4Box("scintbig_box",num*fScint_x/2.+(num-1)*d/2.,fScint_y/2.,num*fScint_z/2.+(num-1)*d/2.);
  fWood_box 		= new G4Box("wood_box",wood_height/2.,wood_width/2.,wood_height/2.);
//  fWood1_box 		= new G4Box("wood_box",wood_width/2.,wood_height/2. - wood_width/2.,wood_height/2. - wood_width/2.); //追加
//  fWood2_box 		= new G4Box("wood_box",wood_height/2. - wood_width/2.,wood_height/2. - wood_width/2.,wood_width/2.); //追加
  fFe_box 		= new G4Box("fe_box",fe_height/2.,fe_width/2.*mm,120./2.);
//  fGlice_box 		= new G4Box("glice_box",num*fScint_x/2.+(num-1.)*d/2.,d/2.+fScint_y/2.,num*fScint_z/2.+(num-1.)*d/2.);
  fGlice_box 		= new G4Box("glice_box",fScint_x/2.,d/2.+fScint_y/2.,num*fScint_z/2.+(num-1.)*d/2.);
  fGlicebig_box		= new G4Box("glicebig_box",num*fScint_x/2.+(num-1.)*d/2.,d/2.+fScint_y/2.,num*fScint_z/2.+(num-1.)*d/2.);
  fAir_box 	= new G4Box("air_box",air_x/2.,air_y/2.,air_z/2.);
  fAir_in_box 	= new G4Box("air_in_box",air_in_x/2.,air_in_y/2.,air_in_z/2.);
  fMPPC_box	= new G4Box("mppc_box",mppc_height/2.,mppc_height/2.,mppc_width/2.);  //rm_x1
//  fMPPC_box	= new G4Box("mppc_box",mppc_height/2.,mppc_width/2.,mppc_height/2.);
  fWin_box	= new G4Box("Win_box",win_longside/2.,win_shortside/2.,win_width/2.);  //rm_x1
//  fWin_box	= new G4Box("Win_box",win_longside/2.,win_shortside/2.,win_width/2.);
  fEnve_box	= new G4Box("Enve_box",enve_longside/2.,enve_shortside/2.,enve_width/2.);
  fCollimator	= new G4Cons("colli_con",colli_s_radius,colli_radius,colli_l_radius,colli_radius,colli_width/2.,0.*degree,360.*degree);
  fSup_box	= new G4Box("Sup_box",sup_yoko/2.,sup_height/2.,sup_width/2.);
  fRadiation	= new G4Tubs("radi_tub",0,radi_radius,radi_width/2.,0.*degree,360.*degree);
  fRadi_case	= new G4Tubs("radi_tub",0,radi_case_radius,radi_case_width/2.,0.*degree,360.*degree);

	fAir_log = new G4LogicalVolume(fAir_box,	G4Material::GetMaterial("Air"),	"scint_log", 0,0,0);
	fAir_in_log = new G4LogicalVolume(fAir_in_box,	G4Material::GetMaterial("Air"),	"scint_log", 0,0,0);
	fScint_log = new G4LogicalVolume(fScint_box,	G4Material::GetMaterial("GAGG"),	"scint_log", 0,0,0);
	fScintbig_log = new G4LogicalVolume(fScintbig_box,    G4Material::GetMaterial("GAGG"),        "scintbig_log", 0,0,0);
	fWood_log = new G4LogicalVolume(fWood_box,	G4Material::GetMaterial("Wood"),	"wood_log", 0,0,0);
//	fWood1_log = new G4LogicalVolume(fWood1_box,	G4Material::GetMaterial("Wood"),	"wood_log", 0,0,0);
//	fWood2_log = new G4LogicalVolume(fWood2_box,	G4Material::GetMaterial("Wood"),	"wood_log", 0,0,0);
	fFe_log = new G4LogicalVolume(fFe_box,	G4Material::GetMaterial("Al2"),	"fe_log", 0,0,0);
	fGlice_log = new G4LogicalVolume(fGlice_box,	G4Material::GetMaterial("Silica"),	"glice_log", 0,0,0);
	fGlicebig_log = new G4LogicalVolume(fGlicebig_box,	G4Material::GetMaterial("Silica"),	"glicebig_log", 0,0,0);
	fMPPC_log = new G4LogicalVolume(fMPPC_box,	G4Material::GetMaterial("Al"),"mppc_log", 0,0,0);
	fWin_log = new G4LogicalVolume(fWin_box,	G4Material::GetMaterial("Silicone"),"win_log",0,0,0);
	fEnve_log = new G4LogicalVolume(fEnve_box,	G4Material::GetMaterial("ABS"),	"enve_log",0,0,0);
	fCollimator_log = new G4LogicalVolume(fCollimator,	G4Material::GetMaterial("Pb"),	"colli_log", 0,0,0);
	fSup_log = new G4LogicalVolume(fSup_box,	G4Material::GetMaterial("Al2"),	"sup_log",0,0,0);
	fSource_log = new G4LogicalVolume(fRadi_case,	G4Material::GetMaterial("PMMA"),	"source_log", 0,0,0);
	fRadiation_log = new G4LogicalVolume(fRadiation,	G4Material::GetMaterial("Air"),	"radiation_log", 0,0,0);


  new G4PVPlacement(0,radi_posi,fRadiation_log,"radiation",fSource_log,false,0);
	
	//scintillator(複数個)
	//for(int i=0;i<scint_num;i++){
	//  fScint_vol[i] = new G4PVPlacement(0,{scint_posi_x[i],scint_posi_y[i],scint_posi_z[i]},fScint_log,"scintillator",fGlice_log,false,0);
	//}  

	//scintillator(大きいの)
	fScintbig_vol = new G4PVPlacement(0,{scintbig_posi_x,scintbig_posi_y,scintbig_posi_z},fScintbig_log,"scintillatorbig",fGlice_log,false,0);


	//glice(複数個)
	//fGlice_vol = new G4PVPlacement(0,glice_posi,fGlice_log,"glice",fAir_in_log,false,0);
//	//glice(大きいの)
	fGlicebig_vol = new G4PVPlacement(0,glice_posi,fGlicebig_log,"glice",fAir_in_log,false,0);
	
	
	//collimator
	//if(colli)
//	new G4PVPlacement(rm_x1,colli_posi,fCollimator_log,"collimator",fAir_in_log,false,0);
	new G4PVPlacement(0,wood1_posi,fWood_log,"Wood",fAir_in_log,false,0);
	new G4PVPlacement(0,wood2_posi,fWood_log,"wood",fAir_in_log,false,0);

	new G4PVPlacement(0,fe1_posi,fFe_log,"fe",fAir_in_log,false,0);  //no colli
	new G4PVPlacement(0,fe2_posi,fFe_log,"fe",fAir_in_log,false,0);


	//Supporter
////	new G4PVPlacement(0,sup_posi,fSup_log,"supporter",fAir_in_log,false,0);


	//radiation source
	new G4PVPlacement(rm_x1,source_posi,fSource_log,"source",fAir_in_log,false,0);
  
  //****************** Build PMTs
  //the "photocathode" is a metal slab at the back of the epoxy that
  //is only a very rough approximation of the real thing since it only
  //absorbs or detects the photons based on the efficiency set below
  
	fMPPC_vol = new G4PVPlacement(0,{mppc_posi_x[0],mppc_posi_y[0],mppc_posi_z[0]},fMPPC_log,"MPPC",fWin_log,false,0);
	for(int i=0;i<mppc_num;i++){	
		fEnve_vol[i] = new G4PVPlacement(rm_x1,{enve_posi_x[i],enve_posi_y[i],enve_posi_z[i]},fEnve_log,"Envelope",fAir_in_log,false,0);
		fWin_vol[i] = new G4PVPlacement(rm_x1,{win_posi_x[i],win_posi_y[i],win_posi_z[i]},fWin_log,"Window",fAir_in_log,false,0);
	}
      		// rm_x1
	fAir_vol = new G4PVPlacement(0,center_posi,fAir_in_log,"airin",fAir_log,false,8);
 

 VisAttributes();
  SurfaceProperties();

  SetLogicalVolume(fAir_log);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeMainVolume::CopyValues(){
  fScint_x=fConstructor->GetScintX();
  fScint_y=fConstructor->GetScintY();
  fScint_z=fConstructor->GetScintZ();
	fD_mtl=fConstructor->GetHousingThickness();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeMainVolume::VisAttributes(){
	//G4VisAttributes* scint_va = new G4VisAttributes(G4Colour(0.8,0.8,0.2));
	//fScint_log->SetVisAttributes(scint_va);
	G4VisAttributes* scintbig_va = new G4VisAttributes(G4Colour(0.8,0.8,0.2));
	fScintbig_log->SetVisAttributes(scintbig_va);
	G4VisAttributes* mppc_va = new G4VisAttributes(G4Colour(0.2,0.4,0.2));
	fMPPC_log->SetVisAttributes(mppc_va);
	G4VisAttributes* win_va = new G4VisAttributes(G4Colour(0.4,0.4,0.6));
	fWin_log->SetVisAttributes(win_va);
	G4VisAttributes* enve_va = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
	fEnve_log->SetVisAttributes(enve_va);
	G4VisAttributes* colli_va = new G4VisAttributes(G4Colour(0.4,0.4,0.6));
	fCollimator_log->SetVisAttributes(colli_va);
	G4VisAttributes* radi_va = new G4VisAttributes(G4Colour(0.8,0.8,0.8));
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

void LXeMainVolume::SurfaceProperties(){
  G4double ephoton[] = {2.76*eV, 3.06*eV};
  const G4int num = sizeof(ephoton)/sizeof(G4double);



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

  //**GAGGの反射設定 (GAGG周り)
	G4double GAGG_RIND[] = {1.0,1.0};
	G4double GAGG_REF[] = {1.0,1.0};
	G4double GAGG2_RIND[] = {1.0,1.0};
	G4double GAGG2_REF[] = {1.0,1.0};

	G4MaterialPropertiesTable* GAGG_PT = new G4MaterialPropertiesTable();
	G4MaterialPropertiesTable* GAGG2_PT = new G4MaterialPropertiesTable();
	
	GAGG_PT->AddProperty("RINDEX",ephoton,GAGG_RIND,num);
	GAGG_PT->AddProperty("REFLECTIVITY",ephoton,GAGG_REF,num);
	GAGG2_PT->AddProperty("RINDEX",ephoton,GAGG2_RIND,num);
	GAGG2_PT->AddProperty("REFLECTIVITY",ephoton,GAGG2_REF,num);

  G4OpticalSurface* GAGG_Surface = new G4OpticalSurface("GAGGSurface",glisur,groundbackpainted,dielectric_dielectric);
	GAGG_Surface->SetMaterialPropertiesTable(GAGG_PT);
  G4OpticalSurface* GAGG2_Surface = new G4OpticalSurface("GAGG2Surface",glisur,polished,dielectric_dielectric);
	GAGG2_Surface->SetMaterialPropertiesTable(GAGG2_PT);

  for(int i=0;i<scint_num;i++){
  new G4LogicalBorderSurface("GAGG_surf",fScint_vol[i],fAir_vol,GAGG_Surface);
	}
  
  new G4LogicalBorderSurface("GAGG_surf",fScintbig_vol,fAir_vol,GAGG_Surface);

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
 
  //**Gliceの反射設定 (GSO周り)
	G4double Glice_RIND[] = {1.0,1.0};
	G4double Glice_REF[] = {1.0,1.0};
	G4MaterialPropertiesTable* Glice_PT = new G4MaterialPropertiesTable();
	Glice_PT->AddProperty("RINDEX",ephoton,Glice_RIND,num);
	Glice_PT->AddProperty("REFLECTIVITY",ephoton,Glice_REF,num);
  G4OpticalSurface* Glice_Surface =
    new G4OpticalSurface("GliceSurface",glisur,groundbackpainted,dielectric_dielectric);
	Glice_Surface->SetMaterialPropertiesTable(Glice_PT);

  new G4LogicalBorderSurface("Glice_surf",fGlice_vol,fAir_vol,Glice_Surface);

  //**Windowの反射設定 (Envelope周り)
         G4double Window_RIND[] = {1.0,1.0};
         G4double Window_REF[] = {1.0,1.0};
         G4MaterialPropertiesTable* Window_PT = new G4MaterialPropertiesTable();
         Window_PT->AddProperty("RINDEX",ephoton,Window_RIND,num);
         Window_PT->AddProperty("REFLECTIVITY",ephoton,Window_REF,num);
   G4OpticalSurface* Window_Surface =
     new G4OpticalSurface("WindowSurface",glisur,groundbackpainted,dielectric_dielectric);
         Window_Surface->SetMaterialPropertiesTable(Window_PT);

	for(int i=0;i<mppc_num;i++){
	   new G4LogicalBorderSurface("Window_surf",fWin_vol[i],fEnve_vol[i],Window_Surface);
	   new G4LogicalBorderSurface("Window_surf",fWin_vol[i],fAir_vol,Window_Surface);
//		for(int j=0;j<mppc_num;j++){
//	   		new G4LogicalBorderSurface("Window_surf",fWin_vol[i],fWin_vol[j],Window_Surface);
//		}
	}




  //**mppcの反射設定(光子感度) 
	G4double MPPC_ReR[]={1.92,1.92};
	G4double MPPC_ImR[]={1.69,1.69};
	G4double MPPC_EFF[]={1.0,1.0}; //enables 'detection' of photons
	G4MaterialPropertiesTable* MPPC_mt = new G4MaterialPropertiesTable();
	MPPC_mt->AddProperty("REALRINDEX",ephoton,MPPC_ReR,num);
	MPPC_mt->AddProperty("IMAGINARYRINDEX",ephoton,MPPC_ImR,num);
	MPPC_mt->AddProperty("EFFICIENCY",ephoton,MPPC_EFF,num);
	G4OpticalSurface* MPPC_Surface = 
	  new G4OpticalSurface("MPPCSurface",glisur,polished,dielectric_metal);
	MPPC_Surface->SetMaterialPropertiesTable(MPPC_mt);



  //**Create logical skin surfaces
	//new G4LogicalSkinSurface("GSO_surf",fScint_log,GSO_Surface);
	new G4LogicalSkinSurface("mppc_surf",fMPPC_log,MPPC_Surface);
}
