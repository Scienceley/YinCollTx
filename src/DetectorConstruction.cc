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
//
/// \file B1/src/DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4Material.hh"
#include "G4SubtractionSolid.hh"
#include "G4MultiUnion.hh"
#include "iostream"


namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters
  //
  G4double env_sizeXY = 50 * cm, env_sizeZ = 140 * cm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_AIR");

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //
  // World
  //
  G4double world_sizeXY = 1.2 * env_sizeXY;
  G4double world_sizeZ = 1.2 * env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

  auto solidWorld =
    new G4Box("World",  // its name
              0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ);  // its size

  auto logicWorld = new G4LogicalVolume(solidWorld,  // its solid
                                        world_mat,  // its material
                                        "World");  // its name

  auto physWorld = new G4PVPlacement(nullptr,  // no rotation
                                     G4ThreeVector(),  // at (0,0,0)
                                     logicWorld,  // its logical volume
                                     "World",  // its name
                                     nullptr,  // its mother  volume
                                     false,  // no boolean operation
                                     0,  // copy number
                                     checkOverlaps);  // overlaps checking

  //
  // Envelope
  //
  auto solidEnv = new G4Box("Envelope",  // its name
                            0.5 * env_sizeXY, 0.5 * env_sizeXY, 0.5 * env_sizeZ);  // its size

  auto logicEnv = new G4LogicalVolume(solidEnv,  // its solid
                                      env_mat,  // its material
                                      "Envelope");  // its name

  new G4PVPlacement(nullptr,  // no rotation
                    G4ThreeVector(),  // at (0,0,0)
                    logicEnv,  // its logical volume
                    "Envelope",  // its name
                    logicWorld,  // its mother  volume
                    false,  // no boolean operation
                    0,  // copy number
                    checkOverlaps);  // overlaps checking
 
  
  // Shape 1
  // Square
  //
  G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_A-150_TISSUE");
  G4ThreeVector pos1 = G4ThreeVector(0, 0 * cm, 50 * cm);
  
  // Cube section shape
  G4double shape1_x = 16.0 *cm;
  G4double shape1_y = 16.0 *cm;
  G4double shape1_z = 16.0 *cm;
    
  auto solidShape1 = new G4Box("Shape1", shape1_x *0.5, shape1_y *0.5, shape1_z *0.5);

  auto logicShape1 = new G4LogicalVolume(solidShape1,  // its solid
                                         shape1_mat,  // its material
                                         "Shape1");  // its name

  new G4PVPlacement(nullptr,  // no rotation
                    pos1,  // at position
                    logicShape1,  // its logical volume
                    "Shape1",  // its name
                    logicEnv,  // its mother  volume
                    false,  // no boolean operation
                    0,  // copy number
                    checkOverlaps);  // overlaps checking
  fScoringVolume = logicShape1;
  
  
  // Collimator
  G4double CollimatorxPos = 0  *cm;
  G4double CollimatoryPos = 0 *cm;
  G4double CollimatorzPos = 0 *cm;
  G4Material* shapeC_mat = nist->FindOrBuildMaterial("G4_Cu"); 
  G4ThreeVector posC = G4ThreeVector(CollimatorxPos, CollimatoryPos, CollimatorzPos);
  
  G4double NoSlits = 5;  //Use low numbers as a test
  G4double SlitxSize = 0.04  *cm;	//0.04
  G4double SlitySize = 5 *cm;
  G4double SlitzSize = 18 *cm;
  
  G4double BoxxSize = 20  *cm;
  G4double BoxySize = 16 *cm;
  G4double BoxzSize = 6.5 *cm;
  
  G4double c2t2cDist = 0.4*cm;//0.4
  
  G4Box* box = new G4Box("Box",BoxxSize*0.5,BoxySize*0.5,BoxzSize*0.5); 				//create box  
  G4Box* Slit = new G4Box("Slit",SlitxSize*0.5,SlitySize*0.5,SlitzSize*0.5);  				//Create the slit object
  
  G4ThreeVector First_slit_pos = G4ThreeVector(-(NoSlits+1)/2 *c2t2cDist, 0, 0); 				//Position of first slit
  G4double First_slit_xpos = G4double((NoSlits+1)/2 *c2t2cDist); 						//x Position of first slit
  //G4double DistBetween =  (BoxxSize - 2*cm) / (NoSlits - 1);						//distance between each slit	

  G4MultiUnion* munion_solid = new G4MultiUnion("Boxes_Union");						// Initialise a MultiUnion structure
    
  std::cout << "Test";  
    
  if (  c2t2cDist > NoSlits / SlitxSize){
  
  std::cout <<"Error, too many slits of this size. \n "; 
  std::cout << "Make slits smaller or decrease amount";
  
  }
    
    
  for (int i=0; i<NoSlits; i++){ 						//For loop
  
  //
  
  // Translation matrix
  G4RotationMatrix rotm  = G4RotationMatrix();  
  G4ThreeVector position1 = G4ThreeVector(First_slit_xpos-((i+1)*c2t2cDist),0*cm,0*cm);
  G4Transform3D tr1 = G4Transform3D(rotm,position1);
  
  // Add the shapes to the structure
  munion_solid->AddNode(Slit,tr1);
  
  // Finally close the structure
  munion_solid->Voxelize();             
       

  } //end loop
  
  
  // Associate it to a solid volume, then a logical volume 

  //Create subtraction shape
  auto solidShapeC = new G4SubtractionSolid("Box+Slit", box, munion_solid); 	
	
	  
	auto logicShapeC = new G4LogicalVolume(solidShapeC,  // its solid
                                         shapeC_mat,  // its material
                                         "Box+Slit");  // its name

  	auto placement = new G4PVPlacement(nullptr,  // no rotation
                    posC,  // at position
                    logicShapeC,  // its logical volume
                    "Box+Slit",  // its name
                    logicEnv,  // its mother  volume
                    false,  // no boolean operation
                    0,  // copy number
                    checkOverlaps);  // overlaps checking
  
                    
  
  
  // always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B1




/*  Cylinders below
  // Shape 1
  // Inner Cylinder. Tumour
  //
  G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_A-150_TISSUE");
  G4ThreeVector pos1 = G4ThreeVector(0, 0 * cm, 0 * cm);
  G4RotationMatrix* xRot1 = new G4RotationMatrix;  // Rotates y and Z axes only
  xRot1->rotateX(M_PI/2*rad);                     // Rotates 90 degrees

  // Inner Cylinder section shape
  G4double shape1_rmin = 0.0 *cm;
  G4double shape1_rmax = 3.0 *cm;
  G4double shape1_pDz = 16.0 *cm;
  G4double shape1_sphi = 0.0;
  G4double shape1_phi = 2.0 * CLHEP::pi;
  
  auto solidShape1 = new G4Tubs("Shape1", shape1_rmin *0.5, shape1_rmax *0.5, shape1_pDz *0.5, shape1_sphi = 0,
                                shape1_phi);

  auto logicShape1 = new G4LogicalVolume(solidShape1,  // its solid
                                         shape1_mat,  // its material
                                         "Shape1");  // its name

  new G4PVPlacement(xRot1,  // no rotation
                    pos1,  // at position
                    logicShape1,  // its logical volume
                    "Shape1",  // its name
                    logicEnv,  // its mother  volume
                    false,  // no boolean operation
                    0,  // copy number
                    checkOverlaps);  // overlaps checking
  
  // Set Shape1 (Inner Cylinder) as scoring volume
  //
  fScoringVolume = logicShape1;
  
  
  // Second layer Cylinder. Brain 
  // 
  //
  G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_A-150_TISSUE");
  G4ThreeVector pos2 = G4ThreeVector(0, 0 * cm, 0 * cm);

  // 2nd Cylinder section shape
  G4double shape2_rmin = 3.0;
  G4double shape2_rmax = 7.0 *cm;
  G4double shape2_pDz = 14.0 *cm;
  G4double shape2_sphi = 0.0;
  G4double shape2_phi = 2.0 * CLHEP::pi;
  
  // 2nd cylinder rotation
  G4RotationMatrix* xRot2 = new G4RotationMatrix;  // Rotates y and Z axes only
  xRot2->rotateX(M_PI/2*rad);                     // Rotates 90 degrees
  
  auto solidShape2 = new G4Tubs("Shape2", shape2_rmin, shape2_rmax *0.5, shape2_pDz *0.5, shape2_sphi *0.5,
                                shape2_phi);

  auto logicShape2 = new G4LogicalVolume(solidShape2,  // its solid
                                         shape2_mat,  // its material
                                         "Shape2");  // its name

  new G4PVPlacement(xRot2,  // no rotation
                    pos1,  // at position
                    logicShape2,  // its logical volume
                    "Shape2",  // its name
                    logicEnv,  // its mother  volume
                    false,  // no boolean operation
                    0,  // copy number
                    checkOverlaps);  // overlaps checking
                    
  
                    
                    
  // Third layer Cylinder. Bone / skull 
  // 
  //
  G4Material* shape3_mat = nist->FindOrBuildMaterial("G4_BONE_COMPACT_ICRU");
  G4ThreeVector pos3 = G4ThreeVector(0, 0 * cm, 0 * cm);

  // 3rd Cylinder section shape
  G4double shape3_rmin = 7.0;
  G4double shape3_rmax = 8.0 *cm;
  G4double shape3_pDz = 12.0 *cm;
  G4double shape3_sphi = 0.0;
  G4double shape3_phi = 2.0 * CLHEP::pi;
  
  // 3rd cylinder rotation
  G4RotationMatrix* xRot3 = new G4RotationMatrix;  // Rotates y and Z axes only
  xRot3->rotateX(M_PI/2*rad);                     // Rotates 90 degrees
  
  auto solidShape3 = new G4Tubs("Shape3", shape3_rmin, shape3_rmax *0.5, shape3_pDz *0.5, shape3_sphi *0.5,
                                shape3_phi);

  auto logicShape3 = new G4LogicalVolume(solidShape3,  // its solid
                                         shape3_mat,  // its material
                                         "Shape3");  // its name

  new G4PVPlacement(xRot3,  // no rotation
                    pos1,  // at position
                    logicShape3,  // its logical volume
                    "Shape3",  // its name
                    logicEnv,  // its mother  volume
                    false,  // no boolean operation
                    0,  // copy number
                    checkOverlaps);  // overlaps checking               
  

  // Fourth layer Cylinder. Skin
  // 
  //
  G4Material* shape4_mat = nist->FindOrBuildMaterial("G4_A-150_TISSUE");
  G4ThreeVector pos4 = G4ThreeVector(0, 0 * cm, 0 * cm);

  // 4th cylinder section shape
  G4double shape4_rmin = 8.0;
  G4double shape4_rmax = 8.5 *cm;
  G4double shape4_pDz = 10.0 *cm;
  G4double shape4_sphi = 0.0;
  G4double shape4_phi = 2.0 * CLHEP::pi;
  
  // 4th cylinder rotation
  G4RotationMatrix* xRot4 = new G4RotationMatrix;  // Rotates y and Z axes only
  xRot4->rotateX(M_PI/2*rad);                     // Rotates 90 degrees
  
  auto solidShape4 = new G4Tubs("Shape4", shape4_rmin, shape4_rmax *0.5, shape4_pDz *0.5, shape4_sphi *0.5,
                                shape4_phi);

  auto logicShape4 = new G4LogicalVolume(solidShape4,  // its solid
                                         shape4_mat,  // its material
                                         "Shape4");  // its name

  new G4PVPlacement(xRot4,  // no rotation
                    pos1,  // at position
                    logicShape4,  // its logical volume
                    "Shape4",  // its name
                    logicEnv,  // its mother  volume
                    false,  // no boolean operation
                    0,  // copy number
                    checkOverlaps);  // overlaps checking
  */
