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
/// \file B1/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the B1::PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);

  // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle(particleName = "gamma");
  fParticleGun->SetParticleDefinition(particle);
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
  //fParticleGun->SetParticleEnergy(150 * MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  // this function is called at the begining of ecah event
  //

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get Envelope volume
  // from G4LogicalVolumeStore.

  G4double envSizeXY = 0;
  G4double envSizeZ = 0;

  if (!fEnvelopeBox) {
    G4LogicalVolume* envLV = G4LogicalVolumeStore::GetInstance()->GetVolume("Envelope");
    if (envLV) fEnvelopeBox = dynamic_cast<G4Box*>(envLV->GetSolid());
  }

  if (fEnvelopeBox) {
    envSizeXY = fEnvelopeBox->GetXHalfLength() * 2.;
    envSizeZ = fEnvelopeBox->GetZHalfLength() * 2.;
  }
  else {
    G4ExceptionDescription msg;
    msg << "Envelope volume of box shape not found.\n";
    msg << "Perhaps you have changed geometry.\n";
    msg << "The gun will be place at the center.";
    G4Exception("PrimaryGeneratorAction::GeneratePrimaries()", "MyCode0002", JustWarning, msg);
  }
  
  /*
 
  
  //Angle dist
  G4double theta, phi;
  G4double px, py, pz;
  //G4double sigmaAngle = 1.;//degrees
  //theta = G4RandGauss::shoot(0.0, (sigmaAngle*CLHEP::pi) / 180.); //for degrees
  G4double sigmaAngle = 0.007; //radians
  theta = G4RandGauss::shoot(0.0,sigmaAngle); //for radians
  phi = CLHEP::twopi * G4UniformRand();
  px = std::sin(theta) * std::cos(phi);
  py = std::sin(theta) * std::sin(phi);
  pz = std::cos(theta);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px, py, pz));
  //fParticleGun->SetParticlePosition(G4ThreeVector(0, 0, 0));
  //Energy
  G4double beamEnergy = 50*MeV;
  G4double energySigma = 0.1;
  G4double finalEnergy = G4RandGauss::shoot(beamEnergy, 
  beamEnergy*energySigma);
  fParticleGun->SetParticleEnergy(finalEnergy);


//  fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
  fParticleGun->SetParticlePosition(G4ThreeVector(0 *cm, 0 *cm, -50 *cm));
  fParticleGun->GeneratePrimaryVertex(event);

  //fParticleGun->GeneratePrimaryVertex(event);
  
  */
  
  G4double startX = 0.;
  G4double startY = 0.;
  G4double startZ = -50.*cm;
  G4double sizeX = 15.*cm;
  G4double sizeY = 15.*cm;
  //Generate distribution from -0.5 and 0.5 and scale by size of beam
  startX = (sizeX * (G4UniformRand() - 0.5)); 
  startY = (sizeY * (G4UniformRand() - 0.5));
  fParticleGun->SetParticlePosition(G4ThreeVector(startX, startY, startZ));
  
  
  G4double theta, phi;
  G4double px, py, pz;
  G4double sigmaAngle = 0.007; //radians
  theta = G4RandGauss::shoot(0.0,sigmaAngle); //for radians
  phi = CLHEP::twopi * G4UniformRand();
  px = std::sin(theta) * std::cos(phi);
  py = std::sin(theta) * std::sin(phi);
  pz = std::cos(theta);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px, py, pz));
  
  //Energy
  G4double beamEnergy = 100*keV;
  G4double energySigma = 0.01;
  G4double finalEnergy = G4RandGauss::shoot(beamEnergy, 
  beamEnergy*energySigma);
  fParticleGun->SetParticleEnergy(finalEnergy);

  fParticleGun->GeneratePrimaryVertex(event);
  
  //  /control/execute run2.mac

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B1
