/**
 * Scintillator Logical Volume
 **/

#include "G4PVPlacement.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Material.hh"

class VolumeScint
{
public:
    VolumeScint();
    //LogicalVolumeを返す
    G4LogicalVolume *getLogicalVolume();
    //可視化
    void VisAttributes();
    //境界の反射
    void SurfaceProperties(G4VPhysicalVolume *fAir_phy,G4VPhysicalVolume *fScint_phy);
    //測定部分のLogicalVolumeを返す
    G4LogicalVolume *getScintLogical();

private:
    //Solid
    G4Box *fScint_Sol;
    G4Box *fmother_Sol;

    //LogicalVolume
    G4LogicalVolume *fScint_log;
    G4LogicalVolume *fmother_log;

    //PhysicalVolume
    G4VPhysicalVolume *fScint_phy;

    //size
    G4double fScint_x;
    G4double fScint_y;
    G4double fScint_z;

    G4VisAttributes *Scint_va;
};
