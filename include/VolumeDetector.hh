/**
 * Constructed Detector Logical Volume
 **/

#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Material.hh"
#include "VolumeMPPC.hh"
#include "VolumeScint.hh"

#include <memory>
#include <vector>

class VolumeDetector
{
public:
    VolumeDetector();
    //LogicalVolumeを返す
    G4LogicalVolume *getLogicalVolume();
    //可視化
    void VisAttributes();
    //境界の反射
    void SurfaceProperties(G4VPhysicalVolume *fAir_phy);
    //測定部分のLogicalVolumeを返す
    G4LogicalVolume *getScintLogical();
    G4LogicalVolume *getMPPCLogical();

private:
    G4double d; //グリスの厚さ
    G4Box *fmother_Sol;
    std::unique_ptr<VolumeMPPC> mppc;
    std::unique_ptr<VolumeScint> scint;

    //LogicalVolume
    G4LogicalVolume *fmother_log;

    //PhysicalVolume
    std::vector<G4PVPlacement*> mppc_phy;
    std::vector<G4PVPlacement*> scint_phy;

    G4VisAttributes *Scint_va;
};
