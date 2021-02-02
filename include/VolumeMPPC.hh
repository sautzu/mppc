/**
 * MPPC Logical Volume
 **/

#include "G4PVPlacement.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Material.hh"

class VolumeMPPC
{
public:
    VolumeMPPC();
    //LogicalVolumeを返す
    G4LogicalVolume *getLogicalVolume();
    //可視化
    void VisAttributes();
    //境界の反射
    void SurfaceProperties(G4VPhysicalVolume*);
    //測定部分のLogicalVolumeを返す
    G4LogicalVolume *getMPPCLogical();

private:
    //Solid
    G4Box *fMPPC_box;
    G4Box *fWin_box;
    G4Box *fEnve_box;

    //LogicalVolume
    G4LogicalVolume *fMPPC_log;
    G4LogicalVolume *fWin_log;
    G4LogicalVolume *fEnve_log;

    //PhysicalVolume
    G4PVPlacement *fMPPC_phy;
    G4PVPlacement *fWin_phy;
    G4PVPlacement *fEnve_phy;

    //MPPC
    G4double mppc_height;
    G4double mppc_width;
    // window
    G4double win_longside;
    G4double win_shortside;
    G4double win_width;
    //Envelop
    G4double enve_longside;
    G4double enve_shortside;
    G4double enve_width;

    G4VisAttributes *mppc_va;
    G4VisAttributes *win_va;
    G4VisAttributes *enve_va;
};
