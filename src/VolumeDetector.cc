#include "VolumeDetector.hh"
#include "G4VisAttributes.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4SystemOfUnits.hh"

VolumeDetector::VolumeDetector()
{
    d = 0.2 * mm;
    fmother_Sol = new G4Box("Detector_Box",10. * mm, 5. * mm, 10. * mm);

    fmother_log = new G4LogicalVolume(fmother_Sol, G4Material::GetMaterial("Air"), "Detector_mother", 0, 0, 0);
    mppc.reset(new VolumeMPPC());
    scint.reset(new VolumeScint());

    for (auto i = 0; i < 4;i++){
        mppc_phy.push_back(new G4PVPlacement(0, {0, 0, i * 10. * mm}, mppc->getLogicalVolume(), "MPPC_in_Detector", fmother_log, false, 1 + i));
    }

    for (auto i = 0; i < 4;i++){
        scint_phy.push_back(new G4PVPlacement(0, {0, -1.5, i * 10. * mm}, scint->getLogicalVolume(), "MPPC_in_Detector", fmother_log, false, 1 + i));
    }
}

G4LogicalVolume *VolumeDetector::getLogicalVolume()
{
    return fmother_log;
}

void VolumeDetector::VisAttributes()
{
    mppc->VisAttributes();
    scint->VisAttributes();
}

void VolumeDetector::SurfaceProperties(G4VPhysicalVolume *fAir_phy)
{
    mppc->SurfaceProperties();
    scint->SurfaceProperties(fAir_phy);
}

G4LogicalVolume *VolumeDetector::getScintLogical()
{
    return scint->getScintLogical();
}
G4LogicalVolume *VolumeDetector::getMPPCLogical()
{
    return mppc->getMPPCLogical();
}
