#include "VolumeScint.hh"
#include "G4VisAttributes.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4SystemOfUnits.hh"

VolumeScint::VolumeScint()
{

    fScint_x = 20.0 * mm;
    fScint_y = 5.0 * mm;
    fScint_z = 20.0 * mm;

    fScint_Sol = new G4Box("Scint_box", fScint_x / 2, fScint_y / 2, fScint_z / 2);
    fmother_Sol = new G4Box("Scint_mother", fScint_x / 2, fScint_y / 2, fScint_z / 2);

    fScint_log = new G4LogicalVolume(fScint_Sol, G4Material::GetMaterial("GAGG"), "Scint", 0, 0, 0);
    fmother_log = new G4LogicalVolume(fmother_Sol, G4Material::GetMaterial("Air"), "Scint_mother", 0, 0, 0);

    fScint_phy = new G4PVPlacement(0, G4ThreeVector(), "Scint", fmother_log, 0, false, 0);
}

G4LogicalVolume *VolumeScint::getLogicalVolume()
{
    return fmother_log;
}

void VolumeScint::VisAttributes()
{
    Scint_va = new G4VisAttributes(G4Colour(0.8, 0.8, 0.2));
    fScint_log->SetVisAttributes(Scint_va);
}

void VolumeScint::SurfaceProperties(G4VPhysicalVolume *fAir_phy)
{
    //エネルギーごとに光学的なパラメータを指定できる
    G4double ephoton[] = {2.76 * eV, 3.06 * eV};
    const G4int num = sizeof(ephoton) / sizeof(G4double);

    //屈折率
    G4double GAGG_RIND[] = {1.0, 1.0};
    //反射率
    G4double GAGG_REF[] = {1.0, 1.0};

    G4MaterialPropertiesTable *GAGG_PT = new G4MaterialPropertiesTable();

    GAGG_PT->AddProperty("RINDEX", ephoton, GAGG_RIND, num);
    GAGG_PT->AddProperty("REFLECTIVITY", ephoton, GAGG_REF, num);

    G4OpticalSurface *GAGG_Surface = new G4OpticalSurface("GAGGSurface", glisur, groundbackpainted, dielectric_dielectric);
    GAGG_Surface->SetMaterialPropertiesTable(GAGG_PT);

    //空気とシンチレータの境界面
    new G4LogicalBorderSurface("GAGG_surf", fScint_phy, fAir_phy, GAGG_Surface);
}

G4LogicalVolume *VolumeScint::getScintLogical()
{
    return fScint_log;
}
