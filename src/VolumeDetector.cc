#include "VolumeDetector.hh"
#include "G4VisAttributes.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4SystemOfUnits.hh"

VolumeDetector::VolumeDetector()
{
    d = 0.2 * mm;
    guide_width = 2.0 * mm;
    fmother_Sol = new G4Box("Detector_Box", 100. * mm / 2., 80. * mm / 2., 100. * mm / 2.);
    fGlice_Sol1 = new G4Box("Glice_Box1", 20. * mm / 2., 5. * mm / 2., d / 2.);
    fGlice_Sol2 = new G4Box("Glice_Box2", 6. * mm / 2., d / 2., 6. * mm / 2.);
    fGuide_Sol = new G4Box("Light_Guide", 20. * mm / 2, 5. * mm / 2, guide_width / 2);

    fmother_log = new G4LogicalVolume(fmother_Sol, G4Material::GetMaterial("Air"), "Detector_mother", 0, 0, 0);
    fAir_log = new G4LogicalVolume(fmother_Sol, G4Material::GetMaterial("Air"), "Detector_in_Air", 0, 0, 0);
    fGlice_log1 = new G4LogicalVolume(fGlice_Sol1, G4Material::GetMaterial("Silica"), "Glice1", 0, 0, 0);
    fGlice_log2 = new G4LogicalVolume(fGlice_Sol2, G4Material::GetMaterial("Silica"), "Glice2", 0, 0, 0);
    fGuide_log = new G4LogicalVolume(fGuide_Sol, G4Material::GetMaterial("PMMA"), "Light_Guide", 0, 0, 0);
    mppc.reset(new VolumeMPPC());
    scint.reset(new VolumeScint());

    G4RotationMatrix *rm_mppc = new G4RotationMatrix();
    rm_mppc->rotateX(-90. * deg);
    //MPPCの配置
    mppc_phy.push_back(new G4PVPlacement(rm_mppc, {0, -33.95 * mm + d, -20. * mm - 2 * d - guide_width}, mppc->getLogicalVolume(), "MPPC_in_Detector1", fAir_log, false, 1));
    mppc_phy.push_back(new G4PVPlacement(rm_mppc, {0, -33.95 * mm + d, 0. * mm}, mppc->getLogicalVolume(), "MPPC_in_Detector2", fAir_log, false, 2));
    mppc_phy.push_back(new G4PVPlacement(rm_mppc, {0, -33.95 * mm + d, 20. * mm + 2 * d + guide_width}, mppc->getLogicalVolume(), "MPPC_in_Detector3", fAir_log, false, 3));

    //グリスの配置
    glice_phy1.push_back(new G4PVPlacement(0, {0, -37.5 * mm, -10. * mm - d / 2.}, fGlice_log1, "glice_between_scint1", fAir_log, false, 1));
    glice_phy1.push_back(new G4PVPlacement(0, {0, -37.5 * mm, -10. * mm - guide_width - d - d / 2.}, fGlice_log1, "glice_between_scint2", fAir_log, false, 2));
    glice_phy1.push_back(new G4PVPlacement(0, {0, -37.5 * mm, 10. * mm + d / 2.}, fGlice_log1, "glice_between_scint3", fAir_log, false, 3));
    glice_phy1.push_back(new G4PVPlacement(0, {0, -37.5 * mm, 10. * mm + guide_width + d + d / 2}, fGlice_log1, "glice_between_scint4", fAir_log, false, 4));
    glice_phy2.push_back(new G4PVPlacement(0, {0, -35. * mm + d / 2, -20. * mm - 2 * d - guide_width}, fGlice_log2, "glice_mppc1", fAir_log, false, 1));
    glice_phy2.push_back(new G4PVPlacement(0, {0, -35. * mm + d / 2, 0. * mm}, fGlice_log2, "glice_mppc2", fAir_log, false, 2));
    glice_phy2.push_back(new G4PVPlacement(0, {0, -35. * mm + d / 2, 20. * mm + 2 * d + guide_width}, fGlice_log2, "glice_mppc3", fAir_log, false, 3));

    //ライトガイドの配置
    guide_phy.push_back(new G4PVPlacement(0, {0, -37.5 * mm, -10. * mm - d - guide_width / 2}, fGuide_log, "LightGUide1", fAir_log, false, 1));
    guide_phy.push_back(new G4PVPlacement(0, {0, -37.5 * mm, 10. * mm + d + guide_width / 2}, fGuide_log, "LightGUide2", fAir_log, false, 1));

    //シンチレーターの配置
    scint_phy.push_back(new G4PVPlacement(0, {0, -37.5 * mm, -20. * mm - 2 * d - guide_width}, scint->getLogicalVolume(), "MPPC_in_Detector1", fAir_log, false, 1));
    scint_phy.push_back(new G4PVPlacement(0, {0, -37.5 * mm, 0 * mm}, scint->getLogicalVolume(), "MPPC_in_Detector2", fAir_log, false, 2));
    scint_phy.push_back(new G4PVPlacement(0, {0, -37.5 * mm, 20. * mm + 2 * d + guide_width}, scint->getLogicalVolume(), "MPPC_in_Detector3", fAir_log, false, 3));

    fAir_phy = new G4PVPlacement(0, G4ThreeVector(), fAir_log, "Detector_in_Air", fmother_log, false, 0);
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
    mppc->SurfaceProperties(this->fAir_phy);
    for (auto i : scint_phy)
    {
        scint->SurfaceProperties(fAir_phy, i);
        scint->SurfaceProperties(this->fAir_phy, i);
    }

    G4double ephoton[] = {2.76 * eV, 3.06 * eV};
    const G4int num = sizeof(ephoton) / sizeof(G4double);

    //**Gliceの反射設定
    G4double Glice_RIND[] = {1.0, 1.0};
    G4double Glice_REF[] = {1.0, 1.0};
    G4double Glice_LOBE[] = {0., 0.};
    G4double Glice_SPIKE[] = {0.95, 0.95};
    G4double Glice_SCATTER[] = {0., 0.};
    G4MaterialPropertiesTable *Glice_PT = new G4MaterialPropertiesTable();
    Glice_PT->AddProperty("RINDEX", ephoton, Glice_RIND, num);
    Glice_PT->AddProperty("REFLECTIVITY", ephoton, Glice_REF, num);
    Glice_PT->AddProperty("SPECULARLOBECONSTANT", ephoton, Glice_LOBE, num);
    Glice_PT->AddProperty("SPECULARSPIKECONSTANT", ephoton, Glice_SPIKE, num);
    Glice_PT->AddProperty("BACKSCATTERCONSTANT", ephoton, Glice_SCATTER, num);
    G4OpticalSurface *Glice_Surface = new G4OpticalSurface("GliceSurface", unified, groundbackpainted, dielectric_dielectric);
    Glice_Surface->SetMaterialPropertiesTable(Glice_PT);
    Glice_Surface->SetSigmaAlpha(0.);
    for (auto i : glice_phy1)
    {
        new G4LogicalBorderSurface("Glice_surf", i, fAir_phy, Glice_Surface);
        new G4LogicalBorderSurface("Glice_surf", i, this->fAir_phy, Glice_Surface);
    }
    for (auto i : glice_phy2)
    {
        new G4LogicalBorderSurface("Glice_surf", i, fAir_phy, Glice_Surface);
        new G4LogicalBorderSurface("Glice_surf", i, this->fAir_phy, Glice_Surface);
    }

    //LightGuideの反射設定
    for (auto i:guide_phy)
    {
        new G4LogicalBorderSurface("Glice_surf", i, fAir_phy, Glice_Surface);
        new G4LogicalBorderSurface("Glice_surf", i, this->fAir_phy, Glice_Surface);
    }
}

G4LogicalVolume *VolumeDetector::getScintLogical()
{
    return scint->getScintLogical();
}
G4LogicalVolume *VolumeDetector::getMPPCLogical()
{
    return mppc->getMPPCLogical();
}
