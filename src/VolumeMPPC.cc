#include "VolumeMPPC.hh"
#include "G4VisAttributes.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4SystemOfUnits.hh"

VolumeMPPC::VolumeMPPC()
{
    mppc_height = 3.0 * mm;
    mppc_width = 0.05 * mm;

    win_longside = 5.0 * mm;
    win_shortside = 5.0 * mm;
    win_width = 0.1 * mm;

    enve_longside = 6.0 * mm;
    enve_shortside = 6.0 * mm;
    enve_width = 2.1 * mm;

    //形の定義
    fMPPC_box = new G4Box("mppc_box", mppc_height / 2., mppc_height / 2., mppc_width / 2.);
    fWin_box = new G4Box("Win_box", win_longside / 2., win_shortside / 2., win_width / 2.);
    fEnve_box = new G4Box("Enve_box", enve_longside / 2., enve_shortside / 2., enve_width / 2.);

    //物質の定義
    fMPPC_log = new G4LogicalVolume(fMPPC_box, G4Material::GetMaterial("Al"), "mppc_log", 0, 0, 0);
    fWin_log = new G4LogicalVolume(fWin_box, G4Material::GetMaterial("Silicone"), "win_log", 0, 0, 0);
    fEnve_log = new G4LogicalVolume(fEnve_box, G4Material::GetMaterial("ABS"), "enve_log", 0, 0, 0);

    //MPPC全体の定義
    G4ThreeVector MPPC_posi(0, 0, 1.025 * mm);
    G4ThreeVector win_posi(0, 0, 1.0 * mm);

    fWin_phy = new G4PVPlacement(0, win_posi, fWin_log, "Window", fEnve_log, false, 0);
    fMPPC_phy = new G4PVPlacement(0, MPPC_posi, fMPPC_log, "MPPC", fEnve_log, false, 0);
    fEnve_phy = new G4PVPlacement(0, G4ThreeVector(), fEnve_log, "Enve", 0, false, 0);
}

G4LogicalVolume *VolumeMPPC::getLogicalVolume()
{
    return fEnve_log;
}

G4LogicalVolume *VolumeMPPC::getMPPCLogical()
{
    return fMPPC_log;
}

void VolumeMPPC::VisAttributes()
{
    mppc_va = new G4VisAttributes(G4Colour(0.2, 0.4, 0.2));
    fMPPC_log->SetVisAttributes(mppc_va);
    win_va = new G4VisAttributes(G4Colour(0.4, 0.4, 0.6));
    fWin_log->SetVisAttributes(win_va);
    enve_va = new G4VisAttributes(G4Colour(0.8, 0.8, 0.8));
    fEnve_log->SetVisAttributes(enve_va);
}

void VolumeMPPC::SurfaceProperties()
{
    G4double ephoton[] = {2.76 * eV, 3.06 * eV};
    const G4int num = sizeof(ephoton) / sizeof(G4double);

    //Windowの反射設定 (Envelope周り)
    G4double Window_RIND[] = {1.0, 1.0};
    G4double Window_REF[] = {1.0, 1.0};
    G4MaterialPropertiesTable *Window_PT = new G4MaterialPropertiesTable();
    Window_PT->AddProperty("RINDEX", ephoton, Window_RIND, num);
    Window_PT->AddProperty("REFLECTIVITY", ephoton, Window_REF, num);
    G4OpticalSurface *Window_Surface = new G4OpticalSurface("WindowSurface", glisur, groundbackpainted, dielectric_dielectric);
    Window_Surface->SetMaterialPropertiesTable(Window_PT);

    G4LogicalBorderSurface *surface1 = new G4LogicalBorderSurface("Window_surf", fWin_phy, fEnve_phy, Window_Surface);

    //mppcの反射設定(光子感度)
    G4double MPPC_ReR[] = {1.92, 1.92};
    G4double MPPC_ImR[] = {1.69, 1.69};
    G4double MPPC_EFF[] = {1.0, 1.0}; //enables 'detection' of photons
    G4MaterialPropertiesTable *MPPC_mt = new G4MaterialPropertiesTable();
    MPPC_mt->AddProperty("REALRINDEX", ephoton, MPPC_ReR, num);
    MPPC_mt->AddProperty("IMAGINARYRINDEX", ephoton, MPPC_ImR, num);
    MPPC_mt->AddProperty("EFFICIENCY", ephoton, MPPC_EFF, num);
    G4OpticalSurface *MPPC_Surface = new G4OpticalSurface("MPPCSurface", glisur, polished, dielectric_metal);
    MPPC_Surface->SetMaterialPropertiesTable(MPPC_mt);

    //Create logical skin surfaces
    G4LogicalSkinSurface *mppc_surf = new G4LogicalSkinSurface("mppc_surf", fMPPC_log, MPPC_Surface);
}
