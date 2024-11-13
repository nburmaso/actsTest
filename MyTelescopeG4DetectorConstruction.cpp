#include <utility>

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4VSolid.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

#include "MyTelescope.h"
#include "MyTelescopeG4DetectorConstruction.h"

MyTelescope::MyTelescopeG4DetectorConstruction::MyTelescopeG4DetectorConstruction(
  const MyTelescopeDetector::Config& cfg,
  std::vector<std::shared_ptr<ActsExamples::RegionCreator>> regionCreators)
  : m_cfg(cfg), m_regionCreators(std::move(regionCreators))
{
}

G4VPhysicalVolume* MyTelescope::MyTelescopeG4DetectorConstruction::Construct()
{
  if (m_world != nullptr) {
    return m_world;
  }

  G4double center = (m_cfg.positions.back() + m_cfg.positions.front()) * 0.5 * cm;
  G4double length = (m_cfg.positions.back() - m_cfg.positions.front()) * cm;

  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // World
  G4double worldSize = std::max({std::abs(m_cfg.offsets[0]) + m_cfg.bounds[0] * 0.5,
                                 std::abs(m_cfg.offsets[1]) + m_cfg.bounds[1] * 0.5,
                                 m_cfg.positions.back() + m_cfg.thickness});

  // Envelope parameters
  G4double envSizeX = m_cfg.bounds[0] * cm;
  G4double envSizeY = m_cfg.bounds[1] * cm;
  G4double envSizeZ = length + m_cfg.thickness * cm;

  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = true;

  // Materials
  G4Material* galactic = nist->FindOrBuildMaterial("G4_Galactic");
  G4Material* silicon = nist->FindOrBuildMaterial("G4_Si");

  // Construct the rotation
  // This assumes the binValue is BinningValue::binX, BinningValue::binY or
  // BinningValue::binZ. No reset is necessary in case of BinningValue::binZ
  G4RotationMatrix* rotation = nullptr;
  if (static_cast<Acts::BinningValue>(m_cfg.binValue) == Acts::BinningValue::binX) {
    rotation = new G4RotationMatrix({0, 0, 1}, {0, 1, 0}, {-1, 0, 0});
  } else if (static_cast<Acts::BinningValue>(m_cfg.binValue) == Acts::BinningValue::binY) {
    rotation = new G4RotationMatrix({1, 0, 0}, {0, 0, 1}, {0, -1, 0});
  }

  // World
  auto* solidWorld = new G4Box("World Solid", worldSize * cm, worldSize * cm, worldSize * cm);

  auto* logicWorld = new G4LogicalVolume(solidWorld, galactic, "World Logic");

  m_world = new G4PVPlacement(nullptr,         // no rotation
                              G4ThreeVector(), // position
                              logicWorld,      // its logical volume
                              "World Phys",    // its name
                              nullptr,         // its mother volume
                              false,           // no boolean operation
                              0,               // copy number
                              checkOverlaps);  // overlaps checking

  // Envelope 1
  auto* solidEnv = new G4Box("Envelope Solid",                                               // its name
                             0.5 * envSizeX * cm, 0.5 * envSizeY * cm, 0.5 * envSizeZ * cm); // its size

  auto* logicEnv1 = new G4LogicalVolume(solidEnv,             // its solid
                                        galactic,             // its material
                                        "Envelope #1 Logic"); // its name

  G4VPhysicalVolume* physEnv1 = new G4PVPlacement(rotation,           // rotation
                                                  G4ThreeVector(),    // at detector center
                                                  logicEnv1,          // its logical volume
                                                  "Envelope #1 Phys", // its name
                                                  logicWorld,         // its mother volume
                                                  false,              // no boolean operation
                                                  0,                  // copy number
                                                  checkOverlaps);     // overlaps checking

  // Envelope 2
  auto* logicEnv2 = new G4LogicalVolume(solidEnv,             // its solid
                                        galactic,             // its material
                                        "Envelope #2 Logic"); // its name

  G4VPhysicalVolume* physEnv2 = new G4PVPlacement(nullptr, // no rotation
                                                  G4ThreeVector(m_cfg.offsets[0] * cm, m_cfg.offsets[1] * cm,
                                                                center), // at detector center
                                                  "Envelope #2 Phys",    // its name
                                                  logicEnv2,             // its logical volume
                                                  physEnv1,              // its mother volume
                                                  false,                 // no boolean operation
                                                  0,                     // copy number
                                                  checkOverlaps);        // overlaps checking

  for (std::size_t i = 0; i < m_cfg.positions.size(); ++i) {
    double rIn = m_cfg.rIns[i];
    double rOut = m_cfg.rOuts[i];

    // Layer
    G4VSolid* solidLayer = nullptr;
    if (m_cfg.surfaceType == static_cast<int>(TelescopeSurfaceType::Plane)) {
      solidLayer = new G4Box("Layer Solid", 0.5 * m_cfg.bounds[0] * cm, 0.5 * m_cfg.bounds[1] * cm, 0.5 * m_cfg.thickness * cm);
    } else {
      solidLayer = new G4Tubs("Layer Solid", rIn * cm, rOut * cm, 0.5 * m_cfg.thickness * cm, 0., 2. * M_PI);
    }

    auto* logicLayer = new G4LogicalVolume(solidLayer,     // its solid
                                           silicon,        // its material
                                           "Layer Logic"); // its name

    new G4PVPlacement(
      nullptr,                                               // no rotation
      G4ThreeVector(0, 0, m_cfg.positions[i] * cm - center), // at position
      "Layer #" + std::to_string(i) + " Phys",               // its name
      logicLayer,                                            // its logical volume
      physEnv2,                                              // its mother volume
      false,                                                 // no boolean operation
      0,                                                     // copy number
      checkOverlaps);                                        // overlaps checking
  }

  // Create regions
  for (const auto& regionCreator : m_regionCreators) {
    regionCreator->Construct();
  }

  return m_world;
}

MyTelescope::MyTelescopeG4DetectorConstructionFactory::
  MyTelescopeG4DetectorConstructionFactory(const MyTelescopeDetector::Config& cfg,
                                           std::vector<std::shared_ptr<ActsExamples::RegionCreator>> regionCreators)
  : m_cfg(cfg), m_regionCreators(std::move(regionCreators)) {}

std::unique_ptr<G4VUserDetectorConstruction> MyTelescope::MyTelescopeG4DetectorConstructionFactory::factorize() const
{
  return std::make_unique<MyTelescopeG4DetectorConstruction>(m_cfg, m_regionCreators);
}
