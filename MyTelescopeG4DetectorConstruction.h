#pragma once

#include "G4VUserDetectorConstruction.hh"
#include "ActsExamples/Geant4/RegionCreator.hpp"
#include "ActsExamples/Geant4/DetectorConstructionFactory.hpp"

#include "MyTelescope.h"

namespace MyTelescope
{

class MyTelescopeG4DetectorConstruction final
  : public G4VUserDetectorConstruction
{
 public:
  MyTelescopeG4DetectorConstruction(
    const MyTelescopeDetector::Config& cfg,
    std::vector<std::shared_ptr<ActsExamples::RegionCreator>> regionCreators = {});

  G4VPhysicalVolume* Construct() final;

 private:
  /// The configuration of the telescope detector
  MyTelescopeDetector::Config m_cfg;
  /// Region creators
  std::vector<std::shared_ptr<ActsExamples::RegionCreator>> m_regionCreators;
  /// The world volume
  G4VPhysicalVolume* m_world{};
};

class MyTelescopeG4DetectorConstructionFactory final
  : public ActsExamples::DetectorConstructionFactory
{
 public:
  explicit MyTelescopeG4DetectorConstructionFactory(const MyTelescopeDetector::Config& cfg,
                                                    std::vector<std::shared_ptr<ActsExamples::RegionCreator>> regionCreators = {});

  std::unique_ptr<G4VUserDetectorConstruction> factorize() const override;

 private:
  /// The configuration of the telescope detector
  MyTelescopeDetector::Config m_cfg;
  /// Region creators
  std::vector<std::shared_ptr<ActsExamples::RegionCreator>> m_regionCreators;
};

} // namespace MyTelescope