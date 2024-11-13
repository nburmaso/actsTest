#include "BuildMyTelescopeDetector.h"

#include <memory>

#include "Acts/Geometry/NavigationLayer.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/DiscLayer.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/ILayerArrayCreator.hpp"
#include "Acts/Geometry/ITrackingVolumeHelper.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"

#include "Acts/Plugins/Geant4/Geant4Converters.hpp"

#include "Geant4/G4Material.hh"
#include "Geant4/G4NistManager.hh"

using Acts::UnitConstants::cm;

std::unique_ptr<const Acts::TrackingGeometry>
  MyTelescope::buildDetector(
    const typename MyTelescope::MyDetectorElement::ContextType& gctx,
    std::vector<std::shared_ptr<MyTelescope::MyDetectorElement>>& detectorStore,
    const std::vector<double>& positions,
    const std::vector<double>& rIns,
    const std::vector<double>& rOuts,
    const std::vector<double>& stereoAngles,
    const std::array<double, 2>& offsets, const std::array<double, 2>& bounds,
    double thickness, MyTelescope::TelescopeSurfaceType surfaceType,
    Acts::BinningValue binValue)
{
  using namespace Acts::UnitLiterals;

  // Create materials
  auto* nist = G4NistManager::Instance();
  G4Material* siMat = nist->FindOrBuildMaterial("G4_Si");
  G4Material* alMat = nist->FindOrBuildMaterial("G4_Al");
  // G4Material* pbMat = nist->FindOrBuildMaterial("G4_Pb");
  G4Material* worldMat = nist->FindOrBuildMaterial("G4_Galactic");

  Acts::Geant4MaterialConverter converter;
  Acts::Material silicon = converter.material(*siMat);
  Acts::Material aluminium = converter.material(*alMat);
  Acts::Material vacuum = converter.material(*worldMat);

  Acts::MaterialSlab matProp(silicon, thickness * cm);
  const auto surfaceMaterial = std::make_shared<Acts::HomogeneousSurfaceMaterial>(matProp);
  const auto volumeMaterial = std::make_shared<Acts::HomogeneousVolumeMaterial>(vacuum);

  // Construct the rotation
  // This assumes the binValue is binX, binY or binZ. No reset is necessary in
  // case of binZ
  Acts::RotationMatrix3 rotation = Acts::RotationMatrix3::Identity();
  if (binValue == Acts::BinningValue::binX) {
    rotation.col(0) = Acts::Vector3(0, 0, -1);
    rotation.col(1) = Acts::Vector3(0, 1, 0);
    rotation.col(2) = Acts::Vector3(1, 0, 0);
  } else if (binValue == Acts::BinningValue::binY) {
    rotation.col(0) = Acts::Vector3(1, 0, 0);
    rotation.col(1) = Acts::Vector3(0, 0, -1);
    rotation.col(2) = Acts::Vector3(0, 1, 0);
  }

  // Construct the surfaces and layers
  std::size_t nLayers = positions.size();
  std::vector<Acts::LayerPtr> layers(nLayers);
  for (unsigned int i = 0; i < nLayers; ++i) {
    double rIn = rIns[i];
    double rOut = rOuts[i];
    const auto rBounds = std::make_shared<const Acts::RadialBounds>(rIn * cm, rOut * cm);     // <- for disk-like layers
    const auto pBounds = std::make_shared<const Acts::RectangleBounds>(rOut * cm, rOut * cm); // <- for square-like layers
    // The translation without rotation yet
    Acts::Translation3 trans(offsets[0], offsets[1], positions[i] * cm);
    // The entire transformation (the coordinate system, whose center is defined
    // by trans, will be rotated as well)
    Acts::Transform3 trafo(rotation * trans);

    // rotate around local z axis by stereo angle
    auto stereo = stereoAngles[i];
    trafo *= Acts::AngleAxis3(stereo, Acts::Vector3::UnitZ());

    // Create the detector element
    std::shared_ptr<MyDetectorElement> detElement = nullptr;
    std::shared_ptr<Acts::Surface> surface = nullptr;

    if (surfaceType == TelescopeSurfaceType::Plane) {
      surface = Acts::Surface::makeShared<Acts::PlaneSurface>(trafo, pBounds);
      auto surArray = std::make_unique<Acts::SurfaceArray>(surface);
      layers[i] = Acts::PlaneLayer::create(trafo, pBounds, std::move(surArray), 1._mm); // <- for square-like layers
    } else {
      surface = Acts::Surface::makeShared<Acts::DiscSurface>(trafo, rBounds);
      auto surArray = std::make_unique<Acts::SurfaceArray>(surface);
      layers[i] = Acts::DiscLayer::create(trafo, rBounds, std::move(surArray), 1._mm); // <- for disc-like layers
    }

    // Associate the layer to the surface
    detElement = std::make_shared<MyDetectorElement>(std::make_shared<const Acts::Transform3>(trafo), surface, thickness * cm, surfaceMaterial);

    surface->assignSurfaceMaterial(surfaceMaterial);
    surface->associateLayer(*layers[i]);

    surface->assignDetectorElement(*detElement);
    // detElement->surface() = *surface;
    detectorStore.push_back(std::move(detElement));
  }

  // The volume transform
  Acts::Translation3 transVol(offsets[0] * cm, offsets[1] * cm, (positions.front() + positions.back()) * 0.5 * cm);
  Acts::Transform3 trafoVol(rotation * transVol);

  // The volume bounds is set to be a bit larger than either cubic with planes
  // or cylinder with discs
  auto length = positions.back() - positions.front();
  std::shared_ptr<Acts::VolumeBounds> boundsVol = nullptr;
  if (surfaceType == TelescopeSurfaceType::Plane) {
    boundsVol = std::make_shared<Acts::CuboidVolumeBounds>(
      bounds[0] * cm + 5._mm, bounds[1] * cm + 5._mm, length * cm + 10._mm);
  } else {
    boundsVol = std::make_shared<Acts::CylinderVolumeBounds>(
      std::max(bounds[0] * cm - 5.0_mm, 0.), bounds[1] * cm + 5._mm, length * cm + 10._mm);
  }

  // Create layer array
  Acts::LayerVector layVec;
  for (unsigned int i = 0; i < nLayers; i++) {
    layVec.push_back(layers[i]);
  }
  Acts::GeometryContext genGctx{gctx};
  Acts::LayerArrayCreator::Config lacConfig;
  Acts::LayerArrayCreator layArrCreator(lacConfig, Acts::getDefaultLogger("LayerArrayCreator", Acts::Logging::INFO));
  std::unique_ptr<const Acts::LayerArray> layArr(layArrCreator.layerArray(genGctx, layVec, positions.front() * cm - 2._mm, positions.back() * cm + 2._mm, Acts::BinningType::arbitrary, binValue));

  // Build mother tracking volume
  auto trackVolume = std::make_shared<Acts::TrackingVolume>(trafoVol, boundsVol, volumeMaterial, std::move(layArr), nullptr, Acts::MutableTrackingVolumeVector{}, "Telescope");

  // Build tracking geometry
  auto trackingGeometry = std::make_unique<Acts::TrackingGeometry>(trackVolume);

  // debug
  const Acts::TrackingVolume* highestTrackingVolume = trackingGeometry->highestTrackingVolume();
  printf("volumeId = %lu\n", highestTrackingVolume->geometryId().value());
  const Acts::LayerArray* confinedLayers = highestTrackingVolume->confinedLayers();
  for (const auto& layer : confinedLayers->arrayObjects()) {
    printf("layerId = %lu, thickness = %f type = %d\n", layer->geometryId().value(), layer->thickness(), layer->layerType());
    if (layer->layerType() == -1) {
      const Acts::NavigationLayer* nlayer = dynamic_cast<const Acts::NavigationLayer*>(layer.get());
      printf(" navigation %p\n", nlayer);
    } else if (layer->layerType() == 0) {
      const Acts::DiscLayer* player = dynamic_cast<const Acts::DiscLayer*>(layer.get());
      double rmin = player->surfaceArray()->surfaces()[0]->bounds().values()[0];
      double rmax = player->surfaceArray()->surfaces()[0]->bounds().values()[1];
      double zpos = player->surfaceArray()->surfaces()[0]->center(genGctx).z();
      printf(" disc %p, rIn:%f, rOut:%f, zPos:%f\n", player, rmin, rmax, zpos);
    } else {
      const Acts::PlaneLayer* player = dynamic_cast<const Acts::PlaneLayer*>(layer.get());
      printf(" plane %p\n", player);
    }
    if (!layer->surfaceArray())
      continue;
    for (const auto& surface : layer->surfaceArray()->surfaces()) {
      printf(" surfaceId = %lu\n", surface->geometryId().value());
    }
  }

  return trackingGeometry;
}