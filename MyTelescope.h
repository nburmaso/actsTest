#pragma once

#include <array>
#include <memory>
#include <utility>
#include <vector>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Plugins/Geant4/Geant4Converters.hpp"

using namespace Acts::UnitLiterals;

namespace Acts
{
class TrackingGeometry;
class IMaterialDecorator;
class Surface;
class PlanarBounds;
class DiscBounds;
class ISurfaceMaterial;
} // namespace Acts

namespace ActsExamples
{
class IContextDecorator;
} // namespace ActsExamples

namespace ActsExamples::Telescope
{
class TelescopeDetectorElement;
class TelescopeG4DetectorConstruction;
} // namespace ActsExamples::Telescope

namespace MyTelescope
{

/// The telescope detector surface type
enum class TelescopeSurfaceType {
  Plane = 0,
  Disc = 1,
};

class MyDetectorElement : public Acts::DetectorElementBase
{
 public:
  MyDetectorElement(std::shared_ptr<const Acts::Transform3> transform,
                    std::shared_ptr<Acts::Surface> surface,
                    double thickness,
                    std::shared_ptr<Acts::HomogeneousSurfaceMaterial> material)
    : Acts::DetectorElementBase(), mTransform(std::move(transform)), mSurface(std::move(surface)), mThickness(thickness), mMaterial(std::move(material))
  {
  }

  struct ContextType {
    /// The current interval of validity
    unsigned int iov = 0;
  };

  const Acts::Transform3& transform(const Acts::GeometryContext& gctx) const override { return *mTransform; }

  const Acts::Surface& surface() const override { return *mSurface; }

  Acts::Surface& surface() override { return *mSurface; }

  double thickness() const override { return mThickness; }

 private:
  std::shared_ptr<const Acts::Transform3> mTransform = nullptr;
  std::shared_ptr<Acts::Surface> mSurface = nullptr;
  std::shared_ptr<Acts::HomogeneousSurfaceMaterial> mMaterial = nullptr;
  double mThickness = 0.;
};

struct MyTelescopeDetector {
  using DetectorElementPtr = std::shared_ptr<MyDetectorElement>;
  using DetectorStore = std::vector<DetectorElementPtr>;

  using ContextDecorators =
    std::vector<std::shared_ptr<ActsExamples::IContextDecorator>>;
  using TrackingGeometryPtr = std::shared_ptr<const Acts::TrackingGeometry>;

  struct Config {
    std::vector<double> positions{{210, 232.5, 255, 277.5, 300}};
    std::vector<double> rIns{{35.7, 35.7, 35.7, 35.7, 35.7}};
    std::vector<double> rOuts{{130., 130., 130., 130., 130.}};
    std::vector<double> stereos{{0., 0., 0., 0., 0.}};
    std::array<double, 2> offsets{{0, 0}};
    std::array<double, 2> bounds{{35.7, 130.}};
    double thickness{0.01};
    int surfaceType{0}; // 0 for plane surface or 1 for disc
    int binValue{2};
  };

  Config config;
  /// The store of the detector elements (lifetime: job)
  DetectorStore detectorStore;

  std::pair<TrackingGeometryPtr, ContextDecorators> finalize(const Config& cfg);
};

} // namespace MyTelescope
