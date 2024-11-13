#include "MyTelescope.h"

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/BinningType.hpp"

#include "BuildMyTelescopeDetector.h"
#include "ActsExamples/TelescopeDetector/TelescopeDetectorElement.hpp"

#include <algorithm>
#include <stdexcept>

auto MyTelescope::MyTelescopeDetector::finalize(const Config& cfg) -> std::pair<TrackingGeometryPtr, ContextDecorators>
{
  MyDetectorElement::ContextType nominalContext;

  if (cfg.surfaceType > 1) {
    throw std::invalid_argument(
      "The surface type could either be 0 for plane surface or 1 for disc "
      "surface.");
  }
  if (cfg.binValue > 2) {
    throw std::invalid_argument("The axis value could only be 0, 1, or 2.");
  }
  // Check if the bounds values are valid
  if (cfg.surfaceType == 1 && cfg.bounds[0] >= cfg.bounds[1]) {
    throw std::invalid_argument(
      "The minR should be smaller than the maxR for disc surface bounds.");
  }

  if (cfg.positions.size() != cfg.stereos.size()) {
    throw std::invalid_argument(
      "The number of provided positions must match the number of "
      "provided stereo angles.");
  }

  config = cfg;

  // Sort the provided distances
  std::vector<double> positions = cfg.positions;
  std::vector<double> stereos = cfg.stereos;
  std::vector<double> rIns = cfg.rIns;
  std::vector<double> rOuts = cfg.rOuts;

  /// Return the telescope detector
  TrackingGeometryPtr gGeometry = MyTelescope::buildDetector(
    nominalContext, detectorStore, positions, rIns, rOuts,
    stereos, cfg.offsets,
    cfg.bounds, cfg.thickness,
    static_cast<MyTelescope::TelescopeSurfaceType>(cfg.surfaceType),
    static_cast<Acts::BinningValue>(cfg.binValue));

  ContextDecorators gContextDecorators = {};

  // return the pair of geometry and empty decorators
  return std::make_pair<TrackingGeometryPtr, ContextDecorators>(
    std::move(gGeometry), std::move(gContextDecorators));
}