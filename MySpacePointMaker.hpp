#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/SpacePointFormation/SpacePointBuilder.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <memory>
#include <string>
#include <vector>

#include "MyFtdDetector.h"

namespace Acts {
class TrackingGeometry;
}

namespace ActsExamples {
struct AlgorithmContext;

class MySpacePointMaker final : public IAlgorithm {
 public:
  struct Config {
    std::string inputMeasurements;
    std::string outputSpacePoints;
    std::shared_ptr<MyFtdDetector> detector{nullptr};
    std::vector<Acts::GeometryIdentifier> geometrySelection;
    int maxDeltaStrawId{6};
    int minMeasPerCand{3};
  };

  struct Candidate {
    std::vector<ActsExamples::IndexSourceLink> sourceLinks;
    int station = -1;
    double chi2 = -1;
    double tx = 0;
    double ty = 0;
    double k = 0;
    int sharedMeasurements = 0;
    double x = 0;
    double y = 0;
    double z = 0;
  };

  MySpacePointMaker(Config cfg, Acts::Logging::Level lvl);
  ProcessCode execute(const AlgorithmContext& ctx) const override;
  const Config& config() const { return m_cfg; }
  double analytic(std::vector<double> &a, std::vector<double> &b, std::vector<double> &g, std::vector<double> &s, 
                  double &t, double &k, double &dt, double &dk, bool debug = 0) const;
  
  double helix(std::vector<double> &z, std::vector<double> &s, std::vector<double> &c, std::vector<double> &g, 
               std::vector<double> &d, double &tx, double &ty, double &p, bool debug = 0) const;

  // double parabolic(const std::vector<double>& z, const std::vector<double>& s,  const std::vector<double>& c,  
  //                  const std::vector<double>& g, const std::vector<double>& d, double &tx, double &ty, double &k) const;

 private:
  Config m_cfg;
  std::optional<IndexSourceLink::SurfaceAccessor> m_slSurfaceAccessor;
  Acts::SpacePointBuilder<SimSpacePoint> m_spacePointBuilder;
  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this, "InputMeasurements"};
  WriteDataHandle<SimSpacePointContainer> m_outputSpacePoints{this, "OutputSpacePoints"};
};
}
