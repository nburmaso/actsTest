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
    std::string inputMeasurementParticlesMap;
    std::string outputSpacePoints;
    std::shared_ptr<MyFtdDetector> detector{nullptr};
    std::vector<Acts::GeometryIdentifier> geometrySelection;
    int maxLoDeltaPStrawId1{2};
    int maxLoDeltaPStrawId2{2};
    int maxLoDeltaPStrawId3{2};
    int maxHiDeltaPStrawId1{6};
    int maxHiDeltaPStrawId2{6};
    int maxHiDeltaPStrawId3{12};
    int maxLoDeltaMStrawId1{2};
    int maxLoDeltaMStrawId2{2};
    int maxLoDeltaMStrawId3{2};
    int maxHiDeltaMStrawId1{6};
    int maxHiDeltaMStrawId2{6};
    int maxHiDeltaMStrawId3{12};
    int minMeasPerCand{3};
    double maxChi2{10.};
    int maximumIterations{10000};
    int maximumSharedStraws{1};
  };

  struct PreCandidate {
    std::vector<ActsExamples::IndexSourceLink> sourceLinks;
    int refTypeA = -1;
    int refTypeB = -1;
    int refIdxA  = -1;
    int refIdxB  = -1;
  };

  struct Candidate {
    std::vector<ActsExamples::IndexSourceLink> sourceLinks;
    int station = -1;
    std::vector<int> straws;
    double chi2 = -1;
    double chi2ndf = -1;
    double tx = 0;
    double ty = 0;
    double k = 0;
    double varxx = 0;
    double varyy = 0;
    double varxy = 0;
    int sharedStraws = 0;
    double x = 0;
    double y = 0;
    double z = 0;
  };

  MySpacePointMaker(Config cfg, Acts::Logging::Level lvl);
  ProcessCode execute(const AlgorithmContext& ctx) const override;
  const Config& config() const { return m_cfg; }

 private:

  std::tuple<double,double,double,double,double,double> linear(PreCandidate& cand, const std::vector<std::array<double, 5>>& cacheZSCGD, bool debug = 0) const;

  std::tuple<double,double,double,double> parabolic(PreCandidate& cand, const std::vector<std::array<double, 5>>& cacheZSCGD, bool debug = 0) const;

  Config m_cfg;
  std::optional<IndexSourceLink::SurfaceAccessor> m_slSurfaceAccessor;
  Acts::SpacePointBuilder<SimSpacePoint> m_spacePointBuilder;
  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this, "InputMeasurements"};
  ReadDataHandle<MeasurementParticlesMap> m_inputMeasurementParticlesMap{this, "InputMeasurementParticlesMap"};
  WriteDataHandle<SimSpacePointContainer> m_outputSpacePoints{this, "OutputSpacePoints"};
};
}
