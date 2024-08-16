// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"

using Stepper = Acts::EigenStepper<>;
using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;
using Fitter = Acts::KalmanFitter<Propagator, Acts::VectorMultiTrajectory>;


#include <memory>
#include <string>

namespace ActsExamples {

class MyCalibrator {
 public:
  MyCalibrator(){}
  void calibrate(const Acts::GeometryContext& gctx,
                 const Acts::CalibrationContext& cctx,
                 const Acts::SourceLink& sourceLink,
                 Acts::VectorMultiTrajectory::TrackStateProxy trackState) const {}
};

struct AlgorithmContext;

class MyRefittingAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    std::string inputTracks;
    std::string outputTracks;
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField;
  };

  // struct SimpleReverseFilteringLogic {
  //   double momentumThreshold = 0;
  //   bool doBackwardFiltering(Acts::VectorMultiTrajectory::ConstTrackStateProxy trackState) const {
  //     auto momentum = fabs(1 / trackState.filtered()[Acts::eBoundQOverP]);
  //     return (momentum <= momentumThreshold);
  //   }
  // };

  MyRefittingAlgorithm(Config config, Acts::Logging::Level level):ActsExamples::IAlgorithm("TrackFittingAlgorithm", level), m_cfg(std::move(config)) {
    m_inputTracks.initialize(m_cfg.inputTracks);
    m_outputTracks.initialize(m_cfg.outputTracks);
  }

  ActsExamples::ProcessCode execute(const AlgorithmContext& ctx) const final {
    const Stepper stepper(std::move(m_cfg.magneticField));
    Acts::Navigator::Config cfg{std::move(m_cfg.trackingGeometry)};
    Acts::Navigator navigator(cfg, logger().cloneWithSuffix("Navigator"));
    Propagator propagator(stepper, std::move(navigator),logger().cloneWithSuffix("Propagator"));
    Fitter fitter(std::move(propagator), logger().cloneWithSuffix("Fitter"));
    IndexSourceLink::SurfaceAccessor slSurfaceAccessor(*m_cfg.trackingGeometry);
    MyCalibrator calibrator;

    Acts::KalmanFitterExtensions<Acts::VectorMultiTrajectory> extensions;
    extensions.surfaceAccessor.connect<&IndexSourceLink::SurfaceAccessor::operator()>(&slSurfaceAccessor);
    extensions.calibrator.connect<&MyCalibrator::calibrate>(&calibrator);

    // Acts::GainMatrixUpdater kfUpdater;
    // Acts::GainMatrixSmoother kfSmoother;
    // SimpleReverseFilteringLogic reverseFilteringLogic;
    //extensions.updater.connect<&Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(&kfUpdater);
    // extensions.smoother.connect<&Acts::GainMatrixSmoother::operator()<Acts::VectorMultiTrajectory>>(&kfSmoother);
    //extensions.reverseFilteringLogic.connect<&SimpleReverseFilteringLogic::doBackwardFiltering>(&reverseFilteringLogic);

    Acts::KalmanFitterOptions<Acts::VectorMultiTrajectory> kfOptions(ctx.geoContext, ctx.magFieldContext, ctx.calibContext, extensions, Acts::PropagatorPlainOptions());
    
    auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
    auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
    TrackContainer tracks(trackContainer, trackStateContainer);

    std::vector<Acts::SourceLink> sourceLinks;
    std::vector<const Acts::Surface*> surfSequence;
    for (const auto& track : m_inputTracks(ctx)) {
      const Acts::BoundTrackParameters initialParams(track.referenceSurface().getSharedPtr(), track.parameters(), track.covariance(), track.particleHypothesis());
      sourceLinks.clear();
      surfSequence.clear();
      for (auto state : track.trackStatesReversed()) {
        surfSequence.push_back(&state.referenceSurface());
        sourceLinks.push_back(state.getUncalibratedSourceLink());
      }
      kfOptions.referenceSurface = &track.referenceSurface();
      fitter.fit(sourceLinks.begin(), sourceLinks.end(), initialParams, kfOptions, tracks);
    }

    ConstTrackContainer constTracks{std::make_shared<Acts::ConstVectorTrackContainer>(std::move(*trackContainer)),
                                    std::make_shared<Acts::ConstVectorMultiTrajectory>(std::move(*trackStateContainer))};

    m_outputTracks(ctx, std::move(constTracks));
    return ActsExamples::ProcessCode::SUCCESS;
  }

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  ReadDataHandle<ConstTrackContainer> m_inputTracks{this, "InputTracks"};
  WriteDataHandle<ConstTrackContainer> m_outputTracks{this, "OutputTracks"};
};

}  // namespace ActsExamples
