#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
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

#include "Acts/TrackFitting/GlobalChiSquareFitter.hpp"

#include "ActsExamples/TrackFitting/RefittingCalibrator.hpp"


using Stepper = Acts::EigenStepper<>;
using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;
using Fitter = Acts::KalmanFitter<Propagator, Acts::VectorMultiTrajectory>;
using Gx2Fitter = Acts::Experimental::Gx2Fitter<Propagator, Acts::VectorMultiTrajectory>;

#include <memory>
#include <string>

namespace ActsExamples {

class MyCalibrator {
 public:
  MyCalibrator() = default;
  void calibrate( const Acts::GeometryContext& gctx,
                  const Acts::CalibrationContext& cctx,
                  const Acts::SourceLink& sourceLink,
                  Acts::VectorMultiTrajectory::TrackStateProxy trackState) const {
    const auto &idxSourceLink = sourceLink.get<IndexSourceLink>();
    const ConstVariableBoundMeasurementProxy measurement = measurements.getMeasurement(idxSourceLink.index());
    Acts::visit_measurement(measurement.size(), [&](auto N) -> void {
        constexpr std::size_t kMeasurementSize = decltype(N)::value;
        const auto fixedMeasurement = static_cast<ConstFixedBoundMeasurementProxy<kMeasurementSize>>(measurement);
        Acts::ActsVector<kMeasurementSize> calibratedParameters = fixedMeasurement.parameters();
        Acts::ActsSquareMatrix<kMeasurementSize> calibratedCovariance = fixedMeasurement.covariance();
        trackState.allocateCalibrated(kMeasurementSize);
        trackState.calibrated<kMeasurementSize>() = calibratedParameters;
        trackState.calibratedCovariance<kMeasurementSize>() = calibratedCovariance;
        trackState.setSubspaceIndices(fixedMeasurement.subspaceIndices());
    });
  }

  ActsExamples::MeasurementContainer measurements;
};

class MyRefittingAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    std::string inputTracks;
    std::string outputTracks;
    std::string inputMeasurements;
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField;
  };

  MyRefittingAlgorithm(Config config, Acts::Logging::Level level):ActsExamples::IAlgorithm("TrackFittingAlgorithm", level), m_cfg(std::move(config)) {
    m_inputMeasurements.initialize(m_cfg.inputMeasurements);
    m_inputTracks.initialize(m_cfg.inputTracks);
    m_outputTracks.initialize(m_cfg.outputTracks);
  }

  ActsExamples::ProcessCode execute(const AlgorithmContext& ctx) const final {
    const Stepper stepper(std::move(m_cfg.magneticField));
    Acts::Navigator::Config cfg{std::move(m_cfg.trackingGeometry)};
    Acts::Navigator navigator(cfg, logger().cloneWithSuffix("Navigator"));
    Propagator propagator(stepper, std::move(navigator),logger().cloneWithSuffix("Propagator"));

    IndexSourceLink::SurfaceAccessor slSurfaceAccessor(*m_cfg.trackingGeometry);
    MyCalibrator calibrator;
    calibrator.measurements = m_inputMeasurements(ctx);

    // Fitter fitter(std::move(propagator), logger().cloneWithSuffix("Fitter"));
    // Acts::KalmanFitterExtensions<Acts::VectorMultiTrajectory> extensions;
    // extensions.surfaceAccessor.connect<&IndexSourceLink::SurfaceAccessor::operator()>(&slSurfaceAccessor);
    // extensions.calibrator.connect<&MyCalibrator::calibrate>(&calibrator);
    // Acts::KalmanFitterOptions<Acts::VectorMultiTrajectory> kfOptions(ctx.geoContext, ctx.magFieldContext, ctx.calibContext, extensions, Acts::PropagatorPlainOptions());

    // gx2
    Gx2Fitter fitter2(std::move(propagator), logger().cloneWithSuffix("Gx2Fitter"));
    Acts::Experimental::Gx2FitterExtensions<Acts::VectorMultiTrajectory> gx2extensions;
    gx2extensions.surfaceAccessor.connect<&IndexSourceLink::SurfaceAccessor::operator()>(&slSurfaceAccessor);
    gx2extensions.calibrator.connect<&MyCalibrator::calibrate>(&calibrator);
    Acts::Experimental::Gx2FitterOptions<Acts::VectorMultiTrajectory> gx2fOptions(ctx.geoContext, ctx.magFieldContext, ctx.calibContext, gx2extensions, Acts::PropagatorPlainOptions(ctx.geoContext, ctx.magFieldContext));

    auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
    auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
    TrackContainer tracks(trackContainer, trackStateContainer);

    std::vector<Acts::SourceLink> sourceLinks;
    std::vector<const Acts::Surface*> surfSequence;
    for (const auto& track : m_inputTracks(ctx)) {
      const Acts::BoundTrackParameters initialParams(track.referenceSurface().getSharedPtr(), track.parameters(), track.covariance(), Acts::ParticleHypothesis::pion());
      //const Acts::BoundTrackParameters initialParams(track.referenceSurface().getSharedPtr(), track.parameters(), track.covariance(), track.particleHypothesis());
      sourceLinks.clear();
      surfSequence.clear();
      for (auto state : track.trackStatesReversed()) {
        surfSequence.push_back(&state.referenceSurface());
        // auto sl = RefittingCalibrator::RefittingSourceLink{state};
        // sourceLinks.push_back(Acts::SourceLink{sl});
        sourceLinks.push_back(state.getUncalibratedSourceLink());
      }

      //ACTS_VERBOSE("MyRefit  " << track.parameters().transpose());
      //kfOptions.referenceSurface = &track.referenceSurface();
      gx2fOptions.referenceSurface = &track.referenceSurface();

      //auto result = fitter.fit(sourceLinks.begin(), sourceLinks.end(), initialParams, kfOptions, tracks);
      auto result2 = fitter2.fit(sourceLinks.begin(), sourceLinks.end(), initialParams, gx2fOptions, tracks);

      // const auto& refittedTrack = result.value();
      // bool hasFittedParams = refittedTrack.hasReferenceSurface();
      // printf("%d\n",hasFittedParams);

      //ACTS_VERBOSE("MyRefit  " << refittedTrack.parameters().transpose());
    }

    ConstTrackContainer constTracks{std::make_shared<Acts::ConstVectorTrackContainer>(std::move(*trackContainer)),
                                    std::make_shared<Acts::ConstVectorMultiTrajectory>(std::move(*trackStateContainer))};

    m_outputTracks(ctx, std::move(constTracks));
    return ActsExamples::ProcessCode::SUCCESS;
  }

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this, "InputMeasurements"};
  ReadDataHandle<ConstTrackContainer> m_inputTracks{this, "InputTracks"};
  WriteDataHandle<ConstTrackContainer> m_outputTracks{this, "OutputTracks"};
};

}  // namespace ActsExamples
