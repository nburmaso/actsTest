#include "MyTrackFindingAlgorithm.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/MeasurementCalibration.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/SympyStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "ActsFatras/EventData/Barcode.hpp"

namespace My{

MyTrackFindingAlgorithm::MyTrackFindingAlgorithm(Config config, Acts::Logging::Level level)
    : IAlgorithm("MyTrackFindingAlgorithm", level), m_cfg(std::move(config)) {
  ACTS_INFO("Constructor for my track finder");
  m_inputMeasurements.initialize(m_cfg.inputMeasurements);
  m_inputMeasurementParticlesMap.initialize(m_cfg.inputMeasurementParticlesMap);
  m_inputInitialTrackParameters.initialize(m_cfg.inputInitialTrackParameters);
  m_outputTracks.initialize(m_cfg.outputTracks);
}

ProcessCode MyTrackFindingAlgorithm::execute(const AlgorithmContext& ctx) const {
  ACTS_INFO("Starting my track finder for event " << ctx.eventNumber);
  const auto& measurements = m_inputMeasurements(ctx);
  const auto& initialParameters = m_inputInitialTrackParameters(ctx);
  const auto& measurementParticlesMap = m_inputMeasurementParticlesMap(ctx);

  auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3{0., 0., 0.});

  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
  TrackContainer tracks(trackContainer, trackStateContainer);

  auto trackContainerTemp = std::make_shared<Acts::VectorTrackContainer>();
  auto trackStateContainerTemp = std::make_shared<Acts::VectorMultiTrajectory>();
  TrackContainer tracksTemp(trackContainerTemp, trackStateContainerTemp);

  // TODO move to constructor
  using Propagator = Acts::Propagator<Acts::EigenStepper<>, Acts::Navigator>;
  using PropagatorOptions = Propagator::template Options<Acts::ActorList<MyTrackFindingActor>>;
  Propagator m_propagator(Acts::EigenStepper<>(m_cfg.magneticField),
                          Acts::Navigator({m_cfg.trackingGeometry}, logger().cloneWithSuffix("Navigator")),
                          logger().clone("MyPropagator",Acts::Logging::INFO));
  PropagatorOptions propOptions(ctx.geoContext, ctx.magFieldContext);
  
  // Extrapolator used to extrapolate tracks to reference surface (sympy stepper instead of eigen stepper)
  using Extrapolator = Acts::Propagator<Acts::SympyStepper, Acts::Navigator>;
  using ExtrapolatorOptions = Extrapolator::template Options<Acts::ActorList<Acts::MaterialInteractor, Acts::EndOfWorldReached>>;
  Extrapolator extrapolator(Acts::SympyStepper(m_cfg.magneticField), 
                            Acts::Navigator({m_cfg.trackingGeometry},logger().cloneWithSuffix("Navigator")),
                            logger().cloneWithSuffix("Propagator"));
  ExtrapolatorOptions extrapolationOptions(ctx.geoContext, ctx.magFieldContext);
  auto extrapolationStrategy = Acts::TrackExtrapolationStrategy::last; // options: first, last, firstOrLast

  ActsExamples::PassThroughCalibrator pcalibrator;
  ActsExamples::MeasurementCalibratorAdapter calibrator(pcalibrator, measurements);
  Acts::GainMatrixUpdater kfUpdater;

  // Define actor properties
  auto& myActor = propOptions.actorList.template get<MyTrackFindingActor>();
  myActor.slAccessor.container = &measurements.orderedIndices();
  myActor.measurements = &measurements;
  myActor.trackStates = &tracksTemp.trackStateContainer();
  // myActor.measurements = &measurements;
  myActor.calibrationContext = &ctx.calibContext;
  myActor.calibrator.template connect<&ActsExamples::MeasurementCalibratorAdapter::calibrate>(&calibrator);
  myActor.updater.connect<&Acts::GainMatrixUpdater::operator()<typename TrackContainer::TrackStateContainerBackend>>(&kfUpdater);

  std::vector<int> particleIds(measurements.size());
  for (auto meas : measurements) {
    for (const auto& [_, barcode] : ActsExamples::makeRange(measurementParticlesMap.equal_range(meas.index()))) {
      particleIds[meas.index()] = barcode.particle();
    }
  }
  myActor.particleIds = &particleIds;

  for (std::size_t iSeed = 0; iSeed < initialParameters.size(); ++iSeed) {
    ACTS_DEBUG("seed " << iSeed);
    tracksTemp.clear();
    auto rootBranch = tracksTemp.makeTrack();
    auto propState = m_propagator.makeState(propOptions);
    auto& r = propState.template get<MyTrackFindingResult>();
    r.collectedTracks.clear();
    r.activeBranches.clear();
    r.activeBranches.push_back(rootBranch);

    auto initResult = m_propagator.initialize(propState, initialParameters.at(iSeed));        
    auto propResult = m_propagator.propagate(propState);
    
    ACTS_VERBOSE("Collected tracks: " << r.collectedTracks.size());
    for (auto collectedTrack : r.collectedTracks) {
      auto track = tracks.makeTrack();
      track.copyFrom(collectedTrack);
      Acts::smoothTrack(ctx.geoContext, track, logger());
      Acts::extrapolateTrackToReferenceSurface(track, *pSurface, extrapolator, extrapolationOptions, 
                                               extrapolationStrategy, logger());
      Acts::calculateTrackQuantities(track);
    }
  }

  // output containers
  auto constTrackStateContainer = std::make_shared<Acts::ConstVectorMultiTrajectory>(std::move(*trackStateContainer));
  auto constTrackContainer = std::make_shared<Acts::ConstVectorTrackContainer>(std::move(*trackContainer));
  ConstTrackContainer constTracks{constTrackContainer, constTrackStateContainer};
  m_outputTracks(ctx, std::move(constTracks));

  return ProcessCode::SUCCESS;
}

}