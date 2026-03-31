#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"

#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Index.hpp"
//#include "ActsExamples/EventData/IndexSourceLink.hpp"

namespace My{

using ActsExamples::ReadDataHandle;
using ActsExamples::WriteDataHandle;
using ActsExamples::IAlgorithm;
using ActsExamples::ProcessCode;
using ActsExamples::AlgorithmContext;
using ActsExamples::ConstTrackContainer;
using ActsExamples::TrackParametersContainer;
using ActsExamples::MeasurementContainer;
using ActsExamples::IndexSourceLinkAccessor;
using ActsExamples::TrackContainer;
using ActsExamples::IndexMultimap;
using ActsExamples::Index;
using ActsExamples::MeasurementParticlesMap;

using TrackProxy = TrackContainer::TrackProxy;
using TrackStateProxy = TrackContainer::TrackStateProxy;
using TrackStateContainerBackend = TrackContainer::TrackStateContainerBackend;
using Updater = typename Acts::KalmanFitterExtensions<TrackStateContainerBackend>::Updater;
using Calibrator = typename Acts::KalmanFitterExtensions<TrackStateContainerBackend>::Calibrator;
using PM = Acts::TrackStatePropMask;

struct MyTrackFindingResult {
  std::vector<TrackStateProxy> trackStateCandidates;
  std::vector<TrackProxy> activeBranches;
  std::vector<TrackProxy> collectedTracks;
};

struct MyTrackFindingActor {
  using result_type = MyTrackFindingResult;
  IndexSourceLinkAccessor slAccessor;
  const MeasurementContainer* measurements{nullptr};
  TrackStateContainerBackend* trackStates{nullptr};
  const Acts::CalibrationContext* calibrationContext{nullptr};
  const std::vector<int>* particleIds;
  Calibrator calibrator{Acts::DelegateFuncTag<Acts::detail::voidFitterCalibrator<TrackStateContainerBackend>>{}};
  Updater updater{Acts::DelegateFuncTag<Acts::detail::voidFitterUpdater<TrackStateContainerBackend>>{}};
  
  template <typename propagator_state_t, typename stepper_t, typename navigator_t>
  void act(propagator_state_t& state, const stepper_t& stepper, const navigator_t& navigator, 
    result_type& result, const Acts::Logger& logger) const {
    if (result.activeBranches.empty()) return;

    auto surface = navigator.currentSurface(state.navigation);
    if (surface != nullptr) {
      // ACTS_DEBUG("On surface: " << surface->geometryId());
      filter(surface, state, stepper, navigator, result, logger);    
    }
    
    if (navigator.endOfWorldReached(state.navigation)){
      // store current branch and remove it from active branches
      auto currentBranch = result.activeBranches.back();
      result.collectedTracks.push_back(currentBranch);
      result.activeBranches.pop_back();
      if (result.activeBranches.empty()) return;
      // set stepper state to the next branch
      currentBranch = result.activeBranches.back();
      ACTS_VERBOSE("Switching navigator to track " << currentBranch.tipIndex());
      auto cs = currentBranch.outermostTrackState(); // current state
      stepper.initialize(state.stepping, cs.filtered(), cs.filteredCovariance(), stepper.particleHypothesis(state.stepping), cs.referenceSurface());
      state.navigation.options.startSurface = &cs.referenceSurface();
      state.navigation.options.targetSurface = nullptr;
      auto navInitRes = navigator.initialize(state.navigation, stepper.position(state.stepping), stepper.direction(state.stepping), state.options.direction);
    }
  }

  template <typename propagator_state_t, typename stepper_t, typename navigator_t>
  void filter(const Acts::Surface* surface, propagator_state_t& state, const stepper_t& stepper, 
    const navigator_t& navigator, result_type& result, const Acts::Logger& logger) const {

    if (surface->associatedDetectorElement() == nullptr) return;

    int layer = surface->geometryId().layer();
    // TODO skip pixel-like surfaces using geo class
    // if (layer== 2) return;
    // if (layer== 6) return;
    // if (layer==13) return;
    // if (layer==20) return;
    // if (layer==27) return;
    // if (layer==34) return;
    // if (layer==38) return;
    auto boundStateRes = stepper.boundState(state.stepping, *surface);
    auto& boundState = *boundStateRes;
    auto& [boundParams, jacobian, pathLength] = boundState;

    // TODO add material

    // create trackState candidates
    result.trackStateCandidates.clear();
    auto [slBegin, slEnd] = slAccessor.range(*surface);
    int passedCandidates = 0;
    double chi2max = 10;

    auto currentBranch = result.activeBranches.back();
    auto tipIndex = currentBranch.tipIndex();
    
    if (slBegin - slEnd>0) result.trackStateCandidates.reserve(slBegin - slEnd);

    char buffer[200]; // Define a fixed-size buffer
    for (auto it = slBegin; it != slEnd; ++it) {
      // if (it != slBegin) continue; // only one hit per surface allowed
      const auto sl = *it; // source link
      auto meas = measurements->getMeasurement(sl.get<ActsExamples::IndexSourceLink>().index());
      auto particleId = particleIds->at(meas.index());
      
      PM mask = PM::Predicted | PM::Jacobian | PM::Calibrated;
      auto ts = trackStates->makeTrackState(mask, tipIndex);
      ts.predicted() = boundParams.parameters();
      ts.predictedCovariance() = *boundParams.covariance();
      ts.jacobian() = jacobian;
      ts.pathLength() = pathLength;
      ts.setReferenceSurface(boundParams.referenceSurface().getSharedPtr());
      calibrator(state.geoContext, *calibrationContext, sl, ts);
      if (ts.calibratedSize()!=1) continue;
      // chi2 for 1-D measurements
      double calib = ts.effectiveCalibrated().data()[0];
      double pred = ts.predicted().data()[0];
      double calibCov = ts.effectiveCalibratedCovariance().data()[0];
      double predCov = ts.predictedCovariance().data()[0];
      double res = calib - pred;
      double chi2 = res*res/(calibCov + predCov);
      std::sprintf(buffer, "track=%2d layer=%2d partId=%2d calib=%4.1f pred=%4.1f cov=%6.4f chi2=%9.5f", 
                   tipIndex, surface->geometryId().layer(), particleId, calib, pred, predCov, chi2);
      ACTS_VERBOSE(buffer);
      // ACTS_DEBUG("track=" << tipIndex << "layer=" << surface->geometryId().layer() << 
      //            "particleId=" << particleId << " calib=" << calib <<" pred=" << pred << " predCov=" << predCov << " chi2=" << chi2);
      ts.chi2() = chi2;
      result.trackStateCandidates.push_back(ts);
      if (ts.chi2() > chi2max) continue;
      passedCandidates++;
    }
    // ACTS_INFO("Track " << tipIndex << " on surface: " << surface->geometryId() << " candidates=" << result.trackStateCandidates.size());
    if (passedCandidates==0) return;
    bool isOutlier = (passedCandidates==0);
    // copy new track state indices from trackStateCandidates to the trackStates container
    // TODO sort track candidates according to chi2
    std::vector<int> trackStateList;
    for (size_t i=0; i<result.trackStateCandidates.size(); i++){
      auto& ts = result.trackStateCandidates[i];
      if (ts.chi2() > chi2max) continue;
      PM mask = PM::Predicted | PM::Filtered | PM::Jacobian | PM::Calibrated;
      auto trackState = trackStates->makeTrackState(mask, ts.previous());
      trackState.copyFrom(ts, mask, false);
      trackState.typeFlags().set(Acts::TrackStateFlag::ParameterFlag);
      trackState.typeFlags().set(Acts::TrackStateFlag::MeasurementFlag);
      trackStateList.push_back(trackState.index());
    }

    auto rootBranch = result.activeBranches.back();
    std::vector<TrackProxy> newBranches;
    for (int i=0;i<trackStateList.size();i++){
      auto shallowCopy = [&] {
        auto sc = rootBranch.container().makeTrack();
        sc.copyFromShallow(rootBranch);
        return sc;
      };
      auto newBranch = (i==0) ? rootBranch : shallowCopy();
      newBranch.tipIndex() = trackStateList[i];
      newBranches.push_back(newBranch);
    }
    result.activeBranches.pop_back();
    
    for (TrackProxy newBranch : newBranches) {
      auto trackState = newBranch.outermostTrackState();
      if (trackState.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
        auto updateRes = updater(state.geoContext, trackState, logger);
        // TODO: add branch stopper
      }
      result.activeBranches.push_back(newBranch);
    }

    currentBranch = result.activeBranches.back();
    auto cs = currentBranch.outermostTrackState(); // current state
    auto freePar = Acts::MultiTrajectoryHelpers::freeFiltered(state.geoContext, cs);
    stepper.update(state.stepping, freePar, cs.filtered(), cs.filteredCovariance(), *surface);
    // TODO add material effects
  }
};


class MyTrackFindingAlgorithm final : public IAlgorithm {
 public:
  class Config {
   public:
    std::string inputMeasurements;
    std::string inputMeasurementParticlesMap;
    std::string inputInitialTrackParameters;
    std::string outputTracks;
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField;
  };

  MyTrackFindingAlgorithm(Config config, Acts::Logging::Level level);
  ProcessCode execute(const AlgorithmContext& ctx) const override;
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  
  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this, "InputMeasurements"};
  ReadDataHandle<MeasurementParticlesMap> m_inputMeasurementParticlesMap{this, "InputMeasurementParticlesMap"};
  ReadDataHandle<TrackParametersContainer> m_inputInitialTrackParameters{this, "InputInitialTrackParameters"};
  WriteDataHandle<ConstTrackContainer> m_outputTracks{this, "OutputTracks"};
};

} // namespace My