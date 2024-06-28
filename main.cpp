#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Seeding/SeedFinderConfig.hpp"
#include "ActsExamples/Io/Json/JsonDigitizationConfig.hpp"
#include "ActsExamples/Generators/EventGenerator.hpp"
#include "ActsExamples/Generators/VertexGenerators.hpp"
#include "ActsExamples/Generators/MultiplicityGenerators.hpp"
#include "ActsExamples/Generators/ParametricParticleGenerator.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Fatras/FatrasSimulation.hpp"
#include "ActsExamples/TelescopeDetector/TelescopeDetector.hpp"
// #include "ActsExamples/Geometry/MaterialWiper.hpp"
#include "ActsExamples/Digitization/DigitizationAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"
#include "ActsExamples/Io/Root/RootSimHitReader.hpp"
#include "ActsExamples/Io/Root/RootSimHitWriter.hpp"
#include "ActsExamples/Io/Root/RootSpacepointWriter.hpp"
#include "ActsExamples/Io/Csv/CsvSeedWriter.hpp"
#include "ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp"
#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"
#include "ActsExamples/Io/Root/RootTrackStatesWriter.hpp"
#include "ActsExamples/Io/Root/RootTrackSummaryWriter.hpp"
#include "ActsExamples/Io/Root/RootParticleWriter.hpp"
#include "ActsExamples/Io/Root/RootParticleReader.hpp"
#include "ActsExamples/Io/Csv/CsvMeasurementWriter.hpp"
#include "ActsExamples/TruthTracking/TrackTruthMatcher.hpp"
// Temporary
#include "ActsExamples/TruthTracking/TruthSeedingAlgorithm.hpp"

// Needed for detector construction
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Geometry/DiscLayer.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/NavigationLayer.hpp"
//  native units: mm, GeV, c=1, e=1

#include "tracker_config.h"
#include "tracker.h"

int main(){
  // Logger
  Acts::Logging::Level logLevel = Acts::Logging::INFO;
  Acts::Logging::Level logLevelFatras = Acts::Logging::INFO;
  Acts::Logging::Level logLevelDigi = Acts::Logging::INFO;
  Acts::Logging::Level logLevelSeed = Acts::Logging::INFO;
  Acts::Logging::Level logLevelFinder = Acts::Logging::INFO;
  Acts::Logging::Level logLevelSequencer = Acts::Logging::ERROR;
  Acts::Logging::Level logLevelMatcher = Acts::Logging::VERBOSE;
 
  // collection names
  std::string particles = "particles";
  std::string vertices = "vertices";
  std::string tracks = "tracks";
  std::string simhits = "simhits";
  std::string seeds = "seeds";

  // Random number generator config
  auto rnd = std::make_shared<ActsExamples::RandomNumbers>(ActsExamples::RandomNumbers::Config({42}));

  int nEvents = 10000;

  // Particle gun generator config
  double etaMin = 2.0;
  double etaMax = 2.0+1.e-100;

  ActsExamples::ParametricParticleGenerator::Config genCfg;
  genCfg.etaUniform = true;
  genCfg.thetaMin = 2 * atan(exp(-etaMax));
  genCfg.thetaMax = 2 * atan(exp(-etaMin));
  genCfg.pMin = 0.2_GeV;
  genCfg.pMax = 0.2_GeV+1.e-100_GeV;
  genCfg.pTransverse = true;
  genCfg.pdg = Acts::PdgParticle::eAntiMuon;
  ActsExamples::EventGenerator::Generator gen{
      std::make_shared<ActsExamples::FixedMultiplicityGenerator>(1),
      std::make_shared<ActsExamples::FixedPrimaryVertexPositionGenerator>(),
      std::make_shared<ActsExamples::ParametricParticleGenerator>(genCfg)};
  ActsExamples::EventGenerator::Config evgenCfg;
  evgenCfg.outputParticles = particles;
  evgenCfg.outputVertices = vertices;
  evgenCfg.generators = {gen};
  evgenCfg.randomNumbers = rnd;

  // Particle reader config
  ActsExamples::RootParticleReader::Config particleReaderCfg;
  particleReaderCfg.outputParticles = particles;
  particleReaderCfg.treeName = "particles";
  particleReaderCfg.filePath = "particles.root";

  auto trackingGeometryPtr = CreateTrackingGeometry();
  auto trackingGeometry = std::make_shared<Acts::TrackingGeometry>(*trackingGeometryPtr);

  // Fatras config
  ActsExamples::FatrasSimulation::Config fatrasCfg;
  fatrasCfg.inputParticles = particles;
  fatrasCfg.outputParticlesInitial = "initial";
  fatrasCfg.outputParticlesFinal = "final";
  fatrasCfg.outputSimHits = simhits;
  fatrasCfg.trackingGeometry = trackingGeometry;
  fatrasCfg.magneticField = std::make_shared<Acts::ConstantBField>(Acts::Vector3(0.0, 0.0, bz));
  fatrasCfg.randomNumbers = rnd;

  // Digitization config
  bool doMerge = 0;
  double mergeSigma = 0.001; // mm
  bool mergeCommonCorner = 0;
  ActsExamples::DigitizationConfig digiCfg(doMerge, mergeSigma, mergeCommonCorner);
  digiCfg.inputSimHits = simhits;
  digiCfg.randomNumbers = rnd;
  digiCfg.trackingGeometry = trackingGeometry;
  digiCfg.digitizationConfigs = ActsExamples::readDigiConfigFromJson("digi-smearing-config.json");

  // Create space points
  ActsExamples::SpacePointMaker::Config spCfg;
  spCfg.inputSourceLinks = digiCfg.outputSourceLinks;
  spCfg.inputMeasurements = digiCfg.outputMeasurements;
  spCfg.trackingGeometry = trackingGeometry;
  spCfg.outputSpacePoints = "spacepoints";
  spCfg.geometrySelection = {Acts::GeometryIdentifier{}};

  // ActsExamples::TruthSeedingAlgorithm::Config truthSeedingCfg;
  // truthSeedingCfg.inputParticles = evgenCfg.outputParticles;
  // truthSeedingCfg.inputSpacePoints = {spCfg.outputSpacePoints};
  // truthSeedingCfg.inputMeasurementParticlesMap = digiCfg.outputMeasurementParticlesMap;
  // truthSeedingCfg.outputSeeds = seeds;
  // truthSeedingCfg.outputParticles = "seeded_particles";
  // truthSeedingCfg.outputProtoTracks = "prototracks";
  // truthSeedingCfg.deltaRMin = 1._mm;
  // truthSeedingCfg.deltaRMax = 100._mm;

  // Seeding
  ActsExamples::SeedingAlgorithm::Config seedingCfg;
  seedingCfg.inputSpacePoints = {spCfg.outputSpacePoints};
  seedingCfg.outputSeeds = seeds;
  seedingCfg.gridOptions.bFieldInZ = bz;
  seedingCfg.gridConfig.minPt = 0.1_GeV;
  seedingCfg.gridConfig.zMin = positions[0] - 1_mm;
  seedingCfg.gridConfig.zMax = positions[2] + 1_mm;
  seedingCfg.gridConfig.rMax = rMax;
  seedingCfg.gridConfig.deltaRMax = rMax - rMin;
  seedingCfg.gridConfig.cotThetaMax = 7.0;
  seedingCfg.gridConfig.impactMax = 300_mm;

  seedingCfg.seedFinderOptions.bFieldInZ = seedingCfg.gridOptions.bFieldInZ;
  seedingCfg.seedFinderConfig.impactMax = seedingCfg.gridConfig.impactMax;
  seedingCfg.seedFinderConfig.rMax = seedingCfg.gridConfig.rMax;
  seedingCfg.seedFinderConfig.deltaRMax = seedingCfg.gridConfig.deltaRMax;
  seedingCfg.seedFinderConfig.deltaRMaxTopSP = seedingCfg.gridConfig.deltaRMax;
  seedingCfg.seedFinderConfig.deltaRMaxBottomSP = seedingCfg.gridConfig.deltaRMax;
  seedingCfg.seedFinderConfig.zMin = seedingCfg.gridConfig.zMin;
  seedingCfg.seedFinderConfig.zMax = seedingCfg.gridConfig.zMax;
  seedingCfg.seedFinderConfig.cotThetaMax = seedingCfg.gridConfig.cotThetaMax;
  seedingCfg.seedFinderConfig.minPt = seedingCfg.gridConfig.minPt;
  seedingCfg.seedFinderConfig.maxSeedsPerSpM = seedingCfg.seedFilterConfig.maxSeedsPerSpM;
  seedingCfg.seedFinderConfig.deltaZMax = positions[2] - positions[0];
  seedingCfg.seedFinderConfig.rMinMiddle = rMin;
  seedingCfg.seedFinderConfig.rMaxMiddle = rMax;

  // Parameter Estimation from seeds
  ActsExamples::TrackParamsEstimationAlgorithm::Config paramsEstimationCfg;
  paramsEstimationCfg.inputSeeds = seeds;
  paramsEstimationCfg.outputTrackParameters = "estimatedparameters";
  paramsEstimationCfg.trackingGeometry = trackingGeometry;
  paramsEstimationCfg.magneticField = fatrasCfg.magneticField;

  // Track finding
  ActsExamples::TrackFindingAlgorithm::Config trackFindingCfg;
  trackFindingCfg.inputMeasurements = digiCfg.outputMeasurements;
  trackFindingCfg.inputSourceLinks = digiCfg.outputSourceLinks;
  trackFindingCfg.inputInitialTrackParameters = paramsEstimationCfg.outputTrackParameters;
  trackFindingCfg.outputTracks = tracks;
  trackFindingCfg.trackingGeometry = trackingGeometry;
  trackFindingCfg.magneticField = fatrasCfg.magneticField;
  trackFindingCfg.measurementSelectorCfg = {{Acts::GeometryIdentifier(), {{}, {std::numeric_limits<double>::max()}, {1u}}}}; // chi2cut, numberOfMeasurementsPerSurface
  trackFindingCfg.findTracks = ActsExamples::TrackFindingAlgorithm::makeTrackFinderFunction(
      trackingGeometry, fatrasCfg.magneticField, *Acts::getDefaultLogger("TrackFinder", logLevelFinder));

  // Track truth matcher
  ActsExamples::TrackTruthMatcher::Config trackTruthMatcherCfg;
  trackTruthMatcherCfg.inputTracks = tracks;
  trackTruthMatcherCfg.inputParticles = particles;
  trackTruthMatcherCfg.inputMeasurementParticlesMap = digiCfg.outputMeasurementParticlesMap;
  trackTruthMatcherCfg.outputTrackParticleMatching = "track_particle_matching";
  trackTruthMatcherCfg.outputParticleTrackMatching = "particle_track_matching";

  // Particle writer config
  ActsExamples::RootParticleWriter::Config particleWriterCfg;
  particleWriterCfg.inputParticles = particles;
  particleWriterCfg.treeName = "particles";
  particleWriterCfg.filePath = "particles.root";

  // SimhitWriter config
  ActsExamples::RootSimHitWriter::Config simhitWriterCfg;
  simhitWriterCfg.inputSimHits = simhits;
  simhitWriterCfg.filePath = "hits.root";

  ActsExamples::RootSimHitReader::Config simhitReaderCfg;
  simhitReaderCfg.outputSimHits = simhits;
  simhitReaderCfg.filePath = "hits.root";

  // Measurement writer
  ActsExamples::CsvMeasurementWriter::Config measWriterCfg;
  measWriterCfg.inputMeasurements = digiCfg.outputMeasurements;
  measWriterCfg.inputMeasurementSimHitsMap = digiCfg.outputMeasurementSimHitsMap;
  measWriterCfg.outputDir = "measurements";

  // SpacepointWriter config
  ActsExamples::RootSpacepointWriter::Config spWriterCfg;
  spWriterCfg.inputSpacepoints = spCfg.outputSpacePoints;
  spWriterCfg.filePath = "spacepoints.root";

  // Seed writer config
  ActsExamples::CsvSeedWriter::Config seedWriterCfg;
  seedWriterCfg.inputSimHits = simhits;
  seedWriterCfg.inputSimSeeds = seeds;
  seedWriterCfg.inputTrackParameters = paramsEstimationCfg.outputTrackParameters;
  seedWriterCfg.inputMeasurementParticlesMap = digiCfg.outputMeasurementParticlesMap;
  seedWriterCfg.inputMeasurementSimHitsMap = digiCfg.outputMeasurementSimHitsMap;
  seedWriterCfg.outputDir = "seeds";

  // write track states from CKF
  ActsExamples::RootTrackStatesWriter::Config trackStatesWriterCfg;
  trackStatesWriterCfg.inputTracks = trackFindingCfg.outputTracks;
  trackStatesWriterCfg.inputParticles = particles;
  trackStatesWriterCfg.inputSimHits = simhits;
  trackStatesWriterCfg.inputTrackParticleMatching = trackTruthMatcherCfg.outputTrackParticleMatching;
  trackStatesWriterCfg.inputMeasurementSimHitsMap = digiCfg.outputMeasurementSimHitsMap;
  trackStatesWriterCfg.filePath = "trackstates_ckf.root";
  trackStatesWriterCfg.treeName = "trackstates";

  ActsExamples::RootTrackSummaryWriter::Config trackSummaryWriterCfg;
  trackSummaryWriterCfg.inputTracks = trackFindingCfg.outputTracks;
  trackSummaryWriterCfg.inputParticles = particles;
  trackSummaryWriterCfg.inputTrackParticleMatching = trackTruthMatcherCfg.outputTrackParticleMatching;

  // Sequencer config
  ActsExamples::Sequencer::Config sequencerCfg;
  sequencerCfg.numThreads = 1;
  sequencerCfg.events = nEvents;
  sequencerCfg.logLevel = logLevelSequencer;

  // Start sequencer
  ActsExamples::Sequencer sequencer(sequencerCfg);
  
  sequencer.addReader(std::make_shared<ActsExamples::EventGenerator>(evgenCfg, logLevel));
  sequencer.addWriter(std::make_shared<ActsExamples::RootParticleWriter>(particleWriterCfg, logLevel));
//  sequencer.addReader(std::make_shared<ActsExamples::RootParticleReader>(particleReaderCfg, logLevel));
  
  sequencer.addElement(std::make_shared<ActsExamples::FatrasSimulation>(fatrasCfg, logLevelFatras));
  sequencer.addWriter(std::make_shared<ActsExamples::RootSimHitWriter>(simhitWriterCfg, logLevel));
//  sequencer.addReader(std::make_shared<ActsExamples::RootSimHitReader>(simhitReaderCfg, logLevel));

  sequencer.addAlgorithm(std::make_shared<ActsExamples::DigitizationAlgorithm>(digiCfg, logLevelDigi));
  sequencer.addWriter(std::make_shared<ActsExamples::CsvMeasurementWriter>(measWriterCfg, logLevel));
  sequencer.addAlgorithm(std::make_shared<ActsExamples::SpacePointMaker>(spCfg, logLevel));
  sequencer.addWriter(std::make_shared<ActsExamples::RootSpacepointWriter>(spWriterCfg, logLevel));
  sequencer.addAlgorithm(std::make_shared<ActsExamples::SeedingAlgorithm>(seedingCfg, logLevelSeed));

  // sequencer.addAlgorithm(std::make_shared<ActsExamples::TruthSeedingAlgorithm>(truthSeedingCfg, logLevelSeed));

  sequencer.addAlgorithm(std::make_shared<ActsExamples::TrackParamsEstimationAlgorithm>(paramsEstimationCfg, logLevel));
  sequencer.addWriter(std::make_shared<ActsExamples::CsvSeedWriter>(seedWriterCfg, logLevel));
  sequencer.addAlgorithm(std::make_shared<ActsExamples::TrackFindingAlgorithm>(trackFindingCfg, logLevelFinder));
  sequencer.addAlgorithm(std::make_shared<ActsExamples::TrackTruthMatcher>(trackTruthMatcherCfg, logLevelMatcher));
  sequencer.addWriter(std::make_shared<ActsExamples::RootTrackStatesWriter>(trackStatesWriterCfg, logLevel));
  sequencer.addWriter(std::make_shared<ActsExamples::RootTrackSummaryWriter>(trackSummaryWriterCfg, logLevel));
  sequencer.run();

  return 0;
}
