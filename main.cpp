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
#include "ActsExamples/Io/Root/RootSeedWriter.hpp"
#include "ActsExamples/Io/Root/RootMeasurementWriter.hpp"
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

#include <filesystem>

#include "MyRefittingAlgorithm.hpp"

using Acts::UnitConstants::cm;

int main(int argc, char *argv[]){
  // default parameters
  TString inputDir = "none";
  TString outputDir = "test";
  int nEvents = 10;
  Acts::PdgParticle pdgCode = Acts::eProton;
  double etaMin = 1.6;
  double vzMax = 50;
  double radLengthPerSeed = 0.01;

  // double ptMin = 0.100_GeV;
  // double ptMax = 1.000_GeV+1.e-100_GeV;

  double ptMin = 0.350_GeV;
  double ptMax = 0.350_GeV +1e-100_GeV;
  // non-default setup
  if (argc>=4) {
    inputDir = argv[1];
    outputDir = argv[2];
    nEvents = TString(argv[3]).Atoi();
  }
  
  if (argc>=6) {
    TString pdg = TString(argv[4]);
    if (pdg.Contains("pi")) pdgCode = Acts::ePionPlus;
    if (pdg.Contains("pr")) pdgCode = Acts::eProton;
    if (pdg.Contains("mu")) pdgCode = Acts::eAntiMuon;
    etaMin = TString(argv[5]).Atof();
  }
  if (argc>=7) vzMax = TString(argv[6]).Atof();
  if (argc>=8) radLengthPerSeed = TString(argv[7]).Atof();

  bool istpc = !outputDir.Contains("notpc");
  bool isframe = !outputDir.Contains("noframe");

  printf("Running acts: %d eta=%.1f vzMax=%.0f\n", pdgCode, etaMin, vzMax);

  inputDir.Append("/");
  outputDir.Append("/");
  std::filesystem::create_directory(outputDir.Data());

  // Logger
  Acts::Logging::Level logLevel = Acts::Logging::INFO;
  Acts::Logging::Level logLevelFatras = Acts::Logging::INFO;
  Acts::Logging::Level logLevelDigi = Acts::Logging::INFO;
  Acts::Logging::Level logLevelSeed = Acts::Logging::INFO;
  Acts::Logging::Level logLevelFinder = Acts::Logging::INFO;
  Acts::Logging::Level logLevelSequencer = Acts::Logging::ERROR;
  Acts::Logging::Level logLevelMatcher = Acts::Logging::INFO;
  Acts::Logging::Level logLevelMeasWriter = Acts::Logging::INFO;
  Acts::Logging::Level logLevelMyRefit = Acts::Logging::VERBOSE;
  // collection names
  std::string particles = "particles";
  std::string vertices = "vertices";
  std::string tracks = "tracks";
  std::string simhits = "simhits";
  std::string spacepoints = "spacepoints";
  std::string measurements = "measurements";
  std::string seeds = "seeds";

  // Random number generator config
  auto rnd = std::make_shared<ActsExamples::RandomNumbers>(ActsExamples::RandomNumbers::Config({42}));

  // Particle gun generator config
  double etaMax = etaMin + 1.e-100;

  ActsExamples::ParametricParticleGenerator::Config genCfg;
  genCfg.etaUniform = true;
  genCfg.thetaMin = 2 * atan(exp(-etaMax));
  genCfg.thetaMax = 2 * atan(exp(-etaMin));
  genCfg.pMin = ptMin;
  genCfg.pMax = ptMax;
  genCfg.pTransverse = true;

  genCfg.pdg = pdgCode;
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
  particleReaderCfg.filePath = TString(inputDir+"particles.root").Data();

  auto trackingGeometryPtr = CreateTrackingGeometry(istpc, isframe);
  auto trackingGeometry = std::make_shared<Acts::TrackingGeometry>(*trackingGeometryPtr);

  // Fatras config
  ActsExamples::FatrasSimulation::Config fatrasCfg;
  fatrasCfg.inputParticles = particles;
  fatrasCfg.outputParticlesInitial = "initial";
  fatrasCfg.outputParticlesFinal = "final";
  fatrasCfg.outputSimHits = simhits;
  fatrasCfg.trackingGeometry = trackingGeometry;
  fatrasCfg.pMin = 0.1_GeV;
  fatrasCfg.magneticField = std::make_shared<Acts::ConstantBField>(Acts::Vector3(0.0, 0.0, bz*Acts::UnitConstants::T));
  fatrasCfg.randomNumbers = rnd;
  //fatrasCfg.generateHitsOnMaterial = true;

  // Digitization config
  ActsExamples::DigitizationConfig digiCfg(0, 0.001, 0);  //doMerge (bool), mergeSigma (mm), mergeCommonCorner (bool)
  digiCfg.inputSimHits = simhits;
  digiCfg.randomNumbers = rnd;
  digiCfg.outputMeasurements = measurements;
  digiCfg.trackingGeometry = trackingGeometry;
  digiCfg.digitizationConfigs = ActsExamples::readDigiConfigFromJson("digi-smearing-config.json");

  // Create space points
  ActsExamples::SpacePointMaker::Config spCfg;
  spCfg.inputSourceLinks = digiCfg.outputSourceLinks;
  spCfg.inputMeasurements = measurements;
  spCfg.trackingGeometry = trackingGeometry;
  spCfg.outputSpacePoints = spacepoints;
  spCfg.geometrySelection = {Acts::GeometryIdentifier{}};

  // ActsExamples::TruthSeedingAlgorithm::Config truthSeedingCfg;
  // truthSeedingCfg.inputParticles = evgenCfg.outputParticles;
  // truthSeedingCfg.inputSpacePoints = {spacepoints};
  // truthSeedingCfg.inputMeasurementParticlesMap = digiCfg.outputMeasurementParticlesMap;
  // truthSeedingCfg.outputSeeds = seeds;
  // truthSeedingCfg.outputParticles = "seeded_particles";
  // truthSeedingCfg.outputProtoTracks = "prototracks";
  // truthSeedingCfg.deltaRMin = 1._mm;
  // truthSeedingCfg.deltaRMax = 100._mm;

  // Seeding
  ActsExamples::SeedingAlgorithm::Config seedingCfg;
  seedingCfg.inputSpacePoints = {spacepoints};
  seedingCfg.outputSeeds = seeds;
  seedingCfg.seedFinderOptions.bFieldInZ = bz*Acts::UnitConstants::T;
  seedingCfg.seedFinderConfig.minPt              = 0.1_GeV;
  seedingCfg.seedFinderConfig.deltaZMax          = (positions[1] - positions[0])*cm + 1_mm;
  seedingCfg.seedFinderConfig.deltaRMax          = (1-positions[1]/positions[2])*rMaxStation*cm + 1_mm;
  seedingCfg.seedFinderConfig.rMin               = rMinStation*cm;
  seedingCfg.seedFinderConfig.rMax               = rMaxStation*cm;
  seedingCfg.seedFinderConfig.rMinMiddle         = rMinStation*cm;
  seedingCfg.seedFinderConfig.rMaxMiddle         = rMaxStation*cm;
  seedingCfg.seedFinderConfig.zMin               = positions[0]*cm - 1_mm;
  seedingCfg.seedFinderConfig.zMax               = positions[2]*cm + 1_mm;
  seedingCfg.seedFinderConfig.cotThetaMax        = 7.0;  
  seedingCfg.seedFinderConfig.impactMax          = rMinStation*cm - 10_mm;
  seedingCfg.seedFinderConfig.collisionRegionMin = -vzMax*cm; // important at low momenta due to mult scattering effects
  seedingCfg.seedFinderConfig.collisionRegionMax = +vzMax*cm; // important at low momenta due to mult scattering effects
  seedingCfg.seedFinderConfig.radLengthPerSeed   = radLengthPerSeed;

  // copy relevant options to grid config
  seedingCfg.gridConfig.rMax        = seedingCfg.seedFinderConfig.rMax;
  seedingCfg.gridConfig.zMin        = seedingCfg.seedFinderConfig.zMin;
  seedingCfg.gridConfig.zMax        = seedingCfg.seedFinderConfig.zMax;
  seedingCfg.gridConfig.cotThetaMax = seedingCfg.seedFinderConfig.cotThetaMax;
  seedingCfg.gridConfig.deltaRMax   = seedingCfg.seedFinderConfig.deltaRMax;
  seedingCfg.gridConfig.impactMax   = seedingCfg.seedFinderConfig.impactMax;
  seedingCfg.gridConfig.minPt       = seedingCfg.seedFinderConfig.minPt;
  seedingCfg.gridOptions.bFieldInZ  = seedingCfg.seedFinderOptions.bFieldInZ;

  seedingCfg.seedFilterConfig.maxSeedsPerSpM = seedingCfg.seedFinderConfig.maxSeedsPerSpM;

  // Parameter Estimation from seeds
  ActsExamples::TrackParamsEstimationAlgorithm::Config paramsEstimationCfg;
  paramsEstimationCfg.inputSeeds = seeds;
  paramsEstimationCfg.outputTrackParameters = "estimatedparameters";
  paramsEstimationCfg.trackingGeometry = trackingGeometry;
  paramsEstimationCfg.magneticField = fatrasCfg.magneticField;

  // Track finding
  ActsExamples::TrackFindingAlgorithm::Config trackFindingCfg;
  trackFindingCfg.inputMeasurements = measurements;
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
  particleWriterCfg.filePath = TString(outputDir+"particles.root").Data();

  // SimhitReader config
  ActsExamples::RootSimHitReader::Config simhitReaderCfg;
  simhitReaderCfg.outputSimHits = simhits;
  simhitReaderCfg.filePath = TString(inputDir+"hits.root").Data();

  // SimhitWriter config
  ActsExamples::RootSimHitWriter::Config simhitWriterCfg;
  simhitWriterCfg.inputSimHits = simhits;
  simhitWriterCfg.filePath = TString(outputDir+"hits.root").Data();

  ActsExamples::RootMeasurementWriter::Config measWriterCfg;
  measWriterCfg.inputMeasurements = measurements;
  measWriterCfg.inputSimHits = simhits;
  measWriterCfg.inputMeasurementSimHitsMap = digiCfg.outputMeasurementSimHitsMap;
  measWriterCfg.filePath = TString(outputDir+"measurements.root").Data();
  measWriterCfg.trackingGeometry = trackingGeometry;
  measWriterCfg.boundIndices = Acts::GeometryHierarchyMap<std::vector<Acts::BoundIndices>>(digiCfg.getBoundIndices());

  // SpacepointWriter config
  ActsExamples::RootSpacepointWriter::Config spWriterCfg;
  spWriterCfg.inputSpacepoints = spacepoints;
  spWriterCfg.filePath = TString(outputDir+"spacepoints.root").Data();

  ActsExamples::RootSeedWriter::Config seedWriterCfg;
  seedWriterCfg.inputSeeds = seeds;
  seedWriterCfg.filePath = TString(outputDir+"seeds.root").Data();

/*
  // Measurement writer
  ActsExamples::CsvMeasurementWriter::Config measCsvWriterCfg;
  measCsvWriterCfg.inputMeasurements = digiCfg.outputMeasurements;
  measCsvWriterCfg.inputMeasurementSimHitsMap = digiCfg.outputMeasurementSimHitsMap;
  measCsvWriterCfg.outputDir = "measurements";
*/

/*
  // Seed csv writer config
  ActsExamples::CsvSeedWriter::Config seedCsvWriterCfg;
  seedCsvWriterCfg.inputSimHits = simhits;
  seedCsvWriterCfg.inputSimSeeds = seeds;
  seedCsvWriterCfg.inputTrackParameters = paramsEstimationCfg.outputTrackParameters;
  seedCsvWriterCfg.inputMeasurementParticlesMap = digiCfg.outputMeasurementParticlesMap;
  seedCsvWriterCfg.inputMeasurementSimHitsMap = digiCfg.outputMeasurementSimHitsMap;
  seedCsvWriterCfg.outputDir = "seeds";
*/

  // write track states from CKF
  ActsExamples::RootTrackStatesWriter::Config trackStatesWriterCfg;
  trackStatesWriterCfg.inputTracks = tracks;
  trackStatesWriterCfg.inputParticles = particles;
  trackStatesWriterCfg.inputSimHits = simhits;
  trackStatesWriterCfg.inputTrackParticleMatching = trackTruthMatcherCfg.outputTrackParticleMatching;
  trackStatesWriterCfg.inputMeasurementSimHitsMap = digiCfg.outputMeasurementSimHitsMap;
  trackStatesWriterCfg.treeName = "trackstates";
  trackStatesWriterCfg.filePath = TString(outputDir+"trackstates.root").Data();

  ActsExamples::RootTrackSummaryWriter::Config trackSummaryWriterCfg;
  trackSummaryWriterCfg.inputTracks = tracks;
  trackSummaryWriterCfg.inputParticles = particles;
  trackSummaryWriterCfg.inputTrackParticleMatching = trackTruthMatcherCfg.outputTrackParticleMatching;
  trackSummaryWriterCfg.filePath = TString(outputDir+"tracksummary.root").Data();

  ActsExamples::MyRefittingAlgorithm::Config refitCfg;
  refitCfg.inputTracks = tracks;
  refitCfg.outputTracks = "refitted_tracks";
  refitCfg.trackingGeometry = trackingGeometry;
  refitCfg.magneticField = fatrasCfg.magneticField;
  //refitCfg.fit = ActsExamples::makeKalmanFitterFunction(trackingGeometry, fatrasCfg.magneticField);

  refitCfg.fit = ActsExamples::makeKalmanFitterFunction(trackingGeometry, fatrasCfg.magneticField, true, true, 0.0, Acts::FreeToBoundCorrection(), *Acts::getDefaultLogger("Kalman", logLevelMyRefit));

  // Acts::FreeToBoundCorrection
  // refitCfg.fit = ActsExamples::makeGsfFitterFunction(trackingGeometry, fatrasCfg.magneticField);
  //   /*BetheHeitlerApprox betheHeitlerApprox*/ , // acts.examples.AtlasBetheHeitlerApprox.makeDefault(),
  //   /*std::size_t maxComponents*/ 4,
  //   /*double weightCutoff*/ 1.0e-4, 
  //   /*Acts::ComponentMergeMethod componentMergeMethod*/,
  //   /*MixtureReductionAlgorithm mixtureReductionAlgorithm*/,
  //   /*const Acts::Logger& logger*/
  
  // Sequencer config
  ActsExamples::Sequencer::Config sequencerCfg;
  sequencerCfg.numThreads = 1;
  sequencerCfg.events = nEvents;
  sequencerCfg.logLevel = logLevelSequencer;

  // Start sequencer
  ActsExamples::Sequencer sequencer(sequencerCfg);
  
  if (inputDir.Contains("none")){ // particle gun + fartras simulation
    sequencer.addReader(std::make_shared<ActsExamples::EventGenerator>(evgenCfg, logLevel));
    sequencer.addElement(std::make_shared<ActsExamples::FatrasSimulation>(fatrasCfg, logLevelFatras));
  } else { // read particles and hits from input file
    sequencer.addReader(std::make_shared<ActsExamples::RootParticleReader>(particleReaderCfg, logLevel));
    sequencer.addReader(std::make_shared<ActsExamples::RootSimHitReader>(simhitReaderCfg, logLevel));
  }
  sequencer.addAlgorithm(std::make_shared<ActsExamples::DigitizationAlgorithm>(digiCfg, logLevelDigi));
  sequencer.addAlgorithm(std::make_shared<ActsExamples::SpacePointMaker>(spCfg, logLevel));
  sequencer.addAlgorithm(std::make_shared<ActsExamples::SeedingAlgorithm>(seedingCfg, logLevelSeed));
  sequencer.addAlgorithm(std::make_shared<ActsExamples::TrackParamsEstimationAlgorithm>(paramsEstimationCfg, logLevel));
  sequencer.addAlgorithm(std::make_shared<ActsExamples::TrackFindingAlgorithm>(trackFindingCfg, logLevelFinder));
  sequencer.addAlgorithm(std::make_shared<ActsExamples::MyRefittingAlgorithm>(refitCfg, logLevelMyRefit));

  sequencer.addAlgorithm(std::make_shared<ActsExamples::TrackTruthMatcher>(trackTruthMatcherCfg, logLevelMatcher));


  sequencer.addWriter(std::make_shared<ActsExamples::RootParticleWriter>(particleWriterCfg, logLevel));
  sequencer.addWriter(std::make_shared<ActsExamples::RootSimHitWriter>(simhitWriterCfg, logLevel));
  sequencer.addWriter(std::make_shared<ActsExamples::RootMeasurementWriter>(measWriterCfg, logLevelMeasWriter));
  sequencer.addWriter(std::make_shared<ActsExamples::RootSpacepointWriter>(spWriterCfg, logLevel));
  sequencer.addWriter(std::make_shared<ActsExamples::RootSeedWriter>(seedWriterCfg, logLevel));
  sequencer.addWriter(std::make_shared<ActsExamples::RootTrackStatesWriter>(trackStatesWriterCfg, logLevel));
  sequencer.addWriter(std::make_shared<ActsExamples::RootTrackSummaryWriter>(trackSummaryWriterCfg, logLevel));

  //  obsolete
  //  sequencer.addAlgorithm(std::make_shared<ActsExamples::TruthSeedingAlgorithm>(truthSeedingCfg, logLevelSeed));
  //  sequencer.addWriter(std::make_shared<ActsExamples::CsvMeasurementWriter>(measCsvWriterCfg, logLevel));
  //  sequencer.addWriter(std::make_shared<ActsExamples::CsvSeedWriter>(seedCsvWriterCfg, logLevel));

  sequencer.run();

  return 0;
}
