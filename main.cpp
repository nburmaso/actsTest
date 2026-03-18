#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/Generators/EventGenerator.hpp"
#include "ActsExamples/Utilities/VertexGenerators.hpp"
#include "ActsExamples/Utilities/MultiplicityGenerators.hpp"
#include "ActsExamples/Utilities/ParametricParticleGenerator.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Fatras/FatrasSimulation.hpp"
#include "ActsExamples/Digitization/DigitizationAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"
#include "ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp"
#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"
#include "ActsExamples/Io/HepMC3/HepMC3InputConverter.hpp"
#include "ActsExamples/Io/Root/RootParticleReader.hpp"
#include "ActsExamples/Io/Root/RootParticleWriter.hpp"
#include "ActsExamples/Io/Root/RootSimHitReader.hpp"
#include "ActsExamples/Io/Root/RootSimHitWriter.hpp"
#include "ActsExamples/Io/Root/RootMeasurementWriter.hpp"
#include "ActsExamples/Io/Root/RootSpacepointWriter.hpp"
#include "ActsExamples/Io/Root/RootSeedWriter.hpp"
#include "ActsExamples/Io/Root/RootTrackStatesWriter.hpp"
#include "ActsExamples/Io/Root/RootTrackSummaryWriter.hpp"
#include "ActsExamples/TruthTracking/TrackTruthMatcher.hpp"
#include "ActsExamples/TruthTracking/TruthSeedingAlgorithm.hpp"
#include "ActsExamples/AmbiguityResolution/GreedyAmbiguityResolutionAlgorithm.hpp"

// #include "MyTrackFindingAlgorithm.hpp"
#include "MyTrackWriter.hpp"
#include "TString.h"
#include <filesystem>

#include "MySpacePointMaker.hpp"
#include "MyRefittingAlgorithm.hpp"
#include "MyDigitizationAlgorithm.hpp"
// #include "tracker.h"

#include "MyFtdDetector.h"

using Acts::UnitConstants::cm;
using namespace Acts::UnitConstants;
using namespace Acts::UnitLiterals;

int main(int argc, char *argv[]){
  // default parameters
  TString inputDir = "none";
  TString outputDir = "test";
//  TString outputDir = "roc_pi_16_7deg";
  int nEvents = 1000;
  //Acts::PdgParticle pdgCode = Acts::eProton;
  Acts::PdgParticle pdgCode = Acts::ePionPlus;
  double etaMin = 1.6;
  double vzMax = 50;
  double radLengthPerSeed = 0.01;

  // double ptMin = 0.950_GeV;
  // double ptMax = 1.000_GeV+1.e-100_GeV;
  // double ptMin = 0.100_GeV;
  // double ptMax = 1.000_GeV +1e-100_GeV;
  // double ptMin = 0.100_GeV;
  // double ptMax = 1.000_GeV +1e-100_GeV;
  double ptMin = 0.100_GeV;
  double ptMax = 0.101_GeV;

  // double ptMin = 0.350_GeV;
  // double ptMax = 0.350_GeV +1e-100_GeV;
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

  bool isroc = !outputDir.Contains("noroc");
  bool isframe = !outputDir.Contains("noframe");
  isroc = 0;
  isframe = 0;
  
  printf("Running acts: %d eta=%.1f vzMax=%.0f\n", pdgCode, etaMin, vzMax);

  inputDir.Append("/");
  outputDir.Append("/");
  std::filesystem::create_directory(outputDir.Data());

  // Logger
  Acts::Logging::Level logLevelF = Acts::Logging::FATAL;
  Acts::Logging::Level logLevelD = Acts::Logging::DEBUG;
  Acts::Logging::Level logLevelV = Acts::Logging::VERBOSE;
  Acts::Logging::Level logLevel = Acts::Logging::INFO;
  // collection names
  auto events = "hepmc3_event";
  auto particles = "particles";
  auto vertices = "vertices";
  auto ckf_tracks = "ckf_tracks";
  auto tracks = "tracks";
  auto simhits = "simhits";
  auto measurements = "measurements";
  auto spacepoints = "spacepoints";
  auto seeds = "seeds";
  auto estimatedparameters = "estimatedparameters";
  auto track_particle_matching = "track_particle_matching";
  auto particle_track_matching = "particle_track_matching";
  auto measurement_particles_map = "measurement_particles_map";
  auto measurement_simhits_map = "measurement_simhits_map";
  // Random number generator config
  auto rnd = std::make_shared<ActsExamples::RandomNumbers>(ActsExamples::RandomNumbers::Config({42}));

  // Particle gun generator config
  double etaMax = etaMin + 1.e-4;

  ActsExamples::ParametricParticleGenerator::Config genCfg;
  genCfg.etaUniform = true;
  genCfg.pTransverse = true;
  genCfg.pdg = pdgCode;
  // genCfg.thetaMin = 2 * atan(exp(-etaMax));
  // genCfg.thetaMax = 2 * atan(exp(-etaMin));
  // genCfg.pMin = 0.2;
  // genCfg.pMax = 1.0;
  genCfg.pMin = 0.5;
  genCfg.pMax = 0.501;
  // genCfg.phiMin = M_PI/180.*45.;
  // genCfg.phiMax = M_PI/180.*45.00001;

  // central UrQMD-like occupancy
  // genCfg.thetaMin = 2 * atan(exp(-1.95));
  // genCfg.thetaMax = 2 * atan(exp(-1.55));
  genCfg.thetaMin = 2 * atan(exp(-1.76));
  genCfg.thetaMax = 2 * atan(exp(-1.74));
  // genCfg.thetaMin = 2 * atan(exp(-1.93));
  // genCfg.thetaMax = 2 * atan(exp(-1.91));
  // genCfg.thetaMin = 2 * atan(exp(-1.59));
  // genCfg.thetaMax = 2 * atan(exp(-1.57));
  // genCfg.thetaMin = 0.35;
  // genCfg.thetaMax = 0.350001;
  // genCfg.pMin = 1.2;
  // genCfg.pMax = 1.21;
//  genCfg.randomizeCharge = true;
//  genCfg.numParticles = 90;
//  genCfg.numParticles = 10;
  genCfg.numParticles = 1;
  ActsExamples::EventGenerator::Generator gen{
      std::make_shared<ActsExamples::FixedMultiplicityGenerator>(1),
      std::make_shared<ActsExamples::FixedPrimaryVertexPositionGenerator>(),
      std::make_shared<ActsExamples::ParametricParticleGenerator>(genCfg)};

  ActsExamples::EventGenerator::Config evgenCfg;
  evgenCfg.outputEvent = events;
  evgenCfg.generators = {gen};
  evgenCfg.randomNumbers = rnd;

  ActsExamples::HepMC3InputConverter::Config hepMC3ConverterCfg;
  hepMC3ConverterCfg.inputEvent = events;
  hepMC3ConverterCfg.outputParticles = particles;
  hepMC3ConverterCfg.outputVertices = vertices;

  // Particle reader config
  ActsExamples::RootParticleReader::Config particleReaderCfg;
  particleReaderCfg.outputParticles = particles;
  particleReaderCfg.treeName = "particles";
  particleReaderCfg.filePath = TString(inputDir+"particles.root").Data();

  auto detector = std::make_shared<MyFtdDetector>();
  std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry = detector->GetTrackingGeometry(true, isroc, isframe, false);
  std::shared_ptr<MyFtdGeo> ftdGeo = detector->FtdGeo();
  const auto& surfaceByIdentifier = trackingGeometry->geoIdSurfaceMap();

  // Magnetic field
  double bz = ftdGeo->GetField() * Acts::UnitConstants::T;
  auto magField = std::make_shared<Acts::ConstantBField>(Acts::Vector3(0.0, 0.0, bz));

  // Fatras config
  ActsExamples::FatrasSimulation::Config fatrasCfg;
  fatrasCfg.inputParticles = particles;
  fatrasCfg.outputParticles = "final";
  fatrasCfg.outputSimHits = simhits;
  fatrasCfg.trackingGeometry = trackingGeometry;
  fatrasCfg.pMin = 0.1_GeV;
  fatrasCfg.magneticField = magField;
  fatrasCfg.randomNumbers = rnd;

  int shift = 2;
  if (isroc) shift++;
  if (isframe) shift++;

  // Digitization component for strips
  ActsExamples::DigiComponentsConfig stripConfig;
  stripConfig.smearingDigiConfig.params.push_back(ActsExamples::ParameterSmearingConfig{Acts::eBoundLoc0, ActsExamples::Digitization::Gauss(0.1)});
  // Digitization component for 2D measurements
  ActsExamples::DigiComponentsConfig digiConfig;
  digiConfig.smearingDigiConfig.params.push_back(ActsExamples::ParameterSmearingConfig{Acts::eBoundLoc0, ActsExamples::Digitization::Gauss(0.1)});
  digiConfig.smearingDigiConfig.params.push_back(ActsExamples::ParameterSmearingConfig{Acts::eBoundLoc1, ActsExamples::Digitization::Gauss(0.1)});

  std::vector<std::pair<Acts::GeometryIdentifier, ActsExamples::DigiComponentsConfig>> elements;

  for (auto surfId : surfaceByIdentifier) {
    auto id = surfId.first;
    const auto* detEl = dynamic_cast<const MyDetectorElement*>(surfId.second->associatedDetectorElement());
    if (detEl == nullptr)
      continue;
    if (detEl->name().find("FTD") != std::string::npos) {
      auto layer = dynamic_cast<const MyFtdDetectorElement*>(detEl)->layer();
      if (ftdGeo->GetLayerType(layer)==MyFtdGeo::FtdLayerTypes::kPixel) {
        elements.emplace_back(id, digiConfig);
      } else {
        elements.emplace_back(id, stripConfig);
      }
    }
  }

  // std::vector<std::pair<Acts::GeometryIdentifier, ActsExamples::DigiComponentsConfig>> elements = { {Acts::GeometryIdentifier{}, digiConfig} };

  // Digitization config
  // ActsExamples::MyDigitizationAlgorithm::Config digiCfg;
  ActsExamples::DigitizationAlgorithm::Config digiCfg;
  digiCfg.inputSimHits = simhits;
  digiCfg.randomNumbers = rnd;
  digiCfg.outputMeasurements = measurements;
  digiCfg.outputMeasurementParticlesMap = measurement_particles_map;
  digiCfg.outputMeasurementSimHitsMap = measurement_simhits_map;
  digiCfg.surfaceByIdentifier = surfaceByIdentifier;
  digiCfg.digitizationConfigs = Acts::GeometryHierarchyMap<ActsExamples::DigiComponentsConfig>(elements);

  // // Create space points
  // ActsExamples::SpacePointMaker::Config spCfg;
  // spCfg.inputMeasurements = measurements;
  // spCfg.trackingGeometry = trackingGeometry;
  // spCfg.outputSpacePoints = spacepoints;
  // spCfg.geometrySelection = {
  //   Acts::GeometryIdentifier().withVolume(1).withLayer( 0+shift), // 3 + shift
  //   Acts::GeometryIdentifier().withVolume(1).withLayer(18+shift), // 17 + shift
  //   Acts::GeometryIdentifier().withVolume(1).withLayer(36+shift)  // 31 + shift
  // };

  // Create space points
  ActsExamples::MySpacePointMaker::Config spCfg;
  spCfg.inputMeasurements = measurements;
  spCfg.detector = detector;
  spCfg.maxDeltaStrawId = 6;
  spCfg.minMeasPerCand = 3;
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
  int iB = 0, iM = 18, iF = 36;
  //int iB = 3, iM = 17, iF = 31;
  //int iB = 4, iM = 6, iF = 8;
  //int iB = 0, iM = 2, iF = 4;
  double rMaxStation = ftdGeo->GetLayerRMax(ftdGeo->GetNumberOfLayers()-1);
  auto positions = ftdGeo->GetLayerPositions();
  ActsExamples::SeedingAlgorithm::Config seedingCfg;
  seedingCfg.inputSpacePoints = {spacepoints};
  seedingCfg.outputSeeds = seeds;
  seedingCfg.seedFinderOptions.bFieldInZ = bz*Acts::UnitConstants::T;
  seedingCfg.seedFinderConfig.useDetailedDoubleMeasurementInfo = true;
  seedingCfg.seedFinderConfig.minPt              = 0.11_GeV;
  seedingCfg.seedFinderConfig.deltaZMax          = (positions[iM] - positions[iB])*cm + 1_mm;
  seedingCfg.seedFinderConfig.deltaRMax          = (1-positions[iM]/positions[iF])*rMaxStation*cm + 1_mm;
  seedingCfg.seedFinderConfig.zMin               = positions[iB]*cm - 1_mm;
  seedingCfg.seedFinderConfig.zMax               = positions[iF]*cm + 1_mm;
  seedingCfg.seedFinderConfig.rMin               = ftdGeo->GetLayerRMin(iB)*cm;
  seedingCfg.seedFinderConfig.rMax               = ftdGeo->GetLayerRMax(iF)*cm;
  seedingCfg.seedFinderConfig.rMinMiddle         = ftdGeo->GetLayerRMin(iM)*cm;
  seedingCfg.seedFinderConfig.rMaxMiddle         = ftdGeo->GetLayerRMax(iM)*cm;
  seedingCfg.seedFinderConfig.cotThetaMax        = 7.0;  
  seedingCfg.seedFinderConfig.impactMax          = ftdGeo->GetLayerRMin(iB)*cm - 10_mm;
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
  paramsEstimationCfg.outputTrackParameters = estimatedparameters;
  paramsEstimationCfg.trackingGeometry = trackingGeometry;
  paramsEstimationCfg.magneticField = fatrasCfg.magneticField;
  // paramsEstimationCfg.initialSigmas = {
  //       0.01 * Acts::UnitConstants::mm,
  //       0.01 * Acts::UnitConstants::mm,
  //       0.01 * Acts::UnitConstants::degree,
  //       0.01 * Acts::UnitConstants::degree,
  //       0.01 * Acts::UnitConstants::e / Acts::UnitConstants::GeV,
  //       1 * Acts::UnitConstants::ns};
  // paramsEstimationCfg.initialSigmaQoverPt = 0; // 0.1 * Acts::UnitConstants::e / Acts::UnitConstants::GeV;
  // paramsEstimationCfg.initialSigmaPtRel = 0; //0.1;

  // Measurement selector
  std::vector<std::pair<Acts::GeometryIdentifier, Acts::MeasurementSelectorCuts>> measSel;

  for (int l = 0; l < positions.size(); ++l) {
    double chi2 = ftdGeo->GetLayerType(l)==MyFtdGeo::FtdLayerTypes::kPixel ? -1 : 10; // std::numeric_limits<double>::max();
    measSel.emplace_back(detector->FtdLayerToGeoId(l).withSensitive(0), Acts::MeasurementSelectorCuts({}, {chi2}, {2u}));
  }
  // trackFindingCfg.measurementSelectorCfg = {{Acts::GeometryIdentifier(), {{}, {std::numeric_limits<double>::max()}, {1u}}}}; // chi2cut, numberOfMeasurementsPerSurface

  // Track finding
  ActsExamples::TrackFindingAlgorithm::Config trackFindingCfg;
  //trackFindingCfg.reverseSearch = true;
  trackFindingCfg.inputMeasurements = measurements;
  trackFindingCfg.inputInitialTrackParameters = paramsEstimationCfg.outputTrackParameters;
  trackFindingCfg.outputTracks = tracks;
  //trackFindingCfg.outputTracks = ckf_tracks;
  trackFindingCfg.trackingGeometry = trackingGeometry;
  trackFindingCfg.magneticField = fatrasCfg.magneticField;
  trackFindingCfg.measurementSelectorCfg = Acts::GeometryHierarchyMap(measSel);
  trackFindingCfg.findTracks = ActsExamples::TrackFindingAlgorithm::makeTrackFinderFunction(
      trackingGeometry, fatrasCfg.magneticField, *Acts::getDefaultLogger("TrackFinder", logLevelF));

  // MyTrackFindingAlgorithm::Config myTrackFindingCfg;
  // myTrackFindingCfg.inputMeasurementParticlesMap = measurement_particles_map;
  // myTrackFindingCfg.inputMeasurements = measurements;
  // myTrackFindingCfg.inputInitialTrackParameters = paramsEstimationCfg.outputTrackParameters;
  // myTrackFindingCfg.outputTracks = ckf_tracks;
  // myTrackFindingCfg.trackingGeometry = trackingGeometry;
  // myTrackFindingCfg.magneticField = fatrasCfg.magneticField;

  // ambiguity resolution
  ActsExamples::GreedyAmbiguityResolutionAlgorithm::Config ambigResCfg;
  ambigResCfg.inputTracks = ckf_tracks;
  ambigResCfg.outputTracks = tracks;
  ambigResCfg.nMeasurementsMin = 5;
  ambigResCfg.maximumSharedHits = 2;
  ambigResCfg.maximumIterations = 1000;

  // Track truth matcher
  ActsExamples::TrackTruthMatcher::Config trackTruthMatcherCfg;
  trackTruthMatcherCfg.inputTracks = tracks;
  trackTruthMatcherCfg.inputParticles = particles;
  trackTruthMatcherCfg.inputMeasurementParticlesMap = measurement_particles_map;
  trackTruthMatcherCfg.outputTrackParticleMatching = track_particle_matching;
  trackTruthMatcherCfg.outputParticleTrackMatching = particle_track_matching;
  trackTruthMatcherCfg.doubleMatching = false;
  trackTruthMatcherCfg.matchingRatio = 0.8;

  // Particle writer config
  ActsExamples::RootParticleWriter::Config particleWriterCfg;
  particleWriterCfg.inputParticles = particles;
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
  measWriterCfg.inputMeasurementSimHitsMap = measurement_simhits_map;
  measWriterCfg.surfaceByIdentifier = trackingGeometry->geoIdSurfaceMap();
  measWriterCfg.filePath = TString(outputDir+"measurements.root").Data();

  // SpacepointWriter config
  ActsExamples::RootSpacepointWriter::Config spWriterCfg;
  spWriterCfg.inputSpacepoints = spacepoints;
  spWriterCfg.filePath = TString(outputDir+"spacepoints.root").Data();

  ActsExamples::RootSeedWriter::Config seedWriterCfg;
  seedWriterCfg.inputSeeds = seeds;
  seedWriterCfg.writingMode = "";
  seedWriterCfg.filePath = TString(outputDir+"seeds.root").Data();

  // write track states from CKF
  ActsExamples::RootTrackStatesWriter::Config trackStatesWriterCfg;
  trackStatesWriterCfg.inputTracks = tracks;
  trackStatesWriterCfg.inputParticles = particles;
  trackStatesWriterCfg.inputSimHits = simhits;
  trackStatesWriterCfg.inputTrackParticleMatching = track_particle_matching;
  trackStatesWriterCfg.inputMeasurementSimHitsMap = measurement_simhits_map;
  trackStatesWriterCfg.filePath = TString(outputDir+"trackstates.root").Data();

  ActsExamples::RootTrackSummaryWriter::Config trackSummaryWriterCfg;
  trackSummaryWriterCfg.inputTracks = tracks;
  trackSummaryWriterCfg.inputParticles = particles;
  trackSummaryWriterCfg.inputTrackParticleMatching = track_particle_matching;
  trackSummaryWriterCfg.filePath = TString(outputDir+"tracksummary.root").Data();

  MyTrackWriter::Config trackWriterCfg;
  trackWriterCfg.inputTracks = tracks;
  trackWriterCfg.filePath = TString(outputDir+"tracks.root").Data();

  // Refitting algorithm config
  ActsExamples::MyRefittingAlgorithm::Config refitCfg;
  refitCfg.inputTracks = tracks;
  refitCfg.inputMeasurements = measurements;
  refitCfg.outputTracks = "refitted_tracks";
  refitCfg.trackingGeometry = trackingGeometry;
  refitCfg.magneticField = fatrasCfg.magneticField;

  ActsExamples::RootTrackSummaryWriter::Config trackRefitSummaryWriterCfg;
  trackRefitSummaryWriterCfg.inputTracks = refitCfg.outputTracks;
  trackRefitSummaryWriterCfg.inputParticles = particles;
  trackRefitSummaryWriterCfg.inputTrackParticleMatching = trackTruthMatcherCfg.outputTrackParticleMatching;
  trackRefitSummaryWriterCfg.filePath = TString(outputDir+"trackrefit.root").Data();

  // Sequencer config
  ActsExamples::Sequencer::Config sequencerCfg;
  sequencerCfg.numThreads = 1;
  sequencerCfg.events = nEvents;
  sequencerCfg.logLevel = logLevelF;

  // Start sequencer
  ActsExamples::Sequencer sequencer(sequencerCfg);
  
  if (inputDir.Contains("none")){ // particle gun + fatras simulation
    sequencer.addReader(std::make_shared<ActsExamples::EventGenerator>(evgenCfg, logLevel));
    sequencer.addAlgorithm(std::make_shared<ActsExamples::HepMC3InputConverter>(hepMC3ConverterCfg, logLevel));
    sequencer.addElement(std::make_shared<ActsExamples::FatrasSimulation>(fatrasCfg, logLevel));
  } else { // read particles and hits from input file
    sequencer.addReader(std::make_shared<ActsExamples::RootParticleReader>(particleReaderCfg, logLevel));
    sequencer.addReader(std::make_shared<ActsExamples::RootSimHitReader>(simhitReaderCfg, logLevel));
  }

  sequencer.addAlgorithm(std::make_shared<ActsExamples::DigitizationAlgorithm>(digiCfg, logLevel));
//  sequencer.addAlgorithm(std::make_shared<ActsExamples::MyDigitizationAlgorithm>(digiCfg, logLevel));
//  sequencer.addAlgorithm(std::make_shared<ActsExamples::SpacePointMaker>(spCfg, logLevel));
  sequencer.addAlgorithm(std::make_shared<ActsExamples::MySpacePointMaker>(spCfg, logLevel));

  sequencer.addAlgorithm(std::make_shared<ActsExamples::SeedingAlgorithm>(seedingCfg, logLevel));
  sequencer.addAlgorithm(std::make_shared<ActsExamples::TrackParamsEstimationAlgorithm>(paramsEstimationCfg, logLevel));

  sequencer.addAlgorithm(std::make_shared<ActsExamples::TrackFindingAlgorithm>(trackFindingCfg, logLevelF));
  // // //sequencer.addAlgorithm(std::make_shared<MyTrackFindingAlgorithm>(myTrackFindingCfg, logLevel));
  // sequencer.addAlgorithm(std::make_shared<ActsExamples::GreedyAmbiguityResolutionAlgorithm>(ambigResCfg, logLevel));  
  // // // //   // sequencer.addAlgorithm(std::make_shared<ActsExamples::MyRefittingAlgorithm>(refitCfg, logLevel));
  sequencer.addAlgorithm(std::make_shared<ActsExamples::TrackTruthMatcher>(trackTruthMatcherCfg, logLevel));

  sequencer.addWriter(std::make_shared<ActsExamples::RootParticleWriter>(particleWriterCfg, logLevel));
  sequencer.addWriter(std::make_shared<ActsExamples::RootSimHitWriter>(simhitWriterCfg, logLevel));
  sequencer.addWriter(std::make_shared<ActsExamples::RootMeasurementWriter>(measWriterCfg, logLevel));
  sequencer.addWriter(std::make_shared<ActsExamples::RootSpacepointWriter>(spWriterCfg, logLevel));

  sequencer.addWriter(std::make_shared<ActsExamples::RootSeedWriter>(seedWriterCfg, logLevel));
  sequencer.addWriter(std::make_shared<ActsExamples::RootTrackStatesWriter>(trackStatesWriterCfg, logLevel));
  sequencer.addWriter(std::make_shared<ActsExamples::RootTrackSummaryWriter>(trackSummaryWriterCfg, logLevel));
  // sequencer.addWriter(std::make_shared<MyTrackWriter>(trackWriterCfg, logLevel));
  // // // sequencer.addWriter(std::make_shared<ActsExamples::RootTrackSummaryWriter>(trackRefitSummaryWriterCfg, logLevel));

  sequencer.run();

  return 0;
}
